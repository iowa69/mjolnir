"""Coverage: mean and median depth, breadth at thresholds, evenness, GC.

The whole module exists to answer one question honestly — *was this genome
sequenced well enough for the rest of the report to mean anything* — and the
answer has to survive being wrong in either direction. A mean depth of 40x over
a genome where a third of the positions have none is not a 40x genome, which is
why breadth and evenness are computed beside the mean rather than in place of
it, and why the evenness definition travels with the number: "evenness" is not a
standard quantity and a reader who is not told the band cannot check the claim.

Depth is accumulated as a histogram rather than as a list of per-position
depths. A 4.4 Mb genome is 4.4 million integers, and the histogram gives the
same mean, the same exact median, every breadth level and the evenness fraction
from a few hundred keys. It also makes the summary functions pure: a test can
hand ``summarise`` a histogram and check the arithmetic without a BAM.

One footgun is worth naming because it silently changes every number here:
``samtools depth`` takes ``-q`` as the *base* quality threshold and ``-Q`` as the
*mapping* quality threshold, which is the reverse of ``samtools mpileup``.
Swapping them produces a plausible coverage profile computed with the wrong
filters.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from ..config import (
    COVERAGE_EVENNESS_BAND,
    Config,
    DEGRADED_DEPTH_FLOOR,
    EVENNESS_DEFINITION,
    GC_TOLERANCE,
    MIN_BREADTH,
    MIN_COVERAGE_EVENNESS,
    MIN_DEPTH,
    MIN_MAPPED_FRACTION,
    MTBSEQ_MINBQUAL,
    source_for,
)
from ..records import (
    PLATFORM_FASTA,
    PLATFORM_ILLUMINA,
    Check,
    QCMetrics,
    normalise_platform,
)
from ..utils import (
    LOG,
    MjolnirError,
    PathLike,
    require,
    require_file,
    safe_fraction,
)
from .map import (
    MIN_BASE_QUALITY,
    MIN_MAPPING_QUALITY,
    iter_output,
)

#: SOURCE: conventional coverage-breadth reporting levels, and the two the
#: record itself names (``QCMetrics.breadth_1x``, ``QCMetrics.breadth_10x``).
#: The clinically load-bearing level is neither of these: it is
#: ``config.DEGRADED_DEPTH_FLOOR``, reported as ``breadth_min_depth``.
BREADTH_LEVELS: Tuple[int, ...] = (1, 10)

#: Read flags excluded from the depth calculation — the same set the pileup
#: excludes, so breadth and allele fractions describe the same reads.
DEPTH_EXCLUDE_FLAGS = "UNMAP,SECONDARY,QCFAIL,DUP"

#: SOURCE: `samtools stats` output format. The SN keys Mjolnir reads. Kept as a
#: table so that a samtools release renaming one of them produces a missing
#: metric — reported as unmeasured — rather than a silently wrong one.
_STATS_KEYS: Dict[str, str] = {
    "raw total sequences:": "total_reads",
    "reads mapped:": "mapped_reads",
    "reads duplicated:": "duplicate_reads",
    "average length:": "mean_read_length",
    "average quality:": "mean_base_quality",
    "bases mapped (cigar):": "bases_mapped",
    "error rate:": "error_rate",
}


# ---------------------------------------------------------------------------
# The depth histogram
# ---------------------------------------------------------------------------

@dataclass
class DepthHistogram:
    """Depth -> number of reference positions at that depth.

    ``reference_length`` is what makes breadth honest. ``samtools depth`` reports
    positions; the denominator for "fraction of the genome covered" is the
    genome, and if the two differ the difference is uncovered sequence that would
    otherwise vanish from the fraction entirely.
    """

    counts: Dict[int, int] = field(default_factory=dict)
    reference_length: Optional[int] = None

    def add(self, depth: int, positions: int = 1) -> None:
        depth = int(depth)
        self.counts[depth] = self.counts.get(depth, 0) + int(positions)

    @property
    def positions_reported(self) -> int:
        return sum(self.counts.values())

    @property
    def total_positions(self) -> int:
        """Positions the fractions are computed over: the reference, when known."""
        if self.reference_length is not None:
            return max(int(self.reference_length), self.positions_reported)
        return self.positions_reported

    @property
    def _implied_zero_positions(self) -> int:
        return self.total_positions - self.positions_reported

    @property
    def total_depth(self) -> int:
        return sum(depth * count for depth, count in self.counts.items())

    def mean(self) -> Optional[float]:
        return safe_fraction(self.total_depth, self.total_positions)

    def median(self) -> Optional[float]:
        """Exact median over every reference position, uncovered ones included."""
        total = self.total_positions
        if total == 0:
            return None
        counts = dict(self.counts)
        # Positions the depth output never mentioned are at depth zero, and
        # leaving them out would report the median of the covered part of the
        # genome as if it were the median of the genome.
        zeros = self._implied_zero_positions
        if zeros:
            counts[0] = counts.get(0, 0) + zeros
        lower_index = (total - 1) // 2
        upper_index = total // 2
        lower_value: Optional[int] = None
        upper_value: Optional[int] = None
        seen = 0
        for depth, count in sorted(counts.items()):
            seen += count
            if lower_value is None and seen > lower_index:
                lower_value = depth
            if seen > upper_index:
                upper_value = depth
                break
        if lower_value is None or upper_value is None:
            return None
        return (float(lower_value) + float(upper_value)) / 2.0

    def positions_at_least(self, min_depth: int) -> int:
        return sum(count for depth, count in self.counts.items() if depth >= int(min_depth))

    def breadth(self, min_depth: int = 1) -> Optional[float]:
        """Fraction of the reference at or above *min_depth*."""
        return safe_fraction(self.positions_at_least(min_depth), self.total_positions)

    def evenness(self, band: Tuple[float, float] = COVERAGE_EVENNESS_BAND
                 ) -> Optional[float]:
        """Fraction of positions inside a band around the mean.

        Definition, which must be printed with the number: the fraction of
        reference positions whose depth lies between ``band[0]`` and ``band[1]``
        times the mean depth. There is no standard evenness statistic, so the
        number is meaningless without this sentence.
        """
        mean = self.mean()
        if not mean:
            return None
        low, high = float(band[0]) * mean, float(band[1]) * mean
        inside = sum(count for depth, count in self.counts.items()
                     if low <= depth <= high)
        if low <= 0 <= high and self._implied_zero_positions:
            inside += self._implied_zero_positions
        return safe_fraction(inside, self.total_positions)


@dataclass
class ReadStats:
    """What ``samtools stats`` says about the library, with absence preserved."""

    total_reads: Optional[int] = None
    mapped_reads: Optional[int] = None
    duplicate_reads: Optional[int] = None
    mean_read_length: Optional[float] = None
    mean_base_quality: Optional[float] = None
    bases_mapped: Optional[int] = None
    error_rate: Optional[float] = None
    gc_content: Optional[float] = None

    @property
    def mapped_fraction(self) -> Optional[float]:
        if self.total_reads is None or self.mapped_reads is None:
            return None
        return safe_fraction(self.mapped_reads, self.total_reads)

    @property
    def duplicate_fraction(self) -> Optional[float]:
        if self.total_reads is None or self.duplicate_reads is None:
            return None
        return safe_fraction(self.duplicate_reads, self.total_reads)


# ---------------------------------------------------------------------------
# Pure command builders
# ---------------------------------------------------------------------------

def samtools_depth_argv(bam: PathLike, *, platform: str = PLATFORM_ILLUMINA,
                        min_base_quality: int = MIN_BASE_QUALITY,
                        min_mapping_quality: int = MIN_MAPPING_QUALITY,
                        all_positions: bool = True,
                        regions_bed: Optional[PathLike] = None,
                        mtbseq_compat: bool = False) -> List[str]:
    """``samtools depth`` over every reference position.

    ``-q`` is the base-quality threshold and ``-Q`` the mapping-quality
    threshold: the opposite way round from ``samtools mpileup``, verified against
    samtools 1.20 and 1.24 on this machine.

    ``-s`` is passed for paired-end data so that the two mates of an overlapping
    pair are counted once. ``samtools mpileup`` removes overlaps by default, so
    without ``-s`` the depth reported here would exceed the depth the allele
    fractions were computed over, in exactly the short-insert libraries where it
    matters most. It is not passed on ONT, where there are no templates to
    overlap, nor under ``mtbseq_compat``, which reproduces MTBseq's ``-x``.

    No depth ceiling is set: modern ``samtools depth`` has none to disable. On
    samtools older than 1.13 the ``-d`` cap still existed and defaulted to 8000,
    which truncates the mean, the median and the evenness band on a deeply
    sequenced sample — a reason to record the samtools version in the report.
    """
    plat = normalise_platform(platform)
    argv = ["samtools", "depth"]
    if mtbseq_compat:
        argv += ["-q", str(MTBSEQ_MINBQUAL), "-Q", "0"]
    else:
        argv += ["-q", str(min_base_quality), "-Q", str(min_mapping_quality)]
        if plat == PLATFORM_ILLUMINA:
            argv.append("-s")
    # Restates samtools' own default filter-out list rather than changing it, so
    # that the command in the methods annex shows which reads were counted.
    argv += ["-G", DEPTH_EXCLUDE_FLAGS]
    if all_positions:
        argv.append("-a")
    if regions_bed is not None:
        argv += ["-b", str(regions_bed)]
    argv.append(str(bam))
    return argv


def samtools_stats_argv(bam: PathLike, *, reference: Optional[PathLike] = None,
                        threads: int = 1) -> List[str]:
    """``samtools stats`` — read counts, lengths, qualities and GC in one pass."""
    argv = ["samtools", "stats", "-@", str(max(1, threads))]
    if reference is not None:
        argv += ["--reference", str(reference)]
    argv.append(str(bam))
    return argv


# ---------------------------------------------------------------------------
# Pure parsing
# ---------------------------------------------------------------------------

def parse_depth_stream(lines: Iterable[str]) -> Dict[str, DepthHistogram]:
    """Accumulate ``samtools depth`` output into one histogram per contig."""
    per_contig: Dict[str, DepthHistogram] = {}
    for line in lines:
        if not line:
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        try:
            depth = int(fields[2])
        except ValueError:
            continue
        histogram = per_contig.get(fields[0])
        if histogram is None:
            histogram = DepthHistogram()
            per_contig[fields[0]] = histogram
        histogram.add(depth)
    return per_contig


def merge_histograms(per_contig: Dict[str, DepthHistogram],
                     reference_length: Optional[int] = None) -> DepthHistogram:
    """One genome-wide histogram from the per-contig ones."""
    overall = DepthHistogram(reference_length=reference_length)
    for histogram in per_contig.values():
        for depth, count in histogram.counts.items():
            overall.counts[depth] = overall.counts.get(depth, 0) + count
    return overall


def _stats_number(text: str) -> Optional[float]:
    try:
        return float(str(text).strip())
    except (TypeError, ValueError):
        return None


def parse_samtools_stats(lines: Iterable[str]) -> ReadStats:
    """Read the ``SN`` summary and the GC distributions out of ``samtools stats``.

    GC comes from the ``GCF``/``GCL`` histograms — GC percentage against read
    count, for first and last fragments — weighted by count. There is no single
    SN line carrying it, and computing it from the per-cycle nucleotide counts
    instead would weight every cycle equally regardless of how many reads reached
    it. If neither histogram is present the GC content stays ``None``: this
    parser was written against the documented format and not verified against a
    samtools binary on this machine, and a missing section must therefore read as
    unmeasured rather than as zero.
    """
    stats = ReadStats()
    gc_weighted = 0.0
    gc_reads = 0.0
    for raw in lines:
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        marker = fields[0]
        if marker == "SN" and len(fields) >= 3:
            attribute = _STATS_KEYS.get(fields[1].strip())
            if attribute is None:
                continue
            value = _stats_number(fields[2])
            if value is None:
                continue
            if attribute in ("total_reads", "mapped_reads", "duplicate_reads",
                             "bases_mapped"):
                setattr(stats, attribute, int(value))
            else:
                setattr(stats, attribute, float(value))
        elif marker in ("GCF", "GCL") and len(fields) >= 3:
            percent = _stats_number(fields[1])
            count = _stats_number(fields[2])
            if percent is None or count is None:
                continue
            gc_weighted += percent * count
            gc_reads += count
    if gc_reads:
        stats.gc_content = (gc_weighted / gc_reads) / 100.0
    return stats


def reference_length_from_fai(reference: PathLike) -> Optional[int]:
    """Total reference length from the ``.fai``, or None when it is not there.

    None rather than a guess: the breadth denominator is the whole point of the
    file, and inventing it from the depth output would make breadth 1.0 for any
    sample, however little of the genome it covered.
    """
    fai = Path(str(reference) + ".fai")
    if not fai.exists():
        return None
    total = 0
    with open(str(fai)) as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            try:
                total += int(fields[1])
            except ValueError:
                continue
    return total or None


# ---------------------------------------------------------------------------
# Rule-derived verdicts
# ---------------------------------------------------------------------------

def build_qc_checks(metrics: QCMetrics, *, config: Optional[Config] = None,
                    reference_gc: Optional[float] = None,
                    duplicates_note: str = "") -> List[Check]:
    """Turn coverage metrics into checks, each carrying the source of its threshold.

    Nothing here invents a threshold. Where config has no published number for a
    metric — the duplicate fraction is the case — no check is manufactured and
    the number is reported on its own, because a pass computed against a bound
    nobody published is exactly the kind of reassurance this project exists to
    avoid.
    """
    min_depth = config.min_depth if config is not None else MIN_DEPTH
    floor = config.degraded_depth_floor if config is not None else DEGRADED_DEPTH_FLOOR
    min_breadth = config.min_breadth if config is not None else MIN_BREADTH

    checks: List[Check] = [
        Check.numeric(
            "mean_depth", metrics.mean_depth,
            minimum=floor, warn_minimum=min_depth,
            source=source_for("min_depth"), unit="x", category="qc",
            reading="target depth {0}x; below {1}x the sample is not callable"
                    .format(min_depth, floor),
            not_measured_why="no alignment was produced, so depth does not exist"),
        Check.numeric(
            "breadth_at_{0}x".format(floor), metrics.breadth_min_depth,
            warn_minimum=min_breadth,
            source=source_for("min_breadth"), unit="fraction", category="qc",
            reading="fraction of the reference reaching the {0}x floor".format(floor),
            not_measured_why="no alignment was produced, so breadth does not exist"),
        Check.numeric(
            "mapped_fraction", metrics.mapped_fraction,
            warn_minimum=MIN_MAPPED_FRACTION,
            source=source_for("min_mapped_fraction"), unit="fraction",
            category="qc",
            reading="reads mapping to the chosen reference; a low value is a "
                    "flag on the library, not a species call",
            not_measured_why="the read total was not available, so the mapped "
                             "fraction has no denominator"),
        Check.numeric(
            "coverage_evenness", metrics.coverage_evenness,
            warn_minimum=MIN_COVERAGE_EVENNESS,
            source=source_for("min_coverage_evenness"), unit="fraction",
            category="qc",
            reading=metrics.evenness_definition or EVENNESS_DEFINITION,
            not_measured_why="evenness needs a mean depth to define its band"),
    ]

    if metrics.gc_content is None:
        checks.append(Check.not_measured(
            "gc_content", "GC content was not measured for this sample",
            source=source_for("gc_tolerance"), category="qc"))
    elif reference_gc is None:
        checks.append(Check.not_measured(
            "gc_content_vs_reference",
            "observed GC is {0:.3f}, but this reference's own GC content is not "
            "recorded in the database registry, so there is no band to compare "
            "it against".format(metrics.gc_content),
            source=source_for("gc_tolerance"), category="qc"))
    else:
        checks.append(Check.numeric(
            "gc_content", metrics.gc_content,
            warn_minimum=float(reference_gc) - GC_TOLERANCE,
            warn_maximum=float(reference_gc) + GC_TOLERANCE,
            source=source_for("gc_tolerance"), unit="fraction", category="qc",
            reading="within {0} of the reference's {1:.3f}".format(
                GC_TOLERANCE, float(reference_gc))))

    if metrics.duplicate_fraction is None:
        checks.append(Check.not_measured(
            "duplicate_fraction",
            duplicates_note or "duplicates were not marked, so the duplicate "
                               "fraction is unknown rather than zero",
            category="qc"))
    return checks


# ---------------------------------------------------------------------------
# The wrapper that runs it
# ---------------------------------------------------------------------------

def measure_coverage(bam: PathLike, reference: PathLike, *,
                     platform: str = PLATFORM_ILLUMINA,
                     config: Optional[Config] = None,
                     reference_length: Optional[int] = None,
                     reference_name: str = "",
                     reference_gc: Optional[float] = None,
                     duplicates_marked: bool = False,
                     duplicates_note: str = "",
                     threads: int = 1) -> QCMetrics:
    """Measure coverage and library composition for one alignment.

    Two samtools passes: ``depth`` for the coverage profile and ``stats`` for the
    read counts, lengths, qualities and GC. Both are streamed, so a 4.4 Mb genome
    never materialises as a file of per-position depths.

    *duplicates_marked* has to be passed rather than inferred. ``samtools stats``
    counts reads carrying the duplicate flag, and on a BAM that was never through
    ``markdup`` that count is zero — a zero duplicate fraction that means "nobody
    looked", printed in the same column as one that means "a clean library". The
    fraction is therefore reported only when something actually marked them.
    """
    plat = normalise_platform(platform)
    if plat == PLATFORM_FASTA:
        raise MjolnirError(
            "an assembly has no read coverage. FASTA input carries no depth, no "
            "mapped fraction and no allele fractions, and the report states that "
            "as a capability loss (design §7)")
    require("samtools", "coverage and library metrics")
    require_file(bam, "alignment BAM")
    require_file(reference, "reference FASTA")
    if config is not None and threads == 1:
        threads = config.threads

    if reference_length is None:
        reference_length = reference_length_from_fai(reference)

    mtbseq_compat = bool(config.mtbseq_compat) if config is not None else False
    per_contig = parse_depth_stream(iter_output(
        samtools_depth_argv(bam, platform=plat, mtbseq_compat=mtbseq_compat)))
    histogram = merge_histograms(per_contig, reference_length=reference_length)
    stats = parse_samtools_stats(
        iter_output(samtools_stats_argv(bam, reference=reference, threads=threads)))

    floor = config.degraded_depth_floor if config is not None else DEGRADED_DEPTH_FLOOR
    metrics = QCMetrics(
        mean_depth=histogram.mean(),
        median_depth=histogram.median(),
        breadth_1x=histogram.breadth(BREADTH_LEVELS[0]),
        breadth_10x=histogram.breadth(BREADTH_LEVELS[1]),
        breadth_min_depth=histogram.breadth(floor),
        coverage_evenness=histogram.evenness(),
        evenness_definition=EVENNESS_DEFINITION,
        mapped_fraction=stats.mapped_fraction,
        gc_content=stats.gc_content,
        total_reads=stats.total_reads,
        mapped_reads=stats.mapped_reads,
        duplicate_fraction=stats.duplicate_fraction if duplicates_marked else None,
        mean_read_length=stats.mean_read_length,
        mean_base_quality=stats.mean_base_quality,
        reference=reference_name or Path(str(reference)).name,
        reference_length=reference_length,
    )
    metrics.checks = build_qc_checks(
        metrics, config=config, reference_gc=reference_gc,
        duplicates_note=duplicates_note or (
            "" if duplicates_marked else
            "duplicates were never marked on this alignment, so the duplicate "
            "fraction is unknown; samtools would report it as zero"))
    if reference_length is None:
        LOG.warning("no .fai beside %s, so breadth is computed over the positions "
                    "samtools reported rather than over the reference length",
                    reference)
    return metrics
