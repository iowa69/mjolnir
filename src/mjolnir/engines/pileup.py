"""Direct pileup at named positions — per-allele depths and allele fractions.

This is deliberately not the variant caller, and the difference is the point.

A caller answers "is there a variant here", and it answers it by first deciding
whether the evidence clears its own model. Two of the questions Mjolnir has to
answer are not of that shape. Lineage barcode genotyping asks "which base is at
this position, and how well is it seen" — pathogen-profiler's approach, adopted
here (design §6) — and a catalogue lookup asks "what fraction of the reads carry
the catalogued allele at this exact coordinate". Both of those have answers at
positions where no caller would emit a record, and on ONT the gap between the
two is measurable: 26 of 27 Illumina-only minor variants in the 508-isolate
comparison were *visible in the ONT pileup* and not called.

So the pileup is read directly, allele by allele, and a position that no read
reached comes back as a site with ``covered=False`` — not as a reference call.
An uncovered barcode site is missing evidence about the lineage, and reporting
it as the reference base would silently manufacture support for whichever
lineage the reference happens to belong to.

Two counting conventions are implemented. ``acgt`` is Mjolnir's: the denominator
is the ACGT depth. ``mtbseq`` reproduces design §9b — MTBseq's frequency
denominator includes N and GAP counts, and GAP wins every tie in the order
A < C < G < T < N < GAP — so that a lab reconciling the two tools can see the
same numbers rather than argue about them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from ..config import (
    Config,
    MIN_MINOR_VARIANT_FRACTION,
    MTBSEQ_MINBQUAL,
    MTBSEQ_UNAMBIG,
    source_for,
)
from ..records import (
    PLATFORM_FASTA,
    PLATFORM_ILLUMINA,
    PLATFORM_ONT,
    normalise_platform,
)
from ..utils import (
    LOG,
    MjolnirError,
    PathLike,
    ensure_dir,
    require,
    require_file,
    safe_fraction,
    tempdir,
)
from .map import (
    MAX_PILEUP_DEPTH,
    MIN_BASE_QUALITY,
    MIN_MAPPING_QUALITY,
    iter_output,
)

#: The two counting conventions. ``acgt`` is what Mjolnir reports; ``mtbseq``
#: exists so ``--compat mtbseq`` can reproduce the legacy numbers exactly.
CONVENTION_ACGT = "acgt"
CONVENTION_MTBSEQ = "mtbseq"
CONVENTIONS: Tuple[str, ...] = (CONVENTION_ACGT, CONVENTION_MTBSEQ)

#: Symbol for a position deleted in a read. samtools writes ``*`` for it, and
#: ``#`` when ``--reverse-del`` is in force; both are the same event.
GAP = "*"

#: SOURCE: MTBseq v1.1.0 source, ``TBtools::call_variants`` (design §9b). The
#: allele ordering MTBseq breaks ties with — GAP beats every base, N beats every
#: real base. Reproduced rather than approved: it is why an MTBseq call at a
#: position with equal A and GAP support is a gap and Mjolnir's is ambiguous.
MTBSEQ_TIE_BREAK_ORDER: Tuple[str, ...] = ("A", "C", "G", "T", "N", GAP)

#: Read flags excluded from every pileup. These are samtools' own defaults,
#: written out because a pileup that silently counted duplicates or secondary
#: alignments would inflate every allele fraction in the report.
PILEUP_EXCLUDE_FLAGS = "UNMAP,SECONDARY,QCFAIL,DUP"

_BASES: Tuple[str, ...] = ("A", "C", "G", "T", "N", GAP)


# ---------------------------------------------------------------------------
# Records
# ---------------------------------------------------------------------------

@dataclass
class PileupSite:
    """What the reads say at one reference position.

    ``covered`` distinguishes "no read reached this position" from "reads reached
    it and none carried the allele asked about". Those are different pieces of
    evidence and the fifth house rule turns on keeping them apart: every fraction
    on an uncovered site is ``None``, never zero.
    """

    chrom: str
    pos: int
    ref_base: str = "N"
    counts: Dict[str, int] = field(default_factory=dict)
    #: Inserted sequence -> reads carrying that insertion, starting after this base.
    insertions: Dict[str, int] = field(default_factory=dict)
    #: Deleted sequence -> reads carrying that deletion, starting after this base.
    deletions: Dict[str, int] = field(default_factory=dict)
    #: samtools' own depth column, before Mjolnir's per-base tokenising. Kept
    #: because it counts things the base column does not, and a disagreement
    #: between the two is worth being able to see.
    raw_depth: Optional[int] = None
    covered: bool = True

    def count(self, allele: str) -> int:
        return int(self.counts.get(str(allele).upper(), 0))

    @property
    def acgt_depth(self) -> int:
        return sum(int(self.counts.get(base, 0)) for base in ("A", "C", "G", "T"))

    @property
    def gap_depth(self) -> int:
        return int(self.counts.get(GAP, 0))

    @property
    def n_depth(self) -> int:
        return int(self.counts.get("N", 0))

    def denominator(self, convention: str = CONVENTION_ACGT) -> int:
        """Total depth under the requested counting convention."""
        if convention == CONVENTION_ACGT:
            return self.acgt_depth
        if convention == CONVENTION_MTBSEQ:
            return self.acgt_depth + self.n_depth + self.gap_depth
        raise MjolnirError(
            "unknown pileup convention {0!r}; expected one of {1}".format(
                convention, ", ".join(CONVENTIONS)))

    def fraction(self, allele: str, convention: str = CONVENTION_ACGT) -> Optional[float]:
        """Fraction of reads carrying *allele*, or None where nothing was seen."""
        return safe_fraction(self.count(allele), self.denominator(convention))

    def insertion_count(self, inserted: str) -> int:
        return int(self.insertions.get(str(inserted).upper(), 0))

    def deletion_count(self, deleted: str) -> int:
        return int(self.deletions.get(str(deleted).upper(), 0))

    def major_allele(self, convention: str = CONVENTION_ACGT) -> Optional[str]:
        """The best-supported allele, or None when there is no single best one.

        Under the ``acgt`` convention an exact tie returns None: two alleles at
        50% each is a real observation and picking one of them by alphabet would
        turn a mixture into a confident genotype. Under ``mtbseq`` the legacy
        tie-break is applied instead, because reproducing that tool means
        reproducing its arbitrariness too.
        """
        pool = _BASES if convention == CONVENTION_MTBSEQ else ("A", "C", "G", "T")
        best: Optional[str] = None
        best_count = 0
        tied = False
        for base in pool:
            count = int(self.counts.get(base, 0))
            if count > best_count:
                best, best_count, tied = base, count, False
            elif count == best_count and count > 0 and base != best:
                if convention == CONVENTION_MTBSEQ:
                    # GAP > N > T > G > C > A, so the later symbol wins.
                    if MTBSEQ_TIE_BREAK_ORDER.index(base) > MTBSEQ_TIE_BREAK_ORDER.index(best):
                        best = base
                else:
                    tied = True
        if best_count == 0 or tied:
            return None
        return best

    def unambiguous_fraction(self,
                             convention: str = CONVENTION_ACGT) -> Optional[float]:
        """Support for the majority allele — MTBseq's de-facto heterozygosity filter.

        Surfaced rather than applied: MTBseq discards a position below 95% and
        says nothing about it, which throws away exactly the minority signal a
        mixed-infection question is asking about (design §8.4).
        """
        best = self.major_allele(convention)
        if best is None:
            return None
        return self.fraction(best, convention)

    def is_unambiguous(self, threshold_percent: float = MTBSEQ_UNAMBIG,
                       convention: str = CONVENTION_ACGT) -> Optional[bool]:
        fraction = self.unambiguous_fraction(convention)
        if fraction is None:
            return None
        return fraction >= (float(threshold_percent) / 100.0)

    @property
    def key(self) -> Tuple[str, int]:
        return (self.chrom, self.pos)

    def to_dict(self) -> Dict[str, object]:
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref_base": self.ref_base,
            "counts": dict(self.counts),
            "insertions": dict(self.insertions),
            "deletions": dict(self.deletions),
            "depth": self.acgt_depth,
            "raw_depth": self.raw_depth,
            "covered": self.covered,
        }


@dataclass
class SiteGenotype:
    """One genotyping decision at one position, with why it was or was not made."""

    chrom: str
    pos: int
    allele: Optional[str] = None
    fraction: Optional[float] = None
    depth: Optional[int] = None
    supported: bool = False
    method: str = ""
    reason: str = ""

    @property
    def key(self) -> Tuple[str, int]:
        return (self.chrom, self.pos)


def uncovered_site(chrom: str, pos: int, ref_base: str = "N") -> PileupSite:
    """A site no read reached. Every count zero, every fraction ``None``."""
    return PileupSite(chrom=chrom, pos=int(pos), ref_base=ref_base,
                      counts=dict((base, 0) for base in _BASES),
                      raw_depth=0, covered=False)


# ---------------------------------------------------------------------------
# Pure parsing
# ---------------------------------------------------------------------------

def parse_pileup_bases(bases: str, ref_base: str = "N"
                       ) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]:
    """Tokenise samtools' base column into per-allele counts.

    The column is a small language, not a string of bases: ``^`` introduces a
    read and swallows the next character as its mapping quality, ``$`` ends one,
    ``.``/``,`` mean the reference base on either strand, ``*``/``#`` mean a
    position deleted by an indel that started earlier, and ``+n<seq>``/``-n<seq>``
    hang an indel off the base just counted. Reading it with a per-character loop
    that does not know about ``^`` mis-counts the mapping-quality character as an
    allele, which shows up as a sprinkle of spurious minority bases at exactly
    the positions where reads start — i.e. everywhere, at low coverage.
    """
    counts: Dict[str, int] = dict((base, 0) for base in _BASES)
    insertions: Dict[str, int] = {}
    deletions: Dict[str, int] = {}
    reference = (str(ref_base) or "N").upper()
    if reference not in counts:
        counts[reference] = 0

    text = str(bases or "")
    index = 0
    length = len(text)
    while index < length:
        char = text[index]
        if char == "^":
            index += 2          # the read-start marker plus its mapping quality
            continue
        if char == "$":
            index += 1
            continue
        if char in "+-":
            cursor = index + 1
            digits = ""
            while cursor < length and text[cursor].isdigit():
                digits += text[cursor]
                cursor += 1
            if not digits:
                # Not an indel marker at all; nothing sane to count.
                index += 1
                continue
            size = int(digits)
            sequence = text[cursor:cursor + size].upper()
            target = insertions if char == "+" else deletions
            target[sequence] = target.get(sequence, 0) + 1
            index = cursor + size
            continue
        if char in ".,":
            counts[reference] = counts.get(reference, 0) + 1
            index += 1
            continue
        upper = char.upper()
        if upper in ("A", "C", "G", "T", "N"):
            counts[upper] = counts.get(upper, 0) + 1
            index += 1
            continue
        if char in "*#":
            counts[GAP] = counts.get(GAP, 0) + 1
            index += 1
            continue
        # '<' and '>' are reference skips, and anything else is a symbol this
        # version of samtools invented after this code was written. Neither is
        # an allele, so neither is counted — and neither is silently treated as
        # one.
        index += 1
    return counts, insertions, deletions


def parse_pileup_line(line: str) -> Optional[PileupSite]:
    """One ``samtools mpileup`` row. Returns None for a blank or malformed row."""
    if not line or not line.strip():
        return None
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 4:
        return None
    chrom, pos_text, ref_base, depth_text = fields[:4]
    try:
        pos = int(pos_text)
    except ValueError:
        return None
    try:
        raw_depth = int(depth_text)
    except ValueError:
        raw_depth = None

    if raw_depth == 0:
        # ``samtools mpileup -a`` writes a literal '*' in the base and quality
        # columns of a zero-depth row as a placeholder. Tokenising it would count
        # a deletion at every uncovered position, which under the MTBseq
        # convention — where GAP is part of the denominator and wins ties — turns
        # an uncovered site into a confident gap call.
        return uncovered_site(chrom, pos, str(ref_base).upper())

    bases = fields[4] if len(fields) > 4 else ""
    counts, insertions, deletions = parse_pileup_bases(bases, ref_base)
    site = PileupSite(chrom=chrom, pos=pos, ref_base=str(ref_base).upper(),
                      counts=counts, insertions=insertions, deletions=deletions,
                      raw_depth=raw_depth,
                      covered=bool(raw_depth) or bool(bases.strip("*#")))
    return site


def parse_pileup(lines: Iterable[str]) -> Iterator[PileupSite]:
    for line in lines:
        site = parse_pileup_line(line)
        if site is not None:
            yield site


# ---------------------------------------------------------------------------
# Pure command building
# ---------------------------------------------------------------------------

def write_positions_file(positions: Iterable[Tuple[str, int]], path: PathLike) -> Path:
    """Write samtools' two-column position list, 1-based.

    Two columns on purpose. ``samtools mpileup -l`` reads a file with three or
    more columns as a BED — half-open and 0-based — and a two-column file as a
    1-based position list. A barcode site written into the wrong one of those is
    off by one, which does not fail: it genotypes the neighbouring base and
    reports a confident lineage from it.
    """
    target = Path(path)
    ensure_dir(target.parent)
    with open(str(target), "w") as handle:
        for chrom, pos in positions:
            handle.write("{0}\t{1}\n".format(chrom, int(pos)))
    return target


def samtools_mpileup_argv(bam: PathLike, reference: PathLike, *,
                          positions_file: Optional[PathLike] = None,
                          platform: str = PLATFORM_ILLUMINA,
                          min_base_quality: int = MIN_BASE_QUALITY,
                          min_mapping_quality: int = MIN_MAPPING_QUALITY,
                          max_depth: int = MAX_PILEUP_DEPTH,
                          all_positions: bool = True,
                          mtbseq_compat: bool = False) -> List[str]:
    """``samtools mpileup`` for the text pileup, not for calling.

    ``-a`` keeps positions with zero depth in the output, which is what makes an
    uncovered barcode site visible as an uncovered site instead of as a missing
    line somebody downstream has to notice.

    ``mtbseq_compat`` is design §9b's legacy stack: ``-B -A -x``, no MAPQ filter,
    base quality 13, and the 250x cap MTBseq inherits by never passing ``-d``.
    """
    plat = normalise_platform(platform)
    if plat == PLATFORM_FASTA:
        raise MjolnirError(
            "an assembly has no reads to pile up; allele fractions do not exist "
            "for FASTA input (design §7)")

    argv = ["samtools", "mpileup", "-f", str(reference),
            "--ff", PILEUP_EXCLUDE_FLAGS]
    if all_positions:
        argv.append("-a")
    if positions_file is not None:
        argv += ["-l", str(positions_file)]

    if mtbseq_compat:
        argv += ["-B", "-A", "-x", "-q", "0", "-Q", str(MTBSEQ_MINBQUAL), "-d", "250"]
    else:
        argv += ["-q", str(min_mapping_quality), "-Q", str(min_base_quality),
                 "-d", str(max_depth)]
        if plat == PLATFORM_ONT:
            # BAQ is an Illumina-shaped correction; on ONT it removes real
            # coverage around every homopolymer, which is where the pncA and
            # rrs questions live.
            argv.append("-B")

    argv.append(str(bam))
    return argv


# ---------------------------------------------------------------------------
# Genotyping decisions
# ---------------------------------------------------------------------------

def genotype_site(site: PileupSite, *, platform: str,
                  min_reads: int, min_fraction: Optional[float] = None,
                  convention: str = CONVENTION_ACGT) -> SiteGenotype:
    """Decide the allele at one site, by the platform's own rule.

    On ONT the highest-depth allele is taken (design §6), with no fraction
    threshold, because ONT's per-read error rate spreads support across
    neighbouring alleles and a fraction cut-off there would refuse to genotype
    perfectly good sites. On Illumina the majority allele must also clear
    *min_fraction*, so a genuinely mixed position is reported as unsupported
    rather than resolved to whichever strain happened to be commoner.

    Either way, *min_reads* is the platform threshold from config, and a site
    that fails it comes back ``supported=False`` with the reason attached — never
    as a reference call.
    """
    plat = normalise_platform(platform)
    depth = site.denominator(convention)
    if not site.covered or depth == 0:
        return SiteGenotype(chrom=site.chrom, pos=site.pos, depth=0, supported=False,
                            method="pileup", allele=None, fraction=None,
                            reason="no reads covered this position")

    allele = site.major_allele(convention)
    if allele is None:
        return SiteGenotype(chrom=site.chrom, pos=site.pos, depth=depth,
                            supported=False, method="pileup",
                            reason="no single best-supported allele: the position "
                                   "is an exact tie between alleles")

    fraction = site.fraction(allele, convention)
    counts = site.count(allele)
    method = ("highest-depth allele (ONT; design §6)" if plat == PLATFORM_ONT
              else "majority allele")

    if counts < min_reads:
        return SiteGenotype(chrom=site.chrom, pos=site.pos, allele=allele,
                            fraction=fraction, depth=depth, supported=False,
                            method=method,
                            reason="allele seen on {0} reads, below the {1} "
                                   "threshold of {2} ({3})".format(
                                       counts, plat, min_reads,
                                       source_for("min_reads_illumina"
                                                  if plat == PLATFORM_ILLUMINA
                                                  else "min_reads_ont")))

    if plat != PLATFORM_ONT and min_fraction is not None and fraction is not None \
            and fraction < min_fraction:
        return SiteGenotype(chrom=site.chrom, pos=site.pos, allele=allele,
                            fraction=fraction, depth=depth, supported=False,
                            method=method,
                            reason="majority allele at {0:.2f} of reads, below the "
                                   "{1:.2f} required to genotype a site".format(
                                       fraction, min_fraction))

    return SiteGenotype(chrom=site.chrom, pos=site.pos, allele=allele,
                        fraction=fraction, depth=depth, supported=True,
                        method=method, reason="")


def alleles_present(site: PileupSite, *, min_fraction: float = MIN_MINOR_VARIANT_FRACTION,
                    min_reads: int = 1,
                    convention: str = CONVENTION_ACGT) -> List[Tuple[str, float, int]]:
    """Every allele above a floor, commonest first — the mixed-lineage input.

    A barcode site carrying two lineage alleles at 60/40 is the observation that
    a mixed infection is made of, and it is invisible to any function that
    returns only the winner.
    """
    depth = site.denominator(convention)
    out: List[Tuple[str, float, int]] = []
    if depth == 0:
        return out
    pool = _BASES if convention == CONVENTION_MTBSEQ else ("A", "C", "G", "T")
    for base in pool:
        count = site.count(base)
        if count < min_reads:
            continue
        fraction = safe_fraction(count, depth)
        if fraction is None or fraction < min_fraction:
            continue
        out.append((base, fraction, count))
    out.sort(key=lambda item: (-item[2], item[0]))
    return out


def allele_fraction_at(site: PileupSite, ref: str, alt: str, *,
                       convention: str = CONVENTION_ACGT) -> Optional[float]:
    """The fraction of reads carrying a catalogued (REF, ALT) allele.

    Handles the three shapes a catalogue entry takes at a coordinate: a
    substitution, an insertion (ALT extends REF) and a deletion (REF extends
    ALT). Anything more complex — a genuine MNV or a replacement — returns None
    rather than a number computed from the first base, because a fraction that
    describes part of an allele is worse than no fraction at all.
    """
    ref = str(ref or "").upper()
    alt = str(alt or "").upper()
    if not site.covered:
        return None
    if len(ref) == 1 and len(alt) == 1:
        return site.fraction(alt, convention)
    denominator = site.denominator(convention)
    if len(alt) > len(ref) and alt.startswith(ref):
        return safe_fraction(site.insertion_count(alt[len(ref):]), denominator)
    if len(ref) > len(alt) and ref.startswith(alt):
        return safe_fraction(site.deletion_count(ref[len(alt):]), denominator)
    return None


# ---------------------------------------------------------------------------
# The wrapper that runs it
# ---------------------------------------------------------------------------

def pileup_at(bam: PathLike, reference: PathLike,
              positions: Sequence[Tuple[str, int]], *,
              platform: str = PLATFORM_ILLUMINA,
              config: Optional[Config] = None,
              min_base_quality: Optional[int] = None,
              min_mapping_quality: Optional[int] = None,
              max_depth: int = MAX_PILEUP_DEPTH,
              scratch_dir: Optional[PathLike] = None) -> Dict[Tuple[str, int], PileupSite]:
    """Pile up *bam* at exactly *positions*, returning one site per position.

    Every requested position is present in the result. Positions samtools did not
    report come back as :func:`uncovered_site`, so a caller iterating the barcode
    scheme cannot mistake a missing key for a reference call.
    """
    plat = normalise_platform(platform)
    wanted = [(str(chrom), int(pos)) for chrom, pos in positions]
    if not wanted:
        return {}
    require("samtools", "direct pileup at catalogue and barcode positions")
    require_file(bam, "alignment BAM")
    require_file(reference, "reference FASTA")

    mtbseq_compat = bool(config.mtbseq_compat) if config is not None else False
    argv_kwargs = dict(
        platform=plat,
        max_depth=max_depth,
        mtbseq_compat=mtbseq_compat,
    )
    if min_base_quality is not None:
        argv_kwargs["min_base_quality"] = min_base_quality
    if min_mapping_quality is not None:
        argv_kwargs["min_mapping_quality"] = min_mapping_quality

    found: Dict[Tuple[str, int], PileupSite] = {}
    with tempdir(prefix="mjolnir.pileup.", parent=scratch_dir) as scratch:
        positions_file = write_positions_file(wanted, scratch / "positions.txt")
        argv = samtools_mpileup_argv(bam, reference, positions_file=positions_file,
                                     **argv_kwargs)
        for site in parse_pileup(iter_output(argv)):
            found[site.key] = site

    missing = 0
    result: Dict[Tuple[str, int], PileupSite] = {}
    for key in wanted:
        site = found.get(key)
        if site is None:
            missing += 1
            site = uncovered_site(key[0], key[1])
        result[key] = site
    if missing:
        LOG.debug("%d of %d requested positions were absent from the pileup and "
                  "are reported as uncovered", missing, len(wanted))
    return result


def genotype_positions(bam: PathLike, reference: PathLike,
                       positions: Sequence[Tuple[str, int]], *,
                       platform: str, config: Optional[Config] = None,
                       min_reads: Optional[int] = None,
                       min_fraction: Optional[float] = None,
                       scratch_dir: Optional[PathLike] = None
                       ) -> Dict[Tuple[str, int], SiteGenotype]:
    """Pile up and genotype in one pass — what ``typing/lineage.py`` calls."""
    plat = normalise_platform(platform)
    if min_reads is None:
        min_reads = config.min_reads(plat) if config is not None else 1
    sites = pileup_at(bam, reference, positions, platform=plat, config=config,
                      scratch_dir=scratch_dir)
    convention = (CONVENTION_MTBSEQ
                  if config is not None and config.mtbseq_compat
                  else CONVENTION_ACGT)
    return dict(
        (key, genotype_site(site, platform=plat, min_reads=min_reads,
                            min_fraction=min_fraction, convention=convention))
        for key, site in sites.items()
    )
