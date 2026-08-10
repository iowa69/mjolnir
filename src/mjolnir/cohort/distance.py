"""Masked pairwise SNP distance, and the denominator that makes it a statement.

Two rules from the design are enforced here rather than recommended.

**Masking is mandatory.** The 508-isolate ONT-vs-Illumina comparison masked
264,525 loci — about 6% of H37Rv — covering repetitive, low-complexity and
error-prone regions, and counted only SNPs with no other SNP within 12 bases.
Distances computed without that mask are inflated in exactly the regions where
mapping is least trustworthy, so :func:`pairwise_distance` will not accept a
missing mask. It takes a :class:`Mask`, and the only way to compare without one
is :meth:`Mask.absent`, which is a named constructor that stamps its reason onto
every distance it produces and onto the checks. A caller cannot end up unmasked
by omission — only on purpose, in writing.

The mask is also not a constant. tbdb's changed twice in three years, the
candidate schemes differ substantially, and Marin et al. found that raising the
mapping-quality threshold outperformed blanket masking. So the mask is a named,
versioned, swappable input, and its name and size travel with every number
computed under it.

**A distance without its denominator is not returned.** Twelve differences over
4.1 Mb of shared callable sequence and twelve over 400 kb are not the same
statement. Every function here returns
:class:`~mjolnir.records.PairwiseDistance` objects, which carry
``shared_callable_sites`` beside ``snps``; there is no function in this module
that hands back a bare integer, and :class:`DistanceMatrix` deliberately has no
``distance(a, b) -> int`` accessor for a caller to reach for. Where the
denominator cannot be established, ``snps`` is ``None`` and the note says what
is missing — never a number that looks computed.

Nothing here runs an external tool. It is arithmetic over the joint table.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ..config import (MASKED_FRACTION_H37RV, MASKED_LOCI_H37RV,
                      MIN_SHARED_CALLABLE_SITES, SNP_PROXIMITY_WINDOW,
                      source_for)
from ..records import (Check, PairwiseDistance, STATUS_FAIL, STATUS_PASS,
                       STATUS_WARN, pair_key)
from ..utils import LOG, MjolnirError, PathLike, sha256sum
from .joint import AMBIGUOUS, JointTable, Regions

#: What the report must say when a comparison was made without a mask. It names
#: the size of the mask that was not applied, because "unmasked" on its own does
#: not convey that ~6% of the genome — the least trustworthy 6% — is back in the
#: count.
MASK_ABSENT_TEXT = (
    "these distances were computed WITHOUT a repeat/low-complexity mask. The "
    "508-isolate MTBC comparison masked {0:,} loci (~{1:.0%} of H37Rv) covering "
    "repetitive, low-complexity and error-prone regions before counting SNPs; "
    "unmasked distances are inflated in precisely those regions and are not "
    "comparable to the published 5-SNP and 12-SNP thresholds."
).format(MASKED_LOCI_H37RV, MASKED_FRACTION_H37RV)

#: How a pair is described when one of its members has no callable regions.
NO_DENOMINATOR_TEXT = (
    "not computed: the shared callable region of this pair is unknown, so a "
    "distance would have no denominator. Supply callable regions for {samples} "
    "(the depth engine writes them) and the pair becomes comparable."
)


# ---------------------------------------------------------------------------
# The mask
# ---------------------------------------------------------------------------

@dataclass
class Mask:
    """A named, versioned set of positions excluded before counting.

    ``applied`` False is the deliberate opt-out built by :meth:`absent`. It is a
    separate state from an empty mask: a mask file that happens to contain no
    intervals is a broken file and :meth:`Regions.from_bed` raises on it, while
    an acknowledged unmasked comparison is a decision somebody made and the
    report prints the sentence that says so.
    """

    name: str
    regions: Optional[Regions] = None
    path: str = ""
    checksum: str = ""
    source: str = ""
    applied: bool = True
    #: Why there is no mask, when ``applied`` is False. Never empty in that state.
    reason: str = ""

    def __post_init__(self) -> None:
        if self.applied and self.regions is None:
            raise MjolnirError(
                "Mask {0!r} claims to be applied but carries no intervals; build an "
                "acknowledged unmasked comparison with Mask.absent(reason=...) "
                "instead".format(self.name))
        if not self.applied and not self.reason:
            raise MjolnirError(
                "an unmasked comparison must state why: Mask.absent(reason=...)")

    @classmethod
    def from_bed(cls, path: PathLike, name: str = "", source: str = "",
                 checksum: bool = True) -> "Mask":
        """Load ``mask.bed`` in tbdb's format: BED3+, 0-based half-open.

        tbdb ships three columns and, in later releases, a fourth naming the
        region — the loader keeps the first three and ignores the rest, so both
        shapes load. The checksum is taken by default because two installations
        with differently-versioned masks produce different clusters from
        identical reads, and the report has to be able to tell them apart.
        """
        regions = Regions.from_bed(
            path, name=name,
            fetch_hint="fetch it with: mjolnir db fetch tbdb  (supplies mask.bed and "
                       "barcode.bed); for a non-MTBC reference, supply --mask-bed or "
                       "state explicitly that no mask exists")
        digest = sha256sum(path) if checksum else ""
        mask = cls(name=name or str(path), regions=regions, path=str(path),
                   checksum=digest, source=source or source_for("masked_loci_h37rv"))
        LOG.info("mask %s: %d intervals, %d bp masked",
                 mask.name, sum(len(v) for v in regions.intervals.values()),
                 regions.length())
        return mask

    @classmethod
    def absent(cls, reason: str, name: str = "none") -> "Mask":
        """An acknowledged unmasked comparison.

        The only route to counting SNPs without a mask, and it demands a reason
        that the report prints. Use it when a non-MTBC reference has no repeat
        mask yet — and expect the distances to read as provisional, because they
        are.
        """
        return cls(name=name, regions=None, applied=False, reason=reason)

    def masked_bases(self) -> Optional[int]:
        return None if self.regions is None else self.regions.length()

    def fraction_of(self, reference_length: Optional[int]) -> Optional[float]:
        if self.regions is None or not reference_length:
            return None
        return float(self.regions.length()) / float(reference_length)

    def contains(self, chrom: str, pos: int) -> bool:
        """Whether a 1-based position is masked. Always False when not applied."""
        if not self.applied or self.regions is None:
            return False
        return self.regions.contains(chrom, pos)

    def describe(self) -> str:
        if not self.applied:
            return "no mask applied ({0})".format(self.reason)
        return "{0} ({1:,} bp masked{2})".format(
            self.name, self.regions.length() if self.regions else 0,
            ", sha256 {0}".format(self.checksum[:12]) if self.checksum else "")

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "applied": self.applied,
            "reason": self.reason,
            "path": self.path,
            "checksum": self.checksum,
            "source": self.source,
            "masked_bases": self.masked_bases(),
            "intervals": (None if self.regions is None
                          else sum(len(v) for v in self.regions.intervals.values())),
        }


def load_mask(path: PathLike, name: str = "", source: str = "") -> Mask:
    """Load a mask BED, or raise saying exactly where to get one."""
    return Mask.from_bed(path, name=name, source=source)


# ---------------------------------------------------------------------------
# Pairwise counting
# ---------------------------------------------------------------------------

def _proximity_survivors(candidates: Sequence[Tuple[str, int]],
                         neighbours: Sequence[Tuple[str, int]],
                         window: int) -> Tuple[List[Tuple[str, int]], int]:
    """Drop candidate differences that have another variable position too close.

    The rule counted in the 508-isolate study is "no other SNP within 12 bases".
    Mjolnir reads *other SNP* as any comparable position at which either member
    of the pair departs from the reference — so the other sample's private
    variants count, and so do indels, since a dense patch of variation is the
    signature of misalignment whichever genome carries it and whatever shape it
    takes. Positions already removed by the mask or with no comparable allele
    are not neighbours: they have been excluded from the comparison and cannot
    then be used to exclude something else. That reading is stricter than
    counting only the differences themselves, and it is spelled out here because
    the published sentence does not settle it.
    """
    if window <= 0:
        return list(candidates), 0
    by_chrom: Dict[str, List[int]] = {}
    for chrom, pos in neighbours:
        by_chrom.setdefault(chrom, []).append(pos)
    for positions in by_chrom.values():
        positions.sort()

    kept: List[Tuple[str, int]] = []
    dropped = 0
    for chrom, pos in candidates:
        positions = by_chrom.get(chrom, [])
        # The nearest neighbouring position on either side decides it, so the
        # two bisects replace a scan of every variable position in the pair —
        # which is the difference between seconds and hours on a 159-sample set.
        left = bisect.bisect_left(positions, pos) - 1
        right = bisect.bisect_right(positions, pos)
        close = ((left >= 0 and pos - positions[left] <= window)
                 or (right < len(positions) and positions[right] - pos <= window))
        if close:
            dropped += 1
        else:
            kept.append((chrom, pos))
    return kept, dropped


def pairwise_distance(table: JointTable, sample_a: str, sample_b: str, mask: Mask,
                      proximity_window: int = SNP_PROXIMITY_WINDOW,
                      ) -> PairwiseDistance:
    """Masked SNP distance between two samples, with its shared denominator.

    The denominator is the callable region of *a* intersected with the callable
    region of *b*, minus the mask, minus the individual positions where either
    sample's allele is ambiguous — so it is the sequence over which a difference
    could actually have been seen, not the length of the reference.

    Returns ``snps=None`` when that cannot be established. A pair with no
    denominator has no distance; it does not have a distance of zero.
    """
    if sample_a not in table.samples or sample_b not in table.samples:
        missing = [s for s in (sample_a, sample_b) if s not in table.samples]
        raise MjolnirError(
            "sample(s) {0} are not in the joint table (it holds {1})".format(
                ", ".join(missing), ", ".join(table.samples)))
    if sample_a == sample_b:
        raise MjolnirError(
            "pairwise_distance was asked for the distance from {0!r} to itself; a "
            "self-comparison is not a measurement".format(sample_a))

    first, second = pair_key(sample_a, sample_b)
    regions_a = table.callable_regions.get(first)
    regions_b = table.callable_regions.get(second)
    if regions_a is None or regions_b is None:
        absent = [name for name, regions in ((first, regions_a), (second, regions_b))
                  if regions is None]
        return PairwiseDistance(
            sample_a=first, sample_b=second, snps=None, shared_callable_sites=None,
            masked_sites=None,
            note=NO_DENOMINATOR_TEXT.format(samples=" and ".join(absent)))

    shared = regions_a.intersect(regions_b, name="{0}&{1}".format(first, second))
    masked_in_shared = 0
    if mask.applied and mask.regions is not None:
        comparable = shared.subtract(mask.regions)
        masked_in_shared = shared.length() - comparable.length()
    else:
        comparable = shared

    candidates: List[Tuple[str, int]] = []
    neighbours: List[Tuple[str, int]] = []
    masked_variant_sites = 0
    ambiguous_sites = 0
    indel_sites = 0
    outside_shared = 0

    index = table.site_index()
    union = sorted(set(table.positions_for(first)) | set(table.positions_for(second)))
    for chrom, pos in union:
        site = index.get((chrom, pos))
        if site is None:  # pragma: no cover - positions come from the table itself
            continue
        if mask.contains(chrom, pos):
            masked_variant_sites += 1
            continue
        if not comparable.contains(chrom, pos):
            outside_shared += 1
            continue
        allele_a = table.allele(first, site)
        allele_b = table.allele(second, site)
        if allele_a is None or allele_b is None:
            # Inside the shared callable region but with no allele: possible when
            # callable regions and calls come from different stages. Counted, not
            # guessed at.
            ambiguous_sites += 1
            continue
        if allele_a == AMBIGUOUS or allele_b == AMBIGUOUS:
            ambiguous_sites += 1
            continue
        neighbours.append((chrom, pos))
        if allele_a == allele_b:
            continue
        if not site.is_snp:
            indel_sites += 1
            continue
        candidates.append((chrom, pos))

    kept, dropped_by_proximity = _proximity_survivors(candidates, neighbours,
                                                      proximity_window)
    shared_sites = comparable.length() - ambiguous_sites
    if shared_sites < 0:  # pragma: no cover - only reachable on inconsistent input
        shared_sites = 0

    detail = [
        "mask: {0}".format(mask.describe()),
        "{0} variant positions inside the mask".format(masked_variant_sites),
        "{0} positions with no comparable allele in one member".format(ambiguous_sites),
        "{0} indel positions not counted as SNPs".format(indel_sites),
        "{0} differences dropped by the {1} bp proximity rule".format(
            dropped_by_proximity, proximity_window),
    ]
    if outside_shared:
        detail.append("{0} variant positions outside the shared callable region"
                      .format(outside_shared))
    if not mask.applied:
        detail.append(MASK_ABSENT_TEXT)

    return PairwiseDistance(
        sample_a=first, sample_b=second, snps=len(kept),
        shared_callable_sites=shared_sites, masked_sites=masked_in_shared,
        note="; ".join(detail))


# ---------------------------------------------------------------------------
# The matrix
# ---------------------------------------------------------------------------

@dataclass
class DistanceMatrix:
    """Every pair's distance, each one carrying its own denominator.

    There is intentionally no accessor here that returns a bare SNP count. The
    unit of this API is :class:`~mjolnir.records.PairwiseDistance`, so a caller
    who wants the number has the shared-callable figure in the same object and
    the report cannot print one without having had the other in its hand.
    """

    samples: List[str] = field(default_factory=list)
    pairs: List[PairwiseDistance] = field(default_factory=list)
    mask: Optional[Mask] = None
    reference: str = ""
    proximity_window: int = SNP_PROXIMITY_WINDOW
    min_shared_callable_sites: int = MIN_SHARED_CALLABLE_SITES
    checks: List[Check] = field(default_factory=list)
    caveats: List[str] = field(default_factory=list)

    def pair(self, a: str, b: str) -> Optional[PairwiseDistance]:
        wanted = pair_key(a, b)
        for entry in self.pairs:
            if entry.key == wanted:
                return entry
        return None

    def comparable_pairs(self) -> List[PairwiseDistance]:
        """Pairs with both a distance and a denominator large enough to mean it."""
        return [p for p in self.pairs
                if p.snps is not None
                and (p.shared_callable_sites or 0) >= self.min_shared_callable_sites]

    def uncomputed_pairs(self) -> List[PairwiseDistance]:
        return [p for p in self.pairs if p.snps is None]

    def thin_pairs(self) -> List[PairwiseDistance]:
        """Pairs whose denominator is too small to compare to the SNP thresholds."""
        return [p for p in self.pairs
                if p.snps is not None
                and (p.shared_callable_sites or 0) < self.min_shared_callable_sites]

    def describe_pair(self, a: str, b: str) -> str:
        """The one line a report prints: the count and what it was counted over."""
        found = self.pair(a, b)
        if found is None:
            return "{0} vs {1}: not compared".format(a, b)
        return "{0} vs {1}: {2}, mask {3}".format(
            a, b, format_distance(found),
            self.mask.describe() if self.mask else "unstated")

    def to_dict(self) -> Dict[str, Any]:
        return {
            "samples": list(self.samples),
            "reference": self.reference,
            "mask": self.mask.to_dict() if self.mask else None,
            "proximity_window": self.proximity_window,
            "min_shared_callable_sites": self.min_shared_callable_sites,
            "pairs": [
                {
                    "sample_a": p.sample_a,
                    "sample_b": p.sample_b,
                    "snps": p.snps,
                    "shared_callable_sites": p.shared_callable_sites,
                    "masked_sites": p.masked_sites,
                    "snps_per_mb": p.snps_per_mb,
                    "note": p.note,
                }
                for p in self.pairs
            ],
            "n_comparable": len(self.comparable_pairs()),
            "n_uncomputed": len(self.uncomputed_pairs()),
            "n_thin": len(self.thin_pairs()),
            "checks": [c.to_dict() for c in self.checks],
            "caveats": list(self.caveats),
        }


def distance_matrix(table: JointTable, mask: Mask,
                    proximity_window: int = SNP_PROXIMITY_WINDOW,
                    min_shared_callable_sites: int = MIN_SHARED_CALLABLE_SITES,
                    samples: Optional[Sequence[str]] = None) -> DistanceMatrix:
    """All pairwise distances over a joint table.

    *mask* is required. Passing ``None`` raises rather than quietly counting
    unmasked, and the message names the two ways forward — a mask file, or
    :meth:`Mask.absent` with a reason that the report will print.
    """
    if mask is None:
        raise MjolnirError(
            "a mask is required before SNP distances are counted: the 508-isolate "
            "MTBC comparison masked {0:,} loci (~{1:.0%} of H37Rv) of repetitive, "
            "low-complexity and error-prone sequence first. Supply one with "
            "load_mask('mask.bed') — 'mjolnir db fetch tbdb' provides the MTBC mask "
            "— or, for a reference with no mask yet, state the decision explicitly "
            "with Mask.absent(reason=...).".format(MASKED_LOCI_H37RV,
                                                   MASKED_FRACTION_H37RV))

    names = list(samples) if samples is not None else list(table.samples)
    unknown = [name for name in names if name not in table.samples]
    if unknown:
        raise MjolnirError(
            "sample(s) {0} are not in the joint table".format(", ".join(unknown)))

    _require_matching_contigs(table, mask)

    matrix = DistanceMatrix(samples=names, mask=mask, reference=table.reference,
                            proximity_window=proximity_window,
                            min_shared_callable_sites=min_shared_callable_sites)

    for i, first in enumerate(names):
        for second in names[i + 1:]:
            matrix.pairs.append(
                pairwise_distance(table, first, second, mask,
                                  proximity_window=proximity_window))

    matrix.checks.extend(_matrix_checks(matrix))
    if not mask.applied:
        matrix.caveats.append(MASK_ABSENT_TEXT)
    matrix.caveats.extend(table.caveats)
    LOG.debug("distance matrix: %d pairs, %d comparable",
              len(matrix.pairs), len(matrix.comparable_pairs()))
    return matrix


def _require_matching_contigs(table: JointTable, mask: Mask) -> None:
    """Refuse a mask whose contig names match nothing in the cohort.

    tbdb's BED files name the H37Rv contig ``Chromosome`` while the FASTA most
    pipelines map against calls it ``NC_000962.3``. A mask under the wrong name
    masks nothing at all, silently, and every distance comes back inflated by
    the repetitive regions it was supposed to remove — a wrong number that looks
    exactly like a right one. So the mismatch is fatal here rather than
    discovered later in a cluster that should not exist.
    """
    if not mask.applied or mask.regions is None:
        return
    cohort_contigs = set(site.chrom for site in table.sites)
    for regions in table.callable_regions.values():
        cohort_contigs.update(regions.chroms())
    if not cohort_contigs:
        return
    mask_contigs = set(mask.regions.chroms())
    if cohort_contigs & mask_contigs:
        return
    raise MjolnirError(
        "the mask {0!r} names contig(s) {1} but this cohort was called against "
        "{2}; under those names the mask would exclude nothing and every distance "
        "would silently include the repetitive regions it exists to remove. "
        "Rename the mask's contigs to match the reference, or supply the mask "
        "built for this reference.".format(
            mask.name, ", ".join(sorted(mask_contigs)[:5]),
            ", ".join(sorted(cohort_contigs)[:5])))


def _matrix_checks(matrix: DistanceMatrix) -> List[Check]:
    """The checks a cohort report prints above the matrix."""
    checks: List[Check] = []

    mask = matrix.mask
    if mask is not None and mask.applied:
        checks.append(Check(
            name="mask_applied", value=mask.name, threshold="a named mask",
            comparison="==", status=STATUS_PASS, category="cohort",
            source=source_for("masked_loci_h37rv"),
            reading="distances were counted after masking {0}. The mask is a named, "
                    "versioned input and not a constant compiled into the tool: tbdb's "
                    "changed twice in three years.".format(mask.describe())))
    else:
        checks.append(Check(
            name="mask_applied", value=False, threshold=True, comparison="==",
            status=STATUS_FAIL, category="cohort",
            source=source_for("masked_loci_h37rv"),
            reading=MASK_ABSENT_TEXT + " Reason recorded: {0}".format(
                mask.reason if mask else "no mask object supplied")))

    checks.append(Check(
        name="snp_proximity_window", value=matrix.proximity_window,
        threshold=SNP_PROXIMITY_WINDOW, comparison="==", unit="bp",
        status=STATUS_PASS if matrix.proximity_window == SNP_PROXIMITY_WINDOW
        else STATUS_WARN,
        category="cohort", source=source_for("snp_proximity_window"),
        reading="a difference was counted only where no other variable position in "
                "either genome lay within {0} bp of it.".format(matrix.proximity_window)))

    if not matrix.pairs:
        checks.append(Check.not_measured(
            "pairwise_distances",
            "no pairs were compared: a cohort of {0} sample(s) has no pair. This is "
            "an absence of comparison, not a finding of relatedness or of "
            "difference.".format(len(matrix.samples)),
            category="cohort"))
        return checks

    denominators = [p.shared_callable_sites for p in matrix.pairs
                    if p.shared_callable_sites is not None]
    smallest = min(denominators) if denominators else None
    checks.append(Check.numeric(
        "min_shared_callable_sites", float(smallest) if smallest is not None else None,
        warn_minimum=float(matrix.min_shared_callable_sites), unit="bp", category="cohort",
        source=source_for("min_shared_callable_sites"),
        reading="the thinnest pair in this cohort was compared over {0} bp of shared "
                "callable sequence; below {1:,} bp a SNP count is not comparable to "
                "the published 5-SNP and 12-SNP thresholds.".format(
                    "{0:,}".format(smallest) if smallest is not None else "an unknown number of",
                    matrix.min_shared_callable_sites),
        not_measured_why="no pair had a shared callable denominator, so no distance "
                         "in this cohort has been computed"))

    uncomputed = matrix.uncomputed_pairs()
    checks.append(Check(
        name="pairs_computed", value=len(matrix.pairs) - len(uncomputed),
        threshold=len(matrix.pairs), comparison="==", unit="pairs", category="cohort",
        status=STATUS_PASS if not uncomputed else STATUS_WARN,
        reading="{0} of {1} pairs were compared; {2} had no shared callable "
                "denominator and are reported as not computed rather than as "
                "distant or identical.".format(
                    len(matrix.pairs) - len(uncomputed), len(matrix.pairs),
                    len(uncomputed))))
    return checks


def pairs_for_cohort(matrix: DistanceMatrix) -> List[PairwiseDistance]:
    """The pair records as :class:`~mjolnir.records.CohortResult` wants them."""
    return list(matrix.pairs)


def masked_fraction(mask: Mask, reference_length: Optional[int]) -> Optional[float]:
    """How much of the reference this mask removes, or None when unknowable."""
    return mask.fraction_of(reference_length)


def format_distance(entry: Optional[PairwiseDistance]) -> str:
    """One human line for a pair, denominator included, never a naked number."""
    if entry is None:
        return "not compared"
    if entry.snps is None:
        return "not computed ({0})".format(entry.note or "no denominator")
    return "{0} SNPs / {1:,} bp shared callable".format(
        entry.snps, entry.shared_callable_sites or 0)


def iter_comparable(pairs: Iterable[PairwiseDistance],
                    min_shared_callable_sites: int = MIN_SHARED_CALLABLE_SITES,
                    ) -> List[PairwiseDistance]:
    """Pairs a cluster may be built from: computed, and over enough sequence."""
    return [p for p in pairs
            if p.snps is not None
            and (p.shared_callable_sites or 0) >= min_shared_callable_sites]
