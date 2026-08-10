"""Heterozygosity: the contamination signal that actually works on mycobacteria.

Everything in this file measures the same physical thing — reads at a position
that disagree with each other — and the design's §8 list of what may honestly be
measured is almost entirely made of it. That is not an accident. A taxonomic
read classifier cannot separate MTBC members and loses ~93% of true target reads
on a standard index, a marker-gene completeness statistic cannot see a
same-species mixture at all, and neither can say anything about the sample this
report is about. Allele fractions can.

Three quantities, and they are deliberately not one number.

**F2 and F47** are the minor-allele fraction summarised over lineage-defining
positions: the operationally validated route to separating a genuine mixed
infection from cross-contamination, because a mixture of two lineages shows up
at exactly the positions that define them, and a batch of cross-contaminated
samples shows the same minor lineage across the batch. What the two names mean
is stated in :func:`f_statistic` — F2 is the mean of the two highest minor-allele
fractions and F47 of the 47 highest — and the honest position on that definition
is recorded there too, because the primary citation was not available on this
machine.

**The genome-wide heterozygous-SNP fraction** under the MixInfect filters
(Q >= 20, DP >= 10) is the broader signal, and it is reported as a two-tier
classification rather than a single cutoff. MixInfect itself fits a model; the
two tiers are Mjolnir's, are labelled as policy in ``config.py``, and exist
because the quantity does not separate cleanly enough for one line to be honest.

**The unambiguous-base fraction** is MTBseq's de-facto heterozygosity filter and
is surfaced here rather than applied silently. MTBseq's 75% majority rule
decides a position and throws the minority reads away; the reads it threw away
are the entire contamination signal, so Mjolnir computes the same number and
prints it instead of consuming it.

Two structural refusals live in this module rather than in its documentation.
On ONT a *low* heterozygosity cannot establish a single-strain sample — 26 of 27
Illumina-only minor SNPs in the 508-isolate comparison were visible in the ONT
pileup and not called — so :func:`classify_mixture` will not return
"single-strain" on that platform; it returns "not-assessed" and says why. On
FASTA there are no allele fractions at all, so every quantity here is ``None``
and the class is "not-assessed". Neither case is allowed to read as clean.

Allele fractions here use a conventional ACGT denominator. MTBseq's includes its
N and GAP counts and lets GAP win every tie, so the same pileup yields different
fractions under the two tools (design §9b); ``--mtbseq-compat`` is where that
convention is reproduced, and it is not the default.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ..config import (
    Config,
    F2_MIXTURE_THRESHOLD,
    F47_MIXTURE_THRESHOLD,
    FASTA_CAPABILITY_LOSS,
    HET_MIN_DEPTH,
    HET_MIN_QUAL,
    HET_SNP_FRACTION_MIXED,
    HET_SNP_FRACTION_WARN,
    MAJOR_VARIANT_FRACTION,
    MIN_BARCODE_CALLABLE_FRACTION,
    MIN_MINOR_VARIANT_FRACTION,
    MIN_UNAMBIGUOUS_FRACTION,
    MTBSEQ_MINFREQ,
    MTBSEQ_UNAMBIG,
    ONT_MINOR_VARIANT_CAVEAT,
    source_for,
    threshold,
)
from ..records import (
    MIXTURE_MIXED,
    MIXTURE_NOT_ASSESSED,
    MIXTURE_POSSIBLE,
    MIXTURE_SINGLE,
    PLATFORM_FASTA,
    PLATFORM_ONT,
    STATUS_FAIL,
    STATUS_PASS,
    STATUS_WARN,
    Check,
    normalise_platform,
)
from ..utils import MjolnirError, mean, safe_fraction

#: The four bases a fraction may be computed over. Anything else a pileup
#: reports — N, deletions, reference skips — is counted in ``raw_depth`` and
#: excluded from the denominator, which is the conventional ACGT depth and is
#: the documented divergence from MTBseq (design §9b).
BASES: Tuple[str, ...] = ("A", "C", "G", "T")

#: Report-facing sentence for the filters, so the number and the filters that
#: produced it are never separated.
MIXINFECT_FILTER_TEXT = (
    "heterozygous sites counted under the MixInfect filters: base quality at or "
    "above {0} and depth at or above {1}x".format(HET_MIN_QUAL, HET_MIN_DEPTH)
)

#: Sentence printed with the unambiguous-base fraction. MTBseq applies this
#: filter and discards what fails it; Mjolnir reports what would have been
#: discarded, because that is the contamination signal.
UNAMBIGUOUS_SURFACED_TEXT = (
    "the unambiguous-base fraction is MTBseq's de-facto heterozygosity filter "
    "({0}% majority at a position, {1}% of positions). MTBseq applies it and "
    "discards the minority reads; Mjolnir reports it, because those reads are "
    "the mixture signal".format(MTBSEQ_MINFREQ, MTBSEQ_UNAMBIG)
)

#: Sentence printed when a single-strain conclusion is withheld on ONT.
ONT_SINGLE_STRAIN_REFUSAL = (
    "a low heterozygosity on ONT cannot establish a single-strain sample, so no "
    "single-strain conclusion is drawn: " + ONT_MINOR_VARIANT_CAVEAT
)

#: How alarming each mixture class is, for folding several signals into one.
#: ``not-assessed`` is deliberately absent: it is not a point on this scale, it
#: is the statement that the scale was never reached.
_MIXTURE_SEVERITY: Dict[str, int] = {
    MIXTURE_SINGLE: 0,
    MIXTURE_POSSIBLE: 1,
    MIXTURE_MIXED: 2,
}

#: The check status each mixture class produces. A possible mixture warns rather
#: than fails because the tier it sits in is Mjolnir policy, not a published
#: cutoff; a called mixture fails because a second strain invalidates the
#: quantities the rest of the report is built from.
_MIXTURE_STATUS: Dict[str, str] = {
    MIXTURE_SINGLE: STATUS_PASS,
    MIXTURE_POSSIBLE: STATUS_WARN,
    MIXTURE_MIXED: STATUS_FAIL,
}

#: Panel sizes for the two F statistics. They are the numbers in the names, not
#: thresholds applied to data, so they are stated here beside the definition
#: rather than registered in ``config.py``.
F2_TOP_SITES = 2
F47_TOP_SITES = 47

F_STATISTIC_DEFINITION = (
    "F{n} is the mean of the {n} highest minor-allele fractions across callable "
    "lineage-defining positions"
)


# ---------------------------------------------------------------------------
# What a caller hands in
# ---------------------------------------------------------------------------

@dataclass
class SiteObservation:
    """One pileup position: how many reads carried each base.

    A pileup row rather than a variant call, because the whole point of this
    module is the signal a majority-rule caller discards. ``engines/pileup.py``
    produces these; nothing here runs a tool.

    ``raw_depth`` is whatever the pileup reported as total depth, including Ns
    and gaps. It is kept because MTBseq's frequency denominator includes them
    and a lab reconciling the two tools needs both numbers, but every fraction
    computed here uses :attr:`depth`, the ACGT total.
    """

    pos: int
    counts: Dict[str, int] = field(default_factory=dict)
    qual: Optional[float] = None
    chrom: str = ""
    ref_base: str = ""
    raw_depth: Optional[int] = None

    def __post_init__(self) -> None:
        self.pos = int(self.pos)
        cleaned: Dict[str, int] = {}
        for base, count in dict(self.counts).items():
            key = str(base).upper()
            if count is None:
                continue
            cleaned[key] = cleaned.get(key, 0) + int(count)
        self.counts = cleaned

    @property
    def depth(self) -> int:
        """ACGT depth. The conventional denominator, not MTBseq's."""
        return sum(self.counts.get(base, 0) for base in BASES)

    @property
    def ranked(self) -> List[Tuple[str, int]]:
        """(base, count) over ACGT, most-supported first, ties broken by base.

        The tie-break is alphabetical and therefore arbitrary; it is stated
        because MTBseq's is not arbitrary — it orders A < C < G < T < N < GAP and
        lets GAP win — and a position split exactly evenly is resolved
        differently by the two tools.
        """
        return sorted(((b, self.counts.get(b, 0)) for b in BASES),
                      key=lambda item: (-item[1], item[0]))

    @property
    def major_base(self) -> Optional[str]:
        top = self.ranked[0]
        return top[0] if top[1] > 0 else None

    @property
    def major_fraction(self) -> Optional[float]:
        return safe_fraction(self.ranked[0][1], self.depth)

    @property
    def minor_base(self) -> Optional[str]:
        second = self.ranked[1]
        return second[0] if second[1] > 0 else None

    @property
    def minor_fraction(self) -> Optional[float]:
        """Second-most-supported ACGT allele over ACGT depth, or None at 0x.

        None rather than 0.0 at zero depth, all the way down: a position with no
        reads is not a homozygous position, and the difference has to survive as
        far as the report.
        """
        return safe_fraction(self.ranked[1][1], self.depth)

    def fraction_of(self, base: str) -> Optional[float]:
        return safe_fraction(self.counts.get(str(base).upper(), 0), self.depth)

    def passes_filters(self, min_depth: int = HET_MIN_DEPTH,
                       min_qual: int = HET_MIN_QUAL) -> bool:
        """MixInfect's filters: depth at or above 10x, quality at or above Q20.

        A site with no quality recorded fails on quality alone, rather than
        being waved through: an unfiltered site counted as heterozygous is a
        false mixture signal, and this module's output gates a sample's validity.
        """
        if self.depth < int(min_depth):
            return False
        if self.qual is None:
            return False
        return float(self.qual) >= float(min_qual)


@dataclass
class LineageSite(SiteObservation):
    """A lineage-defining position, with the allele that defines the lineage.

    ``derived_allele`` is what makes the F statistics interpretable rather than
    merely numeric: with it, a lineage present as a minority is the lineage whose
    defining allele sits between the minor-variant floor and the major-variant
    threshold, which is a statement about strains. Without it all that can be
    said is that the position is heterozygous, and :class:`LineageMixtureStats`
    records which of the two it had.
    """

    lineage: str = ""
    derived_allele: str = ""
    ancestral_allele: str = ""

    def __post_init__(self) -> None:
        SiteObservation.__post_init__(self)
        self.derived_allele = str(self.derived_allele or "").upper()
        self.ancestral_allele = str(self.ancestral_allele or "").upper()

    @property
    def derived_fraction(self) -> Optional[float]:
        """Read support for the lineage-defining allele, or None if unknown."""
        if not self.derived_allele:
            return None
        return self.fraction_of(self.derived_allele)


# ---------------------------------------------------------------------------
# What this module hands back
# ---------------------------------------------------------------------------

@dataclass
class LineageMixtureStats:
    """F2, F47 and which lineages carry minority evidence."""

    f2: Optional[float] = None
    f47: Optional[float] = None
    sites_total: int = 0
    sites_examined: int = 0
    het_sites: int = 0
    mixed_lineages: List[str] = field(default_factory=list)
    per_lineage: List[Dict[str, Any]] = field(default_factory=list)
    #: "derived-allele" when defining alleles were supplied, "minor-allele"
    #: when only the pileup was, since the two support different statements.
    basis: str = ""
    note: str = ""

    @property
    def het_fraction(self) -> Optional[float]:
        return safe_fraction(self.het_sites, self.sites_examined)


@dataclass
class HeterozygosityResult:
    """Every heterozygosity quantity, its class, and the checks it produced."""

    f2: Optional[float] = None
    f47: Optional[float] = None
    lineage_het_sites: Optional[int] = None
    lineage_sites_examined: Optional[int] = None
    het_snp_fraction: Optional[float] = None
    het_snp_count: Optional[int] = None
    snp_sites_examined: Optional[int] = None
    unambiguous_fraction: Optional[float] = None
    mixture_class: str = MIXTURE_NOT_ASSESSED
    mixture_reason: str = ""
    mixed_lineages: List[str] = field(default_factory=list)
    method: str = ""
    caveats: List[str] = field(default_factory=list)
    checks: List[Check] = field(default_factory=list)

    @property
    def assessed(self) -> bool:
        """Whether any heterozygosity quantity was actually computed.

        The report asks this before it says anything about mixture at all, so
        that "not-assessed" cannot be rendered as a quiet pass.
        """
        return any(value is not None for value in
                   (self.f2, self.f47, self.het_snp_fraction, self.unambiguous_fraction))


# ---------------------------------------------------------------------------
# Measurements
# ---------------------------------------------------------------------------

def callable_sites(sites: Sequence[SiteObservation],
                   min_depth: int = HET_MIN_DEPTH,
                   min_qual: int = HET_MIN_QUAL) -> List[SiteObservation]:
    """The subset of *sites* that passes the MixInfect filters."""
    return [site for site in sites if site.passes_filters(min_depth, min_qual)]


def is_heterozygous(site: SiteObservation,
                    min_minor_fraction: float = MIN_MINOR_VARIANT_FRACTION) -> Optional[bool]:
    """Whether a site carries a real second allele, or None if it was not callable.

    None rather than False for an uncallable site. A position nobody could see
    is not a position where nothing was found, and counting it as homozygous
    would dilute the genome-wide fraction towards zero exactly when coverage is
    worst — which is when a mixture is hardest to see and most worth reporting.
    """
    minor = site.minor_fraction
    if minor is None:
        return None
    return minor >= float(min_minor_fraction)


def unambiguous_fraction(sites: Sequence[SiteObservation],
                         majority: float = MTBSEQ_UNAMBIG / 100.0) -> Optional[float]:
    """Fraction of callable positions with an unambiguous majority allele.

    MTBseq computes this to decide which positions to keep. Mjolnir computes it
    to report what MTBseq would have dropped: at 95% of positions unambiguous,
    one position in twenty carries reads that disagree, and whether that is
    library noise or a second organism is the question the rest of this module
    exists to answer.
    """
    considered = [site for site in sites if site.depth > 0]
    if not considered:
        return None
    clean = 0
    for site in considered:
        top = site.major_fraction
        if top is not None and top >= float(majority):
            clean += 1
    return safe_fraction(clean, len(considered))


def genome_wide_heterozygosity(
        snp_sites: Sequence[SiteObservation],
        *,
        min_depth: int = HET_MIN_DEPTH,
        min_qual: int = HET_MIN_QUAL,
        min_minor_fraction: float = MIN_MINOR_VARIANT_FRACTION,
) -> Tuple[Optional[float], int, int]:
    """(fraction, heterozygous sites, sites examined) over the SNP sites given.

    The denominator is the number of *callable SNP sites*, matching MixInfect,
    whose statistic is heterozygous SNPs over called SNPs — which is also why
    :class:`~mjolnir.records.ContaminationResult` names the field
    ``snp_sites_examined``. Handing this function every callable position in the
    genome instead would produce a much smaller number that looks like the same
    statistic and is not, so the caller decides the denominator and the record
    prints it.
    """
    usable = callable_sites(snp_sites, min_depth, min_qual)
    if not usable:
        return None, 0, 0
    het = 0
    for site in usable:
        if is_heterozygous(site, min_minor_fraction):
            het += 1
    return safe_fraction(het, len(usable)), het, len(usable)


def f_statistic(fractions: Sequence[Optional[float]], top_n: int) -> Optional[float]:
    """Mean of the *top_n* highest minor-allele fractions, or None.

    None — never a mean over fewer sites — when fewer than *top_n* callable
    positions exist. The mean of the top 47 of 47 is a different statistic from
    the mean of the top 47 of 1,100, and the mean of "the top 47" computed over
    five sites is not that statistic at all; returning it would make a
    thinly-covered sample look decisively clean or decisively mixed on almost no
    evidence.

    On the definition itself: F2 and F47 are named for their panel sizes, and
    the primary paper defining them was not available on this machine (the
    thresholds in ``config.py`` are marked unverified for the same reason). This
    is Mjolnir's operationalisation of the names, it travels with the site count
    it was computed over, and it is stated in the report as a definition rather
    than assumed to be the reader's.
    """
    ordered = sorted(f for f in fractions if f is not None)
    if len(ordered) < int(top_n):
        return None
    return mean(ordered[-int(top_n):])


def lineage_mixture_statistics(
        sites: Sequence[LineageSite],
        *,
        min_depth: int = HET_MIN_DEPTH,
        min_qual: int = HET_MIN_QUAL,
        min_minor_fraction: float = MIN_MINOR_VARIANT_FRACTION,
        major_fraction: float = MAJOR_VARIANT_FRACTION,
        min_callable_fraction: float = MIN_BARCODE_CALLABLE_FRACTION,
) -> LineageMixtureStats:
    """F2/F47 and per-lineage minority evidence over a lineage-defining panel.

    A lineage is listed in ``mixed_lineages`` when enough of its defining sites
    were callable and its defining allele is present as a minority — above the
    minor-variant floor and below the major-variant threshold. The sample's own
    lineage generally appears in that list alongside the contaminant's, because
    at its defining sites the second strain contributes the ancestral allele as
    the minority. The list is therefore evidence that more than one lineage is
    present, not a resolved pair of strains, and the report must not read it as
    a genotype.
    """
    stats = LineageMixtureStats(sites_total=len(sites))
    usable = [s for s in sites if s.passes_filters(min_depth, min_qual)]
    stats.sites_examined = len(usable)
    if not usable:
        stats.note = (
            "no lineage-defining position passed the MixInfect filters, so F2 "
            "and F47 were not computed"
        )
        return stats

    fractions = [site.minor_fraction for site in usable]
    stats.f2 = f_statistic(fractions, F2_TOP_SITES)
    stats.f47 = f_statistic(fractions, F47_TOP_SITES)
    stats.het_sites = sum(1 for site in usable
                          if is_heterozygous(site, min_minor_fraction))

    have_derived = any(site.derived_allele for site in usable)
    stats.basis = "derived-allele" if have_derived else "minor-allele"
    if not have_derived:
        stats.note = (
            "no lineage-defining alleles were supplied, so minority evidence is "
            "reported per position rather than per lineage"
        )

    by_lineage: Dict[str, List[LineageSite]] = {}
    for site in sites:
        if site.lineage:
            by_lineage.setdefault(site.lineage, []).append(site)

    for lineage in sorted(by_lineage):
        panel = by_lineage[lineage]
        panel_usable = [s for s in panel if s.passes_filters(min_depth, min_qual)]
        callable_share = safe_fraction(len(panel_usable), len(panel))
        if have_derived:
            values = [s.derived_fraction for s in panel_usable if s.derived_fraction is not None]
        else:
            values = [s.minor_fraction for s in panel_usable if s.minor_fraction is not None]
        average = mean(values)
        entry: Dict[str, Any] = {
            "lineage": lineage,
            "sites_total": len(panel),
            "sites_callable": len(panel_usable),
            "callable_fraction": callable_share,
            "mean_fraction": average,
            "basis": stats.basis,
        }
        minority = (
            average is not None
            and callable_share is not None
            and callable_share >= float(min_callable_fraction)
            and float(min_minor_fraction) <= average < float(major_fraction)
        )
        entry["minority_evidence"] = bool(minority)
        stats.per_lineage.append(entry)
        if minority:
            stats.mixed_lineages.append(lineage)
    return stats


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify_mixture(het_snp_fraction: Optional[float],
                     f2: Optional[float] = None,
                     f47: Optional[float] = None,
                     *,
                     platform: str,
                     config: Optional[Config] = None) -> Tuple[str, str]:
    """Two tiers, not one cutoff, and no single-strain conclusion on ONT or FASTA.

    The tiers come from ``config``: below the warn tier the sample reads as a
    single strain, between the tiers a mixture is possible, above the mixed tier
    it is called. Both are Mjolnir policy rather than published cutoffs — the
    underlying quantity does not separate cleanly enough for one line — and the
    report prints them as policy.

    The two platform refusals are the reason this returns a reason string as
    well as a class. On FASTA there are no allele fractions to classify. On ONT a
    low fraction is uninformative in the direction that matters, so the class
    becomes "not-assessed" rather than "single-strain": absence of a minor
    variant on ONT is not absence of a subpopulation, and the one thing this
    module must never emit is a quiet all-clear.
    """
    plat = normalise_platform(platform)
    if plat == PLATFORM_FASTA:
        return MIXTURE_NOT_ASSESSED, FASTA_CAPABILITY_LOSS

    warn_tier = HET_SNP_FRACTION_WARN
    mixed_tier = HET_SNP_FRACTION_MIXED
    minor_floor = config.min_minor_variant_fraction if config else MIN_MINOR_VARIANT_FRACTION

    classes: List[str] = []
    reasons: List[str] = []

    if het_snp_fraction is not None:
        if het_snp_fraction >= mixed_tier:
            classes.append(MIXTURE_MIXED)
            reasons.append(
                "genome-wide heterozygous-SNP fraction {0:.4f} is at or above the "
                "mixed tier {1}".format(het_snp_fraction, mixed_tier))
        elif het_snp_fraction >= warn_tier:
            classes.append(MIXTURE_POSSIBLE)
            reasons.append(
                "genome-wide heterozygous-SNP fraction {0:.4f} is between the "
                "possible-mixture tier {1} and the mixed tier {2}".format(
                    het_snp_fraction, warn_tier, mixed_tier))
        else:
            classes.append(MIXTURE_SINGLE)
            reasons.append(
                "genome-wide heterozygous-SNP fraction {0:.4f} is below the "
                "possible-mixture tier {1}, over sites filtered at Q{2}/{3}x".format(
                    het_snp_fraction, warn_tier, HET_MIN_QUAL, HET_MIN_DEPTH))

    for name, value, limit in (("F2", f2, F2_MIXTURE_THRESHOLD),
                               ("F47", f47, F47_MIXTURE_THRESHOLD)):
        if value is None:
            continue
        if value >= limit:
            classes.append(MIXTURE_POSSIBLE)
            reasons.append(
                "{0} = {1:.4f} at lineage-defining positions is at or above {2}"
                .format(name, value, limit))
        else:
            classes.append(MIXTURE_SINGLE)
            reasons.append(
                "{0} = {1:.4f} at lineage-defining positions is below {2}, and "
                "minor alleles below {3} are not reported".format(
                    name, value, limit, minor_floor))

    if not classes:
        return MIXTURE_NOT_ASSESSED, (
            "no heterozygosity metric could be computed for this sample, so no "
            "statement is made about mixture"
        )

    worst = MIXTURE_SINGLE
    for entry in classes:
        if _MIXTURE_SEVERITY[entry] > _MIXTURE_SEVERITY[worst]:
            worst = entry

    if worst == MIXTURE_SINGLE and plat == PLATFORM_ONT:
        return MIXTURE_NOT_ASSESSED, "; ".join(reasons + [ONT_SINGLE_STRAIN_REFUSAL])
    return worst, "; ".join(reasons)


# ---------------------------------------------------------------------------
# The panel
# ---------------------------------------------------------------------------

def _unverified_mark(name: str) -> str:
    """" [citation unverified]" when the threshold's source has not been checked."""
    return "" if threshold(name).verified else " [citation unverified]"


def heterozygosity_checks(result: "HeterozygosityResult",
                          platform: str) -> List[Check]:
    """Rule-derived verdicts for the contamination panel.

    Every quantity that could not be computed becomes a ``not_measured`` check
    with the reason attached, never an absent row: a panel that silently omits
    the metric it could not compute reads as a panel that passed.
    """
    plat = normalise_platform(platform)
    unavailable = ""
    if plat == PLATFORM_FASTA:
        unavailable = FASTA_CAPABILITY_LOSS

    checks: List[Check] = [
        Check.numeric(
            "heterozygous_snp_fraction", result.het_snp_fraction,
            warn_maximum=HET_SNP_FRACTION_WARN,
            maximum=HET_SNP_FRACTION_MIXED,
            source=source_for("het_snp_fraction_warn"),
            unit="fraction of callable SNP sites", category="contamination",
            reading=MIXINFECT_FILTER_TEXT,
            not_measured_why=unavailable or (
                "no SNP site passed the MixInfect filters (Q>={0}, DP>={1}x), so "
                "the genome-wide heterozygous-SNP fraction was not computed"
                .format(HET_MIN_QUAL, HET_MIN_DEPTH)),
        ),
        Check.numeric(
            "f2_lineage_heterozygosity", result.f2,
            warn_maximum=F2_MIXTURE_THRESHOLD,
            source=source_for("f2_mixture_threshold"),
            unit="minor-allele fraction", category="contamination",
            reading=F_STATISTIC_DEFINITION.format(n=F2_TOP_SITES)
                    + _unverified_mark("f2_mixture_threshold"),
            not_measured_why=unavailable or (
                "fewer than {0} lineage-defining positions were callable, so F2 "
                "was not computed".format(F2_TOP_SITES)),
        ),
        Check.numeric(
            "f47_lineage_heterozygosity", result.f47,
            warn_maximum=F47_MIXTURE_THRESHOLD,
            source=source_for("f47_mixture_threshold"),
            unit="minor-allele fraction", category="contamination",
            reading=F_STATISTIC_DEFINITION.format(n=F47_TOP_SITES)
                    + _unverified_mark("f47_mixture_threshold"),
            not_measured_why=unavailable or (
                "fewer than {0} lineage-defining positions were callable, so F47 "
                "was not computed".format(F47_TOP_SITES)),
        ),
        Check.numeric(
            "unambiguous_base_fraction", result.unambiguous_fraction,
            warn_minimum=MIN_UNAMBIGUOUS_FRACTION,
            source=source_for("min_unambiguous_fraction"),
            unit="fraction of positions", category="contamination",
            reading=UNAMBIGUOUS_SURFACED_TEXT,
            not_measured_why=unavailable or (
                "no position carried reads, so the unambiguous-base fraction was "
                "not computed"),
        ),
    ]

    if result.mixture_class == MIXTURE_NOT_ASSESSED:
        checks.append(Check.not_measured(
            "mixture_classification",
            result.mixture_reason or "mixture was not assessed",
            source=source_for("het_snp_fraction_warn"), category="contamination"))
    else:
        checks.append(Check(
            name="mixture_classification",
            value=result.mixture_class,
            threshold=MIXTURE_SINGLE,
            source=source_for("het_snp_fraction_warn"),
            status=_MIXTURE_STATUS[result.mixture_class],
            reading=result.mixture_reason,
            comparison="==", category="contamination"))

    checks.append(Check.boolean(
        "lineage_minority_evidence",
        bool(result.mixed_lineages) if result.lineage_sites_examined else None,
        expected=False,
        source=source_for("min_minor_variant_fraction"),
        category="contamination",
        reading=("minority evidence at the defining positions of: {0}".format(
            ", ".join(result.mixed_lineages)) if result.mixed_lineages else
            "no lineage carried minority allele evidence above {0}".format(
                MIN_MINOR_VARIANT_FRACTION)),
        fail_status=STATUS_WARN))
    return checks


def assess_heterozygosity(
        *,
        platform: str,
        snp_sites: Optional[Sequence[SiteObservation]] = None,
        lineage_sites: Optional[Sequence[LineageSite]] = None,
        unambiguous_sites: Optional[Sequence[SiteObservation]] = None,
        unambiguous_fraction_value: Optional[float] = None,
        config: Optional[Config] = None,
) -> HeterozygosityResult:
    """Measure every heterozygosity quantity this platform supports.

    ``unambiguous_fraction_value`` exists so a caller that already has MTBseq's
    number — from a compat run, or from ``TBstats`` output being reconciled —
    can pass it instead of the positions. Passing both is an error rather than a
    silent preference, because the two would be computed over different
    denominators and the report would print one while naming the other.

    FASTA short-circuits to "not assessed": an assembly has no allele fractions,
    and computing something else and calling it heterozygosity would turn a
    capability loss into a clean-looking result.
    """
    plat = normalise_platform(platform)
    result = HeterozygosityResult(method="allele fractions at pileup positions; "
                                         + MIXINFECT_FILTER_TEXT)

    if unambiguous_sites is not None and unambiguous_fraction_value is not None:
        raise MjolnirError(
            "assess_heterozygosity was given both unambiguous_sites and "
            "unambiguous_fraction_value; they have different denominators, so "
            "pass exactly one"
        )

    if plat == PLATFORM_FASTA:
        result.mixture_class = MIXTURE_NOT_ASSESSED
        result.mixture_reason = FASTA_CAPABILITY_LOSS
        result.caveats.append(FASTA_CAPABILITY_LOSS)
        result.checks = heterozygosity_checks(result, plat)
        return result

    min_depth = HET_MIN_DEPTH
    min_qual = HET_MIN_QUAL
    minor_floor = config.min_minor_variant_fraction if config else MIN_MINOR_VARIANT_FRACTION
    major = config.major_variant_fraction if config else MAJOR_VARIANT_FRACTION

    if snp_sites:
        fraction, het, examined = genome_wide_heterozygosity(
            snp_sites, min_depth=min_depth, min_qual=min_qual,
            min_minor_fraction=minor_floor)
        result.het_snp_fraction = fraction
        result.het_snp_count = het
        result.snp_sites_examined = examined

    if lineage_sites:
        stats = lineage_mixture_statistics(
            lineage_sites, min_depth=min_depth, min_qual=min_qual,
            min_minor_fraction=minor_floor, major_fraction=major)
        result.f2 = stats.f2
        result.f47 = stats.f47
        result.lineage_het_sites = stats.het_sites
        result.lineage_sites_examined = stats.sites_examined
        result.mixed_lineages = list(stats.mixed_lineages)
        if stats.note:
            result.caveats.append(stats.note)

    if unambiguous_sites is not None:
        result.unambiguous_fraction = unambiguous_fraction(unambiguous_sites)
    elif unambiguous_fraction_value is not None:
        result.unambiguous_fraction = float(unambiguous_fraction_value)

    result.mixture_class, result.mixture_reason = classify_mixture(
        result.het_snp_fraction, result.f2, result.f47,
        platform=plat, config=config)

    if plat == PLATFORM_ONT:
        result.caveats.append(ONT_MINOR_VARIANT_CAVEAT)

    result.checks = heterozygosity_checks(result, plat)
    return result
