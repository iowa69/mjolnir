"""Read-level composition, honestly bounded — and the sample-validity verdict.

The headline this module produces is a verdict on whether the sample can carry
the conclusions the run is about to draw from it, not a purity percentage. The
design says why in one sentence: a sample that was 99.84% *M. tuberculosis*
still produced 13 false-positive SNPs across 12 genes, and 5% *M. avium*
contamination produced 3,325. A percentage invites the reader to compare against
a gate, and every gate at 1% or 5% is a coarse instrument.

The verdict therefore asks *valid for what*. Resistance calling reads a few
hundred catalogued positions and tolerates a great deal — a contaminant that
adds thousands of spurious SNPs across the genome may touch none of them, and
the drug-level ``target_covered`` flag catches the case where it does. An
outbreak SNP distance reads the whole masked genome and tolerates almost
nothing, because 3,325 false-positive variants against a 5- or 12-SNP clustering
threshold is not noise, it is the entire signal. So :func:`sample_validity`
returns a verdict per intended use, and the headline is the worst of the uses
the run actually asked for.

What is measured here is the read-level side: mapped fraction, coverage breadth,
coverage evenness, GC, and — only when a mycobacterial ANI reference set is
present — the non-target read fraction, labelled by what the assignment could
actually resolve. The allele-fraction side lives in ``heterozygosity.py``.

**Four refusals are enforced by the shape of this module, not by its prose.**

1. A taxonomic classifier row is never a species identification for this genus.
   Every label leaving this module passes through :func:`taxon_label_for_report`,
   which collapses MTBC members to the complex: in NCBI taxonomy *M. bovis*
   (taxid 1765) has rank ``no rank`` under *M. tuberculosis*, so "M. bovis 3.2%"
   is not a finding, it is a k-mer artefact with a species name attached.
2. A Kraken2 run against a standard, PlusPF or capped index is not a
   contamination screen for mycobacteria. :class:`TaxonomicScreen` carries
   ``informative`` and its :meth:`~TaxonomicScreen.reportable_rows` returns
   nothing when the index cannot support the statement — so the only thing an
   uninformative screen can contribute to a report is the sentence saying it
   could not be run meaningfully. Measured sensitivity for *M. tuberculosis*
   reads on a standard index is 0.0731 on real Illumina data: about 93% of true
   target reads lost, against ~0.97 with a mycobacterial pangenome database.
3. Kraken2's ``--confidence`` default of 0.0 is refused in ``config.py``, at the
   only point every caller passes through, so no path reaches a 0.0 screen.
4. CheckM, CheckM2 and ConFindr are refused as same-species mixture detectors by
   :func:`assert_mixture_method_supported`, which raises rather than warns.
   CheckM's contamination figure is a multi-copy marker-gene statistic and
   ConFindr has no mycobacterial scheme at all.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ..config import (
    KRAKEN2_REPORTING_FLOOR,
    MAX_NON_TARGET_FRACTION as _MAX_NON_TARGET_FRACTION,
    ANI_MIN_ALIGNED_FRACTION,
    ANI_SPECIES_FLOOR,
    COMPLEX_MAC,
    COMPLEX_MTBC,
    CONTAMINATION_EVIDENCE,
    Config,
    EVENNESS_DEFINITION,
    FASTA_CAPABILITY_LOSS,
    GC_TOLERANCE,
    H37RV_GC,
    KRAKEN2_MTB_SENSITIVITY_PANGENOME,
    KRAKEN2_MTB_SENSITIVITY_STANDARD,
    MAC_SPECIES_ANI_FLOOR,
    MIN_BREADTH,
    MIN_COVERAGE_EVENNESS,
    MIN_MAPPED_FRACTION,
    MIN_UNAMBIGUOUS_FRACTION,
    MTBC_UNRESOLVED_TEXT,
    SPECIES_METHOD_REFUSAL,
    kraken2_confidence,
    kraken2_index_informative,
    source_for,
)
from ..records import (
    MIXTURE_MIXED,
    MIXTURE_NOT_ASSESSED,
    MIXTURE_POSSIBLE,
    MIXTURE_SINGLE,
    PLATFORM_FASTA,
    STATUS_FAIL,
    STATUS_PASS,
    STATUS_WARN,
    VALIDITY_INVALID,
    VALIDITY_NOT_ASSESSED,
    VALIDITY_SUSPECT,
    VALIDITY_VALID,
    Check,
    ContaminationResult,
    QCMetrics,
    normalise_platform,
)
from ..utils import LOG, MjolnirError, PathLike, safe_fraction
from .heterozygosity import (
    HeterozygosityResult,
    LineageSite,
    SiteObservation,
    assess_heterozygosity,
)

# ---------------------------------------------------------------------------
# Intended use — the axis the verdict turns on
# ---------------------------------------------------------------------------

#: Reading a few hundred catalogued positions. Tolerant: a contaminant that
#: never touches a resistance locus does not change a drug call, and the ones
#: that do are caught by the per-drug target-coverage flag.
USE_RESISTANCE = "resistance"

#: Reading the whole masked genome for pairwise SNP distances. Intolerant: the
#: measured consequence of 5% M. avium was 3,325 false-positive variant SNPs,
#: against clustering thresholds of 5 and 12.
USE_TRANSMISSION = "transmission"

INTENDED_USES: Tuple[str, ...] = (USE_RESISTANCE, USE_TRANSMISSION)

USE_LABELS: Dict[str, str] = {
    USE_RESISTANCE: "resistance calling",
    USE_TRANSMISSION: "outbreak SNP distances and clustering",
}

#: Validity ordering, worst last. ``not-assessed`` sits above ``valid`` and not
#: beside it: a sample nobody assessed is not a sample that passed, and the
#: whole point of the ordering is that ``worst`` can never fold an unmeasured
#: dimension into a clean headline.
_VALIDITY_SEVERITY: Dict[str, int] = {
    VALIDITY_VALID: 0,
    VALIDITY_NOT_ASSESSED: 1,
    VALIDITY_SUSPECT: 2,
    VALIDITY_INVALID: 3,
}

#: How a validity verdict renders as a check status. ``not-assessed`` warns and
#: is emitted as an unmeasured check, so it can never appear as a pass.
_VALIDITY_CHECK_STATUS: Dict[str, str] = {
    VALIDITY_VALID: STATUS_PASS,
    VALIDITY_NOT_ASSESSED: STATUS_WARN,
    VALIDITY_SUSPECT: STATUS_WARN,
    VALIDITY_INVALID: STATUS_FAIL,
}


def worst_validity(verdicts: Sequence[str]) -> str:
    """The most severe validity verdict in a set; ``not-assessed`` when empty."""
    chosen = VALIDITY_NOT_ASSESSED
    if not verdicts:
        return chosen
    chosen = VALIDITY_VALID
    for verdict in verdicts:
        if verdict not in _VALIDITY_SEVERITY:
            raise MjolnirError(
                "unknown sample-validity verdict {0!r}".format(verdict))
        if _VALIDITY_SEVERITY[verdict] > _VALIDITY_SEVERITY[chosen]:
            chosen = verdict
    return chosen


# ---------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------

#: Tools that are sometimes pointed at this problem and cannot answer it. The
#: value is the reason, printed verbatim when one of them is requested, because
#: "unsupported" without the reason is how the same suggestion comes back.
REFUSED_MIXTURE_METHODS: Dict[str, str] = {
    "checkm": (
        "CheckM's contamination figure is a multi-copy marker-gene statistic: it "
        "measures how often single-copy marker genes appear more than once, which "
        "detects a second *genome* and says nothing about a second strain of the "
        "same species. It is not a same-species mixture detector"),
    "checkm2": (
        "CheckM2 reports the same marker-gene-derived contamination quantity as "
        "CheckM, from a machine-learning model, and inherits the same limitation: "
        "it is not a same-species mixture detector"),
    "confindr": (
        "ConFindr has no mycobacterial scheme and is validated only against rMLST "
        "databases for three genera, so it cannot be applied to this genus"),
}

MIXTURE_METHOD_ALTERNATIVE = (
    "use mjolnir.contamination.heterozygosity instead: F2/F47 across "
    "lineage-defining positions and the genome-wide heterozygous-SNP fraction "
    "under the MixInfect filters measure a same-species mixture directly"
)


def assert_mixture_method_supported(name: str) -> None:
    """Raise when asked to detect a mixture with a tool that cannot.

    A raise rather than a warning. Each of these tools returns a number that
    looks exactly like a contamination estimate, and a warning next to a number
    loses to the number every time it reaches a report.
    """
    key = str(name or "").strip().lower().replace("-", "").replace("_", "")
    reason = REFUSED_MIXTURE_METHODS.get(key)
    if reason is not None:
        raise MjolnirError(
            "{0} cannot be used as a mixture detector for mycobacteria: {1}.\n"
            "  {2}".format(name, reason, MIXTURE_METHOD_ALTERNATIVE)
        )


# --- MTBC members must never be printed as a classifier's species call -------
#
# SOURCE: design §6 and NCBI taxonomy. MTBC members are later heterotypic
# synonyms of M. tuberculosis and sit at 99.21-99.92% ANI of one another;
# `M. tuberculosis variant bovis` (taxid 1765) has rank `no rank`. No read
# classifier and no ANI value resolves inside that, so no method in Mjolnir
# prints a member name derived from one.

MTBC_MEMBER_MARKERS: Tuple[str, ...] = (
    "tuberculosis", "bovis", "bcg", "africanum", "canettii", "microti",
    "caprae", "orygis", "pinnipedii", "mungi", "suricattae", "dassie",
    "tuberculosis complex",
)

MTBC_CLASSIFIER_LABEL = "M. tuberculosis complex (member not resolved)"


def is_mtbc_member_name(name: str) -> bool:
    """Whether a taxon label names an MTBC member or the complex itself."""
    lowered = " ".join(str(name or "").split()).lower()
    if not lowered:
        return False
    if "mycobact" not in lowered and "bcg" not in lowered:
        return False
    return any(marker in lowered for marker in MTBC_MEMBER_MARKERS)


def taxon_label_for_report(name: str) -> str:
    """The only form in which a taxon label may leave this module.

    MTBC members collapse to the complex. Everything else passes through
    unchanged — the refusal is specific to the complex whose members are not
    separable by the methods that produced the label, not a general distrust of
    taxon names.
    """
    if is_mtbc_member_name(name):
        return MTBC_CLASSIFIER_LABEL
    return str(name or "").strip()


# ---------------------------------------------------------------------------
# The taxonomic screen, and what it is allowed to say
# ---------------------------------------------------------------------------

SCREEN_INFORMATIVE = "informative"
SCREEN_UNINFORMATIVE = "uninformative"

SCREEN_UNINFORMATIVE_HEADLINE = (
    "the taxonomic contamination screen could not be run meaningfully on this "
    "sample; no statement about mycobacterial purity is made from it. Measured "
    "sensitivity for M. tuberculosis reads is {0} on a standard index against "
    "about {1} on a mycobacterial pangenome index.".format(
        KRAKEN2_MTB_SENSITIVITY_STANDARD, KRAKEN2_MTB_SENSITIVITY_PANGENOME)
)


@dataclass
class TaxonomicScreen:
    """A read-classifier screen, and whether it can support any statement.

    ``rows`` is kept even when the screen is uninformative, because a support
    question about why a run said nothing is answered by looking at what the
    classifier actually returned. It is reached through
    :meth:`reportable_rows`, which returns nothing in that case, so the report
    writer cannot print an uninformative tail of low-abundance NTM by accident.
    """

    method: str = ""
    informative: bool = False
    index: str = ""
    confidence: Optional[float] = None
    note: str = ""
    rows: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def status(self) -> str:
        """``informative`` or ``uninformative`` — never ``pass``."""
        return SCREEN_INFORMATIVE if self.informative else SCREEN_UNINFORMATIVE

    def reportable_rows(self) -> List[Dict[str, Any]]:
        """Classifier rows the report may print, which is none when uninformative."""
        if not self.informative:
            return []
        return list(self.rows)

    def to_check(self) -> Check:
        """The panel row for this screen.

        An uninformative screen is ``not_measured``: it warns, it carries the
        reason, and it can never be rendered as a pass. That is the whole
        mechanism by which "we could not look" stops turning into "we looked and
        it was clean".
        """
        if not self.informative:
            return Check.not_measured(
                "taxonomic_contamination_screen",
                "{0} {1}".format(SCREEN_UNINFORMATIVE_HEADLINE, self.note).strip(),
                source=source_for("kraken2_mtb_sensitivity_standard"),
                category="contamination")
        # The screen having been *possible* is not the screen having been
        # *clean*. This branch used to pass unconditionally without reading a
        # single row; the first attempt to fix that read keys the parser does
        # not emit - is_target, fraction, name - so every lookup returned None,
        # the largest foreign share was always 0.0, and a library that was 99%
        # Staphylococcus aureus still passed. The keys below are the ones
        # parse_kraken2_report actually writes.
        detail = "{0} against {1} at --confidence {2}".format(
            self.method or "taxonomic screen", self.index or "an index",
            self.confidence)
        foreign = _foreign_taxa(self.reportable_rows())
        share = max((frac for _label, frac in foreign), default=0.0)
        if not foreign:
            return Check.boolean(
                "taxonomic_contamination_screen", True, expected=True,
                source=source_for("kraken2_mtb_sensitivity_pangenome"),
                category="contamination",
                reading="{0}: no non-target taxon above {1:.1%} of reads.".format(
                    detail, KRAKEN2_REPORTING_FLOOR))
        names = ", ".join("{0} {1:.1%}".format(label, frac)
                          for label, frac in foreign[:3])
        return Check.boolean(
            "taxonomic_contamination_screen", share <= MAX_NON_TARGET_FRACTION,
            expected=True,
            source=source_for("kraken2_mtb_sensitivity_pangenome"),
            category="contamination",
            reading="{0}: non-target reads present ({1}). The largest is "
                    "{2:.1%} against a {3:.1%} ceiling.".format(
                        detail, names, share, MAX_NON_TARGET_FRACTION))


def _foreign_taxa(rows: Sequence[Dict[str, Any]]) -> List[Tuple[str, float]]:
    """Named non-target species in a Kraken2 report, largest share first.

    Restricted to species-rank rows. A Kraken2 report is hierarchical: its
    ``root`` and ``Bacteria`` lines sit at 100%, and counting those as foreign
    taxa would fail every sample ever screened.

    ``percentage`` is 0-100 and every threshold here is a 0-1 fraction, so the
    conversion is not cosmetic - comparing 99.1 against a ceiling of 0.10 and
    comparing 0.991 against it give opposite answers.
    """
    found: List[Tuple[str, float]] = []
    for row in rows:
        if not str(row.get("rank") or "").startswith("S"):
            continue
        if row.get("collapsed_to_complex") or row.get("label") == MTBC_CLASSIFIER_LABEL:
            continue
        percentage = row.get("percentage")
        if percentage is None:
            continue
        fraction = float(percentage) / 100.0
        if fraction < KRAKEN2_REPORTING_FLOOR:
            continue
        found.append((str(row.get("label") or "unnamed"), fraction))
    found.sort(key=lambda pair: -pair[1])
    return found


def parse_kraken2_report(text: str) -> List[Dict[str, Any]]:
    """Parse a ``kraken2 --report`` table into rows.

    Six tab-separated columns: clade percentage, clade fragments, direct
    fragments, rank code, taxid, indented name. Rows that do not have them are
    counted and reported in one error rather than skipped, since a silently
    half-parsed report is a composition estimate computed over an unknown
    denominator.

    Every name is passed through :func:`taxon_label_for_report` on the way out,
    so an MTBC member cannot leave this function under a member name whatever
    the caller does with the rows afterwards.
    """
    rows: List[Dict[str, Any]] = []
    malformed = 0
    for line in str(text or "").splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            malformed += 1
            continue
        try:
            percentage = float(parts[0].strip())
            clade_fragments = int(parts[1].strip())
            direct_fragments = int(parts[2].strip())
        except ValueError:
            malformed += 1
            continue
        raw_name = parts[5].strip()
        rows.append({
            "percentage": percentage,
            "clade_fragments": clade_fragments,
            "direct_fragments": direct_fragments,
            "rank": parts[3].strip(),
            "taxid": parts[4].strip(),
            "label": taxon_label_for_report(raw_name),
            "collapsed_to_complex": is_mtbc_member_name(raw_name),
        })
    if malformed and not rows:
        raise MjolnirError(
            "could not parse any row from the Kraken2 report ({0} malformed "
            "lines); expected six tab-separated columns per line".format(malformed)
        )
    if malformed:
        LOG.warning("ignored %d malformed Kraken2 report lines", malformed)
    return rows


def evaluate_kraken2_screen(db_dir: Optional[PathLike],
                            confidence: Optional[float] = None,
                            report_text: Optional[str] = None) -> TaxonomicScreen:
    """Decide what, if anything, a Kraken2 run here is allowed to contribute.

    The index is interrogated before the output is: a standard, PlusPF or capped
    index has a measured sensitivity of 0.0731 for *M. tuberculosis* reads on
    real Illumina data, so its rows are not evidence about this genus and the
    screen is marked uninformative regardless of how clean they look. An index
    that declares itself a mycobacterial pangenome database (~0.97) is
    informative, and its rows still pass the MTBC collapse.

    ``confidence`` goes through :func:`~mjolnir.config.kraken2_confidence`,
    which raises on Kraken2's own 0.0 default rather than accepting it.
    """
    conf = kraken2_confidence(confidence)
    informative, note = kraken2_index_informative(db_dir)
    screen = TaxonomicScreen(
        method="kraken2",
        informative=informative,
        index=str(db_dir) if db_dir is not None else "",
        confidence=conf,
        note=note,
    )
    if report_text:
        screen.rows = parse_kraken2_report(report_text)
    if not informative:
        LOG.info("Kraken2 screen treated as uninformative: %s", note)
    return screen


def no_screen(reason: str = "no taxonomic contamination screen was run") -> TaxonomicScreen:
    """The screen object for a run that did not perform one.

    A run without a screen and a run with a useless screen produce the same
    thing here — an object whose ``informative`` is False — so neither can reach
    the report as an absence of contamination.
    """
    return TaxonomicScreen(method="", informative=False, note=reason)


# ---------------------------------------------------------------------------
# Read-level composition
# ---------------------------------------------------------------------------

@dataclass
class PurityPanel:
    """Mapped fraction, breadth, evenness and GC, with their checks."""

    mapped_fraction: Optional[float] = None
    coverage_breadth: Optional[float] = None
    coverage_evenness: Optional[float] = None
    evenness_definition: str = EVENNESS_DEFINITION
    gc_content: Optional[float] = None
    reference_gc: Optional[float] = None
    gc_within_tolerance: Optional[bool] = None
    checks: List[Check] = field(default_factory=list)


def measure_purity(qc: Optional[QCMetrics], platform: str,
                   reference_gc: Optional[float] = None,
                   config: Optional[Config] = None) -> PurityPanel:
    """The read-level signals bearing on purity, mirrored from the QC metrics.

    These are the same signals MTBseq's ``TBstats`` emits, which are its only
    purity-adjacent output, and they are mirrored into the contamination panel
    rather than referenced across so the panel is complete on its own in the
    JSON artefact and in what the agent is shown.

    None of them is a contamination measurement on its own. A mapped fraction of
    0.85 says 15% of the library is something else *or* that the reference is
    the wrong one, and only the non-target assignment can tell those apart — so
    each check states its threshold and the verdict weighs them together.
    """
    plat = normalise_platform(platform)
    panel = PurityPanel(reference_gc=reference_gc)
    min_breadth = config.min_breadth if config else MIN_BREADTH

    if qc is not None:
        panel.mapped_fraction = qc.mapped_fraction
        panel.coverage_breadth = (
            qc.breadth_min_depth if qc.breadth_min_depth is not None else qc.breadth_10x)
        panel.coverage_evenness = qc.coverage_evenness
        panel.gc_content = qc.gc_content
        if qc.evenness_definition:
            panel.evenness_definition = qc.evenness_definition

    fasta_note = FASTA_CAPABILITY_LOSS if plat == PLATFORM_FASTA else ""
    if fasta_note:
        no_mapping = fasta_note
    elif qc is None:
        no_mapping = "no QC metrics were supplied for this sample, so the "\
                     "read-level composition signals were not measured"
    else:
        no_mapping = "no reads were mapped to the reference, so this was not measured"

    panel.checks.append(Check.numeric(
        "mapped_read_fraction", panel.mapped_fraction,
        warn_minimum=MIN_MAPPED_FRACTION,
        source=source_for("min_mapped_fraction"),
        unit="fraction of reads", category="contamination",
        reading="reads mapping to the chosen reference; a low value is either a "
                "non-target library or the wrong reference, and the non-target "
                "assignment is what separates them",
        not_measured_why=no_mapping))

    panel.checks.append(Check.numeric(
        "coverage_breadth", panel.coverage_breadth,
        warn_minimum=min_breadth,
        source=source_for("min_breadth"),
        unit="fraction of reference", category="contamination",
        reading="fraction of the reference reaching the configured depth floor; "
                "genome-wide statements need this before they mean anything",
        not_measured_why=no_mapping))

    panel.checks.append(Check.numeric(
        "coverage_evenness", panel.coverage_evenness,
        warn_minimum=MIN_COVERAGE_EVENNESS,
        source=source_for("min_coverage_evenness"),
        unit="fraction of positions", category="contamination",
        reading=panel.evenness_definition,
        not_measured_why=no_mapping))

    expected_gc = reference_gc
    if expected_gc is None and qc is not None and qc.reference:
        # Only H37Rv's composition is known without loading the reference; an NTM
        # reference carries its own value in the database registry and the caller
        # passes it in. Guessing H37Rv's 65.6% for an NTM sample would fire a GC
        # warning on every M. abscessus run, which teaches a reader to ignore it.
        if str(qc.reference).upper().startswith("NC_000962"):
            expected_gc = H37RV_GC
    panel.reference_gc = expected_gc

    if panel.gc_content is None or expected_gc is None:
        panel.checks.append(Check.not_measured(
            "gc_content",
            fasta_note or (
                "GC was not compared: {0}".format(
                    "no GC content was measured" if panel.gc_content is None
                    else "the reference's own GC content is not recorded, so there "
                         "is nothing to compare against")),
            source=source_for("gc_tolerance"), category="contamination"))
    else:
        panel.gc_within_tolerance = abs(panel.gc_content - expected_gc) <= GC_TOLERANCE
        panel.checks.append(Check.numeric(
            "gc_content", panel.gc_content,
            warn_minimum=expected_gc - GC_TOLERANCE,
            warn_maximum=expected_gc + GC_TOLERANCE,
            source=source_for("gc_tolerance"),
            unit="fraction", category="contamination",
            reading="observed GC against the reference's {0:.3f} +/- {1}; a "
                    "composition shift is a flag, never a species claim".format(
                        expected_gc, GC_TOLERANCE)))
    return panel


# ---------------------------------------------------------------------------
# Non-target reads, by ANI assignment only
# ---------------------------------------------------------------------------

#: Derived from ``config.MIN_MAPPED_FRACTION``, not an independent threshold:
#: the reads that did not map to the chosen reference are exactly the reads a
#: non-target assignment is trying to account for, so the two numbers are one
#: number and must not drift apart.
#: Imported from the registry rather than derived here, so it carries a
#: source and appears in `all_thresholds()` like every other number.
MAX_NON_TARGET_FRACTION = _MAX_NON_TARGET_FRACTION

RESOLUTION_SPECIES = "species"
RESOLUTION_COMPLEX = "complex"
RESOLUTION_GENUS = "genus"

MAC_UNRESOLVED_TEXT = (
    "within MAC, ANI at or above {0}% is necessary but not sufficient to name a "
    "species: M. chimaera and M. intracellulare are separated by marker SNPs, "
    "not by ANI, so a MAC component is labelled at complex level here".format(
        MAC_SPECIES_ANI_FLOOR)
)

NO_ANI_REFERENCE_SET = (
    "no mycobacterial ANI reference set is installed, so the non-target read "
    "fraction was not measured. This is an absent measurement, not an absence "
    "of contamination"
)


@dataclass
class AniAssignment:
    """One reference in the ANI set and how much of the library went to it.

    ``fraction`` and ``reads`` are both optional and at least one must be
    present, because the caller may have either a per-reference read count or a
    precomputed share; ``label`` is the taxon name, which leaves this module
    only through :func:`taxon_label_for_report`.
    """

    label: str
    is_target: bool = False
    reads: Optional[int] = None
    fraction: Optional[float] = None
    ani: Optional[float] = None
    aligned_fraction: Optional[float] = None
    complex_name: str = ""

    def __post_init__(self) -> None:
        if self.reads is None and self.fraction is None:
            raise MjolnirError(
                "ANI assignment for {0!r} carries neither a read count nor a "
                "fraction, so it cannot contribute to a composition".format(self.label)
            )

    @property
    def resolves_to_species(self) -> bool:
        """Whether this assignment is strong enough to name a species.

        Three conditions, and the third is the one that matters here: ANI at or
        above the prokaryotic species boundary, over enough aligned genome to
        mean it, and outside the complexes whose members ANI cannot separate.

        MAC is excluded outright, not gated on the higher MAC floor.
        *M. chimaera* and *M. intracellulare* sit above that floor of each other,
        so clearing it is necessary and not sufficient — the marker SNPs that
        finish the job live in ``typing/species.py``, and a composition table has
        none of them. A contaminant is named "MAC" here or it is not named.
        """
        if self.ani is None or self.ani < ANI_SPECIES_FLOOR:
            return False
        if self.aligned_fraction is not None and self.aligned_fraction < ANI_MIN_ALIGNED_FRACTION:
            return False
        if is_mtbc_member_name(self.label) or self.complex_name.upper() == COMPLEX_MTBC:
            return False
        if self.complex_name.upper() == COMPLEX_MAC:
            return False
        return True


@dataclass
class NonTargetResult:
    fraction: Optional[float] = None
    resolution: str = ""
    labels: List[Dict[str, Any]] = field(default_factory=list)
    measured: bool = False
    note: str = ""
    checks: List[Check] = field(default_factory=list)


def assess_non_target(assignments: Optional[Sequence[AniAssignment]],
                      reference_set_present: bool = False) -> NonTargetResult:
    """Non-target read fraction by ANI assignment, when that is possible at all.

    The design permits this measurement "only when a suitable mycobacterial
    reference set is present", so an absent set produces ``measured=False`` and
    the sentence saying so — not a fraction of zero, which would be a claim that
    nothing foreign was found by a method that was never run.

    ``resolution`` labels what the assignment could actually resolve. Inside the
    MTBC it is ``complex`` and can be nothing else; inside MAC it is ``complex``
    unless ANI clears the higher floor that separates *M. chimaera* from
    *M. intracellulare*, which by itself is still not the marker-SNP evidence
    ``typing/species.py`` requires before naming one. A number labelled
    ``species`` and a number labelled ``complex`` are different statements and
    the report prints the label beside the number.
    """
    result = NonTargetResult()
    if not reference_set_present or not assignments:
        result.note = NO_ANI_REFERENCE_SET if not reference_set_present else (
            "the ANI reference set is installed but produced no assignment for "
            "this sample, so the non-target read fraction was not measured")
        result.checks.append(Check.not_measured(
            "non_target_read_fraction", result.note,
            source=source_for("ani_species_floor"), category="contamination"))
        return result

    total_reads = sum(a.reads or 0 for a in assignments)
    non_target = [a for a in assignments if not a.is_target]

    if total_reads:
        result.fraction = safe_fraction(sum(a.reads or 0 for a in non_target), total_reads)
    else:
        shares = [a.fraction for a in non_target if a.fraction is not None]
        result.fraction = sum(shares) if shares else None

    if result.fraction is None:
        result.note = (
            "ANI assignments carried neither read counts nor shares, so the "
            "non-target read fraction was not computed")
        result.checks.append(Check.not_measured(
            "non_target_read_fraction", result.note,
            source=source_for("ani_species_floor"), category="contamination"))
        return result

    result.measured = True
    resolutions = set()
    for assignment in assignments:
        label = taxon_label_for_report(assignment.label)
        resolved = assignment.resolves_to_species
        if resolved:
            resolutions.add(RESOLUTION_SPECIES)
        elif (is_mtbc_member_name(assignment.label)
              or assignment.complex_name.upper() in (COMPLEX_MTBC, COMPLEX_MAC.upper())):
            resolutions.add(RESOLUTION_COMPLEX)
        else:
            resolutions.add(RESOLUTION_GENUS)
        share = assignment.fraction
        if share is None and total_reads:
            share = safe_fraction(assignment.reads or 0, total_reads)
        result.labels.append({
            "label": label,
            "is_target": bool(assignment.is_target),
            "fraction": share,
            "reads": assignment.reads,
            "ani": assignment.ani,
            "aligned_fraction": assignment.aligned_fraction,
            "resolved_to_species": resolved,
            "complex": assignment.complex_name,
        })

    # The weakest resolution any assignment reached is the resolution of the
    # whole composition: one unresolvable component makes the composition
    # unresolvable, and reporting the best of them would overstate all of it.
    for candidate in (RESOLUTION_GENUS, RESOLUTION_COMPLEX, RESOLUTION_SPECIES):
        if candidate in resolutions:
            result.resolution = candidate
            break

    if result.resolution in (RESOLUTION_COMPLEX, RESOLUTION_GENUS):
        if any(is_mtbc_member_name(a.label) for a in assignments):
            result.note = MTBC_UNRESOLVED_TEXT
        elif any(a.complex_name.upper() == COMPLEX_MAC for a in assignments):
            result.note = MAC_UNRESOLVED_TEXT
        else:
            result.note = "ANI resolved this composition only to {0} level".format(
                result.resolution)

    result.checks.append(Check.numeric(
        "non_target_read_fraction", result.fraction,
        warn_maximum=MAX_NON_TARGET_FRACTION,
        source=source_for("min_mapped_fraction"),
        unit="fraction of reads", category="contamination",
        reading="assigned by ANI, resolved to {0} level. {1}".format(
            result.resolution, CONTAMINATION_EVIDENCE)))
    return result


# ---------------------------------------------------------------------------
# The verdict
# ---------------------------------------------------------------------------

@dataclass
class _Dimension:
    """One line of evidence, and what it means for each intended use."""

    name: str
    reason: str
    resistance: str = VALIDITY_VALID
    transmission: str = VALIDITY_VALID

    def for_use(self, use: str) -> str:
        if use == USE_RESISTANCE:
            return self.resistance
        if use == USE_TRANSMISSION:
            return self.transmission
        raise MjolnirError(
            "unknown intended use {0!r}; expected one of {1}".format(
                use, ", ".join(INTENDED_USES))
        )


@dataclass
class ValidityVerdict:
    """The headline, per intended use and folded.

    ``by_use`` is the part that matters clinically. The same sample is routinely
    valid for one question and not for the other, and a single word cannot say
    that — which is why the report prints "valid for resistance calling; invalid
    for outbreak SNP distances" rather than choosing.
    """

    verdict: str = VALIDITY_NOT_ASSESSED
    reason: str = ""
    by_use: Dict[str, str] = field(default_factory=dict)
    by_use_reason: Dict[str, str] = field(default_factory=dict)
    caveats: List[str] = field(default_factory=list)

    def sentence(self) -> str:
        """One line naming the verdict for each use it was asked about."""
        parts = ["{0} for {1}".format(self.by_use[use], USE_LABELS.get(use, use))
                 for use in INTENDED_USES if use in self.by_use]
        return "; ".join(parts) if parts else self.verdict


def _mixture_dimension(mixture_class: str, reason: str) -> _Dimension:
    """How a mixture class bears on each use.

    A called mixture invalidates both: two strains make a drug call ambiguous
    and a distance meaningless. A *possible* mixture parts them, and this is the
    core of the design's §8 argument — 99.84% pure still yielded 13
    false-positive SNPs, which no resistance call noticed and which is more than
    the entire 12-SNP clustering threshold.
    """
    if mixture_class == MIXTURE_MIXED:
        return _Dimension("mixture", reason, VALIDITY_INVALID, VALIDITY_INVALID)
    if mixture_class == MIXTURE_POSSIBLE:
        return _Dimension("mixture", reason, VALIDITY_SUSPECT, VALIDITY_INVALID)
    if mixture_class == MIXTURE_SINGLE:
        return _Dimension("mixture", reason, VALIDITY_VALID, VALIDITY_VALID)
    return _Dimension("mixture", reason, VALIDITY_NOT_ASSESSED, VALIDITY_NOT_ASSESSED)


def sample_validity(*,
                    platform: str,
                    mixture_class: str = MIXTURE_NOT_ASSESSED,
                    mixture_reason: str = "",
                    purity: Optional[PurityPanel] = None,
                    non_target: Optional[NonTargetResult] = None,
                    screen: Optional[TaxonomicScreen] = None,
                    unambiguous_fraction: Optional[float] = None,
                    intended_use: Sequence[str] = INTENDED_USES,
                    config: Optional[Config] = None) -> ValidityVerdict:
    """The sample-validity verdict, per intended use.

    Every dimension of evidence votes once per use, and the worst vote wins —
    including ``not-assessed``, which outranks ``valid`` so that a sample nobody
    could assess never reads as one that passed. That is why a FASTA input can
    never be reported valid: an assembly carries no allele fractions, so the
    mixture dimension is permanently unassessed, and saying "valid" would be
    saying an unmeasured thing was fine.

    A taxonomic screen that could not be run meaningfully does not lower the
    verdict — it contributed no evidence either way — but it is always added as
    a caveat, so a reader can never mistake the silence for a clean result.
    """
    plat = normalise_platform(platform)
    uses = [use for use in intended_use]
    for use in uses:
        if use not in INTENDED_USES:
            raise MjolnirError(
                "unknown intended use {0!r}; expected one or more of {1}".format(
                    use, ", ".join(INTENDED_USES)))
    if not uses:
        raise MjolnirError(
            "sample_validity needs at least one intended use; a validity verdict "
            "with no question attached cannot be answered")

    min_breadth = config.min_breadth if config else MIN_BREADTH
    dimensions: List[_Dimension] = [_mixture_dimension(mixture_class, mixture_reason)]

    if purity is not None:
        mapped = purity.mapped_fraction
        if mapped is None:
            dimensions.append(_Dimension(
                "mapped_read_fraction",
                "the mapped-read fraction was not measured",
                VALIDITY_NOT_ASSESSED, VALIDITY_NOT_ASSESSED))
        elif mapped < MIN_MAPPED_FRACTION:
            dimensions.append(_Dimension(
                "mapped_read_fraction",
                "only {0:.1%} of reads mapped to the reference, below the {1:.0%} "
                "floor".format(mapped, MIN_MAPPED_FRACTION),
                VALIDITY_SUSPECT, VALIDITY_INVALID))

        breadth = purity.coverage_breadth
        if breadth is None:
            dimensions.append(_Dimension(
                "coverage_breadth",
                "coverage breadth was not measured",
                VALIDITY_VALID, VALIDITY_NOT_ASSESSED))
        elif breadth < min_breadth:
            dimensions.append(_Dimension(
                "coverage_breadth",
                "{0:.1%} of the reference reached the depth floor, below the "
                "{1:.0%} needed for genome-wide statements".format(breadth, min_breadth),
                VALIDITY_SUSPECT, VALIDITY_INVALID))

        if purity.coverage_evenness is not None and \
                purity.coverage_evenness < MIN_COVERAGE_EVENNESS:
            dimensions.append(_Dimension(
                "coverage_evenness",
                "coverage is uneven ({0:.1%} of positions inside the band), which "
                "is what a second organism at partial depth looks like".format(
                    purity.coverage_evenness),
                VALIDITY_SUSPECT, VALIDITY_SUSPECT))

        if purity.gc_within_tolerance is False:
            dimensions.append(_Dimension(
                "gc_content",
                "GC content sits outside {0} of the reference's, a composition "
                "shift consistent with foreign DNA".format(GC_TOLERANCE),
                VALIDITY_SUSPECT, VALIDITY_SUSPECT))

    if non_target is not None and non_target.measured and non_target.fraction is not None:
        if non_target.fraction > MAX_NON_TARGET_FRACTION:
            dimensions.append(_Dimension(
                "non_target_read_fraction",
                "{0:.1%} of reads were assigned to non-target references at {1} "
                "resolution. {2}".format(
                    non_target.fraction, non_target.resolution or "unstated",
                    CONTAMINATION_EVIDENCE),
                VALIDITY_SUSPECT, VALIDITY_INVALID))

    if unambiguous_fraction is not None:
        if unambiguous_fraction < MIN_UNAMBIGUOUS_FRACTION:
            dimensions.append(_Dimension(
                "unambiguous_base_fraction",
                "{0:.1%} of positions carried an unambiguous majority allele, "
                "below {1:.0%}; MTBseq would have discarded the disagreeing "
                "reads".format(unambiguous_fraction, MIN_UNAMBIGUOUS_FRACTION),
                VALIDITY_SUSPECT, VALIDITY_SUSPECT))

    verdict = ValidityVerdict()
    for use in uses:
        votes = [d.for_use(use) for d in dimensions]
        decided = worst_validity(votes)
        verdict.by_use[use] = decided
        reasons = [d.reason for d in dimensions
                   if d.for_use(use) == decided and d.reason]
        verdict.by_use_reason[use] = "; ".join(reasons) if reasons else (
            "every measured contamination signal was within its threshold")

    verdict.verdict = worst_validity(list(verdict.by_use.values()))
    verdict.reason = verdict.sentence()
    detail = "; ".join(
        "{0}: {1}".format(USE_LABELS.get(use, use), verdict.by_use_reason[use])
        for use in uses)
    if detail:
        verdict.reason = "{0}. {1}".format(verdict.reason, detail)

    if screen is not None and not screen.informative:
        verdict.caveats.append(
            "{0} {1}".format(SCREEN_UNINFORMATIVE_HEADLINE, screen.note).strip())
    if plat == PLATFORM_FASTA:
        verdict.caveats.append(FASTA_CAPABILITY_LOSS)
    verdict.caveats.append(CONTAMINATION_EVIDENCE)
    return verdict


# ---------------------------------------------------------------------------
# The whole panel
# ---------------------------------------------------------------------------

def assess_contamination(*,
                         platform: str,
                         qc: Optional[QCMetrics] = None,
                         snp_sites: Optional[Sequence[SiteObservation]] = None,
                         lineage_sites: Optional[Sequence[LineageSite]] = None,
                         unambiguous_sites: Optional[Sequence[SiteObservation]] = None,
                         heterozygosity: Optional[HeterozygosityResult] = None,
                         ani_assignments: Optional[Sequence[AniAssignment]] = None,
                         reference_set_present: bool = False,
                         screen: Optional[TaxonomicScreen] = None,
                         reference_gc: Optional[float] = None,
                         intended_use: Sequence[str] = INTENDED_USES,
                         config: Optional[Config] = None) -> ContaminationResult:
    """Everything §8 permits to be measured, and the verdict it supports.

    The one entry point ``pipeline.py`` needs. It fills the contract record so
    that the PDF, the HTML, the JSON artefact and the agent observation all read
    the same fields — building the panel twice is how two outputs come to
    disagree about the same run.

    ``screen`` defaults to :func:`no_screen`, not to a passing screen. A run that
    was never given a Kraken2 index and a run given a standard one arrive at the
    same place: ``screen_informative`` False, with a note saying which.
    """
    plat = normalise_platform(platform)
    het = heterozygosity if heterozygosity is not None else assess_heterozygosity(
        platform=plat, snp_sites=snp_sites, lineage_sites=lineage_sites,
        unambiguous_sites=unambiguous_sites,
        unambiguous_fraction_value=(qc.unambiguous_fraction
                                    if qc is not None and unambiguous_sites is None
                                    else None),
        config=config)

    purity = measure_purity(qc, plat, reference_gc=reference_gc, config=config)
    non_target = assess_non_target(ani_assignments, reference_set_present)
    used_screen = screen if screen is not None else no_screen()

    verdict = sample_validity(
        platform=plat,
        mixture_class=het.mixture_class,
        mixture_reason=het.mixture_reason,
        purity=purity,
        non_target=non_target,
        screen=used_screen,
        unambiguous_fraction=het.unambiguous_fraction,
        intended_use=intended_use,
        config=config)

    caveats: List[str] = []
    for entry in list(het.caveats) + list(verdict.caveats):
        if entry and entry not in caveats:
            caveats.append(entry)
    if non_target.note and non_target.note not in caveats:
        caveats.append(non_target.note)

    checks: List[Check] = list(het.checks) + list(purity.checks) + list(non_target.checks)
    checks.append(used_screen.to_check())
    if verdict.verdict == VALIDITY_NOT_ASSESSED:
        checks.append(Check.not_measured(
            "sample_validity", verdict.reason,
            source=source_for("contamination_evidence"), category="contamination"))
    else:
        checks.append(Check(
            name="sample_validity",
            value=verdict.verdict,
            threshold=VALIDITY_VALID,
            source=source_for("contamination_evidence"),
            status=_VALIDITY_CHECK_STATUS[verdict.verdict],
            reading=verdict.reason,
            comparison="==", category="contamination"))

    return ContaminationResult(
        f2=het.f2,
        f47=het.f47,
        lineage_het_sites=het.lineage_het_sites,
        lineage_sites_examined=het.lineage_sites_examined,
        het_snp_fraction=het.het_snp_fraction,
        het_snp_count=het.het_snp_count,
        snp_sites_examined=het.snp_sites_examined,
        mixture_class=het.mixture_class,
        unambiguous_fraction=het.unambiguous_fraction,
        non_target_fraction=non_target.fraction,
        non_target_resolution=non_target.resolution,
        non_target_labels=non_target.labels,
        mapped_fraction=purity.mapped_fraction,
        coverage_breadth=purity.coverage_breadth,
        coverage_evenness=purity.coverage_evenness,
        gc_content=purity.gc_content,
        verdict=verdict.verdict,
        verdict_reason=verdict.reason,
        screen_informative=used_screen.informative,
        screen_method=used_screen.method,
        screen_note=used_screen.note,
        caveats=caveats,
        checks=checks,
    )


#: Re-exported for the species and report layers, which must apply the same
#: refusal when they print an ANI hit: no method here names an MTBC member.
SPECIES_REFUSAL_TEXT = SPECIES_METHOD_REFUSAL
