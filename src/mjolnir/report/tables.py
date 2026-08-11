"""The report's data layer: TSV and JSON artefacts, and the rows behind them.

Every panel, table and chart in the PDF and in the HTML is built from a function
in this module. That is not tidiness — it is the only way to keep the promise
:mod:`mjolnir.records` makes in its own docstring, that no two outputs can
disagree about what was found. A number formatted once here appears identically
on page 1 of the PDF, in the HTML, and in the TSV a bioinformatician loads into
R three months later.

Three habits in here are load-bearing.

**Absence is spelled.** A missing value is written ``NA`` in TSV and ``null`` in
JSON, never an empty cell. An empty cell in a tab-separated file is
indistinguishable from an empty string, and one careless reader turns "we could
not measure this" into "this was zero". For resistance the same rule is stricter
still: :data:`~mjolnir.records.CALL_NO_CALL` is rendered with its own glyph and
its own sentence, and the word "susceptible" appears nowhere near it.

**Thresholds travel with their sources.** :func:`qc_panel` and
:func:`contamination_panel` do not read a status off a metric; they apply the
registered threshold from :mod:`mjolnir.config` and hand back a
:class:`~mjolnir.records.Check` that carries the citation. If an engine already
computed the same check, that one wins — but the panel is complete either way,
so a report can never quietly omit a metric because an upstream stage forgot it.

**The row builders are pure.** No timestamps, no random ordering, no filesystem
reads. Two runs over the same :class:`~mjolnir.records.SampleResult` produce
byte-identical rows, which is what makes the golden-file test in the design's
§13 possible at all.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .. import __version__
from ..config import (
    ANI_MIN_ALIGNED_FRACTION,
    ANI_SPECIES_FLOOR,
    BCG_PZA_NOTE,
    CATALOGUES,
    CATALOGUE_MTBSEQ,
    CATALOGUE_WHO,
    COMPLEX_MTBC,
    CONTAMINATION_EVIDENCE,
    DEGRADED_DEPTH_FLOOR,
    DRUGS,
    EVENNESS_DEFINITION,
    F2_MIXTURE_THRESHOLD,
    F47_MIXTURE_THRESHOLD,
    GC_TOLERANCE,
    H37RV_ACCESSION,
    H37RV_GC,
    HET_SNP_FRACTION_MIXED,
    HET_SNP_FRACTION_WARN,
    MAJOR_VARIANT_FRACTION,
    MBOVIS_CAVEAT,
    MIN_BARCODE_CALLABLE_FRACTION,
    MIN_BARCODE_SUPPORT_FRACTION,
    MIN_BREADTH,
    MIN_COVERAGE_EVENNESS,
    MIN_DEPTH,
    MIN_MAPPED_FRACTION,
    MIN_SHARED_CALLABLE_SITES,
    MIN_UNAMBIGUOUS_FRACTION,
    MTBC_UNRESOLVED_TEXT,
    MTBSEQ_ASYMMETRY_NOTE,
    SRC_DESIGN,
    SRC_POLICY,
    all_thresholds,
    drug_code,
    normalise_drug,
    platform_caveats,
    source_for,
    unverified,
)
from ..records import (
    CALL_NO_CALL,
    CALL_R_OUTSIDE_WHO,
    CALL_SEVERITY,
    Check,
    CohortResult,
    DrugCall,
    NO_DETERMINANT_TEXT,
    SampleResult,
    STATUS_FAIL,
    STATUS_PASS,
    STATUS_WARN,
    VALIDITY_INVALID,
    VALIDITY_NOT_ASSESSED,
    VALIDITY_SUSPECT,
    VALIDITY_VALID,
    call_label,
    worst_call,
)
from ..utils import (LOG, MjolnirError, ensure_dir, natural_key, plural,
                     round_or_none, safe_name)

#: What a missing value looks like in a TSV. Not the empty string: an empty cell
#: reads as an empty string to every parser and as "nothing was wrong" to every
#: human, and this tool has one rule it will not bend about that.
TSV_NA = "NA"

#: Columns of the drug grid used by the front-page heatmap: the three catalogues
#: side by side, then the consensus Mjolnir applied over them (design §5.5).
CONSENSUS_COLUMN = "Mjolnir"
GRID_COLUMNS: Tuple[str, ...] = tuple(CATALOGUES) + (CONSENSUS_COLUMN,)

#: Short, ASCII-only glyph for each call. ASCII because the PDF is set in
#: Helvetica under WinAnsi encoding, where a typographic symbol silently becomes
#: a black box, and because the PDF and the HTML must show the same mark.
#:
#: ``ND`` is deliberately two letters rather than a dash: "no resistance
#: determinant detected" is a statement about the catalogue, not about the
#: organism, and a dash in a grid reads as "nothing to see here".
CALL_GLYPH: Dict[str, str] = {
    "R": "R",
    "R-outside-WHO": "R*",
    "R-interim": "Ri",
    "Uncertain": "U",
    "S-interim": "s",
    "S": "S",
    CALL_NO_CALL: "ND",
}

#: Legend text, spelled once, printed by both renderers. The first two entries
#: are the distinction the design refuses to let a reader blur.
CALL_LEGEND: Tuple[Tuple[str, str], ...] = (
    (CALL_NO_CALL, "ND - no resistance determinant detected: nothing catalogued "
                   "was found. This is absence of evidence, not susceptibility."),
    ("S", "S - a variant was found and a catalogue graded it as not associated "
          "with resistance. A positive finding, unlike ND."),
    ("S-interim", "s - graded not associated with resistance, interim grade."),
    ("Uncertain", "U - variant of uncertain significance."),
    ("R-interim", "Ri - resistance determinant, interim grade."),
    ("R-outside-WHO", "R* - resistance called by a catalogue that is not WHO, at a "
                      "variant WHO does not grade. Never equivalent to a WHO Group 1 call."),
    ("R", "R - resistance determinant detected (WHO Group 1)."),
)

#: Marks placed beside a drug on page 1. The disagreement mark is required by the
#: design: a drug whose catalogues conflict must be visible on the front page,
#: not only in the annex.
FLAG_DISAGREEMENT = "!="
FLAG_SUPPRESSED = "(x)"
FLAG_CAVEAT = "*"
FLAG_NOT_EVALUABLE = "NE"


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def fmt_number(value: Optional[float], digits: int = 2, unit: str = "",
               na: str = "not measured") -> str:
    """A number for display, or the words for its absence.

    *na* defaults to "not measured" rather than to a dash, because a dash in a
    clinical table is read as "nothing of note".
    """
    if value is None:
        return na
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int) or float(value).is_integer():
        text = "{0:,}".format(int(value))
    else:
        text = "{0:,.{1}f}".format(float(value), digits)
    return text + unit


def fmt_fraction(value: Optional[float], digits: int = 1,
                 na: str = "not measured") -> str:
    """A 0-1 fraction as a percentage string, preserving absence."""
    if value is None:
        return na
    return "{0:.{1}f}%".format(float(value) * 100.0, digits)


def fmt_bool(value: Optional[bool], yes: str = "yes", no: str = "no",
             na: str = "not established") -> str:
    if value is None:
        return na
    return yes if value else no


def _tsv_cell(value: Any) -> str:
    """One TSV cell: absence as ``NA``, no stray tabs, no ``nan``."""
    if value is None:
        return TSV_NA
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        if value != value:  # NaN
            return TSV_NA
        if value.is_integer() and abs(value) < 1e15:
            return str(int(value))
        text = "{0:.6f}".format(value).rstrip("0").rstrip(".")
        return text or "0"
    if isinstance(value, (list, tuple)):
        return ";".join(_tsv_cell(item) for item in value) if value else ""
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    text = str(value)
    return text.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def drug_order(name: str) -> Tuple[int, Any]:
    """Display order: the catalogue's own drug order first, then alphabetical.

    Unknown drugs sort after the known ones rather than being dropped: a
    catalogue edition that adds a drug must still reach the report.
    """
    canonical = normalise_drug(name)
    try:
        return (0, DRUGS.index(canonical))
    except ValueError:
        return (1, natural_key(canonical))


def catalogue_key(catalogue: str) -> str:
    """Column-safe token for a catalogue name: ``WHO v2`` -> ``who_v2``."""
    return safe_name(catalogue).lower()


# ---------------------------------------------------------------------------
# Metric specifications: a threshold, its source, and how to read it
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MetricSpec:
    """One measurable quantity, the threshold it is judged against, and why.

    The specs exist because three consumers need the same triple — the check
    panel, the TSV, and the coverage/depth strip chart, which draws a tick at
    every bound. Deriving the chart's ticks from the same object that produced
    the check is what stops a figure from disagreeing with the table beside it.
    """

    key: str
    label: str
    #: ``x`` for a depth, ``fraction`` for a 0-1 proportion. The strip chart uses
    #: this to decide which axis a metric belongs on.
    unit: str
    source: str
    minimum: Optional[float] = None
    warn_minimum: Optional[float] = None
    maximum: Optional[float] = None
    warn_maximum: Optional[float] = None
    note: str = ""
    higher_is_better: bool = True
    #: Filled in when the metric cannot exist for this input, e.g. depth on an
    #: assembly. Rendered as the reason it was not measured.
    absent_why: str = ""

    def bounds(self) -> List[Tuple[float, str]]:
        """Threshold ticks for the chart, as (value, label) pairs."""
        out: List[Tuple[float, str]] = []
        for value, label in ((self.minimum, "fail below"), (self.warn_minimum, "target"),
                             (self.warn_maximum, "warn above"), (self.maximum, "fail above")):
            if value is not None:
                out.append((float(value), label))
        return out


QC_METRIC_SPECS: Tuple[MetricSpec, ...] = (
    MetricSpec("mean_depth", "mean depth", "x", source_for("min_depth"),
               minimum=DEGRADED_DEPTH_FLOOR, warn_minimum=MIN_DEPTH,
               note="below {0}x the sample is not callable; {1}x is the target".format(
                   DEGRADED_DEPTH_FLOOR, MIN_DEPTH),
               absent_why="assembly input carries no read depth"),
    MetricSpec("median_depth", "median depth", "x", source_for("min_depth"),
               minimum=DEGRADED_DEPTH_FLOOR, warn_minimum=MIN_DEPTH,
               absent_why="assembly input carries no read depth"),
    MetricSpec("breadth_1x", "breadth at 1x", "fraction", source_for("min_breadth"),
               warn_minimum=MIN_BREADTH,
               absent_why="assembly input carries no per-base coverage"),
    MetricSpec("breadth_10x", "breadth at 10x", "fraction", source_for("min_breadth"),
               warn_minimum=MIN_BREADTH,
               absent_why="assembly input carries no per-base coverage"),
    MetricSpec("breadth_min_depth", "breadth at the depth floor", "fraction",
               source_for("min_breadth"), minimum=MIN_BREADTH,
               note="fraction of the reference reaching the degraded depth floor; "
                    "genome-wide statements are not made below this",
               absent_why="assembly input carries no per-base coverage"),
    MetricSpec("coverage_evenness", "coverage evenness", "fraction",
               source_for("min_coverage_evenness"), warn_minimum=MIN_COVERAGE_EVENNESS,
               note=EVENNESS_DEFINITION,
               absent_why="assembly input carries no per-base coverage"),
    MetricSpec("mapped_fraction", "mapped read fraction", "fraction",
               source_for("min_mapped_fraction"), warn_minimum=MIN_MAPPED_FRACTION,
               note="a flag on library purity, not a species claim",
               absent_why="assembly input has no reads to map"),
    MetricSpec("unambiguous_fraction", "unambiguous base fraction", "fraction",
               source_for("min_unambiguous_fraction"), warn_minimum=MIN_UNAMBIGUOUS_FRACTION,
               note="MTBseq's de-facto heterozygosity filter, reported here rather "
                    "than silently applied",
               absent_why="assembly input carries no allele fractions"),
)

#: Measured, but with no registered threshold. They are reported as observations
#: rather than as checks, because a metric with no published bound cannot pass or
#: fail and printing a green tick beside one would be an invention.
QC_OBSERVATION_SPECS: Tuple[Tuple[str, str, str], ...] = (
    ("total_reads", "total reads", "count"),
    ("mapped_reads", "mapped reads", "count"),
    ("duplicate_fraction", "duplicate fraction", "fraction"),
    ("mean_read_length", "mean read length", "bp"),
    ("mean_base_quality", "mean base quality", "phred"),
    ("gc_content", "GC content", "fraction"),
    ("reference_length", "reference length", "bp"),
)

CONTAMINATION_METRIC_SPECS: Tuple[MetricSpec, ...] = (
    MetricSpec("het_snp_fraction", "genome-wide heterozygous-SNP fraction", "fraction",
               source_for("het_snp_fraction_warn"),
               warn_maximum=HET_SNP_FRACTION_WARN, maximum=HET_SNP_FRACTION_MIXED,
               higher_is_better=False,
               note="two tiers rather than one cutoff: above {0} a mixture is "
                    "possible, above {1} it is called".format(
                        HET_SNP_FRACTION_WARN, HET_SNP_FRACTION_MIXED),
               absent_why="assembly input carries no allele fractions, so "
                          "heterozygosity cannot be measured"),
    MetricSpec("f2", "F2 (lineage-defining minor alleles)", "fraction",
               source_for("f2_mixture_threshold"),
               warn_maximum=F2_MIXTURE_THRESHOLD, higher_is_better=False,
               note="citation unverified on this machine; the value is Mjolnir policy "
                    "pending the primary reference",
               absent_why="assembly input carries no allele fractions"),
    MetricSpec("f47", "F47 (lineage-defining minor alleles)", "fraction",
               source_for("f47_mixture_threshold"),
               warn_maximum=F47_MIXTURE_THRESHOLD, higher_is_better=False,
               note="citation unverified on this machine; the value is Mjolnir policy "
                    "pending the primary reference",
               absent_why="assembly input carries no allele fractions"),
    MetricSpec("unambiguous_fraction", "unambiguous base fraction", "fraction",
               source_for("min_unambiguous_fraction"), warn_minimum=MIN_UNAMBIGUOUS_FRACTION,
               absent_why="assembly input carries no allele fractions"),
    MetricSpec("mapped_fraction", "mapped read fraction", "fraction",
               source_for("min_mapped_fraction"), warn_minimum=MIN_MAPPED_FRACTION,
               absent_why="assembly input has no reads to map"),
    MetricSpec("coverage_breadth", "coverage breadth", "fraction",
               source_for("min_breadth"), warn_minimum=MIN_BREADTH,
               absent_why="assembly input carries no per-base coverage"),
    MetricSpec("coverage_evenness", "coverage evenness", "fraction",
               source_for("min_coverage_evenness"), warn_minimum=MIN_COVERAGE_EVENNESS,
               note=EVENNESS_DEFINITION,
               absent_why="assembly input carries no per-base coverage"),
)

CONTAMINATION_OBSERVATION_SPECS: Tuple[Tuple[str, str, str], ...] = (
    ("het_snp_count", "heterozygous sites", "count"),
    ("snp_sites_examined", "SNP sites examined", "count"),
    ("lineage_het_sites", "heterozygous lineage-defining sites", "count"),
    ("lineage_sites_examined", "lineage-defining sites examined", "count"),
    ("non_target_fraction", "non-target read fraction", "fraction"),
    ("gc_content", "GC content", "fraction"),
)


def find_check(checks: Sequence[Check], key: str, label: str) -> Optional[Check]:
    """An engine-supplied check for this metric, matched on name.

    Matching is deliberately narrow — exact on the field name or on the printed
    label — so that a loosely named check from some other panel cannot be
    mistaken for this metric's verdict.
    """
    wanted = {key.lower(), label.lower(), key.replace("_", " ").lower()}
    for check in checks:
        if str(check.name).strip().lower() in wanted:
            return check
    return None


def _panel(values: Any, specs: Sequence[MetricSpec], existing: Sequence[Check],
           category: str, absent_capability: bool) -> List[Check]:
    """Build a complete check panel, preferring checks an engine already made.

    *absent_capability* is True when the input type cannot produce these metrics
    at all — an assembly has no depth — in which case the "not measured" reading
    names the capability loss rather than implying something went wrong.
    """
    panel: List[Check] = []
    for spec in specs:
        found = find_check(existing, spec.key, spec.label)
        if found is not None:
            panel.append(found)
            continue
        value = getattr(values, spec.key, None)
        why = spec.absent_why if (absent_capability and spec.absent_why) else (
            "{0} was not measured".format(spec.label))
        panel.append(Check.numeric(
            spec.label, value,
            minimum=spec.minimum, maximum=spec.maximum,
            warn_minimum=spec.warn_minimum, warn_maximum=spec.warn_maximum,
            source=spec.source, unit=spec.unit, category=category,
            reading=spec.note, not_measured_why=why))
    return panel


def qc_panel(result: SampleResult) -> List[Check]:
    """Coverage and mapping metrics as checks, each carrying its source.

    The GC check is built separately because its expectation is a property of
    the reference, not a universal constant: H37Rv's composition is registered,
    every other reference's is not, and inventing a band for an NTM reference
    would be a number with no source behind it.
    """
    absent = result.platform == "fasta"
    panel = _panel(result.qc, QC_METRIC_SPECS, result.qc.checks, "qc", absent)
    panel.append(_gc_check(result))
    return panel


def _gc_check(result: SampleResult) -> Check:
    reference = "{0} {1}".format(result.reference or "", result.qc.reference or "").lower()
    is_h37rv = H37RV_ACCESSION.lower() in reference or "h37rv" in reference
    if not is_h37rv:
        value = result.qc.gc_content
        if value is None:
            return Check.not_measured(
                "GC content", "GC content was not measured", source=source_for("h37rv_gc"),
                category="qc")
        return Check(
            name="GC content", value=value, threshold=None,
            source="{0}; no expected GC is registered for this reference, so the "
                   "value is reported without a threshold".format(SRC_POLICY),
            status=STATUS_PASS, measured=True, category="qc", unit="fraction",
            reading="reported without a threshold: an expected GC content is only "
                    "registered for {0}".format(H37RV_ACCESSION))
    return Check.numeric(
        "GC content", result.qc.gc_content,
        warn_minimum=H37RV_GC - GC_TOLERANCE, warn_maximum=H37RV_GC + GC_TOLERANCE,
        source="{0}; tolerance {1} ({2})".format(
            source_for("h37rv_gc"), GC_TOLERANCE, source_for("gc_tolerance")),
        unit="fraction", category="qc",
        reading="a composition warning, not a species call",
        not_measured_why="GC content was not measured")


def contamination_panel(result: SampleResult) -> List[Check]:
    """The contamination metrics as checks, plus the two structural refusals.

    The screen-informative check is the one that matters most here. A Kraken2 run
    against a standard index has a measured sensitivity of 0.0731 for
    *M. tuberculosis* reads, so a clean-looking screen from one of those is not
    evidence of anything, and this panel says so rather than omitting the row.
    """
    contamination = result.contamination
    absent = result.platform == "fasta"
    panel = _panel(contamination, CONTAMINATION_METRIC_SPECS, contamination.checks,
                   "contamination", absent)
    panel.append(Check.boolean(
        "contamination screen informative", contamination.screen_informative or None,
        expected=True, source=source_for("kraken2_mtb_sensitivity_standard"),
        category="contamination", fail_status=STATUS_WARN,
        reading=contamination.screen_note or
        "no read-composition screen capable of resolving mycobacteria was run"))
    verdict_status = {
        VALIDITY_VALID: STATUS_PASS,
        VALIDITY_SUSPECT: STATUS_WARN,
        VALIDITY_INVALID: STATUS_FAIL,
        VALIDITY_NOT_ASSESSED: STATUS_WARN,
    }.get(contamination.verdict, STATUS_WARN)
    panel.append(Check(
        name="sample validity verdict", value=contamination.verdict,
        threshold=None, source=source_for("contamination_evidence"),
        status=verdict_status, category="contamination",
        measured=contamination.verdict != VALIDITY_NOT_ASSESSED,
        reading=contamination.verdict_reason or CONTAMINATION_EVIDENCE))
    return panel


def all_checks(result: SampleResult) -> List[Check]:
    """Every check the report shows, in the order it shows them."""
    return list(result.checks) + qc_panel(result) + contamination_panel(result)


# ---------------------------------------------------------------------------
# Label/value panels shared by the PDF and the HTML
# ---------------------------------------------------------------------------

def identity_pairs(result: SampleResult) -> List[Tuple[str, str]]:
    """The sample's identity block: what was run, on what, with what."""
    inputs = list(result.inputs) or ["not recorded"]
    pairs = [
        ("Sample", result.sample_id),
        ("Platform", result.platform),
        ("Input files", "; ".join(str(p) for p in inputs)),
        ("Reference", result.reference or "not recorded"),
        ("Profile", result.profile),
        ("Mjolnir version", result.mjolnir_version or __version__),
    ]
    if result.started_at or result.finished_at:
        pairs.append(("Run", "{0} to {1}".format(result.started_at or "?",
                                                 result.finished_at or "?")))
    if result.runtime_seconds:
        pairs.append(("Runtime", "{0:.1f} s".format(result.runtime_seconds)))
    return pairs


def species_pairs(result: SampleResult) -> List[Tuple[str, str]]:
    """Species, and how far down the tree the evidence actually reaches."""
    species = result.species
    pairs: List[Tuple[str, str]] = [("Species", species.display)]
    if not species.resolved_to_species:
        detail = "not resolved to species"
        if species.complex.upper() == COMPLEX_MTBC.upper() or species.is_mtbc:
            detail = "not resolved to species - " + MTBC_UNRESOLVED_TEXT
        pairs.append(("Resolution", detail))
    else:
        pairs.append(("Resolution", "resolved to species"))
    pairs.append(("Complex", species.complex or "not assigned"))
    pairs.append(("Confidence", species.confidence))
    pairs.append(("Method", species.method or "not recorded"))
    if species.ani is None:
        pairs.append(("ANI", "not measured (species floor {0}%, {1})".format(
            ANI_SPECIES_FLOOR, source_for("ani_species_floor"))))
    else:
        pairs.append(("ANI", "{0:.2f}% against {1} (species floor {2}%)".format(
            species.ani, species.reference or "the best-matching reference",
            ANI_SPECIES_FLOOR)))
    pairs.append(("Aligned fraction", "{0} (a high ANI over a small aligned "
                                      "fraction is not a species call; floor {1})".format(
        fmt_fraction(species.aligned_fraction), ANI_MIN_ALIGNED_FRACTION)))
    if species.subspecies:
        pairs.append(("Subspecies", species.subspecies))
    if species.candidates:
        runners = "; ".join(
            "{0} {1}".format(c.get("name", "?"), fmt_number(c.get("ani"), 2, "%", na="ANI NA"))
            for c in species.candidates[:3])
        pairs.append(("Runners-up", runners))
    for caveat in species.caveats:
        pairs.append(("Caveat", caveat))
    return pairs


def lineage_pairs(result: SampleResult) -> List[Tuple[str, str]]:
    """Lineage with its barcode evidence, and the flags a clinician must see."""
    lineage = result.lineage
    pairs: List[Tuple[str, str]] = [
        ("Lineage", lineage.lineage or "not determined"),
        ("Sublineage", lineage.sublineage or "not determined"),
        ("Confidence", lineage.confidence),
        ("Scheme", lineage.scheme or "not recorded"),
        ("Method", lineage.method or "not recorded"),
    ]
    support = (
        "{0} of {1} callable defining sites carry the derived allele ({2}); "
        "{1} of {3} defining sites were callable ({4}). Support floor {5}, "
        "callable floor {6} ({7})".format(
            lineage.barcode_sites_supporting, lineage.barcode_sites_callable,
            fmt_fraction(lineage.support_fraction), lineage.barcode_sites_total,
            fmt_fraction(lineage.callable_fraction),
            MIN_BARCODE_SUPPORT_FRACTION, MIN_BARCODE_CALLABLE_FRACTION,
            source_for("min_barcode_support_fraction")))
    pairs.append(("Barcode evidence", support))
    if lineage.is_bcg:
        pairs.append(("BCG", "yes - " + BCG_PZA_NOTE))
    if lineage.is_animal or lineage.animal_variant:
        detail = lineage.animal_variant or "animal-adapted lineage"
        if "bovis" in detail.lower() or "bcg" in detail.lower():
            detail = detail + " - " + MBOVIS_CAVEAT
        pairs.append(("Animal lineage", detail))
    if lineage.mixed_lineages:
        pairs.append(("Mixed lineages", "{0} - more than one lineage's defining SNPs "
                                        "were supported, which is a mixed-infection "
                                        "signal and is reported as one".format(
                          ", ".join(lineage.mixed_lineages))))
    for caveat in lineage.caveats:
        pairs.append(("Caveat", caveat))
    return pairs


def validity_pairs(result: SampleResult) -> List[Tuple[str, str]]:
    """The sample-validity verdict and what stands behind it."""
    contamination = result.contamination
    pairs = [
        ("Sample validity", contamination.verdict),
        ("Reason", contamination.verdict_reason or "no reason recorded"),
        ("Mixture class", contamination.mixture_class),
        ("Screen", "{0}{1}".format(
            contamination.screen_method or "no read-composition screen run",
            "" if contamination.screen_informative else " - not informative for "
            "mycobacterial purity")),
    ]
    if contamination.screen_note:
        pairs.append(("Screen note", contamination.screen_note))
    if contamination.non_target_fraction is not None:
        pairs.append(("Non-target reads", "{0}, resolved to {1}".format(
            fmt_fraction(contamination.non_target_fraction),
            contamination.non_target_resolution or "an unstated rank")))
    pairs.append(("Why a verdict and not a percentage", CONTAMINATION_EVIDENCE))
    for caveat in contamination.caveats:
        pairs.append(("Caveat", caveat))
    return pairs


def headline_sentence(result: SampleResult) -> Tuple[str, str]:
    """The one sentence page 1 leads with, and where it came from.

    The model writes prose over rule-derived verdicts, and when there is no
    model — or when the model's answer was discarded by the discipline rules —
    the rule-derived sentence stands in its place and the provenance says which
    of the two happened. It never silently loses the model and never silently
    keeps a bad answer.
    """
    interpretation = result.interpretation
    rule = rule_headline(result)
    if interpretation is None:
        return rule, "rule-derived (no interpretation was produced)"
    if interpretation.rule_only or not interpretation.headline.strip():
        if interpretation.discarded_reason:
            return rule, "rule-derived (the model's answer was discarded: {0})".format(
                interpretation.discarded_reason)
        return rule, "rule-derived (no model host was reachable)"
    origin = "model {0}".format(interpretation.model or "unnamed")
    if interpretation.host:
        origin += " at {0}".format(interpretation.host)
    if interpretation.playbook:
        origin += ", playbook {0}".format(interpretation.playbook)
    return interpretation.headline.strip(), origin


def rule_headline(result: SampleResult) -> str:
    """The rule-derived headline, built only from calls that were actually made.

    Written so that the no-findings branch cannot be read as reassurance: it
    states the number of drugs evaluated and repeats that absence of a
    catalogued mutation is not susceptibility.
    """
    parts: List[str] = []
    verdict = result.contamination.verdict
    if verdict == VALIDITY_INVALID:
        parts.append("This sample is not valid for reporting ({0}).".format(
            result.contamination.verdict_reason or "see the contamination panel"))
    elif verdict == VALIDITY_SUSPECT:
        parts.append("Sample validity is suspect ({0}); read the calls below with "
                     "that in mind.".format(
                         result.contamination.verdict_reason or "see the contamination panel"))

    resistant = sorted(result.resistant_drugs(), key=lambda d: drug_order(d.drug))
    if resistant:
        named = ", ".join("{0} [{1}]".format(d.drug, CALL_GLYPH.get(d.call, d.call))
                          for d in resistant)
        parts.append("Resistance determinants were detected for {0}.".format(named))
        outside = [d.drug for d in resistant if d.call == CALL_R_OUTSIDE_WHO]
        if outside:
            parts.append("{0} rests on a catalogue other than WHO at a variant WHO "
                         "does not grade, and is not equivalent to a WHO Group 1 "
                         "call.".format(", ".join(outside)))
    elif result.drugs:
        parts.append("{0} was reported for all {1} drugs evaluated; no catalogued "
                     "mutation was found, which is not evidence of "
                     "susceptibility.".format(NO_DETERMINANT_TEXT.capitalize(),
                                              len(result.drugs)))
    else:
        parts.append("No drug was evaluated for this sample.")

    disagreeing = [d.drug for d in result.disagreements()]
    if disagreeing:
        parts.append("Catalogues disagree for {0}; the annex shows all three side "
                     "by side.".format(", ".join(sorted(disagreeing, key=natural_key))))
    unmeasured = result.unmeasured()
    if unmeasured:
        parts.append("{0} could not be measured and {1} reported as absent, not "
                     "as normal.".format(plural(len(unmeasured), "metric"),
                                         "is" if len(unmeasured) == 1 else "are"))
    return " ".join(parts)


def cohort_headline(cohort: CohortResult) -> Tuple[str, str]:
    """The cohort's leading sentence, and its provenance."""
    interpretation = cohort.interpretation
    clustered = sum(len(c.members) for c in cohort.clusters)
    thin = [p for p in cohort.pairs
            if p.shared_callable_sites is not None
            and p.shared_callable_sites < MIN_SHARED_CALLABLE_SITES]
    rule = (
        "{0} samples compared at a threshold of {1} SNPs ({2}); {3} "
        "covering {4}.".format(
            len(cohort.samples), cohort.threshold,
            cohort.threshold_basis or "no basis recorded",
            plural(len(cohort.clusters), "cluster"),
            plural(clustered, "sample")))
    uncompared = _uncompared_pairs(cohort)
    if uncompared:
        # Written out in both numbers rather than assembled from fragments: the
        # verb, the pronoun and the noun all have to agree, and gluing them
        # together is how "1 pair was ... reported as absent distances" happens.
        rule += (" 1 pair was never compared, and is reported as an absent "
                 "distance, not as zero." if uncompared == 1 else
                 " {0} pairs were never compared, and are reported as absent "
                 "distances, not as zero.".format(uncompared))
    if thin:
        rule += (" {0} less than {1:,} callable bases, below which a distance "
                 "is not comparable to the published SNP thresholds.".format(
                     plural(len(thin), "pair shares", "pairs share"),
                     MIN_SHARED_CALLABLE_SITES))
    if interpretation is None or interpretation.rule_only or not interpretation.headline.strip():
        reason = ""
        if interpretation is not None and interpretation.discarded_reason:
            reason = " (the model's answer was discarded: {0})".format(
                interpretation.discarded_reason)
        return rule, "rule-derived" + (reason or "")
    return interpretation.headline.strip(), "model {0}".format(
        interpretation.model or "unnamed")


def _uncompared_pairs(cohort: CohortResult) -> int:
    samples = list(cohort.samples)
    total = len(samples) * (len(samples) - 1) // 2
    measured = len([p for p in cohort.pairs if p.snps is not None])
    return max(0, total - measured)


# ---------------------------------------------------------------------------
# Drug rows and the drug grid
# ---------------------------------------------------------------------------

def _catalogue_view(drug_call: Any, catalogue: str) -> Tuple[str, str, str]:
    """One catalogue's (call, grade, evidence) for a drug.

    Several graded variants can contribute to one drug in one catalogue, so the
    call is the most severe of them — the same reduction :func:`worst_call`
    performs everywhere else — and the grades are kept as a list so the annex can
    show that two variants were involved rather than one.
    """
    calls = drug_call.calls_by_catalogue(catalogue)
    if not calls:
        return CALL_NO_CALL, "", ""
    call = worst_call([c.call for c in calls])
    grades = "; ".join(sorted({c.grade for c in calls if c.grade}, key=natural_key))
    evidence = "; ".join(sorted({c.evidence for c in calls if c.evidence}, key=natural_key))
    return call, grades, evidence


def drug_flags(drug_call: Any) -> List[str]:
    """The marks page 1 puts beside a drug."""
    flags: List[str] = []
    if drug_call.disagreement:
        flags.append("{0} catalogues disagree{1}".format(
            FLAG_DISAGREEMENT,
            " ({0})".format(drug_call.disagreement_kind) if drug_call.disagreement_kind else ""))
    if drug_call.suppressed_by:
        flags.append("{0} suppressed: {1}".format(FLAG_SUPPRESSED, drug_call.suppressed_by))
    if drug_call.target_covered is False:
        flags.append("{0} target regions not callable".format(FLAG_NOT_EVALUABLE))
    if drug_call.caveats:
        flags.append("{0} {1}".format(FLAG_CAVEAT,
                                      plural(len(drug_call.caveats), "platform caveat")))
    return flags


def drug_rows(result: SampleResult) -> List[Dict[str, Any]]:
    """One row per drug, with every catalogue's own answer beside the consensus."""
    rows: List[Dict[str, Any]] = []
    for call in sorted(result.drugs, key=lambda d: drug_order(d.drug)):
        row: Dict[str, Any] = {
            "sample": result.sample_id,
            "drug": call.drug,
            "drug_code": drug_code(call.drug),
            "call": call.call,
            "call_glyph": CALL_GLYPH.get(call.call, call.call),
            "call_label": call.label,
            "confidence": call.confidence,
            "who_graded": call.who_graded,
            "who_grade": call.who_grade or None,
            "disagreement": call.disagreement,
            "disagreement_kind": call.disagreement_kind or None,
            "suppressed_by": call.suppressed_by or None,
            "level": call.level or None,
            "cross_resistance": list(call.cross_resistance),
            "target_covered": call.target_covered,
            "supporting_variants": list(call.supporting_variants),
            "determinants_display": determinants_display(call),
            "caveats": list(call.caveats),
            "note": call.note or None,
        }
        for catalogue in CATALOGUES:
            key = catalogue_key(catalogue)
            cat_call, grade, evidence = _catalogue_view(call, catalogue)
            row[key + "_call"] = cat_call
            row[key + "_grade"] = grade or None
            row[key + "_evidence"] = evidence or None
        rows.append(row)
    return rows


#: How many determinants the front-page table prints for one drug before it
#: summarises. SOURCE: Mjolnir policy, forced by a real run — an *M. chimaera*
#: isolate against H37Rv matched enough graded variants per drug to build a
#: table row 4,076 points tall, which is taller than the page, and reportlab
#: raises LayoutError rather than truncating. So the report did not merely look
#: cluttered, it failed to render at all. The full list is always in
#: ``<sample>.catalogue_calls.tsv`` and in the annex.
MAX_DETERMINANTS_SHOWN = 4


def determinants_display(call: DrugCall) -> str:
    """The determinants worth naming beside a drug on page one.

    The variants that *drive* the call come first — those whose own catalogue
    call matches the drug's final call — because for a resistant drug the
    clinically relevant line is the Group 1 mutation, not the eleven benign
    variants that happened to sit in a graded gene beside it. On a *M. bovis*
    isolate that is ``pncA_p.His57Asp`` rather than a paragraph of PPE35.

    Never silently truncated: what is not shown is counted, so the reader knows
    to look in the annex rather than believing they have seen everything.
    """
    keys = list(call.supporting_variants)
    if not keys:
        return ""
    driving = [key for key, entries in
               ((k, [c for c in call.catalogue_calls if c.variant_key == k]) for k in keys)
               if any(c.call == call.call for c in entries)]
    ordered = driving + [key for key in keys if key not in driving]
    shown = ordered[:MAX_DETERMINANTS_SHOWN]
    hidden = len(ordered) - len(shown)
    text = "; ".join(shown)
    if hidden:
        text += "; +{0} more (annex)".format(hidden)
    return text


def catalogue_call_rows(result: SampleResult) -> List[Dict[str, Any]]:
    """One row per (drug, catalogue, variant): the complete catalogue evidence."""
    rows: List[Dict[str, Any]] = []
    for call in sorted(result.drugs, key=lambda d: drug_order(d.drug)):
        for entry in call.catalogue_calls:
            rows.append({
                "sample": result.sample_id,
                "drug": call.drug,
                "mjolnir_call": call.call,
                "catalogue": entry.catalogue,
                "catalogue_call": entry.call,
                "grade": entry.grade or None,
                "variant": entry.variant_key or None,
                "matched_by": entry.matched_by or None,
                "evidence": entry.evidence or None,
                "comment": entry.comment or None,
                "catalogue_version": entry.catalogue_version or None,
                "catalogue_checksum": entry.catalogue_checksum or None,
                "source": entry.source or None,
            })
    return rows


def disagreement_rows(result: SampleResult) -> List[Dict[str, Any]]:
    """Every catalogue's position on every drug where they conflict.

    One row per drug per catalogue, including the catalogues that said nothing,
    because "MTBseq is silent here" is part of the disagreement and omitting the
    row would hide it.
    """
    rows: List[Dict[str, Any]] = []
    for call in sorted(result.disagreements(), key=lambda d: drug_order(d.drug)):
        for catalogue in CATALOGUES:
            cat_call, grade, evidence = _catalogue_view(call, catalogue)
            note = ""
            if catalogue == CATALOGUE_MTBSEQ:
                note = MTBSEQ_ASYMMETRY_NOTE
            elif catalogue == CATALOGUE_WHO and not call.who_graded:
                note = ("WHO does not grade any of the supporting variants; a "
                        "resistance call from another catalogue here is reported as "
                        "{0} and is not equivalent to a WHO Group 1 call".format(
                            CALL_R_OUTSIDE_WHO))
            rows.append({
                "sample": result.sample_id,
                "drug": call.drug,
                "mjolnir_call": call.call,
                "disagreement_kind": call.disagreement_kind or None,
                "catalogue": catalogue,
                "catalogue_call": cat_call,
                "grade": grade or None,
                "evidence": evidence or None,
                "supporting_variants": list(call.supporting_variants),
                "note": note or None,
            })
    return rows


@dataclass
class GridCell:
    """One cell of a call grid: the call, its glyph, and what it stands on."""

    call: str = CALL_NO_CALL
    glyph: str = CALL_GLYPH[CALL_NO_CALL]
    detail: str = ""
    #: True when the row's target regions were not callable at all, which is a
    #: different statement from a no-call and is drawn differently.
    not_evaluable: bool = False


@dataclass
class Grid:
    """A labelled matrix of :class:`GridCell`, ready for a heatmap."""

    row_labels: List[str] = field(default_factory=list)
    row_sublabels: List[str] = field(default_factory=list)
    column_labels: List[str] = field(default_factory=list)
    cells: List[List[GridCell]] = field(default_factory=list)
    #: Flag strings per row, drawn in the margin of the heatmap.
    row_flags: List[List[str]] = field(default_factory=list)
    caption: str = ""

    @property
    def is_empty(self) -> bool:
        return not self.row_labels or not self.column_labels


def drug_grid(result: SampleResult) -> Grid:
    """Drugs by catalogue, with Mjolnir's consensus as the final column.

    This is the front page's heatmap. Putting the three catalogues beside the
    consensus, rather than showing only the consensus, is what makes a
    disagreement visible at a glance instead of buried in an annex.
    """
    grid = Grid(column_labels=list(GRID_COLUMNS),
                caption="Each catalogue's own answer, then Mjolnir's consensus. "
                        "ND marks the absence of a catalogued determinant, which is "
                        "not susceptibility.")
    for call in sorted(result.drugs, key=lambda d: drug_order(d.drug)):
        grid.row_labels.append(drug_code(call.drug))
        grid.row_sublabels.append(call.drug)
        grid.row_flags.append(drug_flags(call))
        row: List[GridCell] = []
        for catalogue in CATALOGUES:
            cat_call, grade, _evidence = _catalogue_view(call, catalogue)
            detail = "{0}: {1}".format(catalogue, call_label(cat_call))
            if grade:
                detail += " [{0}]".format(grade)
            row.append(GridCell(call=cat_call, glyph=CALL_GLYPH.get(cat_call, cat_call),
                                detail=detail))
        row.append(GridCell(
            call=call.call, glyph=CALL_GLYPH.get(call.call, call.call),
            detail="{0}: {1}".format(CONSENSUS_COLUMN, call.label),
            not_evaluable=call.target_covered is False))
        grid.cells.append(row)
    return grid


def cohort_drug_grid(results: Sequence[SampleResult]) -> Grid:
    """Samples by drug: the cohort view of the same consensus calls."""
    drugs: List[str] = []
    for result in results:
        for call in result.drugs:
            if call.drug not in drugs:
                drugs.append(call.drug)
    drugs.sort(key=drug_order)
    grid = Grid(column_labels=[drug_code(d) for d in drugs],
                caption="Consensus call per sample and drug. ND is the absence of a "
                        "catalogued determinant, not susceptibility.")
    for result in results:
        grid.row_labels.append(result.sample_id)
        grid.row_sublabels.append(result.species.display)
        flags = []
        if result.disagreements():
            flags.append("{0} {1} with catalogue disagreement".format(
                FLAG_DISAGREEMENT, plural(len(result.disagreements()), "drug")))
        grid.row_flags.append(flags)
        row: List[GridCell] = []
        for drug in drugs:
            call = result.drug(drug)
            if call is None:
                row.append(GridCell(call=CALL_NO_CALL, glyph=CALL_GLYPH[CALL_NO_CALL],
                                    detail="{0}: drug not evaluated in this "
                                           "sample".format(drug)))
                continue
            row.append(GridCell(call=call.call, glyph=CALL_GLYPH.get(call.call, call.call),
                                detail="{0}: {1}".format(drug, call.label),
                                not_evaluable=call.target_covered is False))
        grid.cells.append(row)
    return grid


# ---------------------------------------------------------------------------
# Variants
# ---------------------------------------------------------------------------

def variant_rows(result: SampleResult) -> List[Dict[str, Any]]:
    """Every observed variant with its coordinate, evidence and per-catalogue grade."""
    rows: List[Dict[str, Any]] = []
    for variant in sorted(result.variants, key=lambda v: (v.chrom, v.pos, v.alt)):
        row: Dict[str, Any] = {
            "sample": result.sample_id,
            "coordinate": "{0}:{1}{2}>{3}".format(variant.chrom, variant.pos,
                                                  variant.ref.upper(), variant.alt.upper()),
            "chrom": variant.chrom,
            "pos": variant.pos,
            "ref": variant.ref,
            "alt": variant.alt,
            "gene": variant.gene or None,
            "hgvs": variant.hgvs or None,
            "hgvs_key": variant.hgvs_key or None,
            "hgvs_alias": variant.hgvs_alias or None,
            "locus_tag": variant.locus_tag or None,
            "variant_type": variant.variant_type,
            "effect": variant.effect or None,
            "depth": variant.depth,
            "alt_reads": variant.alt_reads,
            "ref_reads": variant.ref_reads,
            "allele_fraction": round_or_none(variant.allele_fraction, 4),
            "is_major": variant.is_major,
            "qual": round_or_none(variant.qual, 2),
            "masked": variant.masked,
            "filters": list(variant.filters),
            "source_caller": variant.source_caller or None,
            "note": variant.note or None,
        }
        drugs = sorted({c.drug for c in variant.catalogue_calls}, key=drug_order)
        row["drugs"] = drugs
        row["worst_call"] = worst_call([c.call for c in variant.catalogue_calls]) \
            if variant.catalogue_calls else CALL_NO_CALL
        for catalogue in CATALOGUES:
            key = catalogue_key(catalogue)
            entries = [c for c in variant.catalogue_calls if c.catalogue == catalogue]
            row[key + "_grade"] = "; ".join(
                "{0}={1}".format(c.drug, c.grade or c.call) for c in entries) or None
        rows.append(row)
    return rows


def variant_catalogue_rows(result: SampleResult) -> List[Dict[str, Any]]:
    """One row per (variant, catalogue, drug) — the fully expanded grading table."""
    rows: List[Dict[str, Any]] = []
    for variant in sorted(result.variants, key=lambda v: (v.chrom, v.pos, v.alt)):
        for entry in variant.catalogue_calls:
            rows.append({
                "sample": result.sample_id,
                "variant": variant.display,
                "coordinate": "{0}:{1}{2}>{3}".format(variant.chrom, variant.pos,
                                                      variant.ref.upper(),
                                                      variant.alt.upper()),
                "gene": variant.gene or None,
                "hgvs": variant.hgvs or None,
                "allele_fraction": round_or_none(variant.allele_fraction, 4),
                "depth": variant.depth,
                "catalogue": entry.catalogue,
                "drug": entry.drug,
                "grade": entry.grade or None,
                "call": entry.call,
                "matched_by": entry.matched_by or None,
                "catalogue_version": entry.catalogue_version or None,
                "catalogue_checksum": entry.catalogue_checksum or None,
                "evidence": entry.evidence or None,
                "comment": entry.comment or None,
            })
    return rows


def catalogue_variant_points(result: SampleResult) -> List[Dict[str, Any]]:
    """Variants for the allele-fraction plot, with the call that colours them.

    Variants whose allele fraction is ``None`` are kept rather than filtered out.
    On assembly input every one of them is None, and a plot that quietly dropped
    them would look like a sample with no variants at all.
    """
    # A variant is coloured by the call it actually produced, not by the raw
    # catalogue grade behind it. The distinction matters for exactly the case the
    # design refuses to blur: a variant WHO does not grade that tbdb calls R
    # yields R-outside-WHO, and the plot must not draw it the same red as a WHO
    # Group 1 determinant.
    consensus: Dict[str, str] = {}
    for drug_call in result.drugs:
        for key in drug_call.supporting_variants:
            previous = consensus.get(key)
            consensus[key] = (worst_call([previous, drug_call.call]) if previous
                              else drug_call.call)
    points: List[Dict[str, Any]] = []
    for variant in result.variants:
        catalogue_call = worst_call([c.call for c in variant.catalogue_calls]) \
            if variant.catalogue_calls else CALL_NO_CALL
        call = consensus.get(variant.hgvs_key) or catalogue_call
        points.append({
            "pos": variant.pos,
            "allele_fraction": variant.allele_fraction,
            "depth": variant.depth,
            "call": call,
            "catalogue_call": catalogue_call,
            "catalogued": bool(variant.catalogue_calls),
            "label": variant.display,
            "gene": variant.gene,
            "masked": variant.masked,
            "drugs": sorted({c.drug for c in variant.catalogue_calls}, key=drug_order),
        })
    points.sort(key=lambda p: (-CALL_SEVERITY.get(p["call"], 0), p["pos"]))
    return points


# ---------------------------------------------------------------------------
# Checks, QC, contamination, methods
# ---------------------------------------------------------------------------

def check_rows(checks: Sequence[Check], sample: str = "") -> List[Dict[str, Any]]:
    """Checks as rows, with the threshold that was applied and who published it."""
    rows: List[Dict[str, Any]] = []
    for check in checks:
        rows.append({
            "sample": sample or None,
            "category": check.category or None,
            "check": check.name,
            "measured": check.measured,
            "value": check.value,
            "comparison": check.comparison or None,
            "threshold": check.threshold,
            "unit": check.unit or None,
            "status": "not measured" if not check.measured else check.status,
            "source": check.source or None,
            "reading": check.reading or None,
        })
    return rows


def qc_rows(result: SampleResult) -> List[Dict[str, Any]]:
    return check_rows(qc_panel(result), result.sample_id)


def contamination_rows(result: SampleResult) -> List[Dict[str, Any]]:
    return check_rows(contamination_panel(result), result.sample_id)


def observation_rows(result: SampleResult) -> List[Dict[str, Any]]:
    """Measured quantities that carry no registered threshold.

    They are shown, because they are evidence, and they are shown without a
    status, because a status implies a bound and Mjolnir does not have one for
    these. Inventing one here would be exactly the bare magic number the first
    house rule forbids.
    """
    rows: List[Dict[str, Any]] = []
    for source_obj, specs, panel in ((result.qc, QC_OBSERVATION_SPECS, "qc"),
                                     (result.contamination, CONTAMINATION_OBSERVATION_SPECS,
                                      "contamination")):
        for key, label, unit in specs:
            value = getattr(source_obj, key, None)
            rows.append({
                "sample": result.sample_id,
                "panel": panel,
                "observation": label,
                "value": value,
                "unit": unit,
                "measured": value is not None,
                "note": None if value is not None else "not measured",
            })
    contamination = result.contamination
    if contamination.non_target_labels:
        for entry in contamination.non_target_labels:
            rows.append({
                "sample": result.sample_id,
                "panel": "contamination",
                "observation": "non-target label: {0}".format(entry.get("label", "?")),
                "value": entry.get("fraction"),
                "unit": "fraction",
                "measured": entry.get("fraction") is not None,
                "note": entry.get("note") or contamination.non_target_resolution or None,
            })
    return rows


def lineage_support_rows(result: SampleResult) -> List[Dict[str, Any]]:
    """Per-site barcode evidence, as the typing stage recorded it."""
    preferred = ("position", "lineage", "sublineage", "expected", "observed", "allele",
                 "depth", "allele_fraction", "supports", "callable", "note")
    rows: List[Dict[str, Any]] = []
    for entry in result.lineage.support:
        row: Dict[str, Any] = {"sample": result.sample_id}
        for key in preferred:
            if key in entry:
                row[key] = entry[key]
        for key in sorted(entry):
            if key not in row:
                row[key] = entry[key]
        rows.append(row)
    return rows


def database_rows(databases: Sequence[Any]) -> List[Dict[str, Any]]:
    """Database name, version, checksum and licence — the reproducibility block.

    The checksum is not decoration: a catalogue-version mismatch between two
    installations changes calls, and this table is how two labs find out that
    they were not running the same catalogue.
    """
    rows: List[Dict[str, Any]] = []
    for database in databases:
        rows.append({
            "database": database.name,
            "version": database.version,
            "checksum": database.checksum or None,
            "path": database.path or None,
            "licence": database.licence or None,
            "citation": database.citation or None,
            "url": database.url or None,
            "fetched": database.fetched or None,
            "note": database.note or None,
        })
    return rows


def tool_version_rows(versions: Dict[str, str]) -> List[Dict[str, Any]]:
    return [{"tool": name, "version": versions[name]} for name in sorted(versions)]


def threshold_rows() -> List[Dict[str, Any]]:
    """Every registered threshold, its value, its source and whether it was checked."""
    return [{
        "threshold": entry.name,
        "value": entry.value,
        "unit": entry.unit or None,
        "source": entry.source,
        "citation_verified": entry.verified,
        "note": entry.note or None,
    } for entry in all_thresholds()]


def unverified_rows() -> List[Dict[str, Any]]:
    """The thresholds whose citation nobody has checked against the primary document."""
    return [{
        "threshold": entry.name,
        "value": entry.value,
        "source": entry.source,
        "note": entry.note or None,
    } for entry in unverified()]


def caveat_lines(result: SampleResult) -> List[str]:
    """Everything the report must state about this platform and this sample.

    Platform caveats come from the config table rather than from the result, so
    an ONT report states the minor-variant under-detection consequence even if
    the pipeline that produced the result forgot to attach it.
    """
    lines: List[str] = []
    for text in platform_caveats(result.platform):
        if text not in lines:
            lines.append(text)
    for text in list(result.caveats) + list(result.warnings):
        if text and text not in lines:
            lines.append(text)
    return lines


def methods_pairs(result: SampleResult) -> List[Tuple[str, str]]:
    """The methods block: what was run, against what, and under which numbers."""
    pairs = [
        ("Platform", result.platform),
        ("Reference", result.reference or "not recorded"),
        ("Variant major-allele threshold", "{0} ({1})".format(
            MAJOR_VARIANT_FRACTION, source_for("major_variant_fraction"))),
        ("Target depth", "{0}x, degraded floor {1}x ({2})".format(
            MIN_DEPTH, DEGRADED_DEPTH_FLOOR, source_for("min_depth"))),
        ("Species method", result.species.method or "not recorded"),
        ("Lineage scheme", result.lineage.scheme or "not recorded"),
        ("Consensus rule", "WHO v2 is the anchor; where WHO does not grade a variant "
                           "and another catalogue calls resistance the drug is "
                           "reported as {0} ({1})".format(CALL_R_OUTSIDE_WHO, SRC_DESIGN)),
        ("Catalogue asymmetry", MTBSEQ_ASYMMETRY_NOTE),
    ]
    return pairs


# ---------------------------------------------------------------------------
# Cohort
# ---------------------------------------------------------------------------

def distance_rows(cohort: CohortResult) -> List[Dict[str, Any]]:
    """Pairwise distances, each beside the shared callable sequence it rests on."""
    clustered = {}
    for cluster in cohort.clusters:
        for member in cluster.members:
            clustered[member] = cluster.cluster_id
    rows: List[Dict[str, Any]] = []
    for pair in sorted(cohort.pairs, key=lambda p: (natural_key(p.sample_a),
                                                    natural_key(p.sample_b))):
        shared = pair.shared_callable_sites
        below = None
        if pair.snps is not None and cohort.threshold is not None:
            below = pair.snps <= cohort.threshold
        rows.append({
            "sample_a": pair.sample_a,
            "sample_b": pair.sample_b,
            "snps": pair.snps,
            "shared_callable_sites": shared,
            "snps_per_mb": round_or_none(pair.snps_per_mb, 3),
            "masked_sites": pair.masked_sites,
            "within_threshold": below,
            "shared_sites_sufficient": None if shared is None
            else shared >= MIN_SHARED_CALLABLE_SITES,
            "same_cluster": (clustered.get(pair.sample_a) is not None
                             and clustered.get(pair.sample_a) == clustered.get(pair.sample_b)),
            "note": pair.note or None,
        })
    return rows


def cluster_rows(cohort: CohortResult) -> List[Dict[str, Any]]:
    return [{
        "cluster": cluster.cluster_id,
        "size": len(cluster.members),
        "members": sorted(cluster.members, key=natural_key),
        "threshold": cluster.threshold if cluster.threshold is not None else cohort.threshold,
        "threshold_basis": cohort.threshold_basis or None,
        "max_distance": cluster.max_distance,
        "min_shared_callable_sites": cluster.min_shared_callable_sites,
        "note": cluster.note or None,
    } for cluster in cohort.clusters]


def distance_matrix_rows(cohort: CohortResult) -> List[Dict[str, Any]]:
    """The square matrix, with ``NA`` wherever a pair was never compared.

    Never zero for an uncompared pair. Two samples that were not compared are not
    identical, and a matrix that fills its gaps with zeros produces clusters that
    do not exist.
    """
    matrix = cohort.distance_matrix()
    rows: List[Dict[str, Any]] = []
    for sample in cohort.samples:
        row: Dict[str, Any] = {"sample": sample}
        row.update(dict((other, matrix[sample][other]) for other in cohort.samples))
        rows.append(row)
    return rows


def cohort_pairs(cohort: CohortResult) -> List[Tuple[str, str]]:
    """The cohort's methods block: mask, threshold and denominators."""
    return [
        ("Samples", str(len(cohort.samples))),
        ("Reference", cohort.reference or "not recorded"),
        ("Clustering threshold", "{0} SNPs".format(cohort.threshold)
         if cohort.threshold is not None else "not set"),
        ("Threshold basis", cohort.threshold_basis or "no basis recorded"),
        ("Mask", cohort.mask_name or "no mask recorded - masking is mandatory before "
                                     "counting distances, so an unnamed mask is a gap"),
        ("Masked positions", "{0} ({1})".format(
            fmt_number(cohort.masked_sites, na="not recorded"),
            fmt_fraction(cohort.masked_fraction, na="fraction not recorded"))),
        ("Joint sites before masking", fmt_number(cohort.joint_sites, na="not recorded")),
        ("Shared-callable floor", "{0:,} bp ({1})".format(
            MIN_SHARED_CALLABLE_SITES, source_for("min_shared_callable_sites"))),
        ("Pairs never compared", str(_uncompared_pairs(cohort))),
    ]


def upgma_tree(labels: Sequence[str],
               distance: Dict[Tuple[str, str], Optional[int]]
               ) -> Tuple[List[Dict[str, Any]], List[str], List[str]]:
    """Average-linkage tree over the samples whose distances are all present.

    Returns ``(merges, leaves, excluded)``. Leaf node ids are indices into
    *leaves*; a merge node's ``left`` and ``right`` are node ids, and its
    ``height`` is half the linkage distance so that the drawn branch length is
    the usual UPGMA half-distance.

    Samples with any missing comparison are excluded and returned separately
    rather than imputed. Filling a gap with a large number would draw a
    confident-looking branch out of an absence of data, and the caption is
    expected to name the excluded samples.
    """
    usable = []
    excluded = []
    for a in labels:
        complete = True
        for b in labels:
            if a == b:
                continue
            if distance.get((a, b)) is None and distance.get((b, a)) is None:
                complete = False
                break
        (usable if complete else excluded).append(a)

    if len(usable) < 2:
        return [], list(usable), list(excluded)

    def get(a: str, b: str) -> float:
        value = distance.get((a, b))
        if value is None:
            value = distance.get((b, a))
        return float(value if value is not None else 0)

    clusters: Dict[int, Dict[str, Any]] = {}
    for index, name in enumerate(usable):
        clusters[index] = {"id": index, "members": [name], "height": 0.0, "leaf": name}
    active = list(clusters)
    current: Dict[Tuple[int, int], float] = {}
    for i, a in enumerate(active):
        for b in active[i + 1:]:
            current[(a, b)] = get(clusters[a]["leaf"], clusters[b]["leaf"])

    merges: List[Dict[str, Any]] = []
    next_id = len(usable)
    while len(active) > 1:
        best_pair = None
        best_value = None
        for i, a in enumerate(active):
            for b in active[i + 1:]:
                value = current.get((a, b), current.get((b, a)))
                if value is None:
                    continue
                if best_value is None or value < best_value:
                    best_value, best_pair = value, (a, b)
        if best_pair is None:
            break
        a, b = best_pair
        members = clusters[a]["members"] + clusters[b]["members"]
        node = {"id": next_id, "members": members, "height": best_value / 2.0,
                "left": a, "right": b, "distance": best_value}
        clusters[next_id] = node
        merges.append(node)
        active = [c for c in active if c not in (a, b)]
        for other in active:
            size_a = len(clusters[a]["members"])
            size_b = len(clusters[b]["members"])
            da = current.get((a, other), current.get((other, a), 0.0))
            db = current.get((b, other), current.get((other, b), 0.0))
            current[(next_id, other)] = (da * size_a + db * size_b) / float(size_a + size_b)
        active.append(next_id)
        next_id += 1
    return merges, list(usable), list(excluded)


# ---------------------------------------------------------------------------
# JSON and TSV artefacts
# ---------------------------------------------------------------------------

def sample_json(result: SampleResult, generated: str = "") -> Dict[str, Any]:
    """The machine-readable form of one sample, with the report's own additions.

    ``generated`` defaults to empty so the payload is deterministic; the writer
    passes a timestamp when one is wanted. A golden-file test compares the
    default form.
    """
    payload = result.to_dict()
    headline, provenance = headline_sentence(result)
    payload["report"] = {
        "generated": generated,
        "report_version": __version__,
        "headline": headline,
        "headline_provenance": provenance,
        "caveats": caveat_lines(result),
        "checks": check_rows(all_checks(result), result.sample_id),
        "observations": observation_rows(result),
        "unverified_thresholds": unverified_rows(),
    }
    return payload


def cohort_json(cohort: CohortResult, generated: str = "") -> Dict[str, Any]:
    payload = cohort.to_dict()
    headline, provenance = cohort_headline(cohort)
    payload["report"] = {
        "generated": generated,
        "report_version": __version__,
        "headline": headline,
        "headline_provenance": provenance,
        "pairs_never_compared": _uncompared_pairs(cohort),
        "unverified_thresholds": unverified_rows(),
    }
    return payload


def _columns_of(rows: Sequence[Dict[str, Any]],
                columns: Optional[Sequence[str]] = None) -> List[str]:
    if columns:
        return list(columns)
    ordered: List[str] = []
    for row in rows:
        for key in row:
            if key not in ordered:
                ordered.append(key)
    return ordered


def write_tsv(path: Any, rows: Sequence[Dict[str, Any]],
              columns: Optional[Sequence[str]] = None) -> Path:
    """Write rows as TSV. An empty table still writes its header.

    A header-only file is a statement that the analysis ran and found nothing; a
    missing file is ambiguous between that and a crash, and the difference
    matters when someone is reading a results directory a year later.
    """
    target = Path(path)
    ensure_dir(target.parent)
    names = _columns_of(rows, columns)
    if not names:
        names = ["empty"]
    with open(str(target), "w", encoding="utf-8") as handle:
        handle.write("\t".join(names) + "\n")
        for row in rows:
            handle.write("\t".join(_tsv_cell(row.get(name)) for name in names) + "\n")
    LOG.debug("wrote %s (%d rows)", target, len(rows))
    return target


def write_json(path: Any, payload: Any) -> Path:
    target = Path(path)
    ensure_dir(target.parent)
    with open(str(target), "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False, default=str)
        handle.write("\n")
    LOG.debug("wrote %s", target)
    return target


def write_sample_tables(out_dir: Any, result: SampleResult,
                        generated: str = "") -> List[Path]:
    """Every per-sample artefact: one JSON and one TSV per view."""
    directory = ensure_dir(out_dir)
    stem = safe_name(result.sample_id)
    written = [
        write_json(directory / "{0}.json".format(stem), sample_json(result, generated)),
        write_tsv(directory / "{0}.summary.tsv".format(stem), [result.summary_row()]),
        write_tsv(directory / "{0}.drugs.tsv".format(stem), drug_rows(result)),
        write_tsv(directory / "{0}.catalogue_calls.tsv".format(stem),
                  catalogue_call_rows(result)),
        write_tsv(directory / "{0}.disagreements.tsv".format(stem),
                  disagreement_rows(result)),
        write_tsv(directory / "{0}.variants.tsv".format(stem), variant_rows(result)),
        write_tsv(directory / "{0}.variant_catalogue.tsv".format(stem),
                  variant_catalogue_rows(result)),
        write_tsv(directory / "{0}.checks.tsv".format(stem),
                  check_rows(all_checks(result), result.sample_id)),
        write_tsv(directory / "{0}.observations.tsv".format(stem), observation_rows(result)),
        write_tsv(directory / "{0}.lineage_support.tsv".format(stem),
                  lineage_support_rows(result)),
        write_tsv(directory / "{0}.databases.tsv".format(stem),
                  database_rows(result.database_versions)),
        write_tsv(directory / "{0}.tools.tsv".format(stem),
                  tool_version_rows(result.tool_versions)),
    ]
    return written


def write_cohort_tables(out_dir: Any, cohort: CohortResult,
                        results: Sequence[SampleResult] = (),
                        generated: str = "") -> List[Path]:
    """Cohort artefacts: distances with denominators, clusters, and the matrix."""
    directory = ensure_dir(out_dir)
    written = [
        write_json(directory / "cohort.json", cohort_json(cohort, generated)),
        write_tsv(directory / "cohort.distances.tsv", distance_rows(cohort)),
        write_tsv(directory / "cohort.clusters.tsv", cluster_rows(cohort)),
        write_tsv(directory / "cohort.matrix.tsv", distance_matrix_rows(cohort),
                  columns=["sample"] + list(cohort.samples)),
        write_tsv(directory / "cohort.databases.tsv", database_rows(cohort.database_versions)),
        write_tsv(directory / "cohort.tools.tsv", tool_version_rows(cohort.tool_versions)),
    ]
    if results:
        written.append(write_tsv(directory / "cohort.samples.tsv",
                                 [r.summary_row() for r in results]))
    return written


def write_tables(out_dir: Any, results: Sequence[SampleResult],
                 cohort: Optional[CohortResult] = None,
                 generated: str = "") -> List[Path]:
    """Write every artefact for a run and hand back what was written.

    Raises rather than returning an empty list when there is nothing to write: a
    report step that silently produced no files is the failure mode that turns
    into "the pipeline finished" in somebody's notes.
    """
    if not results and cohort is None:
        raise MjolnirError(
            "no results to write: report.tables.write_tables was called with an "
            "empty result set and no cohort")
    directory = ensure_dir(out_dir)
    written: List[Path] = []
    for result in results:
        written.extend(write_sample_tables(directory, result, generated))
    if cohort is not None:
        written.extend(write_cohort_tables(directory, cohort, results, generated))
    if len(results) > 1:
        written.append(write_tsv(directory / "samples.tsv",
                                 [r.summary_row() for r in results]))
    written.append(write_tsv(directory / "thresholds.tsv", threshold_rows()))
    LOG.info("wrote %d table artefact(s) to %s", len(written), directory)
    return written
