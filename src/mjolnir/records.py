"""The data model every stage reads and writes.

This module is the contract. Engines fill these records, the rule layer turns
them into :class:`Check` verdicts, the agent is shown a projection of them, and
the PDF, the HTML and the TSV/JSON artefacts are all built from
:meth:`SampleResult.to_dict` and :meth:`CohortResult.to_dict` so that no two
outputs can disagree about what was found.

Three choices here are load-bearing rather than cosmetic.

**Absence is a value.** Every measurement that can fail to exist is typed
``Optional`` and defaults to ``None``, never to zero and never to a cheerful
string. ``DrugCall.call`` defaults to ``no-call``, not to ``S``: the design says
absence of a catalogued mutation is reported as "no resistance determinant
detected", and the only way to keep that promise through five modules and a
report writer is to make "susceptible" impossible to reach by default.

**The vocabularies are shared.** ``Status`` and the resistance-call set are
module constants because ``resistance/consensus.py``, ``report/pdf.py``,
``agent/discipline.py`` and the cohort code all have to agree on the exact
strings, and a typo in any one of them would otherwise silently produce a
category nobody handles.

**Records carry their own provenance.** A :class:`Check` knows which threshold
it applied and who published it; a :class:`CatalogueCall` knows which catalogue
and which version graded it. The report never has to reconstruct why a number
was called good, because the record it came from already says.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .utils import MjolnirError, natural_key, round_or_none, safe_fraction, to_jsonable

# ---------------------------------------------------------------------------
# Vocabularies
#
# Closed sets, spelled once. Anything that is not in one of these tuples is a
# bug, and the constructors below say so rather than passing it downstream.
# ---------------------------------------------------------------------------

#: Input platforms. ``fasta`` is an assembly: no allele fractions, therefore no
#: heteroresistance, no mixed-infection detection and no heterozygosity-based
#: contamination metric (design §7).
PLATFORM_ILLUMINA = "illumina"
PLATFORM_ONT = "ont"
PLATFORM_FASTA = "fasta"
PLATFORMS: Tuple[str, ...] = (PLATFORM_ILLUMINA, PLATFORM_ONT, PLATFORM_FASTA)

#: Accepted spellings on the command line and in manifests.
PLATFORM_ALIASES: Dict[str, str] = {
    "illumina": PLATFORM_ILLUMINA,
    "ilm": PLATFORM_ILLUMINA,
    "short": PLATFORM_ILLUMINA,
    "pe": PLATFORM_ILLUMINA,
    "paired": PLATFORM_ILLUMINA,
    "ont": PLATFORM_ONT,
    "nanopore": PLATFORM_ONT,
    "long": PLATFORM_ONT,
    "minion": PLATFORM_ONT,
    "fasta": PLATFORM_FASTA,
    "assembly": PLATFORM_FASTA,
    "contigs": PLATFORM_FASTA,
    "genome": PLATFORM_FASTA,
}

#: Status vocabulary for every rule-derived verdict.
STATUS_PASS = "pass"
STATUS_WARN = "warn"
STATUS_FAIL = "fail"
STATUSES: Tuple[str, ...] = (STATUS_PASS, STATUS_WARN, STATUS_FAIL)

#: Ordering used when a panel of checks has to be reduced to one headline: the
#: worst status wins, always.
STATUS_SEVERITY: Dict[str, int] = {STATUS_PASS: 0, STATUS_WARN: 1, STATUS_FAIL: 2}

#: Resistance call vocabulary.
#:
#: ``S`` means a catalogue actively graded the observed variant as not
#: associated with resistance. ``no-call`` means nothing catalogued was found —
#: which is the common case and which the report must render as "no resistance
#: determinant detected", never as susceptibility. ``R-outside-WHO`` is design
#: §5.5 rule 3: another catalogue calls resistance at a variant WHO does not
#: grade. It is surfaced, and it is never presented as equivalent to WHO Group 1.
CALL_R = "R"
CALL_R_INTERIM = "R-interim"
CALL_UNCERTAIN = "Uncertain"
CALL_S_INTERIM = "S-interim"
CALL_S = "S"
CALL_NO_CALL = "no-call"

#: A drug that was never searched, because variant calling produced nothing for
#: the sample. Distinct from ``no-call``, which reports a search that found
#: nothing, and from ``S``, which reports a positive graded finding.
CALL_NOT_ASSESSED = "not-assessed"
CALL_R_OUTSIDE_WHO = "R-outside-WHO"
RESISTANCE_CALLS: Tuple[str, ...] = (
    CALL_R, CALL_R_INTERIM, CALL_UNCERTAIN, CALL_S_INTERIM, CALL_S,
    CALL_NO_CALL, CALL_R_OUTSIDE_WHO,
)

#: How alarming each call is. Used to pick the drug's headline call when several
#: variants contribute, and to order the front-page table.
CALL_SEVERITY: Dict[str, int] = {
    CALL_R: 6,
    CALL_R_OUTSIDE_WHO: 5,
    CALL_R_INTERIM: 4,
    CALL_UNCERTAIN: 3,
    CALL_NO_CALL: 2,
    CALL_S_INTERIM: 1,
    CALL_S: 0,
}

#: The exact clinical wording for each call. The report and the agent both take
#: their phrasing from here so that "no-call" cannot be paraphrased into
#: "susceptible" by whoever writes the next template.
CALL_LABELS: Dict[str, str] = {
    CALL_R: "resistance determinant detected",
    CALL_R_INTERIM: "resistance determinant detected (interim grade)",
    CALL_UNCERTAIN: "variant of uncertain significance detected",
    CALL_S_INTERIM: "variant graded not associated with resistance (interim)",
    CALL_S: "variant graded not associated with resistance",
    CALL_NO_CALL: "no resistance determinant detected",
    CALL_R_OUTSIDE_WHO: "resistance determinant detected outside the WHO catalogue",
}

#: The sentence the design forbids replacing with "susceptible".
NO_DETERMINANT_TEXT = CALL_LABELS[CALL_NO_CALL]

#: Confidence vocabulary shared by species, lineage and drug calls.
CONFIDENCE_HIGH = "high"
CONFIDENCE_MODERATE = "moderate"
CONFIDENCE_LOW = "low"
CONFIDENCE_NONE = "none"
CONFIDENCES: Tuple[str, ...] = (CONFIDENCE_HIGH, CONFIDENCE_MODERATE,
                                CONFIDENCE_LOW, CONFIDENCE_NONE)

#: Sample-validity verdicts. The headline of the contamination panel is one of
#: these, not a purity percentage (design §8): 99.84% *M. tuberculosis* still
#: produced 13 false-positive SNPs, so a percentage is not the statement a
#: reader needs.
VALIDITY_VALID = "valid"
VALIDITY_SUSPECT = "suspect"
VALIDITY_INVALID = "invalid"
VALIDITY_NOT_ASSESSED = "not-assessed"
SAMPLE_VALIDITY: Tuple[str, ...] = (VALIDITY_VALID, VALIDITY_SUSPECT,
                                    VALIDITY_INVALID, VALIDITY_NOT_ASSESSED)

#: Two-tier mixture classification from the genome-wide heterozygous-SNP
#: fraction. Two tiers rather than one cutoff, because the underlying quantity
#: does not separate cleanly (design §8.2).
MIXTURE_SINGLE = "single-strain"
MIXTURE_POSSIBLE = "possible-mixture"
MIXTURE_MIXED = "mixed"
MIXTURE_NOT_ASSESSED = "not-assessed"
MIXTURE_CLASSES: Tuple[str, ...] = (MIXTURE_SINGLE, MIXTURE_POSSIBLE,
                                    MIXTURE_MIXED, MIXTURE_NOT_ASSESSED)

#: Why two catalogues disagree. A disagreement caused purely by legacy codon
#: numbering is a nomenclature artefact and is reported as such (design §5.3);
#: calling it a biological disagreement would manufacture doubt that does not
#: exist.
DISAGREEMENT_NONE = ""
DISAGREEMENT_NOMENCLATURE = "nomenclature"
DISAGREEMENT_BIOLOGICAL = "biological"
DISAGREEMENT_COVERAGE = "coverage"
DISAGREEMENT_KINDS: Tuple[str, ...] = (DISAGREEMENT_NONE, DISAGREEMENT_NOMENCLATURE,
                                       DISAGREEMENT_BIOLOGICAL, DISAGREEMENT_COVERAGE)

#: Variant classes Mjolnir distinguishes. ``lof`` is not a variant type as such;
#: it is the rule-derived class the WHO grading rules act on (§5.4).
VARIANT_SNP = "snp"
VARIANT_MNV = "mnv"
VARIANT_INS = "insertion"
VARIANT_DEL = "deletion"
VARIANT_INDEL = "indel"
VARIANT_LOF = "lof"
VARIANT_TYPES: Tuple[str, ...] = (VARIANT_SNP, VARIANT_MNV, VARIANT_INS,
                                  VARIANT_DEL, VARIANT_INDEL, VARIANT_LOF)


def normalise_platform(value: str) -> str:
    """Map a user-supplied platform name onto the closed set, or raise."""
    key = str(value or "").strip().lower()
    if key in PLATFORM_ALIASES:
        return PLATFORM_ALIASES[key]
    raise MjolnirError(
        "unknown platform {0!r}; expected one of {1}".format(value, ", ".join(PLATFORMS))
    )


def worst_status(statuses: Sequence[str]) -> str:
    """The most severe status in a panel, ``pass`` when the panel is empty."""
    worst = STATUS_PASS
    for status in statuses:
        if STATUS_SEVERITY.get(status, 0) > STATUS_SEVERITY.get(worst, 0):
            worst = status
    return worst


def worst_call(calls: Sequence[str]) -> str:
    """The call a drug takes when several variants contribute to it.

    Seeded from the calls themselves rather than from ``no-call``. ``no-call``
    sits above ``S`` and ``S-interim`` in the severity order, so seeding with it
    added a phantom contribution that outranked every susceptible finding: a
    drug whose only evidence was a graded not-associated variant came out as
    "no determinant detected" - an absence - when a catalogue had in fact
    spoken.
    """
    ranked = [call for call in calls if call]
    if not ranked:
        return CALL_NO_CALL
    chosen = ranked[0]
    for call in ranked[1:]:
        if CALL_SEVERITY.get(call, -1) > CALL_SEVERITY.get(chosen, -1):
            chosen = call
    return chosen


def call_label(call: str) -> str:
    """Clinical wording for a call. Unknown calls are echoed, never softened."""
    return CALL_LABELS.get(call, call)


def pair_key(a: str, b: str) -> Tuple[str, str]:
    """Order-independent key for a sample pair."""
    return (a, b) if a <= b else (b, a)


class _Record(object):
    """Uniform, JSON-safe serialisation for every record in this module."""

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable(asdict(self))


# ---------------------------------------------------------------------------
# Checks — the rule-derived verdicts the agent is allowed to write prose over
# ---------------------------------------------------------------------------

@dataclass
class Check(_Record):
    """One threshold, applied to one measurement, with its source attached.

    The design's second house rule lives in this class: pass/warn/fail is
    computed in Python from a stated threshold *before* the model is called, and
    the model receives finished checks. ``reading`` is the one-sentence
    rule-derived statement that stands in if the model's answer is discarded or
    if no model is reachable at all — which is why it is never optional in
    practice, even though it defaults to empty here.

    ``measured`` is the flag that keeps the fifth house rule honest. A check
    that could not be computed has ``measured=False`` and ``value=None``, and
    every consumer is expected to render it as "not measured", never as a pass.
    """

    name: str
    value: Any = None
    threshold: Any = None
    source: str = ""
    status: str = STATUS_WARN
    reading: str = ""
    #: How value and threshold relate: ">=", "<=", "==", "in", "range".
    comparison: str = ""
    unit: str = ""
    measured: bool = True
    #: Free-form category so the report can group checks: "qc", "contamination",
    #: "resistance", "typing", "cohort", "platform".
    category: str = ""

    def __post_init__(self) -> None:
        if self.status not in STATUSES:
            raise MjolnirError(
                "Check {0!r} has status {1!r}; expected one of {2}".format(
                    self.name, self.status, ", ".join(STATUSES))
            )

    @property
    def ok(self) -> bool:
        return self.status == STATUS_PASS

    @classmethod
    def not_measured(cls, name: str, why: str, source: str = "",
                     category: str = "", status: str = STATUS_WARN) -> "Check":
        """A check that could not be computed, stated as such.

        Defaults to ``warn`` rather than ``pass`` on purpose. A measurement that
        did not happen is a gap in the evidence, and a gap is not a pass.
        """
        return cls(name=name, value=None, threshold=None, source=source,
                   status=status, reading=why, measured=False, category=category)

    @classmethod
    def numeric(cls, name: str, value: Optional[float], *,
                minimum: Optional[float] = None,
                maximum: Optional[float] = None,
                warn_minimum: Optional[float] = None,
                warn_maximum: Optional[float] = None,
                source: str = "", unit: str = "", category: str = "",
                reading: str = "", not_measured_why: str = "") -> "Check":
        """Apply a numeric threshold and record what was applied.

        ``minimum``/``maximum`` are the hard bounds — outside them the check
        fails. ``warn_minimum``/``warn_maximum`` are the softer bounds, so mean
        depth is expressed as ``minimum=10`` (the degraded floor, below which
        the result is not usable) and ``warn_minimum=25`` (the target, below
        which precision and recall are known to degrade).
        """
        if value is None:
            return cls.not_measured(
                name, not_measured_why or "{0} was not measured".format(name),
                source=source, category=category)

        status = STATUS_PASS
        threshold: Any = None
        comparison = ""
        if minimum is not None and value < minimum:
            status, threshold, comparison = STATUS_FAIL, minimum, ">="
        elif maximum is not None and value > maximum:
            status, threshold, comparison = STATUS_FAIL, maximum, "<="
        elif warn_minimum is not None and value < warn_minimum:
            status, threshold, comparison = STATUS_WARN, warn_minimum, ">="
        elif warn_maximum is not None and value > warn_maximum:
            status, threshold, comparison = STATUS_WARN, warn_maximum, "<="
        else:
            for bound, sign in ((warn_minimum, ">="), (minimum, ">="),
                                (warn_maximum, "<="), (maximum, "<=")):
                if bound is not None:
                    threshold, comparison = bound, sign
                    break

        return cls(name=name, value=value, threshold=threshold, source=source,
                   status=status, reading=reading, comparison=comparison,
                   unit=unit, measured=True, category=category)

    @classmethod
    def boolean(cls, name: str, value: Optional[bool], *, expected: bool = True,
                source: str = "", category: str = "", reading: str = "",
                fail_status: str = STATUS_FAIL) -> "Check":
        """A yes/no condition, with None meaning it was never established."""
        if value is None:
            return cls.not_measured(name, "{0} was not established".format(name),
                                    source=source, category=category)
        status = STATUS_PASS if bool(value) == expected else fail_status
        return cls(name=name, value=bool(value), threshold=expected, source=source,
                   status=status, reading=reading, comparison="==", category=category)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

@dataclass
class SampleInput(_Record):
    """One sample as the user described it, before anything has been measured.

    Kept separate from :class:`SampleResult` so that ``seqio.py`` can validate
    and report on a whole cohort's inputs — wrong pairing, an ONT file passed as
    Illumina, a FASTA with no sequence — before a single read is mapped.
    """

    sample_id: str
    platform: str
    paths: List[Path] = field(default_factory=list)
    #: Explicit reference override (``--ref``). None means Mjolnir chooses one
    #: from the species call.
    reference: Optional[Path] = None
    #: Where the sample came from: "cli", "manifest", "directory".
    origin: str = "cli"
    note: str = ""

    def __post_init__(self) -> None:
        self.platform = normalise_platform(self.platform)
        self.paths = [Path(p) for p in self.paths]
        if self.reference is not None:
            self.reference = Path(self.reference)
        if not self.sample_id:
            raise MjolnirError("SampleInput requires a sample id")
        if not self.paths:
            raise MjolnirError("sample {0!r} has no input files".format(self.sample_id))
        if self.platform == PLATFORM_ILLUMINA and len(self.paths) > 2:
            raise MjolnirError(
                "sample {0!r}: Illumina input takes one or two FASTQ files, got {1}".format(
                    self.sample_id, len(self.paths))
            )
        if self.platform in (PLATFORM_ONT, PLATFORM_FASTA) and len(self.paths) != 1:
            raise MjolnirError(
                "sample {0!r}: {1} input takes exactly one file, got {2}".format(
                    self.sample_id, self.platform, len(self.paths))
            )

    #: Aliases, because "id" and "sample" are both natural at a call site.
    @property
    def id(self) -> str:  # noqa: A003 - deliberate alias
        return self.sample_id

    @property
    def sample(self) -> str:
        return self.sample_id

    @property
    def is_paired(self) -> bool:
        return self.platform == PLATFORM_ILLUMINA and len(self.paths) == 2

    @property
    def is_reads(self) -> bool:
        return self.platform in (PLATFORM_ILLUMINA, PLATFORM_ONT)

    @property
    def r1(self) -> Optional[Path]:
        return self.paths[0] if self.platform == PLATFORM_ILLUMINA else None

    @property
    def r2(self) -> Optional[Path]:
        return self.paths[1] if self.is_paired else None

    @property
    def assembly(self) -> Optional[Path]:
        return self.paths[0] if self.platform == PLATFORM_FASTA else None


# ---------------------------------------------------------------------------
# Variants and catalogue calls
# ---------------------------------------------------------------------------

@dataclass
class Variant(_Record):
    """One observed difference from the reference, with how well it is seen.

    Both keys the design requires are properties rather than stored fields, so
    they cannot drift from the coordinates they are derived from:
    ``coordinate_key`` is WHO's own matching protocol (exact match on chromosome,
    position, reference and alternative nucleotide against NC_000962.3), and
    ``hgvs_key`` is the cross-catalogue join key, because MTBseq and tbdb do not
    share WHO's coordinate table (§5.3).

    ``allele_fraction`` is a 0-1 fraction and is ``None`` for FASTA input, where
    it does not exist. It is not 1.0: an assembly consensus is not evidence that
    an allele was fixed in the population, and writing 1.0 there would turn a
    capability loss into a confident-looking number.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = ""
    hgvs: str = ""
    depth: Optional[int] = None
    allele_fraction: Optional[float] = None
    is_major: Optional[bool] = None
    source_caller: str = ""
    #: Reads supporting the alternative allele. Kept beside the fraction because
    #: the platform thresholds are expressed in reads (>=3 Illumina, >=5 ONT),
    #: not in fractions.
    alt_reads: Optional[int] = None
    ref_reads: Optional[int] = None
    qual: Optional[float] = None
    variant_type: str = VARIANT_SNP
    #: Predicted consequence: synonymous, missense, stop_gained, frameshift,
    #: inframe_deletion, upstream, non_coding.
    effect: str = ""
    locus_tag: str = ""
    #: HGVS as the catalogue spells it, when a legacy alias was needed to match.
    hgvs_alias: str = ""
    #: VCF FILTER entries and Mjolnir's own filters, kept rather than dropped so
    #: the annex can show why a variant was not used.
    filters: List[str] = field(default_factory=list)
    #: Whether this position falls inside the repeat/low-complexity mask.
    masked: bool = False
    #: Per-catalogue grades for this variant, for the annex table.
    catalogue_calls: List["CatalogueCall"] = field(default_factory=list)
    note: str = ""

    def __post_init__(self) -> None:
        self.pos = int(self.pos)
        if self.variant_type not in VARIANT_TYPES:
            raise MjolnirError(
                "variant at {0}:{1} has type {2!r}; expected one of {3}".format(
                    self.chrom, self.pos, self.variant_type, ", ".join(VARIANT_TYPES))
            )

    @property
    def coordinate_key(self) -> Tuple[str, int, str, str]:
        """WHO's matching key: (chromosome, position, reference, alternative)."""
        return (self.chrom, self.pos, self.ref.upper(), self.alt.upper())

    @property
    def hgvs_key(self) -> str:
        """``<gene>_<hgvs>`` — the cross-catalogue join key.

        Empty when the variant has no gene or no HGVS name, which is a real
        state (an intergenic position outside any catalogued region) and must
        not be faked with a coordinate string, or it would join against nothing
        and look like a missing grade.
        """
        if not self.gene or not self.hgvs:
            return ""
        return "{0}_{1}".format(self.gene, self.hgvs)

    @property
    def is_indel(self) -> bool:
        return self.variant_type in (VARIANT_INS, VARIANT_DEL, VARIANT_INDEL) \
            or len(self.ref) != len(self.alt)

    @property
    def display(self) -> str:
        return self.hgvs_key or "{0}:{1}{2}>{3}".format(self.chrom, self.pos, self.ref, self.alt)

    def to_dict(self) -> Dict[str, Any]:
        data = to_jsonable(asdict(self))
        data["coordinate_key"] = "{0}:{1}{2}>{3}".format(
            self.chrom, self.pos, self.ref.upper(), self.alt.upper())
        data["hgvs_key"] = self.hgvs_key
        data["allele_fraction"] = round_or_none(self.allele_fraction, 4)
        return data


@dataclass
class CatalogueCall(_Record):
    """What one catalogue says about one variant for one drug.

    ``grade`` is the source's own grading string, verbatim — for WHO that is one
    of the five numeric-prefixed strings, spelled with a spaced ASCII hyphen.
    MTBseq has no grading at all, so its ``grade`` is empty and its ``call`` is
    R or no-call; that asymmetry is stated in the report rather than hidden
    behind a manufactured grade (§5.5).
    """

    catalogue: str
    drug: str
    grade: str = ""
    comment: str = ""
    source: str = ""
    call: str = CALL_NO_CALL
    #: The variant this call was made about, as its HGVS join key.
    variant_key: str = ""
    #: Version and checksum of the catalogue file the row came from, so two
    #: installations that disagree can be told apart.
    catalogue_version: str = ""
    catalogue_checksum: str = ""
    #: How the row was matched: "coordinate", "hgvs", "rule", "alias".
    matched_by: str = ""
    #: Evidence text the source supplies (WHO's confidence-grading counts, tbdb's
    #: confidence field).
    evidence: str = ""

    def __post_init__(self) -> None:
        if self.call not in RESISTANCE_CALLS:
            raise MjolnirError(
                "catalogue {0!r} produced call {1!r} for {2}; expected one of {3}".format(
                    self.catalogue, self.call, self.drug, ", ".join(RESISTANCE_CALLS))
            )


@dataclass
class DrugCall(_Record):
    """The consensus statement for one drug in one sample.

    ``call`` defaults to ``no-call`` and ``confidence`` to ``none``, so a drug
    that was never evaluated cannot be mistaken for one that came back clear.
    ``caveats`` is where the platform consequences land (design §7): ONT
    under-detection of minor variants, the suppressed *fbiC*/delamanid call, the
    indel-driven LoF caveat, and the total loss of allele fractions on FASTA.
    """

    drug: str
    call: str = CALL_NO_CALL
    confidence: str = CONFIDENCE_NONE
    catalogue_calls: List[CatalogueCall] = field(default_factory=list)
    caveats: List[str] = field(default_factory=list)
    disagreement: bool = False
    disagreement_kind: str = DISAGREEMENT_NONE
    #: HGVS keys of the variants that produced the call.
    supporting_variants: List[str] = field(default_factory=list)
    #: True when WHO graded at least one of the supporting variants. False with
    #: ``call == "R-outside-WHO"`` is exactly the §5.5 rule-3 situation.
    who_graded: bool = False
    #: WHO grade string that anchored the call, when there was one.
    who_grade: str = ""
    #: Suppressed by an epistasis rule (mmpL5 LoF over Rv0678; eis coding LoF
    #: over the eis promoter). Names the rule, and the report prints it.
    suppressed_by: str = ""
    #: Level of resistance and cross-resistance, which come from the catalogue's
    #: Comment column and not from the grade (§5.4).
    level: str = ""
    cross_resistance: List[str] = field(default_factory=list)
    #: Whether the target regions for this drug were callable at all. False means
    #: the drug is unevaluable in this sample, which is not the same as no-call.
    target_covered: Optional[bool] = None
    note: str = ""

    def __post_init__(self) -> None:
        if self.call not in RESISTANCE_CALLS:
            raise MjolnirError(
                "drug {0!r} has call {1!r}; expected one of {2}".format(
                    self.drug, self.call, ", ".join(RESISTANCE_CALLS))
            )
        if self.confidence not in CONFIDENCES:
            raise MjolnirError(
                "drug {0!r} has confidence {1!r}; expected one of {2}".format(
                    self.drug, self.confidence, ", ".join(CONFIDENCES))
            )
        if self.disagreement_kind not in DISAGREEMENT_KINDS:
            raise MjolnirError(
                "drug {0!r} has disagreement kind {1!r}; expected one of {2}".format(
                    self.drug, self.disagreement_kind,
                    ", ".join(k for k in DISAGREEMENT_KINDS if k))
            )

    @property
    def is_resistant(self) -> bool:
        return self.call in (CALL_R, CALL_R_INTERIM, CALL_R_OUTSIDE_WHO)

    @property
    def label(self) -> str:
        """The clinical wording. Never "susceptible" for a no-call."""
        if self.target_covered is False:
            return "not evaluable: target regions not callable"
        return call_label(self.call)

    def calls_by_catalogue(self, catalogue: str) -> List[CatalogueCall]:
        return [c for c in self.catalogue_calls if c.catalogue == catalogue]


# ---------------------------------------------------------------------------
# Typing
# ---------------------------------------------------------------------------

@dataclass
class SpeciesCall(_Record):
    """What organism this is, and how far down the tree the evidence reaches.

    ``resolved_to_species`` is the field the report keys on. Within the MTBC the
    honest answer is often "M. tuberculosis complex, not resolved below complex
    by ANI" — the members sit at 99.21-99.92% ANI and are heterotypic synonyms
    in NCBI taxonomy — and within MAC the *M. chimaera* / *M. intracellulare*
    boundary needs marker SNPs rather than ANI alone (§6). Printing a species
    name in either case would be an invention.

    ``method`` names the evidence, never a taxonomic read classifier: a Kraken2
    row saying "M. bovis 3.2%" is not a species identification and Mjolnir must
    not print one.
    """

    name: str = "unresolved"
    complex: str = ""
    method: str = ""
    ani: Optional[float] = None
    confidence: str = CONFIDENCE_NONE
    resolved_to_species: bool = False
    #: Accession or filename of the best-matching reference in the ANI set.
    reference: str = ""
    #: Aligned fraction reported alongside ANI. A 99.9% ANI over 2% of the
    #: genome is not a species call, and only this field can say so.
    aligned_fraction: Optional[float] = None
    #: Runner-up matches, so the annex can show the margin.
    candidates: List[Dict[str, Any]] = field(default_factory=list)
    subspecies: str = ""
    caveats: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.confidence not in CONFIDENCES:
            raise MjolnirError(
                "species call has confidence {0!r}; expected one of {1}".format(
                    self.confidence, ", ".join(CONFIDENCES))
            )

    @property
    def is_mtbc(self) -> bool:
        return self.complex.upper() in ("MTBC", "M. TUBERCULOSIS COMPLEX")

    @property
    def display(self) -> str:
        if self.resolved_to_species and self.name:
            return self.name
        return self.complex or self.name or "unresolved"


@dataclass
class LineageCall(_Record):
    """MTBC lineage and sublineage from the SNP barcode, or NTM subspecies.

    ``barcode_sites_supporting`` over ``barcode_sites_callable`` is the evidence
    the report prints beside the lineage, and it is the reason the *M. bovis*
    caveat can be stated quantitatively: *M. bovis* is defined by very few
    phylogenetic SNPs — 23 in SNP-IT — so a coverage gap over a handful of
    positions is enough to lose or invent the call.

    ``mixed_lineages`` non-empty means more than one lineage's defining SNPs
    were supported, which is a mixed-infection signal and is reported as one
    rather than resolved by picking the larger.
    """

    lineage: str = ""
    sublineage: str = ""
    barcode_sites_supporting: int = 0
    barcode_sites_total: int = 0
    barcode_sites_callable: int = 0
    is_bcg: bool = False
    is_animal: bool = False
    #: Which animal-adapted member, when one was called: bovis, caprae, orygis,
    #: pinnipedii, microti, BCG.
    animal_variant: str = ""
    caveats: List[str] = field(default_factory=list)
    #: Barcode scheme and its version, e.g. "tbdb barcode.bed <commit>".
    scheme: str = ""
    method: str = "pileup at barcode positions"
    confidence: str = CONFIDENCE_NONE
    mixed_lineages: List[str] = field(default_factory=list)
    #: Per-site detail for the annex: position, expected allele, observed
    #: allele, depth, allele fraction.
    support: List[Dict[str, Any]] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.confidence not in CONFIDENCES:
            raise MjolnirError(
                "lineage call has confidence {0!r}; expected one of {1}".format(
                    self.confidence, ", ".join(CONFIDENCES))
            )

    @property
    def support_fraction(self) -> Optional[float]:
        """Supported over callable barcode sites, or None when none were callable."""
        return safe_fraction(self.barcode_sites_supporting, self.barcode_sites_callable)

    @property
    def callable_fraction(self) -> Optional[float]:
        return safe_fraction(self.barcode_sites_callable, self.barcode_sites_total)

    @property
    def display(self) -> str:
        if self.sublineage:
            return self.sublineage
        return self.lineage or "not determined"


# ---------------------------------------------------------------------------
# QC and contamination
# ---------------------------------------------------------------------------

@dataclass
class QCMetrics(_Record):
    """Coverage, mapping and composition metrics for one sample.

    Every field is ``Optional`` and every one of them is ``None`` on an input
    that cannot produce it — a FASTA has no depth, no mapped fraction and no
    unambiguous-base fraction, and the report says the capability is absent
    rather than printing zeros.

    ``coverage_evenness`` needs its definition stated wherever it is shown,
    because "evenness" is not a standard quantity: Mjolnir defines it as the
    fraction of reference positions whose depth lies within a stated band around
    the mean, and ``evenness_definition`` carries that sentence with the number.
    """

    mean_depth: Optional[float] = None
    median_depth: Optional[float] = None
    #: Fraction of the reference covered at 1x / 10x / the configured minimum.
    breadth_1x: Optional[float] = None
    breadth_10x: Optional[float] = None
    breadth_min_depth: Optional[float] = None
    coverage_evenness: Optional[float] = None
    evenness_definition: str = ""
    mapped_fraction: Optional[float] = None
    gc_content: Optional[float] = None
    #: MTBseq's de-facto heterozygosity filter, surfaced rather than silently
    #: applied: the fraction of called positions where one allele carries at
    #: least the unambiguity threshold of the reads.
    unambiguous_fraction: Optional[float] = None
    total_reads: Optional[int] = None
    mapped_reads: Optional[int] = None
    duplicate_fraction: Optional[float] = None
    mean_read_length: Optional[float] = None
    mean_base_quality: Optional[float] = None
    #: Reference the metrics are against, and its length.
    reference: str = ""
    reference_length: Optional[int] = None
    checks: List[Check] = field(default_factory=list)

    @property
    def status(self) -> str:
        return worst_status([c.status for c in self.checks])


@dataclass
class ContaminationResult(_Record):
    """What could honestly be measured about purity, and the verdict it supports.

    The headline is ``verdict``, one of the sample-validity strings, and not a
    purity percentage. The design is explicit about why: a sample that was
    99.84% *M. tuberculosis* still produced 13 false-positive SNPs across 12
    genes, and 5% *M. avium* contamination produced 3,325 false-positive variant
    SNPs, so any gate at 1% or 5% is a coarse instrument.

    ``screen_informative`` is the second refusal made structural. A Kraken2 run
    against a standard or capped index has a measured sensitivity of 0.0731 for
    *M. tuberculosis* reads on real Illumina data, so its output is not a
    contamination screen for this genus. When this flag is False the report says
    the screen was uninformative — it does not report a clean result.

    The mapping and composition mirrors (``mapped_fraction`` and below) are
    copied from :class:`QCMetrics` so that the contamination panel is complete
    on its own in the JSON artefact and in the agent observation.
    """

    #: Minor-allele frequency across lineage-defining SNP sets (F2/F47).
    f2: Optional[float] = None
    f47: Optional[float] = None
    lineage_het_sites: Optional[int] = None
    lineage_sites_examined: Optional[int] = None
    #: Genome-wide heterozygous-SNP fraction under the MixInfect-style filters.
    het_snp_fraction: Optional[float] = None
    het_snp_count: Optional[int] = None
    snp_sites_examined: Optional[int] = None
    mixture_class: str = MIXTURE_NOT_ASSESSED
    #: MTBseq's unambiguous-base fraction, surfaced rather than discarded.
    unambiguous_fraction: Optional[float] = None
    #: ANI-assigned non-target read fraction, with a label saying what the
    #: assignment could actually resolve ("genus", "complex", "species").
    non_target_fraction: Optional[float] = None
    non_target_resolution: str = ""
    non_target_labels: List[Dict[str, Any]] = field(default_factory=list)
    #: Mirrors of the QC signals that bear on purity.
    mapped_fraction: Optional[float] = None
    coverage_breadth: Optional[float] = None
    coverage_evenness: Optional[float] = None
    gc_content: Optional[float] = None
    #: The verdict, and the sentence that justifies it.
    verdict: str = VALIDITY_NOT_ASSESSED
    verdict_reason: str = ""
    screen_informative: bool = False
    screen_method: str = ""
    #: Why the screen is or is not informative — the index that was used, and
    #: the confidence it was run at.
    screen_note: str = ""
    caveats: List[str] = field(default_factory=list)
    checks: List[Check] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.verdict not in SAMPLE_VALIDITY:
            raise MjolnirError(
                "contamination verdict {0!r}; expected one of {1}".format(
                    self.verdict, ", ".join(SAMPLE_VALIDITY))
            )
        if self.mixture_class not in MIXTURE_CLASSES:
            raise MjolnirError(
                "mixture class {0!r}; expected one of {1}".format(
                    self.mixture_class, ", ".join(MIXTURE_CLASSES))
            )


# ---------------------------------------------------------------------------
# Provenance and interpretation
# ---------------------------------------------------------------------------

@dataclass
class DatabaseVersion(_Record):
    """One database as it was on disk during this run.

    A catalogue-version mismatch between two installations changes calls, so
    every artefact prints version and checksum for everything it consulted.
    """

    name: str
    version: str = "unknown"
    checksum: str = ""
    path: str = ""
    licence: str = ""
    citation: str = ""
    url: str = ""
    fetched: str = ""
    note: str = ""


@dataclass
class Interpretation(_Record):
    """The prose layer, and whether it survived the discipline rules.

    ``rule_only`` True is the normal state when no model host is reachable, and
    it is also what happens when a model answer violated the discipline rules
    and was discarded. Either way the report prints the rule-derived summary and
    says which of the two it is — the output never silently loses the model and
    never silently keeps a bad answer.
    """

    headline: str = ""
    body: str = ""
    rule_only: bool = True
    discarded_reason: str = ""
    model: str = ""
    host: str = ""
    playbook: str = ""


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

@dataclass
class SampleResult(_Record):
    """Everything Mjolnir determined about one sample.

    :meth:`to_dict` is the single source for the JSON artefact, the TSV rows and
    the agent observation. Building them separately is how three outputs come to
    disagree about the same run, and the design forbids it.
    """

    sample_id: str
    platform: str = PLATFORM_ILLUMINA
    inputs: List[str] = field(default_factory=list)
    reference: str = ""
    species: SpeciesCall = field(default_factory=SpeciesCall)
    lineage: LineageCall = field(default_factory=LineageCall)
    variants: List[Variant] = field(default_factory=list)
    drugs: List[DrugCall] = field(default_factory=list)
    qc: QCMetrics = field(default_factory=QCMetrics)
    contamination: ContaminationResult = field(default_factory=ContaminationResult)
    checks: List[Check] = field(default_factory=list)
    #: Caveats that apply to the whole sample — the per-platform consequences of
    #: design §7 land here as well as on the individual drug calls.
    caveats: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    tool_versions: Dict[str, str] = field(default_factory=dict)
    database_versions: List[DatabaseVersion] = field(default_factory=list)
    interpretation: Optional[Interpretation] = None
    mjolnir_version: str = ""
    profile: str = "clinical"
    started_at: str = ""
    finished_at: str = ""
    runtime_seconds: float = 0.0
    status: str = STATUS_WARN

    def __post_init__(self) -> None:
        self.platform = normalise_platform(self.platform)

    @property
    def id(self) -> str:  # noqa: A003 - deliberate alias
        return self.sample_id

    @property
    def sample(self) -> str:
        return self.sample_id

    def drug(self, name: str) -> Optional[DrugCall]:
        for call in self.drugs:
            if call.drug.lower() == str(name).lower():
                return call
        return None

    def resistant_drugs(self) -> List[DrugCall]:
        return [d for d in self.drugs if d.is_resistant]

    def disagreements(self) -> List[DrugCall]:
        return [d for d in self.drugs if d.disagreement]

    def all_checks(self) -> List[Check]:
        """Every check from every panel, in report order."""
        return list(self.checks) + list(self.qc.checks) + list(self.contamination.checks)

    def unmeasured(self) -> List[str]:
        """Names of the checks that could not be computed.

        The agent is handed this list and is forbidden from describing anything
        in it as absent, normal or fine.
        """
        return [c.name for c in self.all_checks() if not c.measured]

    def overall_status(self) -> str:
        return worst_status([c.status for c in self.all_checks()])

    def variants_measured(self) -> bool:
        """Whether this sample got as far as producing a variant call.

        Read from the checks, not from ``len(self.variants)``: a genuine
        zero-variant sample and a sample whose caller died both have an empty
        list, and only one of them supports the phrase "no determinant detected".
        """
        for check in self.all_checks():
            if check.name in ("variant_calling", "sample_analysed") and not check.measured:
                return False
        return True

    def summary_row(self) -> Dict[str, Any]:
        """One flat row for the cohort TSV."""
        row: Dict[str, Any] = {
            "sample": self.sample_id,
            "platform": self.platform,
            "species": self.species.display,
            "species_confidence": self.species.confidence,
            "resolved_to_species": self.species.resolved_to_species,
            "complex": self.species.complex,
            "lineage": self.lineage.display,
            "bcg": self.lineage.is_bcg,
            "animal_lineage": self.lineage.animal_variant or "",
            "mean_depth": round_or_none(self.qc.mean_depth, 1),
            "breadth_min_depth": round_or_none(self.qc.breadth_min_depth, 4),
            "mapped_fraction": round_or_none(self.qc.mapped_fraction, 4),
            "sample_validity": self.contamination.verdict,
            "mixture_class": self.contamination.mixture_class,
            "contamination_screen_informative": self.contamination.screen_informative,
            "n_variants": len(self.variants),
            "status": self.status or self.overall_status(),
        }
        # A sample whose caller failed must not be written as "no-call", the
        # token this codebase defines as "searched and found nothing". The flat
        # TSV is what gets loaded into a spreadsheet and read without the
        # caveats, so the distinction has to survive into the cell itself.
        measured = self.variants_measured()
        for call in sorted(self.drugs, key=lambda d: natural_key(d.drug)):
            row["drug_" + call.drug.lower().replace(" ", "_")] = (
                call.call if measured else CALL_NOT_ASSESSED)
        row["disagreements"] = ";".join(d.drug for d in self.disagreements())
        row["warnings"] = "; ".join(self.warnings)
        return row

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable({
            "mjolnir_version": self.mjolnir_version,
            "sample": self.sample_id,
            "platform": self.platform,
            "profile": self.profile,
            "inputs": self.inputs,
            "reference": self.reference,
            "status": self.status or self.overall_status(),
            "species": self.species.to_dict(),
            "lineage": self.lineage.to_dict(),
            "qc": self.qc.to_dict(),
            "contamination": self.contamination.to_dict(),
            "drugs": [d.to_dict() for d in self.drugs],
            "variants": [v.to_dict() for v in self.variants],
            "checks": [c.to_dict() for c in self.checks],
            "unmeasured": self.unmeasured(),
            "caveats": self.caveats,
            "warnings": self.warnings,
            "interpretation": self.interpretation.to_dict() if self.interpretation else None,
            "tool_versions": self.tool_versions,
            "database_versions": [d.to_dict() for d in self.database_versions],
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "runtime_seconds": round(self.runtime_seconds, 2),
        })


@dataclass
class PairwiseDistance(_Record):
    """A masked SNP distance between two samples, with its denominator.

    The denominator is not optional decoration. As in ``tesseract-ai``'s cgMLST
    output, 12 differences over 4.1 Mb of shared callable sequence and 12 over
    400 kb are not the same statement, and a matrix that prints only the
    numerator invites the reader to treat them as if they were.
    """

    sample_a: str
    sample_b: str
    snps: Optional[int] = None
    shared_callable_sites: Optional[int] = None
    #: Masked positions excluded from this pair's comparison.
    masked_sites: Optional[int] = None
    note: str = ""

    @property
    def key(self) -> Tuple[str, str]:
        return pair_key(self.sample_a, self.sample_b)

    @property
    def snps_per_mb(self) -> Optional[float]:
        if self.snps is None or not self.shared_callable_sites:
            return None
        return self.snps * 1e6 / self.shared_callable_sites


@dataclass
class Cluster(_Record):
    """A set of samples within the clustering threshold of one another."""

    cluster_id: str
    members: List[str] = field(default_factory=list)
    threshold: Optional[int] = None
    max_distance: Optional[int] = None
    min_shared_callable_sites: Optional[int] = None
    note: str = ""

    @property
    def size(self) -> int:
        return len(self.members)


@dataclass
class CohortResult(_Record):
    """Distances, clusters and the shared denominator behind each one.

    The clustering threshold is a field, not a constant: the TB conventions are
    5 and 12 SNPs, and the prior *M. chimaera* run on this machine used 6. The
    basis for whichever value was used is printed beside the clusters, because a
    cluster is a claim about transmission and the threshold is most of the claim.
    """

    samples: List[str] = field(default_factory=list)
    pairs: List[PairwiseDistance] = field(default_factory=list)
    clusters: List[Cluster] = field(default_factory=list)
    threshold: Optional[int] = None
    threshold_basis: str = ""
    #: Mask applied before counting, and how much of the reference it removed.
    mask_name: str = ""
    masked_sites: Optional[int] = None
    masked_fraction: Optional[float] = None
    #: Positions in the joint variant table before masking.
    joint_sites: Optional[int] = None
    reference: str = ""
    checks: List[Check] = field(default_factory=list)
    caveats: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    tool_versions: Dict[str, str] = field(default_factory=dict)
    database_versions: List[DatabaseVersion] = field(default_factory=list)
    interpretation: Optional[Interpretation] = None
    mjolnir_version: str = ""
    runtime_seconds: float = 0.0

    def _index(self) -> Dict[Tuple[str, str], PairwiseDistance]:
        return dict((p.key, p) for p in self.pairs)

    def pair(self, a: str, b: str) -> Optional[PairwiseDistance]:
        return self._index().get(pair_key(a, b))

    def distance(self, a: str, b: str) -> Optional[int]:
        """SNP distance between two samples, or None when it was not computed.

        None rather than zero. Two samples that were never compared are not
        identical, and a distance matrix that fills its gaps with zeros produces
        clusters that do not exist.
        """
        if a == b:
            return 0
        found = self.pair(a, b)
        return None if found is None else found.snps

    def shared_callable_sites(self, a: str, b: str) -> Optional[int]:
        found = self.pair(a, b)
        return None if found is None else found.shared_callable_sites

    def distance_matrix(self) -> Dict[str, Dict[str, Optional[int]]]:
        """Square dict-of-dicts, with None where a pair was not compared."""
        index = self._index()
        matrix: Dict[str, Dict[str, Optional[int]]] = {}
        for a in self.samples:
            row: Dict[str, Optional[int]] = {}
            for b in self.samples:
                if a == b:
                    row[b] = 0
                    continue
                found = index.get(pair_key(a, b))
                row[b] = None if found is None else found.snps
            matrix[a] = row
        return matrix

    def cluster_of(self, sample: str) -> Optional[Cluster]:
        for cluster in self.clusters:
            if sample in cluster.members:
                return cluster
        return None

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable({
            "mjolnir_version": self.mjolnir_version,
            "samples": self.samples,
            "reference": self.reference,
            "threshold": self.threshold,
            "threshold_basis": self.threshold_basis,
            "mask": {
                "name": self.mask_name,
                "masked_sites": self.masked_sites,
                "masked_fraction": round_or_none(self.masked_fraction, 4),
            },
            "joint_sites": self.joint_sites,
            "distance_matrix": self.distance_matrix(),
            "pairs": [p.to_dict() for p in self.pairs],
            "clusters": [c.to_dict() for c in self.clusters],
            "checks": [c.to_dict() for c in self.checks],
            "caveats": self.caveats,
            "warnings": self.warnings,
            "interpretation": self.interpretation.to_dict() if self.interpretation else None,
            "tool_versions": self.tool_versions,
            "database_versions": [d.to_dict() for d in self.database_versions],
            "runtime_seconds": round(self.runtime_seconds, 2),
        })
