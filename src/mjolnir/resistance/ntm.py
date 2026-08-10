"""NTM resistance: the rules none of the three MTBC catalogues contains.

WHO v2, MTBseq's ResSeq lists and tbdb are all *M. tuberculosis* complex
artefacts. Point a MAC or an *M. abscessus* genome at them and the honest answer
is that they have nothing to say — which is exactly the hole the design was
written against, because MTBseq answers it with ``--resilist NONE`` and a silence
that reads like susceptibility.

So the NTM evidence base is implemented here, directly from the primary
literature, as a table rather than as code: organism (or complex) -> drug ->
determinants -> interpretation -> citation. Adding *M. fortuitum* or *M. kansasii*
means adding rows to :data:`NTM_EVIDENCE`, not writing a new branch, and the
report can print the citation beside every call because the citation is a field
on the rule that produced it.

Three properties of this module are load-bearing.

**A pair that is not in the table returns an explicit refusal.** There is no
default branch that falls through to "susceptible". :func:`evidence_for` returns
None, :func:`no_evidence_call` builds a :class:`~mjolnir.records.DrugCall`
carrying :data:`~mjolnir.config.NTM_NO_EVIDENCE_TEXT`, and the assessment lists
the pair in ``no_evidence_base`` so the report prints the absence as an absence.

**Only a measured genotype can produce a negative statement.** An *erm(41)* C28
sequevar, or the truncated *erm(41)* of *M. abscessus* subsp. *massiliense*, is a
positive finding about the isolate and is reported as one. The absence of an
*rrl* mutation is not: it is ``no-call``, "no resistance determinant detected".
And if the gene that would have carried the determinant was never callable, even
a C28 sequevar is downgraded to ``no-call`` with the missing gene named, because
susceptibility to a macrolide cannot be inferred from a gene nobody looked at.

**rRNA positions are in *E. coli* numbering and nothing else.** The 23S of
*M. avium* calls the macrolide site 2274; *E. coli* calls it 2058; the reference
genome calls it whatever coordinate the operon happens to start at. Matching on a
genomic coordinate across four NTM references would silently mis-call every one
of them, so this module matches on the gene-relative HGVS position that
``normalise.py`` is required to supply in *E. coli* numbering, checks the
reference base while it is there, and raises a visible warning when a variant in
a target gene carries no position it can read — never a quiet miss.

Every numeric threshold used here — the *erm(41)* sequevar position, the *rrl*
and *rrs* positions, the major-variant fraction — comes from ``config.py`` with
its source attached. What this module adds are the per-organism primary
citations that ``config.NTM_TARGETS`` did not need to carry, and the rule table
itself is cross-checked against ``config.NTM_TARGETS`` at import so the two
cannot drift apart unnoticed.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from ..config import (
    COMPLEX_ABSCESSUS,
    COMPLEX_MAC,
    ERM41_INDUCIBLE_ALLELE,
    ERM41_SEQUEVAR_POSITION,
    ERM41_SUSCEPTIBLE_ALLELE,
    ERM41_TRUNCATION_NOTE,
    FASTA_CAPABILITY_LOSS,
    MAJOR_VARIANT_FRACTION,
    NTM_NO_EVIDENCE_TEXT,
    NTM_SPECIES_ALIASES,
    NTM_TARGETS,
    RRL_MACROLIDE_POSITIONS,
    RRS_AMIKACIN_NEIGHBOURS,
    RRS_AMIKACIN_POSITIONS,
    SRC_BASTIAN_2011,
    SRC_DESIGN,
    SRC_NASH_2009,
    SRC_PRAMMANANAN_1998,
    SRC_WALLACE_1996,
    is_major_variant,
    normalise_drug,
    source_for,
)
from ..records import (
    CALL_NO_CALL,
    CALL_R,
    CALL_S,
    CALL_SEVERITY,
    CALL_S_INTERIM,
    CALL_UNCERTAIN,
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MODERATE,
    CONFIDENCE_NONE,
    CONFIDENCES,
    NO_DETERMINANT_TEXT,
    PLATFORM_FASTA,
    PLATFORM_ONT,
    RESISTANCE_CALLS,
    STATUS_PASS,
    STATUS_WARN,
    CatalogueCall,
    Check,
    DrugCall,
    SpeciesCall,
    Variant,
    normalise_platform,
    worst_call,
)
from ..utils import LOG, MjolnirError, percentage, to_jsonable

# ---------------------------------------------------------------------------
# Citations
#
# config.py carries the four NTM sources its thresholds needed. These are the
# additional primary citations for the organisms whose rules live here. No new
# numeric threshold is defined in this module: every position below is imported
# from config.py, where it is registered with its source.
# ---------------------------------------------------------------------------

#: Confirmed 2026-08-10: PMID 8192472, AAC 38(2):381-384. A2058G/C/U (E. coli
#: equivalent) in the peptidyltransferase region of the single-copy 23S rRNA of
#: clarithromycin-resistant M. intracellulare, from paired pre/post-monotherapy
#: isolates. M. chimaera sits inside M. intracellulare sensu lato, which is why
#: it inherits this row.
SRC_MEIER_1994 = (
    "Meier et al. 1994, Antimicrob Agents Chemother 38:381 - identification of "
    "23S rRNA (rrl) mutations in clarithromycin-resistant M. intracellulare; "
    "A2058 (E. coli numbering) to G, C or U"
)

#: Confirmed 2026-08-10: PMID 8592991, AAC 39:2625-2630. The macrolide site is
#: numbered 2274 in the M. avium 23S rRNA gene and 2058 in E. coli; this module
#: works exclusively in the E. coli numbering, which is why that is stated on
#: every rrl determinant.
SRC_NASH_1995 = (
    "Nash & Inderlied 1995, Antimicrob Agents Chemother 39:2625 - genetic basis "
    "of macrolide resistance in M. avium; 23S rRNA position 2274 in M. avium "
    "numbering, equivalent to E. coli 2058"
)

#: Confirmed 2026-08-10: J Clin Microbiol 2013, doi:10.1128/JCM.01612-13.
#: Seven of 462 consecutive MAC isolates with high-level amikacin resistance
#: carried rrs A1408G, and they included M. avium, M. intracellulare and
#: M. chimaera by name — which is what makes the MAC amikacin row citable rather
#: than extrapolated from M. abscessus.
SRC_BROWN_ELLIOTT_2013 = (
    "Brown-Elliott et al. 2013, J Clin Microbiol (doi:10.1128/JCM.01612-13) - "
    "amikacin MICs for 462 MAC isolates; high-level resistance carried rrs "
    "A1408G in M. avium, M. intracellulare and M. chimaera"
)

#: Confirmed 2026-08-10: PMID 20536733, doi:10.1111/j.1348-0421.2010.00221.x.
#: The 274-bp deletion (nucleotides 159-432) that truncates erm(41) in
#: M. abscessus subsp. massiliense, leaving an 81-residue product; 89.8% of
#: massiliense isolates were clarithromycin-susceptible, and those that were not
#: carried an rrl A2058/A2059 substitution instead.
SRC_KIM_2010 = (
    "Kim et al. 2010, Microbiol Immunol (doi:10.1111/j.1348-0421.2010.00221.x) - "
    "M. abscessus subsp. massiliense is differentiated by a 274-bp deletion "
    "(nt 159-432) truncating erm(41); resistant massiliense carry rrl A2058/A2059"
)

#: The name this module signs its calls with. It is not one of config.CATALOGUES:
#: those three are MTBC catalogues, and a report that listed "NTM literature
#: rules" beside them as a fourth catalogue would imply a consensus that does not
#: exist. Here there is one source per rule, and the rule names it.
NTM_RULE_SOURCE = "NTM literature rules (Mjolnir)"

# ---------------------------------------------------------------------------
# Platform statements the design requires (design §7)
# ---------------------------------------------------------------------------

#: SOURCE: design §7, "Two gaps found and not papered over". There is no
#: published R10.4.1-era ONT validation of NTM genotypic DST, so every NTM call
#: made from ONT reads carries this sentence. It is not a hedge: erm(41), rrl and
#: rrs have simply never been benchmarked on this chemistry in this genus.
ONT_NTM_NOT_VALIDATED = (
    "no published R10.4.1-era ONT validation exists for NTM genotypic DST "
    "(erm(41), rrl, rrs), so this call is outside any validated envelope on this "
    "platform"
)

#: SOURCE: design §7. The homopolymer concern for rrs and rrl is an inference
#: from ONT's general error mode, not a measured mycobacterial result, and it is
#: labelled as an inference wherever it is printed.
ONT_RRNA_HOMOPOLYMER_CAVEAT = (
    "no ONT-specific homopolymer-error study exists for rrs or rrl in "
    "mycobacteria; the homopolymer concern here is a reasonable inference from "
    "the general ONT error mode, not a measured result"
)

#: SOURCE: Nash et al. 2009 / Meier et al. 1994. Both M. abscessus and the MAC
#: members carry a single rrn operon, which is why one point substitution in rrl
#: or rrs is enough to confer resistance without a wild-type copy masking it.
#: Stated in the report because it is the reason these rules work at all.
SINGLE_RRN_OPERON_NOTE = (
    "M. abscessus and the MAC members carry a single rrn operon, so one "
    "substitution in rrl or rrs is not diluted by a wild-type copy"
)

#: SOURCE: Nash et al. 2009; Bastian et al. 2011. Inducible resistance is
#: invisible to a short-incubation MIC — this is the whole clinical point of
#: sequevar typing, and the report says it beside every T28 call.
INDUCIBLE_MIC_NOTE = (
    "inducible macrolide resistance is not seen on a 3-day MIC; CLSI M24 "
    "requires extended (14-day) incubation, and a functional erm(41) predicts "
    "the resistance the short read-out misses"
)

#: SOURCE: Mjolnir policy, following the design's fifth house rule. What the
#: report prints when a genotype that would support a negative statement was
#: found, but a gene that could have overturned it was never callable.
UNASSESSED_GENE_TEXT = (
    "{genes} was not assessed in this sample, so no negative statement is made: "
    "an acquired mutation there would change this call"
)

#: The macrolide evidence base is clarithromycin. NTM guidelines treat the class
#: together, so an azithromycin question is answered from the clarithromycin
#: rules with this sentence attached rather than with a refusal — but the drug
#: the evidence was measured on is named.
MACROLIDE_CLASS_NOTE = (
    "the erm(41) and rrl evidence base was established with clarithromycin; NTM "
    "guidelines treat the macrolides as a class, and this call is extended to "
    "azithromycin on that basis"
)

#: Drugs answered from another drug's evidence rows, and the sentence that says
#: so. Nothing else is class-extended: aminoglycosides are not interchangeable
#: here, and rrs 1408 is an amikacin statement.
DRUG_CLASS_EQUIVALENTS: Dict[str, Tuple[str, str]] = {
    "Azithromycin": ("Clarithromycin", MACROLIDE_CLASS_NOTE),
}

#: The drugs an NTM sample is assessed for when the caller does not say. Both
#: have an implemented evidence base for every organism in the table; anything
#: else asked for gets the no-evidence-base answer.
NTM_DEFAULT_DRUGS: Tuple[str, ...] = ("Clarithromycin", "Amikacin")

# ---------------------------------------------------------------------------
# The rule vocabulary
# ---------------------------------------------------------------------------

#: An allele at a fixed position that is a property of the isolate rather than a
#: mutation acquired under treatment — the erm(41) T28C polymorphism.
KIND_SEQUEVAR = "sequevar"
#: A gene rendered non-functional by truncation or deletion.
KIND_TRUNCATION = "truncation"
#: An acquired substitution at a numbered rRNA position.
KIND_SUBSTITUTION = "substitution"
DETERMINANT_KINDS: Tuple[str, ...] = (KIND_SEQUEVAR, KIND_TRUNCATION, KIND_SUBSTITUTION)

#: Numbering the rRNA positions are expressed in. Printed with every rrl/rrs
#: determinant, because the native numbering of the same site differs by
#: hundreds of bases between organisms (M. avium 23S calls it 2274).
ECOLI_NUMBERING = "E. coli numbering"


@dataclass(frozen=True)
class Determinant:
    """One genotype, what it means, and who published that meaning.

    ``call`` is deliberately a value from the shared resistance vocabulary rather
    than a local enum: the NTM path has to produce
    :class:`~mjolnir.records.DrugCall` objects the report and the agent already
    know how to render, and a private vocabulary here would be a second place for
    "no-call" to be paraphrased into "susceptible".

    ``alleles`` empty means *any* substitution at ``positions`` satisfies the
    rule — which is the correct reading for *rrl* 2058/2059, where A→G, A→C and
    A→U all confer macrolide resistance. Where one specific substitution carries
    the evidence and the others do not, the allele is named and the remainder
    fall to :data:`UNLISTED_ALLELE_CALL`.
    """

    gene: str
    kind: str
    label: str
    call: str
    confidence: str
    interpretation: str
    citation: str
    positions: Tuple[int, ...] = ()
    numbering: str = ""
    #: Reference base expected at ``positions``. Checked, not assumed: a variant
    #: whose REF disagrees is evidence that the position was not converted to the
    #: numbering this rule uses, and that is reported rather than matched anyway.
    expected_ref: str = ""
    alleles: Tuple[str, ...] = ()
    caveats: Tuple[str, ...] = ()
    #: False when the citation for this specific rule has not been checked
    #: against the primary document. :func:`unverified_determinants` lists them
    #: and the report marks them, on the same principle as config.unverified().
    verified: bool = True

    def __post_init__(self) -> None:
        if self.kind not in DETERMINANT_KINDS:
            raise MjolnirError(
                "determinant {0!r} has kind {1!r}; expected one of {2}".format(
                    self.label, self.kind, ", ".join(DETERMINANT_KINDS)))
        if self.call not in RESISTANCE_CALLS:
            raise MjolnirError(
                "determinant {0!r} has call {1!r}; expected one of {2}".format(
                    self.label, self.call, ", ".join(RESISTANCE_CALLS)))
        if self.confidence not in CONFIDENCES:
            raise MjolnirError(
                "determinant {0!r} has confidence {1!r}; expected one of {2}".format(
                    self.label, self.confidence, ", ".join(CONFIDENCES)))
        if self.kind == KIND_SUBSTITUTION and not self.positions:
            raise MjolnirError(
                "determinant {0!r} is a substitution rule with no positions".format(
                    self.label))
        if self.positions and not self.numbering:
            raise MjolnirError(
                "determinant {0!r} names positions but no numbering system; a "
                "bare rRNA coordinate is ambiguous between organisms".format(self.label))

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable(asdict(self))


@dataclass(frozen=True)
class DrugEvidence:
    """Everything Mjolnir knows about one organism-drug pair.

    ``genes`` is the set that must have been callable before any negative
    statement is made about the drug. It is stated separately from the
    determinants because a gene with no mutation in it still has to have been
    looked at, and "no determinant detected in a gene nobody sequenced" is not a
    result.
    """

    organism: str
    drug: str
    genes: Tuple[str, ...]
    determinants: Tuple[Determinant, ...]
    citation: str
    method: str = "literature rule table"
    note: str = ""
    caveats: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        data = to_jsonable(asdict(self))
        data["determinants"] = [d.to_dict() for d in self.determinants]
        return data


#: What a substitution at a rule's position but with an allele the rule does not
#: list is called. Uncertain, never resistant and never susceptible: the position
#: is implicated, the specific base is not in the evidence.
UNLISTED_ALLELE_CALL = CALL_UNCERTAIN

# ---------------------------------------------------------------------------
# The determinants
# ---------------------------------------------------------------------------

_ERM41_INDUCIBLE = Determinant(
    gene="erm(41)",
    kind=KIND_SEQUEVAR,
    label="erm(41) T{0} (full-length, functional)".format(ERM41_SEQUEVAR_POSITION),
    call=CALL_R,
    confidence=CONFIDENCE_HIGH,
    positions=(ERM41_SEQUEVAR_POSITION,),
    numbering="erm(41) gene numbering",
    expected_ref=ERM41_INDUCIBLE_ALLELE,
    alleles=(ERM41_INDUCIBLE_ALLELE,),
    interpretation=(
        "a functional erm(41) with T at position {0}: inducible macrolide "
        "resistance. Treat as macrolide-resistant.".format(ERM41_SEQUEVAR_POSITION)
    ),
    citation=SRC_NASH_2009,
    caveats=(INDUCIBLE_MIC_NOTE,),
)

_ERM41_SUSCEPTIBLE = Determinant(
    gene="erm(41)",
    kind=KIND_SEQUEVAR,
    label="erm(41) C{0} (T{0}C polymorphism)".format(ERM41_SEQUEVAR_POSITION),
    call=CALL_S,
    confidence=CONFIDENCE_MODERATE,
    positions=(ERM41_SEQUEVAR_POSITION,),
    numbering="erm(41) gene numbering",
    expected_ref=ERM41_INDUCIBLE_ALLELE,
    alleles=(ERM41_SUSCEPTIBLE_ALLELE,),
    interpretation=(
        "the T{0}C sequevar: erm(41) is non-functional and confers no inducible "
        "macrolide resistance.".format(ERM41_SEQUEVAR_POSITION)
    ),
    citation=SRC_BASTIAN_2011,
    caveats=(
        "this is a statement about inducible resistance only; acquired rrl "
        "2058/2059 resistance is assessed separately and must also be absent",
    ),
)

_ERM41_TRUNCATED = Determinant(
    gene="erm(41)",
    kind=KIND_TRUNCATION,
    label="erm(41) truncated or deleted",
    call=CALL_S,
    confidence=CONFIDENCE_MODERATE,
    numbering="",
    interpretation=ERM41_TRUNCATION_NOTE + (
        " — the 274-bp deletion (nt 159-432) leaves an 81-residue product, and "
        "89.8% of such isolates were clarithromycin-susceptible."
    ),
    citation=SRC_KIM_2010,
    caveats=(
        "massiliense isolates that are macrolide-resistant carry an acquired rrl "
        "A2058/A2059 substitution instead, so rrl must have been assessed",
    ),
)

_RRL_MACROLIDE_ABSCESSUS = Determinant(
    gene="rrl",
    kind=KIND_SUBSTITUTION,
    label="rrl {0} substitution".format("/".join(str(p) for p in RRL_MACROLIDE_POSITIONS)),
    call=CALL_R,
    confidence=CONFIDENCE_HIGH,
    positions=tuple(RRL_MACROLIDE_POSITIONS),
    numbering=ECOLI_NUMBERING,
    expected_ref="A",
    alleles=(),  # any substitution: A->G, A->C and A->U all confer resistance
    interpretation=(
        "an acquired substitution in the 23S rRNA peptidyltransferase region: "
        "constitutive macrolide resistance, independent of erm(41)."
    ),
    citation=SRC_WALLACE_1996,
    caveats=(SINGLE_RRN_OPERON_NOTE,),
)

_RRL_MACROLIDE_MAC = Determinant(
    gene="rrl",
    kind=KIND_SUBSTITUTION,
    label="rrl {0} substitution".format("/".join(str(p) for p in RRL_MACROLIDE_POSITIONS)),
    call=CALL_R,
    confidence=CONFIDENCE_HIGH,
    positions=tuple(RRL_MACROLIDE_POSITIONS),
    numbering=ECOLI_NUMBERING,
    expected_ref="A",
    alleles=(),
    interpretation=(
        "an acquired substitution at the macrolide binding site of the "
        "single-copy 23S rRNA: constitutive macrolide resistance. MAC carries no "
        "erm(41), so this is the whole macrolide evidence base for these species."
    ),
    citation=SRC_MEIER_1994,
    caveats=(SINGLE_RRN_OPERON_NOTE,),
)

_RRS_AMIKACIN_1408 = Determinant(
    gene="rrs",
    kind=KIND_SUBSTITUTION,
    label="rrs {0} substitution".format("/".join(str(p) for p in RRS_AMIKACIN_POSITIONS)),
    call=CALL_R,
    confidence=CONFIDENCE_HIGH,
    positions=tuple(RRS_AMIKACIN_POSITIONS),
    numbering=ECOLI_NUMBERING,
    expected_ref="A",
    alleles=("G",),
    interpretation=(
        "the A1408G substitution in the aminoglycoside binding site of the "
        "single-copy 16S rRNA: high-level amikacin resistance, with cross-"
        "resistance across the 2-deoxystreptamine aminoglycosides."
    ),
    citation=SRC_PRAMMANANAN_1998,
    caveats=(SINGLE_RRN_OPERON_NOTE,),
)

_RRS_AMIKACIN_NEIGHBOURS = Determinant(
    gene="rrs",
    kind=KIND_SUBSTITUTION,
    label="rrs {0} substitution".format(
        "/".join(str(p) for p in RRS_AMIKACIN_NEIGHBOURS)),
    call=CALL_UNCERTAIN,
    confidence=CONFIDENCE_LOW,
    positions=tuple(RRS_AMIKACIN_NEIGHBOURS),
    numbering=ECOLI_NUMBERING,
    expected_ref="",
    alleles=(),
    interpretation=(
        "a substitution neighbouring the amikacin binding site. The evidence "
        "here is thinner than for 1408 and this is graded uncertain by design, "
        "not resistant."
    ),
    citation=SRC_PRAMMANANAN_1998,
    verified=False,  # config registers rrs_amikacin_neighbours as unverified too
)

_RRS_AMIKACIN_1408_MAC = Determinant(
    gene="rrs",
    kind=KIND_SUBSTITUTION,
    label=_RRS_AMIKACIN_1408.label,
    call=CALL_R,
    confidence=CONFIDENCE_HIGH,
    positions=tuple(RRS_AMIKACIN_POSITIONS),
    numbering=ECOLI_NUMBERING,
    expected_ref="A",
    alleles=("G",),
    interpretation=(
        "rrs A1408G in MAC: high-level amikacin resistance. Seven of 462 "
        "consecutive MAC isolates with high-level resistance carried it, across "
        "M. avium, M. intracellulare and M. chimaera."
    ),
    citation=SRC_BROWN_ELLIOTT_2013,
    caveats=(SINGLE_RRN_OPERON_NOTE,),
)

# ---------------------------------------------------------------------------
# The table
#
# Organism or complex -> drug -> evidence. Adding a species is adding rows here.
# The species keys are spelled exactly as config.NTM_TARGETS spells them, and
# _validate_evidence_table() below fails at import if the two ever diverge.
# ---------------------------------------------------------------------------

_ABSCESSUS_MACROLIDE = DrugEvidence(
    organism="Mycobacteroides abscessus",
    drug="Clarithromycin",
    genes=("erm(41)", "rrl"),
    determinants=(_ERM41_INDUCIBLE, _ERM41_SUSCEPTIBLE, _ERM41_TRUNCATED,
                  _RRL_MACROLIDE_ABSCESSUS),
    citation=SRC_BASTIAN_2011,
    note=(
        "two independent mechanisms: inducible resistance from a functional "
        "erm(41), read from the T28C sequevar and from the massiliense "
        "truncation, and constitutive resistance from an acquired rrl "
        "substitution. Both must be assessed before anything negative is said."
    ),
)

_ABSCESSUS_AMIKACIN = DrugEvidence(
    organism="Mycobacteroides abscessus",
    drug="Amikacin",
    genes=("rrs",),
    determinants=(_RRS_AMIKACIN_1408, _RRS_AMIKACIN_NEIGHBOURS),
    citation=SRC_PRAMMANANAN_1998,
    note="a single 16S rRNA substitution accounts for acquired amikacin resistance",
)


def _mac_macrolide(organism: str, citation: str) -> DrugEvidence:
    """Macrolide evidence for one MAC member.

    A function rather than four copies, because the rule is the same rule and the
    only thing that differs is which paper measured it in which species — and
    that difference is exactly what the report has to print.
    """
    determinant = _RRL_MACROLIDE_MAC if citation == SRC_MEIER_1994 else Determinant(
        gene=_RRL_MACROLIDE_MAC.gene,
        kind=_RRL_MACROLIDE_MAC.kind,
        label=_RRL_MACROLIDE_MAC.label,
        call=_RRL_MACROLIDE_MAC.call,
        confidence=_RRL_MACROLIDE_MAC.confidence,
        positions=_RRL_MACROLIDE_MAC.positions,
        numbering=_RRL_MACROLIDE_MAC.numbering,
        expected_ref=_RRL_MACROLIDE_MAC.expected_ref,
        alleles=_RRL_MACROLIDE_MAC.alleles,
        interpretation=_RRL_MACROLIDE_MAC.interpretation,
        citation=citation,
        caveats=_RRL_MACROLIDE_MAC.caveats,
    )
    return DrugEvidence(
        organism=organism,
        drug="Clarithromycin",
        genes=("rrl",),
        determinants=(determinant,),
        citation=citation,
        note=(
            "MAC has no erm(41): there is no inducible macrolide mechanism to "
            "type, and the acquired rrl substitution is the entire evidence base"
        ),
    )


def _mac_amikacin(organism: str) -> DrugEvidence:
    return DrugEvidence(
        organism=organism,
        drug="Amikacin",
        genes=("rrs",),
        determinants=(_RRS_AMIKACIN_1408_MAC, _RRS_AMIKACIN_NEIGHBOURS),
        citation=SRC_BROWN_ELLIOTT_2013,
        note="rrs A1408G is the documented high-level amikacin mechanism in MAC",
    )


NTM_EVIDENCE: Dict[str, Dict[str, DrugEvidence]] = {
    "Mycobacteroides abscessus": {
        "Clarithromycin": _ABSCESSUS_MACROLIDE,
        "Amikacin": _ABSCESSUS_AMIKACIN,
    },
    "Mycobacterium chimaera": {
        # M. chimaera sits within M. intracellulare sensu lato, so it inherits
        # the M. intracellulare rrl evidence; the amikacin paper names it
        # explicitly.
        "Clarithromycin": _mac_macrolide("Mycobacterium chimaera", SRC_MEIER_1994),
        "Amikacin": _mac_amikacin("Mycobacterium chimaera"),
    },
    "Mycobacterium intracellulare": {
        "Clarithromycin": _mac_macrolide("Mycobacterium intracellulare", SRC_MEIER_1994),
        "Amikacin": _mac_amikacin("Mycobacterium intracellulare"),
    },
    "Mycobacterium avium": {
        "Clarithromycin": _mac_macrolide("Mycobacterium avium", SRC_NASH_1995),
        "Amikacin": _mac_amikacin("Mycobacterium avium"),
    },
    # Complex-level rows. A sample that ANI could place in MAC but not below it
    # still has an assessable rrl and rrs — refusing to look because the species
    # is unresolved would discard a real finding. The species-level caveat is
    # attached by the caller path, not by pretending the species is known.
    COMPLEX_MAC: {
        "Clarithromycin": _mac_macrolide(COMPLEX_MAC, SRC_MEIER_1994),
        "Amikacin": _mac_amikacin(COMPLEX_MAC),
    },
    COMPLEX_ABSCESSUS: {
        "Clarithromycin": DrugEvidence(
            organism=COMPLEX_ABSCESSUS,
            drug="Clarithromycin",
            genes=_ABSCESSUS_MACROLIDE.genes,
            determinants=_ABSCESSUS_MACROLIDE.determinants,
            citation=_ABSCESSUS_MACROLIDE.citation,
            note=_ABSCESSUS_MACROLIDE.note,
            caveats=(
                "the subspecies was not resolved; erm(41) is read from the gene "
                "itself, so the call stands, but subspecies-level expectations "
                "(massiliense truncation) could not be cross-checked",
            ),
        ),
        "Amikacin": DrugEvidence(
            organism=COMPLEX_ABSCESSUS,
            drug="Amikacin",
            genes=_ABSCESSUS_AMIKACIN.genes,
            determinants=_ABSCESSUS_AMIKACIN.determinants,
            citation=_ABSCESSUS_AMIKACIN.citation,
            note=_ABSCESSUS_AMIKACIN.note,
        ),
    },
}

#: Complex-level keys, which config.NTM_TARGETS does not carry and which the
#: import-time cross-check therefore skips.
_COMPLEX_KEYS: Tuple[str, ...] = (COMPLEX_MAC, COMPLEX_ABSCESSUS)

#: Organism spellings this module accepts, on top of config.NTM_SPECIES_ALIASES.
#: The subspecies trinomials matter: a species call of "Mycobacteroides
#: abscessus subsp. massiliense" must find the abscessus rows, not fall through
#: to "no evidence base".
_ORGANISM_ALIASES: Dict[str, str] = {
    "mac": COMPLEX_MAC,
    "mycobacterium avium complex": COMPLEX_MAC,
    "m. avium complex": COMPLEX_MAC,
    "m. abscessus complex": COMPLEX_ABSCESSUS,
    "mycobacterium abscessus complex": COMPLEX_ABSCESSUS,
    "mycobacteroides abscessus complex": COMPLEX_ABSCESSUS,
    "m. abscessus": "Mycobacteroides abscessus",
    "m. chimaera": "Mycobacterium chimaera",
    "m. intracellulare": "Mycobacterium intracellulare",
    "m. avium": "Mycobacterium avium",
    "mycobacterium intracellulare subsp. chimaera": "Mycobacterium chimaera",
}
for _key, _value in NTM_SPECIES_ALIASES.items():
    _ORGANISM_ALIASES.setdefault(_key, _value)
for _canonical in list(NTM_EVIDENCE):
    _ORGANISM_ALIASES.setdefault(_canonical.lower(), _canonical)

#: Subspecies of M. abscessus, and what each one implies before a single base is
#: read. These are expectations used to cross-check the erm(41) observation, not
#: substitutes for it: massiliense with an intact functional erm(41) is reported
#: as a discordance rather than resolved in either direction.
ABSCESSUS_SUBSPECIES_EXPECTATION: Dict[str, str] = {
    "massiliense": "erm(41) is expected to be truncated (274-bp deletion)",
    "bolletii": "erm(41) is expected to be full-length and functional (T28)",
    "abscessus": "erm(41) may be T28 (functional) or C28 (T28C, non-functional)",
}


def _validate_evidence_table() -> None:
    """Fail at import if this table and ``config.NTM_TARGETS`` have drifted apart.

    Two tables describing the same evidence base is one table too many, but
    config has to own the gene lists because ``ntm_targets()`` is what the rest
    of the pipeline asks "is there anything to say here at all". So they are kept
    in step by a check that runs every time the module is imported, rather than
    by a comment asking the next person to remember.
    """
    for organism, drugs in NTM_EVIDENCE.items():
        if organism in _COMPLEX_KEYS:
            continue
        expected = NTM_TARGETS.get(organism)
        if expected is None:
            raise MjolnirError(
                "resistance/ntm.py has rules for {0!r} but config.NTM_TARGETS "
                "does not list it; the two tables must agree".format(organism))
        for drug, evidence in drugs.items():
            entry = expected.get(drug)
            if entry is None:
                raise MjolnirError(
                    "resistance/ntm.py has {0} rules for {1!r} but "
                    "config.NTM_TARGETS does not".format(drug, organism))
            if tuple(entry.get("genes", ())) != tuple(evidence.genes):
                raise MjolnirError(
                    "gene targets for {0} / {1} disagree: config says {2}, "
                    "resistance/ntm.py says {3}".format(
                        organism, drug, tuple(entry.get("genes", ())), evidence.genes))
    for organism, drugs in NTM_TARGETS.items():
        for drug in drugs:
            if drug not in NTM_EVIDENCE.get(organism, {}):
                raise MjolnirError(
                    "config.NTM_TARGETS lists {0} for {1!r} but "
                    "resistance/ntm.py implements no rule for it".format(drug, organism))


_validate_evidence_table()


# ---------------------------------------------------------------------------
# Lookup
# ---------------------------------------------------------------------------

def known_organisms() -> List[str]:
    """Organisms and complexes with an implemented evidence base."""
    return sorted(NTM_EVIDENCE)


def organism_key(species: Any) -> str:
    """Resolve a species name, a complex name or a :class:`SpeciesCall` to a table key.

    A :class:`~mjolnir.records.SpeciesCall` resolves through its species name
    first and its complex second, which is the order the evidence supports: a
    named species selects the paper that measured that species, and an
    unresolved MAC or *M. abscessus* complex call still selects the shared rRNA
    rules. Returns "" when nothing matches — the caller's cue to emit the
    no-evidence-base result, not to guess.
    """
    candidates: List[str] = []
    if isinstance(species, SpeciesCall):
        if species.name:
            candidates.append(species.name)
        if species.subspecies and species.name:
            candidates.append("{0} subsp. {1}".format(species.name, species.subspecies))
        if species.complex:
            candidates.append(species.complex)
    elif species:
        candidates.append(str(species))

    for candidate in candidates:
        text = " ".join(str(candidate).split())
        if not text:
            continue
        if text in NTM_EVIDENCE:
            return text
        lowered = text.lower()
        found = _ORGANISM_ALIASES.get(lowered)
        if found:
            return found
        # "Mycobacteroides abscessus subsp. massiliense" -> the species rows.
        trimmed = re.split(r"\s+subsp\.?\s+|\s+var\.?\s+", lowered)[0].strip()
        found = _ORGANISM_ALIASES.get(trimmed)
        if found:
            return found
    return ""


def resolve_drug(drug: str) -> Tuple[str, str]:
    """Canonical drug name and the note explaining any class extension.

    Azithromycin is answered from the clarithromycin rows because NTM guidelines
    treat the macrolides as a class, and the returned note says so. Everything
    else maps to itself with no note: aminoglycosides are not interchangeable
    here, and *rrs* 1408 is an amikacin result.
    """
    canonical = normalise_drug(drug)
    mapped = DRUG_CLASS_EQUIVALENTS.get(canonical)
    if mapped is None:
        return canonical, ""
    return mapped[0], mapped[1]


def evidence_for(species: Any, drug: str) -> Optional[DrugEvidence]:
    """The rule rows for one organism-drug pair, or None when there are none.

    None is the whole point of this function: it is the only honest answer for a
    pair nobody has published on, and every caller is required to turn it into
    :func:`no_evidence_call` rather than into a susceptible result.
    """
    key = organism_key(species)
    if not key:
        return None
    canonical, _note = resolve_drug(drug)
    return NTM_EVIDENCE.get(key, {}).get(canonical)


def supported_drugs(species: Any) -> List[str]:
    """Drugs with an implemented evidence base for this organism."""
    key = organism_key(species)
    return sorted(NTM_EVIDENCE.get(key, {}))


def no_evidence_text(species: Any, drug: str) -> str:
    """The sentence the report prints for a pair with no implemented evidence."""
    name = organism_key(species) or _display_name(species)
    return NTM_NO_EVIDENCE_TEXT.format(species=name or "this organism",
                                       drug=normalise_drug(drug) or "this drug")


def unverified_determinants() -> List[Determinant]:
    """Rules whose citation has not been checked against the primary document."""
    seen: List[Determinant] = []
    for drugs in NTM_EVIDENCE.values():
        for evidence in drugs.values():
            for determinant in evidence.determinants:
                if not determinant.verified and determinant not in seen:
                    seen.append(determinant)
    return seen


def evidence_rows() -> List[Dict[str, Any]]:
    """The whole table, flattened, for the report's methods annex.

    One row per determinant, each carrying the citation that justifies it, so
    the annex can print the evidence base itself rather than only the calls it
    produced.
    """
    rows: List[Dict[str, Any]] = []
    for organism in sorted(NTM_EVIDENCE):
        for drug in sorted(NTM_EVIDENCE[organism]):
            evidence = NTM_EVIDENCE[organism][drug]
            for determinant in evidence.determinants:
                rows.append({
                    "organism": organism,
                    "drug": drug,
                    "gene": determinant.gene,
                    "determinant": determinant.label,
                    "kind": determinant.kind,
                    "numbering": determinant.numbering,
                    "call": determinant.call,
                    "confidence": determinant.confidence,
                    "interpretation": determinant.interpretation,
                    "citation": determinant.citation,
                    "citation_verified": determinant.verified,
                })
    return rows


def _display_name(species: Any) -> str:
    if isinstance(species, SpeciesCall):
        return species.display
    return str(species or "")


# ---------------------------------------------------------------------------
# Observations
# ---------------------------------------------------------------------------

#: erm(41) states this module distinguishes. ``STATE_NOT_ASSESSED`` is not a
#: failure mode to be tidied away: it is the common case for a sample where the
#: gene was not covered, and it must never become "susceptible".
STATE_INDUCIBLE = "functional-inducible"
STATE_SEQUEVAR_SUSCEPTIBLE = "non-functional-sequevar"
STATE_TRUNCATED = "non-functional-truncated"
STATE_ABSENT = "absent"
STATE_NOT_ASSESSED = "not-assessed"


@dataclass
class Erm41Observation:
    """What was seen at *erm(41)*, or an explicit record that nothing was.

    This is a separate input rather than another :class:`Variant` because the
    T28C sequevar is a polymorphism, not a difference from whichever reference
    the run happened to map against: an *M. abscessus* subsp. *abscessus*
    reference makes C28 a variant and T28 invisible, and an *M. abscessus* subsp.
    *massiliense* reference makes the truncation invisible instead. The pileup
    engine reads the base at position 28 and the gene's length directly, and
    hands both here.

    Every field defaults to None, which produces :data:`STATE_NOT_ASSESSED`. A
    default of "present and full length" would be an invention, and it is the
    exact invention this project exists to prevent.
    """

    #: Base observed at config.ERM41_SEQUEVAR_POSITION, upper case.
    sequevar_base: Optional[str] = None
    #: Gene detected at all. False is a finding (some isolates lack erm(41)
    #: entirely); None means nobody looked.
    present: Optional[bool] = None
    #: Truncated or interrupted, by deletion or frameshift.
    truncated: Optional[bool] = None
    #: Size of the deletion when one was measured; the massiliense deletion is
    #: 274 bp, and printing the observed size lets a reader see whether this is
    #: that deletion or another one.
    deletion_bp: Optional[int] = None
    depth: Optional[int] = None
    #: Read support for the sequevar base, when the pileup provides it.
    allele_fraction: Optional[float] = None
    method: str = ""
    note: str = ""

    def __post_init__(self) -> None:
        if self.sequevar_base is not None:
            base = str(self.sequevar_base).strip().upper()
            if len(base) != 1 or base not in "ACGTUN":
                raise MjolnirError(
                    "erm(41) sequevar base {0!r} is not a single nucleotide; the "
                    "pileup at position {1} must supply one base or None".format(
                        self.sequevar_base, ERM41_SEQUEVAR_POSITION))
            self.sequevar_base = "T" if base == "U" else base

    @property
    def state(self) -> str:
        """Which erm(41) determinant this observation satisfies, if any."""
        if self.present is False:
            return STATE_ABSENT
        if self.truncated:
            return STATE_TRUNCATED
        base = self.sequevar_base
        if base == ERM41_INDUCIBLE_ALLELE:
            return STATE_INDUCIBLE
        if base == ERM41_SUSCEPTIBLE_ALLELE:
            return STATE_SEQUEVAR_SUSCEPTIBLE
        if base in ("A", "G"):
            # Neither sequevar. Not in the evidence base, and not silently
            # rounded to the nearest one.
            return STATE_NOT_ASSESSED
        return STATE_NOT_ASSESSED

    @property
    def assessed(self) -> bool:
        return self.state != STATE_NOT_ASSESSED

    def to_dict(self) -> Dict[str, Any]:
        data = to_jsonable(asdict(self))
        data["state"] = self.state
        return data


@dataclass
class DeterminantHit:
    """One determinant, matched, with the observation that matched it."""

    determinant: Determinant
    evidence: str
    variant: Optional[Variant] = None
    allele_fraction: Optional[float] = None
    depth: Optional[int] = None
    caveats: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable({
            "gene": self.determinant.gene,
            "determinant": self.determinant.label,
            "call": self.determinant.call,
            "confidence": self.determinant.confidence,
            "interpretation": self.determinant.interpretation,
            "citation": self.determinant.citation,
            "evidence": self.evidence,
            "variant": self.variant.display if self.variant else "",
            "allele_fraction": self.allele_fraction,
            "depth": self.depth,
            "caveats": self.caveats,
        })


@dataclass
class NTMAssessment:
    """Every NTM drug call for one sample, plus what could not be answered.

    ``no_evidence_base`` is a first-class field rather than a note buried in a
    caveat list, because "we have no evidence base for *M. fortuitum* and
    linezolid" is a statement the report has to make in its own right — the
    design's §5.6 requirement that the tool says so instead of guessing.
    """

    organism: str = ""
    display_name: str = ""
    complex: str = ""
    platform: str = ""
    drugs: List[DrugCall] = field(default_factory=list)
    checks: List[Check] = field(default_factory=list)
    caveats: List[str] = field(default_factory=list)
    hits: List[DeterminantHit] = field(default_factory=list)
    #: (organism, drug, sentence) for every pair with no implemented evidence.
    no_evidence_base: List[Dict[str, str]] = field(default_factory=list)
    #: Citations used, so the report can print them without walking the calls.
    citations: List[str] = field(default_factory=list)

    def drug(self, name: str) -> Optional[DrugCall]:
        canonical = normalise_drug(name)
        for call in self.drugs:
            if call.drug == canonical:
                return call
        return None

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable({
            "organism": self.organism,
            "display_name": self.display_name,
            "complex": self.complex,
            "platform": self.platform,
            "drugs": [d.to_dict() for d in self.drugs],
            "checks": [c.to_dict() for c in self.checks],
            "caveats": self.caveats,
            "determinants": [h.to_dict() for h in self.hits],
            "no_evidence_base": self.no_evidence_base,
            "citations": self.citations,
            "method": NTM_RULE_SOURCE,
        })


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

#: Gene-relative HGVS, in any of the spellings the callers produce:
#: ``r.2058a>g``, ``n.2058A>G``, ``2058A>G``. The leading letter is the HGVS
#: reference type; ``c.`` is accepted too, since some callers emit it for rRNA
#: genes even though it is not strictly correct for a non-coding transcript.
_HGVS_POSITION = re.compile(r"^(?:[a-z]\.)?\(?(-?\d+)")


def rrna_position(variant: Variant) -> Optional[int]:
    """The rRNA position a variant sits at, in the rule's numbering, or None.

    Read from the HGVS name, never from ``variant.pos``. The genomic coordinate
    of *rrl* 2058 differs between every NTM reference Mjolnir can map against,
    so a rule matched on genome coordinates would be right for at most one
    organism and silently wrong for the rest. ``normalise.py`` owns the
    conversion into *E. coli* numbering; this function only reads what it wrote.

    None means the position could not be read, and the caller raises that as a
    visible warning rather than treating it as "no mutation here".
    """
    for text in (variant.hgvs, variant.hgvs_alias):
        if not text:
            continue
        match = _HGVS_POSITION.match(str(text).strip().lower())
        if match:
            return int(match.group(1))
    return None


def _gene_matches(variant: Variant, gene: str) -> bool:
    """Whether a variant belongs to a target gene, tolerating naming variants.

    ``erm(41)``, ``erm41`` and ``ermA_41`` all appear in the wild, and ``rrl``
    turns up as ``23S`` and ``rrlA``. Matching is on a stripped, lower-cased
    token so a naming difference cannot silently drop a resistance call.
    """
    observed = re.sub(r"[^a-z0-9]", "", str(variant.gene or "").lower())
    if not observed:
        return False
    wanted = re.sub(r"[^a-z0-9]", "", gene.lower())
    if observed == wanted:
        return True
    synonyms = _GENE_SYNONYMS.get(gene, ())
    return observed in synonyms


#: Spellings of each target gene seen in NTM annotation, normalised to
#: alphanumerics. Deliberately explicit: a regex loose enough to catch every
#: spelling would also catch genes these rules do not apply to.
_GENE_SYNONYMS: Dict[str, Tuple[str, ...]] = {
    "erm(41)": ("erm41", "erm", "ermabs", "erm41mab"),
    "rrl": ("rrl", "rrla", "23s", "23srrna", "rrnal", "mycrrl"),
    "rrs": ("rrs", "rrsa", "16s", "16srrna", "rrnas", "mycrrs"),
}


def _match_substitution(determinant: Determinant, variant: Variant,
                        position: int) -> Optional[DeterminantHit]:
    """Apply one substitution rule to one variant, or return None."""
    if position not in determinant.positions:
        return None
    alt = str(variant.alt or "").upper().replace("U", "T")
    ref = str(variant.ref or "").upper().replace("U", "T")
    caveats: List[str] = list(determinant.caveats)

    if determinant.expected_ref and ref and ref != determinant.expected_ref:
        # The reference base at 2058/2059/1408 is A in every organism here. A
        # different one means the coordinate was not converted into the
        # numbering this rule uses, and a match made anyway would be a coin flip.
        caveats.append(
            "reference base at {0} {1} is {2}, not the expected {3}: the "
            "position may not have been converted to {4}, and this match is "
            "reported for review rather than relied on".format(
                determinant.gene, position, ref, determinant.expected_ref,
                determinant.numbering))

    call_determinant = determinant
    if determinant.alleles and alt not in determinant.alleles:
        call_determinant = Determinant(
            gene=determinant.gene,
            kind=determinant.kind,
            label="{0} (allele {1}, not in the evidence base)".format(
                determinant.label, alt or "?"),
            call=UNLISTED_ALLELE_CALL,
            confidence=CONFIDENCE_LOW,
            positions=determinant.positions,
            numbering=determinant.numbering,
            expected_ref=determinant.expected_ref,
            alleles=(alt,) if alt else (),
            interpretation=(
                "a substitution at {0} {1}, but to {2} rather than the "
                "substitution the evidence describes ({3}); the position is "
                "implicated, this specific allele is not graded".format(
                    determinant.gene, position, alt or "an unrecorded base",
                    ", ".join(determinant.alleles))),
            citation=determinant.citation,
            caveats=determinant.caveats,
            verified=determinant.verified,
        )

    fraction = variant.allele_fraction
    if fraction is not None and is_major_variant(fraction) is False:
        caveats.append(
            "detected as a minority population at {0}% of reads, below the "
            "{1}% major-variant threshold ({2}); a resistant subpopulation is "
            "still clinically relevant and is reported".format(
                percentage(fraction, 1), percentage(MAJOR_VARIANT_FRACTION, 0),
                source_for("major_variant_fraction")))

    return DeterminantHit(
        determinant=call_determinant,
        evidence="{0} {1}{2}>{3}".format(determinant.gene, position, ref or "?",
                                         alt or "?"),
        variant=variant,
        allele_fraction=fraction,
        depth=variant.depth,
        caveats=caveats,
    )


def _erm41_hit(evidence: DrugEvidence,
               observation: Optional[Erm41Observation]) -> Optional[DeterminantHit]:
    """Turn an erm(41) observation into a determinant hit, or None if unassessed."""
    if observation is None:
        return None
    state = observation.state
    if state == STATE_NOT_ASSESSED:
        return None

    by_state = {
        STATE_INDUCIBLE: _ERM41_INDUCIBLE,
        STATE_SEQUEVAR_SUSCEPTIBLE: _ERM41_SUSCEPTIBLE,
        STATE_TRUNCATED: _ERM41_TRUNCATED,
        STATE_ABSENT: _ERM41_TRUNCATED,
    }
    determinant = by_state[state]
    if determinant not in evidence.determinants:
        return None

    if state == STATE_ABSENT:
        evidence_text = "erm(41) not detected in this genome"
    elif state == STATE_TRUNCATED:
        evidence_text = "erm(41) truncated{0}".format(
            " by a {0} bp deletion".format(observation.deletion_bp)
            if observation.deletion_bp else "")
    else:
        evidence_text = "erm(41) position {0} = {1}".format(
            ERM41_SEQUEVAR_POSITION, observation.sequevar_base)

    caveats = list(determinant.caveats)
    if observation.deletion_bp and observation.deletion_bp != 274:
        caveats.append(
            "the deletion measured here is {0} bp, not the 274 bp described for "
            "M. abscessus subsp. massiliense; the gene is truncated either way, "
            "but this is not that deletion".format(observation.deletion_bp))
    if observation.method:
        caveats.append("erm(41) typed by {0}".format(observation.method))

    return DeterminantHit(
        determinant=determinant,
        evidence=evidence_text,
        variant=None,
        allele_fraction=observation.allele_fraction,
        depth=observation.depth,
        caveats=caveats,
    )


# ---------------------------------------------------------------------------
# Calling
# ---------------------------------------------------------------------------

#: The stem of config.NTM_NO_EVIDENCE_TEXT, derived from it rather than retyped,
#: so :func:`is_no_evidence_base` cannot drift from the sentence it looks for.
NO_EVIDENCE_MARKER = NTM_NO_EVIDENCE_TEXT.split("{")[0].strip()


def no_evidence_call(species: Any, drug: str) -> DrugCall:
    """The explicit refusal for a species-drug pair with no implemented evidence.

    ``target_covered`` is False deliberately. The ``no-call`` label is "no
    resistance determinant detected", which would claim that something was
    looked for and not found; ``target_covered=False`` makes
    :attr:`DrugCall.label` read "not evaluable" instead, which is true. The
    precise reason — that no evidence base exists, not that coverage failed — is
    in ``note`` and in the caveat, and :func:`is_no_evidence_base` lets the
    report recognise these calls and print the exact sentence instead of the
    generic label.
    """
    text = no_evidence_text(species, drug)
    return DrugCall(
        drug=normalise_drug(drug),
        call=CALL_NO_CALL,
        confidence=CONFIDENCE_NONE,
        caveats=[text],
        note=text,
        target_covered=False,
    )


def is_no_evidence_base(call: DrugCall) -> bool:
    """Whether a drug call is the "no evidence base" refusal rather than a result."""
    return NO_EVIDENCE_MARKER in (call.note or "")


def _resolve_call(calls: Sequence[str]) -> str:
    """The drug's call given the determinants that matched.

    Not simply :func:`~mjolnir.records.worst_call`. That function ranks
    ``no-call`` above ``S`` on purpose — for the MTBC path, where a catalogue
    grading one variant "not associated with resistance" must not outrank a gene
    that was never examined. Here the situation is inverted: an ``S`` in this
    module comes from a *typed genotype* (erm(41) C28, or the massiliense
    truncation), which is a positive finding and strictly more informative than
    the absence it would otherwise be collapsed into.

    So anything alarming wins first, exactly as elsewhere; only when nothing
    alarming matched does a negative genotype get to speak, and then the weaker
    of the negative statements wins, because ``S-interim`` claims less than
    ``S``. With no determinants at all the answer is ``no-call`` — never ``S``.
    """
    if not calls:
        return CALL_NO_CALL
    alarming = [c for c in calls if CALL_SEVERITY.get(c, 0) > CALL_SEVERITY[CALL_NO_CALL]]
    if alarming:
        return worst_call(alarming)
    negative = [c for c in calls if c in (CALL_S, CALL_S_INTERIM)]
    if negative:
        return max(negative, key=lambda c: CALL_SEVERITY.get(c, 0))
    return CALL_NO_CALL


def _platform_caveats(platform: str, genes: Sequence[str]) -> List[str]:
    """Platform statements that apply to an NTM call (design §7)."""
    plat = normalise_platform(platform) if platform else ""
    caveats: List[str] = []
    if plat == PLATFORM_ONT:
        caveats.append(ONT_NTM_NOT_VALIDATED)
        if any(gene in ("rrl", "rrs") for gene in genes):
            caveats.append(ONT_RRNA_HOMOPOLYMER_CAVEAT)
    elif plat == PLATFORM_FASTA:
        caveats.append(FASTA_CAPABILITY_LOSS)
    return caveats


def call_drug(species: Any, drug: str, variants: Sequence[Variant] = (), *,
              platform: str = "", erm41: Optional[Erm41Observation] = None,
              callable_genes: Optional[Iterable[str]] = None,
              ) -> Tuple[DrugCall, List[DeterminantHit], List[Check]]:
    """Call one drug for one NTM organism from the rule table.

    Returns the call, the determinants that matched, and the checks that record
    what was and was not assessed. The three are returned together because the
    checks are how the report proves the negative: a macrolide ``no-call`` is
    only meaningful beside a check saying *rrl* was covered, and the same
    ``no-call`` beside a check saying it was not is a different statement
    entirely.
    """
    canonical_drug, class_note = resolve_drug(drug)
    evidence = evidence_for(species, canonical_drug)
    if evidence is None:
        return no_evidence_call(species, drug), [], [Check.not_measured(
            "NTM {0} evidence base".format(normalise_drug(drug) or drug),
            no_evidence_text(species, drug),
            source=SRC_DESIGN, category="resistance")]

    assessed: Set[str] = set()
    if callable_genes is not None:
        assessed = set(str(g) for g in callable_genes)
    checks: List[Check] = []
    hits: List[DeterminantHit] = []
    caveats: List[str] = list(evidence.caveats)
    unreadable: List[str] = []

    # erm(41): a typed genotype, not a difference from a reference.
    if "erm(41)" in evidence.genes:
        hit = _erm41_hit(evidence, erm41)
        if hit is not None:
            hits.append(hit)
            assessed.add("erm(41)")
            checks.append(Check(
                name="erm(41) sequevar typed",
                value=erm41.state if erm41 else STATE_NOT_ASSESSED,
                source=source_for("erm41_sequevar_position"),
                status=STATUS_PASS, category="resistance",
                reading=hit.determinant.interpretation))
        else:
            checks.append(Check.not_measured(
                "erm(41) sequevar typed",
                "erm(41) was not typed in this sample, so inducible macrolide "
                "resistance was neither detected nor excluded",
                source=source_for("erm41_sequevar_position"), category="resistance"))

    # rrl / rrs: acquired substitutions, matched on gene-relative position.
    rrna_genes = [g for g in evidence.genes if g != "erm(41)"]
    for gene in rrna_genes:
        gene_variants = [v for v in variants if _gene_matches(v, gene)]
        for variant in gene_variants:
            position = rrna_position(variant)
            if position is None:
                unreadable.append(variant.display)
                continue
            for determinant in evidence.determinants:
                if determinant.gene != gene or determinant.kind != KIND_SUBSTITUTION:
                    continue
                hit = _match_substitution(determinant, variant, position)
                if hit is not None:
                    hits.append(hit)
                    break
        # A variant called in the gene is itself proof the gene was read, so it
        # counts as coverage even when the caller passed no gene list. Absence of
        # both a variant and a declaration is not coverage, and stays None.
        covered: Optional[bool]
        if gene_variants:
            covered = True
        elif callable_genes is not None:
            covered = gene in assessed
        else:
            covered = None
        checks.append(Check.boolean(
            name="{0} assessed for {1}".format(gene, canonical_drug),
            value=covered,
            source=source_for("rrl_macrolide_positions" if gene == "rrl"
                              else "rrs_amikacin_positions"),
            category="resistance",
            reading=(
                "{0} positions {1} ({2}) were examined".format(
                    gene,
                    "/".join(str(p) for p in (
                        RRL_MACROLIDE_POSITIONS if gene == "rrl"
                        else RRS_AMIKACIN_POSITIONS)),
                    ECOLI_NUMBERING)
                if covered else
                "{0} was not established as callable; no negative statement is "
                "made about {1}".format(gene, canonical_drug)),
            fail_status=STATUS_WARN))
        if covered:
            assessed.add(gene)

    if unreadable:
        # Silence here would look exactly like "no mutation found", which is the
        # one thing this module must never fake.
        message = (
            "{0} variant(s) in the target genes carry no gene-relative HGVS "
            "position ({1}), so they could not be tested against the {2} rules: "
            "{3}".format(len(unreadable), ECOLI_NUMBERING, canonical_drug,
                         ", ".join(unreadable[:5])))
        LOG.warning("%s", message)
        caveats.append(message)
        checks.append(Check.not_measured(
            "{0} positions readable".format(canonical_drug), message,
            source=SRC_DESIGN, category="resistance"))

    missing = [g for g in evidence.genes if g not in assessed]
    call = _resolve_call([h.determinant.call for h in hits])
    confidence = CONFIDENCE_NONE
    for hit in hits:
        if hit.determinant.call == call:
            confidence = hit.determinant.confidence
            break

    # Nothing unmeasured is described as fine (house rule 5). A negative
    # statement stands only when every gene that could overturn it was read.
    if call == CALL_S and missing:
        caveats.append(UNASSESSED_GENE_TEXT.format(genes=" and ".join(missing)))
        call = CALL_NO_CALL
        confidence = CONFIDENCE_NONE

    target_covered: Optional[bool]
    if callable_genes is None and not hits:
        target_covered = None
    elif len(missing) == len(evidence.genes):
        target_covered = False
    else:
        target_covered = True

    for hit in hits:
        caveats.extend(hit.caveats)
    caveats.extend(_platform_caveats(platform, evidence.genes))
    if class_note:
        caveats.append(class_note)
    if not hits and target_covered:
        caveats.append(
            "{0}: {1} examined and no listed determinant found. {2}".format(
                NO_DETERMINANT_TEXT, " and ".join(evidence.genes),
                "This is an absence of a known determinant, not a phenotypic "
                "susceptibility result: NTM genotype-phenotype concordance is "
                "not complete and unlisted mechanisms exist."))

    catalogue_calls = [
        CatalogueCall(
            catalogue=NTM_RULE_SOURCE,
            drug=canonical_drug,
            grade="",  # these rules carry no grading scheme; the citation is the evidence
            comment=hit.determinant.interpretation,
            source=hit.determinant.citation,
            call=hit.determinant.call,
            variant_key=hit.variant.hgvs_key if hit.variant else hit.determinant.label,
            catalogue_version=NTM_RULE_SOURCE,
            matched_by="rule",
            evidence=hit.evidence,
        )
        for hit in hits
    ]

    note = evidence.note
    if missing and call != CALL_NO_CALL:
        note = (note + " " if note else "") + UNASSESSED_GENE_TEXT.format(
            genes=" and ".join(missing))

    drug_call = DrugCall(
        drug=normalise_drug(drug),
        call=call,
        confidence=confidence,
        catalogue_calls=catalogue_calls,
        caveats=_dedupe(caveats),
        supporting_variants=[h.variant.hgvs_key or h.variant.display
                             for h in hits if h.variant],
        who_graded=False,  # WHO v2 is an MTBC catalogue; it grades none of this
        target_covered=target_covered,
        note=note,
    )
    return drug_call, hits, checks


def call_ntm_resistance(species: Any, variants: Sequence[Variant] = (), *,
                        platform: str = "", drugs: Optional[Sequence[str]] = None,
                        erm41: Optional[Erm41Observation] = None,
                        callable_genes: Optional[Iterable[str]] = None,
                        ) -> NTMAssessment:
    """Assess every requested drug for one NTM sample.

    ``species`` may be a :class:`~mjolnir.records.SpeciesCall` or a bare name.
    An MTBC sample raises: the MTBC path is ``resistance/consensus.py`` and its
    three catalogues, and quietly running NTM rules over an *M. tuberculosis*
    genome would produce *rrl*/*rrs* calls with no evidence base behind them.

    Drugs with no implemented evidence are not omitted. They come back as
    ``no-call`` with the refusal sentence attached and are listed in
    ``no_evidence_base``, because a drug missing from a report reads as a drug
    that was fine.
    """
    if isinstance(species, SpeciesCall) and species.is_mtbc:
        raise MjolnirError(
            "NTM resistance rules were applied to an MTBC sample ({0}); MTBC "
            "resistance is called by resistance/consensus.py from the WHO, "
            "MTBseq and tbdb catalogues".format(species.display))

    key = organism_key(species)
    assessment = NTMAssessment(
        organism=key,
        display_name=_display_name(species),
        complex=species.complex if isinstance(species, SpeciesCall) else "",
        platform=normalise_platform(platform) if platform else "",
    )

    if not key:
        assessment.caveats.append(
            "no NTM evidence base is implemented for {0}; every drug below is "
            "reported as an absence of evidence, not as susceptibility".format(
                assessment.display_name or "this organism"))

    wanted = list(drugs) if drugs else list(NTM_DEFAULT_DRUGS)
    for drug in wanted:
        call, hits, checks = call_drug(
            species, drug, variants, platform=platform, erm41=erm41,
            callable_genes=callable_genes)
        assessment.drugs.append(call)
        assessment.hits.extend(hits)
        assessment.checks.extend(checks)
        if evidence_for(species, drug) is None:
            assessment.no_evidence_base.append({
                "organism": key or assessment.display_name,
                "drug": normalise_drug(drug),
                "text": no_evidence_text(species, drug),
            })

    if isinstance(species, SpeciesCall) and species.subspecies:
        expectation = ABSCESSUS_SUBSPECIES_EXPECTATION.get(
            species.subspecies.strip().lower())
        if expectation:
            assessment.caveats.append(
                "subspecies {0}: {1}".format(species.subspecies, expectation))
            if (species.subspecies.strip().lower() == "massiliense"
                    and erm41 is not None
                    and erm41.state == STATE_INDUCIBLE):
                assessment.caveats.append(
                    "discordance: this isolate was called M. abscessus subsp. "
                    "massiliense but carries a full-length functional erm(41). "
                    "Either the subspecies call or the erm(41) typing is wrong, "
                    "and the macrolide call above follows the observed gene")

    assessment.caveats.extend(_platform_caveats(
        assessment.platform, ("rrl", "rrs", "erm(41)")))
    assessment.caveats = _dedupe(assessment.caveats)
    assessment.citations = _dedupe(
        [h.determinant.citation for h in assessment.hits]
        + [e.citation for e in _evidence_used(key, wanted)])
    return assessment


def _evidence_used(key: str, drugs: Sequence[str]) -> List[DrugEvidence]:
    rows: List[DrugEvidence] = []
    for drug in drugs:
        canonical, _note = resolve_drug(drug)
        found = NTM_EVIDENCE.get(key, {}).get(canonical)
        if found is not None:
            rows.append(found)
    return rows


def _dedupe(items: Iterable[str]) -> List[str]:
    """Order-preserving de-duplication; the same caveat twice reads as noise."""
    seen: List[str] = []
    for item in items:
        if item and item not in seen:
            seen.append(item)
    return seen
