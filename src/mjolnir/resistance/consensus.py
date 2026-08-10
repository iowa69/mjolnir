"""Three catalogues, reconciled into one statement per drug.

This is design §5.5, and it is short because the rule is short. WHO is the
anchor: where WHO grades the variant, the WHO grade is the Mjolnir call, because
WHO's is the only source with a published, systematically-derived grading. Where
WHO does not grade it and another catalogue calls resistance, the drug is
reported as ``R (outside WHO catalogue)`` — surfaced, never silently dropped,
and never rendered as though it were a Group 1 call. Where the catalogues
conflict, the drug carries a disagreement flag and the annex gets all three side
by side. And where nothing catalogued was found at all, the answer is "no
resistance determinant detected", which is not the same sentence as
"susceptible" and is never allowed to become it.

That last one is the reason ``DrugCall.call`` defaults to ``no-call`` rather
than to ``S`` and the reason this module has no code path that produces ``S``
from an absence. ``S`` here means a catalogue actively graded an observed
variant as not associated with resistance. Nothing else may say it.

Two suppression steps run before the anchor rule, and both leave a record.

*Epistasis* comes from ``rules.py``: *mmpL5* loss of function abrogates *Rv0678*
for bedaquiline and clofazimine, and a coding loss of function in *eis*
abrogates *eis* promoter mutations for amikacin and kanamycin. The abrogated
catalogue calls stay attached to the drug so the annex can show them; they are
simply not allowed to decide the call, and ``DrugCall.suppressed_by`` names the
rule that took them out.

*Platform* comes from ``config.is_suppressed_on_platform``: on ONT an *fbiC*
tandem-repeat deletion driving a delamanid call is suppressed, because such
calls were 47.2% of every discordant drug classification in the 508-isolate
ONT-vs-Illumina comparison. That is not a small correction — it is the single
largest source of platform disagreement in the study the thresholds come from.

A suppressed drug does not fall through to ``no-call``. A determinant *was*
detected; what the rule says is that it is not expected to produce the
phenotype, or that it is probably an artefact of the sequencing platform.
"No resistance determinant detected" would be false, and the closed vocabulary's
only honest landing place is ``Uncertain`` — a variant was seen and Mjolnir is
not calling resistance on it. The caveat printed beside it says which of the two
reasons applies.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..config import (
    ANCHOR_CATALOGUE,
    CATALOGUE_MTBSEQ,
    CATALOGUE_WHO,
    DRUGS,
    FASTA_CAPABILITY_LOSS,
    MTBSEQ_ASYMMETRY_NOTE,
    ONT_INDEL_CAVEAT,
    ONT_MINOR_VARIANT_CAVEAT,
    WHO_GRADE_1,
    WHO_GRADE_2,
    WHO_GRADE_3,
    WHO_GRADE_4,
    WHO_GRADE_5,
    is_suppressed_on_platform,
    normalise_drug,
    normalise_grade,
)
from ..records import (
    CALL_NO_CALL,
    CALL_R,
    CALL_R_INTERIM,
    CALL_R_OUTSIDE_WHO,
    CALL_S,
    CALL_S_INTERIM,
    CALL_SEVERITY,
    CALL_UNCERTAIN,
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MODERATE,
    CONFIDENCE_NONE,
    DISAGREEMENT_BIOLOGICAL,
    DISAGREEMENT_COVERAGE,
    DISAGREEMENT_NOMENCLATURE,
    DISAGREEMENT_NONE,
    NO_DETERMINANT_TEXT,
    PLATFORM_FASTA,
    PLATFORM_ILLUMINA,
    PLATFORM_ONT,
    CatalogueCall,
    DrugCall,
    Variant,
    normalise_platform,
)
from ..utils import natural_key
from . import rules

# ---------------------------------------------------------------------------
# Policy constants
#
# These are mappings and vocabulary choices rather than measured thresholds, so
# they are not numbers in config.py's sense — but each still names where it came
# from, because the report prints the consequence of every one of them.
# ---------------------------------------------------------------------------

#: SOURCE: Mjolnir policy (config.SRC_POLICY), from design §5.5. Calls that
#: assert resistance. ``R-outside-WHO`` is here because it is a resistance
#: statement; it is separated from ``R`` everywhere it is *displayed*, never in
#: whether it counts as a finding.
RESISTANT_CALLS: Tuple[str, ...] = (CALL_R, CALL_R_INTERIM, CALL_R_OUTSIDE_WHO)

#: SOURCE: Mjolnir policy (config.SRC_POLICY). The call a drug takes when the
#: only evidence for resistance was suppressed by an epistasis rule or by a
#: platform artefact rule. Not ``no-call``: a determinant was detected, and
#: ``no-call`` renders as "no resistance determinant detected", which would be
#: false. Not ``S``: nothing graded this sample susceptible.
SUPPRESSED_CALL = CALL_UNCERTAIN

#: SOURCE: Mjolnir policy (config.SRC_POLICY), from design §5.5. How much weight
#: a WHO grade carries. The two definitive grades are high confidence, the two
#: interim grades are moderate because WHO itself marks them provisional, and
#: "Uncertain significance" is low by definition.
CONFIDENCE_FOR_GRADE: Dict[str, str] = {
    WHO_GRADE_1: CONFIDENCE_HIGH,
    WHO_GRADE_2: CONFIDENCE_MODERATE,
    WHO_GRADE_3: CONFIDENCE_LOW,
    WHO_GRADE_4: CONFIDENCE_MODERATE,
    WHO_GRADE_5: CONFIDENCE_HIGH,
}

#: Confidence, most to least. Used only by :func:`downgrade`.
CONFIDENCE_ORDER: Tuple[str, ...] = (
    CONFIDENCE_HIGH, CONFIDENCE_MODERATE, CONFIDENCE_LOW, CONFIDENCE_NONE)

#: SOURCE: design §5.5 rule 3. The sentence that must accompany every
#: ``R-outside-WHO`` call, so it can never be read as a WHO Group 1 finding.
OUTSIDE_WHO_TEXT = (
    "this resistance call comes from a catalogue other than {0}, which does not "
    "grade the variant. It is reported because dropping it would hide a finding, "
    "and it is not equivalent to a WHO Group 1 call: it has no published, "
    "systematically-derived grading behind it".format(ANCHOR_CATALOGUE)
)

#: SOURCE: design §5.5 rule 5, and records.NO_DETERMINANT_TEXT. Spelled out here
#: because the front page prints the short label and the annex prints this.
NO_DETERMINANT_EXPLANATION = (
    "{0}: no variant in this sample matched any of the catalogues for this drug. "
    "That is an absence of evidence, not evidence of susceptibility, and it is "
    "not a phenotypic result".format(NO_DETERMINANT_TEXT)
)

#: SOURCE: design §5.5 and house rule 5. Emitted when nobody established whether
#: the drug's target regions were callable at all.
COVERAGE_UNKNOWN_TEXT = (
    "coverage of this drug's target regions was not established, so an absence "
    "of determinants here has not been shown to be an absence rather than a gap"
)

#: SOURCE: design §5.5 known failure modes. Printed once per sample by the
#: report; exposed here so the sentence has one home.
INHERITED_BLIND_SPOTS_TEXT = (
    "anchoring on {0} inherits its blind spots by construction: a genuinely "
    "novel mechanism cannot be adjudicated by any of the three catalogues, and a "
    "catalogue-version mismatch between installations changes calls, which is "
    "why every catalogue's version and checksum is printed".format(ANCHOR_CATALOGUE)
)


def downgrade(confidence: str, steps: int = 1) -> str:
    """One step less confident, floored at ``none``."""
    try:
        index = CONFIDENCE_ORDER.index(confidence)
    except ValueError:
        return CONFIDENCE_NONE
    return CONFIDENCE_ORDER[min(index + max(0, steps), len(CONFIDENCE_ORDER) - 1)]


def is_resistant_call(call: str) -> bool:
    return call in RESISTANT_CALLS


def no_determinant_statement(drug: str) -> str:
    """The exact sentence for a drug with nothing catalogued found."""
    return "{0}: {1}".format(normalise_drug(drug), NO_DETERMINANT_EXPLANATION)


# ---------------------------------------------------------------------------
# Contributions
# ---------------------------------------------------------------------------

@dataclass
class Contribution:
    """One catalogue's statement about one variant, for one drug.

    Kept as a pair rather than flattened because both halves are needed later:
    the :class:`CatalogueCall` decides the call, and the :class:`Variant` decides
    whether a platform caveat applies, whether the position was masked, and
    whether two catalogues are arguing about a genome or about a spelling.
    """

    variant: Variant
    call: CatalogueCall
    #: Empty while the contribution counts. Otherwise the rule that removed it.
    suppressed_by: str = ""

    @property
    def active(self) -> bool:
        return not self.suppressed_by

    @property
    def catalogue(self) -> str:
        return self.call.catalogue

    @property
    def variant_key(self) -> str:
        return self.call.variant_key or self.variant.hgvs_key or self.variant.display


def collect_contributions(variants: Sequence[Variant], drug: str) -> List[Contribution]:
    """Every catalogue call about *drug*, across a sample's variants.

    Matching is on the normalised drug name because the three catalogues do not
    spell drugs the same way, and an unnormalised "rifampin" would look like a
    drug nobody else called — which reads in the output as a catalogue that had
    no opinion rather than as a join that failed.
    """
    canonical = normalise_drug(drug)
    found: List[Contribution] = []
    for variant in variants:
        for call in variant.catalogue_calls:
            if normalise_drug(call.drug) == canonical:
                found.append(Contribution(variant=variant, call=call))
    return found


def _apply_suppressions(contributions: Sequence[Contribution], drug: str,
                        platform: str,
                        suppressions: Sequence["rules.Suppression"]) -> List[str]:
    """Mark suppressed contributions in place; return the caveats to print."""
    canonical = normalise_drug(drug)
    caveats: List[str] = []

    for suppression in suppressions:
        if normalise_drug(suppression.drug) != canonical:
            continue
        if not suppression.confident:
            # Reported, not applied: the two variants may be in different
            # subpopulations, and a suppression that assumed otherwise would
            # silently remove a real resistance call.
            caveats.append("{0} ({1}); {2}".format(
                suppression.why or suppression.rule, suppression.rule,
                suppression.caveat or rules.MIXED_SUBPOPULATION_CAVEAT))
            continue
        hit = False
        for contribution in contributions:
            if contribution.suppressed_by:
                continue
            if suppression.suppresses(contribution.variant_key):
                contribution.suppressed_by = suppression.rule
                hit = True
        if hit:
            caveats.append("{0} ({1}); the abrogated variant is reported in the "
                           "annex and is not counted towards the call".format(
                               suppression.why or suppression.rule, suppression.rule))

    for contribution in contributions:
        if contribution.suppressed_by:
            continue
        reason = is_suppressed_on_platform(
            contribution.variant.gene or "", canonical, platform)
        if reason:
            contribution.suppressed_by = "platform:{0}".format(normalise_platform(platform))
            if reason not in caveats:
                caveats.append(reason)

    return caveats


# ---------------------------------------------------------------------------
# The anchor rule
# ---------------------------------------------------------------------------

def _worst(contributions: Sequence[Contribution]) -> Optional[Contribution]:
    """The contribution carrying the most alarming call, or None for an empty set."""
    chosen: Optional[Contribution] = None
    for contribution in contributions:
        if chosen is None:
            chosen = contribution
            continue
        if CALL_SEVERITY.get(contribution.call.call, -1) > \
                CALL_SEVERITY.get(chosen.call.call, -1):
            chosen = contribution
    return chosen


def speaking(contributions: Sequence[Contribution]) -> List[Contribution]:
    """The contributions that actually say something.

    A ``no-call`` row is a catalogue declining to comment, not a statement, and
    it must not be reduced together with the graded rows: ``no-call`` outranks
    ``S`` in :data:`records.CALL_SEVERITY` — correctly, since an unexamined drug
    is more alarming than a graded-susceptible one — so a catalogue holding one
    graded row and one blank would otherwise have its grade masked by its blank.
    """
    return [c for c in contributions if c.active and c.call.call != CALL_NO_CALL]


def catalogue_headline(contributions: Sequence[Contribution],
                       catalogue: str) -> Optional[Contribution]:
    """What one catalogue says about this drug, reduced to its worst call.

    None when that catalogue said nothing, which is the state MTBseq is in for
    every variant absent from its flat list.
    """
    return _worst([c for c in speaking(contributions) if c.catalogue == catalogue])


def anchor_call(contributions: Sequence[Contribution]) -> Tuple[str, Optional[Contribution]]:
    """Apply design §5.5 rules 2 and 3 to the active contributions.

    Returns the call and the contribution that produced it. WHO first, always;
    then, only if WHO is silent for this drug, anything another catalogue says —
    with any resistance statement rewritten to ``R-outside-WHO`` so that it can
    never be mistaken for a graded call.
    """
    active = speaking(contributions)
    who = _worst([c for c in active if c.catalogue == CATALOGUE_WHO])
    if who is not None:
        return who.call.call, who

    others = [c for c in active if c.catalogue != CATALOGUE_WHO]
    if not others:
        return CALL_NO_CALL, None

    resistant = [c for c in others if is_resistant_call(c.call.call)]
    if resistant:
        return CALL_R_OUTSIDE_WHO, _worst(resistant)

    chosen = _worst(others)
    return (chosen.call.call if chosen is not None else CALL_NO_CALL), chosen


# ---------------------------------------------------------------------------
# Disagreement
# ---------------------------------------------------------------------------

def _stance(call: str) -> str:
    """Reduce a call to the three positions two catalogues can take."""
    if is_resistant_call(call):
        return "resistant"
    if call in (CALL_S, CALL_S_INTERIM):
        return "not-resistant"
    if call == CALL_UNCERTAIN:
        return "uncertain"
    return ""


def catalogue_stances(contributions: Sequence[Contribution]) -> Dict[str, str]:
    """Per-catalogue stance, skipping catalogues that made no call.

    A catalogue with no row is not disagreeing with anything. That matters most
    for MTBseq, whose list is flat: it can only produce R or no-call, so reading
    its silence as "not resistant" would invent a dissent it never expressed —
    which is precisely the asymmetry ``config.MTBSEQ_ASYMMETRY_NOTE`` describes.
    """
    stances: Dict[str, str] = {}
    for catalogue in sorted(set(c.catalogue for c in speaking(contributions))):
        headline = catalogue_headline(contributions, catalogue)
        if headline is None:
            continue
        stance = _stance(headline.call.call)
        if stance:
            stances[catalogue] = stance
    return stances


def has_disagreement(contributions: Sequence[Contribution]) -> bool:
    """Whether two catalogues that both spoke took different positions."""
    return len(set(catalogue_stances(contributions).values())) > 1


def disagreement_kind(contributions: Sequence[Contribution],
                      target_covered: Optional[bool] = None) -> str:
    """Why two catalogues disagree — and specifically, whether it is real.

    A disagreement caused purely by legacy codon numbering is a nomenclature
    artefact (design §5.3). Reporting it as a biological disagreement would
    manufacture doubt that does not exist, and a clinician reading "the
    catalogues disagree about rifampicin" has no way to tell the two apart
    unless the record does it for them.

    The nomenclature test is coordinate identity: if every contribution that
    disagrees sits on the same ``(chrom, pos, ref, alt)`` and the catalogues
    reached it under different HGVS spellings — or reached it through an alias —
    then the argument is about a name.
    """
    active = [c for c in contributions if c.active and _stance(c.call.call)]
    if not active:
        return DISAGREEMENT_NONE
    if not has_disagreement(active):
        return DISAGREEMENT_NONE
    if target_covered is False:
        return DISAGREEMENT_COVERAGE
    if all(c.variant.masked for c in active):
        return DISAGREEMENT_COVERAGE

    coordinates = set(c.variant.coordinate_key for c in active)
    if len(coordinates) == 1:
        spellings = set(c.variant_key for c in active)
        aliased = any(c.call.matched_by == "alias" or c.variant.hgvs_alias
                      for c in active)
        if len(spellings) > 1 or aliased:
            return DISAGREEMENT_NOMENCLATURE
    return DISAGREEMENT_BIOLOGICAL


def side_by_side(contributions: Sequence[Contribution]) -> List[Dict[str, Any]]:
    """One annex row per catalogue call, for the full three-way comparison.

    Suppressed contributions are included and flagged rather than filtered, so
    the annex shows what was found as well as what was counted.
    """
    rows: List[Dict[str, Any]] = []
    for contribution in contributions:
        call = contribution.call
        rows.append({
            "catalogue": call.catalogue,
            "catalogue_version": call.catalogue_version,
            "catalogue_checksum": call.catalogue_checksum,
            "drug": normalise_drug(call.drug),
            "variant": contribution.variant_key,
            "coordinate": "{0}:{1}{2}>{3}".format(*contribution.variant.coordinate_key),
            "grade": call.grade,
            "call": call.call,
            "matched_by": call.matched_by,
            "comment": call.comment,
            "evidence": call.evidence,
            "suppressed_by": contribution.suppressed_by,
            "counted": contribution.active,
        })
    rows.sort(key=lambda row: (natural_key(row["catalogue"]), natural_key(row["variant"])))
    return rows


# ---------------------------------------------------------------------------
# Platform consequences (design §7)
# ---------------------------------------------------------------------------

def platform_caveats_for(call: str, contributions: Sequence[Contribution],
                         platform: str) -> List[str]:
    """The per-drug consequences of the sequencing platform.

    Three of them, all from design §7 and all stated rather than silently
    applied: ONT indel calls are about 16.6% uncorroborated, so an indel-driven
    call — which is most of what the loss-of-function rules produce — carries a
    caveat; ONT under-detects minor variants, so an ONT result that asserts
    absence must say that absence of a minor variant is not absence of a
    subpopulation; and an assembly has no allele fractions at all, which is a
    capability loss and not a clean result.
    """
    plat = normalise_platform(platform)
    caveats: List[str] = []
    active = [c for c in contributions if c.active]

    if plat == PLATFORM_ONT:
        if any(c.variant.is_indel for c in active):
            caveats.append(ONT_INDEL_CAVEAT)
        rests_on_minor = any(c.variant.is_major is False for c in active)
        if rests_on_minor or not is_resistant_call(call):
            caveats.append(ONT_MINOR_VARIANT_CAVEAT)
    elif plat == PLATFORM_FASTA:
        caveats.append(FASTA_CAPABILITY_LOSS)

    return caveats


# ---------------------------------------------------------------------------
# The consensus
# ---------------------------------------------------------------------------

def consensus_for_drug(drug: str, variants: Sequence[Variant], *,
                       platform: str = PLATFORM_ILLUMINA,
                       suppressions: Sequence["rules.Suppression"] = (),
                       target_covered: Optional[bool] = None) -> DrugCall:
    """One drug's consensus statement, from every catalogue call in the sample.

    The order of operations is the whole rule and it is not interchangeable:
    collect every catalogue's calls, remove the ones an epistasis or platform
    rule abrogates *while recording that it did*, then let WHO anchor whatever is
    left, then decide whether the catalogues that did speak actually disagreed.
    Suppressing after anchoring would let a suppressed variant set the call;
    anchoring before collecting would hide the other catalogues from the annex.
    """
    canonical = normalise_drug(drug)
    plat = normalise_platform(platform)
    contributions = collect_contributions(variants, canonical)

    caveats = _apply_suppressions(contributions, canonical, plat, suppressions)
    suppressed = [c for c in contributions if not c.active]

    call, source = anchor_call(contributions)
    who_graded = source is not None and source.catalogue == CATALOGUE_WHO
    who_grade = normalise_grade(source.call.grade) if who_graded and source else ""

    # A suppressed determinant must not read as an absent one.
    if suppressed and not is_resistant_call(call):
        had_resistance = any(is_resistant_call(c.call.call) for c in suppressed)
        if had_resistance:
            call = SUPPRESSED_CALL

    active = [c for c in contributions if c.active]
    disagreement = has_disagreement(active)
    kind = disagreement_kind(contributions, target_covered) if disagreement \
        else DISAGREEMENT_NONE

    # Confidence.
    if target_covered is False:
        confidence = CONFIDENCE_NONE
    elif call == CALL_NO_CALL:
        confidence = CONFIDENCE_NONE
    elif call == CALL_R_OUTSIDE_WHO:
        confidence = CONFIDENCE_LOW
    elif call == SUPPRESSED_CALL and suppressed:
        confidence = CONFIDENCE_LOW
    elif who_grade:
        confidence = CONFIDENCE_FOR_GRADE.get(who_grade, CONFIDENCE_LOW)
    else:
        confidence = CONFIDENCE_LOW
    if disagreement:
        confidence = downgrade(confidence)
    if plat == PLATFORM_ONT and any(c.variant.is_indel for c in active) \
            and is_resistant_call(call):
        confidence = downgrade(confidence)

    # Caveats, in the order a reader needs them.
    if call == CALL_R_OUTSIDE_WHO:
        caveats.append(OUTSIDE_WHO_TEXT)
    if any(c.catalogue == CATALOGUE_MTBSEQ for c in active) and \
            (disagreement or call == CALL_R_OUTSIDE_WHO):
        caveats.append(MTBSEQ_ASYMMETRY_NOTE)
    if target_covered is None and call == CALL_NO_CALL:
        caveats.append(COVERAGE_UNKNOWN_TEXT)
    caveats.extend(platform_caveats_for(call, contributions, plat))

    # Level and cross-resistance come from the Comment column, not the grade.
    level = ""
    cross: List[str] = []
    for contribution in active:
        comment = contribution.call.comment
        level = level or rules.level_from_comment(comment)
        cross.extend(rules.cross_resistance_from_comment(comment, canonical))

    note = ""
    if call == CALL_NO_CALL:
        note = NO_DETERMINANT_EXPLANATION
    elif suppressed and call == SUPPRESSED_CALL:
        note = ("a catalogued determinant was detected and then suppressed by "
                "{0}; it is reported here rather than dropped, and no resistance "
                "is predicted from it".format(
                    "; ".join(sorted(set(c.suppressed_by for c in suppressed)))))

    return DrugCall(
        drug=canonical,
        call=call,
        confidence=confidence,
        catalogue_calls=[c.call for c in contributions],
        caveats=_dedupe(caveats),
        disagreement=disagreement,
        disagreement_kind=kind,
        supporting_variants=sorted(
            set(c.variant_key for c in active if c.call.call != CALL_NO_CALL),
            key=natural_key),
        who_graded=who_graded,
        who_grade=who_grade,
        suppressed_by="; ".join(sorted(set(c.suppressed_by for c in suppressed))),
        level=level,
        cross_resistance=sorted(set(cross), key=natural_key),
        target_covered=target_covered,
        note=note,
    )


def _dedupe(items: Iterable[str]) -> List[str]:
    """Preserve order, drop repeats and empties."""
    seen = set()
    out: List[str] = []
    for item in items:
        if not item or item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out


def drugs_present(variants: Sequence[Variant]) -> List[str]:
    """Every drug any catalogue call in the sample names, in name order."""
    names = set()
    for variant in variants:
        for call in variant.catalogue_calls:
            canonical = normalise_drug(call.drug)
            if canonical:
                names.add(canonical)
    return sorted(names, key=natural_key)


def consensus(variants: Sequence[Variant], *,
              platform: str = PLATFORM_ILLUMINA,
              drugs: Optional[Sequence[str]] = None,
              target_covered: Optional[Mapping[str, Optional[bool]]] = None,
              suppressions: Optional[Sequence["rules.Suppression"]] = None
              ) -> List[DrugCall]:
    """A :class:`DrugCall` for every drug in the panel.

    *drugs* is the panel to report. It should be the drug set of the catalogue
    that was actually loaded — ``config.DRUGS`` is a display order, not an
    authority, and a third edition of the WHO catalogue must be able to change
    which drugs exist without this module being edited. Any drug a catalogue
    call names but the panel omits is appended rather than dropped, because
    silently discarding a graded row is the one failure mode that produces a
    confident-looking report with a resistance finding missing from it.

    *target_covered* maps drug to whether its target regions were callable.
    A drug absent from the mapping gets ``None``, which is "nobody established
    it" and is reported as such — not as coverage.

    *suppressions* defaults to running the epistasis rules over *variants*.
    Passing an explicit sequence is how a test pins one rule, and how the
    pipeline reuses one computation across the whole panel.
    """
    panel = [normalise_drug(d) for d in (drugs if drugs is not None else DRUGS)]
    for extra in drugs_present(variants):
        if extra not in panel:
            panel.append(extra)

    if suppressions is None:
        suppressions = rules.epistasis_suppressions(variants, drugs=panel)

    coverage: Mapping[str, Optional[bool]] = target_covered or {}
    normalised_coverage = dict(
        (normalise_drug(k), v) for k, v in coverage.items())

    return [
        consensus_for_drug(
            drug, variants, platform=platform, suppressions=suppressions,
            target_covered=normalised_coverage.get(drug))
        for drug in panel
    ]


def annex_rows(variants: Sequence[Variant], drug_call: DrugCall, *,
               platform: str = PLATFORM_ILLUMINA,
               suppressions: Sequence["rules.Suppression"] = ()) -> List[Dict[str, Any]]:
    """The side-by-side catalogue comparison for one drug, for the annex.

    Recomputed from the variants rather than stored on the :class:`DrugCall`,
    because the record deliberately keeps only the calls: an annex row is a view,
    and a view that is serialised alongside its source is a second copy waiting
    to disagree with it.
    """
    contributions = collect_contributions(variants, drug_call.drug)
    _apply_suppressions(contributions, drug_call.drug,
                        normalise_platform(platform), suppressions)
    return side_by_side(contributions)
