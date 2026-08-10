"""The discipline rules: what the model is allowed to have written.

The verdict is a rule and the prose is the model, so the only question here is
whether an answer contradicts, embellishes or invents. Two failure modes are
symmetric and both are tested: a guard that lets a fabricated number through is
useless, and a guard that rejects a correct concessive sentence is worse than
useless — it discards good prose and teaches whoever tunes it to loosen the rule
that was working.

The concessive exemption is not a nicety. In the tesseract-ai corpus it is what
saved 86 of 96 otherwise-correct readings, because "while X was fine, Y was not
measured" is how careful people actually write.
"""

from __future__ import annotations

import json

import pytest

from mjolnir.agent import discipline
from mjolnir.records import (CALL_NO_CALL, CALL_R, DrugCall, SampleResult,
                             SpeciesCall, VALIDITY_SUSPECT, VALIDITY_NOT_ASSESSED)

CONTEXT = (
    "sample: 226-18 (illumina)\n"
    "mean_depth = 42.0 x (threshold 25 x, Hall et al. 2024)\n"
    "mapped_read_fraction = 0.97\n"
    "breadth_min_depth: not measured\n"
    "taxonomic_contamination_screen: not measured - the index is not a "
    "mycobacterial pangenome database\n"
    "species: Mycobacterium tuberculosis complex (not resolved below complex)\n"
    "Rifampicin: R (WHO group 1), rpoB p.Ser450Leu\n"
    "Isoniazid: no resistance determinant detected\n"
    "sample validity: suspect\n"
)

FACTS = discipline.RuleFacts(
    unmeasured=["breadth_min_depth", "taxonomic_contamination_screen"],
    resistant_drugs=["Rifampicin"],
    no_determinant_drugs=["Isoniazid"],
    validity=VALIDITY_SUSPECT,
    platform="illumina",
)


# ------------------------------------------------------------ invented numbers

def test_an_answer_inventing_a_number_is_rejected():
    verdict = discipline.review(
        "Mean depth reached 96.4x, comfortably above the threshold.", CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "unsupported-number"
    assert "96.4" in verdict.reason


def test_a_percentage_the_model_computed_itself_is_still_invented():
    """The observation carries the percentages the report speaks in, for this reason."""
    verdict = discipline.review(
        "Only 3.0% of reads failed to map.", CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "unsupported-number"


def test_a_number_present_in_the_evidence_is_accepted():
    verdict = discipline.review(
        "Mean depth was 42.0x, above the 25x threshold.", CONTEXT, FACTS)
    assert verdict


def test_a_rounding_of_a_supplied_number_is_accepted():
    verdict = discipline.review("Mapping was 0.97 of reads.", CONTEXT, FACTS)
    assert verdict


def test_a_number_that_gains_digits_is_not_a_rounding():
    """18 in the evidence must not license 180 in the answer."""
    verdict = discipline.review("Depth was 420x.", CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "unsupported-number"


def test_small_ordinals_are_allowed_without_support():
    verdict = discipline.review(
        "Only 1 of the 3 catalogues graded this variant.", CONTEXT, FACTS)
    assert verdict


def test_numbers_in_json_scaffolding_are_not_claims():
    """A confidence the system demanded is not an assertion about the evidence."""
    answer = json.dumps({"choice": "warn", "confidence": 0.82,
                         "headline": "Depth was 42.0x."})
    assert discipline.review(answer, CONTEXT, FACTS)


# ---------------------------------------------------- the concessive exemption

def test_a_concessive_clause_is_not_falsely_rejected():
    answer = ("While mapping and depth are adequate, coverage breadth was not "
              "measured, so genome-wide statements are not supported here.")
    verdict = discipline.review(answer, CONTEXT, FACTS)
    assert verdict, verdict.reason


@pytest.mark.parametrize("answer", [
    "Although the rifampicin call is firm, the contamination screen was not informative.",
    "The rifampicin determinant is present; however, breadth was never measured.",
    "Aside from the resistance finding, nothing else was established.",
    "Isoniazid showed no resistance determinant, which is not the same as susceptibility.",
    "Absence of an isoniazid determinant is an absence of evidence, not a susceptible result.",
])
def test_correct_careful_prose_survives(answer):
    verdict = discipline.review(answer, CONTEXT, FACTS)
    assert verdict, "{0!r} rejected: {1}".format(answer, verdict.reason)


def test_a_hedged_sentence_about_an_unmeasured_check_passes():
    answer = ("The contamination screen was not informative, so purity is unknown.")
    assert discipline.review(answer, CONTEXT, FACTS)


# ------------------------------------------------- unmeasured described as fine

def test_calling_an_unmeasured_check_clean_is_rejected():
    verdict = discipline.review(
        "The contamination screen came back clean.", CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "unmeasured-as-settled"
    assert "taxonomic_contamination_screen" in verdict.reason


def test_no_contamination_detected_is_rejected_when_nothing_could_detect_it():
    verdict = discipline.review(
        "No contamination was detected in this library.", CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "unmeasured-as-settled"


# ------------------------------------------------- susceptibility for a no-call

def test_calling_a_no_determinant_drug_susceptible_is_rejected():
    verdict = discipline.review(
        "The isolate is susceptible to isoniazid.", CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "susceptible-for-no-determinant"


def test_an_unqualified_susceptibility_claim_is_rejected_when_no_drug_is_named():
    """With no drug list to check against, any susceptibility predicate falls."""
    verdict = discipline.review("The isolate appears to be susceptible.",
                                CONTEXT, discipline.RuleFacts())
    assert not verdict
    assert verdict.rule == "susceptible-for-no-determinant"


@pytest.mark.xfail(strict=True, reason=(
    "discipline.claims_susceptibility only checks the drug-free form when the "
    "fact set names no drugs: `if not names or any(mentions_drug(...))`. With a "
    "no-determinant drug list present, an unqualified 'the isolate is "
    "susceptible' names no drug and passes, which is the exact sentence §5.5 "
    "rule 5 forbids. Its own docstring says both forms are caught. Reported to "
    "the owner of agent/discipline.py; remove this xfail when it is fixed."))
def test_an_unqualified_susceptibility_claim_is_rejected_alongside_a_drug_list():
    verdict = discipline.review("The isolate appears to be susceptible.",
                                CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "susceptible-for-no-determinant"


def test_recommending_phenotypic_testing_is_not_a_susceptibility_claim():
    answer = "Phenotypic susceptibility testing is required to settle isoniazid."
    assert discipline.review(answer, CONTEXT, FACTS)


# -------------------------------------------------------- contradicting a call

def test_clearing_a_drug_the_rules_called_resistant_is_rejected():
    verdict = discipline.review(
        "Rifampicin remains an option for this patient.", CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "contradicts-drug-call"
    assert "Rifampicin" in verdict.reason


def test_asserting_resistance_the_rules_did_not_call_is_rejected():
    verdict = discipline.review(
        "The isolate is resistant to isoniazid as well.", CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "contradicts-drug-call"


def test_qualifying_a_call_without_reversing_it_is_allowed():
    answer = ("The rifampicin call rests on a single catalogued substitution, so "
              "it should be read with the platform caveats.")
    assert discipline.review(answer, CONTEXT, FACTS)


# -------------------------------------------------------------- the whole sample

def test_declaring_a_suspect_sample_clean_is_rejected():
    verdict = discipline.review("The sample is clean.", CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "contradicts-validity"


def test_no_problems_anywhere_is_rejected_on_a_suspect_sample():
    verdict = discipline.review("There were no major issues with this run.",
                                CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "contradicts-validity"


def test_an_unassessed_sample_may_not_be_endorsed():
    facts = discipline.RuleFacts(validity=VALIDITY_NOT_ASSESSED)
    verdict = discipline.review("The sample is fit for purpose.", CONTEXT, facts)
    assert not verdict
    assert verdict.rule == "contradicts-validity"


# ---------------------------------------------------------- invented organisms

def test_naming_an_organism_the_evidence_never_contained_is_rejected():
    verdict = discipline.review(
        "This is Mycobacterium bovis, so pyrazinamide is intrinsically inactive.",
        CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "invented-organism"


def test_repeating_the_organism_the_evidence_gave_is_allowed():
    assert discipline.review(
        "The isolate is a member of the Mycobacterium tuberculosis complex.",
        CONTEXT, FACTS)


# ------------------------------------------------------------ the closed set

def test_a_gate_answer_outside_the_closed_set_is_rejected():
    answer = json.dumps({"choice": "probably fine", "headline": "All good."})
    verdict = discipline.review(answer, CONTEXT, discipline.RuleFacts(),
                                options=["pass", "warn", "fail"])
    assert not verdict
    assert verdict.rule == "outside-option-set"


def test_a_model_may_not_override_the_rule_derived_choice():
    answer = json.dumps({"choice": "pass", "headline": "Depth was 42.0x."})
    verdict = discipline.review(answer, CONTEXT, discipline.RuleFacts(),
                                options=["pass", "warn", "fail"],
                                expected_choice="warn")
    assert not verdict
    assert verdict.rule == "overrides-rule"


def test_an_empty_answer_is_rejected():
    verdict = discipline.review("   ", CONTEXT, FACTS)
    assert not verdict
    assert verdict.rule == "empty"


# ------------------------------------------- the facts come from the result itself

def test_facts_are_derived_from_the_sample_result_not_passed_by_hand():
    """A drug call added anywhere in the pipeline is policed without rewiring."""
    result = SampleResult(sample_id="226-18", platform="illumina")
    result.species = SpeciesCall(name="unresolved", complex="MTBC")
    result.drugs = [
        DrugCall(drug="Rifampicin", call=CALL_R),
        DrugCall(drug="Isoniazid", call=CALL_NO_CALL),
        DrugCall(drug="Pyrazinamide", call=CALL_NO_CALL, target_covered=False),
    ]
    facts = discipline.facts_from_sample(result)
    assert facts.resistant_drugs == ["Rifampicin"]
    assert facts.no_determinant_drugs == ["Isoniazid"]
    assert facts.unevaluable_drugs == ["Pyrazinamide"]
    assert "Pyrazinamide" not in facts.no_determinant_drugs, \
        "an unevaluable drug is not a drug with nothing found"


def test_an_unevaluable_drug_may_not_be_called_susceptible_either():
    result = SampleResult(sample_id="226-18", platform="illumina")
    result.drugs = [DrugCall(drug="Pyrazinamide", call=CALL_NO_CALL,
                             target_covered=False)]
    facts = discipline.facts_from_sample(result)
    verdict = discipline.review("Pyrazinamide is susceptible.", CONTEXT, facts)
    assert not verdict
    assert verdict.rule == "susceptible-for-no-determinant"
