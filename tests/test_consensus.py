"""The three-catalogue consensus rule (§5.5), including the two paths that lie.

The whole engine exists to stop two specific sentences reaching a clinician:
"susceptible", when what happened is that nothing catalogued was found; and a
bare "R", when the call rests on a catalogue WHO does not grade. Both are tested
here directly, as are the suppressions — an epistasis or platform suppression
that quietly removed a determinant would produce the first sentence by another
route.
"""

from __future__ import annotations

import pytest

from conftest import make_variant, other_call, who_call
from mjolnir import config
from mjolnir.records import (CALL_NO_CALL, CALL_R, CALL_R_OUTSIDE_WHO,
                             CALL_S, CALL_UNCERTAIN, NO_DETERMINANT_TEXT)
from mjolnir.resistance import consensus, rules


# ------------------------------------------------------------- WHO is the anchor

def test_who_grade_is_the_call_even_when_another_catalogue_disagrees():
    """§5.5 rule 2. WHO is the only source with a published, derived grading."""
    variant = make_variant(gene="rpoB", hgvs="p.Leu430Arg", calls=[
        who_call("Rifampicin", config.WHO_GRADE_5, "rpoB_p.Leu430Arg"),
        other_call(config.CATALOGUE_MTBSEQ, "Rifampicin", CALL_R, "rpoB_p.Leu430Arg"),
    ])
    call = consensus.consensus_for_drug("Rifampicin", [variant])
    assert call.call == CALL_S
    assert call.who_graded is True
    assert call.who_grade == config.WHO_GRADE_5


def test_a_who_group_one_call_is_reported_as_r_with_its_grade():
    variant = make_variant(calls=[
        who_call("Rifampicin", config.WHO_GRADE_1, "rpoB_p.Ser450Leu",
                 comment="High-level resistance"),
    ])
    call = consensus.consensus_for_drug("Rifampicin", [variant])
    assert call.call == CALL_R
    assert call.confidence == "high"
    assert call.level == "high-level"
    assert call.who_grade == config.WHO_GRADE_1


def test_disagreement_between_catalogues_is_flagged_not_averaged():
    variant = make_variant(gene="katG", hgvs="p.Ser315Thr", calls=[
        who_call("Isoniazid", config.WHO_GRADE_5, "katG_p.Ser315Thr"),
        other_call(config.CATALOGUE_TBDB, "Isoniazid", CALL_R, "katG_p.Ser315Thr"),
    ])
    call = consensus.consensus_for_drug("Isoniazid", [variant])
    assert call.disagreement is True
    assert call.call == CALL_S, "WHO still anchors a flagged disagreement"
    assert call.confidence != "high", "a disagreement must cost confidence"


def test_a_catalogue_that_said_nothing_is_not_dissenting():
    """MTBseq's list is flat: silence from it is not a susceptible vote."""
    variant = make_variant(calls=[
        who_call("Rifampicin", config.WHO_GRADE_1, "rpoB_p.Ser450Leu"),
        other_call(config.CATALOGUE_MTBSEQ, "Rifampicin", CALL_NO_CALL,
                   "rpoB_p.Ser450Leu"),
    ])
    call = consensus.consensus_for_drug("Rifampicin", [variant])
    assert call.disagreement is False
    assert call.call == CALL_R


# ------------------------------------------------- §5.5 rule 3: outside the catalogue

def test_a_non_who_resistance_call_is_r_outside_who_never_plain_r():
    variant = make_variant(gene="Rv1258c", hgvs="p.Val219Ala", calls=[
        other_call(config.CATALOGUE_MTBSEQ, "Rifampicin", CALL_R, "Rv1258c_p.Val219Ala"),
    ])
    call = consensus.consensus_for_drug("Rifampicin", [variant])
    assert call.call == CALL_R_OUTSIDE_WHO
    assert call.call != CALL_R
    assert call.who_graded is False
    assert call.confidence == "low"


def test_the_outside_who_call_carries_the_sentence_that_bounds_it():
    variant = make_variant(gene="Rv1258c", hgvs="p.Val219Ala", calls=[
        other_call(config.CATALOGUE_TBDB, "Rifampicin", CALL_R, "Rv1258c_p.Val219Ala"),
    ])
    call = consensus.consensus_for_drug("Rifampicin", [variant])
    assert consensus.OUTSIDE_WHO_TEXT in call.caveats
    assert "not equivalent to a WHO Group 1 call" in consensus.OUTSIDE_WHO_TEXT


def test_the_mtbseq_asymmetry_is_stated_when_mtbseq_supplies_the_call():
    variant = make_variant(gene="Rv1258c", hgvs="p.Val219Ala", calls=[
        other_call(config.CATALOGUE_MTBSEQ, "Rifampicin", CALL_R, "Rv1258c_p.Val219Ala"),
    ])
    call = consensus.consensus_for_drug("Rifampicin", [variant])
    assert config.MTBSEQ_ASYMMETRY_NOTE in call.caveats


def test_the_label_of_an_outside_who_call_says_so():
    variant = make_variant(gene="Rv1258c", hgvs="p.Val219Ala", calls=[
        other_call(config.CATALOGUE_TBDB, "Rifampicin", CALL_R, "Rv1258c_p.Val219Ala"),
    ])
    call = consensus.consensus_for_drug("Rifampicin", [variant])
    assert call.label == "resistance determinant detected outside the WHO catalogue"


# ------------------------------------ §5.5 rule 5: absence is absence, not susceptibility

def test_nothing_catalogued_yields_no_determinant_detected():
    call = consensus.consensus_for_drug("Linezolid", [])
    assert call.call == CALL_NO_CALL
    assert call.label == NO_DETERMINANT_TEXT == "no resistance determinant detected"
    assert "susceptible" not in call.label.lower()


def test_the_no_determinant_note_says_it_is_not_a_phenotype():
    call = consensus.consensus_for_drug("Linezolid", [])
    assert "absence of evidence" in call.note
    assert "not a phenotypic result" in call.note


def test_no_call_is_never_rendered_as_susceptible_anywhere_in_the_panel():
    variants = [make_variant(calls=[who_call("Rifampicin", config.WHO_GRADE_1,
                                             "rpoB_p.Ser450Leu")])]
    panel = consensus.consensus(variants, drugs=config.DRUGS)
    for call in panel:
        if call.call == CALL_NO_CALL:
            assert "susceptible" not in call.label.lower()
            assert call.confidence == "none"


def test_uncovered_targets_are_not_evaluable_rather_than_negative():
    call = consensus.consensus_for_drug("Pyrazinamide", [], target_covered=False)
    assert call.label == "not evaluable: target regions not callable"
    assert call.confidence == "none"


def test_unknown_coverage_is_stated_when_nothing_was_found():
    call = consensus.consensus_for_drug("Pyrazinamide", [], target_covered=None)
    assert consensus.COVERAGE_UNKNOWN_TEXT in call.caveats


# ------------------------------------------------------- epistasis, end to end

def _bdq_sample():
    mmpl5 = make_variant(gene="mmpL5", hgvs="p.Ala100fs", pos=778990, ref="G",
                         alt="GA", effect="frameshift_variant",
                         variant_type="insertion")
    rv0678 = make_variant(gene="Rv0678", hgvs="p.Ser53Leu", pos=779000,
                          effect="missense_variant", calls=[
                              who_call("Bedaquiline", config.WHO_GRADE_1,
                                       "Rv0678_p.Ser53Leu"),
                              who_call("Clofazimine", config.WHO_GRADE_1,
                                       "Rv0678_p.Ser53Leu"),
                          ])
    return [mmpl5, rv0678]


@pytest.mark.parametrize("drug", ["Bedaquiline", "Clofazimine"])
def test_mmpl5_knockout_removes_the_rv0678_resistance_call(drug):
    calls = dict((c.drug, c) for c in consensus.consensus(
        _bdq_sample(), drugs=["Bedaquiline", "Clofazimine"]))
    call = calls[drug]
    assert call.call == consensus.SUPPRESSED_CALL == CALL_UNCERTAIN
    assert call.call != CALL_R


@pytest.mark.parametrize("drug", ["Bedaquiline", "Clofazimine"])
def test_the_suppression_is_recorded_rather_than_applied_silently(drug):
    """A suppressed determinant that vanished would read as one never found."""
    calls = dict((c.drug, c) for c in consensus.consensus(
        _bdq_sample(), drugs=["Bedaquiline", "Clofazimine"]))
    call = calls[drug]
    assert call.suppressed_by.startswith("epistasis:mmpL5")
    assert "suppressed" in call.note
    assert any("mmpL5" in caveat for caveat in call.caveats)
    # The abrogated catalogue row is still attached to the drug call.
    assert any(c.variant_key == "Rv0678_p.Ser53Leu" for c in call.catalogue_calls)


def test_a_suppressed_drug_is_not_reported_as_no_determinant_detected():
    calls = dict((c.drug, c) for c in consensus.consensus(
        _bdq_sample(), drugs=["Bedaquiline"]))
    assert calls["Bedaquiline"].label != NO_DETERMINANT_TEXT


def test_the_annex_shows_the_suppressed_row_flagged():
    variants = _bdq_sample()
    suppressions = rules.epistasis_suppressions(variants, drugs=["Bedaquiline"])
    call = consensus.consensus_for_drug("Bedaquiline", variants,
                                        suppressions=suppressions)
    rows = consensus.annex_rows(variants, call, suppressions=suppressions)
    suppressed = [r for r in rows if r["variant"] == "Rv0678_p.Ser53Leu"]
    assert suppressed and all(r["counted"] is False for r in suppressed)
    assert all(r["suppressed_by"] for r in suppressed)


def test_without_the_knockout_the_same_variant_is_resistant():
    """The control: the suppression, not the fixture, is doing the work."""
    _, rv0678 = _bdq_sample()
    calls = dict((c.drug, c) for c in consensus.consensus([rv0678],
                                                          drugs=["Bedaquiline"]))
    assert calls["Bedaquiline"].call == CALL_R


# ------------------------------------------------- ONT: the fbiC delamanid artefact

def _fbic_variant():
    return make_variant(gene="fbiC", hgvs="p.Gly100fs", pos=1303000, ref="GCCT",
                        alt="G", effect="frameshift_variant", variant_type="deletion",
                        allele_fraction=0.99, is_major=True, calls=[
                            who_call("Delamanid", config.WHO_GRADE_2, "fbiC_LoF"),
                        ])


def test_fbic_delamanid_is_suppressed_on_ont():
    """47.2% of all discordant drug classifications in the 508-isolate study."""
    call = consensus.consensus_for_drug("Delamanid", [_fbic_variant()],
                                        platform="ont")
    assert call.call != CALL_R
    assert call.call == consensus.SUPPRESSED_CALL
    assert call.suppressed_by == "platform:ont"
    assert config.ONT_FBIC_CAVEAT in call.caveats


def test_the_same_fbic_variant_is_called_on_illumina():
    call = consensus.consensus_for_drug("Delamanid", [_fbic_variant()],
                                        platform="illumina")
    assert call.call == "R-interim"
    assert call.suppressed_by == ""


def test_the_fbic_suppression_is_specific_to_that_gene_drug_pair():
    assert config.is_suppressed_on_platform("fbiC", "Delamanid", "ont")
    assert config.is_suppressed_on_platform("fbiC", "Pretomanid", "ont") is None
    assert config.is_suppressed_on_platform("ddn", "Delamanid", "ont") is None
    assert config.is_suppressed_on_platform("fbiC", "Delamanid", "illumina") is None


# --------------------------------------------------- ONT: the platform consequences

def test_ont_read_thresholds_are_higher_than_illumina():
    assert config.min_reads_for("ont") == 5
    assert config.min_reads_for("illumina") == 3
    assert config.MIN_READS_ONT > config.MIN_READS_ILLUMINA


def test_an_assembly_has_no_read_threshold_at_all():
    from mjolnir.utils import MjolnirError

    with pytest.raises(MjolnirError) as excinfo:
        config.min_reads_for("fasta")
    assert "no read evidence" in str(excinfo.value)


def test_an_ont_indel_call_carries_the_uncorroborated_caveat():
    variant = make_variant(gene="pncA", hgvs="p.Val130fs", ref="CG", alt="C",
                           effect="frameshift_variant", variant_type="deletion",
                           calls=[who_call("Pyrazinamide", config.WHO_GRADE_2,
                                           "pncA_LoF")])
    call = consensus.consensus_for_drug("Pyrazinamide", [variant], platform="ont")
    assert config.ONT_INDEL_CAVEAT in call.caveats
    assert "16.6%" in config.ONT_INDEL_CAVEAT


def test_an_ont_absence_says_absence_of_a_minor_variant_is_not_absence():
    call = consensus.consensus_for_drug("Linezolid", [], platform="ont")
    assert config.ONT_MINOR_VARIANT_CAVEAT in call.caveats


def test_an_assembly_states_its_capability_loss_on_every_drug():
    call = consensus.consensus_for_drug("Linezolid", [], platform="fasta")
    assert config.FASTA_CAPABILITY_LOSS in call.caveats
    assert "capability loss" in config.FASTA_CAPABILITY_LOSS


def test_illumina_adds_no_platform_caveat():
    call = consensus.consensus_for_drug("Linezolid", [], platform="illumina")
    assert config.ONT_MINOR_VARIANT_CAVEAT not in call.caveats
    assert config.FASTA_CAPABILITY_LOSS not in call.caveats


# ------------------------------------------------------------- the whole panel

def test_the_panel_covers_every_drug_and_appends_unlisted_ones():
    variant = make_variant(gene="rrl", hgvs="n.2814G>T", calls=[
        other_call(config.CATALOGUE_TBDB, "Clarithromycin", CALL_R, "rrl_n.2814G>T"),
    ])
    panel = consensus.consensus([variant], drugs=["Rifampicin"])
    names = [c.drug for c in panel]
    assert names[0] == "Rifampicin"
    assert "Clarithromycin" in names, "a graded row must never be silently dropped"


def test_drug_names_are_normalised_before_they_are_joined():
    """"rifampin" from one catalogue must meet "Rifampicin" from another."""
    variant = make_variant(calls=[
        who_call("rifampin", config.WHO_GRADE_1, "rpoB_p.Ser450Leu"),
    ])
    call = consensus.consensus_for_drug("RIF", [variant])
    assert call.drug == "Rifampicin"
    assert call.call == CALL_R
