"""Thresholds carry their sources, and absence is a value rather than a default.

House rule 1 is that no number appears in logic without a document behind it, and
the mechanism is ``config._define`` — a constant that skipped it would print in
the PDF with no attribution at all. House rule 5 is that nothing unmeasured is
described as fine, and the mechanism is the record layer: ``DrugCall.call``
defaults to ``no-call`` rather than ``S``, ``Check.not_measured`` warns rather
than passes, and every measurement that can fail to exist is ``Optional`` and
``None``.

Both are tested from the outside here, so that a later edit that reintroduces a
cheerful default has to break something.
"""

from __future__ import annotations

import pytest

from mjolnir import config
from mjolnir.records import (CALL_LABELS, CALL_NO_CALL, CALL_R, CALL_S,
                             CONFIDENCE_NONE, MIXTURE_NOT_ASSESSED,
                             NO_DETERMINANT_TEXT, RESISTANCE_CALLS,
                             VALIDITY_NOT_ASSESSED, Check, ContaminationResult,
                             DrugCall, LineageCall, QCMetrics, SampleResult,
                             SpeciesCall, worst_call, worst_status)
from mjolnir.utils import MjolnirError


# ---------------------------------------------------------------- provenance

def test_every_threshold_names_a_source():
    for threshold in config.all_thresholds():
        assert threshold.source.strip(), "{0} has no source".format(threshold.name)
        assert len(threshold.source) > 10, \
            "{0}'s source is not a citation: {1!r}".format(threshold.name,
                                                           threshold.source)


def test_every_threshold_describes_itself_for_the_report():
    for threshold in config.all_thresholds():
        text = threshold.describe()
        assert threshold.name in text
        assert threshold.source in text


def test_an_unverified_citation_is_marked_rather_than_hidden():
    """An unverified citation is worse than none: it looks settled."""
    unverified = config.unverified()
    assert unverified, "the design records several citations as unchecked"
    for threshold in unverified:
        assert "[citation unverified]" in threshold.describe()


def test_an_unregistered_threshold_is_loud_rather_than_empty():
    assert "unregistered" in config.source_for("not_a_real_threshold")
    with pytest.raises(MjolnirError) as excinfo:
        config.threshold("not_a_real_threshold")
    assert "config.py" in str(excinfo.value)


@pytest.mark.parametrize("name", [
    "min_reads_illumina", "min_reads_ont", "min_depth", "major_variant_fraction",
    "cluster_snp_strict", "cluster_snp_relaxed", "snp_proximity_window",
    "min_shared_callable_sites", "who_xlsx_header_row", "rpob_rrdr_codons",
    "rpob_borderline", "epistasis_rules", "kraken2_mtb_sensitivity_standard",
    "ani_species_floor", "erm41_sequevar_position",
])
def test_the_numbers_the_report_prints_are_all_registered(name):
    assert config.threshold(name).source


def test_a_threshold_cannot_be_defined_twice():
    with pytest.raises(MjolnirError):
        config._define("min_depth", 25, "a second definition")


def test_a_cluster_threshold_prints_its_basis():
    assert "Walker" in config.cluster_threshold_basis(5)
    assert "Walker" in config.cluster_threshold_basis(12)
    assert "prior local MTBseq" in config.cluster_threshold_basis(6)
    assert "no published basis" in config.cluster_threshold_basis(7)


def test_the_local_chimaera_distance_is_recorded_as_precedent_not_as_a_default():
    assert config.CHIMAERA_LOCAL_DISTANCE == 6
    assert config.DEFAULT_CLUSTER_DISTANCE == config.CLUSTER_SNP_RELAXED == 12


# --------------------------------------------------------------- the platforms

def test_platform_aliases_resolve_and_nonsense_raises():
    from mjolnir.records import normalise_platform

    assert normalise_platform("nanopore") == "ont"
    assert normalise_platform("PE") == "illumina"
    assert normalise_platform("assembly") == "fasta"
    with pytest.raises(MjolnirError):
        normalise_platform("pacbio")


def test_a_missing_allele_fraction_does_not_default_to_major():
    assert config.is_major_variant(None) is None
    assert config.is_major_variant(0.95) is True
    assert config.is_major_variant(0.5) is False


def test_ont_carries_all_three_platform_consequences():
    caveats = config.platform_caveats("ont")
    assert config.ONT_MINOR_VARIANT_CAVEAT in caveats
    assert config.ONT_FBIC_CAVEAT in caveats
    assert config.ONT_INDEL_CAVEAT in caveats


def test_a_fasta_states_its_capability_loss():
    assert config.platform_caveats("fasta") == (config.FASTA_CAPABILITY_LOSS,)


# ------------------------------------------------------------ absence is a value

def test_a_drug_call_defaults_to_no_call_and_never_to_susceptible():
    call = DrugCall(drug="Linezolid")
    assert call.call == CALL_NO_CALL
    assert call.confidence == CONFIDENCE_NONE
    assert call.label == NO_DETERMINANT_TEXT
    assert "susceptible" not in call.label.lower()


def test_no_call_wording_is_spelled_once_and_is_not_susceptibility():
    assert CALL_LABELS[CALL_NO_CALL] == "no resistance determinant detected"
    assert "susceptible" not in CALL_LABELS[CALL_NO_CALL]
    assert CALL_LABELS[CALL_S] == "variant graded not associated with resistance"


def test_a_call_outside_the_closed_set_raises():
    with pytest.raises(MjolnirError):
        DrugCall(drug="Linezolid", call="probably fine")
    assert "R-outside-WHO" in RESISTANCE_CALLS


def test_an_unexamined_drug_outranks_a_graded_susceptible_one():
    """no-call is more alarming than S, and the ordering has to say so."""
    assert worst_call([CALL_S, CALL_NO_CALL]) == CALL_NO_CALL
    assert worst_call([CALL_NO_CALL, CALL_R]) == CALL_R
    assert worst_call([]) == CALL_NO_CALL


def test_an_unmeasured_check_warns_rather_than_passing():
    check = Check.not_measured("breadth", "no reads were mapped")
    assert check.measured is False
    assert check.status == "warn"
    assert check.ok is False


def test_a_numeric_check_with_no_value_becomes_an_unmeasured_check():
    check = Check.numeric("mean_depth", None, warn_minimum=25.0,
                          not_measured_why="no BAM was produced")
    assert check.measured is False
    assert check.reading == "no BAM was produced"


def test_a_boolean_check_with_no_value_is_not_established():
    check = Check.boolean("reference_matched", None)
    assert check.measured is False


def test_worst_status_folds_to_the_most_severe():
    assert worst_status(["pass", "warn", "fail"]) == "fail"
    assert worst_status([]) == "pass"


def test_a_fresh_sample_result_asserts_nothing():
    result = SampleResult(sample_id="226-18")
    assert result.species.resolved_to_species is False
    assert result.species.display == "unresolved"
    assert result.lineage.display == "not determined"
    assert result.contamination.verdict == VALIDITY_NOT_ASSESSED
    assert result.contamination.mixture_class == MIXTURE_NOT_ASSESSED
    assert result.contamination.screen_informative is False
    assert result.qc.mean_depth is None


def test_qc_metrics_are_none_rather_than_zero_when_absent():
    qc = QCMetrics()
    for field in ("mean_depth", "breadth_1x", "mapped_fraction", "gc_content",
                  "unambiguous_fraction", "coverage_evenness"):
        assert getattr(qc, field) is None, field


def test_a_lineage_support_fraction_is_none_when_nothing_was_callable():
    lineage = LineageCall(lineage="lineage2", barcode_sites_total=1111)
    assert lineage.support_fraction is None
    assert lineage.callable_fraction == 0.0


def test_unmeasured_checks_are_listed_for_the_agent():
    result = SampleResult(sample_id="226-18")
    result.checks = [Check.not_measured("f2", "no lineage sites were genotyped"),
                     Check(name="mean_depth", value=60.0, status="pass")]
    assert result.unmeasured() == ["f2"]


def test_an_unknown_validity_verdict_cannot_be_constructed():
    with pytest.raises(MjolnirError):
        ContaminationResult(verdict="probably ok")


def test_a_species_call_needs_a_known_confidence():
    with pytest.raises(MjolnirError):
        SpeciesCall(name="Mycobacterium chimaera", confidence="quite sure")


# --------------------------------------------------------- the run configuration

def test_a_changed_threshold_is_remembered_as_the_operators_choice():
    cfg = config.Config()
    cfg.set_explicit("min_depth", 8)
    cfg.set_explicit("degraded_depth_floor", 5)
    assert cfg.was_set("min_depth")
    assert cfg.overridden_thresholds() == {"degraded_depth_floor": 5, "min_depth": 8}


def test_a_target_depth_below_the_floor_is_refused():
    cfg = config.Config(min_depth=5, degraded_depth_floor=10)
    with pytest.raises(MjolnirError) as excinfo:
        cfg.validate()
    assert "cannot be lower than the floor" in str(excinfo.value)


def test_the_run_configuration_refuses_krakens_zero_confidence():
    cfg = config.Config(kraken2_confidence=0.0)
    with pytest.raises(MjolnirError):
        cfg.validate()


def test_an_unknown_profile_is_refused_before_the_run():
    cfg = config.Config(profile="pretty")
    with pytest.raises(MjolnirError):
        cfg.validate()


def test_a_valid_default_configuration_validates():
    cfg = config.Config()
    cfg.validate()
    assert cfg.profile in config.PROFILES
    assert cfg.use_llm is True, "the model is on by default and gates nothing"


# ------------------------------------------------------------- the NTM evidence

def test_an_ntm_pair_with_no_evidence_base_is_answered_with_absence():
    assert config.ntm_targets("Mycobacterium chimaera", "Clarithromycin")
    assert config.ntm_targets("Mycobacterium chimaera", "Linezolid") is None
    assert config.ntm_targets("Mycobacterium gordonae", "Amikacin") is None
    assert "not a prediction of susceptibility" in config.NTM_NO_EVIDENCE_TEXT


def test_the_abscessus_genus_split_does_not_lose_the_evidence():
    assert config.ntm_targets("Mycobacterium abscessus", "Clarithromycin") == \
        config.ntm_targets("Mycobacteroides abscessus", "Clarithromycin")
