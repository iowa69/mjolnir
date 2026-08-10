"""What may honestly be said about purity, and what may not (§8).

The measured facts this module is built on: Kraken2's sensitivity for
*M. tuberculosis* reads with a standard database is 0.0731 on real Illumina data
— about 93% of true target reads unclassified or misassigned — against ~0.97 with
a mycobacterial pangenome index. So a screen run against a standard or capped
index cannot support a statement about mycobacterial purity, and the failure mode
that matters is not a wrong number: it is a clean-looking report produced by an
instrument that could not have seen the problem.

Every test here is about that distinction. "We could not look" must never turn
into "we looked and it was clean", and an unassessed dimension must never fold
into a valid headline.
"""

from __future__ import annotations

import json

import pytest

from mjolnir import config
from mjolnir.contamination import purity
from mjolnir.records import (MIXTURE_NOT_ASSESSED, PLATFORM_FASTA,
                             VALIDITY_NOT_ASSESSED, VALIDITY_VALID, QCMetrics)
from mjolnir.utils import MjolnirError


def _index(tmp_path, name, pangenome=None):
    """A directory standing in for a Kraken2 index, optionally self-declaring."""
    path = tmp_path / name
    path.mkdir()
    (path / "hash.k2d").write_text("")
    if pangenome is not None:
        (path / "mjolnir_index.json").write_text(
            json.dumps({"mycobacterial_pangenome": pangenome}))
    return path


CLEAN_LOOKING_REPORT = (
    "99.84\t99840\t99840\tS\t1773\t  Mycobacterium tuberculosis\n"
    "0.16\t160\t160\tS\t1764\t  Mycobacterium avium\n"
)


# ------------------------------------------------------- the index is the gate

@pytest.mark.parametrize("name", [
    "k2_standard_20240112", "standard-8", "capped_16gb", "minikraken2_v2",
    "k2_pluspf_08gb", "core_nt_20250101",
])
def test_a_standard_or_capped_index_is_uninformative(tmp_path, name):
    screen = purity.evaluate_kraken2_screen(_index(tmp_path, name), confidence=0.1)
    assert screen.informative is False
    assert screen.status == purity.SCREEN_UNINFORMATIVE


def test_an_undeclared_index_is_uninformative_by_default(tmp_path):
    """Being wrong in this direction produces a false clean bill of health."""
    screen = purity.evaluate_kraken2_screen(_index(tmp_path, "someones_db"),
                                            confidence=0.1)
    assert screen.informative is False
    assert "mjolnir_index.json" in screen.note


def test_an_index_declaring_itself_not_a_pangenome_is_uninformative(tmp_path):
    screen = purity.evaluate_kraken2_screen(
        _index(tmp_path, "myco_ish", pangenome=False), confidence=0.1)
    assert screen.informative is False


def test_a_declared_mycobacterial_pangenome_index_is_informative(tmp_path):
    screen = purity.evaluate_kraken2_screen(
        _index(tmp_path, "myco_pangenome", pangenome=True), confidence=0.2)
    assert screen.informative is True
    assert screen.status == purity.SCREEN_INFORMATIVE


def test_a_missing_index_is_uninformative_rather_than_an_absence_of_contamination(tmp_path):
    screen = purity.evaluate_kraken2_screen(tmp_path / "not_here", confidence=0.1)
    assert screen.informative is False
    assert "not found" in screen.note


# ------------------------------- an uninformative screen prints nothing reassuring

def test_a_clean_looking_report_from_a_standard_index_is_not_reported(tmp_path):
    screen = purity.evaluate_kraken2_screen(
        _index(tmp_path, "k2_standard_20240112"), confidence=0.1,
        report_text=CLEAN_LOOKING_REPORT)
    assert screen.rows, "the rows are kept for support questions"
    assert screen.reportable_rows() == [], "but none of them may be printed"


def test_the_uninformative_screen_is_an_unmeasured_check_never_a_pass(tmp_path):
    screen = purity.evaluate_kraken2_screen(
        _index(tmp_path, "capped_8gb"), confidence=0.1,
        report_text=CLEAN_LOOKING_REPORT)
    check = screen.to_check()
    assert check.measured is False
    assert check.status == "warn"
    assert check.status != "pass"
    assert "0.0731" in check.reading
    assert "no statement about mycobacterial purity is made" in check.reading


def test_an_informative_screen_may_print_its_rows(tmp_path):
    screen = purity.evaluate_kraken2_screen(
        _index(tmp_path, "myco_pangenome", pangenome=True), confidence=0.2,
        report_text=CLEAN_LOOKING_REPORT)
    assert len(screen.reportable_rows()) == 2
    assert screen.to_check().status == "pass"


def test_a_run_with_no_screen_at_all_lands_in_the_same_place():
    screen = purity.no_screen()
    assert screen.informative is False
    assert screen.reportable_rows() == []
    assert screen.to_check().measured is False


# --------------------------------------------------------- the confidence refusal

def test_krakens_own_confidence_default_is_refused():
    """At 0.0 the screen prints a tail of low-abundance NTM as co-infections."""
    with pytest.raises(MjolnirError) as excinfo:
        config.kraken2_confidence(0.0)
    assert "refused" in str(excinfo.value)
    assert "co-infections" in str(excinfo.value)


def test_no_screen_can_be_evaluated_at_the_refused_default(tmp_path):
    with pytest.raises(MjolnirError):
        purity.evaluate_kraken2_screen(
            _index(tmp_path, "myco_pangenome", pangenome=True), confidence=0.0)


def test_an_unspecified_confidence_takes_mjolnirs_floor():
    assert config.kraken2_confidence(None) == config.KRAKEN2_MIN_CONFIDENCE
    assert config.KRAKEN2_MIN_CONFIDENCE > config.KRAKEN2_REFUSED_DEFAULT


def test_a_confidence_above_one_is_refused():
    with pytest.raises(MjolnirError):
        config.kraken2_confidence(1.5)


# --------------------------------------------------- tools that cannot answer this

@pytest.mark.parametrize("tool", ["CheckM", "checkm2", "check-m", "ConFindr"])
def test_marker_gene_and_rmlst_tools_are_refused_as_mixture_detectors(tool):
    with pytest.raises(MjolnirError) as excinfo:
        purity.assert_mixture_method_supported(tool)
    message = str(excinfo.value)
    assert "not a same-species mixture detector" in message or \
        "cannot be applied to this genus" in message
    assert "heterozygosity" in message, "the refusal must name the alternative"


def test_the_supported_method_is_not_refused():
    assert purity.assert_mixture_method_supported("mjolnir-heterozygosity") is None


# ------------------------------------------------------------- the verdict itself

def test_an_unassessed_dimension_never_folds_into_a_valid_headline():
    assert purity.worst_validity([VALIDITY_VALID, VALIDITY_NOT_ASSESSED]) == \
        VALIDITY_NOT_ASSESSED


def test_an_empty_set_of_verdicts_is_not_assessed():
    assert purity.worst_validity([]) == VALIDITY_NOT_ASSESSED


def test_an_assembly_can_never_be_reported_valid():
    """No allele fractions means the mixture dimension is permanently unassessed."""
    verdict = purity.sample_validity(platform=PLATFORM_FASTA,
                                     mixture_class=MIXTURE_NOT_ASSESSED)
    assert verdict.verdict == VALIDITY_NOT_ASSESSED
    assert config.FASTA_CAPABILITY_LOSS in verdict.caveats


def test_an_uninformative_screen_is_always_carried_into_the_caveats(tmp_path):
    screen = purity.evaluate_kraken2_screen(_index(tmp_path, "k2_standard"),
                                            confidence=0.1)
    verdict = purity.sample_validity(platform="illumina", screen=screen)
    assert any("mycobacterial purity" in caveat for caveat in verdict.caveats)


def test_the_verdict_is_per_intended_use_not_one_word():
    """99.84% pure produced 13 false-positive SNPs: enough to invent a cluster."""
    verdict = purity.sample_validity(
        platform="illumina", mixture_class="possible-mixture",
        mixture_reason="heterozygous-SNP fraction above the warn tier")
    assert verdict.by_use[purity.USE_RESISTANCE] == "suspect"
    assert verdict.by_use[purity.USE_TRANSMISSION] == "invalid"
    assert "resistance calling" in verdict.sentence()


def test_a_verdict_needs_a_question_to_answer():
    with pytest.raises(MjolnirError):
        purity.sample_validity(platform="illumina", intended_use=[])


def test_the_full_panel_reports_an_uninformative_screen_and_no_clean_result(tmp_path):
    screen = purity.evaluate_kraken2_screen(
        _index(tmp_path, "k2_standard_16gb"), confidence=0.1,
        report_text=CLEAN_LOOKING_REPORT)
    result = purity.assess_contamination(
        platform="illumina",
        qc=QCMetrics(mean_depth=60.0, mapped_fraction=0.98, breadth_min_depth=0.99),
        screen=screen)
    assert result.screen_informative is False
    assert result.verdict != VALIDITY_VALID, \
        "no mixture evidence was measured, so 'valid' would be unmeasured"
    assert result.mixture_class == MIXTURE_NOT_ASSESSED
    named = dict((c.name, c) for c in result.checks)
    assert named["taxonomic_contamination_screen"].measured is False


def test_the_evidence_sentence_behind_the_verdict_is_registered():
    assert config.source_for("contamination_evidence") == config.SRC_DESIGN
    assert "3,325" in config.CONTAMINATION_EVIDENCE
    assert "coarse instrument" in config.CONTAMINATION_EVIDENCE


# --------------------------------------------------------------- report parsing

def test_a_wholly_malformed_report_raises_rather_than_estimating_over_nothing():
    with pytest.raises(MjolnirError) as excinfo:
        purity.parse_kraken2_report("not\ta\treport\n")
    assert "six tab-separated columns" in str(excinfo.value)


def test_an_empty_report_is_empty_not_an_error():
    assert purity.parse_kraken2_report("") == []
