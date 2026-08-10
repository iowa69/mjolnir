"""The report's data layer: the last place a call can quietly change meaning.

Everything the PDF, the HTML and the TSV/JSON artefacts print is built from one
projection of :class:`~mjolnir.records.SampleResult`, precisely so that three
outputs cannot come to disagree about one run. These tests hold that projection
to the same rules the engine obeys: a ``no-call`` renders as "no resistance
determinant detected", an ``R-outside-WHO`` never renders as ``R``, a distance
never appears without its denominator, and an unmeasured check never renders as
a pass.

The JSON payload is deterministic by construction — ``generated`` defaults to
empty — so the golden-file comparison is a comparison of content, not of clocks.
"""

from __future__ import annotations

import json

import pytest

from conftest import who_call
from mjolnir import config
from mjolnir.records import (CALL_NO_CALL, CALL_R, CALL_R_OUTSIDE_WHO, Check,
                             Cluster, CohortResult, ContaminationResult,
                             DatabaseVersion, DrugCall, Interpretation,
                             LineageCall, NO_DETERMINANT_TEXT, PairwiseDistance,
                             QCMetrics, SampleResult, SpeciesCall, Variant)
from mjolnir.report import html, tables
from mjolnir.utils import MjolnirError


@pytest.fixture
def result():
    """One sample carrying every shape the front page has to render."""
    variant = Variant(chrom="NC_000962.3", pos=761155, ref="C", alt="T",
                      gene="rpoB", hgvs="p.Ser450Leu", depth=64,
                      allele_fraction=0.98, is_major=True,
                      effect="missense_variant")
    variant.catalogue_calls.append(
        who_call("Rifampicin", config.WHO_GRADE_1, "rpoB_p.Ser450Leu"))
    sample = SampleResult(
        sample_id="226-18", platform="illumina", reference="NC_000962.3",
        variants=[variant],
        species=SpeciesCall(name="unresolved", complex="MTBC", method="skani ANI",
                            ani=99.8, confidence="moderate"),
        lineage=LineageCall(lineage="lineage4", sublineage="lineage4.9",
                            barcode_sites_supporting=18, barcode_sites_callable=18,
                            barcode_sites_total=23, confidence="high",
                            scheme="tbdb barcode.bed"),
        qc=QCMetrics(mean_depth=64.0, breadth_min_depth=0.99, mapped_fraction=0.98,
                     reference="NC_000962.3"),
        contamination=ContaminationResult(verdict="valid", mixture_class="single-strain",
                                          screen_informative=False,
                                          screen_note="no index configured"),
        database_versions=[DatabaseVersion(name=config.CATALOGUE_WHO,
                                           version="WHO-UCN-TB-2023.7",
                                           checksum="abc123", licence="ODC-By v1.0")],
        tool_versions={"bwa-mem2": "2.2.1"},
    )
    sample.drugs = [
        DrugCall(drug="Rifampicin", call=CALL_R, confidence="high", who_graded=True,
                 who_grade=config.WHO_GRADE_1,
                 supporting_variants=["rpoB_p.Ser450Leu"],
                 catalogue_calls=list(variant.catalogue_calls)),
        DrugCall(drug="Isoniazid", call=CALL_NO_CALL),
        DrugCall(drug="Bedaquiline", call=CALL_R_OUTSIDE_WHO, confidence="low",
                 who_graded=False, supporting_variants=["Rv1258c_p.Val219Ala"]),
        DrugCall(drug="Pyrazinamide", call=CALL_NO_CALL, target_covered=False),
    ]
    sample.checks = [
        Check(name="mean_depth", value=64.0, threshold=25, status="pass",
              comparison=">=", unit="x", category="qc"),
        Check.not_measured("taxonomic_contamination_screen",
                           "the index is not a mycobacterial pangenome database",
                           category="contamination"),
    ]
    return sample


@pytest.fixture
def cohort():
    return CohortResult(
        samples=["226-18", "30-20", "41-19"],
        pairs=[
            PairwiseDistance("226-18", "30-20", snps=3,
                             shared_callable_sites=4_100_000, masked_sites=264_525),
            PairwiseDistance("226-18", "41-19", snps=None,
                             shared_callable_sites=None,
                             note="not computed: no shared callable denominator"),
        ],
        clusters=[Cluster(cluster_id="C1", members=["226-18", "30-20"], threshold=12,
                          max_distance=3, min_shared_callable_sites=4_100_000)],
        threshold=12, threshold_basis=config.cluster_threshold_basis(12),
        mask_name="tbdb mask.bed", masked_sites=264_525, masked_fraction=0.06,
        reference="NC_000962.3")


# ------------------------------------------------------- the drug table's wording

def test_a_no_call_renders_as_no_determinant_never_as_susceptible(result):
    rows = dict((row["drug"], row) for row in tables.drug_rows(result))
    assert rows["Isoniazid"]["call"] == CALL_NO_CALL
    assert rows["Isoniazid"]["call_label"] == NO_DETERMINANT_TEXT
    for row in tables.drug_rows(result):
        assert "susceptible" not in str(row["call_label"]).lower()


def test_an_outside_who_call_is_distinguishable_from_a_who_group_one(result):
    rows = dict((row["drug"], row) for row in tables.drug_rows(result))
    assert rows["Bedaquiline"]["call"] == CALL_R_OUTSIDE_WHO
    assert rows["Bedaquiline"]["call"] != rows["Rifampicin"]["call"]
    assert rows["Bedaquiline"]["who_graded"] is False
    assert rows["Rifampicin"]["who_grade"] == config.WHO_GRADE_1
    assert rows["Bedaquiline"]["call_glyph"] != rows["Rifampicin"]["call_glyph"]


def test_an_uncovered_target_renders_as_not_evaluable(result):
    rows = dict((row["drug"], row) for row in tables.drug_rows(result))
    assert rows["Pyrazinamide"]["target_covered"] is False
    assert "not evaluable" in rows["Pyrazinamide"]["call_label"]


def test_the_catalogue_evidence_carries_version_and_checksum(result):
    rows = tables.catalogue_call_rows(result)
    assert rows
    assert rows[0]["catalogue_version"] == "WHO-UCN-TB-2023.7"
    assert rows[0]["matched_by"] == "coordinate"


# ------------------------------------------------------------ the headline

def test_the_headline_is_rule_derived_when_there_is_no_model(result):
    headline, provenance = tables.headline_sentence(result)
    assert "Rifampicin" in headline
    assert "not equivalent to a WHO Group 1 call" in headline
    assert provenance.startswith("rule-derived")


def test_a_discarded_model_answer_is_named_in_the_provenance(result):
    result.interpretation = Interpretation(
        headline="", rule_only=True,
        discarded_reason="states numbers that are not in the evidence: 96.4")
    _headline, provenance = tables.headline_sentence(result)
    assert "discarded" in provenance
    assert "96.4" in provenance


def test_a_kept_model_answer_names_the_model(result):
    result.interpretation = Interpretation(headline="Rifampicin resistance is predicted.",
                                           rule_only=False, model="qwen3:32b",
                                           host="http://localhost:11434")
    headline, provenance = tables.headline_sentence(result)
    assert headline == "Rifampicin resistance is predicted."
    assert "qwen3:32b" in provenance


def test_a_sample_with_no_findings_is_not_written_as_reassurance():
    sample = SampleResult(sample_id="30-20", platform="illumina")
    sample.drugs = [DrugCall(drug=name) for name in ("Rifampicin", "Isoniazid")]
    headline = tables.rule_headline(sample)
    assert "not evidence of susceptibility" in headline
    assert "2 drugs evaluated" in headline


def test_unmeasured_checks_are_named_in_the_headline(result):
    assert "could not be measured" in tables.rule_headline(result)


# --------------------------------------------------- distances keep their denominator

def test_every_distance_row_carries_its_shared_callable_sites(cohort):
    rows = tables.distance_rows(cohort)
    for row in rows:
        assert "shared_callable_sites" in row
    computed = [r for r in rows if r["snps"] not in (None, tables.TSV_NA)]
    assert computed and all(r["shared_callable_sites"] for r in computed)


def test_an_uncomputed_pair_is_na_rather_than_zero(cohort):
    rows = dict(((r["sample_a"], r["sample_b"]), r) for r in tables.distance_rows(cohort))
    uncomputed = rows[("226-18", "41-19")]
    assert uncomputed["snps"] in (None, tables.TSV_NA)
    assert uncomputed["snps"] != 0


def test_the_matrix_writes_na_where_a_pair_was_never_compared(cohort, tmp_path):
    tables.write_cohort_tables(tmp_path, cohort)
    text = (tmp_path / "cohort.matrix.tsv").read_text()
    assert tables.TSV_NA in text
    header = text.splitlines()[0].split("\t")
    assert header == ["sample"] + cohort.samples


def test_the_cohort_methods_block_names_the_mask_and_the_threshold(cohort):
    pairs = dict(tables.cohort_pairs(cohort))
    joined = " ".join("{0}: {1}".format(k, v) for k, v in pairs.items())
    assert "tbdb mask.bed" in joined
    assert "Walker" in joined


# ----------------------------------------------------------------- the artefacts

def test_the_json_payload_is_deterministic(result):
    assert tables.sample_json(result) == tables.sample_json(result)
    assert tables.sample_json(result)["report"]["generated"] == ""


def test_the_json_payload_repeats_the_no_determinant_wording(result):
    payload = tables.sample_json(result)
    isoniazid = [d for d in payload["drugs"] if d["drug"] == "Isoniazid"][0]
    assert isoniazid["call"] == CALL_NO_CALL
    assert "susceptible" not in json.dumps(payload).lower()


def test_the_json_payload_lists_what_could_not_be_measured(result):
    payload = tables.sample_json(result)
    assert "taxonomic_contamination_screen" in payload["unmeasured"]


def test_writing_the_sample_tables_produces_every_view(result, tmp_path):
    written = tables.write_sample_tables(tmp_path, result)
    names = sorted(p.name for p in written)
    assert "226-18.json" in names
    assert "226-18.drugs.tsv" in names
    assert "226-18.checks.tsv" in names
    for path in written:
        assert path.exists() and path.read_text().strip()


def test_an_empty_table_still_writes_its_header(tmp_path):
    """A header-only file says the analysis ran and found nothing."""
    path = tables.write_tsv(tmp_path / "empty.tsv", [], columns=["a", "b"])
    assert path.read_text() == "a\tb\n"


def test_writing_nothing_at_all_raises(tmp_path):
    with pytest.raises(MjolnirError) as excinfo:
        tables.write_tables(tmp_path, [], None)
    assert "empty result set" in str(excinfo.value)


def test_the_thresholds_artefact_carries_every_source(tmp_path, result):
    tables.write_tables(tmp_path, [result])
    rows = (tmp_path / "thresholds.tsv").read_text().splitlines()
    assert len(rows) > 20
    assert "source" in rows[0]


def test_the_check_rows_never_render_an_unmeasured_check_as_a_pass(result):
    rows = tables.check_rows(tables.all_checks(result), result.sample_id)
    unmeasured = [r for r in rows if r.get("measured") is False]
    assert unmeasured
    assert all(r["status"] != "pass" for r in unmeasured)


# -------------------------------------------------------------- the HTML mirror

def test_the_html_is_self_contained_and_carries_the_same_wording(result):
    document = html.render_html(result)
    assert document.lstrip().lower().startswith("<!doctype html>")
    assert "<script src=" not in document, "no script may be fetched at view time"
    assert "<link rel=\"stylesheet\"" not in document
    assert NO_DETERMINANT_TEXT in document
    assert "226-18" in document


def test_the_html_states_the_outside_who_finding(result):
    document = html.render_html(result)
    assert "WHO" in document
    assert "Bedaquiline" in document
