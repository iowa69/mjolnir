"""External tools: how they are found, how their absence reads, and the argv built.

The suite must pass on a machine with no bioinformatics tool installed, so what
is tested here is everything that happens *around* the tools: the argv strings
(pure functions), the error a missing executable produces, and doctor's report of
an environment that has nothing in it. The two tests that genuinely need a binary
skip when it is not there rather than being deleted, because on a machine that
does have samtools they are the only check that the version probe works at all.

The argv assertions are not decoration. ``samtools depth`` takes ``-q`` as base
quality and ``-Q`` as mapping quality, and ``samtools mpileup`` takes them the
other way round; swapping them silently changes every depth in the report.
"""

from __future__ import annotations

import pytest

from mjolnir import config, doctor, shell
from mjolnir.engines import call as call_engine
from mjolnir.engines import depth as depth_engine
from mjolnir.engines import map as map_engine
from mjolnir.utils import MjolnirError

HAVE_SAMTOOLS = shell.first_tool("samtools") is not None


# -------------------------------------------------------------- finding a tool

def test_a_missing_tool_names_the_package_to_install():
    with pytest.raises(MjolnirError) as excinfo:
        shell.tool_path("definitely-not-a-real-tool-xyz", why="testing")
    message = str(excinfo.value)
    assert "not found on PATH" in message
    assert "conda install" in message


def test_an_optional_lookup_returns_none_rather_than_raising():
    assert shell.tool_path("definitely-not-a-real-tool-xyz", required=False) is None
    assert shell.first_tool("definitely-not-a-real-tool-xyz") is None


def test_a_version_string_is_pulled_out_of_a_banner():
    assert shell.parse_version("samtools 1.20\nUsing htslib 1.20") == "1.20"
    assert shell.parse_version("bwa-mem2-2.2.1") == "2.2.1"
    assert shell.parse_version("") == ""


@pytest.mark.skipif(not HAVE_SAMTOOLS, reason="samtools is not installed here")
def test_a_present_tool_yields_a_version_for_the_methods_annex():
    version = shell.tool_version("samtools")
    assert version and any(ch.isdigit() for ch in version)


@pytest.mark.skipif(not HAVE_SAMTOOLS, reason="samtools is not installed here")
def test_captured_versions_accumulate_for_the_report():
    shell.reset_captured()
    shell.record_tool("samtools")
    assert "samtools" in shell.captured_versions()
    shell.reset_captured()


# --------------------------------------------------------------- the argv built

def test_samtools_depth_uses_q_for_base_and_capital_q_for_mapping_quality():
    argv = depth_engine.samtools_depth_argv("s.bam", platform="illumina")
    assert argv[:2] == ["samtools", "depth"]
    assert argv[argv.index("-q") + 1] == str(depth_engine.MIN_BASE_QUALITY)
    assert argv[argv.index("-Q") + 1] == str(depth_engine.MIN_MAPPING_QUALITY)
    assert "-a" in argv, "every position, not only the covered ones"
    assert "-s" in argv, "overlapping mates of a pair are counted once"


def test_ont_depth_does_not_deduplicate_overlapping_mates():
    argv = depth_engine.samtools_depth_argv("s.bam", platform="ont")
    assert "-s" not in argv, "there are no read pairs to overlap on ONT"


def test_mtbseq_compat_reproduces_the_legacy_filters():
    """MTBseq applies no mapping-quality filter at any stage (design §9b)."""
    argv = depth_engine.samtools_depth_argv("s.bam", mtbseq_compat=True)
    assert argv[argv.index("-Q") + 1] == "0"
    assert argv[argv.index("-q") + 1] == str(config.MTBSEQ_MINBQUAL)


def test_the_ont_mapper_uses_the_documented_preset():
    argv = map_engine.minimap2_argv("ref.fa", "reads.fastq.gz")
    assert argv[argv.index("-x") + 1] == config.MINIMAP2_ONT_PRESET == "map-ont"
    assert "--secondary=no" in argv, "a read counted twice invents allele fractions"


def test_the_clair3_model_is_never_defaulted_to_whatever_is_there(tmp_path):
    """An R9 model over R10.4.1 reads produces calls that look ordinary."""
    with pytest.raises(MjolnirError) as excinfo:
        call_engine.clair3_model_path(tmp_path, "r1041_e82_400bps_sup_v500")
    message = str(excinfo.value)
    assert "not found" in message
    assert "R10.4.1" in message and "sup" in message


def test_choosing_an_ont_caller_without_clair3_demands_an_explicit_decision():
    if shell.first_tool("run_clair3.sh"):
        pytest.skip("Clair3 is installed here, so the refusal path is not reachable")
    with pytest.raises(MjolnirError) as excinfo:
        call_engine.choose_caller("ont")
    assert "--allow-degraded-ont-calling" in str(excinfo.value)


# ------------------------------------------------------------------- doctor

def test_doctor_reports_every_tool_with_its_level_and_purpose():
    reports = doctor.check_tools()
    assert reports
    labels = [r.label for r in reports]
    for expected in ("samtools", "bcftools", "minimap2", "kraken2"):
        assert expected in labels
    for report in reports:
        assert report.purpose
        assert report.level in (doctor.LEVEL_REQUIRED, doctor.LEVEL_ALTERNATIVE,
                                doctor.LEVEL_OPTIONAL)


def test_a_missing_required_tool_is_a_problem_with_a_remedy():
    reports = doctor.check_tools()
    problems = doctor.tool_problems(reports)
    for problem in problems:
        assert "install" in problem.lower() or "conda" in problem.lower()


def test_an_alternative_pair_is_satisfied_by_either_member():
    """bwa-mem2 or bwa; skani or mash. Missing one of a pair is not a problem."""
    groups = set(r.group for r in doctor.check_tools() if r.group)
    assert "illumina-mapper" in groups
    assert "ani" in groups


def test_doctor_says_an_unconfigured_screen_was_not_performed():
    report = doctor.check_kraken2_index(None)
    assert report.present is False
    assert report.ready is False
    assert "rather than that the sample was clean" in report.note


def test_doctor_reports_a_standard_index_as_present_but_not_usable(tmp_path):
    """Present is not fit for purpose, and the two must not collapse into a tick."""
    index = tmp_path / "k2_standard_16gb"
    index.mkdir()
    (index / "hash.k2d").write_text("")
    report = doctor.check_kraken2_index(index)
    assert report.present is True
    assert report.usable is False
    assert report.ready is False
    assert "0.0731" in report.note


def test_doctor_accepts_a_declared_pangenome_index(tmp_path):
    index = tmp_path / "myco_pangenome"
    index.mkdir()
    (index / "mjolnir_index.json").write_text('{"mycobacterial_pangenome": true}')
    report = doctor.check_kraken2_index(index)
    assert report.usable is True
    assert report.ready is True


def test_the_python_dependency_report_names_what_each_one_is_for():
    reports = doctor.check_python_deps()
    assert reports
    for report in reports:
        assert report.purpose
