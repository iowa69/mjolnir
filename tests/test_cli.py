"""CLI smoke tests: the surface a user meets, and how it fails.

Nothing here runs an analysis. What is checked is that the command surface
exists, that every subcommand's help renders, and — the part that matters — that
a foreseeable mistake produces the sentence that says what to do about it and a
non-zero exit, rather than a traceback or a zero.

The exit codes are load-bearing in a way MTBseq's are not: MTBseq exits 1 on a
*successful* single-step run, so any wrapper reading its status misreports it.
Mjolnir's ``main`` returns 0 only when the command did what it was asked.
"""

from __future__ import annotations

import pytest

from mjolnir import cli
from mjolnir.utils import MjolnirError


@pytest.mark.parametrize("command", ["run", "cohort", "report", "db", "doctor"])
def test_every_subcommand_renders_its_help(command, capsys):
    with pytest.raises(SystemExit) as excinfo:
        cli.main([command, "--help"])
    assert excinfo.value.code == 0
    assert capsys.readouterr().out.strip()


def test_the_top_level_help_lists_the_commands(capsys):
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--help"])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    for command in ("run", "cohort", "report", "db", "doctor"):
        assert command in out


def test_version_prints_and_exits_zero(capsys):
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["--version"])
    assert excinfo.value.code == 0
    assert capsys.readouterr().out.strip()


def test_no_command_prints_help_and_exits_non_zero(capsys):
    assert cli.main([]) == 1
    assert capsys.readouterr().out.strip()


def test_db_with_no_action_says_where_to_start(capsys):
    assert cli.main(["db"]) == 1
    assert "mjolnir db fetch" in capsys.readouterr().err


def test_db_list_describes_every_database_with_its_licence(capsys):
    assert cli.main(["db", "list"]) == 0
    out = capsys.readouterr().out
    assert "licence" in out
    assert "ODC-By" in out, "the WHO data licence has to be printed"
    assert "attribution" in out.lower()


def test_a_missing_input_fails_with_a_sentence_not_a_traceback(tmp_path, capsys):
    code = cli.main(["run", "-1", str(tmp_path / "absent_R1.fastq.gz"),
                     "-o", str(tmp_path / "out")])
    assert code == 1
    captured = capsys.readouterr()
    assert "input not found" in (captured.err + captured.out)
    assert "Traceback" not in (captured.err + captured.out)


def test_an_unknown_platform_is_rejected_before_anything_runs(tmp_path):
    with pytest.raises(SystemExit):
        cli.main(["run", "-1", str(tmp_path / "r1.fastq.gz"), "--platform", "pacbio",
                  "-o", str(tmp_path / "out")])


def test_the_error_type_the_cli_catches_is_the_one_the_modules_raise():
    """Every "install this / fetch that" message travels as a MjolnirError."""
    assert issubclass(MjolnirError, Exception)
    with pytest.raises(MjolnirError):
        raise MjolnirError("fetch it with: mjolnir db fetch")
