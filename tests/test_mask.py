"""Repeat masks, and the path by which a cohort finds the right one.

A SNP distance is only comparable to a published threshold if the repetitive and
low-complexity regions were excluded first. tbdb ships a curated mask for H37Rv;
nothing ships one for any NTM reference, so :mod:`mjolnir.cohort.mask` computes
one from the reference itself.

The failure this file exists to prevent is not a wrong mask — it is the *right*
mask never being looked for. The joint table was built without a reference, so
the mask lookup had nothing to search beside, fell back to tbdb's H37Rv mask,
found its contigs did not match, and correctly refused to compare. The refusal
was visible; the reason it happened was not.
"""

from __future__ import annotations

import pytest

from mjolnir.cohort import mask as M
from mjolnir.cohort.joint import SampleVariants, build_joint_table
from mjolnir.utils import MjolnirError


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------

def test_overlapping_intervals_merge():
    merged = M.merge([M.Interval("c", 10, 20), M.Interval("c", 15, 30)])
    assert merged == [M.Interval("c", 10, 30)]


def test_adjacent_intervals_merge():
    """Touching intervals are one region; leaving a seam would unmask a base."""
    merged = M.merge([M.Interval("c", 10, 20), M.Interval("c", 20, 30)])
    assert merged == [M.Interval("c", 10, 30)]


def test_intervals_on_different_contigs_do_not_merge():
    merged = M.merge([M.Interval("a", 10, 20), M.Interval("b", 15, 30)])
    assert len(merged) == 2


def test_an_empty_interval_is_dropped():
    assert M.merge([M.Interval("c", 10, 10)]) == []


# ---------------------------------------------------------------------------
# What is found from the sequence alone
# ---------------------------------------------------------------------------

def test_a_homopolymer_run_is_masked():
    """An aligner cannot place a gap inside one, so a difference there is noise."""
    sequence = "ACGT" * 5 + "A" * 20 + "ACGT" * 5
    found = M.homopolymer_intervals(sequence, "c")
    assert found, "the 20-base A run was not found"
    assert found[0].length >= 20


def test_a_short_run_is_not_masked():
    """Masking every 4-mer would swallow the genome."""
    assert M.homopolymer_intervals("ACGTAAAACGT", "c") == []


def test_a_two_base_region_is_low_complexity():
    sequence = "ATATATATAT" * 12
    assert M.low_complexity_intervals(sequence, "c")


def test_ordinary_sequence_is_not_low_complexity():
    sequence = "ACGTTGCAGGATCCAGTCAGGCATTACGGATCCAAGT" * 4
    assert M.low_complexity_intervals(sequence, "c") == []


# ---------------------------------------------------------------------------
# Repeats from the self-alignment
# ---------------------------------------------------------------------------

def _paf(q, qs, qe, t, ts, te, matches, block):
    return "\t".join(str(x) for x in
                     [q, 1000, qs, qe, "+", t, 1000, ts, te, matches, block, 60])


def test_the_self_diagonal_is_dropped():
    """A genome aligned to itself matches itself end to end; that is not a repeat."""
    text = _paf("c", 0, 1000, "c", 0, 1000, 1000, 1000)
    assert M.parse_paf(text) == []


def test_a_repeat_on_the_same_contig_is_kept():
    """Dropping every same-contig hit would miss the IS elements a mask is for."""
    text = _paf("c", 0, 900, "c", 5000, 5900, 880, 900)
    found = M.parse_paf(text)
    assert len(found) == 2, found


def test_a_short_repeat_is_ignored():
    """Below a read length, a repeat does not stop unique placement."""
    text = _paf("c", 0, 50, "c", 5000, 5050, 50, 50)
    assert M.parse_paf(text) == []


def test_a_diverged_repeat_is_ignored():
    """Two copies a read can tell apart do not need masking."""
    text = _paf("c", 0, 900, "c", 5000, 5900, 500, 900)
    assert M.parse_paf(text) == []


def test_a_malformed_paf_line_is_skipped_not_fatal():
    assert M.parse_paf("rubbish\tline\n" + _paf("c", 0, 900, "c", 5000, 5900, 880, 900))


# ---------------------------------------------------------------------------
# Building
# ---------------------------------------------------------------------------

@pytest.fixture
def reference(tmp_path):
    # Long enough that a single 900 bp repeat is a small fraction of it. A
    # short fixture trips the MAX_MASKED_FRACTION ceiling for the right reason
    # and tests nothing about the assembly of the result.
    body = ("ACGTTGCAGGATCCAGTCAGGCATTACGGATCCAAGT" * 400
            + "A" * 30
            + "ACGTTGCAGGATCCAGTCAGGCATTACGGATCCAAGT" * 400)
    path = tmp_path / "ref.fna"
    path.write_text(">c1 test\n{0}\n".format(body))
    return path


def test_build_mask_uses_the_injected_runner(reference):
    """No minimap2 needed to test the assembly of the result."""
    calls = []

    def runner(command):
        calls.append(command)
        return _paf("c1", 0, 900, "c1", 1200, 2100, 880, 900)

    intervals, lengths = M.build_mask(reference, runner=runner)
    assert calls and calls[0][0] == "minimap2"
    assert lengths["c1"] > 0
    assert intervals


def test_an_absurd_mask_is_refused(reference):
    """A mask over a third of the genome is a broken reference or a broken run.

    Distances through it would not mean anything, so it raises rather than
    quietly shrinking every denominator.
    """
    length = len(reference.read_text().split("\n", 1)[1].replace("\n", ""))

    def runner(_command):
        return _paf("c1", 0, length, "c1", 0, length - 1, length, length)

    with pytest.raises(MjolnirError) as excinfo:
        M.build_mask(reference, runner=runner)
    assert "%" in str(excinfo.value)


def test_the_written_mask_says_what_it_is_not(tmp_path):
    """A computed mask is a floor, not an equal of a curated one.

    Someone comparing two cohorts has to be able to tell which kind each had, so
    the method and its limits are in the file rather than in a README.
    """
    target = tmp_path / "out.bed"
    M.write_mask(target, [M.Interval("c1", 10, 40)],
                 reference="ref.fna", lengths={"c1": 1000})
    text = target.read_text()
    assert "NOT equivalent to a curated mask" in text
    assert "PE/PPE" in text
    assert "c1\t10\t40" in text


def test_the_default_mask_path_sits_beside_the_reference():
    assert M.default_mask_path("/db/ani/x.fna").name == "x.fna.mask.bed"


# ---------------------------------------------------------------------------
# The lookup path
# ---------------------------------------------------------------------------

def test_the_joint_table_carries_the_reference_the_mask_is_chosen_by():
    """Built without it, every NTM cohort silently fell back to tbdb's H37Rv mask.

    The mask is looked for beside the reference. A joint table with an empty
    reference has nothing to look beside, so the computed mask sitting right
    next to the genome was never found and the cohort refused to compare.
    """
    entries = [
        SampleVariants(sample_id="a", reference="/db/ani/chimaera.fna"),
        SampleVariants(sample_id="b", reference="/db/ani/chimaera.fna"),
    ]
    table = build_joint_table(entries, reference="/db/ani/chimaera.fna")
    assert table.reference == "/db/ani/chimaera.fna"
    assert M.default_mask_path(table.reference).name == "chimaera.fna.mask.bed"
