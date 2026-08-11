"""How ``mash dist`` is invoked, which decides whether a species call means anything.

``mash dist a b c query`` does not compare *query* against a, b and c. It names
only the first path as the reference and treats every later path — the query
included — as a query against it. Passing the reference set as a list therefore
produces a table of the reference genomes' distances *to each other*, in which
the sample is one row among many and is almost never the top hit.

That is not a subtle inaccuracy. On real data an *M. chimaera* isolate and an
*M. bovis* isolate both reported an identical "99.1268% against H37Rv", because
that number is the distance from H37Rv to *M. bovis* and neither sample took any
part in producing it. The chimaera isolate was called MTBC.

Nothing about the output looked wrong: the hits were well-formed, sorted, and
carried plausible ANI values. Only the sample was missing. So these tests assert
the *shape of the invocation* rather than the shape of the result, because a
result-shaped test passed throughout.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mjolnir.typing import species


@pytest.fixture
def reference_set(tmp_path):
    """Three reference genomes on disk, with a manifest naming each."""
    names = [("a.fna", "Mycobacterium tuberculosis"),
             ("b.fna", "Mycobacterium intracellulare subsp. chimaera"),
             ("c.fna", "Mycobacteroides abscessus")]
    refs = []
    for filename, name in names:
        path = tmp_path / filename
        path.write_text(">seq\nACGTACGTACGT\n")
        refs.append(species.ReferenceGenome(path=path, name=name))
    return refs


class RecordingRunner:
    """Captures every command, and answers a dist call with one row per reference."""

    def __init__(self, refs, top=None):
        self.commands = []
        self.refs = refs
        self.top = top or refs[0]

    def __call__(self, command):
        self.commands.append([str(part) for part in command])
        if "dist" not in command:
            return ""
        lines = []
        for index, ref in enumerate(self.refs):
            distance = 0.005 if ref is self.top else 0.18 + index * 0.01
            lines.append("{0}\tquery.fq\t{1}\t0\t900/1000".format(ref.path, distance))
        return "\n".join(lines) + "\n"

    @property
    def dist_command(self):
        return next(c for c in self.commands if "dist" in c)

    @property
    def sketch_command(self):
        return next((c for c in self.commands if "sketch" in c), None)


def test_the_query_is_compared_against_the_references_not_the_references_against_each_other(
        tmp_path, reference_set):
    """``mash dist`` must receive exactly one reference argument and one query.

    This is the whole bug. With more than one path before the query, mash
    silently reinterprets the command and the sample stops being measured.
    """
    query = tmp_path / "query.fq"
    query.write_text("@r\nACGT\n+\n!!!!\n")
    runner = RecordingRunner(reference_set)

    species.run_mash(query, reference_set, from_reads=True, runner=runner)

    command = runner.dist_command
    paths = [part for part in command
             if part.endswith((".fna", ".fq", ".fastq", ".msh", ".fasta"))]
    assert len(paths) == 2, "mash dist got {0} paths: {1}".format(len(paths), paths)
    assert paths[-1] == str(query), "the query must be the last argument"
    assert paths[0].endswith(".msh"), "the references must be given as one sketch"


def test_the_references_are_sketched_into_a_single_file_first(tmp_path, reference_set):
    """Every reference must reach the sketch, or it cannot be a candidate."""
    query = tmp_path / "query.fq"
    query.write_text("@r\nACGT\n+\n!!!!\n")
    runner = RecordingRunner(reference_set)

    species.run_mash(query, reference_set, from_reads=True, runner=runner)

    sketch = runner.sketch_command
    assert sketch is not None, "the reference set was never sketched"
    for reference in reference_set:
        assert str(reference.path) in sketch


def test_no_reference_genome_is_ever_passed_as_a_query(tmp_path, reference_set):
    """A reference appearing after the reference argument becomes a query.

    That is what produced a candidate list of reference-to-reference distances.
    """
    query = tmp_path / "query.fq"
    query.write_text("@r\nACGT\n+\n!!!!\n")
    runner = RecordingRunner(reference_set)

    species.run_mash(query, reference_set, from_reads=True, runner=runner)

    command = runner.dist_command
    for reference in reference_set:
        assert str(reference.path) not in command


def test_the_top_hit_is_the_reference_the_query_is_closest_to(tmp_path, reference_set):
    """With the invocation right, the nearest reference is the one that is named."""
    query = tmp_path / "query.fq"
    query.write_text("@r\nACGT\n+\n!!!!\n")
    chimaera = reference_set[1]
    runner = RecordingRunner(reference_set, top=chimaera)

    matches = species.run_mash(query, reference_set, from_reads=True, runner=runner)

    assert matches, "no matches parsed"
    assert matches[0].name == chimaera.name
    assert matches[0].ani > 99.0


def test_a_supplied_sketch_is_used_directly(tmp_path, reference_set):
    """An operator-built sketch skips the sketching step but keeps the shape."""
    query = tmp_path / "query.fq"
    query.write_text("@r\nACGT\n+\n!!!!\n")
    sketch = tmp_path / "mycobacteria.msh"
    sketch.write_bytes(b"sketch")
    runner = RecordingRunner(reference_set)

    species.run_mash(query, reference_set, from_reads=True, sketch=sketch, runner=runner)

    assert runner.sketch_command is None, "a supplied sketch must not be rebuilt"
    command = runner.dist_command
    assert str(sketch) in command
    assert command[-1] == str(query)


def test_read_input_discards_single_copy_kmers(tmp_path, reference_set):
    """``-r -m 2`` on reads; without it a read set's distance to everything inflates."""
    query = tmp_path / "query.fq"
    query.write_text("@r\nACGT\n+\n!!!!\n")
    runner = RecordingRunner(reference_set)

    species.run_mash(query, reference_set, from_reads=True, runner=runner)

    command = runner.dist_command
    assert "-r" in command and "-m" in command
    assert command[command.index("-m") + 1] == "2"
