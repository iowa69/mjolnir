"""Repeat masks for references that do not have one.

A pairwise SNP distance is only comparable to a published threshold if the
repetitive and low-complexity regions have been excluded, because a read that
maps ambiguously produces a difference that is a mapping artefact rather than a
mutation. For the M. tuberculosis complex tbdb ships a mask and the matter is
settled. For everything else there is nothing, and design §9 says NTM references
get their own mask computed at database build time — which nothing implemented.

The consequence was not a wrong number, because :mod:`mjolnir.cohort.distance`
refuses a mask whose contigs do not match the reference. It was no number at
all: a four-isolate *M. chimaera* cohort produced no distances, no clusters and
no outbreak analysis, which is the entire reason such a collection exists.

**What this computes, and what it cannot.** Two things are found from the
reference alone:

*Repeats*, by aligning the genome to itself and keeping every alignment that is
not a sequence against itself. Anything covered twice is somewhere a short read
cannot be placed uniquely.

*Low-complexity*, by a sliding window over the sequence: homopolymers and runs
whose base composition collapses to one or two letters. These are where an
aligner's gap placement becomes arbitrary.

What it cannot find is the third category tbdb's mask carries — regions that are
*empirically* error-prone across many samples, like the PE/PPE families, which
were identified by looking at thousands of isolates rather than at one genome.
So a computed mask is a floor, not an equal of a curated one, and it says so:
:func:`write_mask` records the method in the BED header and the report prints it.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..config import source_for
from ..utils import LOG, MjolnirError, PathLike, ensure_dir, require, require_file

#: SOURCE: Mjolnir policy. Short repeats do not defeat a 150 bp read pair; the
#: shortest IS element in these genomes is around 700 bp and is what a mask is
#: for. Set below the read length and the mask swallows genuine variation.
MIN_REPEAT_LENGTH = 200

#: SOURCE: Mjolnir policy. Two copies at this identity or above cannot be told
#: apart by a short read, so both are excluded.
MIN_REPEAT_IDENTITY = 0.90

#: SOURCE: Mjolnir policy. Window and threshold for the low-complexity scan: a
#: window in which two bases account for this much of the sequence is one where
#: an aligner's gap placement is arbitrary.
COMPLEXITY_WINDOW = 48
COMPLEXITY_MAX_TWO_BASE_FRACTION = 0.92

#: What a mask may exclude before it stops being a mask and becomes a decision
#: to not measure. SOURCE: Mjolnir policy, informed by the ~6% tbdb excludes
#: from H37Rv; an order of magnitude more than that is a broken reference or a
#: broken computation, and either way the distances would not mean anything.
MAX_MASKED_FRACTION = 0.35


@dataclass(frozen=True)
class Interval:
    """A half-open BED interval: ``start`` is 0-based, ``end`` exclusive."""

    contig: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


def merge(intervals: Iterable[Interval]) -> List[Interval]:
    """Sorted, non-overlapping intervals. Adjacent ones are joined."""
    ordered = sorted(intervals, key=lambda i: (i.contig, i.start, i.end))
    out: List[Interval] = []
    for interval in ordered:
        if interval.end <= interval.start:
            continue
        if out and out[-1].contig == interval.contig and interval.start <= out[-1].end:
            if interval.end > out[-1].end:
                out[-1] = Interval(out[-1].contig, out[-1].start, interval.end)
            continue
        out.append(interval)
    return out


# ---------------------------------------------------------------------------
# Low complexity
# ---------------------------------------------------------------------------

def low_complexity_intervals(sequence: str, contig: str, *,
                             window: int = COMPLEXITY_WINDOW,
                             limit: float = COMPLEXITY_MAX_TWO_BASE_FRACTION
                             ) -> List[Interval]:
    """Windows dominated by one or two bases.

    Deliberately cheap and deliberately blunt: this is a mask, and a region
    wrongly included costs a little sensitivity while a region wrongly excluded
    costs a false SNP difference between two isolates.
    """
    upper = sequence.upper()
    found: List[Interval] = []
    if len(upper) < window:
        return found
    for start in range(0, len(upper) - window + 1, window // 2):
        chunk = upper[start:start + window]
        counts = sorted((chunk.count(base) for base in "ACGT"), reverse=True)
        if not counts or counts[0] == 0:
            continue
        if (counts[0] + counts[1]) / float(window) >= limit:
            found.append(Interval(contig, start, start + window))
    return merge(found)


def homopolymer_intervals(sequence: str, contig: str, *,
                          min_length: int = 12) -> List[Interval]:
    """Runs of one base long enough that an aligner cannot place a gap in them."""
    found: List[Interval] = []
    for match in re.finditer(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})" % (
            min_length, min_length, min_length, min_length), sequence.upper()):
        found.append(Interval(contig, match.start(), match.end()))
    return merge(found)


# ---------------------------------------------------------------------------
# Repeats
# ---------------------------------------------------------------------------

def self_alignment_argv(reference: PathLike, *, threads: int = 1) -> List[str]:
    """``minimap2`` aligning a genome to itself.

    ``-DP`` skips self-hits at the seeding stage and ``-c`` asks for base-level
    alignment so identity can be computed. ``asm20`` tolerates the divergence
    between real repeat copies, which are rarely identical.
    """
    return ["minimap2", "-DP", "-c", "-x", "asm20",
            "-t", str(max(1, int(threads))),
            str(reference), str(reference)]


def parse_paf(text: str, *, min_length: int = MIN_REPEAT_LENGTH,
              min_identity: float = MIN_REPEAT_IDENTITY) -> List[Interval]:
    """Repeat intervals from PAF, dropping the diagonal.

    A genome aligned to itself matches itself end to end; that alignment is
    discarded by coordinate rather than by name, because a repeat can and does
    occur on the same contig.
    """
    found: List[Interval] = []
    for line in (text or "").splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 12:
            continue
        try:
            q_name, q_start, q_end = fields[0], int(fields[2]), int(fields[3])
            t_name, t_start, t_end = fields[5], int(fields[7]), int(fields[8])
            matches, block = int(fields[9]), int(fields[10])
        except ValueError:
            continue
        if q_name == t_name and q_start == t_start and q_end == t_end:
            continue
        if block <= 0 or (q_end - q_start) < min_length:
            continue
        if matches / float(block) < min_identity:
            continue
        found.append(Interval(q_name, q_start, q_end))
        found.append(Interval(t_name, t_start, t_end))
    return merge(found)


# ---------------------------------------------------------------------------
# Building and writing
# ---------------------------------------------------------------------------

def build_mask(reference: PathLike, *, threads: int = 1,
               runner: Optional[object] = None) -> Tuple[List[Interval], Dict[str, int]]:
    """``(intervals, contig lengths)`` for a reference with no curated mask."""
    from ..engines.annotate import load_fasta

    resolved = require_file(reference, "the reference to build a mask for")
    sequences = load_fasta(resolved)
    if not sequences:
        raise MjolnirError("{0} contains no sequence".format(resolved))
    lengths = {name: len(seq) for name, seq in sequences.items()}

    intervals: List[Interval] = []
    for name, sequence in sequences.items():
        intervals.extend(low_complexity_intervals(sequence, name))
        intervals.extend(homopolymer_intervals(sequence, name))

    command = self_alignment_argv(resolved, threads=threads)
    if runner is None:
        import subprocess

        require("minimap2", why="computing a repeat mask for a reference")
        LOG.info("self-aligning %s to find repeats", resolved.name)
        completed = subprocess.run(command, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        if completed.returncode != 0:
            raise MjolnirError(
                "minimap2 failed while self-aligning {0} (exit {1}):\n  {2}".format(
                    resolved, completed.returncode,
                    completed.stderr.decode("utf-8", "replace")[-500:]))
        paf = completed.stdout.decode("utf-8", "replace")
    else:
        paf = runner(command)

    intervals.extend(parse_paf(paf))
    merged = merge(intervals)

    total = sum(lengths.values())
    masked = sum(i.length for i in merged)
    if total and masked / float(total) > MAX_MASKED_FRACTION:
        raise MjolnirError(
            "a computed mask for {0} would exclude {1:.1%} of it, above the {2:.0%} "
            "ceiling. A mask that large is a broken reference or a broken "
            "computation, and distances through it would not mean anything.\n"
            "  Inspect the reference, or supply a curated mask with --mask.".format(
                resolved.name, masked / float(total), MAX_MASKED_FRACTION))
    return merged, lengths


def write_mask(path: PathLike, intervals: Sequence[Interval],
               *, reference: str, lengths: Optional[Dict[str, int]] = None) -> Path:
    """Write a BED, with a header recording how it was made.

    The header matters as much as the intervals. A computed mask finds repeats
    and low-complexity sequence from one genome; it does not find the
    empirically error-prone regions a curated mask carries, which were derived
    from thousands of isolates. A reader comparing two cohorts needs to know
    which kind they had.
    """
    target = Path(path)
    ensure_dir(target.parent)
    total = sum((lengths or {}).values())
    masked = sum(i.length for i in intervals)
    with open(str(target), "w") as handle:
        handle.write("# Mjolnir computed repeat mask\n")
        handle.write("# reference: {0}\n".format(reference))
        handle.write("# method: minimap2 self-alignment (repeats >= {0} bp at >= "
                     "{1:.0%} identity), plus low-complexity windows and "
                     "homopolymers\n".format(MIN_REPEAT_LENGTH, MIN_REPEAT_IDENTITY))
        handle.write("# NOT equivalent to a curated mask: empirically error-prone "
                     "regions such as PE/PPE were identified across thousands of "
                     "isolates and cannot be found from one genome\n")
        if total:
            handle.write("# masked: {0} bp of {1} ({2:.2%})\n".format(
                masked, total, masked / float(total)))
        handle.write("# source: {0}\n".format(source_for("cluster_snp_strict")))
        for interval in intervals:
            handle.write("{0}\t{1}\t{2}\n".format(
                interval.contig, interval.start, interval.end))
    LOG.info("wrote %s: %d intervals, %d bp masked", target, len(intervals), masked)
    return target


def default_mask_path(reference: PathLike) -> Path:
    """Where a computed mask for *reference* lives: beside it, ``.mask.bed``."""
    return Path(str(reference) + ".mask.bed")


__all__ = [
    "Interval", "build_mask", "default_mask_path", "homopolymer_intervals",
    "low_complexity_intervals", "merge", "parse_paf", "self_alignment_argv",
    "write_mask", "MAX_MASKED_FRACTION", "MIN_REPEAT_IDENTITY",
    "MIN_REPEAT_LENGTH",
]
