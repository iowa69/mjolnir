"""Read alignment: bwa-mem2 (bwa fallback) for Illumina, minimap2 for ONT.

Every command this module issues is built by a pure function that returns an
argv list and executes nothing. That split is not decoration: a mapping command
is the one part of the pipeline that cannot be unit-tested by running it — the
machine that runs the tests is not the machine with 40 GB of reads on it — so
the commands are tested as data and only the thin wrappers at the bottom of the
file spawn a process.

Three choices here change the numbers downstream and are therefore stated
rather than buried.

**Unmapped reads stay in the BAM.** No ``--sam-hit-only``, no ``-F 4`` filter at
sort time. Mapped fraction is a contamination signal the design requires
(§8.3), and it cannot be computed from a file that has already discarded its
denominator.

**The aligner is piped straight into ``samtools sort``.** A 4.4 Mb genome at
100x is a few gigabytes of SAM text if it is written out first, and that text is
never read by anything else.

**Duplicate marking is Illumina-only.** ``samtools markdup`` decides duplication
from identical outer coordinates, which for ONT reads of variable length from a
single molecule is not a duplicate but a coincidence. Marking them would remove
real coverage, so ONT alignments are left unmarked and ``duplicate_fraction``
is reported as unmeasured rather than as zero.

``shell.py`` owns single-command execution and version capture for the rest of
the tool. What lives here instead is a *pipeline* runner, because a chain of
processes joined by pipes is not a single command, and because the stderr of a
mapper must go to a file: minimap2 emits a warning per malformed read, and on a
pipe those fill the 64 KiB buffer, block the mapper, starve the sorter of input
and hang the run with no error at all.
"""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence

from ..config import (
    Config,
    HET_MIN_QUAL,
    MAPPERS,
    MINIMAP2_ONT_PRESET,
)
from ..records import (
    PLATFORM_FASTA,
    PLATFORM_ILLUMINA,
    PLATFORM_ONT,
    SampleInput,
    normalise_platform,
)
from ..utils import (
    LOG,
    MjolnirError,
    PathLike,
    conda_package,
    ensure_dir,
    first_available,
    require,
    require_file,
    tempdir,
)

# ---------------------------------------------------------------------------
# Alignment-filter parameters
#
# These are tool-invocation parameters, not clinical thresholds — the clinical
# thresholds all live in config.py and are imported, never re-spelled. Each one
# still carries the source that fixes its value, because a MAPQ floor changes
# which variants exist and a depth cap changes every allele fraction.
# ---------------------------------------------------------------------------

#: SOURCE: Marin et al. 2022 (as cited in the Mjolnir design §9): raising the
#: mapping-quality threshold to MQ >= 40 gave precision 99.1% / recall 85.8% and
#: outperformed blanket masking of repetitive regions. MTBseq applies no MAPQ
#: filter at any stage (design §9b), which is precisely why a Mjolnir number and
#: an MTBseq number at the same position can differ; `--compat mtbseq` drops
#: this floor to 0 so the two can be reconciled.
SRC_MARIN_2022 = (
    "Marin et al. 2022, mapping-quality thresholding versus masking in M. "
    "tuberculosis (MQ>=40: precision 99.1%, recall 85.8%), as cited in the "
    "Mjolnir design §9"
)
MIN_MAPPING_QUALITY = 40

#: SOURCE: Sobkowiak et al. 2018 (MixInfect), the same Q>=20 floor config.py
#: registers as ``het_min_qual``. Deliberately the same number: a base the
#: heterozygosity analysis would refuse to count must not be counted into the
#: allele fraction that feeds it, or the two would disagree about the same site.
MIN_BASE_QUALITY = HET_MIN_QUAL

#: SOURCE: Mjolnir policy, forced by design §9b. ``samtools mpileup`` caps depth
#: at 250 by default and MTBseq never passes ``-d``, so on a deep sample MTBseq's
#: allele frequencies are computed from downsampled counts. Mjolnir raises the
#: cap far above any depth a clinical isolate reaches, so the fraction it prints
#: is the fraction that was sequenced. A cap is still needed: one pathological
#: pile-up column must not allocate unbounded memory.
MAX_PILEUP_DEPTH = 8000

#: SOURCE: bwa/bwa-mem2 manual, ``-K INT`` — process INT input bases per batch
#: regardless of thread count. Without it, the number of threads changes the
#: batching and therefore the alignment of a handful of reads, so two runs of the
#: same sample on differently-loaded machines can produce different variant
#: calls. 100 Mbase is the value the GATK best-practice invocation uses.
BWA_DETERMINISTIC_CHUNK = 100000000

#: SOURCE: Mjolnir policy. Sort memory per thread; conservative because several
#: samples may be in flight on one machine.
SORT_MEMORY_PER_THREAD = "768M"

#: SAM platform codes (SAM specification, @RG PL field). Written into the read
#: group so that a BAM handed to a different tool still declares what it is.
_SAM_PLATFORM: Dict[str, str] = {
    PLATFORM_ILLUMINA: "ILLUMINA",
    PLATFORM_ONT: "ONT",
}

#: What each mapper leaves beside the FASTA once it has been indexed. Used to
#: tell "not indexed yet" from "indexed with the other tool", because the two
#: need different error messages.
INDEX_SUFFIXES: Dict[str, Sequence[str]] = {
    "bwa-mem2": (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    "bwa": (".amb", ".ann", ".bwt", ".pac", ".sa"),
}


# ---------------------------------------------------------------------------
# Process execution: one pipeline runner and one streaming reader
# ---------------------------------------------------------------------------

def _tail(path: Path, limit: int = 2000) -> str:
    """Last *limit* characters of a log file, for an error message."""
    try:
        return path.read_text(errors="replace")[-limit:].strip()
    except OSError:
        return ""


def run_pipeline(stages: Sequence[Sequence[str]], *,
                 stdout_path: Optional[PathLike] = None,
                 log_dir: Optional[PathLike] = None,
                 keep_logs: bool = False) -> List[str]:
    """Run ``a | b | c`` without a shell, returning each stage's stderr tail.

    Failure is reported from the downstream end first. When ``samtools sort``
    dies, the aligner upstream of it also exits non-zero — but only because it
    got SIGPIPE, and printing that as the cause sends the reader looking at the
    wrong tool.
    """
    commands = [[str(arg) for arg in stage] for stage in stages]
    if not commands or any(not stage for stage in commands):
        raise MjolnirError("run_pipeline needs at least one non-empty command")
    for stage in commands:
        require(stage[0], "alignment and calling pipeline")

    LOG.debug("exec: %s", " | ".join(" ".join(stage) for stage in commands))

    with tempdir(prefix="mjolnir.pipe.", keep=keep_logs) as scratch:
        where = ensure_dir(log_dir) if log_dir is not None else scratch
        logs = [Path(where) / "{0:02d}.{1}.err".format(index, Path(stage[0]).name)
                for index, stage in enumerate(commands)]
        handles: List = []
        procs: List[subprocess.Popen] = []
        out_handle = None
        try:
            for path in logs:
                handles.append(open(str(path), "wb"))
            if stdout_path is not None:
                ensure_dir(Path(stdout_path).parent)
                out_handle = open(str(stdout_path), "wb")

            previous = None
            for index, stage in enumerate(commands):
                is_last = index == len(commands) - 1
                target = out_handle if is_last else subprocess.PIPE
                try:
                    proc = subprocess.Popen(stage, stdin=previous, stdout=target,
                                            stderr=handles[index])
                except OSError as exc:
                    for started in procs:
                        started.kill()
                    raise MjolnirError(
                        "could not start '{0}': {1}\n"
                        "  conda install -c conda-forge -c bioconda {2}".format(
                            stage[0], exc, conda_package(stage[0]))
                    ) from exc
                if previous is not None:
                    # Drop our copy of the upstream pipe so that this stage is the
                    # only reader; otherwise the upstream process never sees EOF.
                    previous.close()
                previous = proc.stdout
                procs.append(proc)

            for proc in procs:
                proc.wait()
        finally:
            for handle in handles:
                handle.close()
            if out_handle is not None:
                out_handle.close()

        tails = [_tail(path) for path in logs]
        for index in range(len(procs) - 1, -1, -1):
            if procs[index].returncode != 0:
                raise MjolnirError(
                    "{0} failed (exit {1}) in: {2}\n{3}".format(
                        commands[index][0], procs[index].returncode,
                        " | ".join(" ".join(stage) for stage in commands),
                        tails[index] or "(no stderr output)")
                )
        return tails


def iter_output(argv: Sequence[str], *, log_dir: Optional[PathLike] = None) -> Iterator[str]:
    """Yield the stdout lines of one command, raising if it exits non-zero.

    Used by ``depth.py`` and ``pileup.py``, whose inputs are one command each and
    whose outputs are megabytes of text that never need to exist as a file. The
    generator kills the process if the consumer abandons it, so a parser that
    raises part-way through does not leave samtools running.
    """
    command = [str(arg) for arg in argv]
    if not command:
        raise MjolnirError("iter_output needs a command")
    require(command[0], "pileup and coverage")
    LOG.debug("exec: %s", " ".join(command))

    with tempdir(prefix="mjolnir.out.") as scratch:
        log_path = Path(ensure_dir(log_dir) if log_dir is not None else scratch)
        log_path = log_path / "{0}.err".format(Path(command[0]).name)
        with open(str(log_path), "wb") as errfile:
            try:
                proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                                        stderr=errfile, universal_newlines=True)
            except OSError as exc:
                raise MjolnirError(
                    "could not start '{0}': {1}\n"
                    "  conda install -c conda-forge -c bioconda {2}".format(
                        command[0], exc, conda_package(command[0]))
                ) from exc
            try:
                for line in proc.stdout:
                    yield line.rstrip("\n")
            finally:
                if proc.stdout is not None:
                    proc.stdout.close()
                if proc.poll() is None:
                    proc.kill()
                proc.wait()
        if proc.returncode != 0:
            raise MjolnirError(
                "{0} failed (exit {1}):\n{2}".format(
                    command[0], proc.returncode, _tail(log_path) or "(no stderr output)")
            )


# ---------------------------------------------------------------------------
# The result record
# ---------------------------------------------------------------------------

@dataclass
class Alignment:
    """A sorted, indexed BAM and the exact commands that produced it.

    ``commands`` is kept so the methods annex can print what was run rather than
    a prose description of it: the report's claim to reproducibility is the argv,
    not a sentence about bwa.
    """

    sample_id: str
    platform: str
    bam: Path
    reference: Path
    mapper: str
    commands: List[List[str]] = field(default_factory=list)
    duplicates_marked: bool = False
    #: Set when duplicate marking was deliberately not attempted, so QC reports
    #: "not measured" instead of a zero duplicate fraction.
    duplicates_note: str = ""

    @property
    def index(self) -> Path:
        return Path(str(self.bam) + ".bai")


# ---------------------------------------------------------------------------
# Pure command builders
# ---------------------------------------------------------------------------

def read_group(sample_id: str, platform: str, library: str = "") -> str:
    """A SAM ``@RG`` line for the aligner's ``-R``.

    The tabs are literal backslash-t: bwa and minimap2 both unescape the string
    themselves, and since Mjolnir spawns processes without a shell there is
    nothing else in the way to do it.
    """
    plat = normalise_platform(platform)
    code = _SAM_PLATFORM.get(plat)
    if code is None:
        raise MjolnirError(
            "no read group applies to {0} input: an assembly has no reads".format(plat))
    return "@RG\\tID:{0}\\tSM:{0}\\tPL:{1}\\tLB:{2}".format(
        sample_id, code, library or sample_id)


def bwa_mem_argv(tool: str, reference: PathLike, r1: PathLike,
                 r2: Optional[PathLike] = None, *, threads: int = 1,
                 sample_id: str = "", chunk: int = BWA_DETERMINISTIC_CHUNK) -> List[str]:
    """``bwa-mem2 mem`` / ``bwa mem`` for Illumina reads, emitting SAM on stdout.

    ``-M`` is deliberately absent. It flags supplementary alignments as secondary
    for the benefit of Picard-era tools, and a caller that then discards
    secondary alignments loses exactly the split reads that carry the large
    deletions the LoF rules depend on.
    """
    if tool not in ("bwa-mem2", "bwa"):
        raise MjolnirError(
            "unknown Illumina mapper {0!r}; expected one of {1}".format(
                tool, ", ".join(MAPPERS[PLATFORM_ILLUMINA])))
    argv = [tool, "mem", "-t", str(max(1, threads)), "-K", str(chunk)]
    if sample_id:
        argv += ["-R", read_group(sample_id, PLATFORM_ILLUMINA)]
    argv.append(str(reference))
    argv.append(str(r1))
    if r2 is not None:
        argv.append(str(r2))
    return argv


def minimap2_argv(reference: PathLike, reads: PathLike, *, threads: int = 1,
                  sample_id: str = "", preset: str = MINIMAP2_ONT_PRESET) -> List[str]:
    """``minimap2 -a -x map-ont`` for ONT reads, emitting SAM on stdout.

    ``--secondary=no`` because a variant caller counting the same read twice at
    two loci invents allele fractions in every repeat family; ``-L`` because a
    long read spanning a tandem repeat can exceed the 65,535-operation CIGAR
    limit of BAM and would otherwise be written unreadably; ``--MD`` because
    downstream QC and some callers want the tag and computing it later means
    reading the reference again.
    """
    return [
        "minimap2", "-a", "-x", str(preset),
        "-t", str(max(1, threads)),
        "--secondary=no", "-L", "--MD",
    ] + (["-R", read_group(sample_id, PLATFORM_ONT)] if sample_id else []) + [
        str(reference), str(reads),
    ]


def samtools_fixmate_argv(*, threads: int = 1) -> List[str]:
    """``samtools fixmate -m`` on the aligner's name-ordered stream.

    ``-m`` adds the mate-score tag that ``markdup`` needs to choose which read of
    a duplicate pair to keep; without it markdup keeps an arbitrary one.
    """
    return ["samtools", "fixmate", "-m", "-u", "-@", str(max(1, threads)), "-", "-"]


def samtools_sort_argv(out_bam: Optional[PathLike] = None, *, threads: int = 1,
                       tmp_prefix: Optional[PathLike] = None,
                       uncompressed: bool = False,
                       memory: str = SORT_MEMORY_PER_THREAD) -> List[str]:
    """``samtools sort`` reading SAM on stdin.

    With *out_bam* it writes the finished file; with ``uncompressed`` it writes a
    raw BAM stream for the next stage, which saves compressing data that is about
    to be decompressed again.
    """
    argv = ["samtools", "sort", "-@", str(max(1, threads)), "-m", str(memory)]
    if tmp_prefix is not None:
        argv += ["-T", str(tmp_prefix)]
    if uncompressed:
        argv.append("-u")
    if out_bam is not None:
        argv += ["-o", str(out_bam)]
    argv.append("-")
    return argv


def samtools_markdup_argv(out_bam: PathLike, *, threads: int = 1,
                          stats_path: Optional[PathLike] = None) -> List[str]:
    """``samtools markdup`` — marks, never removes.

    A removed duplicate cannot be counted, and ``duplicate_fraction`` is a
    library-quality signal the QC panel reports.
    """
    argv = ["samtools", "markdup", "-@", str(max(1, threads))]
    if stats_path is not None:
        argv += ["-f", str(stats_path)]
    argv += ["-", str(out_bam)]
    return argv


def samtools_index_argv(bam: PathLike, *, threads: int = 1) -> List[str]:
    return ["samtools", "index", "-@", str(max(1, threads)), str(bam)]


def samtools_faidx_argv(reference: PathLike) -> List[str]:
    return ["samtools", "faidx", str(reference)]


def index_reference_argv(reference: PathLike, tool: str) -> List[str]:
    """The command that builds *tool*'s index beside the FASTA."""
    if tool not in INDEX_SUFFIXES:
        raise MjolnirError(
            "no index command known for {0!r}; expected one of {1}".format(
                tool, ", ".join(sorted(INDEX_SUFFIXES))))
    return [tool, "index", str(reference)]


def mapping_pipeline(sample: SampleInput, reference: PathLike, out_bam: PathLike, *,
                     mapper: str, threads: int = 1,
                     mark_duplicates: bool = False,
                     tmp_prefix: Optional[PathLike] = None,
                     markdup_stats: Optional[PathLike] = None) -> List[List[str]]:
    """The whole ``aligner | ... | samtools`` chain for one sample, as argv lists.

    This is the function the tests exercise. ``map_reads`` below is its wrapper
    and adds nothing to the command line, so a test of this function is a test of
    what actually runs.
    """
    platform = normalise_platform(sample.platform)
    if platform == PLATFORM_FASTA:
        raise MjolnirError(
            "sample {0!r} is an assembly: there are no reads to map. FASTA input "
            "is compared directly to the reference and never goes through the "
            "aligner (design §7)".format(sample.sample_id))

    map_threads = max(1, threads - 1) if threads > 1 else 1
    sort_threads = max(1, threads // 4)

    if platform == PLATFORM_ILLUMINA:
        reads = sample.paths
        stages = [bwa_mem_argv(mapper, reference, reads[0],
                               reads[1] if len(reads) > 1 else None,
                               threads=map_threads, sample_id=sample.sample_id)]
    else:
        stages = [minimap2_argv(reference, sample.paths[0], threads=map_threads,
                                sample_id=sample.sample_id)]

    if mark_duplicates:
        if platform != PLATFORM_ILLUMINA:
            raise MjolnirError(
                "duplicate marking was requested for {0} data: samtools markdup "
                "identifies duplicates by identical outer coordinates, which for "
                "long single-molecule reads removes real coverage rather than "
                "PCR artefacts".format(platform))
        stages.append(samtools_fixmate_argv(threads=sort_threads))
        stages.append(samtools_sort_argv(None, threads=sort_threads,
                                         tmp_prefix=tmp_prefix, uncompressed=True))
        stages.append(samtools_markdup_argv(out_bam, threads=sort_threads,
                                            stats_path=markdup_stats))
    else:
        stages.append(samtools_sort_argv(out_bam, threads=sort_threads,
                                         tmp_prefix=tmp_prefix))
    return stages


# ---------------------------------------------------------------------------
# Reference and tool availability
# ---------------------------------------------------------------------------

def choose_mapper(platform: str) -> str:
    """The mapper to use, preferring bwa-mem2 and falling back to bwa.

    The fallback is a genuine equivalent — same algorithm, same output, a
    different SIMD implementation — which is the only kind of fallback this
    project allows. If neither is present the run stops and says what to install,
    because mapping Illumina reads with minimap2 instead would silently change
    every allele fraction in the report.
    """
    plat = normalise_platform(platform)
    candidates = MAPPERS[plat]
    if not candidates:
        raise MjolnirError(
            "no mapper applies to {0} input: an assembly is compared to the "
            "reference directly".format(plat))
    found = first_available(*candidates)
    if found is None:
        raise MjolnirError(
            "no {0} mapper found on PATH (looked for {1}).\n"
            "  conda install -c conda-forge -c bioconda {2}".format(
                plat, ", ".join(candidates),
                " ".join(conda_package(tool) for tool in candidates))
        )
    if found != candidates[0]:
        LOG.info("%s not found; using %s (equivalent aligner, different SIMD "
                 "implementation)", candidates[0], found)
    return found


def reference_index_paths(reference: PathLike, tool: str) -> List[Path]:
    """The files *tool*'s index consists of, in the order it writes them."""
    if tool not in INDEX_SUFFIXES:
        raise MjolnirError("no index layout known for {0!r}".format(tool))
    base = str(reference)
    return [Path(base + suffix) for suffix in INDEX_SUFFIXES[tool]]


def ensure_reference_index(reference: PathLike, tool: str, *, build: bool = False) -> Path:
    """Check that *reference* is indexed for *tool*, or say exactly how to index it.

    The two indexes are treated differently because they cost differently.

    ``samtools faidx`` on a 4.4 Mb genome takes about a second, so it is built on
    demand whenever the directory is writable. Refusing it only sends the user
    away to run a command that takes less time than reading the error did — and
    ``mjolnir db fetch`` followed by ``mjolnir run`` failed on exactly that.

    An aligner index is minutes, and a pipeline that quietly starts one in the
    middle of a batch is indistinguishable from a pipeline that has hung, so that
    one still refuses unless asked with ``--build-index``. A read-only database
    directory also refuses, with the command to run.
    """
    fasta = require_file(reference, "reference FASTA")
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        if build or os.access(str(fasta.parent), os.W_OK):
            LOG.info("indexing %s with samtools faidx", fasta.name)
            run_pipeline([samtools_faidx_argv(fasta)])
        else:
            raise MjolnirError(
                "reference {0} has no .fai index and {1} is not writable.\n"
                "  samtools faidx {0}".format(fasta, fasta.parent))

    if tool in INDEX_SUFFIXES:
        missing = [path for path in reference_index_paths(fasta, tool) if not path.exists()]
        if missing:
            if build:
                LOG.info("building %s index for %s (one-off, minutes)", tool, fasta.name)
                run_pipeline([index_reference_argv(fasta, tool)])
            else:
                raise MjolnirError(
                    "reference {0} is not indexed for {1} (missing {2}).\n"
                    "  {1} index {0}".format(
                        fasta, tool, ", ".join(path.name for path in missing))
                )
    # minimap2 needs no prebuilt index: it indexes a 4.4 Mb genome in about a
    # second, and a .mmi built with the wrong preset silently changes alignment.
    return fasta


# ---------------------------------------------------------------------------
# The wrapper that runs it
# ---------------------------------------------------------------------------

def map_reads(sample: SampleInput, reference: PathLike, out_bam: PathLike, *,
              threads: int = 1, mapper: Optional[str] = None,
              mark_duplicates: Optional[bool] = None,
              build_index: bool = False,
              config: Optional[Config] = None,
              keep_logs: bool = False) -> Alignment:
    """Align one sample's reads and return the sorted, indexed BAM.

    *mark_duplicates* defaults to True on Illumina and False on ONT; passing it
    explicitly for ONT raises rather than being ignored, so a caller that thinks
    it asked for duplicate marking finds out that it did not.
    """
    platform = normalise_platform(sample.platform)
    if platform == PLATFORM_FASTA:
        raise MjolnirError(
            "sample {0!r} is an assembly and has no reads to align".format(sample.sample_id))
    if config is not None and threads == 1:
        threads = config.threads

    tool = mapper or choose_mapper(platform)
    require(tool, "read alignment")
    require("samtools", "BAM sorting and indexing")
    fasta = ensure_reference_index(reference, tool, build=build_index)
    for path in sample.paths:
        require_file(path, "input reads for sample {0!r}".format(sample.sample_id))

    out_bam = Path(out_bam)
    ensure_dir(out_bam.parent)
    if mark_duplicates is None:
        mark_duplicates = platform == PLATFORM_ILLUMINA

    stats_path = out_bam.with_suffix(".markdup.txt") if mark_duplicates else None
    stages = mapping_pipeline(sample, fasta, out_bam, mapper=tool, threads=threads,
                              mark_duplicates=mark_duplicates,
                              tmp_prefix=out_bam.with_suffix(".sorttmp"),
                              markdup_stats=stats_path)
    LOG.info("mapping %s (%s) with %s", sample.sample_id, platform, tool)
    run_pipeline(stages, log_dir=out_bam.parent if keep_logs else None,
                 keep_logs=keep_logs)
    run_pipeline([samtools_index_argv(out_bam, threads=threads)])

    note = ""
    if not mark_duplicates:
        note = (
            "duplicates were not marked: samtools markdup identifies duplicates "
            "from identical outer coordinates, which on ONT data marks independent "
            "molecules. The duplicate fraction is therefore unmeasured, not zero."
            if platform == PLATFORM_ONT else
            "duplicate marking was disabled for this run"
        )
    return Alignment(
        sample_id=sample.sample_id,
        platform=platform,
        bam=out_bam,
        reference=Path(fasta),
        mapper=tool,
        commands=[list(stage) for stage in stages] + [samtools_index_argv(out_bam)],
        duplicates_marked=bool(mark_duplicates),
        duplicates_note=note,
    )
