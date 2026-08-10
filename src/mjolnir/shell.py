"""Running external tools, and remembering exactly which ones ran.

Every mapper, caller and sketcher Mjolnir drives is somebody else's binary, and
the report has to be able to say which build of it produced the numbers. That is
the reason this module exists rather than scattered ``subprocess.run`` calls: a
command that runs through :func:`run` has its version captured on first use and
filed in a run-wide registry, so ``SampleResult.tool_versions`` is populated as a
side effect of doing the work rather than by a separate pass that can drift out
of step with it.

Two behaviours here are deliberate and are not conveniences:

A non-zero exit raises :class:`~mjolnir.utils.MjolnirError` naming the tool, the
exit code, the tail of its stderr and the conda line that installs it. The
alternative — returning a result object the caller may or may not inspect — is
how a failed alignment becomes an empty BAM becomes a sample with no variants
becomes a report that says "no resistance determinant detected". There is no
``check=False`` default; a caller that genuinely tolerates failure has to say so.

A missing binary raises before anything else happens, with the install command.
``FileNotFoundError: 'bwa-mem2'`` three stages into a run is a support ticket;
"bwa-mem2 not found on PATH ... conda install -c bioconda bwa-mem2" is not.

Version probing itself never raises. A tool that is present but whose version
cannot be read is recorded as ``"present, version not reported"`` — which is a
true statement — and never as a version string that was guessed.
"""

from __future__ import annotations

import os
import re
import shlex
import signal
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from .utils import (
    LOG,
    MjolnirError,
    PathLike,
    conda_package,
    to_jsonable,
    which,
)

# ---------------------------------------------------------------------------
# Plumbing constants. These are not scientific thresholds — they do not change
# any call — so they live here beside the code that uses them rather than in the
# threshold registry, but they are still named and justified.
# ---------------------------------------------------------------------------

#: Lines of stderr quoted back in an error message. Enough to carry a real
#: diagnostic ("[E::hts_open] fail to open file") without pasting a mapper's
#: entire progress log into a traceback.
STDERR_TAIL_LINES = 40

#: Seconds allowed for a ``--version`` probe. A version flag that has not
#: answered in ten seconds is a tool that is broken or is trying to do real work,
#: and either way the run should not stall on it.
VERSION_PROBE_TIMEOUT_SECONDS = 10

#: Bytes of a pipeline's stdout read back into memory when the caller did not
#: redirect it to a file. A pipeline is normally used to keep a SAM stream off
#: the disk entirely, so a caller that lets its stdout be captured is either
#: producing something small or has made a mistake; 8 MB is generous for the
#: first case and bounded for the second.
MAX_CAPTURED_PIPE_BYTES = 8 << 20

#: Longest version string kept. Some tools answer ``--version`` with a banner.
MAX_VERSION_CHARS = 120

#: What is recorded for a tool that ran but would not say what it was. Written
#: out in full rather than left empty, because an empty string in the report's
#: version column reads as "not used".
VERSION_UNKNOWN = "present, version not reported"


# ---------------------------------------------------------------------------
# Version probes
# ---------------------------------------------------------------------------

#: Argument vectors that make each tool print its version, in the order to try.
#: Not every tool uses ``--version``: ``bwa`` and ``freebayes`` predate the
#: convention (``bwa`` with no arguments prints a usage banner carrying
#: ``Version:`` and exits 1), ``bwa-mem2`` and ``seqkit`` use a subcommand, and
#: skani answers ``-V``. An entry of ``()`` means "run it with no arguments".
VERSION_ARGS: Dict[str, Tuple[Tuple[str, ...], ...]] = {
    "bwa": ((), ("--version",)),
    "bwa-mem2": (("version",), ("--version",)),
    "minimap2": (("--version",),),
    "samtools": (("--version",),),
    "bcftools": (("--version",),),
    "freebayes": (("--version",),),
    "run_clair3.sh": (("--version",),),
    "kraken2": (("--version",),),
    "kraken2-build": (("--version",),),
    "mash": (("--version",), ()),
    "skani": (("--version",), ("-V",)),
    "seqkit": (("version",),),
    "iqtree": (("--version",),),
    "iqtree2": (("--version",),),
    "iqtree3": (("--version",),),
    "fastp": (("--version",),),
    "nanoq": (("--version",),),
    "seqtk": ((),),
    "bedtools": (("--version",),),
    "tabix": (("--version",),),
    "bgzip": (("--version",),),
    "trimmomatic": (("-version",),),
}

#: Tried when a tool is not in :data:`VERSION_ARGS`.
DEFAULT_VERSION_ARGS: Tuple[Tuple[str, ...], ...] = (
    ("--version",), ("-version",), ("version",), ("-V",), ("-v",),
)

#: ``bwa``-style banners put the number on its own labelled line; samtools and
#: bcftools put it after the program name on line 1. Both forms are extracted so
#: that the registry stores "0.7.17-r1188" and not four lines of usage text.
_VERSION_LABEL = re.compile(r"^\s*[Vv]ersion:?\s+(\S+)", re.MULTILINE)
_VERSION_ANYWHERE = re.compile(r"\b[vV]?(\d+\.\d+(?:\.\d+)*(?:[-+.~][A-Za-z0-9._-]+)?)\b")


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

@dataclass
class CommandResult:
    """One external command that ran, and what it said.

    Carries the command itself rather than only its output: the methods annex
    prints the exact argv of every stage, which is the difference between a
    report a second lab can reproduce and one it has to guess at.
    """

    command: List[str]
    returncode: int
    stdout: str = ""
    stderr: str = ""
    duration_seconds: float = 0.0
    tool: str = ""
    version: str = ""
    timed_out: bool = False

    @property
    def ok(self) -> bool:
        return self.returncode == 0 and not self.timed_out

    def stderr_tail(self, lines: int = STDERR_TAIL_LINES) -> str:
        """The last *lines* of stderr, which is where the diagnostic usually is."""
        text = (self.stderr or "").rstrip()
        if not text:
            return ""
        parts = text.splitlines()
        return "\n".join(parts[-lines:])

    def to_dict(self) -> Dict[str, object]:
        return to_jsonable({
            "command": format_command(self.command),
            "tool": self.tool,
            "version": self.version,
            "returncode": self.returncode,
            "duration_seconds": self.duration_seconds,
            "timed_out": self.timed_out,
            "stderr_tail": self.stderr_tail(),
        })


# ---------------------------------------------------------------------------
# PATH lookup
# ---------------------------------------------------------------------------

def format_command(cmd: Sequence[object]) -> str:
    """The command as a line that can be pasted into a shell and re-run."""
    return " ".join(shlex.quote(str(part)) for part in cmd)


def tool_path(tool: str, why: str = "", required: bool = True) -> Optional[str]:
    """Resolve *tool* on PATH.

    The one lookup helper the rest of the package should use, because it is also
    where the "what do I install" answer is attached. With ``required=False`` it
    returns None instead of raising, for the places that legitimately choose
    between equivalents — ``bwa-mem2`` or ``bwa`` — rather than degrading.
    """
    if os.path.sep in str(tool):
        candidate = Path(tool).expanduser()
        if os.access(str(candidate), os.X_OK) and candidate.is_file():
            return str(candidate)
        found = None
    else:
        found = which(tool)
    if found:
        return found
    if not required:
        return None
    extra = " ({0})".format(why) if why else ""
    raise MjolnirError(
        "required executable '{0}' not found on PATH{1}.\n"
        "  conda install -c conda-forge -c bioconda {2}".format(
            tool, extra, conda_package(tool))
    )


def first_tool(*tools: str) -> Optional[str]:
    """The first of *tools* present on PATH, or None.

    Only for genuine equivalents. Where the alternatives differ in what they
    measure — bcftools versus Clair3 on ONT — the caller must record which one
    it used and state the consequence, so it picks explicitly instead.
    """
    for tool in tools:
        found = which(tool)
        if found:
            return tool
    return None


# ---------------------------------------------------------------------------
# Version capture
# ---------------------------------------------------------------------------

_VERSION_CACHE: Dict[str, str] = {}
_CAPTURED: Dict[str, str] = {}


def _probe_version(tool: str, path: str) -> str:
    """Best-effort version string. Never raises.

    Accepts a non-zero exit, because several of these tools print their banner
    to stderr and exit 1 when given no arguments — ``bwa`` is the canonical
    example — and refusing those would leave the report's most-used aligner
    without a version.
    """
    attempts = VERSION_ARGS.get(Path(tool).name, DEFAULT_VERSION_ARGS)
    for args in attempts:
        try:
            proc = subprocess.run(
                [path] + list(args),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdin=subprocess.DEVNULL,
                timeout=VERSION_PROBE_TIMEOUT_SECONDS,
                universal_newlines=True,
            )
        except (OSError, subprocess.SubprocessError):
            # A probe that will not run tells us nothing about the tool's
            # version; it does not tell us the tool is absent either, so the
            # next argument form is tried and the caller still gets a result.
            continue
        text = "{0}\n{1}".format(proc.stdout or "", proc.stderr or "").strip()
        if not text:
            continue
        parsed = parse_version(text)
        if parsed:
            return parsed
    return VERSION_UNKNOWN


def parse_version(text: str) -> str:
    """Pull a version out of a tool's banner, or return the banner's first line.

    Separated out so it can be unit-tested against real banners without running
    anything: the shapes differ enough (``Version: 0.7.17-r1188``,
    ``samtools 1.19.2``, a bare ``2.24``) that a single regex gets one of them
    wrong.
    """
    if not text:
        return ""
    label = _VERSION_LABEL.search(text)
    if label:
        return label.group(1)[:MAX_VERSION_CHARS]
    first = next((line.strip() for line in text.splitlines() if line.strip()), "")
    number = _VERSION_ANYWHERE.search(first)
    if number:
        return number.group(1)[:MAX_VERSION_CHARS]
    number = _VERSION_ANYWHERE.search(text)
    if number:
        return number.group(1)[:MAX_VERSION_CHARS]
    return first[:MAX_VERSION_CHARS]


def tool_version(tool: str, refresh: bool = False) -> Optional[str]:
    """The version of *tool*, or None when it is not on PATH.

    Cached per process: the probe is a subprocess, and mapping a cohort would
    otherwise spawn one per sample to learn something that cannot change
    mid-run.
    """
    name = Path(tool).name
    if not refresh and name in _VERSION_CACHE:
        return _VERSION_CACHE[name]
    path = which(tool) if os.path.sep not in str(tool) else str(tool)
    if not path or not os.access(path, os.X_OK):
        return None
    version = _probe_version(tool, path)
    _VERSION_CACHE[name] = version
    return version


def record_tool(tool: str) -> Optional[str]:
    """Capture *tool*'s version and file it for the report's methods annex."""
    name = Path(tool).name
    version = tool_version(tool)
    if version is not None:
        _CAPTURED[name] = version
    return version


def captured_versions() -> Dict[str, str]:
    """Every tool that has actually run in this process, with its version.

    This is what goes into ``SampleResult.tool_versions``. It lists what ran,
    not what is installed — ``doctor`` answers the second question — so a report
    never claims to have used a caller it never invoked.
    """
    return dict(sorted(_CAPTURED.items()))


def record_versions(*tools: str) -> Dict[str, str]:
    """Capture several tools up front, for a stage that shells out repeatedly."""
    for tool in tools:
        record_tool(tool)
    return captured_versions()


def reset_captured() -> None:
    """Forget which tools have run. For tests, and for a second sample in-process."""
    _CAPTURED.clear()


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------

def _merged_env(extra_env: Optional[Mapping[str, str]]) -> Optional[Dict[str, str]]:
    if not extra_env:
        return None
    merged = dict(os.environ)
    merged.update(dict((str(k), str(v)) for k, v in extra_env.items()))
    return merged


def _fail(result: CommandResult, why: str) -> MjolnirError:
    tail = result.stderr_tail()
    detail = "\n  stderr:\n{0}".format(_indent(tail)) if tail else ""
    context = "\n  needed for: {0}".format(why) if why else ""
    if result.timed_out:
        headline = "{0} timed out".format(result.tool or result.command[0])
    else:
        headline = "{0} failed with exit code {1}".format(
            result.tool or result.command[0], result.returncode)
    return MjolnirError(
        "{0}{1}\n  command: {2}{3}\n"
        "  if {4} is not installed: conda install -c conda-forge -c bioconda {5}".format(
            headline, context, format_command(result.command), detail,
            result.tool or result.command[0],
            conda_package(Path(result.command[0]).name))
    )


def _indent(text: str, prefix: str = "    ") -> str:
    return "\n".join(prefix + line for line in text.splitlines())


def run(cmd: Sequence[object], *,
        check: bool = True,
        capture: bool = True,
        stdout_path: Optional[PathLike] = None,
        stdin_data: Optional[str] = None,
        cwd: Optional[PathLike] = None,
        extra_env: Optional[Mapping[str, str]] = None,
        timeout: Optional[float] = None,
        why: str = "") -> CommandResult:
    """Run one external command and capture what it produced.

    *why* is folded into the error message when the command fails, so a failure
    reads "samtools failed ... needed for: coverage depth over the reference"
    rather than making the operator infer the stage from an argv.

    ``stdout_path`` sends stdout straight to a file, which is how the alignment
    and pileup stages avoid holding a SAM stream in memory; ``stdout`` on the
    returned result is then empty by construction, not because nothing was
    written.
    """
    argv = [str(part) for part in cmd]
    if not argv:
        raise MjolnirError("run() was given an empty command")
    tool = Path(argv[0]).name

    # Resolve first, so an absent binary is an install instruction rather than a
    # FileNotFoundError from inside subprocess.
    resolved = tool_path(argv[0], why=why)
    argv[0] = resolved or argv[0]
    version = record_tool(argv[0])

    LOG.debug("exec: %s", format_command(argv))
    started = time.time()
    stdout_handle = None
    try:
        if stdout_path is not None:
            try:
                stdout_handle = open(str(stdout_path), "w")
            except OSError as exc:
                raise MjolnirError(
                    "cannot write {0}'s output to {1}: {2}".format(
                        tool, stdout_path, exc)) from exc
            stdout_target = stdout_handle
        else:
            stdout_target = subprocess.PIPE if capture else None
        try:
            proc = subprocess.Popen(
                argv,
                stdin=subprocess.PIPE if stdin_data is not None else subprocess.DEVNULL,
                stdout=stdout_target,
                stderr=subprocess.PIPE if capture else None,
                cwd=str(cwd) if cwd else None,
                env=_merged_env(extra_env),
                universal_newlines=True,
                start_new_session=True,
            )
        except OSError as exc:
            raise MjolnirError(
                "could not start {0}: {1}\n  command: {2}\n"
                "  conda install -c conda-forge -c bioconda {3}".format(
                    tool, exc, format_command(argv), conda_package(tool))
            ) from exc

        timed_out = False
        try:
            out, err = proc.communicate(input=stdin_data, timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            _kill(proc)
            out, err = proc.communicate()
            err = (err or "") + "\ntimed out after {0}s".format(timeout)
    finally:
        if stdout_handle is not None:
            stdout_handle.close()

    result = CommandResult(
        command=argv,
        returncode=proc.returncode if not timed_out else -1,
        stdout=out or "",
        stderr=err or "",
        duration_seconds=time.time() - started,
        tool=tool,
        version=version or "",
        timed_out=timed_out,
    )
    if check and not result.ok:
        raise _fail(result, why)
    return result


def _kill(proc: "subprocess.Popen") -> None:
    """Kill a timed-out process and anything it spawned.

    ``start_new_session`` puts each command in its own process group precisely so
    this can reach the children. Killing only the parent leaves a mapper's
    threads running against the same output file.
    """
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (OSError, ProcessLookupError):
        try:
            proc.kill()
        except OSError:
            pass


@dataclass
class PipelineResult:
    """A shell pipeline: every stage's exit status, not just the last one's.

    ``bwa-mem2 mem ... | samtools sort`` returns 0 from a shell when the aligner
    dies and sort writes an empty BAM, because a shell reports only the last
    exit status. Every stage is checked here for exactly that reason.
    """

    stages: List[CommandResult] = field(default_factory=list)
    stdout: str = ""
    duration_seconds: float = 0.0

    @property
    def ok(self) -> bool:
        return all(stage.ok for stage in self.stages)

    @property
    def failed(self) -> Optional[CommandResult]:
        for stage in self.stages:
            if not stage.ok:
                return stage
        return None

    def to_dict(self) -> Dict[str, object]:
        return to_jsonable({
            "stages": [stage.to_dict() for stage in self.stages],
            "duration_seconds": self.duration_seconds,
        })


def run_pipe(stages: Sequence[Sequence[object]], *,
             check: bool = True,
             stdout_path: Optional[PathLike] = None,
             cwd: Optional[PathLike] = None,
             extra_env: Optional[Mapping[str, str]] = None,
             timeout: Optional[float] = None,
             why: str = "") -> PipelineResult:
    """Run commands joined by pipes, checking every stage.

    Each stage's stderr goes to its own temporary file rather than a pipe: with
    pipes, a stage that writes more stderr than the pipe buffer holds blocks
    forever while the parent waits on the last stage, and a mapper's progress
    log is easily that large. Deadlocking a run to save a temp file is not a
    trade worth making.
    """
    if not stages:
        raise MjolnirError("run_pipe() was given no stages")

    argvs: List[List[str]] = []
    for stage in stages:
        argv = [str(part) for part in stage]
        if not argv:
            raise MjolnirError("run_pipe() was given an empty stage")
        argv[0] = tool_path(argv[0], why=why) or argv[0]
        record_tool(argv[0])
        argvs.append(argv)

    LOG.debug("exec pipeline: %s", " | ".join(format_command(a) for a in argvs))
    started = time.time()
    env = _merged_env(extra_env)
    err_files = [tempfile.TemporaryFile(mode="w+") for _ in argvs]
    procs: List[subprocess.Popen] = []
    out_handle = None
    captured = ""
    try:
        if stdout_path is not None:
            try:
                out_handle = open(str(stdout_path), "w")
            except OSError as exc:
                raise MjolnirError(
                    "cannot write pipeline output to {0}: {1}".format(
                        stdout_path, exc)) from exc
            final_stdout = out_handle
        else:
            final_stdout = tempfile.TemporaryFile(mode="w+")

        previous_stdout = subprocess.DEVNULL
        for index, argv in enumerate(argvs):
            last = index == len(argvs) - 1
            try:
                proc = subprocess.Popen(
                    argv,
                    stdin=previous_stdout,
                    stdout=final_stdout if last else subprocess.PIPE,
                    stderr=err_files[index],
                    cwd=str(cwd) if cwd else None,
                    env=env,
                    universal_newlines=True,
                    start_new_session=True,
                )
            except OSError as exc:
                for running in procs:
                    _kill(running)
                for handle in err_files:
                    handle.close()
                if out_handle is None:
                    final_stdout.close()
                raise MjolnirError(
                    "could not start {0}: {1}\n  pipeline: {2}\n"
                    "  conda install -c conda-forge -c bioconda {3}".format(
                        Path(argv[0]).name, exc,
                        " | ".join(format_command(a) for a in argvs),
                        conda_package(Path(argv[0]).name))
                ) from exc
            # Close our copy of the upstream pipe so the upstream stage sees
            # EPIPE when this one exits early.
            if previous_stdout not in (subprocess.DEVNULL, None):
                previous_stdout.close()
            previous_stdout = proc.stdout
            procs.append(proc)

        timed_out = False
        deadline = None if timeout is None else started + timeout
        for proc in reversed(procs):
            remaining = None if deadline is None else max(0.1, deadline - time.time())
            try:
                proc.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                timed_out = True
                for running in procs:
                    _kill(running)
                break

        if out_handle is None:
            final_stdout.seek(0)
            captured = final_stdout.read(MAX_CAPTURED_PIPE_BYTES)
            final_stdout.close()
    finally:
        if out_handle is not None:
            out_handle.close()

    results: List[CommandResult] = []
    for index, (argv, proc) in enumerate(zip(argvs, procs)):
        err_files[index].seek(0)
        stderr = err_files[index].read()
        err_files[index].close()
        name = Path(argv[0]).name
        results.append(CommandResult(
            command=argv,
            returncode=proc.returncode if proc.returncode is not None else -1,
            stdout="",
            stderr=stderr,
            duration_seconds=time.time() - started,
            tool=name,
            version=_VERSION_CACHE.get(name, ""),
            timed_out=timed_out,
        ))

    pipeline = PipelineResult(stages=results, stdout=captured,
                              duration_seconds=time.time() - started)
    if check and not pipeline.ok:
        failed = pipeline.failed or results[-1]
        raise _fail(failed, why or "pipeline: {0}".format(
            " | ".join(format_command(a) for a in argvs)))
    return pipeline


def run_text(cmd: Sequence[object], *, why: str = "",
             timeout: Optional[float] = None) -> str:
    """Run a command and return its stdout. Convenience over :func:`run`."""
    return run(cmd, why=why, timeout=timeout).stdout
