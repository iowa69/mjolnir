"""Shared low-level helpers: logging, errors, files, small sequence maths.

Nothing in here knows what a resistance call is. It is the layer every other
module may import without creating a cycle, so it holds exactly two things that
matter to the design: the one error type the CLI prints without a traceback,
and the helpers that make a failure loud instead of quiet.

Two of those helpers deserve naming here, because they exist to enforce house
rules rather than to save typing:

``require`` and ``require_database`` refuse to continue when a tool or a
database is absent, and say what to install or fetch. The alternative — falling
back to some degraded path — produces a report that looks finished and is not,
which is the failure mode this whole project was written against.

``safe_fraction`` returns ``None`` when the denominator is zero, never ``0.0``.
A metric that could not be computed must stay distinguishable from a metric
that was computed and came out at zero; collapsing the two is how "unmeasured"
becomes "clean" three layers downstream.

Running external commands lives in ``shell.py``, not here, so that version
capture and command logging have one home.
"""

from __future__ import annotations

import gzip
import hashlib
import logging
import os
import re
import shutil
import sys
import tempfile
from contextlib import contextmanager
from dataclasses import is_dataclass, asdict
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union

LOG = logging.getLogger("mjolnir")

PathLike = Union[str, Path]

_LEVEL_COLOURS = {
    "DEBUG": "\033[36m",
    "INFO": "\033[32m",
    "WARNING": "\033[33m",
    "ERROR": "\033[31m",
    "CRITICAL": "\033[1;31m",
}
_RESET = "\033[0m"


class MjolnirError(RuntimeError):
    """Fatal, user-facing error. The CLI prints these without a traceback.

    Every raise site must name the thing that is missing and the command that
    would fix it. "Could not run" is not an acceptable message; "minimap2 not
    found on PATH; conda install -c bioconda minimap2" is.
    """


class _Formatter(logging.Formatter):
    def __init__(self, colour: bool) -> None:
        super().__init__("[%(asctime)s] %(levelname)s %(message)s", datefmt="%H:%M:%S")
        self.colour = colour

    def format(self, record: logging.LogRecord) -> str:
        text = super().format(record)
        if self.colour:
            prefix = _LEVEL_COLOURS.get(record.levelname, "")
            if prefix:
                text = text.replace(record.levelname, prefix + record.levelname + _RESET, 1)
        return text


def setup_logging(verbosity: int = 0, quiet: bool = False) -> None:
    """Configure the package logger. verbosity: 0=INFO, >=1=DEBUG."""
    level = logging.DEBUG if verbosity >= 1 else logging.INFO
    if quiet:
        level = logging.WARNING
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(_Formatter(colour=sys.stderr.isatty()))
    LOG.handlers[:] = [handler]
    LOG.setLevel(level)
    LOG.propagate = False


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------

#: Executables whose conda package is named differently from the command. Used
#: only to make the error message actionable — nothing dispatches on it.
CONDA_PACKAGES: Dict[str, str] = {
    "bwa-mem2": "bwa-mem2",
    "bwa": "bwa",
    "minimap2": "minimap2",
    "samtools": "samtools",
    "bcftools": "bcftools",
    "freebayes": "freebayes",
    "run_clair3.sh": "clair3",
    "kraken2": "kraken2",
    "kraken2-build": "kraken2",
    "skani": "skani",
    "mash": "mash",
    "seqtk": "seqtk",
    "fastp": "fastp",
    "trimmomatic": "trimmomatic",
    "nanoq": "nanoq",
    "bedtools": "bedtools",
    "tabix": "htslib",
    "bgzip": "htslib",
}


def conda_package(tool: str) -> str:
    """The conda package that provides *tool*."""
    return CONDA_PACKAGES.get(tool, tool)


def have(tool: str) -> bool:
    """True when *tool* is on PATH."""
    return shutil.which(tool) is not None


def which(tool: str) -> Optional[str]:
    return shutil.which(tool)


def require(tool: str, why: str = "") -> str:
    """Absolute path to *tool*, or a MjolnirError naming the install command."""
    path = shutil.which(tool)
    if path is None:
        extra = " ({0})".format(why) if why else ""
        raise MjolnirError(
            "required executable '{0}' not found on PATH{1}.\n"
            "  conda install -c conda-forge -c bioconda {2}".format(
                tool, extra, conda_package(tool))
        )
    return path


def first_available(*tools: str) -> Optional[str]:
    """The first of *tools* present on PATH, or None.

    Used where a fallback is a genuine equivalent rather than a degradation —
    ``bwa-mem2`` then ``bwa`` — and never to paper over an absent capability.
    """
    for tool in tools:
        if have(tool):
            return tool
    return None


def require_database(path: PathLike, name: str, fetch_hint: str) -> Path:
    """A database directory or file that must exist, or a MjolnirError.

    *fetch_hint* is the literal command that would obtain it, because "database
    not found" without it turns into a support question every single time.
    """
    resolved = Path(path).expanduser()
    if not resolved.exists():
        raise MjolnirError(
            "{0} not found at {1}.\n  fetch it with: {2}".format(name, resolved, fetch_hint)
        )
    return resolved


def require_file(path: PathLike, what: str, fetch_hint: str = "") -> Path:
    """An input file that must exist and be non-empty."""
    resolved = Path(path).expanduser()
    if not resolved.exists():
        hint = "\n  {0}".format(fetch_hint) if fetch_hint else ""
        raise MjolnirError("{0} not found: {1}{2}".format(what, resolved, hint))
    if resolved.is_file() and resolved.stat().st_size == 0:
        raise MjolnirError("{0} is empty: {1}".format(what, resolved))
    return resolved


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------

@contextmanager
def tempdir(prefix: str = "mjolnir.", keep: bool = False,
            parent: Optional[PathLike] = None) -> Iterator[Path]:
    path = Path(tempfile.mkdtemp(prefix=prefix, dir=str(parent) if parent else None))
    try:
        yield path
    finally:
        if keep:
            LOG.debug("keeping temp dir %s", path)
        else:
            shutil.rmtree(path, ignore_errors=True)


def ensure_dir(path: PathLike) -> Path:
    resolved = Path(path).expanduser()
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def smart_open(path: PathLike, mode: str = "rt", encoding: Optional[str] = None):
    """Open a plain, gzip, bzip2 or xz file transparently.

    Text reads replace undecodable bytes. Sequence data is ASCII but FASTA
    description lines and catalogue comment columns are free text, and losing a
    character in a comment beats aborting a run on a UnicodeDecodeError.

    *encoding* is for the files that are known not to be UTF-8. MTBseq's
    resistance list is latin-1 and raises at byte 0x98, offset 7324; naming its
    encoding recovers the comment text that ``errors="replace"`` would corrupt
    into replacement characters, and those comments carry the MIC statements the
    report prints.
    """
    name = str(path)
    text = {"errors": "replace"} if "b" not in mode else {}
    if encoding and "b" not in mode:
        text["encoding"] = encoding
    if name.endswith(".gz"):
        return gzip.open(name, mode, **text)
    if name.endswith(".bz2"):
        import bz2

        return bz2.open(name, mode, **text)
    if name.endswith(".xz"):
        import lzma

        return lzma.open(name, mode, **text)
    return open(name, mode, **text)


def is_gzip(path: PathLike) -> bool:
    try:
        with open(str(path), "rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def sha256sum(path: PathLike, chunk: int = 1 << 20) -> str:
    """Checksum of a file, for the database registry.

    §5.5 of the design: a catalogue-version mismatch between two installations
    changes calls, so the report prints the checksum of every catalogue it used.
    """
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def free_bytes(path: PathLike) -> int:
    """Bytes available on the filesystem holding *path* (0 when unknown)."""
    try:
        stats = os.statvfs(str(path))
    except (OSError, AttributeError):
        return 0
    return stats.f_bavail * stats.f_frsize


def cpu_count() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


# ---------------------------------------------------------------------------
# Names
# ---------------------------------------------------------------------------

_READ_SUFFIXES = ("_R1_001", "_R2_001", "_R1", "_R2", "_1", "_2",
                  ".R1", ".R2", ".1", ".2")
_SEQ_EXTENSIONS = (".fasta", ".fa", ".fna", ".ffn", ".fsa", ".seq",
                   ".fastq", ".fq", ".contigs", ".scaffolds")


def safe_name(text: str) -> str:
    """Filesystem/column-safe token."""
    keep: List[str] = []
    for ch in str(text):
        keep.append(ch if (ch.isalnum() or ch in "._-") else "_")
    return "".join(keep).strip("_") or "sample"


def sample_name_from_path(path: PathLike, strip_read_suffix: bool = False) -> str:
    """Derive a sample name from a file path by stripping known extensions."""
    name = Path(path).name
    for ext in (".gz", ".bz2", ".xz"):
        if name.endswith(ext):
            name = name[: -len(ext)]
    for ext in _SEQ_EXTENSIONS:
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    if strip_read_suffix:
        for suffix in _READ_SUFFIXES:
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
    return name or "sample"


def natural_key(text: str) -> Tuple:
    """Sort key that orders embedded numbers numerically.

    ``rpoB_p.Ser450Leu`` beside ``rpoB_p.Ser45Leu``: plain string ordering puts
    codon 450 before codon 45, which makes an annex table read as if it were
    unsorted.
    """
    parts: List[Tuple[int, int, str]] = []
    number = ""
    for ch in str(text):
        if ch.isdigit():
            number += ch
        else:
            if number:
                parts.append((1, int(number), ""))
                number = ""
            parts.append((0, 0, ch))
    if number:
        parts.append((1, int(number), ""))
    return tuple(parts)


# ---------------------------------------------------------------------------
# Numbers
# ---------------------------------------------------------------------------

def safe_fraction(numerator: float, denominator: float) -> Optional[float]:
    """``numerator / denominator``, or None when the denominator is zero.

    Deliberately not 0.0. A breadth of coverage that could not be computed
    because no reads mapped is a different statement from a breadth of zero,
    and only one of them may be printed as a number.
    """
    if not denominator:
        return None
    return float(numerator) / float(denominator)


def percentage(fraction: Optional[float], digits: int = 2) -> Optional[float]:
    """A 0-1 fraction as a percentage, preserving None."""
    if fraction is None:
        return None
    return round(float(fraction) * 100.0, digits)


def mean(values: Sequence[float]) -> Optional[float]:
    if not values:
        return None
    return float(sum(values)) / len(values)


def median(values: Sequence[float]) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[mid])
    return (float(ordered[mid - 1]) + float(ordered[mid])) / 2.0


def round_or_none(value: Optional[float], digits: int = 3) -> Optional[float]:
    return None if value is None else round(float(value), digits)


def human_bytes(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024.0:
            return "{0:.0f}{1}".format(n, unit) if unit == "B" else "{0:.1f}{1}".format(n, unit)
        n /= 1024.0
    return "{0:.1f}PB".format(n)


def human_time(seconds: float) -> str:
    if seconds < 60:
        return "{0:.1f}s".format(seconds)
    minutes, secs = divmod(seconds, 60)
    if minutes < 60:
        return "{0}m{1:04.1f}s".format(int(minutes), secs)
    hours, minutes = divmod(minutes, 60)
    return "{0}h{1}m".format(int(hours), int(minutes))


def chunked(items: Iterable, size: int) -> Iterator[List]:
    """Yield lists of at most *size* items."""
    batch: List = []
    for item in items:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


# ---------------------------------------------------------------------------
# Sequence
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv",
                            "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


_CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

#: One-letter to three-letter amino acid, for HGVS ``p.`` names. The design
#: fixes three-letter form as the cross-catalogue join key (§5.3), so both
#: directions are needed: tbdb writes ``p.Ser450Leu``, some legacy MTBseq rows
#: write ``S450L``.
AA_1_TO_3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter",
    "X": "Xaa",
}
AA_3_TO_1 = dict((three, one) for one, three in AA_1_TO_3.items())


def translate(seq: str, start_is_met: bool = True) -> str:
    """Translate in frame 1 with the standard code.

    Bacterial alternative starts (GTG/TTG/ATT/CTG) are rendered ``M`` when
    *start_is_met*, matching how reference protein records store them — without
    this, every *M. tuberculosis* gene beginning GTG looks like a variant.
    """
    seq = seq.upper().replace("U", "T")
    out: List[str] = []
    for i in range(0, len(seq) - 2, 3):
        out.append(_CODONS.get(seq[i:i + 3], "X"))
    if start_is_met and out:
        out[0] = "M"
    return "".join(out)


def gc_fraction(seq: str) -> Optional[float]:
    """GC as a 0-1 fraction over unambiguous bases, or None for an empty string."""
    upper = seq.upper()
    acgt = sum(upper.count(base) for base in "ACGT")
    if not acgt:
        return None
    return float(upper.count("G") + upper.count("C")) / acgt


#: A run of nucleotides long enough that it can only be sequence. 50 is well
#: above any gene name, allele identifier or accession and far below anything a
#: model could use. ``agent/observation.py`` raises on a match; the pattern
#: lives here so the report writer can apply the same test before it embeds a
#: field it did not generate.
NUCLEOTIDE_RUN = re.compile(r"[ACGTNacgtn]{50,}")


def looks_like_sequence(text: str) -> bool:
    """Whether a string carries a long nucleotide run.

    House rule: the model never sees raw sequence, and that is enforced in code
    rather than asked for in a prompt.
    """
    return bool(NUCLEOTIDE_RUN.search(text or ""))


# ---------------------------------------------------------------------------
# Serialisation
# ---------------------------------------------------------------------------

def to_jsonable(obj: Any, digits: int = 6) -> Any:
    """Recursively convert a record structure into something ``json`` accepts.

    Paths become strings, sets and tuples become lists, dataclasses become
    dicts, and non-string mapping keys are joined with ``|`` — cohort pair keys
    are ``(sample_a, sample_b)`` tuples in memory and ``"a|b"`` on disk. Floats
    are rounded so that two runs of the same data produce byte-identical JSON,
    which is what makes a golden-file test on the report's data layer possible.
    """
    if obj is None or isinstance(obj, (str, int, bool)):
        return obj
    if isinstance(obj, float):
        return round(obj, digits)
    if isinstance(obj, Path):
        return str(obj)
    if is_dataclass(obj) and not isinstance(obj, type):
        return to_jsonable(asdict(obj), digits)
    if isinstance(obj, dict):
        out: Dict[str, Any] = {}
        for key, value in obj.items():
            if isinstance(key, str):
                name = key
            elif isinstance(key, tuple):
                name = "|".join(str(part) for part in key)
            else:
                name = str(key)
            out[name] = to_jsonable(value, digits)
        return out
    if isinstance(obj, (list, tuple, set, frozenset)):
        items = sorted(obj, key=str) if isinstance(obj, (set, frozenset)) else obj
        return [to_jsonable(item, digits) for item in items]
    return str(obj)
