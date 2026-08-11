"""Working out what the user actually handed us, and saying how we know.

Three input shapes reach this module — a pair of Illumina FASTQs, a single ONT
FASTQ, an assembled FASTA — and getting the answer wrong is not a cosmetic
error. Calling ONT reads "Illumina" applies a 3-read variant threshold instead
of 5, drops the whole ONT caveat block from the report, and lets a spurious
``fbiC`` deletion through as delamanid resistance. So detection here is
evidence-producing rather than boolean: :func:`detect_platform` returns the read
lengths it measured, the filename hints it saw, the conflicts between them and
the confidence that follows, and every caller is expected to put that in the
report rather than just the label.

The platform call is made from read lengths, not from filenames. Filenames are
advisory and are frequently wrong — a MinKNOW run copied into a directory called
``illumina_backup`` is a real thing — so a hint that contradicts the measurement
lowers confidence and is recorded as a conflict; it never flips the call. The
one case where a hint decides anything is a file too thin to measure, and then
the evidence says so in as many words.

Nothing here reads a whole file. Length profiling stops at a bounded number of
reads *and* a bounded number of bases, because an ultra-long ONT library would
otherwise pull hundreds of megabytes through gzip to answer a question that a
few hundred reads settle.

Compression is decided by magic bytes, not by the extension. A ``.fastq`` that
is really gzip and a ``.gz`` that is really plain text both exist in the wild,
and reading the first as text produces mojibake that looks like a malformed
FASTQ rather than a misnamed file.
"""

from __future__ import annotations

import bz2
import gzip
import lzma
import re
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple

from . import config
from .records import (
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MODERATE,
    CONFIDENCE_NONE,
    PLATFORM_FASTA,
    PLATFORM_ILLUMINA,
    PLATFORM_ONT,
    STATUS_FAIL,
    STATUS_PASS,
    STATUS_WARN,
    Check,
    SampleInput,
    normalise_platform,
)
from .utils import (
    LOG,
    MjolnirError,
    PathLike,
    natural_key,
    safe_fraction,
    to_jsonable,
    safe_name,
)

# ---------------------------------------------------------------------------
# Thresholds
#
# House rule: no bare number in logic, and every threshold names its source.
# These belong to input detection, so they are declared beside the code that
# applies them — but they are registered in config.py's one registry, which is
# what ``source_for()`` reads and what the report prints. A number registered
# anywhere else would appear in the PDF with no attribution.
# ---------------------------------------------------------------------------

SRC_ILLUMINA_CHEMISTRY = (
    "Illumina sequencing chemistry - the longest routine paired-end read length "
    "offered on any Illumina instrument is 2x300 nt (MiSeq Reagent Kit v3)"
)
SRC_BCL2FASTQ = (
    "Illumina bcl2fastq2 Conversion Software v2.20 Software Guide (document "
    "15051736) - FASTQ naming: SampleName_S#_L00#_R#_001.fastq.gz"
)
SRC_FASTQ_FORMAT = (
    "Cock et al. 2010, Nucleic Acids Res 38:1767 - the Sanger FASTQ format and "
    "the Solexa/Illumina quality-encoding variants"
)


def _register(name: str, value: Any, source: str, unit: str = "", note: str = "",
              verified: bool = True) -> Any:
    """File a seqio threshold in config.py's registry and return its value."""
    return config._define(name, value, source, unit=unit, note=note, verified=verified)


#: SOURCE: Illumina chemistry (above). Above this, no Illumina instrument
#: produced the read. Used as the single hard discriminator between platforms,
#: with headroom over 300 nt so that a 301-nt read with its adapter still
#: retained does not read as ONT.
ILLUMINA_MAX_READ_LENGTH = _register(
    "seqio_illumina_max_read_length", 400, SRC_ILLUMINA_CHEMISTRY, unit="nt",
    note="longest read length attributable to Illumina chemistry, with headroom")

#: SOURCE: Mjolnir policy. A handful of over-length reads in an Illumina file
#: are concatemers and adapter artefacts, not evidence of a long-read platform;
#: a tenth of the library being over-length is not an artefact.
ONT_LONG_READ_FRACTION = _register(
    "seqio_ont_long_read_fraction", 0.10, config.SRC_POLICY, unit="fraction",
    note="fraction of sampled reads above the Illumina length ceiling that "
         "makes the library long-read")

#: SOURCE: Mjolnir policy, consistent with the design's ONT configuration floor
#: (R10.4.1 + Dorado sup). A median above 1 kb is unambiguous long-read data and
#: raises the platform call from moderate to high confidence; a lower median
#: does not argue against ONT, since amplicon and degraded libraries are short.
ONT_INDICATIVE_MEDIAN_LENGTH = _register(
    "seqio_ont_indicative_median_length", 1000, config.SRC_POLICY, unit="nt",
    note="median read length above which a long-read call is high confidence")

#: SOURCE: Mjolnir policy. Reads sampled per file for the platform call.
#: Bounded because this runs before anything else and must not become the
#: expensive step; 2,000 reads settle a length distribution that differs by an
#: order of magnitude between the platforms.
PLATFORM_SAMPLE_READS = _register(
    "seqio_platform_sample_reads", 2000, config.SRC_POLICY, unit="reads",
    note="reads examined per FASTQ when detecting the platform")

#: SOURCE: Mjolnir policy. Second bound on the same sampling, in bases, so that
#: an ultra-long ONT library stops early rather than decompressing 200 MB.
PLATFORM_SAMPLE_BASES = _register(
    "seqio_platform_sample_bases", 20_000_000, config.SRC_POLICY, unit="bases",
    note="bases examined per FASTQ when detecting the platform")

#: SOURCE: Mjolnir policy. Below this many reads the length distribution is not
#: a measurement, and the platform call is reported as low confidence resting on
#: the filename rather than on the data.
MIN_READS_FOR_PLATFORM_CALL = _register(
    "seqio_min_reads_for_platform_call", 50, config.SRC_POLICY, unit="reads",
    note="reads needed before the length distribution decides the platform")

#: SOURCE: Cock et al. 2010. Sanger/Illumina-1.8+ quality characters occupy
#: ASCII 33-74; the legacy Illumina-1.3+ offset-64 encodings occupy 64-104, and
#: 59-104 once Solexa's negative scores are included. A file whose quality
#: characters never drop below 59 and reach above 74 is offset-64, which every
#: base-quality threshold in this tool would otherwise misread by 31 points.
PHRED64_MIN_ASCII = _register(
    "seqio_phred64_min_ascii", 59, SRC_FASTQ_FORMAT, unit="ASCII",
    note="lowest quality character in the Solexa/Illumina-1.3+ offset-64 range")
PHRED33_MAX_ASCII = _register(
    "seqio_phred33_max_ascii", 74, SRC_FASTQ_FORMAT, unit="ASCII",
    note="highest quality character in the Sanger/Illumina-1.8+ offset-33 range")

#: SOURCE: Mjolnir policy. Quality characters sampled before pronouncing on the
#: encoding. Low-quality bases are common enough that a few hundred reads settle
#: it, and an encoding call is only used to warn, never to rescale.
QUALITY_ENCODING_SAMPLE_READS = _register(
    "seqio_quality_encoding_sample_reads", 500, config.SRC_POLICY, unit="reads",
    note="reads examined when sniffing the FASTQ quality encoding")


# ---------------------------------------------------------------------------
# Extensions, magic bytes and filename hints
# ---------------------------------------------------------------------------

FASTQ_EXTENSIONS: Tuple[str, ...] = (".fastq", ".fq")
FASTA_EXTENSIONS: Tuple[str, ...] = (".fasta", ".fa", ".fna", ".ffn", ".fsa",
                                     ".fas", ".seq", ".contigs", ".scaffolds")
COMPRESSION_EXTENSIONS: Tuple[str, ...] = (".gz", ".bz2", ".xz", ".bgz")

#: Magic bytes, longest first so that xz is not shadowed.
_MAGIC: Tuple[Tuple[bytes, str], ...] = (
    (b"\xfd7zXZ\x00", "xz"),
    (b"\x1f\x8b", "gzip"),
    (b"BZh", "bzip2"),
)

#: Substrings that suggest a long-read origin, matched against the file name and
#: its immediate parent directory. MinKNOW writes into ``fastq_pass/barcodeNN``,
#: which is why the parent is looked at at all.
ONT_FILENAME_HINTS: Tuple[str, ...] = (
    "ont", "nanopore", "minion", "gridion", "promethion", "flongle",
    "dorado", "guppy", "fastq_pass", "fastq_fail", "barcode", "minknow",
)

#: Substrings that suggest Illumina. ``_r1_001`` is the bcl2fastq stamp and is
#: the strongest of these by some distance.
ILLUMINA_FILENAME_HINTS: Tuple[str, ...] = (
    "illumina", "miseq", "nextseq", "novaseq", "hiseq", "iseq", "nextera",
    "bcl2fastq", "_r1_001", "_r2_001",
)


# ---------------------------------------------------------------------------
# Compression and opening
# ---------------------------------------------------------------------------

def compression(path: PathLike) -> str:
    """Compression of *path* by magic bytes: gzip, bzip2, xz or none.

    Deliberately ignores the extension. Sniffing costs one 6-byte read and
    catches the two mistakes that would otherwise surface much later as
    unreadable sequence: a gzip file named ``.fastq``, and a plain file named
    ``.fastq.gz``.
    """
    try:
        with open(str(path), "rb") as handle:
            head = handle.read(6)
    except OSError as exc:
        raise MjolnirError("cannot read {0}: {1}".format(path, exc)) from exc
    for magic, name in _MAGIC:
        if head.startswith(magic):
            return name
    return "none"


def open_text(path: PathLike) -> Any:
    """Open a sequence file as text, decompressing by magic bytes.

    Undecodable bytes are replaced rather than fatal: FASTA description lines
    are free text and occasionally carry a stray byte, and losing one character
    of a description is not a reason to abandon a run.
    """
    kind = compression(path)
    name = str(path)
    try:
        if kind == "gzip":
            return gzip.open(name, "rt", errors="replace")
        if kind == "bzip2":
            return bz2.open(name, "rt", errors="replace")
        if kind == "xz":
            return lzma.open(name, "rt", errors="replace")
        return open(name, "rt", errors="replace")
    except OSError as exc:
        raise MjolnirError("cannot open {0}: {1}".format(path, exc)) from exc


def compression_mismatch(path: PathLike) -> str:
    """A description of a name/content disagreement about compression, or "".

    Reported rather than corrected. Mjolnir reads the file correctly either way,
    but a lab whose ``.fastq.gz`` files are not gzipped has a broken step
    upstream and should be told.
    """
    name = str(path).lower()
    actual = compression(path)
    claimed = ""
    if name.endswith((".gz", ".bgz")):
        claimed = "gzip"
    elif name.endswith(".bz2"):
        claimed = "bzip2"
    elif name.endswith(".xz"):
        claimed = "xz"
    if claimed and actual != claimed:
        return "{0} is named as {1} but its first bytes are {2}".format(
            Path(path).name, claimed, actual)
    if not claimed and actual != "none":
        return "{0} has no compression suffix but is {1}-compressed".format(
            Path(path).name, actual)
    return ""


def strip_extensions(path: PathLike) -> str:
    """The file name with its compression and sequence extensions removed."""
    name = Path(path).name
    lowered = name.lower()
    for ext in COMPRESSION_EXTENSIONS:
        if lowered.endswith(ext):
            name = name[: -len(ext)]
            lowered = name.lower()
            break
    for ext in FASTQ_EXTENSIONS + FASTA_EXTENSIONS:
        if lowered.endswith(ext):
            name = name[: -len(ext)]
            break
    return name


# ---------------------------------------------------------------------------
# Format sniffing
# ---------------------------------------------------------------------------

FORMAT_FASTQ = "fastq"
FORMAT_FASTA = "fasta"


def sniff_format(path: PathLike) -> str:
    """"fastq", "fasta", or "" when the first record line settles neither.

    Content wins over the extension, always. A ``.fasta`` holding reads and a
    ``.fastq`` holding contigs both happen, and the consequence of believing the
    name is that an assembly gets mapped as if it were a read set.
    """
    try:
        with open_text(path) as handle:
            for _ in range(100):
                line = handle.readline()
                if not line:
                    break
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith(">"):
                    return FORMAT_FASTA
                if stripped.startswith("@"):
                    return FORMAT_FASTQ
                if stripped.startswith(";"):  # legacy FASTA comment line
                    continue
                return ""
    except (OSError, EOFError, gzip.BadGzipFile, lzma.LZMAError) as exc:
        raise MjolnirError(
            "{0} could not be read: {1}\n"
            "  the file is either truncated or not the format its name "
            "claims".format(path, exc)) from exc
    return ""


def is_fastq(path: PathLike) -> bool:
    return sniff_format(path) == FORMAT_FASTQ


def is_fasta(path: PathLike) -> bool:
    return sniff_format(path) == FORMAT_FASTA


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def read_fastq_head(path: PathLike, max_records: int = PLATFORM_SAMPLE_READS,
                    max_bases: int = PLATFORM_SAMPLE_BASES
                    ) -> Iterator[Tuple[str, str, str]]:
    """Yield up to *max_records* ``(name, sequence, quality)`` from the top.

    Both bounds are honoured: an ultra-long library hits *max_bases* first and
    stops there, which is the whole point of profiling from the head rather than
    from the file.
    """
    bases = 0
    count = 0
    with open_text(path) as handle:
        while count < max_records and bases < max_bases:
            header = handle.readline()
            if not header:
                return
            if not header.startswith("@"):
                raise MjolnirError(
                    "{0} is not a valid FASTQ: record {1} starts with {2!r}, "
                    "not '@'".format(path, count + 1, header[:1]))
            seq = handle.readline().strip()
            plus = handle.readline()
            qual = handle.readline().strip()
            if not plus or not qual:
                # A truncated final record. Everything before it was read
                # cleanly, so the profile stands on what was complete.
                return
            if len(seq) != len(qual):
                raise MjolnirError(
                    "{0} is malformed at record {1}: sequence is {2} nt but "
                    "quality is {3} characters".format(
                        path, count + 1, len(seq), len(qual)))
            yield header[1:].strip(), seq, qual
            bases += len(seq)
            count += 1


def read_fasta(path: PathLike) -> Iterator[Tuple[str, str]]:
    """Yield ``(header, sequence)`` pairs; the header keeps its description."""
    header: Optional[str] = None
    chunks: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------

@dataclass
class LengthProfile:
    """What a bounded sample of a FASTQ's reads looks like.

    ``truncated`` records that the sample hit one of its bounds, so a consumer
    can tell "these are all the reads there are" from "these are the first two
    thousand". ``reads_sampled`` is never presented as a read count for the
    file; estimating that from a compressed head is guesswork and this module
    does not do guesswork.
    """

    reads_sampled: int = 0
    bases_sampled: int = 0
    min_length: Optional[int] = None
    max_length: Optional[int] = None
    mean_length: Optional[float] = None
    median_length: Optional[float] = None
    n50: Optional[int] = None
    #: Fraction of sampled reads longer than the Illumina ceiling.
    long_fraction: Optional[float] = None
    #: Standard deviation over mean. Illumina libraries are near-uniform in
    #: length until they are trimmed; ONT libraries never are.
    length_variability: Optional[float] = None
    mean_quality: Optional[float] = None
    quality_encoding: str = ""
    truncated: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable({
            "reads_sampled": self.reads_sampled,
            "bases_sampled": self.bases_sampled,
            "min_length": self.min_length,
            "max_length": self.max_length,
            "mean_length": self.mean_length,
            "median_length": self.median_length,
            "n50": self.n50,
            "long_fraction": self.long_fraction,
            "length_variability": self.length_variability,
            "mean_quality": self.mean_quality,
            "quality_encoding": self.quality_encoding,
            "truncated": self.truncated,
        })


def _n50(lengths: Sequence[int]) -> Optional[int]:
    if not lengths:
        return None
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2.0
    running = 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return ordered[-1]


def quality_encoding(qualities: Sequence[str]) -> str:
    """"phred33", "phred64" or "" from observed quality characters.

    Only ever used to warn. Rescaling a mis-encoded file silently is how a
    sample ends up with plausible-looking base qualities that are 31 points
    wrong, so Mjolnir reports the encoding and refuses to reinterpret it.
    """
    lowest = None
    highest = None
    for qual in qualities:
        for ch in qual:
            value = ord(ch)
            if lowest is None or value < lowest:
                lowest = value
            if highest is None or value > highest:
                highest = value
    if lowest is None or highest is None:
        return ""
    if lowest < PHRED64_MIN_ASCII:
        return "phred33"
    if highest > PHRED33_MAX_ASCII:
        return "phred64"
    # Everything seen sits in the overlap of the two ranges. Modern data is
    # offset-33, but "probably" is not a measurement, so this says nothing.
    return ""


def length_profile(path: PathLike, max_records: int = PLATFORM_SAMPLE_READS,
                   max_bases: int = PLATFORM_SAMPLE_BASES) -> LengthProfile:
    """Profile the head of a FASTQ. Never reads the whole file."""
    lengths: List[int] = []
    quals: List[str] = []
    qual_sum = 0
    qual_bases = 0
    for _name, seq, qual in read_fastq_head(path, max_records, max_bases):
        lengths.append(len(seq))
        if len(quals) < QUALITY_ENCODING_SAMPLE_READS:
            quals.append(qual)
        qual_sum += sum(ord(ch) - 33 for ch in qual)
        qual_bases += len(qual)

    profile = LengthProfile(reads_sampled=len(lengths), bases_sampled=sum(lengths))
    if not lengths:
        return profile
    profile.min_length = min(lengths)
    profile.max_length = max(lengths)
    profile.mean_length = round(float(sum(lengths)) / len(lengths), 1)
    profile.median_length = float(statistics.median(lengths))
    profile.n50 = _n50(lengths)
    long_reads = sum(1 for length in lengths if length > ILLUMINA_MAX_READ_LENGTH)
    profile.long_fraction = safe_fraction(long_reads, len(lengths))
    if profile.mean_length:
        profile.length_variability = safe_fraction(
            statistics.pstdev(lengths), profile.mean_length)
    profile.mean_quality = round(qual_sum / qual_bases, 1) if qual_bases else None
    profile.quality_encoding = quality_encoding(quals)
    profile.truncated = (len(lengths) >= max_records
                         or profile.bases_sampled >= max_bases)
    return profile


@dataclass
class AssemblyStats:
    """Contiguity and composition of an assembly, for the FASTA input path."""

    contigs: int = 0
    total_length: int = 0
    n50: Optional[int] = None
    largest_contig: Optional[int] = None
    gc_fraction: Optional[float] = None
    ambiguous_bases: int = 0
    #: Fraction of bases that are unambiguous ACGT. Feeds QCMetrics; an assembly
    #: padded with N is not a clean assembly and the report says so.
    unambiguous_fraction: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable({
            "contigs": self.contigs,
            "total_length": self.total_length,
            "n50": self.n50,
            "largest_contig": self.largest_contig,
            "gc_fraction": self.gc_fraction,
            "ambiguous_bases": self.ambiguous_bases,
            "unambiguous_fraction": self.unambiguous_fraction,
        })


def assembly_stats(path: PathLike) -> AssemblyStats:
    """Read a FASTA once and summarise it.

    This is the only place a whole file is read, and it is unavoidable: an
    assembly's length and N50 cannot be sampled. Assemblies are megabytes, not
    gigabytes, so the cost is bounded by the input rather than by a choice made
    here.
    """
    stats = AssemblyStats()
    lengths: List[int] = []
    gc = 0
    acgt = 0
    for _header, seq in read_fasta(path):
        if not seq:
            continue
        upper = seq.upper()
        lengths.append(len(seq))
        stats.total_length += len(seq)
        counts = dict((base, upper.count(base)) for base in "ACGT")
        acgt += sum(counts.values())
        gc += counts["G"] + counts["C"]
        stats.ambiguous_bases += len(seq) - sum(counts.values())
    stats.contigs = len(lengths)
    if not lengths:
        return stats
    stats.largest_contig = max(lengths)
    stats.n50 = _n50(lengths)
    stats.gc_fraction = safe_fraction(gc, acgt)
    stats.unambiguous_fraction = safe_fraction(acgt, stats.total_length)
    return stats


# ---------------------------------------------------------------------------
# Filenames: pairing and sample naming
# ---------------------------------------------------------------------------

#: SOURCE: bcl2fastq2 v2.20 guide. ``SampleName_S#_L00#_R#_001.fastq.gz``, with
#: the lane token optional when the run was demultiplexed with ``--no-lane-
#: splitting``. Matched first and exactly, because it is the only convention
#: that tells us unambiguously which part of the name is the sample: given
#: ``226-18_S8_L001_R1_001``, ``_S8`` is a sample-sheet index and ``226-18`` is
#: what the lab calls the isolate.
_BCL2FASTQ = re.compile(
    r"^(?P<sample>.+?)_S(?P<index>\d+)"
    r"(?:_L(?P<lane>\d{3}))?"
    r"_(?P<read>[RI][12])"
    r"_(?P<chunk>\d{3})$"
)

#: ``sample_R1`` / ``sample.R2`` / ``sample-r1``. Unambiguous enough to strip on
#: sight: an ``R`` immediately before the read digit is not an accident.
_R_SUFFIX = re.compile(r"^(?P<sample>.+?)[._-](?P<read>[Rr][12])$")

#: ``sample_1`` / ``sample.2``. Ambiguous — an isolate genuinely called
#: ``patient_2`` looks identical — so this form is only ever believed when the
#: partner file is present in the same set. See :func:`group_reads`.
_BARE_SUFFIX = re.compile(r"^(?P<sample>.+?)[._-](?P<read>[12])$")

CONVENTION_BCL2FASTQ = "bcl2fastq"
CONVENTION_R = "R1/R2"
CONVENTION_BARE = "_1/_2"
CONVENTION_NONE = "none"


@dataclass(frozen=True)
class ReadFileName:
    """A FASTQ path taken apart into sample, read number and lane."""

    path: Path
    stem: str
    sample: str
    read: Optional[int] = None
    index_read: bool = False
    lane: str = ""
    convention: str = CONVENTION_NONE

    @property
    def ambiguous(self) -> bool:
        """Whether the read number rests on a bare ``_1``/``_2`` suffix."""
        return self.convention == CONVENTION_BARE


def parse_read_name(path: PathLike) -> ReadFileName:
    """Split a FASTQ filename into sample, read number and lane.

    Handles the four conventions that actually occur in the wild, in decreasing
    order of how much they tell us. The bcl2fastq form is the informative one:
    it is the only one that identifies the ``_S8`` sample-sheet index as
    decoration rather than as part of the name.
    """
    resolved = Path(path)
    stem = strip_extensions(resolved)

    match = _BCL2FASTQ.match(stem)
    if match:
        read_token = match.group("read")
        return ReadFileName(
            path=resolved,
            stem=stem,
            sample=match.group("sample"),
            read=int(read_token[1]),
            index_read=read_token[0].upper() == "I",
            lane=match.group("lane") or "",
            convention=CONVENTION_BCL2FASTQ,
        )

    match = _R_SUFFIX.match(stem)
    if match:
        return ReadFileName(path=resolved, stem=stem, sample=match.group("sample"),
                            read=int(match.group("read")[1]),
                            convention=CONVENTION_R)

    match = _BARE_SUFFIX.match(stem)
    if match:
        return ReadFileName(path=resolved, stem=stem, sample=match.group("sample"),
                            read=int(match.group("read")),
                            convention=CONVENTION_BARE)

    return ReadFileName(path=resolved, stem=stem, sample=stem,
                        convention=CONVENTION_NONE)


def sample_name(path: PathLike) -> str:
    """The sample name a FASTQ or FASTA path implies.

    ``30-20_S1_R1_001.fastq.gz`` is sample ``30-20``;
    ``226-18_S8_L001_R1_001.fastq.gz`` is sample ``226-18``. Both are real
    filenames from the data this tool was written for, and both are wrong under
    a naive "strip after the last underscore" rule.
    """
    parsed = parse_read_name(path)
    if parsed.convention == CONVENTION_BARE:
        # Not stripped on its own evidence; see _BARE_SUFFIX.
        return parsed.stem
    return parsed.sample


@dataclass
class ReadGroup:
    """One sample's read files, after pairing."""

    sample_id: str
    files: List[Path] = field(default_factory=list)
    paired: bool = False
    convention: str = CONVENTION_NONE
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable({
            "sample_id": self.sample_id,
            "files": [str(p) for p in self.files],
            "paired": self.paired,
            "convention": self.convention,
            "notes": list(self.notes),
        })


def group_reads(paths: Sequence[PathLike]) -> List[ReadGroup]:
    """Pair FASTQ paths into per-sample groups.

    Two rules keep this honest. A read number is only believed when its partner
    is present: a lone ``patient_2.fastq.gz`` is one single-end file for sample
    ``patient_2``, not read 2 of a pair that lost its mate, because the second
    reading would silently rename the sample. And a sample split across lanes is
    an error rather than a guess — merging lanes is a decision with consequences
    for read groups and duplicate marking, so the user is asked to concatenate
    rather than having it done behind their back.
    """
    parsed = [parse_read_name(p) for p in paths]

    index_reads = [p for p in parsed if p.index_read]
    parsed = [p for p in parsed if not p.index_read]

    buckets: Dict[str, List[ReadFileName]] = {}
    for entry in parsed:
        buckets.setdefault(entry.sample, []).append(entry)

    groups: List[ReadGroup] = []
    for name in sorted(buckets, key=natural_key):
        entries = buckets[name]
        lanes = sorted(set(e.lane for e in entries if e.lane))
        if len(lanes) > 1:
            raise MjolnirError(
                "sample {0!r} is split across lanes {1}: {2}\n"
                "  concatenate the lanes first, keeping R1 and R2 separate, "
                "e.g.\n    cat {0}_*_L*_R1_001.fastq.gz > {0}_R1.fastq.gz\n"
                "    cat {0}_*_L*_R2_001.fastq.gz > {0}_R2.fastq.gz\n"
                "  Mjolnir will not merge them silently: how lanes are merged "
                "changes duplicate marking and read groups.".format(
                    name, ", ".join(lanes),
                    ", ".join(str(e.path.name) for e in entries))
            )

        by_read: Dict[int, List[ReadFileName]] = {}
        unnumbered: List[ReadFileName] = []
        for entry in entries:
            if entry.read is None:
                unnumbered.append(entry)
            else:
                by_read.setdefault(entry.read, []).append(entry)

        for entry in unnumbered:
            groups.append(ReadGroup(sample_id=entry.sample, files=[entry.path],
                                    paired=False, convention=CONVENTION_NONE))

        duplicated = [read for read, items in by_read.items() if len(items) > 1]
        if duplicated:
            raise MjolnirError(
                "sample {0!r} has more than one file for read {1}: {2}\n"
                "  give one file per read, or rename the extras.".format(
                    name, duplicated[0],
                    ", ".join(str(e.path.name) for e in entries))
            )

        if 1 in by_read and 2 in by_read:
            first = by_read[1][0]
            second = by_read[2][0]
            groups.append(ReadGroup(
                sample_id=name,
                files=[first.path, second.path],
                paired=True,
                convention=first.convention,
            ))
        else:
            for read, items in sorted(by_read.items()):
                entry = items[0]
                notes: List[str] = []
                if entry.ambiguous:
                    # The `_1` was not treated as a read number, so the sample
                    # keeps the whole stem; say so, because a user who really
                    # did lose an R2 needs to see why the name looks odd.
                    resolved_name = entry.stem
                    notes.append(
                        "{0} ends in '_{1}' but no partner file was given; "
                        "treated as single-end reads for sample {2!r} rather "
                        "than as read {1} of a pair".format(
                            entry.path.name, read, resolved_name))
                else:
                    resolved_name = entry.sample
                    notes.append(
                        "{0} is read {1} but its mate was not given; "
                        "treated as single-end".format(entry.path.name, read))
                groups.append(ReadGroup(
                    sample_id=resolved_name, files=[entry.path], paired=False,
                    convention=entry.convention, notes=notes))

    for entry in index_reads:
        LOG.warning("ignoring index read %s: barcode reads carry no sample "
                    "sequence", entry.path.name)
    if index_reads and groups:
        groups[0].notes.append(
            "ignored {0} index-read file(s) ({1}): index reads carry barcodes, "
            "not sample sequence".format(
                len(index_reads), ", ".join(e.path.name for e in index_reads)))
    return groups


# ---------------------------------------------------------------------------
# Platform detection
# ---------------------------------------------------------------------------

def filename_hints(path: PathLike) -> Dict[str, List[str]]:
    """Platform hints in the file name and its parent directory."""
    resolved = Path(path)
    haystack = "{0}/{1}".format(resolved.parent.name, resolved.name).lower()
    return {
        PLATFORM_ONT: [h for h in ONT_FILENAME_HINTS if h in haystack],
        PLATFORM_ILLUMINA: [h for h in ILLUMINA_FILENAME_HINTS if h in haystack],
    }


@dataclass
class PlatformEvidence:
    """Why the platform was called what it was called.

    Exists so that no caller has to take the label on trust. ``reasons`` is the
    measurement, ``hints`` is the filename, ``conflicts`` is where the two
    disagree, and ``checks`` is the same information in the form the report and
    the agent consume. A conflict never changes ``platform`` — the reads decide
    that — but it does lower ``confidence`` and it is always printed.
    """

    platform: str
    confidence: str = CONFIDENCE_NONE
    basis: str = ""
    reasons: List[str] = field(default_factory=list)
    hints: Dict[str, List[str]] = field(default_factory=dict)
    conflicts: List[str] = field(default_factory=list)
    profile: Optional[LengthProfile] = None
    stats: Optional[AssemblyStats] = None
    files: List[Path] = field(default_factory=list)
    forced: bool = False
    checks: List[Check] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable({
            "platform": self.platform,
            "confidence": self.confidence,
            "basis": self.basis,
            "reasons": list(self.reasons),
            "hints": dict(self.hints),
            "conflicts": list(self.conflicts),
            "profile": self.profile.to_dict() if self.profile else None,
            "stats": self.stats.to_dict() if self.stats else None,
            "files": [str(p) for p in self.files],
            "forced": self.forced,
            "checks": [c.to_dict() for c in self.checks],
        })


def _platform_check(evidence: PlatformEvidence) -> Check:
    """The platform call as a Check, so the report can print it like any other."""
    if evidence.forced:
        # An operator override is not a measurement and is not a failure; it is
        # a stated assumption, and it is flagged so the report carries it.
        status = STATUS_WARN
    elif evidence.confidence == CONFIDENCE_HIGH:
        status = STATUS_PASS
    elif evidence.confidence == CONFIDENCE_NONE:
        status = STATUS_FAIL
    else:
        status = STATUS_WARN
    return Check(
        name="platform_detection",
        value=evidence.platform,
        threshold=None,
        source=config.source_for("seqio_illumina_max_read_length"),
        status=status,
        reading=evidence.basis,
        comparison="==",
        category="platform",
        measured=not evidence.forced,
    )


def _length_checks(profile: LengthProfile, platform: str) -> List[Check]:
    """The read-length measurement as checks, including whether it agrees.

    The length check does not ask "are these reads good"; it asks whether the
    measured length distribution is consistent with the platform that was
    called. Those are different questions, and only the second one belongs
    here — depth, breadth and quality are ``engines/depth.py``'s to answer.
    """
    checks: List[Check] = []
    checks.append(Check.numeric(
        "reads_sampled_for_platform_call", float(profile.reads_sampled),
        warn_minimum=float(MIN_READS_FOR_PLATFORM_CALL),
        source=config.source_for("seqio_min_reads_for_platform_call"),
        unit="reads", category="platform",
        reading="reads examined from the head of the file to profile read length",
        not_measured_why="no reads could be read from the file"))
    if profile.max_length is None:
        checks.append(Check.not_measured(
            "read_length_distribution",
            "no reads were readable, so read length could not be profiled",
            source=config.source_for("seqio_illumina_max_read_length"),
            category="platform"))
    else:
        over_ceiling = profile.max_length > ILLUMINA_MAX_READ_LENGTH
        agrees = over_ceiling == (platform == PLATFORM_ONT)
        checks.append(Check(
            name="max_read_length",
            value=profile.max_length,
            threshold=ILLUMINA_MAX_READ_LENGTH,
            source=config.source_for("seqio_illumina_max_read_length"),
            status=STATUS_PASS if agrees else STATUS_WARN,
            reading=(
                "longest read in the sampled head, {0} the {1} nt Illumina "
                "ceiling - consistent with the {2} call".format(
                    "above" if over_ceiling else "within",
                    ILLUMINA_MAX_READ_LENGTH, platform)
                if agrees else
                "longest read in the sampled head is {0} the {1} nt Illumina "
                "ceiling, which does not follow from the {2} call".format(
                    "above" if over_ceiling else "within",
                    ILLUMINA_MAX_READ_LENGTH, platform)),
            comparison="<=",
            unit="nt",
            category="platform",
        ))
    if profile.quality_encoding == "phred64":
        checks.append(Check(
            name="fastq_quality_encoding",
            value="phred64",
            threshold="phred33",
            source=config.source_for("seqio_phred64_min_ascii"),
            status=STATUS_FAIL,
            reading="quality characters are in the legacy offset-64 range; "
                    "every base-quality threshold in this tool assumes "
                    "offset-33 and would be misread by 31 points. Re-convert "
                    "the file before running.",
            comparison="==",
            category="platform",
        ))
    elif profile.quality_encoding == "phred33":
        checks.append(Check(
            name="fastq_quality_encoding", value="phred33", threshold="phred33",
            source=config.source_for("seqio_phred33_max_ascii"), status=STATUS_PASS,
            reading="Sanger/Illumina-1.8+ offset-33 quality encoding",
            comparison="==", category="platform"))
    else:
        checks.append(Check.not_measured(
            "fastq_quality_encoding",
            "every sampled quality character sits in the range the offset-33 "
            "and offset-64 encodings share, so the encoding is undetermined",
            source=config.source_for("seqio_phred33_max_ascii"),
            category="platform"))
    return checks


def detect_platform(paths: Sequence[PathLike], paired: bool = False,
                    sample_reads: int = PLATFORM_SAMPLE_READS
                    ) -> PlatformEvidence:
    """Decide Illumina / ONT / FASTA for one sample's files, with the evidence.

    The order of argument is the order of strength: an assembled FASTA is
    settled by its first character; a pair of FASTQs that pair by filename and
    profile as short reads is Illumina beyond argument; a single FASTQ is
    decided by the length distribution of a bounded sample of its reads.
    """
    files = [Path(p) for p in paths]
    if not files:
        raise MjolnirError("detect_platform() was given no files")
    for path in files:
        if not path.exists():
            raise MjolnirError("input not found: {0}".format(path))
        if path.stat().st_size == 0:
            raise MjolnirError("input file is empty: {0}".format(path))

    fmt = sniff_format(files[0])
    if fmt == "":
        raise MjolnirError(
            "{0} is neither FASTA nor FASTQ: its first record line begins with "
            "neither '>' nor '@'.\n  If it is compressed, Mjolnir reads gzip, "
            "bzip2 and xz - check the file is not truncated.".format(files[0]))

    hints = filename_hints(files[0])

    if fmt == FORMAT_FASTA:
        stats = assembly_stats(files[0])
        if stats.contigs == 0:
            raise MjolnirError(
                "{0} is FASTA but contains no sequence records".format(files[0]))
        evidence = PlatformEvidence(
            platform=PLATFORM_FASTA,
            confidence=CONFIDENCE_HIGH,
            basis="assembled sequence: {0} contigs, {1:,} bp".format(
                stats.contigs, stats.total_length),
            reasons=["first record line begins with '>'",
                     "{0} contigs totalling {1:,} bp".format(
                         stats.contigs, stats.total_length)],
            hints=hints,
            stats=stats,
            files=files,
        )
        evidence.reasons.append(config.FASTA_CAPABILITY_LOSS)
        evidence.checks.append(_platform_check(evidence))
        evidence.checks.append(Check(
            name="input_capability",
            value="assembly",
            threshold="reads",
            source=config.source_for("fasta_capability_loss"),
            status=STATUS_WARN,
            reading=config.FASTA_CAPABILITY_LOSS,
            comparison="==",
            category="platform",
        ))
        return evidence

    profile = length_profile(files[0], max_records=sample_reads)
    if profile.reads_sampled == 0:
        raise MjolnirError(
            "{0} is FASTQ but no complete read could be parsed from it; the "
            "file is empty or truncated at its first record".format(files[0]))

    platform, confidence, reasons = _platform_from_profile(profile, paired)
    basis = reasons[0] if reasons else ""

    # A file too thin to profile is the one case where a filename carries any
    # weight, and then the evidence says exactly that. Handled before the
    # conflict check below, so that a hint which has already been acted on is
    # not then reported as a disagreement with itself.
    if profile.reads_sampled < MIN_READS_FOR_PLATFORM_CALL:
        hinted = [p for p in (PLATFORM_ONT, PLATFORM_ILLUMINA) if hints.get(p)]
        if len(hinted) == 1 and hinted[0] != platform:
            platform = hinted[0]
            basis = (
                "only {0} reads were readable - too few to profile - so the "
                "call rests on the filename ({1}), not on the data".format(
                    profile.reads_sampled, ", ".join(hints[platform])))
            reasons.append(basis)
        confidence = CONFIDENCE_LOW

    conflicts: List[str] = []
    other = PLATFORM_ONT if platform == PLATFORM_ILLUMINA else PLATFORM_ILLUMINA
    if hints.get(other):
        conflicts.append(
            "the filename suggests {0} ({1}) but the read lengths say {2}; the "
            "reads decide, and the disagreement is reported rather than "
            "resolved".format(other, ", ".join(hints[other]), platform))
        confidence = _lower_confidence(confidence)

    evidence = PlatformEvidence(
        platform=platform,
        confidence=confidence,
        basis=basis,
        reasons=reasons,
        hints=hints,
        conflicts=conflicts,
        profile=profile,
        files=files,
    )
    evidence.checks.append(_platform_check(evidence))
    evidence.checks.extend(_length_checks(profile, platform))
    for note in (compression_mismatch(f) for f in files):
        if note:
            evidence.conflicts.append(note)
    return evidence


def _lower_confidence(confidence: str) -> str:
    order = [CONFIDENCE_HIGH, CONFIDENCE_MODERATE, CONFIDENCE_LOW, CONFIDENCE_NONE]
    index = order.index(confidence) if confidence in order else len(order) - 1
    return order[min(index + 1, len(order) - 1)]


def _platform_from_profile(profile: LengthProfile, paired: bool
                           ) -> Tuple[str, str, List[str]]:
    """The read-length rule, separated so it can be tested without a file."""
    reasons: List[str] = []
    long_fraction = profile.long_fraction if profile.long_fraction is not None else 0.0
    max_length = profile.max_length or 0
    median = profile.median_length or 0.0

    if long_fraction >= ONT_LONG_READ_FRACTION:
        reasons.append(
            "{0:.1%} of {1} sampled reads are longer than {2} nt, which no "
            "Illumina chemistry produces".format(
                long_fraction, profile.reads_sampled, ILLUMINA_MAX_READ_LENGTH))
        reasons.append("median read length {0:,.0f} nt, longest {1:,} nt".format(
            median, max_length))
        if median >= ONT_INDICATIVE_MEDIAN_LENGTH or long_fraction >= 0.5:
            confidence = CONFIDENCE_HIGH
        else:
            confidence = CONFIDENCE_MODERATE
            reasons.append(
                "median below {0:,} nt, so this is long-read data but a short "
                "or degraded library".format(ONT_INDICATIVE_MEDIAN_LENGTH))
        if paired:
            reasons.append(
                "the two files pair by filename, which contradicts the read "
                "lengths")
        return PLATFORM_ONT, confidence, reasons

    if max_length <= ILLUMINA_MAX_READ_LENGTH:
        reasons.append(
            "all {0} sampled reads are {1} nt or shorter, within Illumina "
            "chemistry".format(profile.reads_sampled, max_length))
        if paired:
            reasons.append("the two files pair as R1/R2")
            confidence = CONFIDENCE_HIGH
        elif profile.reads_sampled >= MIN_READS_FOR_PLATFORM_CALL:
            confidence = CONFIDENCE_HIGH
            reasons.append(
                "single file, so these are unpaired or interleaved short reads")
        else:
            confidence = CONFIDENCE_LOW
        return PLATFORM_ILLUMINA, confidence, reasons

    reasons.append(
        "{0:.1%} of {1} sampled reads exceed {2} nt - below the {3:.0%} that "
        "would make this a long-read library, so these are read through into "
        "adapter or concatemers rather than ONT reads".format(
            long_fraction, profile.reads_sampled, ILLUMINA_MAX_READ_LENGTH,
            ONT_LONG_READ_FRACTION))
    reasons.append("longest sampled read {0:,} nt".format(max_length))
    return PLATFORM_ILLUMINA, CONFIDENCE_MODERATE, reasons


# ---------------------------------------------------------------------------
# The public entry point
# ---------------------------------------------------------------------------

@dataclass
class DetectedSample:
    """One sample as detected, with the evidence that produced it."""

    sample: SampleInput
    evidence: PlatformEvidence
    notes: List[str] = field(default_factory=list)

    @property
    def sample_id(self) -> str:
        return self.sample.sample_id

    @property
    def platform(self) -> str:
        return self.sample.platform

    @property
    def checks(self) -> List[Check]:
        return list(self.evidence.checks)

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable({
            "sample": self.sample.to_dict(),
            "evidence": self.evidence.to_dict(),
            "notes": list(self.notes),
        })


def expand_inputs(paths: Sequence[PathLike], recursive: bool = False) -> List[Path]:
    """Turn a mix of files and directories into a sorted list of sequence files.

    Directories are expanded to the FASTQ and FASTA files inside them, in
    natural order so that ``L2`` sorts after ``L10``'s sibling rather than
    between ``L1`` and ``L10``. A directory containing no sequence files is an
    error, not an empty cohort.
    """
    found: List[Path] = []
    for raw in paths:
        path = Path(raw).expanduser()
        if path.is_dir():
            walker = path.rglob("*") if recursive else path.glob("*")
            inside = [p for p in walker if p.is_file() and _looks_like_sequence_file(p)]
            if not inside:
                raise MjolnirError(
                    "no FASTQ or FASTA files found in {0}{1}".format(
                        path, "" if recursive else " (use --recursive to "
                                                    "descend into subdirectories)"))
            found.extend(inside)
        elif path.exists():
            found.append(path)
        else:
            raise MjolnirError("input not found: {0}".format(path))
    return sorted(set(found), key=lambda p: natural_key(str(p)))


def _looks_like_sequence_file(path: Path) -> bool:
    lowered = path.name.lower()
    for comp in COMPRESSION_EXTENSIONS:
        if lowered.endswith(comp):
            lowered = lowered[: -len(comp)]
            break
    return lowered.endswith(FASTQ_EXTENSIONS + FASTA_EXTENSIONS)


def detect_inputs(paths: Sequence[PathLike],
                  platform: Optional[str] = None,
                  sample_id: Optional[str] = None,
                  reference: Optional[PathLike] = None,
                  origin: str = "cli",
                  recursive: bool = False,
                  sample_reads: int = PLATFORM_SAMPLE_READS
                  ) -> List[DetectedSample]:
    """Turn what the user typed into validated :class:`SampleInput` records.

    *platform* forces the call rather than detecting it. Forcing is honoured —
    the operator may know something the reads do not show — but the measured
    evidence is still collected and, where it disagrees, the disagreement is
    recorded on the evidence and marked in the checks. Forcing a platform is not
    a way to make a contradiction disappear.

    *sample_id* may only be given for a single sample; naming a cohort with one
    name would silently collapse it.
    """
    files = expand_inputs(paths, recursive=recursive)
    if not files:
        raise MjolnirError("no input files given")

    forced = normalise_platform(platform) if platform else None

    fastas: List[Path] = []
    fastqs: List[Path] = []
    for path in files:
        fmt = sniff_format(path)
        if fmt == FORMAT_FASTA:
            fastas.append(path)
        elif fmt == FORMAT_FASTQ:
            fastqs.append(path)
        else:
            raise MjolnirError(
                "{0} is neither FASTA nor FASTQ; Mjolnir reads assemblies and "
                "reads, and cannot tell what this is".format(path))

    detected: List[DetectedSample] = []

    groups: List[ReadGroup] = group_reads(fastqs) if fastqs else []
    for group in groups:
        evidence = detect_platform(group.files, paired=group.paired,
                                   sample_reads=sample_reads)
        evidence = _apply_forced_platform(evidence, forced, group)
        detected.append(DetectedSample(
            sample=SampleInput(
                sample_id=group.sample_id,
                platform=evidence.platform,
                paths=list(group.files),
                reference=Path(reference) if reference else None,
                origin=origin,
                note="; ".join(group.notes),
            ),
            evidence=evidence,
            notes=list(group.notes),
        ))

    for path in fastas:
        evidence = detect_platform([path])
        evidence = _apply_forced_platform(evidence, forced, None)
        detected.append(DetectedSample(
            sample=SampleInput(
                # Not sample_name(): an assembly has no read number, so
                # "SA_R1.fasta" is a sample called SA_R1, not read 1 of a pair.
                sample_id=strip_extensions(path),
                platform=evidence.platform,
                paths=[path],
                reference=Path(reference) if reference else None,
                origin=origin,
            ),
            evidence=evidence,
        ))

    if sample_id:
        if len(detected) != 1:
            raise MjolnirError(
                "--sample names one sample, but {0} were detected in the given "
                "inputs: {1}".format(
                    len(detected), ", ".join(d.sample_id for d in detected)))
        detected[0].sample.sample_id = sample_id

    # Keyed on the name the OUTPUT FILES will carry, not on the raw id. Every
    # artefact is written as safe_name(sample_id), so "sample 1" and "sample/1"
    # are different ids and one set of files: the guard passed and the second
    # sample silently overwrote the first's report.
    seen: Dict[str, List[str]] = {}
    for entry in detected:
        seen.setdefault(safe_name(entry.sample_id), []).append(entry.sample_id)
    collisions = sorted((key, ids) for key, ids in seen.items() if len(ids) > 1)
    if collisions:
        detail = "; ".join(
            "{0} -> {1}".format(", ".join(sorted(set(ids))), key)
            for key, ids in collisions)
        raise MjolnirError(
            "these inputs resolve to the same output name: {0}\n"
            "  results are written as <sample>.json and friends, so two samples "
            "sharing an output name would overwrite each other. Rename the "
            "inputs or pass --sample.".format(detail))

    detected.sort(key=lambda d: natural_key(d.sample_id))
    return detected


def _apply_forced_platform(evidence: PlatformEvidence, forced: Optional[str],
                           group: Optional[ReadGroup]) -> PlatformEvidence:
    """Honour an operator's ``--platform``, keeping the measured disagreement."""
    if not forced or forced == evidence.platform:
        return evidence
    evidence.conflicts.append(
        "platform was forced to {0}; the measured evidence says {1} ({2})".format(
            forced, evidence.platform, evidence.basis))
    if group is not None and group.paired and forced != PLATFORM_ILLUMINA:
        raise MjolnirError(
            "sample {0!r} was given as two paired files but --platform {1} "
            "takes a single file. Pass one file, or drop --platform and let "
            "the reads decide.".format(group.sample_id, forced))
    evidence.platform = forced
    evidence.forced = True
    evidence.confidence = CONFIDENCE_NONE
    evidence.basis = "platform forced to {0} by the operator".format(forced)
    evidence.checks = [c for c in evidence.checks if c.name != "platform_detection"]
    evidence.checks.insert(0, _platform_check(evidence))
    return evidence


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_fastq(path: PathLike) -> None:
    """Raise unless *path* is a readable, non-empty FASTQ."""
    resolved = Path(path)
    if not resolved.exists():
        raise MjolnirError("input not found: {0}".format(path))
    if resolved.stat().st_size == 0:
        raise MjolnirError("input file is empty: {0}".format(path))
    fmt = sniff_format(path)
    if fmt == FORMAT_FASTA:
        raise MjolnirError(
            "{0} is named like reads but contains FASTA. Pass assemblies "
            "without --platform, or as --platform fasta.".format(path))
    if fmt != FORMAT_FASTQ:
        raise MjolnirError("{0} is not a FASTQ file".format(path))
    for _name, seq, _qual in read_fastq_head(path, max_records=1):
        if seq:
            return
    raise MjolnirError("no reads found in FASTQ: {0}".format(path))


def validate_fasta(path: PathLike) -> None:
    """Raise unless *path* is a readable FASTA holding at least one sequence."""
    resolved = Path(path)
    if not resolved.exists():
        raise MjolnirError("input not found: {0}".format(path))
    if resolved.stat().st_size == 0:
        raise MjolnirError("input file is empty: {0}".format(path))
    fmt = sniff_format(path)
    if fmt == FORMAT_FASTQ:
        raise MjolnirError(
            "{0} is named like an assembly but contains FASTQ reads.".format(path))
    if fmt != FORMAT_FASTA:
        raise MjolnirError("{0} is not a FASTA file".format(path))
    for _header, seq in read_fasta(path):
        if seq:
            return
    raise MjolnirError("no sequence records found in FASTA: {0}".format(path))


def describe_inputs(detected: Sequence[DetectedSample]) -> str:
    """A human-readable block for the log and for ``mjolnir doctor``.

    Prints the evidence, not only the verdict, because "sample 226-18: illumina"
    is exactly the kind of unexplained assertion this tool exists not to make.
    """
    lines: List[str] = []
    for entry in detected:
        evidence = entry.evidence
        lines.append("{0}  [{1}, {2} confidence]".format(
            entry.sample_id, evidence.platform, evidence.confidence))
        for path in entry.sample.paths:
            lines.append("    file:     {0}".format(path))
        for reason in evidence.reasons:
            lines.append("    evidence: {0}".format(reason))
        for note in entry.notes:
            lines.append("    note:     {0}".format(note))
        for conflict in evidence.conflicts:
            lines.append("    CONFLICT: {0}".format(conflict))
    return "\n".join(lines)
