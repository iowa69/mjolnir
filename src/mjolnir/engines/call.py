"""Variant calling, dispatched by platform, with the per-platform thresholds applied.

The design fixes the callers (§7): bcftools or freebayes on Illumina, Clair3 on
ONT, direct comparison for an assembly. What it also fixes, and what most of
this module is about, is what happens when the preferred caller is not there.

On ONT the fallback is a genuine degradation, not an equivalent. bcftools is
specifically weak on ONT indels — bcftools' own ``ont`` profile disables indel
calling outright — so a Mjolnir run that falls back to it has no indel calls at
all, which silently removes every loss-of-function resistance rule (§5.4). That
is a capability loss and the tool refuses to take it by accident: the fallback
has to be asked for, it stamps a caveat on the result, and the result carries
``degraded=True`` so the report can say which caller produced the numbers.

Allele fractions from a caller are not the whole story either, and are not meant
to be. Barcode genotyping and the allele fraction at every catalogue position
come from a direct pileup (``pileup.py``), because a caller that did not emit a
record at a position has said nothing about that position — and on ONT, where 26
of 27 Illumina-only minor variants were visible in the pileup but never called,
the difference between "not called" and "not there" is the whole question.

Variant normalisation — left-alignment, parsimony, HGVS naming — belongs to
``resistance/normalise.py``. This module reports what the caller emitted, at the
coordinates the caller emitted it, and marks rather than deletes anything that
fails a threshold, so the annex can show why a variant was not used.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

from ..config import (
    CALLERS,
    Config,
    DEGRADED_DEPTH_FLOOR,
    MAJOR_VARIANT_FRACTION,
    MIN_MINOR_VARIANT_FRACTION,
    MTBSEQ_MINBQUAL,
    ONT_MINOR_VARIANT_CAVEAT,
    min_reads_for,
    source_for,
)
from ..records import (
    PLATFORM_FASTA,
    PLATFORM_ILLUMINA,
    PLATFORM_ONT,
    VARIANT_DEL,
    VARIANT_INDEL,
    VARIANT_INS,
    VARIANT_MNV,
    VARIANT_SNP,
    Variant,
    normalise_platform,
)
from ..utils import (
    LOG,
    MjolnirError,
    PathLike,
    ensure_dir,
    first_available,
    have,
    require,
    require_file,
    safe_fraction,
    smart_open,
)
from .map import (
    MAX_PILEUP_DEPTH,
    MIN_BASE_QUALITY,
    MIN_MAPPING_QUALITY,
    SRC_MARIN_2022,
    run_pipeline,
)

# ---------------------------------------------------------------------------
# Caller parameters
# ---------------------------------------------------------------------------

#: SOURCE: bcftools mpileup ``--config ont`` profile, read from ``bcftools
#: mpileup -X list`` (bcftools 1.21 on this machine), which expands to exactly
#: ``-B -Q5 --max-BQ 30 -I``. Mjolnir spells the flags out rather than passing
#: ``-X ont`` so the invocation does not depend on the bcftools version, and so
#: that the ``-I`` — no indel calling — is visible in the command the report
#: prints rather than hidden inside a profile name.
#:
#: bcftools 1.20 and later also ship an ``ont-sup`` profile
#: (``--indels-cns -B -Q1 --max-BQ 35 -F0.2 -o15 -e1 -h110 …``) which does call
#: indels and is aimed at exactly the R10.4.1 + Dorado `sup` data the design
#: requires. Mjolnir does not select it, for two reasons worth stating: it needs
#: bcftools >= 1.20, and this path exists only as an announced degradation from
#: Clair3 — quietly turning indel calling back on inside the fallback would make
#: the fallback look like a substitute for the caller the design specifies.
ONT_BCFTOOLS_MIN_BASE_QUALITY = 5
ONT_BCFTOOLS_MAX_BASE_QUALITY = 30

#: SOURCE: design §7. The sentence the report must carry when ONT variants came
#: from the fallback rather than from Clair3.
ONT_FALLBACK_CAVEAT = (
    "Clair3 was unavailable, so ONT variants were called with bcftools using "
    "bcftools' own ONT profile, which disables indel calling because BCFtools "
    "indel calls on ONT are unreliable. This run therefore contains NO indel "
    "calls: loss-of-function resistance rules could not fire. That is an absence "
    "of capability, not an absence of indels."
)

#: SOURCE: Clair3 documentation. ``--haploid_precise`` treats only a homozygous
#: call as a variant, which is the correct model for a haploid bacterial genome
#: and is also the reason a Clair3 run reports no minor alleles at all. The
#: minor-variant question is answered from the pileup instead (design §6, §7).
CLAIR3_HAPLOID_MODE = "--haploid_precise"

#: SOURCE: Clair3 defaults for the ONT platform. Kept explicit because a caller
#: reading this file should not have to know what Clair3's built-in default is
#: to know what was run.
CLAIR3_SNP_MIN_AF = 0.08
CLAIR3_INDEL_MIN_AF = 0.15

#: SOURCE: design §7 — R10.4.1 + Dorado `sup` is the minimum credible ONT
#: configuration. The matching Clair3 model is the one trained for that
#: chemistry and basecaller. The exact directory name travels with the model
#: release, so this is the default Mjolnir looks for and not a value it invents:
#: when it is absent the error lists what the models directory actually holds.
CLAIR3_DEFAULT_MODEL = "r1041_e82_400bps_sup_v420"
CLAIR3_MODELS_SUBDIR = "clair3_models"

#: SOURCE: freebayes manual. ``--pooled-continuous`` reports every allele that
#: passes the input filters with its observation counts, instead of forcing a
#: genotype. That is what makes heteroresistance visible on Illumina, which the
#: design requires to be reported (§7).
FREEBAYES_POOLED_CONTINUOUS = True

#: SOURCE: freebayes manual, ``--haplotype-length``. Zero disables the joining of
#: nearby variants into haplotype alleles. WHO's catalogue is matched on
#: decomposed variants — its own MNVs are split and joined by `&` (§5.2) — so
#: emitting components is what the matcher expects. The cost is real and stated:
#: two SNPs in one codon arrive as two records, and the amino-acid consequence of
#: the pair is reconstructed by resistance/normalise.py, not by the caller.
FREEBAYES_HAPLOTYPE_LENGTH = 0

#: Filter labels Mjolnir attaches itself, spelled once so the annex and the
#: consensus engine cannot disagree about them.
FILTER_MIN_READS = "min-reads"
FILTER_MINOR_AF = "below-minor-af"
FILTER_LOW_DEPTH = "low-depth"

_ALT_PLACEHOLDERS = ("<*>", "<NON_REF>", "*", ".", "")


# ---------------------------------------------------------------------------
# The result record
# ---------------------------------------------------------------------------

@dataclass
class VariantCallResult:
    """Variants from one sample, with which caller produced them and how well."""

    sample_id: str
    platform: str
    caller: str
    vcf: Optional[Path] = None
    variants: List[Variant] = field(default_factory=list)
    commands: List[List[str]] = field(default_factory=list)
    #: True when the caller used was not the one the design prefers for this
    #: platform. The report prints this; it never has to infer it from the name.
    degraded: bool = False
    caveats: List[str] = field(default_factory=list)
    reference: str = ""

    @property
    def passing(self) -> List[Variant]:
        """Variants with no filter attached — the ones the rule layer may use."""
        return [variant for variant in self.variants if not variant.filters]

    @property
    def filtered_out(self) -> List[Variant]:
        return [variant for variant in self.variants if variant.filters]


# ---------------------------------------------------------------------------
# Pure command builders
# ---------------------------------------------------------------------------

def _output_type(path: PathLike) -> str:
    """bcftools ``-O`` letter implied by an output filename."""
    name = str(path)
    if name.endswith(".vcf.gz"):
        return "z"
    if name.endswith(".bcf"):
        return "b"
    if name.endswith(".vcf"):
        return "v"
    raise MjolnirError(
        "cannot tell the output format of {0!r}: use .vcf, .vcf.gz or .bcf".format(name))


def bcftools_mpileup_argv(bam: PathLike, reference: PathLike, *,
                          platform: str = PLATFORM_ILLUMINA,
                          threads: int = 1,
                          max_depth: int = MAX_PILEUP_DEPTH,
                          min_mapping_quality: int = MIN_MAPPING_QUALITY,
                          min_base_quality: int = MIN_BASE_QUALITY,
                          regions_bed: Optional[PathLike] = None,
                          skip_indels: Optional[bool] = None,
                          mtbseq_compat: bool = False) -> List[str]:
    """``bcftools mpileup``, emitting uncompressed BCF for the next stage.

    ``mtbseq_compat`` reproduces the legacy filter stack of design §9b: no MAPQ
    filter at all, base quality 13, anomalous read pairs re-included (``-A``),
    overlap detection off (``-x``), BAQ off (``-B``) and the 250x depth cap
    MTBseq inherits by never passing ``-d``. It is not a full reproduction —
    MTBseq's caller is a Perl majority rule over its own position table, whose
    denominator includes N and GAP — and the comparable path for that is
    ``pileup.py``'s ``mtbseq`` convention, not this one.
    """
    plat = normalise_platform(platform)
    argv = ["bcftools", "mpileup", "-f", str(reference),
            "-a", "FORMAT/AD,FORMAT/DP,FORMAT/SP",
            "--threads", str(max(1, threads))]

    if mtbseq_compat:
        argv += ["-B", "-A", "-x",
                 "-q", "0",
                 "-Q", str(MTBSEQ_MINBQUAL),
                 "-d", "250"]
    elif plat == PLATFORM_ONT:
        # bcftools' own ONT profile: no BAQ (it is an Illumina-shaped
        # correction), a low base-quality floor because ONT qualities are not
        # comparable to Illumina's, and a ceiling on the quality any base may
        # claim.
        argv += ["-B",
                 "-q", str(min_mapping_quality),
                 "-Q", str(ONT_BCFTOOLS_MIN_BASE_QUALITY),
                 "--max-BQ", str(ONT_BCFTOOLS_MAX_BASE_QUALITY),
                 "-d", str(max_depth)]
    else:
        argv += ["-q", str(min_mapping_quality),
                 "-Q", str(min_base_quality),
                 "-d", str(max_depth)]

    if skip_indels is None:
        skip_indels = plat == PLATFORM_ONT and not mtbseq_compat
    if skip_indels:
        argv.append("-I")
    if regions_bed is not None:
        argv += ["-R", str(regions_bed)]
    argv += ["-Ou", str(bam)]
    return argv


def bcftools_call_argv(*, threads: int = 1, ploidy: str = "1",
                       variants_only: bool = True,
                       output: Optional[PathLike] = None) -> List[str]:
    """``bcftools call`` in multiallelic mode, haploid.

    ``--ploidy 1`` is one of bcftools' predefined ploidies (confirmed against
    ``bcftools call --ploidy list``) and is not cosmetic: bcftools assumes all
    sites are diploid when none of ``--samples-file``, ``--ploidy`` or
    ``--ploidy-file`` is given, and under that model a fixed variant in a haploid
    genome is genotyped 1/1 while a real minority allele is genotyped 0/1 —
    inviting every downstream reader to treat an ONT error and a mixed infection
    as the same observation.
    """
    argv = ["bcftools", "call", "-m", "--ploidy", str(ploidy),
            "--threads", str(max(1, threads))]
    if variants_only:
        argv.append("-v")
    if output is None:
        argv.append("-Ou")
    else:
        argv += ["-O", _output_type(output), "-o", str(output)]
    return argv


def bcftools_norm_argv(reference: PathLike, output: PathLike, *,
                       threads: int = 1) -> List[str]:
    """``bcftools norm`` — left-align and split multiallelic records.

    Splitting is what makes ``Variant.coordinate_key`` usable at all: WHO's
    matching protocol is an exact match on one (CHROM, POS, REF, ALT), and a
    record carrying ``ALT=A,G`` matches neither.
    """
    return ["bcftools", "norm", "-f", str(reference), "-m", "-any",
            "--threads", str(max(1, threads)),
            "-O", _output_type(output), "-o", str(output)]


def bcftools_pipeline(bam: PathLike, reference: PathLike, output: PathLike, *,
                      platform: str = PLATFORM_ILLUMINA, threads: int = 1,
                      regions_bed: Optional[PathLike] = None,
                      mtbseq_compat: bool = False) -> List[List[str]]:
    """mpileup | call | norm, as argv lists."""
    return [
        bcftools_mpileup_argv(bam, reference, platform=platform, threads=threads,
                              regions_bed=regions_bed, mtbseq_compat=mtbseq_compat),
        bcftools_call_argv(threads=threads),
        bcftools_norm_argv(reference, output, threads=threads),
    ]


def freebayes_argv(bam: PathLike, reference: PathLike, output: PathLike, *,
                   min_alternate_count: int = 3,
                   min_alternate_fraction: float = MIN_MINOR_VARIANT_FRACTION,
                   min_coverage: int = DEGRADED_DEPTH_FLOOR,
                   min_mapping_quality: int = MIN_MAPPING_QUALITY,
                   min_base_quality: int = MIN_BASE_QUALITY,
                   pooled_continuous: bool = FREEBAYES_POOLED_CONTINUOUS,
                   haplotype_length: int = FREEBAYES_HAPLOTYPE_LENGTH,
                   regions_bed: Optional[PathLike] = None) -> List[str]:
    """``freebayes`` configured to report alleles rather than genotypes.

    *min_alternate_count* is the platform read-count threshold from config
    (``min_reads_for``), passed in rather than defaulted here so that the number
    the caller enforces and the number the report cites are the same number.
    """
    argv = ["freebayes", "-f", str(reference),
            "--min-alternate-count", str(min_alternate_count),
            "--min-alternate-fraction", str(min_alternate_fraction),
            "--min-coverage", str(min_coverage),
            "--min-mapping-quality", str(min_mapping_quality),
            "--min-base-quality", str(min_base_quality),
            "--haplotype-length", str(haplotype_length)]
    if pooled_continuous:
        argv.append("--pooled-continuous")
    else:
        argv += ["--ploidy", "1"]
    if regions_bed is not None:
        argv += ["--targets", str(regions_bed)]
    argv += ["--vcf", str(output), str(bam)]
    return argv


def clair3_argv(bam: PathLike, reference: PathLike, out_dir: PathLike, *,
                model_path: PathLike, threads: int = 1, sample_id: str = "",
                snp_min_af: float = CLAIR3_SNP_MIN_AF,
                indel_min_af: float = CLAIR3_INDEL_MIN_AF,
                haploid_mode: str = CLAIR3_HAPLOID_MODE,
                regions_bed: Optional[PathLike] = None) -> List[str]:
    """``run_clair3.sh`` for a haploid, non-human genome.

    Three flags are load-bearing on bacteria. ``--include_all_ctgs`` because
    Clair3 otherwise calls only human-shaped contig names and would return an
    empty VCF for NC_000962.3 without an error. ``--no_phasing_for_fa`` because
    whatshap phasing of a haploid genome is meaningless work. And
    ``--var_pct_full=1.0`` because the full-alignment model is the accurate one
    and running it over every candidate on a 4.4 Mb genome costs minutes, not
    the hours it would cost on a human sample.
    """
    if haploid_mode not in ("--haploid_precise", "--haploid_sensitive", ""):
        raise MjolnirError(
            "unknown Clair3 haploid mode {0!r}; expected --haploid_precise or "
            "--haploid_sensitive".format(haploid_mode))
    argv = ["run_clair3.sh",
            "--bam_fn={0}".format(bam),
            "--ref_fn={0}".format(reference),
            "--output={0}".format(out_dir),
            "--model_path={0}".format(model_path),
            "--threads={0}".format(max(1, threads)),
            "--platform=ont",
            "--include_all_ctgs",
            "--no_phasing_for_fa",
            "--var_pct_full=1.0",
            "--ref_pct_full=1.0",
            "--snp_min_af={0}".format(snp_min_af),
            "--indel_min_af={0}".format(indel_min_af)]
    if haploid_mode:
        argv.append(haploid_mode)
    if sample_id:
        argv.append("--sample_name={0}".format(sample_id))
    if regions_bed is not None:
        argv.append("--bed_fn={0}".format(regions_bed))
    return argv


def assembly_comparison_pipeline(assembly: PathLike, reference: PathLike, *,
                                 threads: int = 1,
                                 preset: str = "asm5") -> List[List[str]]:
    """The FASTA path: whole-genome alignment, then variants from the alignment.

    The VCF arrives on stdout — ``paftools.js call`` has no output flag — so the
    caller redirects it with ``run_pipeline(stdout_path=...)`` rather than being
    handed a filename here that nothing would use.

    ``minimap2 -c --cs`` then ``paftools.js call`` is minimap2's own documented
    route from an assembly to a VCF; the ``sort`` between them is required by
    paftools, which expects the PAF ordered by reference start. ``asm5`` is the
    right preset within a species — MTBC members sit at 99.21-99.92% ANI, well
    inside the 5% divergence the preset is tuned for.

    Everything this path produces has ``allele_fraction=None``. An assembly
    consensus is not evidence that an allele was fixed in the population, and
    writing 1.0 there would turn the design's stated capability loss (§7) into a
    confident-looking number.
    """
    return [
        ["minimap2", "-c", "--cs", "-x", str(preset), "-t", str(max(1, threads)),
         str(reference), str(assembly)],
        ["sort", "-k6,6", "-k8,8n"],
        ["paftools.js", "call", "-f", str(reference), "-"],
    ]


# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------

def classify_variant_type(ref: str, alt: str) -> str:
    """The ``records`` variant class implied by a REF/ALT pair."""
    ref, alt = ref.upper(), alt.upper()
    if len(ref) == len(alt):
        return VARIANT_SNP if len(ref) == 1 else VARIANT_MNV
    if len(alt) > len(ref):
        return VARIANT_INS if alt.startswith(ref) else VARIANT_INDEL
    return VARIANT_DEL if ref.startswith(alt) else VARIANT_INDEL


def _split_field(text: str) -> List[str]:
    return [part for part in str(text).split(",")]


def _as_int(text: Any) -> Optional[int]:
    try:
        return int(str(text))
    except (TypeError, ValueError):
        return None


def _as_float(text: Any) -> Optional[float]:
    try:
        return float(str(text))
    except (TypeError, ValueError):
        return None


def _parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for entry in str(info).split(";"):
        if not entry or entry == ".":
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
            out[key] = value
        else:
            out[entry] = "1"
    return out


def parse_vcf_lines(lines: Iterable[str], *, source_caller: str = "") -> List[Variant]:
    """Turn VCF text into :class:`~mjolnir.records.Variant` records.

    Written by hand rather than through pysam because the VCFs Mjolnir reads are
    small, because pysam is an optional extra, and mostly because the three
    callers spell read support three different ways — bcftools writes ``AD``,
    freebayes writes ``RO``/``AO``, Clair3 writes ``AD`` and ``AF`` — and the
    conversion between them is the part that has to be tested.

    Multiallelic records are split, one :class:`Variant` per ALT, with that ALT's
    own read counts. The allele fraction denominator is the sum of ``AD`` where
    it exists rather than ``DP``: ``DP`` at a site includes reads supporting
    neither allele, and dividing by it reports a fraction lower than the one the
    catalogue thresholds were derived against.
    """
    variants: List[Variant] = []
    for raw in lines:
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        chrom, pos_text, _id, ref, alt_text, qual_text, filter_text, info_text = fields[:8]

        info = _parse_info(info_text)
        sample_values: Dict[str, str] = {}
        if len(fields) >= 10:
            keys = fields[8].split(":")
            values = fields[9].split(":")
            sample_values = dict(zip(keys, values))

        vcf_filters = [entry for entry in str(filter_text).split(";")
                       if entry not in ("PASS", ".", "")]

        alts = _split_field(alt_text)
        # AD is ref-first; AO (freebayes) is alt-only.
        ad = [_as_int(value) for value in _split_field(sample_values.get("AD", ""))]
        ao = [_as_int(value) for value in _split_field(sample_values.get("AO",
                                                                        info.get("AO", "")))]
        ro = _as_int(sample_values.get("RO", info.get("RO", "")))
        af_field = [_as_float(value) for value in _split_field(
            sample_values.get("AF", info.get("AF", "")))]
        depth = _as_int(sample_values.get("DP", info.get("DP", "")))

        for index, alt in enumerate(alts):
            if alt in _ALT_PLACEHOLDERS:
                continue

            ref_reads: Optional[int] = None
            alt_reads: Optional[int] = None
            if len(ad) == len(alts) + 1:
                ref_reads, alt_reads = ad[0], ad[index + 1]
            elif len(ad) == len(alts) and len(alts) == 1:
                # Some writers emit AD with only the alt when the record is
                # already biallelic and split.
                alt_reads = ad[0]
            if alt_reads is None and index < len(ao):
                alt_reads, ref_reads = ao[index], ro

            counted = [value for value in ([ref_reads] + [alt_reads]) if value is not None]
            denominator: Optional[int] = None
            if ref_reads is not None and alt_reads is not None:
                denominator = sum(int(value) for value in ad[:len(alts) + 1]
                                  if value is not None) or sum(counted)
            elif depth is not None:
                denominator = depth

            fraction: Optional[float] = None
            if alt_reads is not None and denominator:
                fraction = safe_fraction(alt_reads, denominator)
            elif index < len(af_field) and af_field[index] is not None:
                fraction = af_field[index]

            variants.append(Variant(
                chrom=chrom,
                pos=int(pos_text),
                ref=ref.upper(),
                alt=alt.upper(),
                depth=depth,
                allele_fraction=fraction,
                is_major=None,          # set by apply_platform_filters, with the threshold
                source_caller=source_caller,
                alt_reads=alt_reads,
                ref_reads=ref_reads,
                qual=_as_float(qual_text),
                variant_type=classify_variant_type(ref, alt),
                filters=list(vcf_filters),
            ))
    return variants


def parse_vcf(path: PathLike, *, source_caller: str = "") -> List[Variant]:
    """Read a VCF from disk, transparently handling ``.vcf.gz``."""
    with smart_open(path, "rt") as handle:
        return parse_vcf_lines(handle, source_caller=source_caller)


# ---------------------------------------------------------------------------
# The platform thresholds
# ---------------------------------------------------------------------------

def apply_platform_filters(variants: Sequence[Variant], platform: str, *,
                           config: Optional[Config] = None) -> List[Variant]:
    """Apply §7's read-count and allele-fraction thresholds, marking not deleting.

    Variants that fail keep their place in the list with a filter label, because
    the annex has to be able to show that a catalogued mutation was seen at two
    reads and rejected — which is a different statement from never having been
    seen, and the one a reviewer will ask about.

    On ONT no minor-allele floor is applied. Minor variants there are not
    quantifiable: the caller under-detects them (26 of 27 in the 508-isolate
    comparison), so a fraction below the major threshold is annotated with the
    platform caveat rather than being filtered as if the number could be trusted
    either way.
    """
    plat = normalise_platform(platform)
    if plat == PLATFORM_FASTA:
        # No reads, so neither threshold has a denominator. Marking these as
        # failing a read-count test would be inventing evidence about them.
        for variant in variants:
            if not variant.note:
                variant.note = ("assembly input: no read support and no allele "
                                "fraction exist for this variant")
        return list(variants)

    min_reads = config.min_reads(plat) if config is not None else min_reads_for(plat)
    major_fraction = (config.major_variant_fraction if config is not None
                      else MAJOR_VARIANT_FRACTION)
    minor_floor = (config.min_minor_variant_fraction if config is not None
                   else MIN_MINOR_VARIANT_FRACTION)
    depth_floor = (config.degraded_depth_floor if config is not None
                   else DEGRADED_DEPTH_FLOOR)

    for variant in variants:
        if variant.allele_fraction is not None:
            variant.is_major = variant.allele_fraction >= major_fraction
        if variant.alt_reads is not None and variant.alt_reads < min_reads:
            variant.filters.append("{0}<{1}".format(FILTER_MIN_READS, min_reads))
        if variant.depth is not None and variant.depth < depth_floor:
            variant.filters.append("{0}<{1}x".format(FILTER_LOW_DEPTH, depth_floor))
        if plat == PLATFORM_ILLUMINA:
            if variant.allele_fraction is not None and variant.allele_fraction < minor_floor:
                variant.filters.append("{0}<{1}".format(FILTER_MINOR_AF, minor_floor))
        elif variant.is_major is False and not variant.note:
            variant.note = ONT_MINOR_VARIANT_CAVEAT
    return list(variants)


def threshold_sources(platform: str) -> Dict[str, str]:
    """Which published number each applied threshold came from, for the annex."""
    plat = normalise_platform(platform)
    key = "min_reads_illumina" if plat == PLATFORM_ILLUMINA else "min_reads_ont"
    return {
        "min_reads": source_for(key),
        "major_variant_fraction": source_for("major_variant_fraction"),
        "min_minor_variant_fraction": source_for("min_minor_variant_fraction"),
        "min_mapping_quality": SRC_MARIN_2022,
    }


# ---------------------------------------------------------------------------
# Clair3 model resolution
# ---------------------------------------------------------------------------

def clair3_model_path(db_dir: Optional[PathLike] = None,
                      model: Optional[str] = None) -> Path:
    """Locate the Clair3 model directory, or say exactly what to fetch.

    A Clair3 model is chemistry- and basecaller-specific. Running an R9 model
    over R10.4.1 ``sup`` reads produces calls that look ordinary and are wrong,
    so there is no default-to-whatever-is-there behaviour here: either the named
    model exists or the run stops, listing what the models directory does hold.
    """
    name = model or CLAIR3_DEFAULT_MODEL
    candidate = Path(str(name)).expanduser()
    if candidate.is_dir():
        return candidate

    if db_dir is None:
        raise MjolnirError(
            "Clair3 model {0!r} not found and no database directory was given.\n"
            "  set MJOLNIR_DB, or pass --clair3-model /path/to/{0}".format(name))

    root = Path(str(db_dir)).expanduser() / CLAIR3_MODELS_SUBDIR
    target = root / name
    if target.is_dir():
        return target

    available = sorted(child.name for child in root.iterdir()) if root.is_dir() else []
    raise MjolnirError(
        "Clair3 model {0!r} not found under {1}.\n"
        "  models present: {2}\n"
        "  download one from https://github.com/HKU-BAL/Clair3 (pre-trained "
        "models) and unpack it into that directory, or pass "
        "--clair3-model /path/to/{0}\n"
        "  (the model must match the chemistry and basecaller: the design "
        "requires R10.4.1 with Dorado sup)".format(
            name, root, ", ".join(available) if available else "none")
    )


# ---------------------------------------------------------------------------
# The wrapper that runs it
# ---------------------------------------------------------------------------

def choose_caller(platform: str, *, allow_degraded_fallback: bool = False) -> str:
    """The caller to use on this platform, given what is installed.

    On Illumina bcftools and freebayes are both first-class and either is
    accepted. On ONT there is one preferred caller and one degraded one, and the
    degraded one has to be asked for.
    """
    plat = normalise_platform(platform)
    if plat == PLATFORM_FASTA:
        return "direct-comparison"
    if plat == PLATFORM_ILLUMINA:
        found = first_available("bcftools", "freebayes")
        if found is None:
            raise MjolnirError(
                "no Illumina variant caller found on PATH (looked for bcftools, "
                "freebayes).\n"
                "  conda install -c conda-forge -c bioconda bcftools")
        return found

    if have("run_clair3.sh"):
        return "clair3"
    if allow_degraded_fallback and have("bcftools"):
        LOG.warning("Clair3 not found; falling back to bcftools for ONT. %s",
                    ONT_FALLBACK_CAVEAT)
        return "bcftools"
    raise MjolnirError(
        "Clair3 is the ONT variant caller (design §7) and run_clair3.sh is not on "
        "PATH.\n"
        "  conda install -c conda-forge -c bioconda clair3\n"
        "  then install a Clair3 model matching your chemistry (R10.4.1 + "
        "Dorado sup, e.g. {0}) and point at it with --clair3-model\n"
        "The bcftools fallback disables indel calling entirely and must be asked "
        "for explicitly with --allow-degraded-ont-calling.".format(CLAIR3_DEFAULT_MODEL)
    )


def call_variants(sample_id: str, bam: PathLike, reference: PathLike, platform: str,
                  out_dir: PathLike, *, config: Optional[Config] = None,
                  caller: Optional[str] = None, threads: int = 1,
                  clair3_model: Optional[str] = None,
                  allow_degraded_fallback: bool = False,
                  regions_bed: Optional[PathLike] = None,
                  assembly: Optional[PathLike] = None) -> VariantCallResult:
    """Call variants for one sample and return them with their provenance.

    *assembly* is required for FASTA input and *bam* is ignored there: an
    assembly is compared to the reference directly and never goes through an
    aligner-plus-caller path (§7).
    """
    plat = normalise_platform(platform)
    if config is not None and threads == 1:
        threads = config.threads
    out_dir = ensure_dir(out_dir)
    reference = require_file(reference, "reference FASTA")
    if not Path(str(reference) + ".fai").exists():
        raise MjolnirError(
            "reference {0} has no .fai index, which every caller here needs.\n"
            "  samtools faidx {0}".format(reference))

    chosen = caller or choose_caller(plat, allow_degraded_fallback=allow_degraded_fallback)
    if chosen not in CALLERS[plat] and chosen != "direct-comparison":
        raise MjolnirError(
            "caller {0!r} is not a {1} caller; expected one of {2}".format(
                chosen, plat, ", ".join(CALLERS[plat])))

    mtbseq_compat = bool(config.mtbseq_compat) if config is not None else False
    commands: List[List[str]] = []
    caveats: List[str] = []
    degraded = False

    if plat == PLATFORM_FASTA:
        if assembly is None:
            raise MjolnirError(
                "sample {0!r} is FASTA input but no assembly path was given to "
                "call_variants".format(sample_id))
        require("minimap2", "assembly-to-reference comparison")
        if not have("paftools.js"):
            # paftools.js ships inside the minimap2 distribution and needs the k8
            # javascript shell to run, so neither name alone is a useful hint.
            raise MjolnirError(
                "paftools.js not found on PATH; it is what turns a whole-genome "
                "alignment into variant calls for FASTA input.\n"
                "  conda install -c conda-forge -c bioconda minimap2 k8")
        vcf = out_dir / "{0}.assembly.vcf".format(sample_id)
        stages = assembly_comparison_pipeline(assembly, reference, threads=threads)
        run_pipeline(stages, stdout_path=vcf)
        commands.extend(stages)
        variants = parse_vcf(vcf, source_caller="minimap2/paftools")

    elif chosen == "clair3":
        require("run_clair3.sh", "ONT variant calling")
        model = clair3_model_path(
            config.db_dir if config is not None else None, clair3_model)
        clair_dir = ensure_dir(out_dir / "clair3")
        argv = clair3_argv(bam, reference, clair_dir, model_path=model,
                           threads=threads, sample_id=sample_id,
                           regions_bed=regions_bed)
        run_pipeline([argv])
        commands.append(argv)
        vcf = clair_dir / "merge_output.vcf.gz"
        if not vcf.exists():
            raise MjolnirError(
                "Clair3 finished but produced no {0}. Check {1} for its logs; the "
                "usual cause on a bacterial genome is a missing --include_all_ctgs, "
                "which Mjolnir does pass.".format(vcf.name, clair_dir))
        variants = parse_vcf(vcf, source_caller="clair3")

    elif chosen == "freebayes":
        require("freebayes", "Illumina variant calling")
        vcf = out_dir / "{0}.freebayes.vcf".format(sample_id)
        min_alt = config.min_reads(plat) if config is not None else min_reads_for(plat)
        minor = (config.min_minor_variant_fraction if config is not None
                 else MIN_MINOR_VARIANT_FRACTION)
        floor = (config.degraded_depth_floor if config is not None
                 else DEGRADED_DEPTH_FLOOR)
        argv = freebayes_argv(bam, reference, vcf, min_alternate_count=min_alt,
                              min_alternate_fraction=minor, min_coverage=floor,
                              regions_bed=regions_bed)
        run_pipeline([argv])
        commands.append(argv)
        variants = parse_vcf(vcf, source_caller="freebayes")

    elif chosen == "bcftools":
        require("bcftools", "variant calling")
        vcf = out_dir / "{0}.bcftools.vcf".format(sample_id)
        stages = bcftools_pipeline(bam, reference, vcf, platform=plat, threads=threads,
                                   regions_bed=regions_bed, mtbseq_compat=mtbseq_compat)
        run_pipeline(stages)
        commands.extend(stages)
        variants = parse_vcf(vcf, source_caller="bcftools")
        if plat == PLATFORM_ONT:
            degraded = True
            caveats.append(ONT_FALLBACK_CAVEAT)

    else:
        raise MjolnirError("no implementation for caller {0!r}".format(chosen))

    apply_platform_filters(variants, plat, config=config)
    if plat == PLATFORM_ONT and not degraded:
        caveats.append(ONT_MINOR_VARIANT_CAVEAT)

    LOG.info("%s: %d variants from %s (%d passing the %s thresholds)",
             sample_id, len(variants), chosen,
             len([v for v in variants if not v.filters]), plat)

    return VariantCallResult(
        sample_id=sample_id,
        platform=plat,
        caller=chosen,
        vcf=vcf,
        variants=variants,
        commands=commands,
        degraded=degraded,
        caveats=caveats,
        reference=str(reference),
    )
