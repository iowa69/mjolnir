"""Every threshold Mjolnir applies, with the source that published it.

The first house rule is that no number appears in logic. A bare ``25`` inside
``depth.py`` is a number nobody can audit; ``MIN_DEPTH`` defined here, with the
paper that measured it, is a number the report can print and a reviewer can
disagree with.

So every constant below is registered as well as assigned. ``_define`` returns
the value — so the constant is an ordinary module-level number and costs nothing
at the call site — and simultaneously files it in :data:`THRESHOLDS` with its
source string, its unit and whether the citation was verified against the
primary document. :func:`source_for` is what the report calls when it prints a
threshold, and :func:`unverified` is what ``doctor`` calls to list the numbers
that still need someone to check them.

Two entries here are refusals rather than thresholds. Kraken2's ``--confidence``
must never default to 0.0 for a contamination screen, and a standard or capped
Kraken2 index is not a mycobacterial contamination screen at all. Both are
enforced by functions in this module rather than left as advice in a docstring,
because advice in a docstring is how a tool ends up printing a tail of
low-abundance NTM as "co-infections".
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .records import (
    CALL_NO_CALL,
    CALL_R,
    CALL_R_INTERIM,
    CALL_S,
    CALL_S_INTERIM,
    CALL_UNCERTAIN,
    PLATFORM_FASTA,
    PLATFORM_ILLUMINA,
    PLATFORM_ONT,
    normalise_platform,
)
from .utils import LOG, MjolnirError, cpu_count

PACKAGE_DIR = Path(__file__).resolve().parent
BUNDLED_DATA = PACKAGE_DIR / "data"

# ---------------------------------------------------------------------------
# Sources. One string per document, so a citation is spelled once and a typo in
# it is visible in every threshold that shares it.
# ---------------------------------------------------------------------------

SRC_COLPUS_2026 = (
    "Colpus et al. 2026, bioRxiv - 508-isolate ONT vs Illumina MTBC comparison "
    "(clinical DST thresholds, masking, platform discordance)"
)
SRC_HALL_2024 = (
    "Hall et al. 2024, eLife - variant calling accuracy across depth in "
    "M. tuberculosis; precision and recall degrade notably below 25x"
)
SRC_MTBSEQ_MANUAL = "MTBseq MANUAL.md (ngs-fzb/MTBseq_source), default parameters"
SRC_MIXINFECT = (
    "Sobkowiak et al. 2018, BMC Genomics 19:613 - MixInfect: detection of mixed "
    "M. tuberculosis infections from WGS"
)
SRC_WHO_V2 = (
    "WHO catalogue of mutations in Mycobacterium tuberculosis complex and their "
    "association with drug resistance, 2nd edition (WHO-UCN-TB-2023.7)"
)
SRC_WALKER_2013 = (
    "Walker et al. 2013, Lancet Infect Dis 13:137 - WGS to delineate "
    "M. tuberculosis outbreaks; 5-SNP and 12-SNP thresholds"
)
SRC_NC_000962 = "GenBank NC_000962.3 - M. tuberculosis H37Rv complete genome"
SRC_TBDB = "jodyphelan/tbdb - TB-Profiler library (mutations.csv, barcode.bed, mask.bed)"
SRC_SNPIT = (
    "Lipworth et al. 2019, Emerg Infect Dis 25:482 - SNP-IT; M. bovis is defined "
    "by 23 phylogenetic SNPs"
)
SRC_ANI_SPECIES = (
    "Richter & Rossello-Mora 2009, PNAS 106:19126 and Jain et al. 2018, Nat "
    "Commun 9:5114 - 95-96% ANI as the prokaryotic species boundary"
)
SRC_NASH_2009 = (
    "Nash, Brown-Elliott & Wallace 2009, Antimicrob Agents Chemother 53:1367 - "
    "erm(41) confers inducible macrolide resistance in M. abscessus"
)
SRC_BASTIAN_2011 = (
    "Bastian et al. 2011, Antimicrob Agents Chemother 55:775 - erm(41) T28 vs "
    "C28 sequevar and rrl sequencing predict clarithromycin susceptibility"
)
SRC_WALLACE_1996 = (
    "Wallace et al. 1996, Antimicrob Agents Chemother 40:1676 - rrl 2058/2059 "
    "mutations and clarithromycin resistance in M. chelonae/M. abscessus"
)
SRC_PRAMMANANAN_1998 = (
    "Prammananan et al. 1998, J Infect Dis 177:1573 - rrs 1408 substitution "
    "confers amikacin resistance in M. abscessus and M. chelonae"
)
SRC_DESIGN = "Mjolnir design, docs/superpowers/specs/2026-08-10-mjolnir-design.md"
SRC_POLICY = (
    "Mjolnir policy - not a published threshold; chosen here, printed as a "
    "policy choice wherever it is applied"
)

# ---------------------------------------------------------------------------
# The threshold registry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Threshold:
    """One number and the document that justifies it."""

    name: str
    value: Any
    source: str
    unit: str = ""
    note: str = ""
    #: False when the citation was written from memory and has not been checked
    #: against the primary document on this machine. ``mjolnir doctor`` lists
    #: these, and the report marks them, because an unverified citation is
    #: worse than none: it looks settled.
    verified: bool = True

    def describe(self) -> str:
        unit = " {0}".format(self.unit) if self.unit else ""
        mark = "" if self.verified else " [citation unverified]"
        return "{0} = {1}{2} ({3}){4}".format(
            self.name, self.value, unit, self.source, mark)


THRESHOLDS: Dict[str, Threshold] = {}


def _define(name: str, value: Any, source: str, unit: str = "", note: str = "",
            verified: bool = True) -> Any:
    """Register a threshold and hand back its value.

    Registration is not optional bookkeeping: ``source_for(name)`` is how the
    report prints the provenance of a number, and a constant that skipped this
    function would appear in the PDF with no attribution at all.
    """
    if name in THRESHOLDS:
        raise MjolnirError("threshold {0!r} defined twice in config.py".format(name))
    THRESHOLDS[name] = Threshold(name=name, value=value, source=source, unit=unit,
                                 note=note, verified=verified)
    return value


def source_for(name: str) -> str:
    """The source string for a registered threshold, or a loud placeholder."""
    entry = THRESHOLDS.get(name)
    if entry is None:
        return "unregistered threshold {0!r} - no source recorded".format(name)
    return entry.source


def threshold(name: str) -> Threshold:
    entry = THRESHOLDS.get(name)
    if entry is None:
        raise MjolnirError(
            "no threshold named {0!r}; every number in logic must be defined in "
            "config.py with its source".format(name)
        )
    return entry


def unverified() -> List[Threshold]:
    """Thresholds whose citation has not been checked against the primary source."""
    return [t for t in THRESHOLDS.values() if not t.verified]


def all_thresholds() -> List[Threshold]:
    return sorted(THRESHOLDS.values(), key=lambda t: t.name)


# ---------------------------------------------------------------------------
# Reference
# ---------------------------------------------------------------------------

#: SOURCE: GenBank NC_000962.3. Every MTBC coordinate in this tool, in the WHO
#: catalogue's Genomic_coordinates sheet and in tbdb's barcode.bed and mask.bed
#: is against this accession. A coordinate against any other assembly is wrong
#: by up to several kilobases and produces silently mis-graded variants.
H37RV_ACCESSION = _define(
    "h37rv_accession", "NC_000962.3", SRC_NC_000962,
    note="H37Rv, the coordinate system for all MTBC catalogues Mjolnir reads")

#: SOURCE: NC_000962.3 record length.
H37RV_LENGTH = _define(
    "h37rv_length", 4411532, SRC_NC_000962, unit="bp")

#: SOURCE: NC_000962.3 composition. Used only as an expectation band for the QC
#: GC check; a sample far outside it is a composition warning, not a species call.
H37RV_GC = _define(
    "h37rv_gc", 0.656, SRC_NC_000962, unit="fraction",
    note="H37Rv GC content; NTM references carry their own value in the registry")

#: SOURCE: Mjolnir policy. GC band around the reference's own GC beyond which the
#: composition check warns. Wide enough not to fire on ordinary library bias.
GC_TOLERANCE = _define(
    "gc_tolerance", 0.03, SRC_POLICY, unit="fraction",
    note="warn when observed GC differs from the reference by more than this")

# ---------------------------------------------------------------------------
# Platform thresholds (design §7)
# ---------------------------------------------------------------------------

#: SOURCE: Colpus et al. 2026 bioRxiv, 508-isolate ONT-vs-Illumina MTBC study.
#: The published clinical-DST thresholds: a variant needs at least 3 supporting
#: reads on Illumina and at least 5 on ONT.
MIN_READS_ILLUMINA = _define(
    "min_reads_illumina", 3, SRC_COLPUS_2026, unit="reads",
    note="minimum reads supporting a variant on Illumina")

#: SOURCE: Colpus et al. 2026 bioRxiv (as above). Higher on ONT because the
#: per-read error rate is higher, not because more depth is available.
MIN_READS_ONT = _define(
    "min_reads_ont", 5, SRC_COLPUS_2026, unit="reads",
    note="minimum reads supporting a variant on ONT")

#: SOURCE: Colpus et al. 2026 bioRxiv. At or above 90% read support the variant
#: is major; below it, it is a minor variant, which is a different clinical
#: statement and is reported as one.
MAJOR_VARIANT_FRACTION = _define(
    "major_variant_fraction", 0.90, SRC_COLPUS_2026, unit="fraction",
    note="read support at or above which a variant is called major")

#: SOURCE: Mjolnir policy. Floor for reporting a minor variant at all on
#: Illumina; below this the signal is indistinguishable from sequencing error at
#: ordinary depths. Not applied to ONT, where minor variants are reported as
#: under-detected rather than quantified (see ONT_MINOR_VARIANT_CAVEAT).
MIN_MINOR_VARIANT_FRACTION = _define(
    "min_minor_variant_fraction", 0.05, SRC_POLICY, unit="fraction",
    note="lowest minor-allele fraction Mjolnir will report on Illumina")

#: SOURCE: Hall et al. 2024, eLife. 25x is the target depth; precision and
#: recall degrade notably below it.
MIN_DEPTH = _define(
    "min_depth", 25, SRC_HALL_2024, unit="x",
    note="target mean depth; below this, variant calling accuracy is known to fall")

#: SOURCE: Hall et al. 2024, eLife, read with the design's §7 table. 10x is the
#: degraded floor: between 10x and 25x the result is reported with a degraded-
#: depth caveat on every call; below 10x it is not a result.
DEGRADED_DEPTH_FLOOR = _define(
    "degraded_depth_floor", 10, SRC_HALL_2024, unit="x",
    note="below this depth the sample is not callable and is reported as such")

#: SOURCE: Mjolnir policy. Fraction of the reference that must reach
#: DEGRADED_DEPTH_FLOOR before genome-wide statements (lineage, cohort distance)
#: are made at all.
MIN_BREADTH = _define(
    "min_breadth", 0.95, SRC_POLICY, unit="fraction",
    note="fraction of the reference covered at the degraded floor")

#: SOURCE: Mjolnir policy. Mapped-read fraction below which the library is
#: flagged as substantially non-target. It is a flag, not a species claim.
MIN_MAPPED_FRACTION = _define(
    "min_mapped_fraction", 0.90, SRC_POLICY, unit="fraction",
    note="reads mapping to the chosen reference; below this, purity is suspect")

#: SOURCE: Mjolnir policy, definition stated wherever the number is printed.
#: Evenness is the fraction of reference positions whose depth lies between 0.5x
#: and 2x the mean. There is no standard definition, so the definition travels
#: with the number.
COVERAGE_EVENNESS_BAND = _define(
    "coverage_evenness_band", (0.5, 2.0), SRC_POLICY, unit="x mean depth",
    note="depth band, as multiples of mean depth, used to define evenness")
MIN_COVERAGE_EVENNESS = _define(
    "min_coverage_evenness", 0.80, SRC_POLICY, unit="fraction",
    note="fraction of positions inside the evenness band")
EVENNESS_DEFINITION = (
    "fraction of reference positions with depth between {0}x and {1}x the mean"
    .format(*COVERAGE_EVENNESS_BAND)
)

#: Mapper and caller per platform (design §7). Fallbacks are genuine
#: equivalents; there is no fallback that silently changes what is measured.
MAPPERS: Dict[str, Tuple[str, ...]] = {
    PLATFORM_ILLUMINA: ("bwa-mem2", "bwa"),
    PLATFORM_ONT: ("minimap2",),
    PLATFORM_FASTA: (),
}
MINIMAP2_ONT_PRESET = _define(
    "minimap2_ont_preset", "map-ont", SRC_DESIGN,
    note="minimap2 preset for ONT reads")
CALLERS: Dict[str, Tuple[str, ...]] = {
    PLATFORM_ILLUMINA: ("bcftools", "freebayes"),
    # Clair3 preferred on ONT: Clair3/DeepVariant lead on ONT bacterial data and
    # BCFtools is specifically weak on ONT indels (design §7).
    PLATFORM_ONT: ("clair3", "bcftools"),
    PLATFORM_FASTA: ("direct-comparison",),
}

#: SOURCE: Colpus et al. 2026 bioRxiv. 26 of 27 Illumina-only minor SNPs were
#: visible in the ONT pileup but not called, so absence of a minor variant on
#: ONT is not evidence of absence of a subpopulation.
ONT_MINOR_VARIANT_CAVEAT = _define(
    "ont_minor_variant_caveat",
    "ONT under-detects minor resistance variants: in the 508-isolate comparison "
    "26 of 27 Illumina-only minor SNPs were visible in the ONT pileup but not "
    "called. Absence of a minor variant here is not absence of a subpopulation.",
    SRC_COLPUS_2026)

#: SOURCE: Colpus et al. 2026 bioRxiv. fbiC tandem-repeat deletions produced
#: 47.2% of all discordant drug classifications, as spurious delamanid
#: resistance. Mjolnir suppresses this specific gene-drug call on ONT and says
#: it did.
ONT_SUPPRESSED_GENE_DRUG: Tuple[Tuple[str, str], ...] = _define(
    "ont_suppressed_gene_drug", (("fbiC", "Delamanid"),), SRC_COLPUS_2026,
    note="gene-drug calls suppressed on ONT because they are platform artefacts")
ONT_FBIC_DISCORDANCE_FRACTION = _define(
    "ont_fbic_discordance_fraction", 0.472, SRC_COLPUS_2026, unit="fraction",
    note="share of all discordant drug classifications attributable to fbiC on ONT")
ONT_FBIC_CAVEAT = (
    "delamanid resistance from an fbiC tandem-repeat deletion is suppressed on "
    "ONT: such calls accounted for 47.2% of all discordant drug classifications "
    "in the 508-isolate comparison"
)

#: SOURCE: Colpus et al. 2026 bioRxiv. About 16.6% of ONT indel calls were not
#: corroborated by Illumina, so any indel-driven call — notably the loss-of-
#: function rules — carries a platform caveat on ONT.
ONT_INDEL_UNCORROBORATED_FRACTION = _define(
    "ont_indel_uncorroborated_fraction", 0.166, SRC_COLPUS_2026, unit="fraction",
    note="ONT indel calls not corroborated by Illumina in the same isolates")
ONT_INDEL_CAVEAT = (
    "this call rests on an indel, and about 16.6% of ONT indel calls were "
    "uncorroborated by Illumina in the 508-isolate comparison"
)

#: SOURCE: design §7. R10.4.1 chemistry with Dorado `sup` basecalling is the
#: minimum credible ONT configuration; `fast` is not acceptable. Mjolnir cannot
#: always detect this from the reads, so it is stated as a requirement in the
#: report rather than silently assumed.
ONT_MINIMUM_CONFIGURATION = _define(
    "ont_minimum_configuration", "R10.4.1 + Dorado sup", SRC_DESIGN,
    note="fast-basecalled ONT data is outside the validated envelope")

#: SOURCE: design §7. FASTA input has no allele fractions at all.
FASTA_CAPABILITY_LOSS = _define(
    "fasta_capability_loss",
    "assembly input carries no allele fractions: heteroresistance, mixed-infection "
    "detection and the heterozygosity-based contamination metrics are unavailable, "
    "which is a capability loss and not a clean result.",
    SRC_DESIGN)

PLATFORM_CAVEATS: Dict[str, Tuple[str, ...]] = {
    PLATFORM_ILLUMINA: (),
    PLATFORM_ONT: (ONT_MINOR_VARIANT_CAVEAT, ONT_FBIC_CAVEAT, ONT_INDEL_CAVEAT),
    PLATFORM_FASTA: (FASTA_CAPABILITY_LOSS,),
}


def min_reads_for(platform: str) -> int:
    """Minimum reads supporting a variant, by platform.

    Raises for FASTA rather than returning zero: an assembly has no reads, and a
    caller that asks this question about a FASTA has a bug that must surface
    here rather than as a threshold of 0 applied to nothing.
    """
    plat = normalise_platform(platform)
    if plat == PLATFORM_ILLUMINA:
        return MIN_READS_ILLUMINA
    if plat == PLATFORM_ONT:
        return MIN_READS_ONT
    raise MjolnirError(
        "no read-support threshold applies to {0} input: an assembly has no "
        "read evidence".format(PLATFORM_FASTA)
    )


def is_major_variant(allele_fraction: Optional[float]) -> Optional[bool]:
    """Whether a variant is major. None in, None out — never a default of True."""
    if allele_fraction is None:
        return None
    return float(allele_fraction) >= MAJOR_VARIANT_FRACTION


def platform_caveats(platform: str) -> Tuple[str, ...]:
    """The consequences the report must state for this platform (design §7)."""
    return PLATFORM_CAVEATS[normalise_platform(platform)]


def is_suppressed_on_platform(gene: str, drug: str, platform: str) -> Optional[str]:
    """The reason a gene-drug call is suppressed on this platform, or None."""
    if normalise_platform(platform) != PLATFORM_ONT:
        return None
    for sup_gene, sup_drug in ONT_SUPPRESSED_GENE_DRUG:
        if gene.lower() == sup_gene.lower() and drug.lower() == sup_drug.lower():
            return ONT_FBIC_CAVEAT
    return None


# ---------------------------------------------------------------------------
# MTBseq legacy defaults, kept so a Mjolnir run can be compared to an MTBseq one
# ---------------------------------------------------------------------------
#
# SOURCE for all eight: MTBseq MANUAL.md. These are not Mjolnir's operating
# thresholds — MIN_READS_* and MIN_DEPTH above are — they exist so that
# `--compat mtbseq` reproduces the legacy filter stack and the two tools can be
# compared on the same reads. The prior chimaera run on this machine shows the
# convention in its own filenames: try_joint_cf5_cr5_fr75_ph4_samples2.tab.

MTBSEQ_MINCOVF = _define(
    "mtbseq_mincovf", 4, SRC_MTBSEQ_MANUAL, unit="reads",
    note="minimum forward-strand reads covering a position")
MTBSEQ_MINCOVR = _define(
    "mtbseq_mincovr", 4, SRC_MTBSEQ_MANUAL, unit="reads",
    note="minimum reverse-strand reads covering a position")
MTBSEQ_MINPHRED = _define(
    "mtbseq_minphred", 4, SRC_MTBSEQ_MANUAL, unit="reads",
    note="minimum reads with phred score at or above the quality cutoff")
MTBSEQ_MINFREQ = _define(
    "mtbseq_minfreq", 75, SRC_MTBSEQ_MANUAL, unit="percent",
    note="majority rule: allele frequency for a variant call. This is the filter "
         "that discards minority signal; Mjolnir surfaces what it would drop")
MTBSEQ_UNAMBIG = _define(
    "mtbseq_unambig", 95, SRC_MTBSEQ_MANUAL, unit="percent",
    note="unambiguous-base threshold; MTBseq's de-facto heterozygosity filter")
MTBSEQ_MINBQUAL = _define(
    "mtbseq_minbqual", 13, SRC_MTBSEQ_MANUAL, unit="phred",
    note="minimum base quality considered in the pileup")
MTBSEQ_WINDOW = _define(
    "mtbseq_window", 12, SRC_MTBSEQ_MANUAL, unit="bp",
    note="window within which a second variant invalidates the first")
MTBSEQ_DISTANCE = _define(
    "mtbseq_distance", 12, SRC_MTBSEQ_MANUAL, unit="SNPs",
    note="TBgroup clustering distance")

MTBSEQ_DEFAULTS: Dict[str, Any] = {
    "mincovf": MTBSEQ_MINCOVF,
    "mincovr": MTBSEQ_MINCOVR,
    "minphred": MTBSEQ_MINPHRED,
    "minfreq": MTBSEQ_MINFREQ,
    "unambig": MTBSEQ_UNAMBIG,
    "minbqual": MTBSEQ_MINBQUAL,
    "window": MTBSEQ_WINDOW,
    "distance": MTBSEQ_DISTANCE,
}

# ---------------------------------------------------------------------------
# Contamination (design §8)
# ---------------------------------------------------------------------------

#: SOURCE: Sobkowiak et al. 2018, BMC Genomics (MixInfect). Heterozygous-site
#: filters: base quality at or above 20, depth at or above 10.
HET_MIN_QUAL = _define(
    "het_min_qual", 20, SRC_MIXINFECT, unit="phred",
    note="minimum quality for a site to count as heterozygous")
HET_MIN_DEPTH = _define(
    "het_min_depth", 10, SRC_MIXINFECT, unit="x",
    note="minimum depth for a site to count as heterozygous")

#: SOURCE: Mjolnir policy, informed by MixInfect. MixInfect itself fits a model
#: rather than applying a cutoff, so these two tiers are Mjolnir's, and the
#: report says so: below the first, single-strain; between them, possible
#: mixture; above the second, mixed.
HET_SNP_FRACTION_WARN = _define(
    "het_snp_fraction_warn", 0.010, SRC_POLICY, unit="fraction",
    note="genome-wide heterozygous-SNP fraction above which a mixture is possible")
HET_SNP_FRACTION_MIXED = _define(
    "het_snp_fraction_mixed", 0.050, SRC_POLICY, unit="fraction",
    note="genome-wide heterozygous-SNP fraction above which a mixture is called")

#: SOURCE: design §8.1 — F2/F47 minor-allele-frequency statistics over
#: lineage-defining SNP sets, the operationally validated route to separating
#: mixed infection from cross-contamination via batch patterns. The primary
#: citation for the F2/F47 definitions has not been checked on this machine, so
#: the thresholds are marked unverified and the report marks them too.
F2_MIXTURE_THRESHOLD = _define(
    "f2_mixture_threshold", 0.05, SRC_DESIGN, unit="fraction",
    note="F2 statistic above which a mixture is suspected", verified=False)
F47_MIXTURE_THRESHOLD = _define(
    "f47_mixture_threshold", 0.05, SRC_DESIGN, unit="fraction",
    note="F47 statistic above which a mixture is suspected", verified=False)

#: SOURCE: design §8, item 4. MTBseq's unambiguous-base fraction, reported
#: rather than silently applied.
MIN_UNAMBIGUOUS_FRACTION = _define(
    "min_unambiguous_fraction", 0.95, SRC_MTBSEQ_MANUAL, unit="fraction",
    note="fraction of called positions with an unambiguous majority allele")

#: SOURCE: design §8. The measured consequence of small amounts of contamination,
#: and the reason the headline is a validity verdict rather than a purity figure.
CONTAMINATION_EVIDENCE = _define(
    "contamination_evidence",
    "a sample that was 99.84% M. tuberculosis still produced 13 false-positive "
    "SNPs across 12 genes, and 5% M. avium contamination produced 3,325 "
    "false-positive variant SNPs; any gate at 1% or 5% is a coarse instrument",
    SRC_DESIGN)

# --- Kraken2: the documented refusal ---------------------------------------
#
# SOURCE: design §8. Kraken2's own default for --confidence is 0.0, and at 0.0 a
# contamination screen prints a tail of low-abundance NTM that is an artefact of
# k-mer promiscuity. Mjolnir refuses that default outright: kraken2_confidence()
# raises rather than accepting it, so no code path can reach a 0.0 screen.
#
# The second refusal is larger. Measured Kraken2 sensitivity for M. tuberculosis
# reads with the standard database is 0.0731 on real Illumina data — about 93% of
# true target reads are unclassified or misassigned — against ~0.97 with a
# mycobacterial pangenome database. A standard or capped index therefore cannot
# support any statement about mycobacterial purity, and a screen run against one
# is reported as uninformative rather than as clean.

KRAKEN2_REFUSED_DEFAULT = _define(
    "kraken2_refused_default", 0.0, SRC_DESIGN,
    note="Kraken2's own --confidence default, which Mjolnir refuses to use")
KRAKEN2_MIN_CONFIDENCE = _define(
    "kraken2_min_confidence", 0.10, SRC_POLICY, unit="confidence",
    note="Mjolnir's --confidence floor for any Kraken2 run it performs")
KRAKEN2_MTB_SENSITIVITY_STANDARD = _define(
    "kraken2_mtb_sensitivity_standard", 0.0731, SRC_DESIGN, unit="fraction",
    note="measured Kraken2 sensitivity for M. tuberculosis reads, standard index",
    verified=False)
KRAKEN2_MTB_SENSITIVITY_PANGENOME = _define(
    "kraken2_mtb_sensitivity_pangenome", 0.97, SRC_DESIGN, unit="fraction",
    note="measured Kraken2 sensitivity with a mycobacterial pangenome index",
    verified=False)

#: Index names that are known not to be mycobacterial pangenome databases. Used
#: by :func:`kraken2_index_informative` when the index carries no manifest.
KRAKEN2_UNINFORMATIVE_INDEX_HINTS: Tuple[str, ...] = (
    "standard", "capped", "minikraken", "pluspf", "pluspfp", "plusp",
    "std8", "std16", "std_8", "std_16", "viral", "core_nt", "nt",
)

KRAKEN2_UNINFORMATIVE_TEXT = (
    "the Kraken2 index available here is not a mycobacterial pangenome database. "
    "Measured sensitivity for M. tuberculosis reads with a standard index is "
    "0.0731 on real Illumina data, so this screen cannot support a statement "
    "about mycobacterial purity and no such statement is made."
)


def kraken2_confidence(value: Optional[float] = None) -> float:
    """Validate a Kraken2 confidence, refusing 0.0.

    This is a refusal, not a preference. Kraken2's default of 0.0 makes a
    contamination screen report a tail of low-abundance NTM as if they were
    co-infections, and the design forbids that output existing at all — so the
    check lives at the only point every caller must pass through.
    """
    if value is None:
        return KRAKEN2_MIN_CONFIDENCE
    conf = float(value)
    if conf <= KRAKEN2_REFUSED_DEFAULT:
        raise MjolnirError(
            "Kraken2 --confidence {0} is refused: at Kraken2's 0.0 default a "
            "contamination screen reports a tail of low-abundance NTM as "
            "co-infections. Use at least {1} (--kraken2-confidence)."
            .format(conf, KRAKEN2_MIN_CONFIDENCE)
        )
    if conf > 1.0:
        raise MjolnirError(
            "Kraken2 --confidence must be between 0 and 1 (got {0})".format(conf))
    return conf


def kraken2_index_informative(db_dir: Optional[Any]) -> Tuple[bool, str]:
    """Whether this Kraken2 index can say anything about mycobacterial purity.

    An index may declare itself by shipping ``mjolnir_index.json`` with
    ``{"mycobacterial_pangenome": true}`` beside its ``hash.k2d`` — that is what
    ``mjolnir db`` writes when it builds one. Without that declaration the name
    is the only evidence available, and the answer defaults to "not informative",
    because the cost of being wrong in that direction is a report that says
    "no contamination detected" about a screen that could not have detected it.
    """
    if db_dir is None:
        return False, "no Kraken2 index configured"
    path = Path(str(db_dir)).expanduser()
    if not path.exists():
        return False, "Kraken2 index not found at {0}".format(path)
    manifest = path / "mjolnir_index.json"
    if manifest.exists():
        try:
            declared = json.loads(manifest.read_text())
        except (OSError, ValueError) as exc:
            raise MjolnirError(
                "Kraken2 index manifest {0} could not be read: {1}".format(manifest, exc)
            )
        if bool(declared.get("mycobacterial_pangenome")):
            return True, "mycobacterial pangenome index declared by {0}".format(manifest)
        return False, "{0} declares this index is not a mycobacterial pangenome".format(manifest)
    lowered = path.name.lower()
    for hint in KRAKEN2_UNINFORMATIVE_INDEX_HINTS:
        if hint in lowered:
            return False, KRAKEN2_UNINFORMATIVE_TEXT
    return False, (
        "the Kraken2 index at {0} carries no mjolnir_index.json declaring it a "
        "mycobacterial pangenome database, so its output is treated as "
        "uninformative for mycobacterial purity. {1}".format(path, KRAKEN2_UNINFORMATIVE_TEXT)
    )


# ---------------------------------------------------------------------------
# Cohort (design §9)
# ---------------------------------------------------------------------------

#: SOURCE: Walker et al. 2013, Lancet Infect Dis. The two TB conventions: 5 SNPs
#: for a recent-transmission cluster, 12 SNPs for the wider epidemiological link.
CLUSTER_SNP_STRICT = _define(
    "cluster_snp_strict", 5, SRC_WALKER_2013, unit="SNPs",
    note="recent transmission")
CLUSTER_SNP_RELAXED = _define(
    "cluster_snp_relaxed", 12, SRC_WALKER_2013, unit="SNPs",
    note="wider epidemiological link; also MTBseq's TBgroup default")

#: SOURCE: prior M. chimaera run on this machine (MTBseq 1.0.3, 2022-11-11,
#: M_Chimaera_TN/samp), which was run with --distance 6. Recorded because it is
#: the comparison baseline for the NTM validation set — it is not a TB threshold
#: and is not a default.
CHIMAERA_LOCAL_DISTANCE = _define(
    "chimaera_local_distance", 6, "prior local MTBseq M. chimaera run (--distance 6)",
    unit="SNPs",
    note="local precedent for the NTM outbreak set, not a published threshold")

#: Default clustering distance. The relaxed TB convention, because reporting a
#: cluster that is not one is a smaller harm than missing a link — and because
#: the threshold is a flag, printed with its basis, not a constant.
DEFAULT_CLUSTER_DISTANCE = _define(
    "default_cluster_distance", CLUSTER_SNP_RELAXED, SRC_WALKER_2013, unit="SNPs")

CLUSTER_THRESHOLD_BASIS: Dict[int, str] = {
    CLUSTER_SNP_STRICT: "5 SNPs - recent transmission ({0})".format(SRC_WALKER_2013),
    CLUSTER_SNP_RELAXED: "12 SNPs - epidemiological link ({0})".format(SRC_WALKER_2013),
    CHIMAERA_LOCAL_DISTANCE: (
        "6 SNPs - the value used by the prior local MTBseq M. chimaera run, kept "
        "for comparability with that baseline; not a published TB threshold"),
}


def cluster_threshold_basis(distance: int) -> str:
    """The sentence printed beside a cluster, explaining where the number came from."""
    known = CLUSTER_THRESHOLD_BASIS.get(int(distance))
    if known:
        return known
    return (
        "{0} SNPs - operator-supplied threshold with no published basis recorded; "
        "the TB conventions are {1} and {2} ({3})".format(
            distance, CLUSTER_SNP_STRICT, CLUSTER_SNP_RELAXED, SRC_WALKER_2013)
    )


#: SOURCE: Colpus et al. 2026 bioRxiv. Masking is mandatory before counting SNP
#: distances: 264,525 loci — about 6% of H37Rv — covering repetitive,
#: low-complexity and error-prone regions were masked, and only SNPs with no
#: other SNP within 12 bases were counted.
MASKED_LOCI_H37RV = _define(
    "masked_loci_h37rv", 264525, SRC_COLPUS_2026, unit="positions",
    note="repetitive, low-complexity and error-prone loci masked before counting")
MASKED_FRACTION_H37RV = _define(
    "masked_fraction_h37rv", 0.06, SRC_COLPUS_2026, unit="fraction")
SNP_PROXIMITY_WINDOW = _define(
    "snp_proximity_window", 12, SRC_COLPUS_2026, unit="bp",
    note="a SNP with another SNP within this distance is not counted")

#: SOURCE: Mjolnir policy, following tesseract-ai's cgMLST output. A distance
#: computed over less shared callable sequence than this is reported with the
#: denominator emphasised and is not used for clustering.
MIN_SHARED_CALLABLE_SITES = _define(
    "min_shared_callable_sites", 3_000_000, SRC_POLICY, unit="bp",
    note="shared callable sequence below which a pairwise distance is not "
         "comparable to the published SNP thresholds")

# ---------------------------------------------------------------------------
# WHO catalogue v2 (design §5)
# ---------------------------------------------------------------------------
#
# SOURCE for the grade strings: WHO-UCN-TB-2023.7 Catalogue_master_file. They are
# numeric-prefixed and use a SPACED ASCII HYPHEN. The published PDF renders an
# en-dash; matching on the PDF form silently matches nothing, which is a trap the
# design records as verified rather than hypothetical. Do not "tidy" these.

WHO_GRADE_1 = "1) Assoc w R"
WHO_GRADE_2 = "2) Assoc w R - Interim"
WHO_GRADE_3 = "3) Uncertain significance"
WHO_GRADE_4 = "4) Not assoc w R - Interim"
WHO_GRADE_5 = "5) Not assoc w R"

WHO_GRADES: Tuple[str, ...] = _define(
    "who_grades",
    (WHO_GRADE_1, WHO_GRADE_2, WHO_GRADE_3, WHO_GRADE_4, WHO_GRADE_5),
    SRC_WHO_V2,
    note="verbatim grade strings; spaced ASCII hyphen, not an en-dash")

#: WHO grade to Mjolnir call. WHO is the anchor: where WHO grades a variant, the
#: WHO grade is the Mjolnir call (§5.5 rule 2).
WHO_GRADE_TO_CALL: Dict[str, str] = {
    WHO_GRADE_1: CALL_R,
    WHO_GRADE_2: CALL_R_INTERIM,
    WHO_GRADE_3: CALL_UNCERTAIN,
    WHO_GRADE_4: CALL_S_INTERIM,
    WHO_GRADE_5: CALL_S,
}

#: Spellings seen in the wild that mean one of the five grades: the en-dash form
#: printed in the PDF, and the bare group numbers used by several downstream
#: tools. Anything not here is not silently coerced.
_GRADE_ALIASES: Dict[str, str] = {}
for _canonical in WHO_GRADES:
    _GRADE_ALIASES[_canonical.lower()] = _canonical
    _GRADE_ALIASES[_canonical.replace(" - ", " – ").lower()] = _canonical
    _GRADE_ALIASES[_canonical.replace(" - ", " — ").lower()] = _canonical
    _GRADE_ALIASES[_canonical.split(")", 1)[0]] = _canonical
    _GRADE_ALIASES["group " + _canonical.split(")", 1)[0]] = _canonical


def normalise_grade(grade: str) -> str:
    """Map a grade string onto one of the five canonical forms.

    Returns "" for anything unrecognised rather than guessing. A grade Mjolnir
    cannot place is not a grade, and a call built on a guessed grade would carry
    WHO's authority without WHO's evidence.
    """
    key = " ".join(str(grade or "").split()).lower()
    if not key:
        return ""
    found = _GRADE_ALIASES.get(key)
    if found:
        return found
    LOG.debug("unrecognised WHO grade string %r", grade)
    return ""


def call_for_grade(grade: str) -> str:
    """The Mjolnir call a WHO grade maps to; ``no-call`` when unrecognised."""
    return WHO_GRADE_TO_CALL.get(normalise_grade(grade), CALL_NO_CALL)


#: SOURCE: WHO-UCN-TB-2023.7 header layout. The xlsx header is on row 3; rows 1-2
#: are a merged banner, so a naive read_excel yields garbage columns.
WHO_XLSX_HEADER_ROW = _define(
    "who_xlsx_header_row", 3, SRC_WHO_V2, unit="1-based row",
    note="pandas read_excel needs header=2 for this")
WHO_XLSX_MASTER_SHEET = _define(
    "who_xlsx_master_sheet", "Catalogue_master_file", SRC_WHO_V2)
WHO_XLSX_COORDINATES_SHEET = _define(
    "who_xlsx_coordinates_sheet", "Genomic_coordinates", SRC_WHO_V2)

#: SOURCE: WHO-UCN-TB-2023.7, verified counts from the design. Loaders assert
#: against these and say so when the file they were given does not match, since
#: the repo's .txt master file is a different and smaller thing (40,178 rows, 14
#: drugs, Streptomycin absent entirely).
WHO_V2_EXPECTED_ROWS = _define(
    "who_v2_expected_rows", 48152, SRC_WHO_V2, unit="variant-drug rows")
WHO_V2_EXPECTED_VARIANTS = _define(
    "who_v2_expected_variants", 30699, SRC_WHO_V2, unit="unique variants")
WHO_V2_EXPECTED_GENES = _define(
    "who_v2_expected_genes", 65, SRC_WHO_V2, unit="genes")
WHO_V2_EXPECTED_DRUGS = _define(
    "who_v2_expected_drugs", 15, SRC_WHO_V2, unit="drugs")

#: SOURCE: WHO-UCN-TB-2023.7. 38,884 of 48,152 rows carry this literal string in
#: the `genomic position` column instead of a coordinate. It is a decoy: the real
#: coordinates are on the Genomic_coordinates sheet.
WHO_COORDINATE_DECOY = _define(
    "who_coordinate_decoy", '(see "Genomic_coordinates" sheet)', SRC_WHO_V2,
    note="literal contents of the decoy `genomic position` column")

#: SOURCE: WHO-UCN-TB-2023.7 VCF INFO encoding. MNVs are decomposed and one
#: genomic variant maps to several graded variants, joined by this separator.
WHO_MNV_SEPARATOR = _define("who_mnv_separator", "&", SRC_WHO_V2)

#: Catalogue identities. Spelled once so three modules cannot disagree about
#: whether the anchor is called "WHO", "WHOv2" or "WHO v2".
CATALOGUE_WHO = "WHO v2"
CATALOGUE_MTBSEQ = "MTBseq"
CATALOGUE_TBDB = "tbdb"
CATALOGUES: Tuple[str, ...] = (CATALOGUE_WHO, CATALOGUE_MTBSEQ, CATALOGUE_TBDB)

#: The anchor catalogue (§5.5 rule 2): the only source with a published,
#: systematically-derived grading.
ANCHOR_CATALOGUE = _define(
    "anchor_catalogue", CATALOGUE_WHO, SRC_DESIGN,
    note="where WHO grades a variant, the WHO grade is the Mjolnir call")

#: MTBseq's list has no grading at all, so it can only ever produce R or
#: no-call. The report states the asymmetry rather than hiding it.
MTBSEQ_ASYMMETRY_NOTE = (
    "MTBseq's resistance list is flat: it has no confidence grading, so it can "
    "only contribute R or no-call, and its agreement with a WHO grade is weaker "
    "evidence than it looks"
)

# --- Drugs -----------------------------------------------------------------
#
# SOURCE: WHO-UCN-TB-2023.7 covers 15 drugs. Capreomycin, graded in the 1st
# edition, is absent from the 2nd. This tuple fixes display order and gives the
# report a stable column layout; catalogues.py must take the authoritative drug
# set from the xlsx it actually loaded, because a 3rd edition was called for on
# 2024-08-26 and this list must not be what decides which drugs exist.

DRUGS: Tuple[str, ...] = _define(
    "drugs",
    ("Rifampicin", "Isoniazid", "Ethambutol", "Pyrazinamide",
     "Levofloxacin", "Moxifloxacin", "Bedaquiline", "Linezolid",
     "Clofazimine", "Delamanid", "Pretomanid", "Amikacin",
     "Kanamycin", "Streptomycin", "Ethionamide"),
    SRC_WHO_V2, note="display order; the loaded catalogue is authoritative")

DRUG_CODES: Dict[str, str] = {
    "Rifampicin": "RIF", "Isoniazid": "INH", "Ethambutol": "EMB",
    "Pyrazinamide": "PZA", "Levofloxacin": "LFX", "Moxifloxacin": "MFX",
    "Bedaquiline": "BDQ", "Linezolid": "LZD", "Clofazimine": "CFZ",
    "Delamanid": "DLM", "Pretomanid": "PMD", "Amikacin": "AMK",
    "Kanamycin": "KAN", "Streptomycin": "STM", "Ethionamide": "ETO",
    "Capreomycin": "CAP", "Clarithromycin": "CLR",
    # Graded by MTBseq but not by the WHO catalogue v2. They are here because
    # dropping them loses 30 real resistance rows: MTBseq carries 22 for
    # para-aminosalicylic acid and 8 for cycloserine, and a determinant that
    # reaches no drug column is indistinguishable from one that was never found.
    "Para-aminosalicylic acid": "PAS", "Cycloserine": "CS",
}

#: Drug *classes* that a catalogue may name where an agent is expected. MTBseq
#: writes ``fluoroquinolones (FQ)`` on 22 rows, which is a statement about the
#: class - gyrA mutations confer it - and WHO grades those same genes for both
#: fluoroquinolones it covers. Expanding is therefore correct biologically, but
#: the report must say the source graded a class rather than these two agents,
#: so the expansion is recorded rather than applied silently.
DRUG_CLASSES: Dict[str, Tuple[str, ...]] = {
    "fluoroquinolones": ("Levofloxacin", "Moxifloxacin"),
    "fq": ("Levofloxacin", "Moxifloxacin"),
    "fluoroquinolone": ("Levofloxacin", "Moxifloxacin"),
}

#: Spellings the three catalogues use for the same drug. Cross-catalogue
#: consensus joins on the drug name, so an unnormalised "rifampin" would look
#: like a drug nobody else called.
DRUG_ALIASES: Dict[str, str] = {}
for _drug, _code in DRUG_CODES.items():
    DRUG_ALIASES[_drug.lower()] = _drug
    DRUG_ALIASES[_code.lower()] = _drug
DRUG_ALIASES.update({
    "rifampin": "Rifampicin", "rmp": "Rifampicin", "rif": "Rifampicin",
    "inh": "Isoniazid", "isoniazide": "Isoniazid",
    "emb": "Ethambutol", "pza": "Pyrazinamide",
    "levofloxacine": "Levofloxacin", "lev": "Levofloxacin", "lvx": "Levofloxacin",
    "moxifloxacine": "Moxifloxacin", "mxf": "Moxifloxacin",
    "bdq": "Bedaquiline", "lzd": "Linezolid", "cfz": "Clofazimine",
    "dlm": "Delamanid", "pmd": "Pretomanid", "pa": "Pretomanid",
    "amk": "Amikacin", "kan": "Kanamycin", "km": "Kanamycin",
    "streptomycine": "Streptomycin", "str": "Streptomycin", "sm": "Streptomycin",
    "ethionamid": "Ethionamide", "eth": "Ethionamide", "eto": "Ethionamide",
    "clarithromycine": "Clarithromycin", "clr": "Clarithromycin",
    "pas": "Para-aminosalicylic acid",
    "para-aminosalicylic acid": "Para-aminosalicylic acid",
    "para aminosalicylic acid": "Para-aminosalicylic acid",
    "p-aminosalicylic acid": "Para-aminosalicylic acid",
    "aminosalicylic acid": "Para-aminosalicylic acid",
    "cs": "Cycloserine", "cycloserine": "Cycloserine",
    "d-cycloserine": "Cycloserine", "dcs": "Cycloserine",
    "terizidone": "Cycloserine",
})


def normalise_drug(name: str) -> str:
    """Canonical drug name, or the input title-cased when it is unknown.

    Unknown drugs are passed through rather than dropped: a catalogue edition
    that adds a drug must still be reportable, and silently discarding its rows
    would be the worse failure.
    """
    key = " ".join(str(name or "").split()).lower()
    if not key:
        return ""
    return DRUG_ALIASES.get(key, str(name).strip())


def drug_code(name: str) -> str:
    canonical = normalise_drug(name)
    return DRUG_CODES.get(canonical, canonical[:3].upper())


# --- Rules beyond the table (design §5.4) ----------------------------------
#
# SOURCE for all of these: WHO-UCN-TB-2023.7, "additional grading rules" — the
# catalogue is three components, not one table, and lookup alone is incorrect.

#: The rifampicin resistance-determining region, inclusive codon range. Any
#: non-synonymous mutation or indel inside it grades as Group 2.
RPOB_RRDR_CODONS: Tuple[int, int] = _define(
    "rpob_rrdr_codons", (426, 452), SRC_WHO_V2, unit="codon",
    note="RRDR; any non-synonymous change or indel inside it is Group 2")

#: The four borderline rpoB mutations that are Group 1 by rule. They sit at or
#: near the edge of phenotypic detectability and are the classic source of
#: "phenotypically susceptible, genotypically resistant" disagreement.
RPOB_BORDERLINE: Tuple[str, ...] = _define(
    "rpob_borderline", ("Leu430Pro", "His445Asn", "His445Ser", "Ile491Phe"),
    SRC_WHO_V2, note="Group 1 by rule; Ile491Phe sits outside the RRDR")

#: Genes whose loss of function grades as Group 2 for the named drug.
LOF_GROUP2_GENES: Dict[str, Tuple[str, ...]] = _define(
    "lof_group2_genes",
    {
        "Isoniazid": ("katG",),
        "Pyrazinamide": ("pncA",),
        "Bedaquiline": ("Rv0678", "pepQ"),
        "Clofazimine": ("Rv0678", "pepQ"),
        "Delamanid": ("ddn", "fbiA", "fbiB", "fbiC", "fbiD", "fgd1"),
        "Pretomanid": ("ddn", "fbiA", "fbiB", "fbiC", "fbiD", "fgd1"),
        "Ethionamide": ("ethA",),
    },
    SRC_WHO_V2, note="loss of function in these genes is Group 2 for the drug")

#: Consequences Mjolnir treats as loss of function.
LOF_EFFECTS: Tuple[str, ...] = _define(
    "lof_effects",
    ("frameshift", "stop_gained", "start_lost", "feature_ablation",
     "transcript_ablation", "gene_deletion"),
    SRC_WHO_V2)

#: Epistasis: a loss-of-function in the first gene abrogates the effect of
#: variants in the second, for the listed drugs. Encoded in WHO's Comment column
#: and applied here as a suppression step — with the suppression stated in the
#: report, never applied silently.
EPISTASIS_RULES: Tuple[Dict[str, Any], ...] = _define(
    "epistasis_rules",
    (
        {"suppressor_gene": "mmpL5", "suppressor_effect": "lof",
         "suppressed_gene": "Rv0678", "drugs": ("Bedaquiline", "Clofazimine"),
         "why": "mmpL5 loss of function abrogates the efflux phenotype that "
                "Rv0678 variants derepress"},
        {"suppressor_gene": "eis", "suppressor_effect": "lof",
         "suppressed_gene": "eis", "suppressed_region": "promoter",
         "drugs": ("Amikacin", "Kanamycin"),
         "why": "a coding loss of function in eis abrogates the effect of eis "
                "promoter mutations, which act by over-expressing the protein"},
    ),
    SRC_WHO_V2)

#: Any novel silent variant grades Group 4 (§5.4). Stated as a constant so the
#: rule and the report share one string.
SILENT_VARIANT_GRADE = _define(
    "silent_variant_grade", WHO_GRADE_4, SRC_WHO_V2,
    note="any novel synonymous variant is Not assoc w R - Interim")

#: Comment-column phrases that carry level of resistance and cross-resistance.
#: These come from the Comment column, not from the grade (§5.4).
COMMENT_HIGH_LEVEL = _define(
    "comment_high_level", "High-level resistance", SRC_WHO_V2)
CROSS_RESISTANCE_PAIRS: Tuple[Tuple[str, str], ...] = _define(
    "cross_resistance_pairs",
    (("Delamanid", "Pretomanid"), ("Bedaquiline", "Clofazimine")),
    SRC_WHO_V2, note="documented cross-resistance, reported on both drugs")

# ---------------------------------------------------------------------------
# Species and lineage (design §6)
# ---------------------------------------------------------------------------

#: SOURCE: Richter & Rossello-Mora 2009; Jain et al. 2018. The prokaryotic
#: species boundary. Below this, no species claim is made.
ANI_SPECIES_FLOOR = _define(
    "ani_species_floor", 95.0, SRC_ANI_SPECIES, unit="percent ANI",
    note="ANI below which Mjolnir will not name a species")

#: SOURCE: Mjolnir policy, informed by the ~80% ANI that separates Mycobacterium
#: from its neighbours. Not a species boundary and never used as one: this is the
#: floor below which a hit is not close enough for its genome to serve as a
#: mapping reference for the query. A query under it is not something the
#: mycobacterial reference set can speak for, and calling variants against a
#: genome that distant would put every coordinate in the wrong organism without
#: saying so.
#: SOURCE: WHO catalogue v2 practice. Its graded promoter variants run to
#: c.-669 (mshA) and the catalogue's own upstream rows sit inside this window,
#: so a narrower one would drop Group 1 determinants - eis c.-14C>T and
#: inhA c.-154G>A among them - and report them as intergenic.
PROMOTER_UPSTREAM_BP = _define(
    "promoter_upstream_bp", 1000, SRC_WHO_V2, unit="bp",
    note="how far before a start codon a variant is still named for that gene")

ANI_GENUS_FLOOR = _define(
    "ani_genus_floor", 80.0, "Mjolnir policy (design section 6)",
    unit="percent ANI",
    note="ANI below which a reference genome is too distant to map the query against")

#: SOURCE: Mjolnir policy. ANI is computed over an aligned fraction, and a very
#: high ANI over a small aligned fraction is not a species identification.
ANI_MIN_ALIGNED_FRACTION = _define(
    "ani_min_aligned_fraction", 0.60, SRC_POLICY, unit="fraction",
    note="aligned fraction below which an ANI value is not used for a species call")

#: SOURCE: design §6. MTBC members sit at 99.21-99.92% ANI of one another and are
#: later heterotypic synonyms of M. tuberculosis in NCBI taxonomy — M. bovis
#: (taxid 1765) has rank `no rank`. ANI therefore cannot resolve inside the
#: complex, and the honest outcome is "MTBC, not resolved below complex".
MTBC_INTRA_ANI_RANGE: Tuple[float, float] = _define(
    "mtbc_intra_ani_range", (99.21, 99.92), SRC_DESIGN, unit="percent ANI",
    note="ANI between MTBC members; inside this range ANI resolves nothing")
MTBC_UNRESOLVED_TEXT = (
    "M. tuberculosis complex members lie at 99.21-99.92% ANI of one another and "
    "are later heterotypic synonyms of M. tuberculosis, so ANI cannot resolve "
    "below complex; the member is called from lineage-defining SNPs instead"
)

#: SOURCE: Mjolnir policy. Within MAC, M. chimaera and M. intracellulare sit
#: above the 95% species boundary, so ANI alone does not separate them and
#: marker SNPs are required. This is the distinction the local outbreak data
#: needs, so it is a policy floor rather than an inherited one.
MAC_SPECIES_ANI_FLOOR = _define(
    "mac_species_ani_floor", 99.0, SRC_POLICY, unit="percent ANI",
    note="within MAC, ANI must reach this AND marker SNPs must agree before a "
         "species is named")

COMPLEX_MTBC = "MTBC"
COMPLEX_MAC = "MAC"
COMPLEX_ABSCESSUS = "M. abscessus complex"

#: SOURCE: design §6. A taxonomic read classifier is not a species identifier
#: for this genus, and Mjolnir never prints one of its rows as a species call.
SPECIES_METHOD_REFUSAL = _define(
    "species_method_refusal",
    "species identification is ANI-based; a taxonomic read-classifier row such as "
    "\"M. bovis 3.2%\" is not a species identification and is never printed as one",
    SRC_DESIGN)

#: SOURCE: Lipworth et al. 2019 (SNP-IT). M. bovis is defined by very few
#: phylogenetic SNPs, so the call is highly sensitive to coverage gaps and to
#: contamination — reported as a confidence caveat on the call itself.
MBOVIS_DEFINING_SNPS = _define(
    "mbovis_defining_snps", 23, SRC_SNPIT, unit="SNPs",
    note="so few that a coverage gap can lose or invent the call")
MBOVIS_CAVEAT = (
    "M. bovis is defined by very few phylogenetic SNPs (23 in SNP-IT), so this "
    "call is highly sensitive to coverage gaps and to contamination"
)

#: SOURCE: Mjolnir policy. Barcode-support thresholds: the fraction of a
#: lineage's defining sites that must be callable, and the fraction of those
#: that must carry the derived allele, before the lineage is called.
MIN_BARCODE_CALLABLE_FRACTION = _define(
    "min_barcode_callable_fraction", 0.80, SRC_POLICY, unit="fraction")
MIN_BARCODE_SUPPORT_FRACTION = _define(
    "min_barcode_support_fraction", 0.90, SRC_POLICY, unit="fraction")

#: SOURCE: design §6, following pathogen-profiler. Barcode genotyping is done
#: from a direct pileup, not from the variant caller, and on ONT the
#: highest-depth allele is taken at each barcode site.
BARCODE_FROM_PILEUP = _define(
    "barcode_from_pileup", True, SRC_DESIGN,
    note="barcode sites are genotyped from the pileup; on ONT the highest-depth "
         "allele wins")

#: SOURCE: design §6. BCG carries intrinsic pyrazinamide resistance, which is a
#: clinical fact and not a catalogue lookup, so it is flagged explicitly.
BCG_PZA_NOTE = _define(
    "bcg_pza_note",
    "BCG is intrinsically resistant to pyrazinamide; this is a property of the "
    "organism and is not derived from the mutation catalogue",
    SRC_DESIGN)

# ---------------------------------------------------------------------------
# NTM resistance (design §5.6)
# ---------------------------------------------------------------------------
#
# Not covered by any of the three MTBC catalogues, so these are implemented from
# the primary literature with explicit gene targets. Every NTM call names its
# supporting reference in the report, and where no evidence base exists for a
# species-drug pair the tool says so instead of guessing.

ERM41_SEQUEVAR_POSITION = _define(
    "erm41_sequevar_position", 28, SRC_BASTIAN_2011, unit="nucleotide position in erm(41)",
    note="T28 gives inducible macrolide resistance; C28 (T28C) is susceptible")
ERM41_INDUCIBLE_ALLELE = _define(
    "erm41_inducible_allele", "T", SRC_NASH_2009,
    note="T at position 28 of erm(41): inducible clarithromycin resistance")
ERM41_SUSCEPTIBLE_ALLELE = _define(
    "erm41_susceptible_allele", "C", SRC_BASTIAN_2011,
    note="C28, the T28C polymorphism: no inducible resistance")
ERM41_TRUNCATION_NOTE = _define(
    "erm41_truncation_note",
    "M. abscessus subsp. massiliense carries a truncated, non-functional erm(41) "
    "and is not inducibly macrolide-resistant",
    SRC_BASTIAN_2011)

RRL_MACROLIDE_POSITIONS: Tuple[int, ...] = _define(
    "rrl_macrolide_positions", (2058, 2059), SRC_WALLACE_1996,
    unit="E. coli numbering",
    note="acquired constitutive macrolide resistance in rrl (23S rRNA)")
RRS_AMIKACIN_POSITIONS: Tuple[int, ...] = _define(
    "rrs_amikacin_positions", (1408,), SRC_PRAMMANANAN_1998,
    unit="E. coli numbering",
    note="rrs 1408 confers amikacin resistance; neighbouring positions 1406 and "
         "1409-1491 are reported as uncertain rather than as resistance")
RRS_AMIKACIN_NEIGHBOURS: Tuple[int, ...] = _define(
    "rrs_amikacin_neighbours", (1406, 1409, 1491), SRC_PRAMMANANAN_1998,
    unit="E. coli numbering",
    note="reported, graded uncertain: the evidence base is thinner than for 1408",
    verified=False)

#: Which NTM species-drug pairs Mjolnir has an evidence base for. A pair absent
#: from this table is answered with "no evidence base", never with a call.
NTM_TARGETS: Dict[str, Dict[str, Dict[str, Any]]] = {
    "Mycobacteroides abscessus": {
        "Clarithromycin": {"genes": ("erm(41)", "rrl"), "source": SRC_BASTIAN_2011},
        "Amikacin": {"genes": ("rrs",), "source": SRC_PRAMMANANAN_1998},
    },
    "Mycobacterium chimaera": {
        "Clarithromycin": {"genes": ("rrl",), "source": SRC_WALLACE_1996},
        "Amikacin": {"genes": ("rrs",), "source": SRC_PRAMMANANAN_1998},
    },
    "Mycobacterium intracellulare": {
        "Clarithromycin": {"genes": ("rrl",), "source": SRC_WALLACE_1996},
        "Amikacin": {"genes": ("rrs",), "source": SRC_PRAMMANANAN_1998},
    },
    "Mycobacterium avium": {
        "Clarithromycin": {"genes": ("rrl",), "source": SRC_WALLACE_1996},
        "Amikacin": {"genes": ("rrs",), "source": SRC_PRAMMANANAN_1998},
    },
}

#: Aliases so a species call spelled either way finds its targets. The genus was
#: split — M. abscessus is now Mycobacteroides abscessus — and both names are in
#: current use in the databases Mjolnir reads.
NTM_SPECIES_ALIASES: Dict[str, str] = {
    "mycobacterium abscessus": "Mycobacteroides abscessus",
    "mycobacteroides abscessus": "Mycobacteroides abscessus",
    "mycobacterium chimaera": "Mycobacterium chimaera",
    "mycobacterium intracellulare": "Mycobacterium intracellulare",
    "mycobacterium avium": "Mycobacterium avium",
}

NTM_NO_EVIDENCE_TEXT = (
    "no published evidence base is implemented for {species} and {drug}, so no "
    "genotypic prediction is made; this is an absence of evidence, not a "
    "prediction of susceptibility"
)


def ntm_targets(species: str, drug: str) -> Optional[Dict[str, Any]]:
    """Gene targets and the citation for one NTM species-drug pair, or None.

    None means Mjolnir has no evidence base for the pair, and the caller is
    expected to emit :data:`NTM_NO_EVIDENCE_TEXT` rather than a susceptible call.
    """
    key = " ".join(str(species or "").split()).lower()
    canonical = NTM_SPECIES_ALIASES.get(key, species)
    entry = NTM_TARGETS.get(canonical)
    if not entry:
        return None
    return entry.get(normalise_drug(drug))


# ---------------------------------------------------------------------------
# Environment and run configuration
# ---------------------------------------------------------------------------

DB_ENV_VAR = "MJOLNIR_DB"
LLM_HOST_ENV_VAR = "MJOLNIR_LLM_HOST"
LLM_MODEL_ENV_VAR = "MJOLNIR_LLM_MODEL"
THREADS_ENV_VAR = "MJOLNIR_THREADS"
KRAKEN2_DB_ENV_VAR = "MJOLNIR_KRAKEN2_DB"

#: Report profiles. `clinical` puts the drug table and validity verdict first;
#: `research` reorders the same content, annexes forward.
PROFILES: Tuple[str, ...] = ("clinical", "research")


def default_db_dir() -> Path:
    """Resolve the database root: ``$MJOLNIR_DB`` -> ``~/.mjolnir/db`` -> bundled."""
    env = os.environ.get(DB_ENV_VAR)
    if env:
        return Path(env).expanduser().resolve()
    user_dir = Path.home() / ".mjolnir" / "db"
    if user_dir.exists():
        return user_dir
    bundled = BUNDLED_DATA / "db"
    if bundled.exists():
        return bundled
    return user_dir


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError:
        raise MjolnirError("{0} must be an integer (got {1!r})".format(name, raw))


@dataclass
class Config:
    """Everything a run needs that is not a per-sample input.

    Thresholds default to the module constants above, so a run that overrides
    nothing is running published numbers. ``explicit`` records which ones the
    operator changed, because a report must be able to say "this run used a
    depth floor of 8, not the published 25" — a changed threshold that is not
    announced is the same problem as a threshold with no source.
    """

    db_dir: Path = field(default_factory=default_db_dir)
    out_dir: Path = field(default_factory=lambda: Path("mjolnir_out"))
    threads: int = field(default_factory=cpu_count)
    #: None means "detect from the input" (seqio.py).
    platform: Optional[str] = None
    reference: Optional[Path] = None
    profile: str = "clinical"

    min_depth: int = MIN_DEPTH
    degraded_depth_floor: int = DEGRADED_DEPTH_FLOOR
    min_breadth: float = MIN_BREADTH
    min_reads_illumina: int = MIN_READS_ILLUMINA
    min_reads_ont: int = MIN_READS_ONT
    major_variant_fraction: float = MAJOR_VARIANT_FRACTION
    min_minor_variant_fraction: float = MIN_MINOR_VARIANT_FRACTION

    cluster_distance: int = DEFAULT_CLUSTER_DISTANCE
    mask_bed: Optional[Path] = None

    kraken2_db: Optional[Path] = None
    kraken2_confidence: float = KRAKEN2_MIN_CONFIDENCE

    llm_host: str = ""
    llm_model: str = ""
    #: Off means the report is rule-only and says so. It always runs without the
    #: model; the model never gates the result.
    use_llm: bool = True

    keep_temp: bool = False
    tmp_dir: Optional[Path] = None
    #: MTBseq-comparable filter stack, for benchmarking against a legacy run.
    mtbseq_compat: bool = False
    explicit: set = field(default_factory=set)

    def __post_init__(self) -> None:
        self.db_dir = Path(self.db_dir)
        self.out_dir = Path(self.out_dir)
        if self.reference is not None:
            self.reference = Path(self.reference)
        if self.mask_bed is not None:
            self.mask_bed = Path(self.mask_bed)
        if self.tmp_dir is not None:
            self.tmp_dir = Path(self.tmp_dir)
        if self.kraken2_db is None:
            env_db = os.environ.get(KRAKEN2_DB_ENV_VAR)
            if env_db:
                self.kraken2_db = Path(env_db).expanduser()
        elif not isinstance(self.kraken2_db, Path):
            self.kraken2_db = Path(str(self.kraken2_db)).expanduser()
        if not self.llm_host:
            self.llm_host = os.environ.get(LLM_HOST_ENV_VAR, "")
        if not self.llm_model:
            self.llm_model = os.environ.get(LLM_MODEL_ENV_VAR, "")
        if THREADS_ENV_VAR in os.environ and "threads" not in self.explicit:
            self.threads = _env_int(THREADS_ENV_VAR, self.threads)
        if self.platform is not None:
            self.platform = normalise_platform(self.platform)
        if self.threads < 1:
            self.threads = 1

    def set_explicit(self, name: str, value: Any) -> None:
        """Set a threshold and remember that the operator, not the paper, chose it."""
        if not hasattr(self, name):
            raise MjolnirError("no configuration field named {0!r}".format(name))
        setattr(self, name, value)
        self.explicit.add(name)

    def was_set(self, name: str) -> bool:
        return name in self.explicit

    def overridden_thresholds(self) -> Dict[str, Any]:
        """What this run changed away from its published defaults."""
        return dict((name, getattr(self, name)) for name in sorted(self.explicit))

    def min_reads(self, platform: str) -> int:
        plat = normalise_platform(platform)
        if plat == PLATFORM_ILLUMINA:
            return self.min_reads_illumina
        if plat == PLATFORM_ONT:
            return self.min_reads_ont
        return min_reads_for(plat)  # raises, with the explanation

    def validate(self) -> None:
        """Fail before the run rather than in the middle of it."""
        if self.profile not in PROFILES:
            raise MjolnirError(
                "--profile must be one of {0} (got {1!r})".format(
                    ", ".join(PROFILES), self.profile))
        if self.degraded_depth_floor < 1:
            raise MjolnirError("--degraded-depth-floor must be at least 1x")
        if self.min_depth < self.degraded_depth_floor:
            raise MjolnirError(
                "--min-depth ({0}) is below --degraded-depth-floor ({1}); the "
                "target depth cannot be lower than the floor at which a sample "
                "stops being callable".format(self.min_depth, self.degraded_depth_floor))
        for name in ("min_breadth", "major_variant_fraction", "min_minor_variant_fraction"):
            value = getattr(self, name)
            if not 0.0 <= float(value) <= 1.0:
                raise MjolnirError(
                    "--{0} must be a fraction between 0 and 1 (got {1})".format(
                        name.replace("_", "-"), value))
        if self.min_minor_variant_fraction >= self.major_variant_fraction:
            raise MjolnirError(
                "--min-minor-variant-fraction ({0}) must be below "
                "--major-variant-fraction ({1})".format(
                    self.min_minor_variant_fraction, self.major_variant_fraction))
        for name in ("min_reads_illumina", "min_reads_ont"):
            if getattr(self, name) < 1:
                raise MjolnirError("--{0} must be at least 1".format(name.replace("_", "-")))
        if self.cluster_distance < 0:
            raise MjolnirError("--distance must not be negative")
        # Refuses 0.0 outright; see kraken2_confidence().
        self.kraken2_confidence = kraken2_confidence(self.kraken2_confidence)

    def as_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        for key in ("db_dir", "out_dir", "reference", "mask_bed", "tmp_dir", "kraken2_db"):
            value = getattr(self, key)
            data[key] = str(value) if value is not None else None
        data["explicit"] = sorted(self.explicit)
        data["platform"] = self.platform
        return data
