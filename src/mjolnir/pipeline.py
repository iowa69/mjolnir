"""End-to-end orchestration: one sample from files to a finished result.

Every other module in Mjolnir answers one question and refuses to answer any
other. This module is where those answers are assembled into a
:class:`~mjolnir.records.SampleResult`, and its whole job is to make sure that
assembly cannot invent an answer nobody gave.

Three properties are deliberate and are what the rest of the file is shaped
around.

**A stage that cannot run removes a capability, it does not lower a standard.**
Every optional stage goes through :meth:`Pipeline._stage`, which catches
:class:`~mjolnir.utils.MjolnirError`, records the reason as a warning *and* as a
``measured=False`` :class:`~mjolnir.records.Check`, and returns ``None``. The
downstream stage then sees ``None`` and says what it could not do. Nothing falls
back to a cheaper method that produces a number of a different kind — a missing
ANI reference set means no species was identified, not a species guessed from a
read classifier, and a missing catalogue means no resistance call, not a
susceptible one.

**Order follows evidence, not convenience.** The species call is made twice from
one ANI run: once before mapping, because the reference has to be chosen from
something, and once afterwards with the MAC marker pileup, because within MAC
ANI alone may not name a species (§6). The two calls come from the same
:func:`~mjolnir.typing.species.ani_matches` result, so the reference a sample was
mapped to and the species printed beside it can never disagree.

**One pileup, not four.** The barcode sites, the MAC markers and the catalogue
coordinates are piled up in a single ``samtools mpileup`` pass over the union of
their positions, and the callers slice what they need out of it. The pileup runs
with ``-Q 20`` — MixInfect's own base-quality filter — which is why the
heterozygosity observations built here record that threshold as their quality:
the filter was applied per base before the counts existed, which is a stronger
statement than a per-site average, and it is stated wherever the number appears.

Reads are analysed as submitted. Mjolnir performs no trimming and no adapter
removal, exactly as MTBseq does not (§9b) — but unlike MTBseq it says so, in a
caveat on every read-derived result, because "we did not trim" and "trimming was
not necessary" are different statements.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass, fields as dataclass_fields
from datetime import datetime
from pathlib import Path
from typing import (Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set,
                    Tuple)

from . import __version__, config, shell
from .config import Config
from .contamination.heterozygosity import LineageSite, SiteObservation
from .contamination.purity import (INTENDED_USES, TaxonomicScreen,
                                   assess_contamination, evaluate_kraken2_screen,
                                   no_screen)
from .cohort.cluster import cluster_samples
from .cohort.distance import Mask, distance_matrix, load_mask
from .cohort.joint import Regions, SampleVariants, build_joint_table, cohort_size_check
from .db import fetch as db_fetch
from .db import registry as db_registry
from .engines.call import call_variants
from .engines.depth import measure_coverage, samtools_depth_argv
from .engines.map import ensure_reference_index, iter_output, map_reads
from .engines.pileup import PileupSite, allele_fraction_at, pileup_at
from .records import (PLATFORM_FASTA, PLATFORM_ILLUMINA, PLATFORM_ONT,
                      STATUS_FAIL, STATUS_WARN, Check, CohortResult,
                      DatabaseVersion, Interpretation, LineageCall, QCMetrics,
                      SampleInput, SampleResult, SpeciesCall, Variant,
                      VARIANT_SNP, normalise_platform)
from .resistance import consensus as consensus_engine
from .resistance import ntm as ntm_engine
from .resistance import rules as rules_engine
from .resistance.catalogues import (Catalogue, calls_for_variant,
                                    database_versions as catalogue_versions,
                                    load_catalogues)
from .seqio import assembly_stats, validate_fasta, validate_fastq
from .typing.lineage import (H37RV_CHROM_ALIASES, BarcodeSite, barcode_path,
                             call_lineage, lineage_checks, lineage_not_applicable,
                             load_barcode, scheme_description)
from .typing.species import (ANI_FETCH_HINT, MarkerSnp, ReferenceGenome,
                             ani_matches, identify_species, load_mac_markers,
                             load_reference_set, species_checks)
from .utils import LOG, MjolnirError, ensure_dir, human_time, safe_name

#: Read handling is not a silent step. MTBseq performs no trimming, adapter
#: removal or read QC at any point (design §9b); Mjolnir does the same and
#: prints this instead of leaving the reader to assume either way.
NO_TRIMMING_NOTE = (
    "reads were analysed as submitted: Mjolnir performs no trimming, adapter "
    "removal or read filtering, so adapter content and quality decay are "
    "reflected in the coverage and allele-fraction numbers rather than removed "
    "from them"
)

#: Every allele fraction, depth and heterozygosity number in a Mjolnir report
#: comes out of one pileup, and that pileup applies MixInfect's base-quality
#: filter per base rather than per site (samtools mpileup -Q, see
#: engines/pileup.py). Stated wherever the heterozygosity numbers are.
PILEUP_QUALITY_NOTE = (
    "heterozygosity was measured from a pileup run with a minimum base quality "
    "of {0} (MixInfect's filter, applied per base before the counts existed), "
    "not from a per-site quality average"
).format(config.HET_MIN_QUAL)

#: The catalogues, and the tbdb barcode, key their coordinates to gene and HGVS
#: names as well as to positions. Mjolnir ships no variant annotator, so the
#: HGVS join key is only present when the caller or an upstream annotation
#: supplied it — and WHO's own protocol is coordinate-based, which is why this
#: is a stated limitation rather than a blocking one (§5.3).
NO_ANNOTATION_NOTE = (
    "no gene-level annotation is attached to the called variants, so "
    "cross-catalogue matching by HGVS name (MTBseq, tbdb) and the NTM "
    "rrl/rrs/erm(41) rules cannot fire; WHO's own coordinate-based matching is "
    "unaffected and is what produced any graded call below"
)

#: Positions Mjolnir piles up for reasons other than variant calling.
_PILEUP_PURPOSES = ("lineage barcode", "MAC marker panel", "catalogue coordinates")

#: Output formats ``write_outputs`` understands. ``tsv`` and ``json`` are one
#: writer — :func:`mjolnir.report.tables.write_sample_tables` emits both, and
#: splitting them would let a run produce a TSV whose JSON says something else.
OUTPUT_FORMATS: Tuple[str, ...] = ("tables", "html", "pdf")
FORMAT_ALIASES: Dict[str, str] = {
    "tsv": "tables", "json": "tables", "tables": "tables",
    "html": "html", "htm": "html", "pdf": "pdf",
}


def _stamp() -> str:
    """Local time with an offset, so two machines' outputs can be ordered."""
    return datetime.now().astimezone().isoformat(timespec="seconds")


# ---------------------------------------------------------------------------
# Run options
# ---------------------------------------------------------------------------

@dataclass
class RunOptions:
    """What to run and with which tool, separately from the thresholds.

    Kept apart from :class:`~mjolnir.config.Config` because these are choices
    about *how* a run is executed — which caller, whether to build a missing
    index — while ``Config`` holds the numbers a reader may need to see beside a
    result. A threshold override belongs in the report; ``--build-index`` does
    not.
    """

    #: Force a mapper or caller instead of choosing by platform and PATH.
    mapper: str = ""
    caller: str = ""
    clair3_model: str = ""
    #: Build a missing reference index instead of refusing. Off by default: the
    #: database directory may be shared and read-only, and a silent multi-minute
    #: index build in the middle of a batch looks exactly like a hang.
    build_index: bool = False
    #: None means "by platform" (marked on Illumina, not on ONT).
    mark_duplicates: Optional[bool] = None
    #: Allow the bcftools fallback on ONT. It disables indel calling entirely,
    #: so it has to be asked for and the report says it was used.
    allow_degraded_ont_calling: bool = False
    #: A Kraken2 report the operator already produced, instead of running it.
    kraken2_report: Optional[Path] = None
    #: Which intended uses the validity verdict is answered for.
    intended_use: Tuple[str, ...] = tuple(INTENDED_USES)
    #: Compute per-sample callable regions. Needed for cohort distances (they
    #: are the denominator) and costs a second pass of ``samtools depth``, so it
    #: is off unless a cohort is being built.
    callable_regions: bool = False
    #: Stages that can be turned off for a partial run.
    typing: bool = True
    resistance: bool = True
    contamination: bool = True
    interpret: bool = True

    def __post_init__(self) -> None:
        if self.kraken2_report is not None:
            self.kraken2_report = Path(self.kraken2_report)
        self.intended_use = tuple(self.intended_use)


# ---------------------------------------------------------------------------
# Small helpers shared by the stages
# ---------------------------------------------------------------------------

def read_fai(reference: Path) -> List[Tuple[str, int]]:
    """Contig names and lengths from ``<reference>.fai``.

    The names matter as much as the lengths: tbdb writes ``Chromosome`` where
    the NCBI FASTA says ``NC_000962.3``, and a pileup requested under the wrong
    name returns nothing at all rather than failing, which would read as a
    genome-wide coverage gap.
    """
    fai = Path(str(reference) + ".fai")
    if not fai.exists():
        # Built on demand when the directory allows it: a .fai on a 4.4 Mb genome
        # is a one-second job, and refusing it only sends the reader away to run
        # a command shorter than the error explaining it.
        ensure_reference_index(reference, tool="", build=False)
    contigs: List[Tuple[str, int]] = []
    with open(str(fai), "rt", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise MjolnirError(
                    "{0} line {1} is not a samtools .fai record".format(fai, line_number))
            contigs.append((fields[0], int(fields[1])))
    if not contigs:
        raise MjolnirError("{0} lists no contigs".format(fai))
    return contigs


def contig_translation(declared: Iterable[str],
                       reference_contigs: Sequence[str]) -> Dict[str, str]:
    """Map catalogue/barcode contig names onto the names this BAM uses.

    Only two translations are allowed, and both are identities rather than
    guesses: a name that is already present, and a known H37Rv alias when the
    reference carries exactly one of them. Anything else maps to itself, which
    means the pileup returns nothing for those positions and the callable count
    is zero — the truthful outcome for a barcode or a catalogue applied to the
    wrong reference, and the one that shows up in the report as a coverage gap
    rather than as a clean absence of variants.
    """
    present = set(reference_contigs)
    aliases = [name for name in reference_contigs if name in H37RV_CHROM_ALIASES]
    mapping: Dict[str, str] = {}
    for name in declared:
        if name in present:
            mapping[name] = name
        elif name in H37RV_CHROM_ALIASES and len(aliases) == 1:
            mapping[name] = aliases[0]
        elif len(present) == 1:
            mapping[name] = next(iter(present))
        else:
            mapping[name] = name
    return mapping


def callable_regions(bam: Path, *, platform: str, min_depth: int,
                     cfg: Optional[Config] = None, name: str = "") -> Regions:
    """Intervals at or above *min_depth*, as the cohort's distance denominator.

    A pairwise SNP distance without this is a numerator with no denominator, and
    the design is explicit that 12 differences over 4.1 Mb and 12 over 400 kb are
    not the same statement (§9). It is a second pass of ``samtools depth`` rather
    than a by-product of :func:`~mjolnir.engines.depth.measure_coverage`, which
    keeps a histogram and not a position list — so it is only computed when a
    cohort is actually being built.
    """
    argv = samtools_depth_argv(
        bam, platform=normalise_platform(platform),
        mtbseq_compat=bool(cfg.mtbseq_compat) if cfg is not None else False)
    regions = Regions(name=name or "callable at >= {0}x".format(min_depth))
    chrom = ""
    start: Optional[int] = None
    previous = 0
    for line in iter_output(argv):
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        this_chrom = fields[0]
        try:
            pos = int(fields[1])
            depth = int(fields[2])
        except ValueError:
            raise MjolnirError(
                "unexpected 'samtools depth' output while measuring callable "
                "regions of {0}: {1!r}".format(bam, line[:120]))
        if depth >= min_depth:
            if start is not None and this_chrom == chrom and pos == previous + 1:
                previous = pos
                continue
            if start is not None:
                regions.add_1based(chrom, start, previous)
            chrom, start, previous = this_chrom, pos, pos
        elif start is not None:
            regions.add_1based(chrom, start, previous)
            start = None
    if start is not None:
        regions.add_1based(chrom, start, previous)
    return regions


def _dedupe(items: Iterable[str]) -> List[str]:
    """Order-preserving de-duplication. The same caveat twice reads as noise,
    and a reader who skims a repeated sentence stops reading the list."""
    seen: Set[str] = set()
    out: List[str] = []
    for item in items:
        if not item or item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out


def _observation_from_site(site: PileupSite) -> SiteObservation:
    """A heterozygosity observation from one pileup column.

    ``qual`` is the pileup's own minimum base quality rather than a measured
    site average: every base in ``site.counts`` already cleared it (the pileup
    runs with ``-Q``), so recording it here applies MixInfect's Q>=20 filter
    exactly once, per base, where it is strictest. :data:`PILEUP_QUALITY_NOTE`
    travels with the resulting numbers so the reader knows which of the two it
    is looking at.
    """
    return SiteObservation(
        pos=site.pos, counts=dict(site.counts), qual=float(config.HET_MIN_QUAL),
        chrom=site.chrom, ref_base=site.ref_base, raw_depth=site.raw_depth)


def _snp_observation(variant: Variant) -> Optional[SiteObservation]:
    """A heterozygosity observation for one called SNP, or None if it cannot be one.

    Built from the caller's own allele depths. Without them there is no minor
    allele fraction to measure and the site is left out of the denominator
    rather than being counted as homozygous, which would push the mixture
    statistic downwards for exactly the samples whose evidence is thinnest.
    """
    if variant.variant_type != VARIANT_SNP or len(variant.ref) != 1 or len(variant.alt) != 1:
        return None
    if variant.alt_reads is None or variant.ref_reads is None:
        return None
    counts = {variant.ref.upper(): int(variant.ref_reads),
              variant.alt.upper(): int(variant.alt_reads)}
    return SiteObservation(pos=variant.pos, counts=counts, qual=variant.qual,
                           chrom=variant.chrom, ref_base=variant.ref.upper(),
                           raw_depth=variant.depth)


def _lineage_observations(sites: Sequence[BarcodeSite],
                          pileup: Mapping[Tuple[str, int], PileupSite],
                          translation: Mapping[str, str]) -> List[LineageSite]:
    """F2/F47 input: the barcode positions, with which allele is derived.

    The derived allele is the one ``barcode.bed`` names for the taxon; the
    ancestral allele is the reference base as the pileup saw it. A site the
    pileup did not cover is skipped here rather than being passed on with zero
    counts, because :mod:`mjolnir.contamination.heterozygosity` would then treat
    it as callable-and-clean.
    """
    out: List[LineageSite] = []
    for site in sites:
        key = (translation.get(site.chrom, site.chrom), site.pos)
        column = pileup.get(key)
        if column is None or not column.covered:
            continue
        alleles = site.alleles
        derived = alleles[0] if alleles else ""
        out.append(LineageSite(
            pos=site.pos, counts=dict(column.counts), qual=float(config.HET_MIN_QUAL),
            chrom=column.chrom, ref_base=column.ref_base,
            raw_depth=column.raw_depth, lineage=site.taxon,
            derived_allele=derived, ancestral_allele=column.ref_base))
    return out


# ---------------------------------------------------------------------------
# The pipeline
# ---------------------------------------------------------------------------

class Pipeline(object):
    """Runs the stages for one sample, and joins a cohort from their results.

    One instance per run. The catalogues, the barcode scheme, the ANI reference
    set and the model client are loaded once and shared: a 159-sample cohort
    parsing the WHO xlsx 159 times would spend most of its wall clock in
    openpyxl, and — more importantly — two samples in one report must be graded
    against byte-identical catalogue files or their calls are not comparable.

    A load failure is cached too. If the WHO catalogue is absent, every sample
    records the same reason and none of them re-raises it, so the operator gets
    one explanation and a complete set of results that say what is missing.
    """

    def __init__(self, cfg: Optional[Config] = None,
                 options: Optional[RunOptions] = None) -> None:
        self.config = cfg if cfg is not None else Config()
        self.options = options if options is not None else RunOptions()
        self._catalogues_cache: Optional[Dict[str, Catalogue]] = None
        self._catalogues_error = ""
        self._barcode_cache: Optional[List[BarcodeSite]] = None
        self._barcode_error = ""
        self._barcode_scheme = ""
        self._references_cache: Optional[List[ReferenceGenome]] = None
        self._references_error = ""
        self._markers_cache: Optional[List[MarkerSnp]] = None
        self._client_resolved = False
        self._client: Any = None
        #: Tools that actually ran, for ``SampleResult.tool_versions``.
        self.tools_used: Set[str] = set()
        #: Databases actually consulted, for the version annex.
        self.databases_used: Set[str] = set()
        #: Notes about how the mapping reference was chosen, raised into
        #: the sample's caveats. A reference picked by ANI proximity rather
        #: than by a species call changes what every coordinate means, so
        #: it has to reach the report rather than only the log.
        self.reference_notes: List[str] = []
        #: Filled by :meth:`run_sample` when ``options.callable_regions`` is on.
        self.sample_variants: Dict[str, SampleVariants] = {}

    # -- shared, lazily loaded inputs ---------------------------------------

    def catalogues(self) -> Dict[str, Catalogue]:
        """The MTBC catalogues, loaded once; ``{}`` when they could not be.

        The paths are resolved from the database registry first and from
        :mod:`mjolnir.resistance.catalogues`' own defaults second, because the
        registry is what writes the files and the two modules were written to
        different layouts. Passing explicit paths is what keeps a fetched
        database findable without either module having to guess.
        """
        if self._catalogues_cache is not None:
            return self._catalogues_cache
        paths = self.catalogue_paths()
        try:
            loaded = load_catalogues(
                db_dir=self.config.db_dir,
                who=paths.get("who"), who_coordinates=paths.get("who_coordinates"),
                mtbseq=paths.get("mtbseq"), tbdb=paths.get("tbdb"),
                require_who=True)
        except MjolnirError as exc:
            self._catalogues_error = str(exc)
            LOG.warning("resistance catalogues unavailable: %s", exc)
            self._catalogues_cache = {}
            return self._catalogues_cache
        for name in ("who-catalogue-v2", "tbdb", "mtbseq-resistance"):
            self.databases_used.add(name)
        self._catalogues_cache = loaded
        return loaded

    def catalogue_paths(self) -> Dict[str, Optional[Path]]:
        """Where each catalogue file actually is, registry layout first."""
        root = Path(self.config.db_dir)
        wanted = (
            ("who", "who-catalogue-v2", "WHO-UCN-TB-2023.7-eng.xlsx",
             "who/WHO-UCN-TB-2023.7-eng.xlsx"),
            ("who_coordinates", "who-catalogue-v2", "Genomic_coordinates_7May2024.vcf.gz",
             "who/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"),
            ("mtbseq", "mtbseq-resistance", "MTB_Resistance_Mediating.txt",
             "mtbseq/MTB_Resistance_Mediating.txt"),
            ("tbdb", "tbdb", "mutations.csv", "tbdb/mutations.csv"),
        )
        found: Dict[str, Optional[Path]] = {}
        for key, database, filename, legacy in wanted:
            candidates: List[Path] = []
            try:
                candidates.append(db_registry.spec_for(database).path_to(root, filename))
            except MjolnirError:
                pass
            candidates.append(root / legacy)
            # When nothing is installed, the registry path is still the right
            # answer to hand on: the error a loader raises then names the place
            # ``mjolnir db fetch`` would have written the file, rather than a
            # legacy layout nothing writes to.
            found[key] = next((p for p in candidates if p.exists()),
                              candidates[0] if candidates else None)
        return found

    def barcode(self) -> List[BarcodeSite]:
        """tbdb's lineage barcode, loaded once; ``[]`` when it is not installed."""
        if self._barcode_cache is not None:
            return self._barcode_cache
        try:
            path = barcode_path(self.config.db_dir)
            sites = load_barcode(path)
        except MjolnirError as exc:
            self._barcode_error = str(exc)
            LOG.warning("MTBC lineage barcode unavailable: %s", exc)
            self._barcode_cache = []
            return self._barcode_cache
        self._barcode_scheme = scheme_description(sites, name="tbdb barcode.bed")
        self.databases_used.add("tbdb")
        self._barcode_cache = sites
        return sites

    def references(self) -> List[ReferenceGenome]:
        """The ANI reference set, loaded once; ``[]`` when it is not installed."""
        if self._references_cache is not None:
            return self._references_cache
        try:
            found = load_reference_set(self.config.db_dir)
        except MjolnirError as exc:
            self._references_error = str(exc)
            LOG.warning("ANI reference set unavailable: %s", exc)
            self._references_cache = []
            return self._references_cache
        self._references_cache = found
        return found

    def markers(self) -> List[MarkerSnp]:
        """The MAC marker panel. Absent is a supported state, not an error."""
        if self._markers_cache is None:
            self._markers_cache = load_mac_markers(self.config.db_dir, required=False)
        return self._markers_cache

    def client(self) -> Any:
        """The model client, or None when the operator turned the model off."""
        if not self._client_resolved:
            self._client_resolved = True
            if self.options.interpret:
                from .agent import client_from_config

                self._client = client_from_config(self.config)
            else:
                self._client = None
        return self._client

    # -- stage plumbing -----------------------------------------------------

    def _stage(self, result: SampleResult, name: str, category: str,
               consequence: str, func: Any, *args: Any, **kwargs: Any) -> Any:
        """Run one optional stage; on failure record what was lost and continue.

        The recorded :class:`~mjolnir.records.Check` is ``measured=False``, which
        every consumer renders as "not measured". That is the whole mechanism
        behind the fifth house rule: there is no path from "this stage did not
        run" to a passing check, so a report cannot end up cleaner because less
        of it worked.
        """
        try:
            return func(*args, **kwargs)
        except MjolnirError as exc:
            LOG.warning("%s: %s not available: %s", result.sample_id, name, exc)
            result.warnings.append("{0}: {1}".format(name, exc))
            result.checks.append(Check.not_measured(
                name, "{0}. Consequence: {1}".format(exc, consequence),
                category=category))
            return None

    def _record_tools(self, *tools: str) -> None:
        for tool in tools:
            if tool:
                self.tools_used.add(tool)

    # -- reference ----------------------------------------------------------

    def resolve_reference(self, sample: SampleInput, species: Optional[SpeciesCall]) -> Path:
        """Pick the reference this sample is called against, or refuse.

        The order is evidence first, convenience never: an explicit ``--ref``, a
        per-sample reference from the manifest, H37Rv for anything in the MTBC,
        the ANI set's own genome for the species that was actually identified,
        and MTBseq's NTM references as the last resort for the three species it
        ships. If none of those apply the run stops — mapping an *M. chimaera*
        isolate to H37Rv "so that something runs" is the failure this project
        was written against, and it is what the 2022 local run did.
        """
        if sample.reference is not None:
            return Path(sample.reference)
        if self.config.reference is not None:
            return Path(self.config.reference)

        if species is not None and species.complex == config.COMPLEX_MTBC:
            try:
                path = db_fetch.database_file(
                    "h37rv", "{0}.fasta".format(config.H37RV_ACCESSION),
                    db_root=self.config.db_dir)
            except MjolnirError as exc:
                raise MjolnirError(
                    "sample {0!r} is in the M. tuberculosis complex and needs the "
                    "H37Rv reference every MTBC catalogue is written against.\n"
                    "  {1}".format(sample.sample_id, exc))
            self.databases_used.add("h37rv")
            return path

        if species is not None and species.resolved_to_species:
            for reference in self.references():
                if reference.name.strip().lower() == species.name.strip().lower() \
                        and reference.exists:
                    return Path(reference.path)
            fallback = self._mtbseq_ntm_reference(species.name)
            if fallback is not None:
                return fallback

        # Not resolved to a species, but ANI still named the genome this isolate
        # is nearest to. Mapping to that is not the failure this method guards
        # against: refusing to name a MAC isolate below the complex is the honest
        # answer (design §6), and the nearest genome is the right thing to call
        # it against. Without this branch every M. chimaera isolate stops here,
        # because "cannot resolve below complex" is the correct species result
        # and it left no reference behind.
        nearest = self._nearest_reference(species)
        if nearest is not None:
            path, name, ani = nearest
            LOG.info("%s: species not resolved below %s; calling against the "
                     "nearest reference by ANI, %s (%.2f%%)",
                     sample.sample_id, species.complex or "genus", name, ani)
            self.reference_notes.append(
                "the species could not be resolved below {0}, so this sample was "
                "called against the nearest genome in the ANI reference set, {1} "
                "at {2:.2f}% — variant positions are relative to that genome and "
                "not to a genome of the sample's own species".format(
                    species.complex or "genus", name, ani))
            return path

        raise MjolnirError(
            "no reference could be chosen for sample {0!r}: the species call is "
            "{1!r} and no genome for it is installed.\n"
            "  pass one with --ref, or install the ANI reference set: {2}".format(
                sample.sample_id,
                species.display if species is not None else "not determined",
                ANI_FETCH_HINT)
        )

    def _nearest_reference(self, species: Optional[SpeciesCall]
                           ) -> Optional[Tuple[Path, str, float]]:
        """The installed genome the sample is closest to, by the ANI candidates.

        Returns ``(path, name, ani)``, or None when there are no candidates or
        none of them clears the genus floor — a query that matches nothing in the
        set is not a mycobacterium the set can speak for, and guessing a
        reference for it would put every downstream coordinate in the wrong
        genome without saying so.
        """
        if species is None or not species.candidates:
            return None
        by_name = {}
        for reference in self.references():
            if reference.exists:
                by_name.setdefault(reference.name.strip().lower(), reference)
        for candidate in species.candidates:
            ani = candidate.get("ani_percent")
            name = str(candidate.get("name") or "").strip()
            if ani is None or not name:
                continue
            if float(ani) < config.ANI_GENUS_FLOOR:
                break
            reference = by_name.get(name.lower())
            if reference is not None:
                return Path(reference.path), reference.name, float(ani)
        return None

    def _mtbseq_ntm_reference(self, species_name: str) -> Optional[Path]:
        """One of the three NTM genomes MTBseq ships, when it matches the call."""
        known = {
            "mycobacteroides abscessus": "M._abscessus_CIP-104536T_2014-02-03.fasta",
            "mycobacterium abscessus": "M._abscessus_CIP-104536T_2014-02-03.fasta",
            "mycobacterium chimaera": "M._chimaera_DSM44623_2016-01-28.fasta",
            "mycobacterium fortuitum": "M._fortuitum_CT6_2016-01-08.fasta",
        }
        filename = known.get(" ".join(str(species_name or "").split()).lower())
        if not filename:
            return None
        try:
            path = db_fetch.database_file("mtbseq-ntm-references", filename,
                                          db_root=self.config.db_dir)
        except MjolnirError:
            return None
        self.databases_used.add("mtbseq-ntm-references")
        return path

    def reference_gc(self, reference: Path) -> Optional[float]:
        """GC of the reference itself, so the composition check has a baseline.

        H37Rv's value is a registered threshold; for any other genome the only
        honest baseline is the file on disk, so it is measured. ``None`` when it
        cannot be read, which makes the GC check "not measured" rather than a
        comparison against the wrong organism.
        """
        name = Path(reference).name
        if name.startswith(config.H37RV_ACCESSION):
            return config.H37RV_GC
        try:
            stats = assembly_stats(reference)
        except MjolnirError as exc:
            LOG.debug("reference GC not measured for %s: %s", reference, exc)
            return None
        return stats.gc_fraction

    # -- the sample ---------------------------------------------------------

    def run_sample(self, sample: SampleInput, workdir: Path,
                   evidence_checks: Sequence[Check] = ()) -> SampleResult:
        """Everything Mjolnir can determine about one sample, honestly bounded.

        Never raises for a missing capability — those become recorded absences.
        It does raise for the two things that make the rest meaningless: input
        files that cannot be read, and a reference that cannot be chosen.
        """
        started = time.time()
        platform = normalise_platform(sample.platform)
        result = SampleResult(
            sample_id=sample.sample_id,
            platform=platform,
            inputs=[str(p) for p in sample.paths],
            mjolnir_version=__version__,
            profile=self.config.profile,
            started_at=_stamp(),
        )
        result.checks.extend(evidence_checks)
        result.caveats.extend(config.platform_caveats(platform))
        result.caveats.extend(self.reference_notes)
        if platform != PLATFORM_FASTA:
            result.caveats.append(NO_TRIMMING_NOTE)
        if platform == PLATFORM_ONT:
            result.caveats.append(
                "ONT data is interpreted on the assumption of {0}; Mjolnir cannot "
                "verify the chemistry or basecaller from the reads".format(
                    config.ONT_MINIMUM_CONFIGURATION))
        if self.config.mtbseq_compat:
            result.caveats.append(
                "--mtbseq-compat: allele fractions use MTBseq's denominator "
                "(N and GAP counted, GAP winning ties) and no mapping-quality "
                "floor is applied, so these numbers are comparable to an MTBseq "
                "run and not to Mjolnir's own defaults")
        overrides = self.config.overridden_thresholds()
        if overrides:
            result.caveats.append(
                "this run overrode published defaults: {0}".format(
                    ", ".join("{0}={1}".format(k, v) for k, v in sorted(overrides.items()))))

        self._validate_inputs(sample)

        # -- species, first pass: it decides the reference -------------------
        query = sample.paths[0]
        matches = self._stage(
            result, "ani_reference_set", "typing",
            "no species can be named, and no reference can be chosen from one",
            self._ani, query, platform)
        provisional = None
        if matches is not None:
            provisional = self._stage(
                result, "species_identification", "typing",
                "the sample is reported as an unidentified mycobacterium",
                identify_species, query, self.config.db_dir,
                is_reads=platform != PLATFORM_FASTA,
                threads=self.config.threads, matches=matches,
                markers=self.markers(), marker_counts={})
        identified = provisional is not None
        if provisional is None:
            # No method is claimed on a call that was never made. Naming one
            # would make ``species_checks`` report the identification as having
            # come from a method other than ANI, which is the one thing the
            # design forbids saying about a species call.
            provisional = SpeciesCall(
                name="unresolved", method="",
                caveats=[self._references_error or
                         "species identification did not run",
                         config.SPECIES_METHOD_REFUSAL])
        result.species = provisional

        # Cleared per sample: a note about how one sample's reference was
        # chosen must not appear in the next sample's caveats.
        self.reference_notes = []
        reference = self.resolve_reference(sample, provisional)
        result.reference = str(reference)
        contigs = read_fai(reference)
        contig_names = [name for name, _length in contigs]
        reference_length = sum(length for _name, length in contigs)
        ref_gc = self.reference_gc(reference)

        # -- mapping and coverage -------------------------------------------
        alignment = None
        bam: Optional[Path] = None
        if platform == PLATFORM_FASTA:
            result.checks.append(Check.not_measured(
                "coverage", config.FASTA_CAPABILITY_LOSS,
                source=config.source_for("fasta_capability_loss"), category="qc"))
        else:
            alignment = self._stage(
                result, "read_mapping", "qc",
                "no coverage, no variants and no pileup for this sample",
                map_reads, sample, reference,
                Path(workdir) / "{0}.bam".format(safe_name(sample.sample_id)),
                threads=self.config.threads,
                mapper=self.options.mapper or None,
                mark_duplicates=self.options.mark_duplicates,
                build_index=self.options.build_index,
                config=self.config)
            if alignment is not None:
                bam = alignment.bam
                self._record_tools(alignment.mapper, "samtools")
                qc = self._stage(
                    result, "coverage", "qc",
                    "depth, breadth and evenness are unknown for this sample",
                    measure_coverage, bam, reference, platform=platform,
                    config=self.config, reference_length=reference_length,
                    reference_name=Path(reference).name, reference_gc=ref_gc,
                    duplicates_marked=alignment.duplicates_marked,
                    duplicates_note=alignment.duplicates_note,
                    threads=self.config.threads)
                if qc is not None:
                    result.qc = qc

        # -- variant calling -------------------------------------------------
        variants: List[Variant] = []
        if platform == PLATFORM_FASTA or bam is not None:
            called = self._stage(
                result, "variant_calling", "resistance",
                "no variants were called, so no resistance determinant could be "
                "detected — which is not the same as none being present",
                call_variants, sample.sample_id, bam or Path("."), reference,
                platform, Path(workdir) / "variants",
                config=self.config, caller=self.options.caller or None,
                threads=self.config.threads,
                clair3_model=self.options.clair3_model or None,
                allow_degraded_fallback=self.options.allow_degraded_ont_calling,
                assembly=sample.assembly)
            if called is not None:
                variants = list(called.variants)
                result.variants = variants
                result.caveats.extend(called.caveats)
                self._record_tools(_caller_binary(called.caller))
                if called.degraded:
                    result.warnings.append(
                        "variants were called by the degraded ONT path ({0})".format(
                            called.caller))
        if variants and not any(v.gene for v in variants):
            result.caveats.append(NO_ANNOTATION_NOTE)
            result.checks.append(Check.not_measured(
                "variant_gene_annotation", NO_ANNOTATION_NOTE, category="resistance"))

        # -- one pileup for the barcode, the markers and the catalogue -------
        barcode_sites = self.barcode() if self.options.typing else []
        marker_snps = self.markers() if self.options.typing else []
        is_mtbc, assumption = self.mtbc_context(result.species, reference, contig_names)
        if assumption:
            result.caveats.append(assumption)
            result.checks.append(Check.not_measured(
                "mtbc_membership", assumption,
                source=config.source_for("species_method_refusal"), category="typing"))
        catalogue_positions: Dict[str, Set[Tuple[str, int]]] = {}
        if self.options.resistance and is_mtbc:
            # Catalogue coordinates are H37Rv coordinates. Piling them up against
            # an NTM reference would read 30,000 unrelated positions and print
            # allele fractions for variants that cannot exist there.
            catalogue_positions = self._catalogue_positions()
        pileup, translation = self._pileup(
            result, bam, reference, contig_names, platform,
            barcode_sites, marker_snps, catalogue_positions)

        # -- species, second pass, with the MAC markers ----------------------
        if matches is not None and marker_snps and pileup:
            final = self._stage(
                result, "species_markers", "typing",
                "MAC is reported at complex level",
                identify_species, query, self.config.db_dir,
                is_reads=platform != PLATFORM_FASTA, threads=self.config.threads,
                matches=matches, markers=marker_snps,
                marker_counts=self._counts_for(marker_snps, pileup, translation))
            if final is not None:
                result.species = final
        if identified:
            result.checks.extend(species_checks(result.species))
        else:
            result.checks.append(Check.not_measured(
                "species_identification",
                "no species identification was attempted or completed, so this "
                "sample is reported as an unidentified mycobacterium. {0}".format(
                    config.SPECIES_METHOD_REFUSAL),
                source=config.source_for("species_method_refusal"),
                category="typing"))

        # -- lineage ----------------------------------------------------------
        result.lineage = self._lineage(result, barcode_sites, pileup, translation,
                                       platform, is_mtbc)
        result.checks.extend(lineage_checks(result.lineage))
        if result.lineage.is_bcg:
            result.caveats.append(config.BCG_PZA_NOTE)

        # -- allele fractions at catalogue positions -------------------------
        self._backfill_allele_fractions(variants, pileup, translation)

        # -- resistance -------------------------------------------------------
        if self.options.resistance:
            self._resistance(result, variants, pileup, translation, catalogue_positions,
                             platform, is_mtbc)

        # -- contamination ----------------------------------------------------
        if self.options.contamination:
            self._contamination(result, sample, barcode_sites, pileup, translation,
                                platform, ref_gc)

        # -- cohort input -----------------------------------------------------
        if self.options.callable_regions:
            self._collect_cohort_input(result, sample, bam, reference, platform, variants)

        # -- provenance and prose ---------------------------------------------
        result.caveats = _dedupe(result.caveats)
        result.warnings = _dedupe(result.warnings)
        result.tool_versions = self._tool_versions()
        result.database_versions = self._database_versions()
        result.status = result.overall_status()
        result.finished_at = _stamp()
        result.runtime_seconds = time.time() - started

        if self.options.interpret:
            result.interpretation = self._interpret(result)
        else:
            result.interpretation = Interpretation(
                headline="", body="", rule_only=True,
                discarded_reason="interpretation was switched off for this run")

        LOG.info("%s: %s | %s | %d variants | %s | %s",
                 result.sample_id, result.species.display, result.lineage.display,
                 len(result.variants), result.contamination.verdict,
                 human_time(result.runtime_seconds))
        return result

    # -- stages -------------------------------------------------------------

    def _validate_inputs(self, sample: SampleInput) -> None:
        """Refuse a sample whose files cannot be read, before anything runs."""
        for path in sample.paths:
            if sample.platform == PLATFORM_FASTA:
                validate_fasta(path)
            else:
                validate_fastq(path)

    def _ani(self, query: Path, platform: str) -> List[Any]:
        """ANI matches against the reference set, used twice and run once."""
        references = self.references()
        if not references:
            raise MjolnirError(
                self._references_error or
                "no ANI reference set is installed, so no species can be named")
        found = ani_matches(query, references, is_reads=platform != PLATFORM_FASTA,
                            threads=self.config.threads, db_dir=self.config.db_dir)
        self._record_tools(found[0].method if found else "")
        return found

    def _pileup(self, result: SampleResult, bam: Optional[Path], reference: Path,
                contig_names: Sequence[str], platform: str,
                barcode_sites: Sequence[BarcodeSite],
                marker_snps: Sequence[MarkerSnp],
                catalogue_positions: Mapping[str, Set[Tuple[str, int]]]
                ) -> Tuple[Dict[Tuple[str, int], PileupSite], Dict[str, str]]:
        """One pileup over the union of every position anything downstream needs.

        Returns the sites keyed by the *reference's* contig names together with
        the translation that got there, so a caller holding barcode or catalogue
        coordinates can find its own rows without either side guessing at contig
        naming.
        """
        declared: Set[str] = set()
        wanted: Set[Tuple[str, int]] = set()
        for site in barcode_sites:
            declared.add(site.chrom)
            wanted.add((site.chrom, site.pos))
        for marker in marker_snps:
            declared.add(marker.chrom)
            wanted.add((marker.chrom, marker.pos))
        for positions in catalogue_positions.values():
            for chrom, pos in positions:
                declared.add(chrom)
                wanted.add((chrom, pos))
        translation = contig_translation(declared, contig_names)
        if bam is None or not wanted:
            if wanted and bam is None and platform != PLATFORM_FASTA:
                result.checks.append(Check.not_measured(
                    "direct_pileup",
                    "no alignment was produced, so the {0} could not be read from "
                    "a pileup".format(", ".join(_PILEUP_PURPOSES)),
                    category="typing"))
            return {}, translation

        translated = sorted(set((translation.get(chrom, chrom), pos)
                                for chrom, pos in wanted))
        sites = self._stage(
            result, "direct_pileup", "typing",
            "the lineage barcode, the MAC markers and the allele fractions at "
            "catalogue positions are all unavailable",
            pileup_at, bam, reference, translated, platform=platform,
            config=self.config, scratch_dir=self.config.tmp_dir)
        if sites is None:
            return {}, translation
        self._record_tools("samtools")
        LOG.debug("%s: piled up %d positions (%s)", result.sample_id, len(sites),
                  ", ".join(_PILEUP_PURPOSES))
        return sites, translation

    def _counts_for(self, markers: Sequence[MarkerSnp],
                    pileup: Mapping[Tuple[str, int], PileupSite],
                    translation: Mapping[str, str]) -> Dict[Tuple[str, int], Dict[str, int]]:
        """Marker-panel counts, keyed as the marker file spells its contig."""
        counts: Dict[Tuple[str, int], Dict[str, int]] = {}
        for marker in markers:
            site = pileup.get((translation.get(marker.chrom, marker.chrom), marker.pos))
            if site is not None:
                counts[(marker.chrom, marker.pos)] = dict(site.counts)
        return counts

    def mtbc_context(self, species: SpeciesCall, reference: Path,
                     contigs: Sequence[str] = ()) -> Tuple[bool, str]:
        """Whether the MTBC catalogues and barcode apply, and on what basis.

        Measured membership is the first answer. The second exists because
        Mjolnir ships no ANI reference set (design §14 leaves its composition
        open), so on a fresh installation no species can be *measured* at all —
        and a tool that then refuses to grade a *M. tuberculosis* genome mapped
        to H37Rv is useless in its own default configuration.

        So a run explicitly pointed at H37Rv is treated as an operator assertion
        that this is an MTBC isolate, and the second return value is the sentence
        that says exactly that. It is recorded as a caveat and as an unmeasured
        check, never as a species call: the report still prints "unresolved" for
        the species, because nothing identified it.
        """
        if species.complex == config.COMPLEX_MTBC or species.is_mtbc:
            return True, ""
        if species.resolved_to_species:
            return False, ""
        name = Path(reference).name
        # The accession is the real signal: every MTBC catalogue coordinate,
        # barcode site and mask interval in this tool is written against
        # NC_000962.3, so a reference carrying that accession is the coordinate
        # system they need. The bare alias "Chromosome" is deliberately not
        # accepted — any single-contig assembly can be called that.
        h37rv_contig = any(str(c).startswith("NC_000962") or str(c).startswith("AL123456")
                           for c in contigs)
        if h37rv_contig or name.startswith(config.H37RV_ACCESSION) \
                or "h37rv" in name.lower():
            return True, (
                "MTBC membership was not established from sequence: no species "
                "identification was made, and the catalogues and lineage barcode "
                "were applied because this run was called against the H37Rv "
                "coordinate system ({0}). That is an operator assertion, not a "
                "measurement.".format(name))
        return False, ""

    def _lineage(self, result: SampleResult, barcode_sites: Sequence[BarcodeSite],
                 pileup: Mapping[Tuple[str, int], PileupSite],
                 translation: Mapping[str, str], platform: str,
                 is_mtbc: bool) -> LineageCall:
        """The MTBC member from the SNP barcode, or a stated non-answer.

        Three distinct outcomes, and none of them is silence: not an MTBC
        isolate, no barcode installed, or no pileup to read it out of.
        """
        species = result.species
        if not is_mtbc:
            return lineage_not_applicable(
                species.display,
                "the lineage barcode defines MTBC lineages only" if species.complex
                else "no MTBC membership was established, so the barcode was not applied")
        if not barcode_sites:
            return lineage_not_applicable(
                species.display,
                self._barcode_error or "the tbdb lineage barcode is not installed")
        if not pileup:
            return lineage_not_applicable(
                species.display,
                "no pileup was available at the barcode positions, so no lineage "
                "was called; this is a coverage gap, not lineage 4")

        counts: Dict[Tuple[str, int], Dict[str, int]] = {}
        for site in barcode_sites:
            column = pileup.get((translation.get(site.chrom, site.chrom), site.pos))
            if column is not None:
                counts[(column.chrom, column.pos)] = dict(column.counts)
        called = self._stage(
            result, "lineage", "typing",
            "no lineage, sublineage, BCG or animal-lineage call",
            call_lineage, barcode_sites, counts, platform,
            scheme=self._barcode_scheme or "tbdb barcode.bed")
        if called is None:
            return lineage_not_applicable(species.display,
                                          "the lineage stage did not complete")
        return called

    def _backfill_allele_fractions(self, variants: Sequence[Variant],
                                   pileup: Mapping[Tuple[str, int], PileupSite],
                                   translation: Mapping[str, str]) -> None:
        """Fill in an allele fraction from the pileup where the caller gave none.

        Only where it is missing, and never on FASTA input: an assembly has no
        allele fractions at all, and writing 1.0 there would turn a capability
        loss into a confident number (§7).
        """
        if not pileup:
            return
        for variant in variants:
            if variant.allele_fraction is not None:
                continue
            site = pileup.get((translation.get(variant.chrom, variant.chrom), variant.pos))
            if site is None or not site.covered:
                continue
            fraction = allele_fraction_at(site, variant.ref, variant.alt)
            if fraction is None:
                continue
            variant.allele_fraction = fraction
            variant.is_major = config.is_major_variant(fraction)
            if variant.depth is None:
                variant.depth = site.acgt_depth
            variant.note = (variant.note + "; " if variant.note else "") + \
                "allele fraction measured from the direct pileup"

    # -- resistance ---------------------------------------------------------

    def _catalogue_positions(self) -> Dict[str, Set[Tuple[str, int]]]:
        """Catalogued coordinates per drug, from the loaded catalogues.

        Used for two things: the pileup that gives an allele fraction at every
        catalogued position whether or not the caller emitted a record there,
        and the per-drug coverage that decides whether a drug was evaluable at
        all. A drug whose positions were never callable is "not evaluable",
        which is a different statement from "no determinant detected".
        """
        positions: Dict[str, Set[Tuple[str, int]]] = {}
        for catalogue in self.catalogues().values():
            for entry in catalogue.entries:
                if not entry.coordinates:
                    continue
                bucket = positions.setdefault(config.normalise_drug(entry.drug), set())
                for chrom, pos, _ref, _alt in entry.coordinates:
                    bucket.add((chrom, int(pos)))
        return positions

    def _target_coverage(self, catalogue_positions: Mapping[str, Set[Tuple[str, int]]],
                         pileup: Mapping[Tuple[str, int], PileupSite],
                         translation: Mapping[str, str], platform: str
                         ) -> Dict[str, Optional[bool]]:
        """Whether each drug's catalogued positions were callable.

        ``None`` means it could not be established — on FASTA input, or with no
        pileup — and stays ``None`` rather than becoming ``True``. The threshold
        is :data:`mjolnir.config.MIN_BREADTH`, the same fraction the genome-wide
        breadth check applies, read here over one drug's target positions rather
        than over the whole reference; no second number is introduced for it.
        """
        if platform == PLATFORM_FASTA or not pileup or not catalogue_positions:
            return dict((drug, None) for drug in catalogue_positions)
        floor = self.config.degraded_depth_floor
        covered: Dict[str, Optional[bool]] = {}
        for drug, positions in catalogue_positions.items():
            if not positions:
                covered[drug] = None
                continue
            callable_count = 0
            for chrom, pos in positions:
                site = pileup.get((translation.get(chrom, chrom), pos))
                if site is not None and site.covered and site.acgt_depth >= floor:
                    callable_count += 1
            covered[drug] = (float(callable_count) / len(positions)) >= self.config.min_breadth
        return covered

    def _gene_drugs(self, catalogues: Mapping[str, Catalogue]) -> Dict[str, Tuple[str, ...]]:
        """Which drugs each gene is catalogued for.

        The silent-variant rule (§5.4) needs this and Mjolnir will not guess it:
        the map is derived from the catalogue that was actually loaded, so a
        catalogue edition that adds a gene brings its drugs with it.
        """
        mapping: Dict[str, Set[str]] = {}
        for catalogue in catalogues.values():
            for entry in catalogue.entries:
                if entry.gene:
                    mapping.setdefault(entry.gene.strip().lower(), set()).add(
                        config.normalise_drug(entry.drug))
        return dict((gene, tuple(sorted(drugs))) for gene, drugs in mapping.items())

    def _resistance(self, result: SampleResult, variants: List[Variant],
                    pileup: Mapping[Tuple[str, int], PileupSite],
                    translation: Mapping[str, str],
                    catalogue_positions: Mapping[str, Set[Tuple[str, int]]],
                    platform: str, is_mtbc: bool) -> None:
        """Drug calls, by whichever route the organism has an evidence base for."""
        if is_mtbc:
            self._resistance_mtbc(result, variants, pileup, translation,
                                  catalogue_positions, platform)
            return
        self._resistance_ntm(result, variants, platform)

    def _resistance_mtbc(self, result: SampleResult, variants: List[Variant],
                         pileup: Mapping[Tuple[str, int], PileupSite],
                         translation: Mapping[str, str],
                         catalogue_positions: Mapping[str, Set[Tuple[str, int]]],
                         platform: str) -> None:
        """The three-catalogue consensus (§5.5), or a stated absence of it."""
        catalogues = self.catalogues()
        if not catalogues:
            result.checks.append(Check.not_measured(
                "resistance_calling",
                "{0} No resistance call was made for any drug; this is an absence "
                "of evidence and not susceptibility.".format(
                    self._catalogues_error or "no catalogue was loaded."),
                source=config.source_for("anchor_catalogue"), category="resistance"))
            result.caveats.append(
                "no drug-resistance call was made: the catalogues were not "
                "available in this run")
            return

        who = catalogues.get(config.CATALOGUE_WHO)
        version = who.version if who is not None else ""
        checksum = who.checksum if who is not None else ""
        gene_drugs = self._gene_drugs(catalogues)

        for variant in variants:
            calls = calls_for_variant(catalogues, variant)
            variant.catalogue_calls = calls
            variant.catalogue_calls.extend(rules_engine.annotate_variant(
                variant,
                gene_drugs=gene_drugs.get((variant.gene or "").strip().lower(), ()),
                catalogue_version=version, catalogue_checksum=checksum))

        suppressions = rules_engine.epistasis_suppressions(variants)
        target_covered = self._target_coverage(catalogue_positions, pileup,
                                               translation, platform)
        drugs = sorted(set(list(catalogue_positions.keys())
                           + [d for c in catalogues.values() for d in c.drugs]),
                       key=lambda name: config.DRUGS.index(name)
                       if name in config.DRUGS else len(config.DRUGS))
        result.drugs = consensus_engine.consensus(
            variants, platform=platform, drugs=drugs,
            target_covered=target_covered, suppressions=suppressions)

        missing = [name for name in config.CATALOGUES if name not in catalogues]
        if missing:
            result.caveats.append(
                "{0} did not contribute to this report because {1} not loaded; "
                "their silence is not agreement".format(
                    " and ".join(missing), "it was" if len(missing) == 1 else "they were"))
        for suppression in suppressions:
            result.caveats.append("{0}: {1}".format(suppression.rule, suppression.why))
        result.checks.append(Check.boolean(
            "who_catalogue_anchor", who is not None,
            source=config.source_for("anchor_catalogue"), category="resistance",
            reading="WHO v2 {0} is the anchor: where it grades a variant its grade "
                    "is the call".format(version) if who is not None else
                    "the WHO catalogue was not loaded"))

    def _resistance_ntm(self, result: SampleResult, variants: List[Variant],
                        platform: str) -> None:
        """NTM resistance from the literature rule table (§5.6)."""
        assessment = self._stage(
            result, "ntm_resistance", "resistance",
            "no genotypic drug prediction was made for this organism",
            ntm_engine.call_ntm_resistance, result.species, variants,
            platform=platform)
        if assessment is None:
            return
        result.drugs = assessment.drugs
        result.checks.extend(assessment.checks)
        result.caveats.extend(assessment.caveats)
        for row in assessment.no_evidence_base:
            LOG.debug("%s: %s", result.sample_id, row["text"])
        if assessment.citations:
            result.caveats.append(
                "NTM calls rest on: {0}".format("; ".join(assessment.citations)))
        # Only where erm(41) is part of this organism's evidence base. Saying
        # "erm(41) was not typed" about M. chimaera would be true and useless;
        # saying it about M. abscessus is the single most important gap in the
        # macrolide call.
        uses_erm41 = any(
            "erm(41)" in (ntm_engine.evidence_for(result.species, call.drug).genes
                          if ntm_engine.evidence_for(result.species, call.drug) else ())
            for call in assessment.drugs)
        if uses_erm41:
            result.checks.append(Check.not_measured(
                "erm41_sequevar",
                "erm(41) sequevar typing needs the gene's coordinates in this "
                "reference, which Mjolnir does not yet ship for NTM genomes; the "
                "macrolide call therefore rests on rrl alone and inducible "
                "resistance was not assessed",
                source=config.source_for("erm41_sequevar_position"),
                category="resistance"))

    # -- contamination ------------------------------------------------------

    def _contamination(self, result: SampleResult, sample: SampleInput,
                       barcode_sites: Sequence[BarcodeSite],
                       pileup: Mapping[Tuple[str, int], PileupSite],
                       translation: Mapping[str, str], platform: str,
                       reference_gc: Optional[float]) -> None:
        """Everything §8 allows to be measured, and the verdict it supports."""
        snp_sites: List[SiteObservation] = []
        if platform != PLATFORM_FASTA:
            for variant in result.variants:
                observation = _snp_observation(variant)
                if observation is not None:
                    snp_sites.append(observation)
        lineage_sites = _lineage_observations(barcode_sites, pileup, translation) \
            if platform != PLATFORM_FASTA else []
        unambiguous_sites = [_observation_from_site(site) for site in pileup.values()
                             if site.covered] if platform != PLATFORM_FASTA else []

        screen = self._screen(result, sample, platform)
        contamination = self._stage(
            result, "contamination", "contamination",
            "no sample-validity verdict was reached",
            assess_contamination, platform=platform, qc=result.qc,
            snp_sites=snp_sites or None, lineage_sites=lineage_sites or None,
            unambiguous_sites=unambiguous_sites or None,
            reference_set_present=bool(self.references()),
            screen=screen, reference_gc=reference_gc,
            intended_use=self.options.intended_use, config=self.config)
        if contamination is None:
            return
        result.contamination = contamination
        if snp_sites or lineage_sites or unambiguous_sites:
            contamination.caveats.append(PILEUP_QUALITY_NOTE)
            contamination.caveats.append(
                "the genome-wide heterozygous-SNP denominator is the {0} SNP "
                "sites this run could examine, not every callable position in "
                "the genome".format(len(snp_sites)))
        result.qc.unambiguous_fraction = contamination.unambiguous_fraction

    def _screen(self, result: SampleResult, sample: SampleInput,
                platform: str) -> TaxonomicScreen:
        """A Kraken2 screen if one was asked for, and its refusal if it was not.

        An index that is not a mycobacterial pangenome yields an *uninformative*
        screen, not a clean one: measured Kraken2 sensitivity for
        *M. tuberculosis* reads against a standard index is 0.0731 (§8).
        """
        if self.options.kraken2_report is not None:
            text = Path(self.options.kraken2_report).read_text(errors="replace")
            return evaluate_kraken2_screen(self.config.kraken2_db,
                                           confidence=self.config.kraken2_confidence,
                                           report_text=text)
        if self.config.kraken2_db is None:
            return no_screen("no Kraken2 index was configured (--kraken2-db)")
        if platform == PLATFORM_FASTA:
            return no_screen("the read composition screen needs reads; this sample "
                             "is an assembly")
        informative, why = config.kraken2_index_informative(self.config.kraken2_db)
        if not informative:
            # Running it anyway would produce a table nobody may use, and its
            # mere presence in the report invites exactly the reading §8 forbids.
            return evaluate_kraken2_screen(self.config.kraken2_db,
                                           confidence=self.config.kraken2_confidence,
                                           report_text=None)
        report_path = Path(self.config.out_dir) / "{0}.kraken2.report".format(
            safe_name(sample.sample_id))
        text = self._stage(
            result, "kraken2_screen", "contamination",
            "the non-target read fraction is unmeasured, which is not the same "
            "as a clean sample",
            self._run_kraken2, sample, report_path, why)
        if text is None:
            return no_screen("the Kraken2 screen did not complete; see the warning above")
        return evaluate_kraken2_screen(self.config.kraken2_db,
                                       confidence=self.config.kraken2_confidence,
                                       report_text=text)

    def _run_kraken2(self, sample: SampleInput, report_path: Path, why: str) -> str:
        """Run Kraken2 at a confidence Mjolnir will accept, and read its report."""
        ensure_dir(report_path.parent)
        confidence = config.kraken2_confidence(self.config.kraken2_confidence)
        argv: List[Any] = [
            "kraken2", "--db", str(self.config.kraken2_db),
            "--threads", str(self.config.threads),
            "--confidence", "{0:g}".format(confidence),
            "--report", str(report_path), "--output", os.devnull,
        ]
        if sample.is_paired:
            argv.append("--paired")
        argv.extend(str(path) for path in sample.paths)
        LOG.info("%s: Kraken2 screen (%s)", sample.sample_id, why)
        shell.run(argv, stdout_path=None, capture=False,
                  why="read composition screen")
        self._record_tools("kraken2")
        return report_path.read_text(errors="replace")

    # -- cohort input -------------------------------------------------------

    def _collect_cohort_input(self, result: SampleResult, sample: SampleInput,
                              bam: Optional[Path], reference: Path, platform: str,
                              variants: Sequence[Variant]) -> None:
        """Keep this sample's variants and callable regions for the joint table."""
        regions = None
        if bam is not None:
            regions = self._stage(
                result, "callable_regions", "cohort",
                "this sample has no denominator, so every pairwise distance "
                "involving it is reported as not computed",
                callable_regions, bam, platform=platform,
                min_depth=self.config.degraded_depth_floor, cfg=self.config,
                name="{0} callable at >= {1}x".format(
                    result.sample_id, self.config.degraded_depth_floor))
        elif platform == PLATFORM_FASTA:
            result.checks.append(Check.not_measured(
                "callable_regions",
                "an assembly carries no depth, so its callable regions cannot be "
                "established and its pairwise distances have no denominator",
                category="cohort"))
        self.sample_variants[result.sample_id] = SampleVariants(
            sample_id=result.sample_id,
            variants=list(variants),
            callable_regions=regions,
            reference=Path(reference).name,
            platform=platform,
            note=result.contamination.verdict)

    # -- provenance ---------------------------------------------------------

    def _tool_versions(self) -> Dict[str, str]:
        """Versions of the tools that actually ran in this run.

        Recorded explicitly rather than harvested, because the engines run their
        pipelines through their own executor: a version table assembled from
        whatever happened to pass through ``shell.run`` would silently omit the
        aligner.
        """
        shell.record_versions(*sorted(self.tools_used))
        return shell.captured_versions()

    def _database_versions(self) -> List[DatabaseVersion]:
        """Version and checksum of every database this run consulted."""
        versions: List[DatabaseVersion] = []
        seen: Set[str] = set()
        for entry in catalogue_versions(self.catalogues()):
            versions.append(entry)
            seen.add(entry.name)
        for name in sorted(self.databases_used):
            try:
                entry = db_fetch.database_version(name, db_root=self.config.db_dir)
            except MjolnirError as exc:
                LOG.debug("no recorded version for database %s: %s", name, exc)
                continue
            if entry.name not in seen:
                versions.append(entry)
                seen.add(entry.name)
        return versions

    def _interpret(self, result: SampleResult) -> Interpretation:
        """The prose layer. Never fatal: a model failure is a rule-only report."""
        from .agent import interpret_sample

        try:
            return interpret_sample(result, client=self.client(),
                                    run_config=self.config)
        except MjolnirError as exc:
            LOG.warning("%s: interpretation not written: %s", result.sample_id, exc)
            return Interpretation(rule_only=True, discarded_reason=str(exc))

    # -- batch and cohort ---------------------------------------------------

    def run(self, samples: Sequence[SampleInput], workdir: Path,
            evidence: Optional[Mapping[str, Sequence[Check]]] = None
            ) -> List[SampleResult]:
        """Run every sample, keeping a failure to one sample.

        Sequential on purpose. Each stage already parallelises internally with
        ``--threads``, and a mycobacterial genome is small enough that running
        samples concurrently mostly multiplies peak memory and interleaves the
        logs of whichever one failed.
        """
        results: List[SampleResult] = []
        for index, sample in enumerate(samples, start=1):
            LOG.info("[%d/%d] %s (%s)", index, len(samples), sample.sample_id,
                     sample.platform)
            checks = list((evidence or {}).get(sample.sample_id, ()))
            try:
                results.append(self.run_sample(sample, workdir, checks))
            except MjolnirError as exc:
                LOG.error("%s: not analysed: %s", sample.sample_id, exc)
                failed = SampleResult(
                    sample_id=sample.sample_id, platform=sample.platform,
                    inputs=[str(p) for p in sample.paths],
                    mjolnir_version=__version__, profile=self.config.profile,
                    started_at=_stamp(), finished_at=_stamp(), status=STATUS_FAIL)
                failed.warnings.append("not analysed: {0}".format(exc))
                failed.checks.append(Check.not_measured(
                    "sample_analysed", "not analysed: {0}".format(exc),
                    category="qc", status=STATUS_FAIL))
                failed.checks.extend(checks)
                results.append(failed)
        return results

    def build_cohort(self, results: Sequence[SampleResult]) -> CohortResult:
        """Joint table, masked distances and clusters over the samples just run.

        Raises only when the joint table itself cannot exist — different
        references, duplicate sample ids. A one-sample cohort is not an error:
        it produces a table, no pairs, and a check saying that the absence of a
        comparison is an absence and not a distance of zero.
        """
        started = time.time()
        entries = [self.sample_variants[r.sample_id] for r in results
                   if r.sample_id in self.sample_variants]
        if not entries:
            raise MjolnirError(
                "no sample contributed variants to the cohort. Cohort mode needs "
                "the per-sample analysis to have produced a variant call and a "
                "callable-region set; check the per-sample warnings above.")

        table = build_joint_table(entries)
        checks: List[Check] = list(table.checks)
        if not any(check.name == "cohort_size" for check in checks):
            checks.append(cohort_size_check(len(entries)))
        caveats: List[str] = list(table.caveats)

        mask = self._mask(table.reference)
        matrix = distance_matrix(
            table, mask, proximity_window=config.SNP_PROXIMITY_WINDOW,
            min_shared_callable_sites=config.MIN_SHARED_CALLABLE_SITES)
        checks.extend(matrix.checks)
        caveats.extend(matrix.caveats)

        assignment = cluster_samples(
            [entry.sample_id for entry in entries], matrix.pairs,
            threshold=self.config.cluster_distance,
            threshold_basis=config.cluster_threshold_basis(self.config.cluster_distance),
            min_shared_callable_sites=config.MIN_SHARED_CALLABLE_SITES)
        checks.extend(assignment.checks)
        caveats.extend(assignment.caveats)

        reference_length = None
        for result in results:
            if result.qc.reference_length:
                reference_length = result.qc.reference_length
                break

        cohort = CohortResult(
            samples=[entry.sample_id for entry in entries],
            pairs=list(matrix.pairs),
            clusters=list(assignment.clusters),
            threshold=self.config.cluster_distance,
            threshold_basis=assignment.threshold_basis,
            mask_name=mask.name,
            masked_sites=mask.masked_bases(),
            masked_fraction=mask.fraction_of(reference_length),
            joint_sites=table.site_count,
            reference=table.reference,
            checks=checks,
            caveats=_dedupe(caveats),
            tool_versions=self._tool_versions(),
            database_versions=self._database_versions(),
            mjolnir_version=__version__,
        )
        cohort.runtime_seconds = time.time() - started
        if self.options.interpret:
            cohort.interpretation = self._interpret_cohort(cohort)
        LOG.info("cohort: %d samples, %d joint sites, %d clusters at %d SNPs (%s)",
                 len(cohort.samples), table.site_count, len(cohort.clusters),
                 self.config.cluster_distance, human_time(cohort.runtime_seconds))
        return cohort

    def _mask(self, reference: str) -> Mask:
        """The repeat/low-complexity mask, or an explicit refusal to pretend.

        Masking is mandatory before counting SNP distances (§9), so the only way
        to get an unmasked matrix is :meth:`Mask.absent`, which stamps the reason
        on every caveat list and fails the ``mask_applied`` check. It is never
        reached silently.
        """
        if self.config.mask_bed is not None:
            self.databases_used.add("operator-supplied mask")
            return load_mask(self.config.mask_bed, name=Path(self.config.mask_bed).name,
                             source="operator-supplied (--mask)")
        try:
            path = db_fetch.database_file("tbdb", "mask.bed", db_root=self.config.db_dir)
        except MjolnirError as exc:
            return Mask.absent(
                "no mask was applied: {0}. An unmasked distance over repetitive "
                "and error-prone regions is not comparable to the published SNP "
                "thresholds.".format(exc))
        self.databases_used.add("tbdb")
        return load_mask(path, name="tbdb mask.bed", source=config.SRC_TBDB)

    def _interpret_cohort(self, cohort: CohortResult) -> Interpretation:
        from .agent import interpret_cohort

        try:
            return interpret_cohort(cohort, client=self.client())
        except MjolnirError as exc:
            LOG.warning("cohort interpretation not written: %s", exc)
            return Interpretation(rule_only=True, discarded_reason=str(exc))


def _caller_binary(caller: str) -> str:
    """The executable behind a caller name, for the version table."""
    return {"clair3": "run_clair3.sh",
            "direct-comparison": "minimap2"}.get(caller, caller)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def resolve_formats(requested: Optional[Sequence[str]]) -> List[str]:
    """Normalise ``--format`` values, refusing anything not implemented."""
    if not requested:
        return ["tables", "html", "pdf"]
    out: List[str] = []
    for item in requested:
        for token in str(item).split(","):
            token = token.strip().lower()
            if not token:
                continue
            resolved = FORMAT_ALIASES.get(token)
            if resolved is None:
                raise MjolnirError(
                    "unknown --format {0!r}; choose from {1} (tsv and json are "
                    "written together)".format(
                        token, ", ".join(sorted(set(FORMAT_ALIASES)))))
            if resolved not in out:
                out.append(resolved)
    return out


def write_outputs(results: Sequence[SampleResult], out_dir: Path,
                  cohort: Optional[CohortResult] = None,
                  formats: Optional[Sequence[str]] = None,
                  profile: str = "clinical") -> List[Path]:
    """Write every requested artefact and return what was written.

    The PDF is optional at runtime — reportlab and matplotlib are extras — and a
    missing one is reported as a missing PDF rather than as a failed run, since
    the HTML report carries the same content and the TSV/JSON carry the same
    numbers. Nothing else is allowed to fail quietly: a report step that wrote
    no files while returning success is the failure this project was written
    against.
    """
    from .report import (write_cohort_html, write_cohort_pdf, write_cohort_tables,
                         write_html, write_pdf, write_sample_tables)
    from .report.pdf import generated_stamp

    wanted = resolve_formats(formats)
    directory = ensure_dir(out_dir)
    generated = generated_stamp()
    written: List[Path] = []

    if "tables" in wanted:
        for result in results:
            written.extend(write_sample_tables(directory, result, generated))
        if cohort is not None:
            written.extend(write_cohort_tables(directory, cohort, results, generated))

    if "html" in wanted:
        for result in results:
            written.append(write_html(
                directory / "{0}.html".format(safe_name(result.sample_id)),
                result, cohort=cohort, profile=profile, generated=generated))
        if cohort is not None:
            written.append(write_cohort_html(directory / "cohort.html", cohort,
                                             results, generated))

    if "pdf" in wanted:
        for result in results:
            path = directory / "{0}.pdf".format(safe_name(result.sample_id))
            try:
                written.append(write_pdf(path, result, cohort=cohort, profile=profile))
            except MjolnirError as exc:
                LOG.warning("no PDF for %s: %s", result.sample_id, exc)
        if cohort is not None:
            try:
                written.append(write_cohort_pdf(directory / "cohort.pdf", cohort, results))
            except MjolnirError as exc:
                LOG.warning("no cohort PDF: %s", exc)

    LOG.info("wrote %d file(s) to %s", len(written), directory)
    return written


# ---------------------------------------------------------------------------
# Reading a finished run back in
# ---------------------------------------------------------------------------

def _fields_of(cls: Any) -> Set[str]:
    return set(f.name for f in dataclass_fields(cls))


def _build(cls: Any, data: Mapping[str, Any]) -> Any:
    """Construct a record from a JSON object, ignoring keys it does not have.

    ``Variant.to_dict`` adds derived keys (``coordinate_key``, ``hgvs_key``) that
    are properties rather than fields, and a future version may add more, so the
    reader drops what it does not recognise instead of refusing to open a file
    written by a slightly different release.
    """
    known = _fields_of(cls)
    return cls(**dict((k, v) for k, v in data.items() if k in known))


def load_sample_result(path: Path) -> SampleResult:
    """Rebuild a :class:`~mjolnir.records.SampleResult` from its JSON artefact.

    This is what ``mjolnir report`` runs on. It reads the artefact rather than
    re-deriving anything: regenerating a report must not be able to change a
    call, and the only way to guarantee that is for the second render to see
    exactly the fields the first one wrote.
    """
    import json

    from .records import (CatalogueCall, ContaminationResult, DrugCall,
                          LineageCall, SpeciesCall)  # noqa: PLC0415 - see below

    raw = json.loads(Path(path).read_text(errors="replace"))
    if not isinstance(raw, dict) or "sample" not in raw:
        raise MjolnirError(
            "{0} is not a Mjolnir sample result: no 'sample' key. The file "
            "'mjolnir run' writes is <sample>.json in the output directory."
            .format(path))

    def checks(items: Any) -> List[Check]:
        return [_build(Check, item) for item in (items or [])]

    variants: List[Variant] = []
    for item in raw.get("variants") or []:
        variant = _build(Variant, item)
        variant.catalogue_calls = [_build(CatalogueCall, c)
                                   for c in (item.get("catalogue_calls") or [])]
        variants.append(variant)

    drugs: List[DrugCall] = []
    for item in raw.get("drugs") or []:
        drug = _build(DrugCall, item)
        drug.catalogue_calls = [_build(CatalogueCall, c)
                                for c in (item.get("catalogue_calls") or [])]
        drugs.append(drug)

    qc = _build(QCMetrics, raw.get("qc") or {})
    qc.checks = checks((raw.get("qc") or {}).get("checks"))
    contamination = _build(ContaminationResult, raw.get("contamination") or {})
    contamination.checks = checks((raw.get("contamination") or {}).get("checks"))

    result = SampleResult(
        sample_id=str(raw.get("sample") or ""),
        platform=str(raw.get("platform") or PLATFORM_ILLUMINA),
        inputs=list(raw.get("inputs") or []),
        reference=str(raw.get("reference") or ""),
        species=_build(SpeciesCall, raw.get("species") or {}),
        lineage=_build(LineageCall, raw.get("lineage") or {}),
        variants=variants,
        drugs=drugs,
        qc=qc,
        contamination=contamination,
        checks=checks(raw.get("checks")),
        caveats=list(raw.get("caveats") or []),
        warnings=list(raw.get("warnings") or []),
        tool_versions=dict(raw.get("tool_versions") or {}),
        database_versions=[_build(DatabaseVersion, d)
                           for d in (raw.get("database_versions") or [])],
        interpretation=(_build(Interpretation, raw["interpretation"])
                        if raw.get("interpretation") else None),
        mjolnir_version=str(raw.get("mjolnir_version") or ""),
        profile=str(raw.get("profile") or "clinical"),
        started_at=str(raw.get("started_at") or ""),
        finished_at=str(raw.get("finished_at") or ""),
        runtime_seconds=float(raw.get("runtime_seconds") or 0.0),
        status=str(raw.get("status") or STATUS_WARN),
    )
    if result.mjolnir_version and result.mjolnir_version != __version__:
        LOG.warning("%s was written by Mjolnir %s and is being re-rendered by %s",
                    path, result.mjolnir_version, __version__)
    return result


def load_cohort_result(path: Path) -> CohortResult:
    """Rebuild a :class:`~mjolnir.records.CohortResult` from ``cohort.json``."""
    import json

    from .records import Cluster, PairwiseDistance

    raw = json.loads(Path(path).read_text(errors="replace"))
    if not isinstance(raw, dict) or "samples" not in raw:
        raise MjolnirError(
            "{0} is not a Mjolnir cohort result: no 'samples' key".format(path))
    mask = raw.get("mask") or {}
    return CohortResult(
        samples=list(raw.get("samples") or []),
        pairs=[_build(PairwiseDistance, p) for p in (raw.get("pairs") or [])],
        clusters=[_build(Cluster, c) for c in (raw.get("clusters") or [])],
        threshold=raw.get("threshold"),
        threshold_basis=str(raw.get("threshold_basis") or ""),
        mask_name=str(mask.get("name") or ""),
        masked_sites=mask.get("masked_sites"),
        masked_fraction=mask.get("masked_fraction"),
        joint_sites=raw.get("joint_sites"),
        reference=str(raw.get("reference") or ""),
        checks=[_build(Check, c) for c in (raw.get("checks") or [])],
        caveats=list(raw.get("caveats") or []),
        warnings=list(raw.get("warnings") or []),
        tool_versions=dict(raw.get("tool_versions") or {}),
        database_versions=[_build(DatabaseVersion, d)
                           for d in (raw.get("database_versions") or [])],
        interpretation=(_build(Interpretation, raw["interpretation"])
                        if raw.get("interpretation") else None),
        mjolnir_version=str(raw.get("mjolnir_version") or ""),
        runtime_seconds=float(raw.get("runtime_seconds") or 0.0),
    )


__all__ = [
    "NO_ANNOTATION_NOTE", "NO_TRIMMING_NOTE", "OUTPUT_FORMATS",
    "PILEUP_QUALITY_NOTE", "Pipeline", "RunOptions", "callable_regions",
    "contig_translation", "load_cohort_result", "load_sample_result",
    "read_fai", "resolve_formats", "write_outputs",
]
