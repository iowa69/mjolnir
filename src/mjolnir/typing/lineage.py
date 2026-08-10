"""MTBC lineage and sublineage from a SNP barcode, genotyped from the pileup.

Species identification stops at the complex for MTBC — the members are later
heterotypic synonyms of *M. tuberculosis* at 99.21-99.92% ANI, so no distance
can separate them (see ``typing/species.py``). The member is a phylogenetic
question, and this module answers it the only way it can be answered: from
lineage-defining SNPs, at positions read out of ``barcode.bed``.

Four decisions here are the ones that matter.

**The genotype comes from a direct pileup, not from the variant caller.**
Following pathogen-profiler, each barcode position is genotyped from base counts
rather than from a VCF. A variant caller's filters are tuned for calling
variants, and a barcode site that a caller declined to emit is indistinguishable
in the VCF from a site that matched the reference — which turns a coverage gap
into a confident wrong lineage. On ONT the highest-depth allele wins outright,
because ONT per-read error makes a fraction threshold at a single site a worse
estimator than the mode.

**Support is reported, not just the label.** Every call carries sites supporting
over sites callable over sites the scheme defines for that taxon. A lineage
called from 4 of 5 defining sites and one called from 58 of 60 are different
claims, and only the denominators say so.

**BCG is flagged explicitly**, because it is intrinsically pyrazinamide-resistant
and that is a property of the organism rather than anything a mutation catalogue
will report. A pyrazinamide result on a BCG isolate that does not say so is
clinically misleading.

**Every *M. bovis* call carries a confidence caveat.** *M. bovis* is defined by
very few phylogenetic SNPs — 23 in SNP-IT, and in the tbdb scheme its taxa carry
four or five defining sites each — so a coverage gap over a handful of positions,
or a contaminating read set, is enough to lose the call or to invent it.

The scheme is also honest about what it cannot do: tbdb's ``barcode.bed`` carries
defining sites for *M. bovis* (La1), *M. caprae* (La2) and *M. orygis* (La3), and
none at all for *M. microti* or *M. pinnipedii*. Absence of a *M. microti* call
from this scheme is absence of a barcode, not absence of the organism, and every
MTBC call says so.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from ..config import (
    BCG_PZA_NOTE,
    FASTA_CAPABILITY_LOSS,
    MBOVIS_CAVEAT,
    MBOVIS_DEFINING_SNPS,
    MIN_BARCODE_CALLABLE_FRACTION,
    MIN_BARCODE_SUPPORT_FRACTION,
    MIN_MINOR_VARIANT_FRACTION,
    ONT_MINOR_VARIANT_CAVEAT,
    min_reads_for,
    source_for,
)
from ..records import (
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MODERATE,
    CONFIDENCE_NONE,
    PLATFORM_FASTA,
    PLATFORM_ONT,
    STATUS_WARN,
    Check,
    LineageCall,
    normalise_platform,
)
from ..utils import LOG, MjolnirError, PathLike, natural_key, require_database, safe_fraction

# ---------------------------------------------------------------------------
# The pileup interface
# ---------------------------------------------------------------------------

#: What ``engines/pileup.py`` hands over: base counts per reference position,
#: keyed by ``(chrom, 1-based position)``. Anything that is not A, C, G or T is
#: ignored when depth is computed — Mjolnir uses a conventional ACGT depth
#: outside ``--mtbseq-compat``, where MTBseq's denominator includes N and GAP
#: (design §9b) and therefore produces different fractions at the same position.
PileupCounts = Mapping[Tuple[str, int], Mapping[str, int]]

BASES: Tuple[str, ...] = ("A", "C", "G", "T")

#: IUPAC expansion for the barcode's allele column. All 1,111 rows of the tbdb
#: file are unambiguous today, but the BED format permits an ambiguity code and
#: pathogen-profiler expands them, so a scheme that starts using one must not
#: silently stop matching.
IUPAC_BASES: Dict[str, Tuple[str, ...]] = {
    "A": ("A",), "C": ("C",), "G": ("G",), "T": ("T",), "U": ("T",),
    "R": ("A", "G"), "Y": ("C", "T"), "S": ("C", "G"), "W": ("A", "T"),
    "K": ("G", "T"), "M": ("A", "C"),
    "B": ("C", "G", "T"), "D": ("A", "G", "T"),
    "H": ("A", "C", "T"), "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T"),
}

#: An assembly contributes exactly one observation per position. This is not a
#: threshold that could be tuned — it is what a consensus base is — so it is
#: named here rather than registered in config.py.
FASTA_SITE_OBSERVATIONS = 1


def major_allele(counts: Mapping[str, int]) -> Tuple[str, int, bool]:
    """The highest-depth ACGT allele, the ACGT depth, and whether it tied.

    A tie returns an empty allele: two bases at equal depth is not an
    observation of either, and resolving it by alphabetical order would make a
    coin toss look like a genotype. The depth is still returned, because the
    caller needs to distinguish "tied at 40x" from "no coverage".
    """
    depth = 0
    best_base = ""
    best_count = 0
    tied = False
    for base in BASES:
        count = int(counts.get(base, 0) or 0)
        depth += count
        if count > best_count:
            best_base, best_count, tied = base, count, False
        elif count == best_count and count > 0:
            tied = True
    if tied or best_count <= 0:
        return "", depth, tied
    return best_base, depth, False


# ---------------------------------------------------------------------------
# barcode.bed
# ---------------------------------------------------------------------------

#: tbdb's ``barcode.bed`` is an 8-column BED, verified against
#: jodyphelan/tbdb@barcode.bed on 2026-08-10: 1,111 single-base rows, 126 taxa,
#: all on the contig named ``Chromosome``. The columns are
#: ``chrom, start(0-based), end(1-based), taxon, derived allele, family,
#: spoligotype families, RD region``, with the literal string ``None`` used for
#: an empty field. pathogen-profiler reads ``end`` as the position and expands
#: the allele through IUPAC, and Mjolnir does the same so the two tools genotype
#: the same base.
BARCODE_MIN_COLUMNS = 5
BARCODE_FULL_COLUMNS = 8
BARCODE_EMPTY_TOKEN = "None"

BARCODE_FILE_NAME = "barcode.bed"
BARCODE_DIR_NAME = "tbdb"
BARCODE_FETCH_HINT = "mjolnir db fetch tbdb"

#: Contig names H37Rv travels under. tbdb's BED says ``Chromosome``; the WHO
#: catalogue, the FASTA a BAM was mapped against and the pileup keys usually say
#: ``NC_000962.3``. The names are aliases for one sequence, and a barcode that
#: silently matched nothing because of the spelling would report every lineage as
#: uncallable.
H37RV_CHROM_ALIASES: Tuple[str, ...] = (
    "Chromosome", "NC_000962.3", "NC_000962", "AL123456.3", "AL123456",
    "H37Rv", "MTB_anc", "NC000962.3",
)


@dataclass(frozen=True)
class BarcodeSite:
    """One lineage-defining position, as ``barcode.bed`` states it."""

    chrom: str
    pos: int
    taxon: str
    allele: str
    family: str = ""
    spoligotype: str = ""
    rd: str = ""
    line_number: int = 0

    @property
    def key(self) -> Tuple[str, int]:
        return (self.chrom, self.pos)

    @property
    def alleles(self) -> Tuple[str, ...]:
        """The allele column expanded through IUPAC."""
        return IUPAC_BASES.get(self.allele.upper(), (self.allele.upper(),))


def barcode_path(db_dir: PathLike) -> Path:
    """``<db>/tbdb/barcode.bed``, or a MjolnirError naming the fetch command."""
    return require_database(
        Path(db_dir) / BARCODE_DIR_NAME / BARCODE_FILE_NAME,
        "the tbdb lineage barcode (barcode.bed)",
        BARCODE_FETCH_HINT,
    )


def load_barcode(path: PathLike) -> List[BarcodeSite]:
    """Parse ``barcode.bed``.

    Strict about the columns it needs and forgiving about the ones it does not:
    a scheme that grows a ninth column keeps working, one that drops the allele
    column fails immediately and names the line. The BED ``end`` coordinate is
    the 1-based position, which is what a pileup is keyed on.
    """
    resolved = Path(path)
    if not resolved.exists():
        raise MjolnirError(
            "lineage barcode not found: {0}\n  fetch it with: {1}".format(
                resolved, BARCODE_FETCH_HINT)
        )
    sites: List[BarcodeSite] = []
    with open(str(resolved), "rt", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip() or line.lstrip().startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < BARCODE_MIN_COLUMNS:
                raise MjolnirError(
                    "{0} line {1}: expected at least {2} tab-separated columns "
                    "(chrom, start, end, taxon, allele), found {3}".format(
                        resolved, line_number, BARCODE_MIN_COLUMNS, len(fields))
                )
            try:
                position = int(fields[2])
            except ValueError:
                raise MjolnirError(
                    "{0} line {1}: BED end coordinate {2!r} is not an "
                    "integer".format(resolved, line_number, fields[2])
                )
            allele = fields[4].strip().upper()
            if allele not in IUPAC_BASES:
                raise MjolnirError(
                    "{0} line {1}: allele {2!r} is not a nucleotide or IUPAC "
                    "code".format(resolved, line_number, fields[4])
                )
            taxon = fields[3].strip()
            if not taxon:
                raise MjolnirError(
                    "{0} line {1}: the taxon column is empty".format(resolved, line_number))
            sites.append(BarcodeSite(
                chrom=fields[0].strip(),
                pos=position,
                taxon=taxon,
                allele=allele,
                family=_clean(fields[5]) if len(fields) > 5 else "",
                spoligotype=_clean(fields[6]) if len(fields) > 6 else "",
                rd=_clean(fields[7]) if len(fields) > 7 else "",
                line_number=line_number,
            ))
    if not sites:
        raise MjolnirError("{0} contains no barcode rows".format(resolved))
    LOG.debug("loaded %d barcode sites covering %d taxa from %s",
              len(sites), len({s.taxon for s in sites}), resolved)
    return sites


def _clean(value: str) -> str:
    """Strip the literal ``None`` tbdb writes for an empty field."""
    text = str(value or "").strip()
    return "" if text in ("", BARCODE_EMPTY_TOKEN, ".") else text


def barcode_positions(sites: Sequence[BarcodeSite]) -> List[Tuple[str, int]]:
    """Sorted, de-duplicated positions for ``engines/pileup.py`` to pile up."""
    return sorted({site.key for site in sites})


def scheme_description(sites: Sequence[BarcodeSite], name: str = "") -> str:
    """One line naming the scheme, its size and its taxon count."""
    taxa = {site.taxon for site in sites}
    return "{0} ({1} defining sites, {2} taxa)".format(
        name or "SNP barcode", len(sites), len(taxa))


# ---------------------------------------------------------------------------
# Taxon hierarchy
# ---------------------------------------------------------------------------

#: ``lineage4.2.2.1`` and ``La1.2.BCG`` are hierarchical; ``M.canetti`` is not.
#: The pattern is what tells them apart, so a taxon whose name merely contains a
#: dot is not chopped into a fake ancestry.
_HIERARCHICAL = re.compile(r"^(lineage|La)([0-9]+)((?:\.[A-Za-z0-9_\-]+)*)$")


def taxon_components(taxon: str) -> List[str]:
    """The dot-separated levels of a hierarchical taxon, or the whole name."""
    if _HIERARCHICAL.match(str(taxon or "")):
        return str(taxon).split(".")
    return [str(taxon or "")]


def taxon_depth(taxon: str) -> int:
    return len(taxon_components(taxon))


def taxon_root(taxon: str) -> str:
    """``lineage4.3.1`` -> ``lineage4``; ``La1.2.BCG`` -> ``La1``."""
    return taxon_components(taxon)[0]


def taxon_ancestors(taxon: str) -> List[str]:
    """Every ancestor of a taxon, shallowest first, excluding the taxon itself."""
    parts = taxon_components(taxon)
    return [".".join(parts[:i]) for i in range(1, len(parts))]


def is_descendant(child: str, ancestor: str) -> bool:
    return ancestor in taxon_ancestors(child)


# ---------------------------------------------------------------------------
# Animal lineages and BCG
# ---------------------------------------------------------------------------

#: The animal-adapted members a clinician expects an MTBC report to distinguish.
#: La1/La2/La3 is the Zwyer et al. 2021 nomenclature tbdb adopts verbatim.
MTBC_ANIMAL_MEMBERS: Tuple[str, ...] = (
    "Mycobacterium bovis",
    "Mycobacterium caprae",
    "Mycobacterium microti",
    "Mycobacterium orygis",
    "Mycobacterium pinnipedii",
)

#: Roots that are animal-adapted in the tbdb scheme, with the member each denotes.
ANIMAL_ROOTS: Dict[str, str] = {
    "La1": "Mycobacterium bovis",
    "La2": "Mycobacterium caprae",
    "La3": "Mycobacterium orygis",
}

#: The BCG marker in the scheme: the taxon id ends ``.BCG`` and the spoligotype
#: column reads ``BCG``. Both are checked, because either alone would break on a
#: scheme revision that changed one of them.
BCG_TAXON_SUFFIX = ".BCG"
BCG_SPOLIGOTYPE = "BCG"

#: The sentence attached to every MTBC call, naming what the loaded scheme
#: cannot call at all. Formatted with the missing members.
SCHEME_GAP_TEXT = (
    "the loaded barcode scheme carries no defining sites for {members}, so "
    "absence of such a call is absence of a barcode rather than absence of the "
    "organism"
)


def expand_member_name(label: str) -> str:
    """``M.bovis`` -> ``Mycobacterium bovis``; anything else is returned as given."""
    text = str(label or "").strip()
    if text.lower().startswith("m.") and len(text) > 2 and text[2] != " ":
        return "Mycobacterium " + text[2:].strip()
    return text


def animal_member_for(taxon: str, family: str = "") -> str:
    """The animal-adapted member a taxon denotes, or "" when it is not one."""
    named = expand_member_name(family)
    if named in MTBC_ANIMAL_MEMBERS:
        return named
    return ANIMAL_ROOTS.get(taxon_root(taxon), "")


def scheme_animal_coverage(
        sites: Sequence[BarcodeSite]) -> Tuple[List[str], List[str]]:
    """Which animal-adapted members the loaded scheme can and cannot call."""
    present: Set[str] = set()
    for site in sites:
        member = animal_member_for(site.taxon, site.family)
        if member:
            present.add(member)
    covered = [m for m in MTBC_ANIMAL_MEMBERS if m in present]
    missing = [m for m in MTBC_ANIMAL_MEMBERS if m not in present]
    return covered, missing


# ---------------------------------------------------------------------------
# Genotyping
# ---------------------------------------------------------------------------

@dataclass
class SiteGenotype:
    """One barcode position as it was actually observed."""

    site: BarcodeSite
    depth: int = 0
    observed: str = ""
    target_reads: int = 0
    target_fraction: Optional[float] = None
    is_callable: bool = False
    supports: bool = False
    #: Target allele present but not the majority, above the minor-variant floor.
    minor_support: bool = False
    tie: bool = False
    counts: Dict[str, int] = field(default_factory=dict)

    def as_row(self) -> Dict[str, Any]:
        """The annex row: position, expected allele, observed allele, depth, AF."""
        return {
            "taxon": self.site.taxon,
            "family": self.site.family,
            "rd": self.site.rd,
            "chrom": self.site.chrom,
            "pos": self.site.pos,
            "expected": self.site.allele,
            "observed": self.observed or "",
            "depth": self.depth,
            "target_reads": self.target_reads,
            "allele_fraction": (None if self.target_fraction is None
                                else round(self.target_fraction, 4)),
            "callable": self.is_callable,
            "supports": self.supports,
            "minor_support": self.minor_support,
            "tie": self.tie,
        }


def _chrom_lookup(sites: Sequence[BarcodeSite],
                  counts: PileupCounts) -> Dict[str, str]:
    """Map each barcode contig name onto the name the pileup uses.

    tbdb writes ``Chromosome``; a BAM built against the NCBI FASTA says
    ``NC_000962.3``. They are the same sequence, so the alias is resolved rather
    than reported as a genome-wide coverage gap. When the pileup uses a name that
    is not a known H37Rv alias and is not the barcode's own name, nothing is
    guessed: the sites are left uncallable and the caller sees a callable
    fraction of zero, which is the truthful outcome for a barcode applied to the
    wrong reference.
    """
    pileup_chroms = {chrom for chrom, _pos in counts.keys()}
    mapping: Dict[str, str] = {}
    for chrom in {site.chrom for site in sites}:
        if chrom in pileup_chroms:
            mapping[chrom] = chrom
            continue
        aliases = [c for c in pileup_chroms if c in H37RV_CHROM_ALIASES]
        if chrom in H37RV_CHROM_ALIASES and len(aliases) == 1:
            mapping[chrom] = aliases[0]
            LOG.debug("barcode contig %r resolved to pileup contig %r", chrom, aliases[0])
        elif len(pileup_chroms) == 1:
            only = next(iter(pileup_chroms))
            mapping[chrom] = only
            LOG.debug("barcode contig %r resolved to the pileup's only contig %r",
                      chrom, only)
        else:
            mapping[chrom] = chrom
    return mapping


def genotype_barcode_sites(sites: Sequence[BarcodeSite], counts: PileupCounts,
                           platform: str, *,
                           min_site_depth: Optional[int] = None,
                           minor_fraction: float = MIN_MINOR_VARIANT_FRACTION
                           ) -> List[SiteGenotype]:
    """Genotype every barcode position from the pileup.

    The rule is the same on both read platforms — the highest-depth ACGT allele
    wins — and the difference is the depth a site must reach before it is
    genotyped at all: the published clinical-DST minima, 3 reads on Illumina and
    5 on ONT. Taking the mode rather than applying a fraction threshold is what
    the design asks for on ONT, and applying it on Illumina too keeps one
    definition of "the allele at this site" in the codebase instead of two.

    ``minor_support`` records the other case: the derived allele is present but
    is not the majority. On Illumina that is a mixed-lineage signal worth
    reporting; on ONT it is recorded and caveated, because minor variants are
    measurably under-detected there.
    """
    plat = normalise_platform(platform)
    if min_site_depth is None:
        min_site_depth = (FASTA_SITE_OBSERVATIONS if plat == PLATFORM_FASTA
                          else min_reads_for(plat))
    chrom_map = _chrom_lookup(sites, counts)
    genotypes: List[SiteGenotype] = []
    for site in sites:
        key = (chrom_map.get(site.chrom, site.chrom), site.pos)
        site_counts = dict(counts.get(key) or {})
        observed, depth, tie = major_allele(site_counts)
        targets = site.alleles
        target_reads = sum(int(site_counts.get(base, 0) or 0) for base in targets)
        fraction = safe_fraction(target_reads, depth)
        is_callable = depth >= min_site_depth and bool(observed)
        supports = bool(is_callable and observed in targets
                        and target_reads >= min_site_depth)
        minor = bool(is_callable and not supports and fraction is not None
                     and fraction >= minor_fraction)
        genotypes.append(SiteGenotype(
            site=site, depth=depth, observed=observed, target_reads=target_reads,
            target_fraction=fraction, is_callable=is_callable, supports=supports,
            minor_support=minor, tie=tie,
            counts=dict((b, int(site_counts.get(b, 0) or 0)) for b in BASES),
        ))
    return genotypes


@dataclass
class TaxonSupport:
    """Everything the barcode says about one taxon in one sample."""

    taxon: str
    family: str = ""
    spoligotype: str = ""
    rd: str = ""
    total: int = 0
    callable_sites: int = 0
    supporting: int = 0
    minor: int = 0
    genotypes: List[SiteGenotype] = field(default_factory=list)

    @property
    def support_fraction(self) -> Optional[float]:
        return safe_fraction(self.supporting, self.callable_sites)

    @property
    def callable_fraction(self) -> Optional[float]:
        return safe_fraction(self.callable_sites, self.total)

    @property
    def minor_fraction(self) -> Optional[float]:
        return safe_fraction(self.minor, self.callable_sites)

    @property
    def is_supported(self) -> bool:
        """Enough sites callable, and enough of those carrying the derived allele.

        Both fractions, not one: a taxon whose sites all carry the derived
        allele but of which only one was callable has not been demonstrated, and
        a taxon fully covered but half-supported has been contradicted.
        """
        callable_fraction = self.callable_fraction
        support_fraction = self.support_fraction
        if callable_fraction is None or support_fraction is None:
            return False
        return (callable_fraction >= MIN_BARCODE_CALLABLE_FRACTION
                and support_fraction >= MIN_BARCODE_SUPPORT_FRACTION)

    @property
    def is_minor_supported(self) -> bool:
        """The derived alleles are there, but as a minority population."""
        minor_fraction = self.minor_fraction
        callable_fraction = self.callable_fraction
        if minor_fraction is None or callable_fraction is None:
            return False
        return (callable_fraction >= MIN_BARCODE_CALLABLE_FRACTION
                and minor_fraction >= MIN_BARCODE_SUPPORT_FRACTION)


def summarise_taxa(genotypes: Sequence[SiteGenotype]) -> Dict[str, TaxonSupport]:
    """Group genotyped sites by the taxon they define."""
    taxa: Dict[str, TaxonSupport] = {}
    for genotype in genotypes:
        site = genotype.site
        entry = taxa.get(site.taxon)
        if entry is None:
            entry = TaxonSupport(taxon=site.taxon, family=site.family,
                                 spoligotype=site.spoligotype, rd=site.rd)
            taxa[site.taxon] = entry
        entry.total += 1
        entry.genotypes.append(genotype)
        if genotype.is_callable:
            entry.callable_sites += 1
        if genotype.supports:
            entry.supporting += 1
        if genotype.minor_support:
            entry.minor += 1
    return taxa


def family_for(taxa: Mapping[str, TaxonSupport], taxon: str) -> str:
    """The named family for a taxon (``lineage2`` -> ``East-Asian``).

    The scheme file carries the mapping, so Mjolnir ships no hand-built lineage
    -> family table. Where the taxon itself has no family the nearest ancestor
    that does supplies it, which is how ``lineage4.3.4.2.1`` inherits
    Euro-American.
    """
    for candidate in [taxon] + list(reversed(taxon_ancestors(taxon))):
        entry = taxa.get(candidate)
        if entry is not None and entry.family:
            return entry.family
    return ""


# ---------------------------------------------------------------------------
# The call
# ---------------------------------------------------------------------------

def _select(taxa: Mapping[str, TaxonSupport]) -> Tuple[str, List[str], List[str]]:
    """Pick the most specific defensible taxon, and note what contradicts it.

    Depth first, then weight of evidence. A deep taxon is only accepted if no
    ancestor of it was callable and contradicted — a sample supporting
    ``lineage4.3.4.2`` while actively refuting ``lineage4.3`` has not been
    typed, it has been confused, and the call falls back to the deepest
    consistent ancestor rather than trusting the specific label.
    """
    supported = sorted(
        [name for name, entry in taxa.items() if entry.is_supported],
        key=lambda n: (-taxon_depth(n), -(taxa[n].supporting), natural_key(n)),
    )
    if not supported:
        return "", [], []
    supported_set = set(supported)
    for candidate in supported:
        chain = taxon_ancestors(candidate) + [candidate]
        contradicted = [
            a for a in taxon_ancestors(candidate)
            if a in taxa and a not in supported_set and taxa[a].callable_sites > 0
        ]
        if contradicted:
            continue
        off_chain = sorted(
            [n for n in supported_set
             if n not in chain and not is_descendant(n, candidate)],
            key=natural_key,
        )
        return candidate, chain, off_chain
    # Every supported taxon is contradicted by one of its own ancestors. Report
    # the shallowest supported taxon and let the conflict list say the rest.
    fallback = sorted(supported, key=lambda n: (taxon_depth(n), natural_key(n)))[0]
    chain = [a for a in taxon_ancestors(fallback) if a in supported_set] + [fallback]
    off_chain = sorted([n for n in supported_set if n not in chain], key=natural_key)
    return fallback, chain, off_chain


def call_lineage(sites: Sequence[BarcodeSite], counts: PileupCounts,
                 platform: str, *, scheme: str = "tbdb barcode.bed",
                 min_site_depth: Optional[int] = None,
                 genotypes: Optional[Sequence[SiteGenotype]] = None) -> LineageCall:
    """Lineage, sublineage, animal lineage and BCG, with the support behind them.

    *genotypes* is accepted so a caller that already piled up the barcode
    positions for the contamination screen — F2/F47 uses the same lineage-defining
    sites — does not pile them up twice.
    """
    plat = normalise_platform(platform)
    if genotypes is None:
        genotypes = genotype_barcode_sites(sites, counts, plat,
                                           min_site_depth=min_site_depth)
    taxa = summarise_taxa(genotypes)
    chosen, chain, off_chain = _select(taxa)

    callable_total = sum(1 for g in genotypes if g.is_callable)
    scheme_text = "{0}; {1} of {2} sites callable in this sample".format(
        scheme_description(sites, scheme), callable_total, len(genotypes))

    caveats: List[str] = []
    _covered, missing = scheme_animal_coverage(sites)
    if missing:
        caveats.append(SCHEME_GAP_TEXT.format(members=", ".join(missing)))
    if plat == PLATFORM_FASTA:
        caveats.append(FASTA_CAPABILITY_LOSS)

    if not chosen:
        caveats.append(
            "no taxon reached the barcode support thresholds ({0:.0%} of its "
            "defining sites callable and {1:.0%} of those carrying the derived "
            "allele); the lineage is not determined, which is not the same as "
            "lineage 4".format(MIN_BARCODE_CALLABLE_FRACTION,
                               MIN_BARCODE_SUPPORT_FRACTION)
        )
        return LineageCall(
            lineage="", sublineage="",
            barcode_sites_supporting=0,
            barcode_sites_total=len(genotypes),
            barcode_sites_callable=callable_total,
            caveats=caveats, scheme=scheme_text,
            method=_method_text(plat), confidence=CONFIDENCE_NONE,
            support=[g.as_row() for g in genotypes if g.is_callable and g.supports],
        )

    chain_entries = [taxa[name] for name in chain if name in taxa]
    total = sum(e.total for e in chain_entries)
    callable_sites = sum(e.callable_sites for e in chain_entries)
    supporting = sum(e.supporting for e in chain_entries)

    root = taxon_root(chosen)
    family = family_for(taxa, chosen)
    entry = taxa[chosen]
    is_bcg = bool(chosen.endswith(BCG_TAXON_SUFFIX)
                  or entry.spoligotype.upper() == BCG_SPOLIGOTYPE)
    member = animal_member_for(chosen, family)
    is_animal = bool(member) or is_bcg
    animal_variant = ""
    if is_bcg:
        animal_variant = "BCG"
    elif member:
        animal_variant = member.split()[-1]

    mixed = list(off_chain)
    minor_mixed = sorted(
        [name for name, taxon in taxa.items()
         if taxon.is_minor_supported and name not in chain and name not in off_chain],
        key=natural_key,
    )
    mixed.extend(minor_mixed)

    confidence = _confidence(chain_entries, mixed)
    if is_animal and confidence == CONFIDENCE_HIGH:
        # Capped deliberately. The animal taxa are the thinnest in the scheme —
        # La2, La3 and La1.2.BCG carry five defining sites each — so even a
        # perfect 5-of-5 is not the same weight of evidence as a human-lineage
        # call resting on tens of positions, and printing "high" beside the
        # coverage caveat below would contradict it.
        confidence = CONFIDENCE_MODERATE

    if mixed:
        caveats.append(
            "defining sites for more than one taxon were supported ({0}); this "
            "is a mixed-infection or contamination signal and the lineage is "
            "reported with it rather than resolved by taking the larger".format(
                ", ".join(mixed[:6]))
        )
    if minor_mixed:
        caveats.append(
            "the additional support for {0} is minority signal (derived allele "
            "present below the majority threshold)".format(", ".join(minor_mixed[:6]))
        )
    if plat == PLATFORM_ONT:
        caveats.append(ONT_MINOR_VARIANT_CAVEAT)

    if is_bcg:
        caveats.append(BCG_PZA_NOTE)
    if member == "Mycobacterium bovis" or is_bcg:
        caveats.append(MBOVIS_CAVEAT)
        caveats.append(
            "in this scheme {0} is defined by {1} site(s) and its whole La1 "
            "chain by {2}; SNP-IT defines M. bovis on {3} phylogenetic SNPs, so "
            "a coverage gap over a handful of positions moves the call".format(
                chosen, entry.total, total, MBOVIS_DEFINING_SNPS)
        )
    elif is_animal:
        caveats.append(
            "{0} is defined by {1} site(s) in this scheme and its whole chain by "
            "{2}; animal-lineage calls rest on far fewer positions than a human "
            "lineage call and are correspondingly sensitive to coverage".format(
                chosen, entry.total, total)
        )
    if callable_sites < total:
        caveats.append(
            "{0} of the {1} sites defining this call were not callable at the "
            "depth required".format(total - callable_sites, total)
        )
    # A supported taxon below the chosen one, reached only by stepping over an
    # ancestor the data refutes. Held back rather than reported, and said so:
    # the more specific label would be the more interesting answer, and that is
    # exactly why it must not be given on inconsistent evidence.
    held_back = sorted(
        [name for name, entry in taxa.items()
         if entry.is_supported and name not in chain and is_descendant(name, chosen)],
        key=natural_key,
    )
    if held_back:
        caveats.append(
            "defining sites for {0} were also supported, but an intervening "
            "ancestor was callable and contradicted, so the call is held at the "
            "deepest consistent level ({1})".format(", ".join(held_back[:4]), chosen)
        )

    support_rows = [g.as_row() for g in genotypes
                    if g.site.taxon in chain or g.site.taxon in mixed]
    support_rows.sort(key=lambda row: (natural_key(row["taxon"]), row["pos"]))

    return LineageCall(
        lineage=root,
        sublineage=chosen if chosen != root else "",
        barcode_sites_supporting=supporting,
        barcode_sites_total=total,
        barcode_sites_callable=callable_sites,
        is_bcg=is_bcg,
        is_animal=is_animal,
        animal_variant=animal_variant,
        caveats=caveats,
        scheme=scheme_text,
        method=_method_text(plat),
        confidence=confidence,
        mixed_lineages=mixed,
        support=support_rows,
    )


def _method_text(platform: str) -> str:
    """What the call was made from, said plainly enough for the methods annex."""
    base = ("direct pileup at barcode positions, highest-depth allele per site "
            "(not from the variant caller)")
    if platform == PLATFORM_ONT:
        return base + "; ONT minimum 5 supporting reads per site"
    if platform == PLATFORM_FASTA:
        return ("consensus base per barcode position from the assembly; no allele "
                "fractions, so no mixed-lineage signal is available")
    return base + "; Illumina minimum 3 supporting reads per site"


def _confidence(chain_entries: Sequence[TaxonSupport], mixed: Sequence[str]) -> str:
    """How much the barcode support justifies claiming.

    Capped at moderate whenever another taxon is also supported: a mixed signal
    is a reason to doubt the label, and a "high confidence" beside a reported
    mixture would contradict the caveat sitting next to it.
    """
    if not chain_entries:
        return CONFIDENCE_NONE
    total = sum(e.total for e in chain_entries)
    callable_sites = sum(e.callable_sites for e in chain_entries)
    supporting = sum(e.supporting for e in chain_entries)
    callable_fraction = safe_fraction(callable_sites, total)
    support_fraction = safe_fraction(supporting, callable_sites)
    if callable_fraction is None or support_fraction is None:
        return CONFIDENCE_NONE
    if mixed:
        return CONFIDENCE_LOW if support_fraction < 1.0 else CONFIDENCE_MODERATE
    if callable_fraction >= 1.0 and support_fraction >= 1.0:
        return CONFIDENCE_HIGH
    if callable_fraction >= MIN_BARCODE_CALLABLE_FRACTION \
            and support_fraction >= MIN_BARCODE_SUPPORT_FRACTION:
        return CONFIDENCE_MODERATE
    return CONFIDENCE_LOW


# ---------------------------------------------------------------------------
# Non-MTBC input
# ---------------------------------------------------------------------------

#: What an NTM isolate gets instead of a lineage. There is no equivalent of the
#: MTBC barcode for the NTM species Mjolnir reports, and inventing a "lineage 1"
#: for M. chimaera would be a fabrication. Subspecies resolution for NTM lives in
#: ``typing/species.py``, where it is done on ANI plus marker SNPs.
NTM_NO_BARCODE_TEXT = (
    "the MTBC lineage barcode does not apply to {species}; no validated "
    "genome-wide lineage scheme is implemented for it, so no lineage is "
    "reported. This is an absence of a scheme, not an absence of structure"
)


def lineage_not_applicable(species: str, reason: str = "") -> LineageCall:
    """The lineage record for an isolate the MTBC barcode does not describe."""
    return LineageCall(
        lineage="", sublineage="", scheme="", method="not applicable",
        confidence=CONFIDENCE_NONE,
        caveats=[reason or NTM_NO_BARCODE_TEXT.format(species=species or "this isolate")],
    )


# ---------------------------------------------------------------------------
# Checks for the report
# ---------------------------------------------------------------------------

def lineage_checks(call: LineageCall) -> List[Check]:
    """Rule-derived verdicts over the barcode evidence, with their sources."""
    checks: List[Check] = [
        Check.numeric(
            "lineage_barcode_callable_fraction",
            call.callable_fraction,
            warn_minimum=float(MIN_BARCODE_CALLABLE_FRACTION),
            source=source_for("min_barcode_callable_fraction"),
            unit="fraction",
            category="typing",
            reading=("{0} of {1} sites defining this lineage were callable".format(
                call.barcode_sites_callable, call.barcode_sites_total)),
            not_measured_why=("the barcode scheme defines no sites for the "
                              "reported taxon, so support cannot be computed"),
        ),
        Check.numeric(
            "lineage_barcode_support_fraction",
            call.support_fraction,
            warn_minimum=float(MIN_BARCODE_SUPPORT_FRACTION),
            source=source_for("min_barcode_support_fraction"),
            unit="fraction",
            category="typing",
            reading=("{0} of {1} callable defining sites carried the derived "
                     "allele".format(call.barcode_sites_supporting,
                                     call.barcode_sites_callable)),
            not_measured_why=("no defining site for the reported taxon was "
                              "callable, so no support fraction exists"),
        ),
        Check.boolean(
            "lineage_single_taxon",
            not call.mixed_lineages,
            expected=True,
            source=source_for("barcode_from_pileup"),
            category="typing",
            reading=("one taxon's defining sites were supported"
                     if not call.mixed_lineages
                     else "defining sites for {0} were supported".format(
                         ", ".join(call.mixed_lineages[:6]))),
            fail_status=STATUS_WARN,
        ),
    ]
    if call.is_bcg:
        checks.append(Check(
            name="bcg_intrinsic_pyrazinamide_resistance",
            value=True,
            threshold=True,
            source=source_for("bcg_pza_note"),
            status=STATUS_WARN,
            reading=BCG_PZA_NOTE,
            comparison="==",
            category="typing",
        ))
    if call.animal_variant == "bovis":
        checks.append(Check(
            name="mbovis_defining_snp_count",
            value=call.barcode_sites_total,
            threshold=MBOVIS_DEFINING_SNPS,
            source=source_for("mbovis_defining_snps"),
            status=STATUS_WARN,
            reading=MBOVIS_CAVEAT,
            comparison=">=",
            unit="sites",
            category="typing",
        ))
    return checks


def describe_call(call: LineageCall,
                  taxa: Optional[Mapping[str, TaxonSupport]] = None) -> str:
    """One line for the report: label, named family and the support behind it."""
    label = call.display
    if not label or label == "not determined":
        return "lineage not determined ({0} of {1} defining sites callable)".format(
            call.barcode_sites_callable, call.barcode_sites_total)
    family = ""
    if taxa is not None:
        family = family_for(taxa, call.sublineage or call.lineage)
    parts = [label]
    if family:
        parts.append("({0})".format(family))
    parts.append("- {0}/{1} defining sites supported, {2}/{3} callable".format(
        call.barcode_sites_supporting, call.barcode_sites_callable,
        call.barcode_sites_callable, call.barcode_sites_total))
    if call.is_bcg:
        parts.append("- BCG")
    elif call.is_animal and call.animal_variant:
        parts.append("- animal lineage: {0}".format(call.animal_variant))
    return " ".join(parts)


def barcode_taxa(sites: Sequence[BarcodeSite]) -> List[str]:
    """Every taxon the loaded scheme can call, in natural order."""
    return sorted({site.taxon for site in sites}, key=natural_key)


def sites_for_taxon(sites: Sequence[BarcodeSite], taxon: str) -> List[BarcodeSite]:
    """The defining sites for one taxon, for the annex and for the doctor report."""
    return [site for site in sites if site.taxon == taxon]


def iter_chain_sites(sites: Sequence[BarcodeSite],
                     taxon: str) -> Iterable[BarcodeSite]:
    """Defining sites for a taxon and every ancestor of it."""
    chain = set(taxon_ancestors(taxon) + [taxon])
    for site in sites:
        if site.taxon in chain:
            yield site
