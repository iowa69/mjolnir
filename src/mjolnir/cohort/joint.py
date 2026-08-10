"""The joint variant table: every sample's allele at every variable position.

This is MTBseq's ``TBjoin`` equivalent, written against the way that stage
actually failed on this machine. The 2022 *M. chimaera* run
(``M_Chimaera_TN/samp/MTBseq_2022-11-11_iowa.log``) ended with::

    <ERROR> No joint variant file try_joint_cf5_cr5_fr75_ph4_samples2.tab to amend!

Three separate assumptions produced that line: that a cohort is large, that the
reference is H37Rv, and that the joint file will be found by reconstructing its
name from the run's parameters. None of them hold here, so none of them is made.
A cohort of one and a cohort of two are ordinary inputs that produce a table and
a stated consequence, not an error; the reference is whatever the samples were
called against and is carried as data; and nothing in this module reaches for a
file by guessing what it was called.

**No cell is filled in by assumption.** The trap in every joint table is the
empty cell: a sample with no variant at a position is written as reference, and
a position with no coverage becomes evidence of identity. Here a cell is the
reference allele only where that sample's callable regions say the position was
callable. Everywhere else it is ``None`` — unknown — and
:mod:`mjolnir.cohort.distance` removes those positions from both the numerator
and the denominator rather than counting them as agreement. A sample supplied
without callable regions therefore has no known cells at all, which is the
honest reading of "we were not told what could be seen", and the checks say so
in words.

**A minority allele is not a genotype.** Where a sample's variant is below the
major-variant fraction the cell is :data:`AMBIGUOUS`, not the alternative
allele and not the reference: a 30% alternative allele is neither, and
collapsing it either way manufactures a difference or hides one.

Coordinates: :class:`~mjolnir.records.Variant` positions are 1-based, as in VCF.
:class:`Regions` stores 0-based half-open intervals, as in BED, and converts at
its own boundary — so BED files load unchanged and variant positions are asked
about unchanged.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from pathlib import Path
from typing import (Any, Dict, Iterable, Iterator, List, Mapping, Optional,
                    Sequence, Set, Tuple)

from ..config import MAJOR_VARIANT_FRACTION, is_major_variant
from ..records import (Check, STATUS_PASS, STATUS_WARN, Variant, VARIANT_SNP)
from ..utils import LOG, MjolnirError, PathLike, require_file, safe_fraction

#: What a cell holds when the sample was callable at the position but no single
#: allele holds the majority — a heterozygous or minor-variant call. It is not a
#: base and must never be compared as one.
AMBIGUOUS = "N"

#: How the TSV renders a cell that is not known. Deliberately not "-" and not
#: the reference base: a reader scanning the column has to be unable to mistake
#: "not callable here" for "same as the reference".
UNKNOWN_SYMBOL = "?"

#: VCF FILTER values that mean "this call passed". Anything else in
#: :attr:`Variant.filters` keeps the variant out of the joint table, and the
#: table counts how many were dropped that way.
PASSING_FILTERS = ("", "PASS", ".")


# ---------------------------------------------------------------------------
# Interval geometry
# ---------------------------------------------------------------------------

class Regions(object):
    """A set of genomic intervals, per contig, in BED coordinates.

    Callable regions, masks and their intersections are all the same shape, so
    they are all this class: the denominator of a pairwise distance is
    ``callable_a.intersect(callable_b).subtract(mask).length()``, computed on
    intervals rather than on a per-base structure, because a 4.4 Mb genome times
    a 159-sample cohort is not a place to allocate one Python object per base.

    Intervals are 0-based half-open ``[start, end)`` exactly as BED writes them.
    :meth:`contains` takes a **1-based** position, because that is what a VCF
    record and a :class:`~mjolnir.records.Variant` carry, and doing the
    conversion here is what keeps it from being done wrongly at nine call sites.
    """

    __slots__ = ("name", "_raw", "_merged", "_index")

    def __init__(self, name: str = "") -> None:
        self.name = name
        self._raw: Dict[str, List[Tuple[int, int]]] = {}
        self._merged: Optional[Dict[str, List[Tuple[int, int]]]] = None
        self._index: Dict[str, Tuple[List[int], List[int]]] = {}

    # -- construction -------------------------------------------------------

    def add(self, chrom: str, start: int, end: int) -> "Regions":
        """Add a 0-based half-open interval. Empty and reversed intervals raise."""
        start, end = int(start), int(end)
        if end <= start:
            raise MjolnirError(
                "interval {0}:{1}-{2} is empty or reversed; BED intervals are "
                "0-based half-open and must satisfy start < end".format(chrom, start, end)
            )
        if start < 0:
            raise MjolnirError(
                "interval {0}:{1}-{2} starts before the contig".format(chrom, start, end))
        self._raw.setdefault(str(chrom), []).append((start, end))
        self._merged = None
        self._index = {}
        return self

    def add_1based(self, chrom: str, start: int, end: int) -> "Regions":
        """Add a 1-based inclusive interval, as a GFF or a pileup would state it."""
        return self.add(chrom, int(start) - 1, int(end))

    def add_position(self, chrom: str, pos: int) -> "Regions":
        """Add a single 1-based position."""
        return self.add(chrom, int(pos) - 1, int(pos))

    @classmethod
    def whole(cls, lengths: Mapping[str, int], name: str = "") -> "Regions":
        """Every base of every contig — the callable set of a complete assembly."""
        regions = cls(name=name)
        for chrom, length in lengths.items():
            if int(length) > 0:
                regions.add(chrom, 0, int(length))
        return regions

    @classmethod
    def from_bed(cls, path: PathLike, name: str = "",
                 fetch_hint: str = "") -> "Regions":
        """Load a BED file: ``chrom``, ``start``, ``end`` and anything after them.

        Malformed lines raise rather than being skipped. A mask file that
        silently loses half its rows to a stray header would produce distances
        that look fine and are not, which is the exact failure this project is
        written against.
        """
        resolved = require_file(path, "BED file", fetch_hint)
        regions = cls(name=name or resolved.name)
        kept = 0
        with resolved.open("rt") as handle:
            for number, raw in enumerate(handle, start=1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                if line.startswith("track") or line.startswith("browser"):
                    continue
                fields = line.split("\t") if "\t" in line else line.split()
                if len(fields) < 3:
                    raise MjolnirError(
                        "{0} line {1}: expected at least 3 tab-separated columns "
                        "(chrom, start, end), got {2}".format(resolved, number, len(fields))
                    )
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    raise MjolnirError(
                        "{0} line {1}: start and end must be integers, got {2!r} and "
                        "{3!r}. BED has no header line; if this file has one, remove "
                        "it or comment it out with '#'.".format(
                            resolved, number, fields[1], fields[2])
                    )
                if end <= start:
                    raise MjolnirError(
                        "{0} line {1}: interval {2}:{3}-{4} is empty or reversed; BED "
                        "is 0-based half-open".format(resolved, number, fields[0], start, end)
                    )
                regions.add(fields[0], start, end)
                kept += 1
        if not kept:
            raise MjolnirError(
                "{0} contains no intervals; an empty mask or callable-region file is "
                "not the same as no masking, and Mjolnir will not treat it as "
                "either".format(resolved)
            )
        LOG.debug("loaded %d intervals from %s", kept, resolved)
        return regions

    # -- access -------------------------------------------------------------

    @property
    def intervals(self) -> Dict[str, List[Tuple[int, int]]]:
        """Merged, sorted intervals per contig."""
        if self._merged is None:
            merged: Dict[str, List[Tuple[int, int]]] = {}
            for chrom, spans in self._raw.items():
                out: List[Tuple[int, int]] = []
                for start, end in sorted(spans):
                    if out and start <= out[-1][1]:
                        if end > out[-1][1]:
                            out[-1] = (out[-1][0], end)
                    else:
                        out.append((start, end))
                merged[chrom] = out
            self._merged = merged
            self._index = {}
        return self._merged

    def chroms(self) -> List[str]:
        return sorted(self.intervals)

    def contains(self, chrom: str, pos: int) -> bool:
        """Whether a **1-based** position falls inside the set."""
        spans = self.intervals.get(str(chrom))
        if not spans:
            return False
        starts, ends = self._index.get(str(chrom), ([], []))
        if not starts:
            starts = [s for s, _ in spans]
            ends = [e for _, e in spans]
            self._index[str(chrom)] = (starts, ends)
        offset = int(pos) - 1
        idx = bisect.bisect_right(starts, offset) - 1
        return idx >= 0 and offset < ends[idx]

    def length(self, chrom: Optional[str] = None) -> int:
        """Total bases covered, over one contig or all of them."""
        spans = self.intervals
        if chrom is not None:
            return sum(end - start for start, end in spans.get(str(chrom), []))
        return sum(end - start for chrom_spans in spans.values() for start, end in chrom_spans)

    def is_empty(self) -> bool:
        return self.length() == 0

    def __bool__(self) -> bool:
        return not self.is_empty()

    def __len__(self) -> int:
        return self.length()

    def __repr__(self) -> str:
        return "Regions({0!r}, {1} contigs, {2} bp)".format(
            self.name, len(self.intervals), self.length())

    # -- set algebra --------------------------------------------------------

    def intersect(self, other: "Regions", name: str = "") -> "Regions":
        """The intervals present in both sets."""
        out = Regions(name=name or "{0} & {1}".format(self.name, other.name))
        mine, theirs = self.intervals, other.intervals
        for chrom in set(mine) & set(theirs):
            a, b = mine[chrom], theirs[chrom]
            i = j = 0
            while i < len(a) and j < len(b):
                start = max(a[i][0], b[j][0])
                end = min(a[i][1], b[j][1])
                if start < end:
                    out.add(chrom, start, end)
                if a[i][1] <= b[j][1]:
                    i += 1
                else:
                    j += 1
        return out

    def subtract(self, other: "Regions", name: str = "") -> "Regions":
        """This set with *other* removed."""
        out = Regions(name=name or "{0} - {1}".format(self.name, other.name))
        theirs = other.intervals
        for chrom, spans in self.intervals.items():
            cuts = theirs.get(chrom, [])
            cut_starts = [s for s, _ in cuts]
            for start, end in spans:
                cursor = start
                idx = max(bisect.bisect_right(cut_starts, start) - 1, 0)
                for cut_start, cut_end in cuts[idx:]:
                    if cut_start >= end:
                        break
                    if cut_end <= cursor:
                        continue
                    if cut_start > cursor:
                        out.add(chrom, cursor, min(cut_start, end))
                    cursor = max(cursor, cut_end)
                    if cursor >= end:
                        break
                if cursor < end:
                    out.add(chrom, cursor, end)
        return out


# ---------------------------------------------------------------------------
# Inputs to the join
# ---------------------------------------------------------------------------

@dataclass
class SampleVariants:
    """One sample's calls, and what the sample could be called at.

    ``callable_regions`` is not optional in spirit even though it is typed
    ``Optional``: without it, nothing about this sample can be compared to
    anything, because every position that is not a variant is unknown rather
    than reference. The pipeline fills it from the depth engine (positions at or
    above the degraded depth floor); an assembly fills it with
    :meth:`Regions.whole` over the contig lengths that aligned.
    """

    sample_id: str
    variants: Sequence[Variant] = field(default_factory=list)
    callable_regions: Optional[Regions] = None
    reference: str = ""
    platform: str = ""
    note: str = ""

    def __post_init__(self) -> None:
        if not self.sample_id:
            raise MjolnirError("SampleVariants requires a sample id")

    @property
    def has_callable_regions(self) -> bool:
        return self.callable_regions is not None and bool(self.callable_regions.length())


# ---------------------------------------------------------------------------
# The table
# ---------------------------------------------------------------------------

@dataclass
class JointSite:
    """One position at which at least one sample carries a non-reference allele.

    ``alleles`` holds only the samples that departed from the reference — the
    non-reference base, or :data:`AMBIGUOUS` where no allele held the majority.
    Every other sample's cell is derived from its callable regions when the
    table is read, so a sample's absence from this dict never by itself means
    "reference".
    """

    chrom: str
    pos: int
    ref: str
    alleles: Dict[str, str] = field(default_factory=dict)
    #: A substitution of one base for one base in every sample that carries it.
    #: Indels and MNVs are kept in the table and excluded from SNP counting.
    is_snp: bool = True
    variant_type: str = VARIANT_SNP
    genes: List[str] = field(default_factory=list)
    note: str = ""

    @property
    def key(self) -> Tuple[str, int]:
        return (self.chrom, self.pos)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "is_snp": self.is_snp,
            "variant_type": self.variant_type,
            "genes": list(self.genes),
            "alleles": dict(self.alleles),
            "note": self.note,
        }


@dataclass
class JointTable:
    """Every variable position across a cohort, with each sample's allele.

    Read it through :meth:`allele`, never through ``site.alleles`` directly:
    the dict holds departures from the reference, and :meth:`allele` is what
    turns "absent from the dict" into either the reference allele or ``None``
    by consulting the sample's callable regions.
    """

    samples: List[str] = field(default_factory=list)
    sites: List[JointSite] = field(default_factory=list)
    reference: str = ""
    callable_regions: Dict[str, Regions] = field(default_factory=dict)
    checks: List[Check] = field(default_factory=list)
    caveats: List[str] = field(default_factory=list)
    #: Variants left out of the table, by reason, for the annex.
    excluded: Dict[str, int] = field(default_factory=dict)
    #: Per-sample sorted variable positions, built once at construction.
    _positions: Dict[str, List[Tuple[str, int]]] = field(default_factory=dict, repr=False)
    #: Position-keyed view of :attr:`sites`, built on first use.
    _index: Dict[Tuple[str, int], JointSite] = field(default_factory=dict, repr=False)

    # -- shape --------------------------------------------------------------

    @property
    def site_count(self) -> int:
        return len(self.sites)

    @property
    def snp_site_count(self) -> int:
        return sum(1 for site in self.sites if site.is_snp)

    def samples_without_callable_regions(self) -> List[str]:
        return [s for s in self.samples
                if s not in self.callable_regions or not self.callable_regions[s].length()]

    # -- reading ------------------------------------------------------------

    def allele(self, sample: str, site: JointSite) -> Optional[str]:
        """This sample's allele at this site, or ``None`` when it is not known.

        ``None`` is returned for a position outside the sample's callable
        regions and for a sample with no callable regions at all. It is not the
        reference allele: an uncovered position is not evidence of agreement,
        and treating it as one is how a coverage gap becomes a transmission
        cluster.
        """
        found = site.alleles.get(sample)
        if found is not None:
            return found
        regions = self.callable_regions.get(sample)
        if regions is None:
            return None
        return site.ref if regions.contains(site.chrom, site.pos) else None

    def positions_for(self, sample: str) -> List[Tuple[str, int]]:
        """Sorted ``(chrom, pos)`` at which this sample departs from the reference."""
        return self._positions.get(sample, [])

    def site_index(self) -> Dict[Tuple[str, int], JointSite]:
        """``(chrom, pos)`` to site, built once and kept.

        The distance code walks pair after pair over the same sites, so the
        lookup is built here rather than rebuilt per pair — and it is rebuilt if
        the site list has changed length under it, so a table amended after a
        first read cannot serve a stale index.
        """
        if len(self._index) != len(self.sites):
            self._index = dict((site.key, site) for site in self.sites)
        return self._index

    def site_at(self, chrom: str, pos: int) -> Optional[JointSite]:
        return self.site_index().get((str(chrom), int(pos)))

    def rows(self) -> Iterator[Dict[str, Any]]:
        """One dict per site, sample columns filled, for TSV or a DataFrame."""
        for site in self.sites:
            row: Dict[str, Any] = {
                "chrom": site.chrom,
                "pos": site.pos,
                "ref": site.ref,
                "type": site.variant_type,
                "gene": ",".join(site.genes),
            }
            for sample in self.samples:
                allele = self.allele(sample, site)
                row[sample] = UNKNOWN_SYMBOL if allele is None else allele
            yield row

    def write_tsv(self, path: PathLike, comments: bool = True) -> Path:
        """Write the table. ``?`` is not callable; ``N`` is called but ambiguous.

        The legend is written into the file as ``#`` comment lines by default,
        because this file outlives the run that produced it and a bare ``?`` in
        a column of bases is exactly the kind of symbol somebody guesses at.
        """
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        columns = ["chrom", "pos", "ref", "type", "gene"] + list(self.samples)
        with target.open("wt") as handle:
            if comments:
                handle.write("# mjolnir joint variant table\n")
                handle.write("# reference: {0}\n".format(self.reference or "unstated"))
                handle.write("# samples: {0}\n".format(len(self.samples)))
                handle.write("# variable positions: {0} ({1} SNP)\n".format(
                    self.site_count, self.snp_site_count))
                handle.write("# {0} = position not callable in that sample "
                             "(absence of evidence, not the reference allele)\n"
                             .format(UNKNOWN_SYMBOL))
                handle.write("# {0} = callable but no allele reached the "
                             "major-variant fraction of {1}\n".format(
                                 AMBIGUOUS, MAJOR_VARIANT_FRACTION))
            handle.write("\t".join(columns) + "\n")
            for row in self.rows():
                handle.write("\t".join(str(row[column]) for column in columns) + "\n")
        return target

    def to_dict(self, include_sites: bool = False) -> Dict[str, Any]:
        data: Dict[str, Any] = {
            "reference": self.reference,
            "samples": list(self.samples),
            "n_samples": len(self.samples),
            "joint_sites": self.site_count,
            "snp_sites": self.snp_site_count,
            "callable_bases": dict(
                (sample, regions.length()) for sample, regions in self.callable_regions.items()),
            "samples_without_callable_regions": self.samples_without_callable_regions(),
            "excluded": dict(self.excluded),
            "checks": [c.to_dict() for c in self.checks],
            "caveats": list(self.caveats),
        }
        if include_sites:
            data["sites"] = [site.to_dict() for site in self.sites]
        return data


# ---------------------------------------------------------------------------
# Building it
# ---------------------------------------------------------------------------

def cohort_size_check(n_samples: int) -> Check:
    """The check that stands where MTBseq's ``TBjoin`` aborted.

    A cohort of one or two is a legitimate thing to ask for — two isolates from
    one patient, a sample against its own repeat — and MTBseq treats it as a
    fatal condition, reconstructing a filename such as
    ``try_joint_cf5_cr5_fr75_ph4_samples2.tab`` and stopping when it is not
    there. Mjolnir states the consequence instead: with one sample there is
    nothing to compare it to, with two there is exactly one comparison, and
    neither is an error.
    """
    if n_samples <= 0:
        raise MjolnirError("a cohort needs at least one sample")
    if n_samples == 1:
        return Check(
            name="cohort_size", value=1, threshold=2, comparison=">=",
            status=STATUS_WARN, unit="samples", category="cohort",
            reading="single-sample cohort: the joint variant table was built, but a "
                    "distance is a statement about a pair and there is no pair here. "
                    "No distances and no clusters follow, which is an absence of "
                    "comparison rather than a finding of no relatedness.")
    if n_samples == 2:
        return Check(
            name="cohort_size", value=2, threshold=2, comparison=">=",
            status=STATUS_PASS, unit="samples", category="cohort",
            reading="two-sample cohort: exactly one pairwise distance. MTBseq's TBjoin "
                    "aborts on cohorts this small; Mjolnir reports the single distance "
                    "with its shared-callable denominator.")
    return Check(
        name="cohort_size", value=n_samples, threshold=2, comparison=">=",
        status=STATUS_PASS, unit="samples", category="cohort",
        reading="{0} samples: {1} pairwise comparisons.".format(
            n_samples, n_samples * (n_samples - 1) // 2))


def _passes_filters(variant: Variant) -> bool:
    for entry in variant.filters or []:
        if str(entry).strip().upper() not in tuple(f.upper() for f in PASSING_FILTERS):
            return False
    return True


def _allele_for(variant: Variant) -> Optional[str]:
    """The genotype this variant implies, or ``AMBIGUOUS``, or None to skip it.

    ``is_major`` decides, falling back to the allele fraction when the caller
    did not set it. When neither exists the input is an assembly, where there is
    no allele fraction to have — the consensus base is the only genotype
    available and it is used, with the capability loss recorded on the table
    rather than pretended away.
    """
    major = variant.is_major
    if major is None:
        major = is_major_variant(variant.allele_fraction)
    if major is False:
        return AMBIGUOUS
    return variant.alt.upper()


def build_joint_table(samples: Sequence[SampleVariants], reference: str = "",
                      include_indels: bool = True) -> JointTable:
    """Join per-sample calls into one table of variable positions.

    Works against any reference. The reference is taken from the samples
    themselves when they agree and raises when they do not, because a joint
    table across two coordinate systems is not a table of anything — this is the
    other half of what broke the *M. chimaera* run, where every MTB-specific
    input was ``NONE`` and the stage carried on regardless.

    Indels are kept by default. They are excluded from SNP distances by
    :mod:`mjolnir.cohort.distance`, but they are real differences and dropping
    them from the table would also drop them from the proximity rule, which
    exists precisely because variants cluster around alignment trouble.
    """
    if not samples:
        raise MjolnirError(
            "cohort mode needs at least one sample; nothing was supplied to join")

    names: List[str] = []
    for entry in samples:
        if entry.sample_id in names:
            raise MjolnirError(
                "sample {0!r} appears twice in the cohort; sample ids must be "
                "unique or the joint table columns are ambiguous".format(entry.sample_id))
        names.append(entry.sample_id)

    references = set(e.reference for e in samples if e.reference)
    if len(references) > 1:
        raise MjolnirError(
            "cohort samples were called against different references ({0}); a joint "
            "variant table needs one coordinate system. Re-run the differing samples "
            "against the same reference.".format(", ".join(sorted(references)))
        )
    resolved_reference = reference or (sorted(references)[0] if references else "")

    sites: Dict[Tuple[str, int], JointSite] = {}
    positions: Dict[str, List[Tuple[str, int]]] = dict((name, []) for name in names)
    excluded: Dict[str, int] = {"filtered": 0, "no_alt": 0, "indel_excluded": 0}
    fasta_like: List[str] = []

    for entry in samples:
        seen: Set[Tuple[str, int]] = set()
        for variant in entry.variants:
            if not _passes_filters(variant):
                excluded["filtered"] += 1
                continue
            if not variant.alt or variant.alt == variant.ref:
                excluded["no_alt"] += 1
                continue
            is_snp = (len(variant.ref) == 1 and len(variant.alt) == 1
                      and variant.variant_type == VARIANT_SNP)
            if not is_snp and not include_indels:
                excluded["indel_excluded"] += 1
                continue
            if variant.is_major is None and variant.allele_fraction is None:
                if entry.sample_id not in fasta_like:
                    fasta_like.append(entry.sample_id)

            key = (variant.chrom, int(variant.pos))
            site = sites.get(key)
            if site is None:
                site = JointSite(chrom=variant.chrom, pos=int(variant.pos),
                                 ref=variant.ref.upper(), is_snp=is_snp,
                                 variant_type=variant.variant_type)
                sites[key] = site
            if not is_snp:
                # One sample's indel makes the position an indel position for the
                # whole cohort: it is not a base-for-base comparison there.
                site.is_snp = False
                site.variant_type = variant.variant_type
            if variant.gene and variant.gene not in site.genes:
                site.genes.append(variant.gene)

            allele = _allele_for(variant)
            if allele is None:
                continue
            if key in seen and site.alleles.get(entry.sample_id) not in (None, allele):
                # Two different alternative alleles at one position in one sample:
                # neither is the genotype, and picking one would invent a call.
                site.alleles[entry.sample_id] = AMBIGUOUS
                site.note = ("more than one alternative allele was called here in at "
                             "least one sample")
                continue
            site.alleles[entry.sample_id] = allele
            seen.add(key)

    ordered = sorted(sites.values(), key=lambda s: (s.chrom, s.pos))
    for site in ordered:
        for sample in site.alleles:
            positions[sample].append((site.chrom, site.pos))

    callable_regions = dict(
        (e.sample_id, e.callable_regions) for e in samples if e.callable_regions is not None)

    table = JointTable(samples=names, sites=ordered, reference=resolved_reference,
                       callable_regions=callable_regions, excluded=excluded)
    table._positions = dict((name, sorted(pos)) for name, pos in positions.items())

    table.checks.append(cohort_size_check(len(names)))
    table.checks.append(Check(
        name="joint_variant_sites", value=table.site_count, status=STATUS_PASS,
        unit="positions", category="cohort",
        reading="{0} positions vary across the cohort, of which {1} are single-base "
                "substitutions.".format(table.site_count, table.snp_site_count)))

    missing = table.samples_without_callable_regions()
    if missing:
        table.checks.append(Check.not_measured(
            "callable_regions",
            "no callable regions were supplied for {0}, so every position where "
            "these samples have no variant is unknown rather than reference. "
            "Distances involving them cannot be computed and are reported as "
            "not computed.".format(", ".join(missing)),
            category="cohort"))
        table.caveats.append(
            "callable regions are missing for {0} of {1} samples; a joint table cannot "
            "distinguish 'same as the reference' from 'not covered' without them"
            .format(len(missing), len(names)))
    else:
        table.checks.append(Check(
            name="callable_regions", value=len(names), threshold=len(names),
            comparison="==", status=STATUS_PASS, unit="samples", category="cohort",
            reading="callable regions are known for every sample, so an empty cell is "
                    "reported as not-callable rather than as agreement with the "
                    "reference."))

    if fasta_like:
        table.caveats.append(
            "{0} carried no allele fractions (assembly input), so minority alleles "
            "could not be distinguished from fixed ones at any position"
            .format(", ".join(sorted(fasta_like))))

    if not resolved_reference:
        table.caveats.append(
            "no reference was named for this cohort; the coordinates in this table "
            "are only meaningful against whatever the samples were called on")

    LOG.debug("joint table: %d samples, %d variable positions (%d SNP)",
              len(names), table.site_count, table.snp_site_count)
    return table


def callable_summary(table: JointTable) -> Dict[str, Optional[float]]:
    """Callable bases per sample as a fraction of the largest callable set.

    A cheap cohort-level look at whether one sample is going to drag every
    distance it takes part in: the denominator of a pair is bounded by the
    smaller of the two, so a sample at 40% of the cohort maximum has already
    decided that its distances are not comparable to the published SNP
    thresholds, before anything is counted.
    """
    lengths = dict((sample, regions.length())
                   for sample, regions in table.callable_regions.items())
    if not lengths:
        return {}
    largest = max(lengths.values())
    return dict((sample, safe_fraction(value, largest)) for sample, value in lengths.items())


def merged_positions(table: JointTable, samples: Iterable[str]) -> List[Tuple[str, int]]:
    """Sorted union of the variable positions of the named samples.

    Used by the distance code to walk only the positions that can possibly
    differ between two genomes, rather than every position in the cohort.
    """
    union: Set[Tuple[str, int]] = set()
    for sample in samples:
        union.update(table.positions_for(sample))
    return sorted(union)
