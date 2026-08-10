"""Mycobacterial species identification, by ANI and by nothing else.

The central prohibition of design §6 is enforced structurally in this module
rather than asked for in a comment. Species identification for this genus must
not come from a taxonomic read classifier, because in current NCBI taxonomy the
MTBC members are not at species rank at all: ``Mycobacterium tuberculosis
variant bovis`` (taxid 1765) has rank ``no rank`` beneath species
*M. tuberculosis*, the MTBC members being later heterotypic synonyms of
*M. tuberculosis* at 99.21-99.92% ANI. A Kraken2 row reading "M. bovis 3.2%" is
therefore not a species identification, and a tool that prints one has invented
a result its evidence cannot support.

Three things make that impossible here rather than merely discouraged:

* The only entry point that produces a :class:`~mjolnir.records.SpeciesCall` is
  :func:`identify_species`, and its evidence argument is a query genome or read
  set plus an ANI reference set. There is no parameter that accepts classifier
  output, so no caller can feed one in.
* :func:`species_from_classifier` exists and always raises. It is the name a
  future caller would reach for, and reaching for it yields the taxonomy
  explanation instead of a call.
* Every call is built through :func:`_species_call`, which refuses a method
  string outside :data:`ANI_METHODS` and refuses to mark any MTBC member name as
  a resolved species. Even a reference manifest that lists *M. bovis* as a
  reference genome cannot produce "M. bovis" as an identification — the call is
  demoted to the complex, with the reason attached.

What the module *will* do is name a species where ANI can carry the claim, and
say "cannot resolve below complex" where it cannot. That outcome is not a
failure: inside the MTBC it is the only honest answer, and inside MAC the
*M. chimaera* / *M. intracellulare* / *M. avium* boundary needs marker SNPs
beside the ANI, which is exactly the distinction the local outbreak data
requires.

Process execution lives in ``shell.py``; this module builds commands and parses
their output, and takes an injectable *runner* so the parsers can be unit-tested
without a mycobacterial genome on disk.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..config import (
    ANI_MIN_ALIGNED_FRACTION,
    ANI_SPECIES_FLOOR,
    COMPLEX_ABSCESSUS,
    COMPLEX_MAC,
    COMPLEX_MTBC,
    MAC_SPECIES_ANI_FLOOR,
    MTBC_INTRA_ANI_RANGE,
    MTBC_UNRESOLVED_TEXT,
    SPECIES_METHOD_REFUSAL,
    default_db_dir,
    source_for,
)
from ..records import (
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MODERATE,
    CONFIDENCE_NONE,
    STATUS_FAIL,
    STATUS_PASS,
    STATUS_WARN,
    Check,
    SpeciesCall,
)
from ..utils import (
    LOG,
    MjolnirError,
    PathLike,
    first_available,
    require,
    require_database,
    require_file,
)
from .lineage import PileupCounts, major_allele

# ---------------------------------------------------------------------------
# Where the reference set lives, and what it is
# ---------------------------------------------------------------------------

#: Layout under the database root, written by the operator (see below);
#: nothing here is vendored, because the reference set is the largest single
#: download Mjolnir needs and its composition is still an open question (design
#: §14).
ANI_DIR_NAME = "ani"
REFERENCE_MANIFEST_NAME = "references.tsv"
MASH_SKETCH_NAME = "mycobacteria.msh"
MAC_MARKER_NAME = "mac_markers.tsv"

#: Mjolnir ships no ANI reference set — design §14 leaves its composition and
#: download budget open — so the hint says how to build one rather than naming a
#: `mjolnir db fetch` target that does not exist in the registry.
ANI_FETCH_HINT = (
    "put the reference genomes under <db>/ani/ and list them in "
    "<db>/ani/references.tsv (columns: file, name, and optionally complex, "
    "accession, subspecies, note, source)"
)

#: Columns the manifest must carry, and the ones it may. ``complex`` is optional
#: because :func:`complex_for` can derive it from the name, but a manifest that
#: states it wins — a curated file knows things a name table does not.
MANIFEST_REQUIRED_COLUMNS: Tuple[str, ...] = ("file", "name")
MANIFEST_OPTIONAL_COLUMNS: Tuple[str, ...] = (
    "complex", "accession", "subspecies", "note", "source",
)

#: The method strings a species call may carry. Anything else is refused by
#: :func:`_species_call`, which is what stops a classifier result from ever
#: being dressed up as an identification.
METHOD_SKANI = "skani ANI against the curated mycobacterial reference set"
METHOD_MASH = "mash distance against the curated mycobacterial reference set"
ANI_METHODS: Tuple[str, ...] = (METHOD_SKANI, METHOD_MASH)

#: Mash reports a distance, not an ANI, and the conversion below is an estimate
#: from the mutation model in Ondov et al. 2016 (Genome Biol 17:132), not a
#: measured alignment. It is good to roughly a tenth of a percent near the
#: species boundary and is not adequate on its own inside a complex — which is
#: why a mash-derived call always carries MASH_ESTIMATE_CAVEAT.
MASH_ESTIMATE_CAVEAT = (
    "ANI here is estimated from a mash k-mer distance rather than measured by "
    "alignment; it is adequate for a genus/species-level statement and is not "
    "used to resolve inside a complex"
)

#: skani is an assembly-to-assembly comparison and is not defined for raw reads,
#: so read input goes to mash. Stated rather than silently switched, because the
#: two methods do not produce the same number.
READS_METHOD_NOTE = (
    "skani compares assembled sequence and is not defined for raw reads, so the "
    "read-level ANI screen uses mash with single-copy k-mers discarded"
)


# ---------------------------------------------------------------------------
# Taxonomy tables
#
# These are names, not thresholds: which binomials belong to which complex. They
# exist so that a reference hit, or a classifier label somebody else's module
# obtained, can be placed in the complex whose members ANI cannot separate.
# ---------------------------------------------------------------------------

#: MTBC members, in every spelling the databases on hand use. Note that NCBI now
#: writes most of these as "Mycobacterium tuberculosis variant <x>"; both forms
#: are here because both are in current use.
MTBC_MEMBERS: Tuple[str, ...] = (
    "mycobacterium tuberculosis",
    "mycobacterium tuberculosis variant bovis",
    "mycobacterium tuberculosis variant bovis bcg",
    "mycobacterium tuberculosis variant africanum",
    "mycobacterium tuberculosis variant caprae",
    "mycobacterium tuberculosis variant microti",
    "mycobacterium tuberculosis variant pinnipedii",
    "mycobacterium tuberculosis variant orygis",
    "mycobacterium bovis",
    "mycobacterium bovis bcg",
    "mycobacterium africanum",
    "mycobacterium canettii",
    "mycobacterium caprae",
    "mycobacterium microti",
    "mycobacterium pinnipedii",
    "mycobacterium orygis",
    "mycobacterium mungi",
    "mycobacterium suricattae",
    "mycobacterium tuberculosis complex",
)

#: *Mycobacterium avium* complex. *M. chimaera* was split out of
#: *M. intracellulare* on ITS sequevar rather than on ANI, which is precisely why
#: ANI alone does not separate them here.
MAC_MEMBERS: Tuple[str, ...] = (
    "mycobacterium avium",
    "mycobacterium avium subsp. avium",
    "mycobacterium avium subsp. hominissuis",
    "mycobacterium avium subsp. paratuberculosis",
    "mycobacterium avium subsp. silvaticum",
    "mycobacterium intracellulare",
    "mycobacterium paraintracellulare",
    "mycobacterium chimaera",
    "mycobacterium colombiense",
    "mycobacterium arosiense",
    "mycobacterium marseillense",
    "mycobacterium timonense",
    "mycobacterium bouchedurhonense",
    "mycobacterium yongonense",
    "mycobacterium vulneris",
    "mycobacterium ituriense",
    "mycobacterium avium complex",
)

#: *Mycobacteroides abscessus* complex. The genus was split, and both the old and
#: the new binomial appear in the databases Mjolnir reads.
ABSCESSUS_MEMBERS: Tuple[str, ...] = (
    "mycobacteroides abscessus",
    "mycobacterium abscessus",
    "mycobacteroides abscessus subsp. abscessus",
    "mycobacteroides abscessus subsp. massiliense",
    "mycobacteroides abscessus subsp. bolletii",
    "mycobacterium abscessus subsp. abscessus",
    "mycobacterium abscessus subsp. massiliense",
    "mycobacterium abscessus subsp. bolletii",
    "mycobacterium massiliense",
    "mycobacterium bolletii",
)

#: The three MAC species the outbreak data needs separated. Kept apart from
#: :data:`MAC_MEMBERS` because these are the only ones the marker table is
#: expected to carry.
MAC_TARGET_SPECIES: Tuple[str, ...] = (
    "Mycobacterium chimaera",
    "Mycobacterium intracellulare",
    "Mycobacterium avium",
)

_COMPLEX_BY_NAME: Dict[str, str] = {}
for _name in MTBC_MEMBERS:
    _COMPLEX_BY_NAME[_name] = COMPLEX_MTBC
for _name in MAC_MEMBERS:
    _COMPLEX_BY_NAME[_name] = COMPLEX_MAC
for _name in ABSCESSUS_MEMBERS:
    _COMPLEX_BY_NAME[_name] = COMPLEX_ABSCESSUS

#: Long form of each complex, for prose. ``SpeciesCall.display`` returns the
#: short code in ``complex``; this is what the report prints beside it.
COMPLEX_LONG_NAME: Dict[str, str] = {
    COMPLEX_MTBC: "Mycobacterium tuberculosis complex",
    COMPLEX_MAC: "Mycobacterium avium complex",
    COMPLEX_ABSCESSUS: "Mycobacteroides abscessus complex",
}

#: What a classifier row about an MTBC or MAC member may honestly be turned into.
#: Used by the contamination screen, which does have a legitimate use for
#: classifier labels — as a description of what is in the library, never as an
#: identification of the isolate.
CLASSIFIER_DEMOTION: Dict[str, str] = {
    COMPLEX_MTBC: (
        "M. tuberculosis complex (member not resolvable by a read classifier)"),
    COMPLEX_MAC: (
        "M. avium complex (member not resolvable by a read classifier)"),
}


def normalise_species_name(name: str) -> str:
    """Lower-case, whitespace-collapsed form used for the taxonomy lookups."""
    return " ".join(str(name or "").split()).lower()


def complex_for(name: str) -> str:
    """The complex a binomial belongs to, or "" when it belongs to none.

    Matching walks from the full name down to the bare binomial, so
    "Mycobacterium avium subsp. hominissuis TH135" places as MAC through
    "mycobacterium avium" even though the strain suffix is not in the table.
    """
    key = normalise_species_name(name)
    if not key:
        return ""
    if key in _COMPLEX_BY_NAME:
        return _COMPLEX_BY_NAME[key]
    tokens = key.split()
    # Longest-prefix match: "mycobacterium avium subsp. avium ATCC 25291" ->
    # "mycobacterium avium subsp. avium" -> "mycobacterium avium".
    for stop in range(len(tokens), 1, -1):
        prefix = " ".join(tokens[:stop])
        if prefix in _COMPLEX_BY_NAME:
            return _COMPLEX_BY_NAME[prefix]
    return ""


def is_mtbc_member(name: str) -> bool:
    """Whether a name denotes an MTBC member, under any of its spellings."""
    return complex_for(name) == COMPLEX_MTBC


def demote_classifier_label(label: str) -> str:
    """What a taxonomic classifier row may be printed as.

    A classifier can legitimately say "this library contains reads that look
    mycobacterial and reads that look like *Cutibacterium*". It cannot say which
    MTBC member an isolate is, and it cannot separate *M. chimaera* from
    *M. intracellulare*. Rather than dropping those rows — which would hide real
    non-target signal — this replaces the member name with the complex and says
    the classifier could not resolve inside it.
    """
    demoted = CLASSIFIER_DEMOTION.get(complex_for(label))
    return demoted if demoted else str(label or "").strip()


def species_from_classifier(*args: Any, **kwargs: Any) -> SpeciesCall:
    """Always raises. There is no such thing, and this is where you find out why.

    This function exists so that the obvious call site — someone wiring Kraken2
    output into the typing layer — fails loudly with the taxonomy explanation
    instead of quietly producing a call. See :func:`demote_classifier_label` for
    what classifier output may legitimately be used for.
    """
    raise MjolnirError(
        "species identification from a taxonomic read classifier is refused. "
        + SPECIES_METHOD_REFUSAL
        + ". In current NCBI taxonomy the MTBC members are not at species rank: "
        "Mycobacterium tuberculosis variant bovis (taxid 1765) has rank 'no "
        "rank' under species M. tuberculosis, the members being later "
        "heterotypic synonyms at 99.21-99.92% ANI. Use identify_species() with "
        "an ANI reference set, and demote_classifier_label() if you need to "
        "report what a classifier saw in the library."
    )


# ---------------------------------------------------------------------------
# The reference set
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ReferenceGenome:
    """One genome in the curated ANI reference set, as the manifest describes it."""

    path: Path
    name: str
    complex: str = ""
    accession: str = ""
    subspecies: str = ""
    note: str = ""
    source: str = ""

    @property
    def exists(self) -> bool:
        return self.path.exists()

    @property
    def label(self) -> str:
        return self.accession or self.path.name


def ani_dir(db_dir: Optional[PathLike] = None) -> Path:
    """The ANI reference directory under the database root."""
    root = Path(db_dir) if db_dir is not None else default_db_dir()
    return Path(root) / ANI_DIR_NAME


def load_reference_set(db_dir: Optional[PathLike] = None) -> List[ReferenceGenome]:
    """Read ``<db>/ani/references.tsv``, or say exactly what to fetch.

    The manifest is authoritative about which genomes are in the set and what
    they are called, because the alternative — inferring a taxon from a filename
    — is how an assembly accession ends up printed as a species. Where the
    manifest does not state a complex, :func:`complex_for` supplies one; where it
    does, the manifest wins.
    """
    directory = ani_dir(db_dir)
    manifest = require_database(
        directory / REFERENCE_MANIFEST_NAME,
        "the mycobacterial ANI reference manifest",
        ANI_FETCH_HINT,
    )
    rows = _read_tsv(manifest, MANIFEST_REQUIRED_COLUMNS)
    references: List[ReferenceGenome] = []
    for line_number, row in rows:
        name = row["name"].strip()
        if not name:
            raise MjolnirError(
                "{0} line {1}: the 'name' column is empty; a reference with no "
                "taxon name cannot be reported as an identification".format(
                    manifest, line_number)
            )
        relative = row["file"].strip()
        if not relative:
            raise MjolnirError(
                "{0} line {1}: the 'file' column is empty".format(manifest, line_number))
        path = Path(relative)
        if not path.is_absolute():
            path = directory / path
        declared = row.get("complex", "").strip()
        references.append(ReferenceGenome(
            path=path,
            name=name,
            complex=declared or complex_for(name),
            accession=row.get("accession", "").strip(),
            subspecies=row.get("subspecies", "").strip(),
            note=row.get("note", "").strip(),
            source=row.get("source", "").strip(),
        ))
    if not references:
        raise MjolnirError(
            "{0} lists no reference genomes.\n  fetch it with: {1}".format(
                manifest, ANI_FETCH_HINT)
        )
    LOG.debug("loaded %d ANI references from %s", len(references), manifest)
    return references


def _read_tsv(path: Path, required: Sequence[str]) -> List[Tuple[int, Dict[str, str]]]:
    """Header-keyed TSV reader that names the file and line in every error.

    Deliberately not ``csv.DictReader``: a missing column there yields ``None``
    for that field and the failure surfaces three functions later as a
    ``TypeError`` about a NoneType, by which point the file that caused it is no
    longer in scope.
    """
    rows: List[Tuple[int, Dict[str, str]]] = []
    header: List[str] = []
    with open(str(path), "rt", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if not header:
                header = [f.strip().lstrip("#").strip().lower() for f in fields]
                missing = [c for c in required if c not in header]
                if missing:
                    raise MjolnirError(
                        "{0} is missing the required column(s) {1}; the header "
                        "found was {2}".format(
                            path, ", ".join(missing), ", ".join(header))
                    )
                continue
            if len(fields) < len(header):
                fields = fields + [""] * (len(header) - len(fields))
            rows.append((line_number, dict(zip(header, fields))))
    if not header:
        raise MjolnirError("{0} has no header row".format(path))
    return rows


# ---------------------------------------------------------------------------
# Running the ANI tools
# ---------------------------------------------------------------------------

#: A runner takes a command and returns its standard output. ``shell.py`` owns
#: process execution; this indirection exists so the parsers below can be tested
#: against captured tool output without either tool installed.
Runner = Callable[[Sequence[str]], str]


def _default_runner(command: Sequence[str]) -> str:
    """Run through ``mjolnir.shell``, and complain precisely if it cannot."""
    try:
        from ..shell import run  # noqa: WPS433 - shell.py owns process execution
    except ImportError as exc:
        raise MjolnirError(
            "mjolnir.shell could not be imported, so external tools cannot be "
            "run ({0}). Pass a runner= callable if you are driving the typing "
            "layer directly.".format(exc)
        )
    proc = run([str(part) for part in command])
    if isinstance(proc, str):
        return proc
    text = getattr(proc, "stdout", None)
    if text is None:
        raise MjolnirError(
            "mjolnir.shell.run returned {0!r}, which carries no stdout; the "
            "typing layer needs the captured output of {1}".format(
                type(proc).__name__, command[0])
        )
    return text


@dataclass
class AniMatch:
    """One query-to-reference comparison, with what the tool actually measured."""

    name: str
    complex: str
    reference: str
    ani: Optional[float] = None
    aligned_fraction: Optional[float] = None
    aligned_fraction_reference: Optional[float] = None
    #: mash only: the shared-hash fraction, which is evidence about the
    #: comparison but is not an aligned fraction and is not used as one.
    shared_hashes: Optional[float] = None
    method: str = ""
    subspecies: str = ""
    source: str = ""

    def as_candidate(self) -> Dict[str, Any]:
        """The annex row, so the reader can see the margin over the runner-up."""
        return {
            "name": self.name,
            "complex": self.complex,
            "reference": self.reference,
            "ani_percent": None if self.ani is None else round(self.ani, 4),
            "aligned_fraction": (None if self.aligned_fraction is None
                                 else round(self.aligned_fraction, 4)),
            "shared_hashes": (None if self.shared_hashes is None
                              else round(self.shared_hashes, 4)),
            "method": self.method,
        }


def _index_references(references: Sequence[ReferenceGenome]) -> Dict[str, ReferenceGenome]:
    """Look-up by every string a tool might echo back for a reference."""
    index: Dict[str, ReferenceGenome] = {}
    for reference in references:
        for key in (str(reference.path), reference.path.name, reference.name,
                    reference.accession):
            if key:
                index.setdefault(key, reference)
                index.setdefault(key.lower(), reference)
    return index


def _resolve_reference(token: str,
                       index: Mapping[str, ReferenceGenome]) -> Optional[ReferenceGenome]:
    """Match a tool's reference id back to a manifest entry, or give up.

    Giving up means returning None and the caller dropping the row with a debug
    line, which is right: a hit Mjolnir cannot name is a hit it cannot report,
    and inventing a name from the filename is the failure mode the manifest
    exists to prevent.
    """
    candidates = [token, token.lower(), Path(token).name, Path(token).name.lower()]
    for key in candidates:
        found = index.get(key)
        if found is not None:
            return found
    return None


def run_skani(query: PathLike, references: Sequence[ReferenceGenome], *,
              threads: int = 1, runner: Optional[Runner] = None) -> List[AniMatch]:
    """ANI by alignment, with skani, against every reference in the set.

    skani prints a header line and one row per pair that clears its own internal
    screening, so an absent row means "below skani's floor", not "identical".
    Columns are read by header name because skani has added columns between
    releases and positional parsing would silently shift.
    """
    # The PATH check belongs to the default runner. An injected runner is the
    # caller taking responsibility for execution — a test replaying captured
    # skani output, or a wrapper that runs it inside a container — and demanding
    # the binary anyway would only make the parser untestable.
    binary = require("skani", why="ANI-based species identification") \
        if runner is None else "skani"
    present = [r for r in references if r.exists]
    if not present:
        raise MjolnirError(
            "none of the {0} genomes in the ANI reference manifest are present "
            "on disk.\n  fetch them with: {1}".format(len(references), ANI_FETCH_HINT)
        )
    query_path = require_file(query, "the query genome for ANI identification")
    command: List[str] = [binary, "dist", "-t", str(max(1, int(threads))),
                          "-q", str(query_path), "-r"]
    command.extend(str(r.path) for r in present)
    text = (runner or _default_runner)(command)
    return _parse_skani(text, _index_references(present))


#: skani's column names. Read by name; a release that renames them fails loudly
#: here rather than producing ANI values taken from the wrong column.
SKANI_COLUMNS: Tuple[str, ...] = ("Ref_file", "ANI", "Align_fraction_query")


def _parse_skani(text: str, index: Mapping[str, ReferenceGenome]) -> List[AniMatch]:
    matches: List[AniMatch] = []
    header: List[str] = []
    for line in (text or "").splitlines():
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if not header:
            header = [f.strip() for f in fields]
            missing = [c for c in SKANI_COLUMNS if c not in header]
            if missing:
                raise MjolnirError(
                    "skani output is missing the column(s) {0}; Mjolnir reads "
                    "skani's columns by name and cannot guess. Header was: "
                    "{1}".format(", ".join(missing), ", ".join(header))
                )
            continue
        row = dict(zip(header, fields))
        reference = _resolve_reference(row.get("Ref_file", ""), index)
        if reference is None:
            LOG.debug("skani hit %r is not in the reference manifest; ignoring",
                      row.get("Ref_file", ""))
            continue
        matches.append(AniMatch(
            name=reference.name,
            complex=reference.complex or complex_for(reference.name),
            reference=reference.label,
            ani=_float_or_none(row.get("ANI")),
            aligned_fraction=_percent_to_fraction(row.get("Align_fraction_query")),
            aligned_fraction_reference=_percent_to_fraction(row.get("Align_fraction_ref")),
            method=METHOD_SKANI,
            subspecies=reference.subspecies,
            source=reference.source,
        ))
    matches.sort(key=lambda m: (-(m.ani or 0.0), m.name))
    return matches


def run_mash(query: PathLike, references: Sequence[ReferenceGenome], *,
             threads: int = 1, from_reads: bool = False,
             sketch: Optional[PathLike] = None,
             runner: Optional[Runner] = None) -> List[AniMatch]:
    """ANI estimated from a mash distance.

    Used for read input, where skani is not defined, and as the fallback when
    skani is not installed. ``-r -m 2`` on reads discards single-copy k-mers,
    which at ordinary depths are almost entirely sequencing error; without it a
    read set's distance to everything is inflated and nothing clears the floor.
    """
    binary = require("mash", why="ANI-based species identification") \
        if runner is None else "mash"
    query_path = require_file(query, "the query genome or read set for ANI identification")
    command: List[str] = [binary, "dist", "-p", str(max(1, int(threads)))]
    if from_reads:
        command += ["-r", "-m", "2"]
    if sketch is not None:
        command.append(str(require_file(
            sketch, "the mash sketch of the mycobacterial reference set", ANI_FETCH_HINT)))
        present: List[ReferenceGenome] = list(references)
    else:
        present = [r for r in references if r.exists]
        if not present:
            raise MjolnirError(
                "none of the {0} genomes in the ANI reference manifest are "
                "present on disk and no {1} sketch was found.\n  fetch them "
                "with: {2}".format(len(references), MASH_SKETCH_NAME, ANI_FETCH_HINT)
            )
        command.extend(str(r.path) for r in present)
    command.append(str(query_path))
    text = (runner or _default_runner)(command)
    return _parse_mash(text, _index_references(present))


def _parse_mash(text: str, index: Mapping[str, ReferenceGenome]) -> List[AniMatch]:
    """Parse ``mash dist``: reference, query, distance, p-value, shared hashes."""
    matches: List[AniMatch] = []
    for line in (text or "").splitlines():
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        reference = _resolve_reference(fields[0], index)
        if reference is None:
            LOG.debug("mash hit %r is not in the reference manifest; ignoring", fields[0])
            continue
        distance = _float_or_none(fields[2])
        if distance is None:
            continue
        matches.append(AniMatch(
            name=reference.name,
            complex=reference.complex or complex_for(reference.name),
            reference=reference.label,
            # Mash's own ANI estimate: 1 - distance, as a percentage.
            ani=(1.0 - distance) * 100.0,
            aligned_fraction=None,
            shared_hashes=_shared_hash_fraction(fields[4]) if len(fields) > 4 else None,
            method=METHOD_MASH,
            subspecies=reference.subspecies,
            source=reference.source,
        ))
    matches.sort(key=lambda m: (-(m.ani or 0.0), m.name))
    return matches


def _float_or_none(text: Optional[str]) -> Optional[float]:
    if text is None:
        return None
    stripped = str(text).strip()
    if not stripped:
        return None
    try:
        return float(stripped)
    except ValueError:
        return None


def _percent_to_fraction(text: Optional[str]) -> Optional[float]:
    """skani reports aligned fraction as a percentage; records store 0-1."""
    value = _float_or_none(text)
    return None if value is None else value / 100.0


def _shared_hash_fraction(text: str) -> Optional[float]:
    """``"421/1000"`` -> 0.421. Reported, never used as an aligned fraction."""
    parts = str(text).split("/")
    if len(parts) != 2:
        return None
    numerator = _float_or_none(parts[0])
    denominator = _float_or_none(parts[1])
    if numerator is None or not denominator:
        return None
    return numerator / denominator


def ani_matches(query: PathLike, references: Sequence[ReferenceGenome], *,
                is_reads: bool = False, threads: int = 1,
                db_dir: Optional[PathLike] = None,
                tool: Optional[str] = None,
                runner: Optional[Runner] = None) -> List[AniMatch]:
    """Run the best available ANI tool, or say which one to install.

    skani is preferred for assemblies because it measures ANI over an alignment
    and reports the aligned fraction, which is the number that distinguishes a
    species call from a 99.9% match over 2% of the genome. Reads go to mash
    regardless, because skani is not defined for them.
    """
    if tool is not None and tool not in ("skani", "mash"):
        raise MjolnirError(
            "unknown ANI tool {0!r}; Mjolnir drives skani or mash".format(tool))
    sketch_path = ani_dir(db_dir) / MASH_SKETCH_NAME
    sketch = sketch_path if sketch_path.exists() else None
    if is_reads:
        if tool == "skani":
            raise MjolnirError(READS_METHOD_NOTE)
        return run_mash(query, references, threads=threads, from_reads=True,
                        sketch=sketch, runner=runner)
    chosen = tool or first_available("skani", "mash")
    if chosen is None:
        raise MjolnirError(
            "no ANI tool found on PATH. Species identification for this genus is "
            "ANI-based and has no substitute — a taxonomic read classifier is "
            "not one.\n  conda install -c conda-forge -c bioconda skani\n"
            "  (or: conda install -c conda-forge -c bioconda mash)"
        )
    if chosen == "skani":
        return run_skani(query, references, threads=threads, runner=runner)
    return run_mash(query, references, threads=threads, from_reads=False,
                    sketch=sketch, runner=runner)


# ---------------------------------------------------------------------------
# MAC marker SNPs
# ---------------------------------------------------------------------------

#: Columns of ``<db>/ani/mac_markers.tsv``. Every row carries its own ``source``
#: because these positions are evidence, and a marker position with no citation
#: is exactly the bare magic number the house rules forbid. Mjolnir does not
#: hard-code them: the marker set is versioned with the database, appears in the
#: database registry, and a run without it says so rather than guessing.
MARKER_REQUIRED_COLUMNS: Tuple[str, ...] = ("chrom", "pos", "allele", "species")
MARKER_OPTIONAL_COLUMNS: Tuple[str, ...] = ("ref", "gene", "note", "source")


@dataclass(frozen=True)
class MarkerSnp:
    """One position whose allele separates two members of a complex."""

    chrom: str
    pos: int
    allele: str
    species: str
    ref: str = ""
    gene: str = ""
    note: str = ""
    source: str = ""

    @property
    def key(self) -> Tuple[str, int]:
        return (self.chrom, self.pos)


def load_mac_markers(db_dir: Optional[PathLike] = None,
                     required: bool = False) -> List[MarkerSnp]:
    """Marker SNPs separating *M. chimaera*, *M. intracellulare* and *M. avium*.

    Returns an empty list when the file is absent and *required* is False, and
    the caller then reports MAC at complex level with the reason attached. That
    is a stated capability loss, not a silent fallback: no species is named on
    ANI alone inside MAC, because ANI does not separate these three.
    """
    path = ani_dir(db_dir) / MAC_MARKER_NAME
    if not path.exists():
        if required:
            raise MjolnirError(
                "the MAC marker-SNP table {0} is absent, so M. chimaera cannot "
                "be separated from M. intracellulare.\n  fetch it with: "
                "{1}".format(path, ANI_FETCH_HINT)
            )
        LOG.debug("no MAC marker table at %s", path)
        return []
    markers: List[MarkerSnp] = []
    for line_number, row in _read_tsv(path, MARKER_REQUIRED_COLUMNS):
        position = _float_or_none(row["pos"])
        if position is None:
            raise MjolnirError(
                "{0} line {1}: 'pos' is not a number ({2!r})".format(
                    path, line_number, row["pos"])
            )
        allele = row["allele"].strip().upper()
        if not allele:
            raise MjolnirError(
                "{0} line {1}: 'allele' is empty".format(path, line_number))
        markers.append(MarkerSnp(
            chrom=row["chrom"].strip(),
            pos=int(position),
            allele=allele,
            species=row["species"].strip(),
            ref=row.get("ref", "").strip().upper(),
            gene=row.get("gene", "").strip(),
            note=row.get("note", "").strip(),
            source=row.get("source", "").strip(),
        ))
    if not markers:
        raise MjolnirError("{0} contains no marker rows".format(path))
    return markers


@dataclass
class MarkerResult:
    """What the marker panel said, including when it said nothing."""

    species: str = ""
    observed: int = 0
    callable_sites: int = 0
    total: int = 0
    per_species: Dict[str, Dict[str, int]] = field(default_factory=dict)
    conflict: bool = False
    reason: str = ""
    sources: List[str] = field(default_factory=list)


def genotype_markers(markers: Sequence[MarkerSnp], counts: PileupCounts,
                     min_site_depth: int = 1) -> MarkerResult:
    """Resolve a complex member from marker alleles observed in the pileup.

    A species is named only when every one of its callable markers carries the
    marker allele and no other species' markers do the same. Two species both
    fully supported is a conflict, and a conflict is reported as one — inside
    MAC that pattern is what a mixture of *M. chimaera* and *M. intracellulare*
    would look like, and picking the larger count would erase it.
    """
    result = MarkerResult(total=len(markers))
    if not markers:
        result.reason = "no marker table was loaded"
        return result
    if not counts:
        result.reason = "no pileup was supplied at the marker positions"
        return result

    per_species: Dict[str, Dict[str, int]] = {}
    sources: List[str] = []
    for marker in markers:
        stats = per_species.setdefault(
            marker.species, {"total": 0, "callable": 0, "observed": 0})
        stats["total"] += 1
        if marker.source and marker.source not in sources:
            sources.append(marker.source)
        site = counts.get(marker.key)
        if not site:
            continue
        observed, depth, _tie = major_allele(site)
        if depth < min_site_depth or not observed:
            continue
        stats["callable"] += 1
        result.callable_sites += 1
        if observed == marker.allele:
            stats["observed"] += 1

    result.per_species = per_species
    result.sources = sources
    complete = [name for name, stats in per_species.items()
                if stats["callable"] > 0 and stats["observed"] == stats["callable"]]
    if not complete:
        result.reason = (
            "no species' markers were fully supported at the callable marker "
            "positions ({0} of {1} markers callable)".format(
                result.callable_sites, result.total)
        )
        return result
    if len(complete) > 1:
        result.conflict = True
        result.reason = (
            "markers for more than one species were fully supported ({0}), which "
            "is what a mixture looks like; no species is named".format(
                ", ".join(sorted(complete)))
        )
        return result
    result.species = complete[0]
    result.observed = per_species[result.species]["observed"]
    result.callable_sites = per_species[result.species]["callable"]
    result.total = per_species[result.species]["total"]
    return result


# ---------------------------------------------------------------------------
# The call
# ---------------------------------------------------------------------------

def _species_call(name: str, complex_name: str, method: str, *,
                  resolved: bool, ani: Optional[float],
                  aligned_fraction: Optional[float], reference: str,
                  confidence: str, candidates: Sequence[Dict[str, Any]],
                  caveats: Sequence[str], subspecies: str = "") -> SpeciesCall:
    """Build the record, refusing the two things that must never be emitted.

    First, a method outside :data:`ANI_METHODS`. Second, an MTBC member name
    marked as a resolved species: even if a reference manifest carries
    *M. bovis* as a genome, ANI cannot make that claim, and the call is demoted
    to the complex with the reason attached rather than printed.
    """
    if method not in ANI_METHODS:
        raise MjolnirError(
            "refusing to build a species call with method {0!r}. {1}".format(
                method, SPECIES_METHOD_REFUSAL)
        )
    caveat_list = list(caveats)
    if resolved and is_mtbc_member(name) and normalise_species_name(name) not in (
            "mycobacterium tuberculosis complex",):
        LOG.debug("demoting resolved species %r to %s", name, COMPLEX_MTBC)
        resolved = False
        complex_name = COMPLEX_MTBC
        name = COMPLEX_LONG_NAME[COMPLEX_MTBC]
        if MTBC_UNRESOLVED_TEXT not in caveat_list:
            caveat_list.append(MTBC_UNRESOLVED_TEXT)
    return SpeciesCall(
        name=name,
        complex=complex_name,
        method=method,
        ani=None if ani is None else round(ani, 4),
        confidence=confidence,
        resolved_to_species=resolved,
        reference=reference,
        aligned_fraction=aligned_fraction,
        candidates=list(candidates),
        subspecies=subspecies,
        caveats=caveat_list,
    )


def _usable(match: AniMatch) -> bool:
    """Whether a hit clears the species boundary at all."""
    return match.ani is not None and match.ani >= ANI_SPECIES_FLOOR


def _aligned_fraction_ok(match: AniMatch) -> Optional[bool]:
    """True/False when the aligned fraction was measured, None when it was not."""
    if match.aligned_fraction is None:
        return None
    return match.aligned_fraction >= ANI_MIN_ALIGNED_FRACTION


def identify_species(query: PathLike, db_dir: Optional[PathLike] = None, *,
                     is_reads: bool = False, threads: int = 1,
                     references: Optional[Sequence[ReferenceGenome]] = None,
                     matches: Optional[Sequence[AniMatch]] = None,
                     markers: Optional[Sequence[MarkerSnp]] = None,
                     marker_counts: Optional[PileupCounts] = None,
                     tool: Optional[str] = None,
                     runner: Optional[Runner] = None) -> SpeciesCall:
    """Identify the isolate from ANI, and refuse to over-resolve.

    The three outcomes, in the order the evidence allows them:

    * a named species, when the best hit clears the ANI floor over an adequate
      aligned fraction and does not sit inside a complex ANI cannot resolve;
    * a complex, for MTBC always and for MAC unless the marker panel resolves
      it, with the reason stated on the call;
    * unresolved, when nothing clears the floor — which is a real result for a
      genus this reference set does not cover, and is not softened into a guess.

    *matches* is accepted so the pipeline can reuse an ANI run it already did
    (the contamination screen wants the same numbers), and so this function is
    testable without either tool installed.
    """
    if matches is None:
        reference_set = list(references) if references is not None \
            else load_reference_set(db_dir)
        matches = ani_matches(query, reference_set, is_reads=is_reads,
                              threads=threads, db_dir=db_dir, tool=tool,
                              runner=runner)
    ordered = sorted(matches, key=lambda m: (-(m.ani or 0.0), m.name))
    candidates = [m.as_candidate() for m in ordered[:10]]
    method = ordered[0].method if ordered else (METHOD_MASH if is_reads else METHOD_SKANI)

    base_caveats: List[str] = []
    if method == METHOD_MASH:
        base_caveats.append(MASH_ESTIMATE_CAVEAT)
    if is_reads:
        base_caveats.append(READS_METHOD_NOTE)

    usable = [m for m in ordered if _usable(m)]
    if not usable:
        best_ani = ordered[0].ani if ordered else None
        reason = (
            "no reference in the mycobacterial set reached the {0}% ANI species "
            "boundary{1}; the isolate is not identified rather than assigned to "
            "the nearest match".format(
                ANI_SPECIES_FLOOR,
                "" if best_ani is None else " (best was {0:.2f}%)".format(best_ani))
        )
        return _species_call(
            "unresolved", "", method, resolved=False, ani=best_ani,
            aligned_fraction=ordered[0].aligned_fraction if ordered else None,
            reference=ordered[0].reference if ordered else "",
            confidence=CONFIDENCE_NONE, candidates=candidates,
            caveats=base_caveats + [reason])

    best = usable[0]
    complex_name = best.complex or complex_for(best.name)
    aligned_ok = _aligned_fraction_ok(best)
    caveats = list(base_caveats)

    if aligned_ok is False:
        caveats.append(
            "the best ANI match aligned only {0:.1f}% of the query, below the "
            "{1:.0%} Mjolnir requires for a species claim; a high ANI over a "
            "small aligned fraction is not an identification".format(
                (best.aligned_fraction or 0.0) * 100.0, ANI_MIN_ALIGNED_FRACTION)
        )
        return _species_call(
            COMPLEX_LONG_NAME.get(complex_name, "unresolved"), complex_name, method,
            resolved=False, ani=best.ani, aligned_fraction=best.aligned_fraction,
            reference=best.reference, confidence=CONFIDENCE_LOW,
            candidates=candidates, caveats=caveats)
    if aligned_ok is None:
        caveats.append(
            "this method reports no aligned fraction, so the ANI value is not "
            "qualified by how much of the query it covers")

    if complex_name == COMPLEX_MTBC:
        return _call_mtbc(best, method, candidates, caveats)
    if complex_name == COMPLEX_MAC:
        return _call_mac(best, usable, method, candidates, caveats,
                         markers=markers, marker_counts=marker_counts, db_dir=db_dir)

    # Outside the two complexes ANI cannot resolve, ANI is the identification.
    competitors = sorted({m.name for m in usable if m.name != best.name})
    if competitors:
        caveats.append(
            "other references also cleared the species boundary ({0}); the call "
            "is the highest ANI and the margin is in the annex".format(
                ", ".join(competitors[:4]))
        )
        confidence = CONFIDENCE_MODERATE
    else:
        confidence = CONFIDENCE_HIGH if aligned_ok else CONFIDENCE_MODERATE
    return _species_call(
        best.name, complex_name, method, resolved=True, ani=best.ani,
        aligned_fraction=best.aligned_fraction, reference=best.reference,
        confidence=confidence, candidates=candidates, caveats=caveats,
        subspecies=best.subspecies)


def _call_mtbc(best: AniMatch, method: str, candidates: Sequence[Dict[str, Any]],
               caveats: Sequence[str]) -> SpeciesCall:
    """MTBC: the complex, always, with the member left to the SNP barcode.

    There is no configuration that makes this resolve to a species. The members
    lie at 99.21-99.92% ANI of one another and are later heterotypic synonyms of
    *M. tuberculosis* in NCBI taxonomy, so the number ANI produces inside the
    complex carries no information about which member this is. The member comes
    from ``typing/lineage.py``, from lineage-defining SNPs.
    """
    notes = list(caveats)
    notes.append(MTBC_UNRESOLVED_TEXT)
    low, high = MTBC_INTRA_ANI_RANGE
    if best.ani is not None and best.ani < low:
        notes.append(
            "the best MTBC reference matched at {0:.2f}% ANI, below the {1}% "
            "floor of the range MTBC members occupy with one another, so "
            "membership of the complex itself is not firm".format(best.ani, low)
        )
        confidence = CONFIDENCE_LOW
    else:
        confidence = CONFIDENCE_HIGH
    notes.append(
        "the MTBC member (lineage, animal lineage or BCG) is called from "
        "lineage-defining SNPs, not from ANI and never from a read classifier")
    return _species_call(
        COMPLEX_LONG_NAME[COMPLEX_MTBC], COMPLEX_MTBC, method, resolved=False,
        ani=best.ani, aligned_fraction=best.aligned_fraction,
        reference=best.reference, confidence=confidence, candidates=candidates,
        caveats=notes)


def _call_mac(best: AniMatch, usable: Sequence[AniMatch], method: str,
              candidates: Sequence[Dict[str, Any]], caveats: Sequence[str],
              markers: Optional[Sequence[MarkerSnp]],
              marker_counts: Optional[PileupCounts],
              db_dir: Optional[PathLike]) -> SpeciesCall:
    """MAC: ANI must clear the policy floor *and* the markers must agree.

    *M. chimaera* was separated from *M. intracellulare* on ITS sequevar, not on
    a genome-wide distance, and the two sit close enough that a 95% boundary
    says nothing about which one an isolate is. So the species is named only
    when the ANI to that species' reference clears the higher MAC floor and a
    marker panel independently agrees. Anything less is MAC at complex level —
    which is a real answer, and the one the *M. chimaera* outbreak work needs
    said out loud rather than guessed.
    """
    notes = list(caveats)
    marker_set = list(markers) if markers is not None else load_mac_markers(db_dir)
    marker_result = genotype_markers(marker_set, marker_counts or {})

    ani_ok = best.ani is not None and best.ani >= MAC_SPECIES_ANI_FLOOR
    rivals = sorted({m.name for m in usable
                     if m.name != best.name and m.ani is not None
                     and m.ani >= MAC_SPECIES_ANI_FLOOR})

    if not ani_ok:
        notes.append(
            "the best MAC reference matched at {0} ANI, below the {1}% Mjolnir "
            "requires before naming a MAC species; within MAC the 95% species "
            "boundary does not separate M. chimaera from M. intracellulare".format(
                "no measured" if best.ani is None else "{0:.2f}%".format(best.ani),
                MAC_SPECIES_ANI_FLOOR)
        )
    if marker_result.species:
        notes.append(
            "MAC marker panel: {0} supported by {1} of {2} callable markers{3}".format(
                marker_result.species, marker_result.observed,
                marker_result.callable_sites,
                " ({0})".format("; ".join(marker_result.sources))
                if marker_result.sources else "")
        )
    else:
        notes.append("MAC marker panel did not resolve a species: "
                     + (marker_result.reason or "no result"))
    if marker_result.species and ani_ok and \
            normalise_species_name(marker_result.species) != normalise_species_name(best.name):
        notes.append(
            "the marker panel names {0} while the highest ANI is to {1}; the two "
            "lines of evidence disagree and no species is named".format(
                marker_result.species, best.name)
        )
        return _species_call(
            COMPLEX_LONG_NAME[COMPLEX_MAC], COMPLEX_MAC, method, resolved=False,
            ani=best.ani, aligned_fraction=best.aligned_fraction,
            reference=best.reference, confidence=CONFIDENCE_LOW,
            candidates=candidates, caveats=notes)

    if ani_ok and marker_result.species:
        if rivals:
            notes.append(
                "other MAC references also cleared {0}% ANI ({1}); the marker "
                "panel is what separates them".format(
                    MAC_SPECIES_ANI_FLOOR, ", ".join(rivals[:4]))
            )
        return _species_call(
            best.name, COMPLEX_MAC, method, resolved=True, ani=best.ani,
            aligned_fraction=best.aligned_fraction, reference=best.reference,
            confidence=CONFIDENCE_HIGH if not rivals else CONFIDENCE_MODERATE,
            candidates=candidates, caveats=notes, subspecies=best.subspecies)

    notes.append(
        "MAC is reported at complex level: ANI alone does not separate "
        "M. chimaera, M. intracellulare and M. avium, and a species is named "
        "only when the ANI floor and an independent marker panel agree")
    return _species_call(
        COMPLEX_LONG_NAME[COMPLEX_MAC], COMPLEX_MAC, method, resolved=False,
        ani=best.ani, aligned_fraction=best.aligned_fraction,
        reference=best.reference, confidence=CONFIDENCE_MODERATE,
        candidates=candidates, caveats=notes)


# ---------------------------------------------------------------------------
# Checks for the report
# ---------------------------------------------------------------------------

def species_checks(call: SpeciesCall) -> List[Check]:
    """Rule-derived verdicts the report and the agent read.

    Written here rather than in the report so that the number, the threshold and
    the source travel together — the report never has to remember that 95% came
    from Richter & Rossello-Mora.
    """
    checks: List[Check] = [
        Check.numeric(
            "species_ani",
            call.ani,
            warn_minimum=float(ANI_SPECIES_FLOOR),
            source=source_for("ani_species_floor"),
            unit="% ANI",
            category="typing",
            reading=("average nucleotide identity to the closest reference in the "
                     "curated mycobacterial set"),
            not_measured_why=("no ANI was measured, so no species-level claim is "
                              "made from sequence identity"),
        ),
        Check.numeric(
            "species_aligned_fraction",
            call.aligned_fraction,
            warn_minimum=float(ANI_MIN_ALIGNED_FRACTION),
            source=source_for("ani_min_aligned_fraction"),
            unit="fraction",
            category="typing",
            reading="fraction of the query covered by the ANI alignment",
            not_measured_why=("the ANI method used reports no aligned fraction, so "
                              "the identity value is not qualified by coverage"),
        ),
    ]
    if call.complex == COMPLEX_MTBC:
        checks.append(Check(
            name="species_resolved_below_complex",
            value=False,
            threshold=False,
            source=source_for("mtbc_intra_ani_range"),
            status=STATUS_PASS,
            reading=MTBC_UNRESOLVED_TEXT,
            comparison="==",
            category="typing",
        ))
    else:
        checks.append(Check.boolean(
            "species_resolved_below_complex",
            call.resolved_to_species,
            expected=True,
            source=source_for("mac_species_ani_floor" if call.complex == COMPLEX_MAC
                              else "ani_species_floor"),
            category="typing",
            reading=("a species was named" if call.resolved_to_species
                     else "the evidence reaches the complex but not the species"),
            fail_status=STATUS_WARN,
        ))
    checks.append(Check(
        name="species_method_is_ani",
        value=call.method in ANI_METHODS,
        threshold=True,
        source=source_for("species_method_refusal"),
        status=STATUS_PASS if call.method in ANI_METHODS else STATUS_FAIL,
        reading=SPECIES_METHOD_REFUSAL,
        comparison="==",
        category="typing",
    ))
    return checks


def describe_reference_set(references: Iterable[ReferenceGenome]) -> str:
    """One line naming what the identification was made against."""
    items = list(references)
    present = sum(1 for r in items if r.exists)
    complexes = sorted({r.complex for r in items if r.complex})
    return "{0} mycobacterial reference genomes ({1} present on disk){2}".format(
        len(items), present,
        "; complexes represented: " + ", ".join(complexes) if complexes else "")
