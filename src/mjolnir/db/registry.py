"""What every database is, where it comes from, and whether we may ship it.

This module is a declaration, not a downloader. It holds the one description of
each database that the fetcher, the doctor and the report all read, so that the
licence printed in an annex and the licence the fetcher enforced cannot drift
apart.

Three things here are load-bearing rather than documentary.

**Redistribution is a field the code acts on.** Mjolnir is MIT. Where a source's
licence does not permit an MIT project to ship its files onward, the database is
fetched on the installing machine instead of vendored into the distribution, and
:func:`check_redistribution` refuses to build a bundle containing it. The WHO
catalogue *data* — the xlsx and the coordinate VCF — is ODC-By v1.0 and may be
redistributed with attribution; the WHO 2nd-edition report PDF is CC BY-NC-SA
3.0 IGO and may not; MTBseq is GPL-3.0 and tbdb is LGPL-3.0, both copyleft. Only
the first of those could be mirrored, and even it carries an attribution
obligation that travels with the file, which is why ``LICENSE.md`` is fetched
beside the catalogue rather than left behind.

**The catalogue edition is data, not a constant.** WHO v2 is the latest edition
as of the snapshot date below, but a 3rd edition was called for on 2024-08-26
and is unreleased. It is registered as one member of the ``who-catalogue``
family, and a v3 becomes available by adding a second :class:`Database` and
setting ``superseded_by`` on this one. Nothing downstream asks "is this v2"; it
asks :func:`latest_in_family`.

**Integrity is anchored to what could actually be checked.** No file here was
downloaded on the machine that wrote this registry, so no sha256 could be
computed, and inventing one would be worse than leaving it empty. What GitHub's
API does give without transferring the file is the git blob object id, and those
are recorded verbatim. ``fetch.py`` recomputes the blob id locally from the
bytes it received, so the pin is real; the sha256 it computes at fetch time is
what the report prints and what a second installation is compared against.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ..config import (
    CATALOGUE_MTBSEQ,
    CATALOGUE_TBDB,
    CATALOGUE_WHO,
    DB_ENV_VAR,
    H37RV_ACCESSION,
    H37RV_LENGTH,
    KRAKEN2_DB_ENV_VAR,
    KRAKEN2_UNINFORMATIVE_TEXT,
)
from ..utils import MjolnirError, PathLike

#: The day the sizes, blob ids and licence facts below were read from upstream.
#: Everything in this file is a claim about upstream on that date; ``fetch.py``
#: says so out loud when what it receives differs.
SNAPSHOT_DATE = "2026-08-10"

#: Mjolnir's own licence, which is what makes redistribution a question at all.
PROJECT_LICENCE = "MIT"

#: SOURCE: design §12. Stated once, printed by ``mjolnir db list``, because a
#: user who does not know why a database is downloaded rather than bundled will
#: assume the download is an oversight and mirror it themselves.
REDISTRIBUTION_POLICY = (
    "Mjolnir is {0}-licensed. A database whose licence does not permit an {0} "
    "project to redistribute it is fetched on the installing machine at install "
    "time rather than vendored into the distribution, and that decision is "
    "recorded per database below.".format(PROJECT_LICENCE)
)


# ---------------------------------------------------------------------------
# Licences
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Licence:
    """One licence, and the single question the code asks of it.

    ``redistributable`` is not a summary of the licence — it is the narrower
    question of whether *this MIT project* may ship the files onward. A licence
    can be entirely permissive for use and analysis and still answer no here;
    GPL-3.0 is the obvious case.
    """

    spdx: str
    name: str
    url: str = ""
    redistributable: bool = False
    #: The attribution that must travel with the data when it is redistributed
    #: or its results are published. Empty when the licence asks for none.
    attribution: str = ""
    #: True when the licence file must be re-read on the installing machine
    #: because upstream may have changed it, or because its application to data
    #: files (as opposed to code) is not settled.
    verify_at_fetch: bool = False
    note: str = ""

    def describe(self) -> str:
        verdict = ("redistributable with attribution" if self.redistributable
                   else "not redistributable by Mjolnir; fetched at install time")
        return "{0} ({1}) - {2}".format(self.spdx, self.name, verdict)


#: Verified 2026-08-10 from `Final Result Files/LICENSE.md` in
#: GTB-tbsequencing/mutation-catalogue-2023, whose first line reads: "All data
#: published here, in excel or VCF format, are licensed under the Open Data
#: Commons Attribution License (ODC-By) v1.0." The repository's root LICENSE is
#: MIT and covers the pipeline code, not the data — two different licences in
#: one repository, which is exactly the trap this field exists to avoid.
ODC_BY_1_0 = Licence(
    spdx="ODC-By-1.0",
    name="Open Data Commons Attribution License v1.0",
    url="https://opendatacommons.org/licenses/by/1-0/",
    redistributable=True,
    attribution=(
        "World Health Organization. Catalogue of mutations in Mycobacterium "
        "tuberculosis complex and their association with drug resistance, 2nd "
        "edition. Geneva: WHO; 2023 (WHO-UCN-TB-2023.7). Data under ODC-By v1.0."
    ),
    note=("covers the xlsx and VCF data files only; the PDFs in the same "
          "directory are WHO publications under a different licence"),
)

#: The WHO 2nd-edition report itself. The NonCommercial and ShareAlike terms are
#: both incompatible with redistribution inside an MIT distribution, and the
#: document is a human-readable report rather than a machine-readable input, so
#: Mjolnir links to it and never downloads it.
CC_BY_NC_SA_3_0_IGO = Licence(
    spdx="CC-BY-NC-SA-3.0-IGO",
    name="Creative Commons Attribution-NonCommercial-ShareAlike 3.0 IGO",
    url="https://creativecommons.org/licenses/by-nc-sa/3.0/igo/",
    redistributable=False,
    attribution=(
        "World Health Organization. Catalogue of mutations in Mycobacterium "
        "tuberculosis complex and their association with drug resistance, 2nd "
        "edition. Geneva: WHO; 2023. Licence: CC BY-NC-SA 3.0 IGO."
    ),
    note="NonCommercial and ShareAlike terms; the PDF is documentation, not an input",
)

#: Verified 2026-08-10 from LICENSE.md in ngs-fzb/MTBseq_source. Note its exact
#: wording: "The code of this software can be redistributed and/or modified
#: under the terms of the GNU General Public License ... version 3". The licence
#: scopes itself to the code and says nothing about the resistance lists or the
#: reference FASTAs, so the safe reading is that the repository licence governs
#: them and that copyleft applies. Either way an MIT project does not vendor them.
GPL_3_0 = Licence(
    spdx="GPL-3.0-or-later",
    name="GNU General Public License v3.0 or later",
    url="https://www.gnu.org/licenses/gpl-3.0.html",
    redistributable=False,
    attribution=(
        "Kohl TA, Utpatel C, Schleusener V, et al. MTBseq: a comprehensive "
        "pipeline for whole genome sequence analysis of Mycobacterium "
        "tuberculosis complex isolates. PeerJ 2018;6:e5895."
    ),
    verify_at_fetch=True,
    note=("the licence text scopes to 'the code'; the data files under var/ are "
          "not separately licensed, so they are fetched rather than vendored"),
)

#: Verified 2026-08-10: jodyphelan/tbdb reports LGPL-3.0 through the GitHub API
#: and its LICENCE file is the LGPL v3 text. Whether "library" semantics mean
#: anything for a BED file or a CSV is unresolved — design §14 lists it as an
#: open question — so the honest answer is to fetch and to re-read the licence
#: on the installing machine rather than to decide the question here.
LGPL_3_0 = Licence(
    spdx="LGPL-3.0-or-later",
    name="GNU Lesser General Public License v3.0 or later",
    url="https://www.gnu.org/licenses/lgpl-3.0.html",
    redistributable=False,
    attribution=(
        "Phelan JE, O'Sullivan DM, Machado D, et al. Integrating informatics "
        "tools and portable sequencing technology for rapid detection of "
        "resistance to anti-tuberculous drugs. Genome Med 2019;11:41 (tbdb / "
        "TB-Profiler library)."
    ),
    verify_at_fetch=True,
    note=("copyleft; application to data files is unsettled (design §14), so "
          "tbdb is fetched and its LICENCE file is fetched beside it"),
)

#: NCBI states that it places no restrictions on the use or distribution of the
#: data in GenBank, while noting that individual records may carry third-party
#: claims. NC_000962.3 is a public reference genome with no such claim.
NCBI_PUBLIC = Licence(
    spdx="NCBI-public",
    name="NCBI GenBank - no restrictions placed by NCBI on use or distribution",
    url="https://www.ncbi.nlm.nih.gov/home/about/policies/",
    redistributable=True,
    attribution="GenBank {0}, Mycobacterium tuberculosis H37Rv".format(H37RV_ACCESSION),
    note="individual GenBank records can carry third-party restrictions; this one does not",
)

#: For entries Mjolnir describes but does not obtain, where the operator's own
#: build or their institution's copy is the source.
OPERATOR_SUPPLIED = Licence(
    spdx="operator-supplied",
    name="built or obtained by the operator; licence is theirs to determine",
    redistributable=False,
    verify_at_fetch=True,
    note="Mjolnir never downloads this; it only records what it needs and why",
)


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------

#: Verification strategies. A file gets the strongest one that could honestly be
#: applied to it at the time this registry was written.
VERIFY_CHECKSUM = "checksum"       # sha256 if pinned, else the git blob id
VERIFY_FASTA_LENGTH = "fasta-length"  # the sequence is the invariant, not the bytes
VERIFY_NON_EMPTY = "non-empty"     # nothing stable upstream to compare against
VERIFY_STRATEGIES = (VERIFY_CHECKSUM, VERIFY_FASTA_LENGTH, VERIFY_NON_EMPTY)


@dataclass(frozen=True)
class DatabaseFile:
    """One file inside a database, with whatever integrity anchor exists for it.

    ``sha256`` is empty throughout this registry and that is deliberate: no file
    was downloaded on the machine that wrote it, so every digest here would have
    been copied from somewhere unverifiable. ``git_blob_sha1`` was read from the
    GitHub API, which serves the object id without serving the object, and
    ``fetch.py`` recomputes it from the received bytes — a pin that was actually
    obtained beats a digest that was asserted.

    Once an installation has fetched a file, its measured sha256 goes into the
    on-disk manifest and into every report, which is what makes two
    installations comparable (design §5.5).
    """

    name: str
    url: str
    #: Size upstream reported on :data:`SNAPSHOT_DATE`. Used to budget disk and
    #: to notice a truncated or redirected download, never as proof of identity.
    size_bytes: int = 0
    sha256: str = ""
    git_blob_sha1: str = ""
    required: bool = True
    verify: str = VERIFY_CHECKSUM
    #: For VERIFY_FASTA_LENGTH: the sequence length the file must contain.
    expect_sequence_length: Optional[int] = None
    #: For VERIFY_FASTA_LENGTH: a string that must appear in the header.
    expect_header: str = ""
    #: True for tar or zip archives that are expanded in place after
    #: verification. Nothing in the registry needs it today — every current
    #: source publishes plain files — but the ANI reference set of design §14
    #: will arrive as an archive, and an unpack step bolted on later is where a
    #: path-traversal check gets forgotten.
    unpack: bool = False
    note: str = ""

    def __post_init__(self) -> None:
        if self.verify not in VERIFY_STRATEGIES:
            raise MjolnirError(
                "database file {0!r} declares unknown verification strategy "
                "{1!r}; expected one of {2}".format(
                    self.name, self.verify, ", ".join(VERIFY_STRATEGIES)))
        if self.verify == VERIFY_FASTA_LENGTH and not self.expect_sequence_length:
            raise MjolnirError(
                "database file {0!r} asks for {1} verification without an "
                "expected sequence length".format(self.name, VERIFY_FASTA_LENGTH))

    @property
    def has_pin(self) -> bool:
        """Whether anything at all can be compared after the download."""
        return bool(self.sha256 or self.git_blob_sha1
                    or self.verify == VERIFY_FASTA_LENGTH)


# ---------------------------------------------------------------------------
# Databases
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Database:
    """One database as Mjolnir knows it: identity, provenance, licence, contents."""

    name: str
    title: str
    #: Prose for ``mjolnir db list``: what this is and what it is used for.
    what: str
    #: The versioned lineage this entry belongs to. Members of a family are
    #: ordered by ``edition`` and resolved through :func:`latest_in_family`, so
    #: a new catalogue edition is a registration, not an edit to call sites.
    family: str
    version: str
    #: ISO date of the version itself, not of this snapshot.
    version_date: str
    provider: str
    homepage: str
    citation: str
    licence: Licence
    files: Tuple[DatabaseFile, ...] = ()
    #: A fetched JSON file that states the database's own version, and the key
    #: holding it. Where a database says what version it is, that statement is a
    #: fact about the files on disk; ``version`` above is only a fact about the
    #: day this registry was written, so the fetcher prefers the file.
    version_file: str = ""
    version_key: str = ""
    #: The fetched file carrying upstream's own licence text, where one is
    #: fetched. ``Licence.verify_at_fetch`` without one of these is a
    #: registration error: the promise to check the licence on the installing
    #: machine is empty if the licence never lands there.
    licence_file: str = ""
    #: Directory under the database root. Defaults to ``name``.
    subdir: str = ""
    #: Ordinal within the family; higher is newer. 0 for unversioned families.
    edition: int = 0
    #: Set to the name of the successor entry once one is registered.
    superseded_by: str = ""
    #: What would supersede this, and where to look for it. Printed by
    #: ``mjolnir db list`` so an operator can tell "current" from "stale".
    successor_watch: str = ""
    #: False when a human must obtain it; ``mjolnir db fetch`` will not try.
    fetchable: bool = True
    #: In the default fetch set.
    auto: bool = True
    #: Which parts of Mjolnir stop working without it. Used by the doctor to
    #: turn "database missing" into "this is what you lose".
    required_for: Tuple[str, ...] = field(default_factory=tuple)
    note: str = ""

    def __post_init__(self) -> None:
        if self.fetchable and not self.files:
            raise MjolnirError(
                "database {0!r} is marked fetchable but declares no files".format(self.name))
        seen = set()
        for item in self.files:
            if item.name in seen:
                raise MjolnirError(
                    "database {0!r} declares {1!r} twice".format(self.name, item.name))
            seen.add(item.name)
        if self.version_file and self.version_file not in seen:
            raise MjolnirError(
                "database {0!r} reads its version from {1!r} but does not fetch "
                "it".format(self.name, self.version_file))
        if self.version_file and not self.version_key:
            raise MjolnirError(
                "database {0!r} names a version file but no key to read from "
                "it".format(self.name))
        if self.licence_file and self.licence_file not in seen:
            raise MjolnirError(
                "database {0!r} names {1!r} as its licence file but does not "
                "fetch it".format(self.name, self.licence_file))
        if self.fetchable and self.licence.verify_at_fetch and not self.licence_file:
            raise MjolnirError(
                "database {0!r} carries a licence marked for verification at "
                "fetch time ({1}) but fetches no licence file to verify".format(
                    self.name, self.licence.spdx))

    # -- derived -----------------------------------------------------------

    @property
    def redistributable(self) -> bool:
        return self.licence.redistributable

    @property
    def must_fetch(self) -> bool:
        """True when this may not be vendored and has to be fetched on install.

        This is the field the code acts on: :func:`check_redistribution` refuses
        to bundle anything for which it is true, and ``mjolnir db list`` prints
        the reason next to the licence.
        """
        return not self.licence.redistributable

    @property
    def total_bytes(self) -> int:
        return sum(item.size_bytes for item in self.files)

    @property
    def required_bytes(self) -> int:
        return sum(item.size_bytes for item in self.files if item.required)

    def directory(self, db_root: PathLike) -> Path:
        return Path(db_root).expanduser() / (self.subdir or self.name)

    def file(self, name: str) -> DatabaseFile:
        for item in self.files:
            if item.name == name:
                return item
        raise MjolnirError(
            "database {0!r} has no file named {1!r} (has: {2})".format(
                self.name, name, ", ".join(item.name for item in self.files)))

    def required_files(self) -> List[DatabaseFile]:
        return [item for item in self.files if item.required]

    def path_to(self, db_root: PathLike, name: str) -> Path:
        """Where *name* lives once fetched. Raises if the database has no such file."""
        return self.directory(db_root) / self.file(name).name

    def describe(self) -> str:
        return "{0} {1} - {2} [{3}]".format(
            self.name, self.version, self.title, self.licence.spdx)


DATABASES: Dict[str, Database] = {}


def _register(spec: Database) -> Database:
    if spec.name in DATABASES:
        raise MjolnirError("database {0!r} registered twice".format(spec.name))
    DATABASES[spec.name] = spec
    return spec


# ---------------------------------------------------------------------------
# WHO catalogue (design §5.1)
# ---------------------------------------------------------------------------

#: The commit that last touched `Final Result Files/` in
#: GTB-tbsequencing/mutation-catalogue-2023 (2024-05-08), read from the API on
#: SNAPSHOT_DATE. The repository publishes no tags, so a commit id is the only
#: stable ref available; fetching from `main` instead would silently change the
#: catalogue under a running installation.
WHO_V2_COMMIT = "e26e535109dcf5562d57ebac4c3e15d00c1f3c94"
WHO_V2_RAW = (
    "https://raw.githubusercontent.com/GTB-tbsequencing/mutation-catalogue-2023/"
    "{0}/Final%20Result%20Files/".format(WHO_V2_COMMIT)
)

#: Design §5.2, first trap. The repository also ships a `.txt` master file, and
#: it is not the same data: 40,178 rows against 48,152, 14 drugs against 15,
#: Streptomycin absent entirely, rpoB/rpoC/rpsL/gid affected. It is deliberately
#: not registered as a file here, and this string exists so that the reason
#: survives the next person who notices the omission.
WHO_TXT_REFUSAL = (
    "the repository's WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt is not "
    "equivalent to the xlsx - fewer rows, one drug missing, four genes affected "
    "- so Mjolnir reads the xlsx and does not fetch the txt"
)

#: Design §5.1 and the WHO announcement of 2024-08-26. Kept as prose rather than
#: as a boolean, because "no third edition yet" is a statement with a date on it.
WHO_EDITION_WATCH = (
    "WHO v2 (WHO-UCN-TB-2023.7) is the latest edition as of {0}. A call for data "
    "for a 3rd edition was issued on 2024-08-26 with a 2024-10-15 deadline and no "
    "edition has been released. When one is, register it as a new member of the "
    "'who-catalogue' family and set superseded_by on this entry; nothing in "
    "Mjolnir asks whether the catalogue is v2.".format(SNAPSHOT_DATE)
)

WHO_CATALOGUE_V2 = _register(Database(
    name="who-catalogue-v2",
    title="WHO catalogue of mutations in MTBC, 2nd edition (data files)",
    what=("48,152 graded variant-drug rows over 30,699 unique variants, 65 genes "
          "and 15 drugs, with the Genomic_coordinates sheet that makes WHO's own "
          "coordinate-based matching protocol possible. The anchor catalogue: "
          "where it grades a variant, its grade is the Mjolnir call."),
    family="who-catalogue",
    version="v2 (WHO-UCN-TB-2023.7)",
    version_date="2023-11-14",
    provider="World Health Organization / GTB-tbsequencing",
    homepage="https://github.com/GTB-tbsequencing/mutation-catalogue-2023",
    citation=ODC_BY_1_0.attribution,
    licence=ODC_BY_1_0,
    licence_file="LICENSE.md",
    edition=2,
    successor_watch=WHO_EDITION_WATCH,
    required_for=("resistance calling for MTBC", "WHO grades in the report"),
    note=WHO_TXT_REFUSAL,
    files=(
        DatabaseFile(
            name="WHO-UCN-TB-2023.7-eng.xlsx",
            url=WHO_V2_RAW + "WHO-UCN-TB-2023.7-eng.xlsx",
            size_bytes=30943740,
            git_blob_sha1="c7cc287d5396fd5243404279d828fe4ee4c78823",
            note=("header on row 3; Catalogue_master_file and Genomic_coordinates "
                  "sheets are both required"),
        ),
        DatabaseFile(
            name="Genomic_coordinates_7May2024.vcf.gz",
            url=WHO_V2_RAW + "Genomic_coordinates_7May2024.vcf.gz",
            size_bytes=850540,
            git_blob_sha1="71e602ffb599598c8b475547b6fd4310fe5d8975",
            required=False,
            note=("the same coordinates as the xlsx sheet, in VCF form, with MNV "
                  "components joined by '&' in INFO; useful for cross-checking a "
                  "coordinate match, not required to make one"),
        ),
        DatabaseFile(
            name="LICENSE.md",
            url=WHO_V2_RAW + "LICENSE.md",
            size_bytes=19918,
            git_blob_sha1="4fb15501d4343ed5aadcd0253d8269d7218dd54f",
            note=("ODC-By v1.0 text; fetched because the attribution obligation "
                  "travels with the data and a mirror without it is a breach"),
        ),
    ),
))

#: The report PDF. Registered so that ``mjolnir db list`` accounts for it and
#: states its licence, and never fetched: it is a human document under
#: CC BY-NC-SA 3.0 IGO, and nothing in the pipeline reads it.
WHO_REPORT_PDF = _register(Database(
    name="who-catalogue-v2-report",
    title="WHO catalogue of mutations in MTBC, 2nd edition (the report PDF)",
    what=("The narrative report behind the catalogue: how the grades were "
          "derived, what each group means, the additional grading rules Mjolnir "
          "implements in resistance/rules.py. Documentation for the human, not "
          "an input to the pipeline."),
    family="who-catalogue",
    version="v2 report",
    version_date="2023-11-14",
    provider="World Health Organization",
    homepage="https://www.who.int/publications/i/item/9789240082410",
    citation=CC_BY_NC_SA_3_0_IGO.attribution,
    licence=CC_BY_NC_SA_3_0_IGO,
    edition=2,
    fetchable=False,
    auto=False,
    successor_watch=WHO_EDITION_WATCH,
    required_for=(),
    note=("not downloaded and not redistributable: NonCommercial and ShareAlike "
          "terms. Read it at the homepage above. Note that the PDF prints the "
          "grade strings with an en-dash while the xlsx uses a spaced ASCII "
          "hyphen (design §5.2) - matching on the PDF's form fails."),
))


# ---------------------------------------------------------------------------
# tbdb (design §5.1, §6, §9)
# ---------------------------------------------------------------------------

#: tbdb publishes tags (v1.0, who-v2-strict) but develops on master, and the
#: barcode and mask files that matter here change on master rather than at tags.
#: The ref is therefore master, and the commit below is what master pointed at
#: on SNAPSHOT_DATE — recorded so a fetch can say "upstream has moved on" rather
#: than pretending the two installations hold the same file.
TBDB_REF = "master"
TBDB_SNAPSHOT_COMMIT = "618cf0ff5f22886971bd437929d2c49defa6c7bf"
TBDB_RAW = "https://raw.githubusercontent.com/jodyphelan/tbdb/{0}/".format(TBDB_REF)

TBDB = _register(Database(
    name="tbdb",
    title="tbdb - the TB-Profiler library",
    what=("Three things Mjolnir needs and one it cross-checks: barcode.bed, the "
          "1,111-SNP / 126-taxon lineage barcode including the La1/La2/La3 animal "
          "nomenclature and BCG; mask.bed, the repeat and low-complexity mask "
          "that pairwise distances are computed against; mutations.csv, the third "
          "catalogue in the consensus; and genome.fasta, tbdb's own H37Rv, kept "
          "only to confirm that its coordinates and ours are the same."),
    family="tbdb",
    version="master@{0}".format(TBDB_SNAPSHOT_COMMIT[:12]),
    version_date="2026-03-11",
    provider="jodyphelan / TB-Profiler",
    homepage="https://github.com/jodyphelan/tbdb",
    citation=LGPL_3_0.attribution,
    licence=LGPL_3_0,
    licence_file="LICENCE",
    version_file="variables.json",
    version_key="db-schema-version",
    required_for=("MTBC lineage barcode", "cohort distance masking",
                  "the tbdb column of the consensus"),
    successor_watch=("tbdb is a moving library, not an edition: its mask changed "
                     "twice in three years (Modlin blind spots 2023-03, merged "
                     "Marin + Modlin 2025-08). The version recorded at fetch time "
                     "is db-schema-version from variables.json plus the file "
                     "checksums, and the report prints which mask it used."),
    note=("licence re-read at fetch time (design §14): LGPL-3.0 is copyleft and "
          "its application to BED and CSV data files is unsettled, so these are "
          "fetched, never vendored."),
    files=(
        DatabaseFile(
            name="barcode.bed", url=TBDB_RAW + "barcode.bed",
            size_bytes=78083, git_blob_sha1="624f5526fb58c537f0fe6899c78276063e2e79e3",
            note=("8 columns, 1,111 SNP rows, 126 taxa; also carries the "
                  "lineage -> named-family mapping, so no hand-built lookup is "
                  "needed. La1.2.BCG rests on 5 SNPs - typing/lineage.py must "
                  "report support, not just the label."),
        ),
        DatabaseFile(
            name="mask.bed", url=TBDB_RAW + "mask.bed",
            size_bytes=60942, git_blob_sha1="8a6f3d8c58c43fcb6ac4e8d130a83c8a27fa854c",
            note="named, versioned and swappable; the report prints which mask ran",
        ),
        DatabaseFile(
            name="mutations.csv", url=TBDB_RAW + "mutations.csv",
            size_bytes=4709614, git_blob_sha1="2d837b43a1a457fe92715334c4a2cd6fc505b11a",
            note="TB-Profiler's variant library, with its own confidence field",
        ),
        DatabaseFile(
            name="genome.gff", url=TBDB_RAW + "genome.gff",
            size_bytes=2346945,
            note=("H37Rv gene models. Without them a variant is a coordinate and "
                  "nothing else: WHO is matched on coordinates and still works, "
                  "but MTBseq and tbdb are keyed on <gene>_<hgvs> and match "
                  "nothing at all, so the consensus quietly becomes WHO alone "
                  "and the NTM rrl/rrs/erm(41) rules cannot fire."),
        ),
        DatabaseFile(
            name="genome.fasta", url=TBDB_RAW + "genome.fasta",
            size_bytes=4485135, git_blob_sha1="17f502d3a1a23b1f0ba339a364aa02bacd1fe10c",
            required=False,
            verify=VERIFY_FASTA_LENGTH,
            expect_sequence_length=H37RV_LENGTH,
            note=("tbdb's H37Rv. Verified by sequence length rather than by bytes: "
                  "what matters is that its coordinate system is the same "
                  "{0} as ours, not that the file is formatted identically."
                  .format(H37RV_ACCESSION)),
        ),
        DatabaseFile(
            name="variables.json", url=TBDB_RAW + "variables.json",
            size_bytes=631, git_blob_sha1="6f7dc72b374f19985759909ce700f9d3fb1c2c3e",
            note=("carries db-schema-version and tbdb's own drug spellings; the "
                  "fetcher reads the version out of it rather than guessing one"),
        ),
        DatabaseFile(
            name="LICENCE", url=TBDB_RAW + "LICENCE",
            size_bytes=7360, git_blob_sha1="40e88e0b45705e60679df9b7288ab5437032edcb",
            note="fetched so the licence can be verified on the installing machine",
        ),
    ),
))


# ---------------------------------------------------------------------------
# MTBseq (design §5.1, §1)
# ---------------------------------------------------------------------------

#: MTBseq tags its releases, so the ref is a tag and the fetch is reproducible.
MTBSEQ_TAG = "v1.1.0"
MTBSEQ_RAW = "https://raw.githubusercontent.com/ngs-fzb/MTBseq_source/{0}/".format(MTBSEQ_TAG)

MTBSEQ_RESISTANCE = _register(Database(
    name="mtbseq-resistance",
    title="MTBseq resistance-mediating lists and gene categories",
    what=("The second catalogue in the consensus. A flat list with no confidence "
          "grading, so it can only ever contribute R or no-call, and the report "
          "states that asymmetry rather than letting agreement with it look like "
          "corroboration of a WHO grade."),
    family="mtbseq",
    version=MTBSEQ_TAG,
    version_date="2021-11-19",
    provider="FZ Borstel (ngs-fzb)",
    homepage="https://github.com/ngs-fzb/MTBseq_source",
    citation=GPL_3_0.attribution,
    licence=GPL_3_0,
    licence_file="LICENSE.md",
    required_for=("the MTBseq column of the consensus", "--mtbseq-compat runs"),
    note=("these lists are MTB-only. A M. chimaera run under MTBseq gets "
          "--resilist NONE --categories NONE --basecalib NONE, which is the gap "
          "Mjolnir exists to close - they are not a fallback for NTM."),
    files=(
        DatabaseFile(
            name="MTB_Resistance_Mediating.txt",
            url=MTBSEQ_RAW + "var/res/MTB_Resistance_Mediating.txt",
            size_bytes=275033, git_blob_sha1="55e6ce4c84a950e46ed006e40d2ca96d3f82a405",
        ),
        DatabaseFile(
            name="MTB_Extended_Resistance_Mediating.txt",
            url=MTBSEQ_RAW + "var/res/MTB_Extended_Resistance_Mediating.txt",
            size_bytes=3204, git_blob_sha1="d851a0081f05ba618627a3e8e1da111ab9106291",
        ),
        DatabaseFile(
            name="MTB_Gene_Categories.txt",
            url=MTBSEQ_RAW + "var/cat/MTB_Gene_Categories.txt",
            size_bytes=78887, git_blob_sha1="e6a2b1a41754cf84a5823d33d790e415710f2a4d",
            note="gene category table, used when reproducing MTBseq's own output",
        ),
        DatabaseFile(
            name="MTB_Base_Calibration_List.vcf",
            url=MTBSEQ_RAW + "var/res/MTB_Base_Calibration_List.vcf",
            size_bytes=52730, git_blob_sha1="b4072bf2aa69cb110794eb7d3253f6a6d05ea65b",
            required=False,
            note="known sites for BQSR; only --mtbseq-compat has a use for it",
        ),
        DatabaseFile(
            name="LICENSE.md", url=MTBSEQ_RAW + "LICENSE.md",
            size_bytes=1071, git_blob_sha1="31cc2550b7f2a94bfafab0161d3904bd6a082f7f",
            note=("MTBseq's own licence statement, fetched so the scoping "
                  "question - it licenses 'the code' and is silent about var/ - "
                  "can be re-read on the installing machine"),
        ),
    ),
))

MTBSEQ_REFERENCES = _register(Database(
    name="mtbseq-ntm-references",
    title="MTBseq NTM reference genomes",
    what=("The three non-tuberculous references MTBseq ships alongside H37Rv: "
          "M. abscessus CIP-104536T, M. chimaera DSM44623 and M. fortuitum CT6, "
          "each with its gene table. These are what an NTM sample is mapped "
          "against, and the M. chimaera reference is the one the outbreak data on "
          "hand needs."),
    family="mtbseq",
    version=MTBSEQ_TAG,
    version_date="2021-11-19",
    provider="FZ Borstel (ngs-fzb)",
    homepage="https://github.com/ngs-fzb/MTBseq_source",
    citation=GPL_3_0.attribution,
    licence=GPL_3_0,
    licence_file="LICENSE.md",
    required_for=("NTM mapping and variant calling",),
    note=("MTBseq's H37Rv copy is deliberately not fetched here: the canonical "
          "one is the NCBI record (see the h37rv entry), and holding two H37Rv "
          "FASTAs with different headers is how a coordinate system quietly "
          "forks. The BWA index files that sit beside these FASTAs upstream are "
          "also skipped - Mjolnir builds its own with bwa-mem2 or minimap2."),
    files=(
        DatabaseFile(
            name="M._abscessus_CIP-104536T_2014-02-03.fasta",
            url=MTBSEQ_RAW + "var/ref/M._abscessus_CIP-104536T_2014-02-03.fasta",
            size_bytes=5067198, git_blob_sha1="0da36f7e238591fac22ad7aa57d4e8e6871a4d04",
        ),
        DatabaseFile(
            name="M._abscessus_CIP-104536T_2014-02-03_genes.txt",
            url=MTBSEQ_RAW + "var/ref/M._abscessus_CIP-104536T_2014-02-03_genes.txt",
            size_bytes=486462, git_blob_sha1="161c433050e5d2d04cd4bcd387bf105764c4bfc4",
            note="gene coordinates; erm(41) and rrl are located through this",
        ),
        DatabaseFile(
            name="M._chimaera_DSM44623_2016-01-28.fasta",
            url=MTBSEQ_RAW + "var/ref/M._chimaera_DSM44623_2016-01-28.fasta",
            size_bytes=5865677, git_blob_sha1="24faace958535ab6a0f0d01aa2520926f5e63fb5",
        ),
        DatabaseFile(
            name="M._chimaera_DSM44623_2016-01-28_genes.txt",
            url=MTBSEQ_RAW + "var/ref/M._chimaera_DSM44623_2016-01-28_genes.txt",
            size_bytes=559609, git_blob_sha1="dbc4b575b8a20206bef79139497b64919882b90f",
        ),
        DatabaseFile(
            name="M._fortuitum_CT6_2016-01-08.fasta",
            url=MTBSEQ_RAW + "var/ref/M._fortuitum_CT6_2016-01-08.fasta",
            size_bytes=6254645, git_blob_sha1="8d4297c362666256d684d514d9b811ebca3dce73",
        ),
        DatabaseFile(
            name="M._fortuitum_CT6_2016-01-08_genes.txt",
            url=MTBSEQ_RAW + "var/ref/M._fortuitum_CT6_2016-01-08_genes.txt",
            size_bytes=524452, git_blob_sha1="b6a24bb5145e7ab74f6c4f23b12425263983d8ea",
        ),
        DatabaseFile(
            name="LICENSE.md", url=MTBSEQ_RAW + "LICENSE.md",
            size_bytes=1071, git_blob_sha1="31cc2550b7f2a94bfafab0161d3904bd6a082f7f",
            note="the same licence statement, kept beside the references it governs",
        ),
    ),
))


# ---------------------------------------------------------------------------
# H37Rv (design §6, and the coordinate system for everything above)
# ---------------------------------------------------------------------------

#: NCBI E-utilities. `rettype=fasta` on nuccore returns the record with an NCBI
#: header, which differs between fetches in whitespace and description text
#: while the sequence does not — hence VERIFY_FASTA_LENGTH for this one file.
H37RV_EFETCH = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id={0}&rettype=fasta&retmode=text".format(H37RV_ACCESSION)
)

H37RV = _register(Database(
    name="h37rv",
    title="M. tuberculosis H37Rv reference genome ({0})".format(H37RV_ACCESSION),
    what=("The coordinate system every MTBC catalogue in this registry is "
          "written against: the WHO Genomic_coordinates sheet, tbdb's barcode.bed "
          "and mask.bed, and MTBseq's lists. A coordinate against any other "
          "assembly is wrong by up to several kilobases and mis-grades variants "
          "silently."),
    family="reference",
    version=H37RV_ACCESSION,
    version_date="2013-06-19",
    provider="NCBI GenBank",
    homepage="https://www.ncbi.nlm.nih.gov/nuccore/{0}".format(H37RV_ACCESSION),
    citation=("Cole ST, Brosch R, Parkhill J, et al. Deciphering the biology of "
              "Mycobacterium tuberculosis from the complete genome sequence. "
              "Nature 1998;393:537 (GenBank {0})".format(H37RV_ACCESSION)),
    licence=NCBI_PUBLIC,
    required_for=("MTBC mapping", "every catalogue coordinate lookup"),
    note=("the accession is pinned including its version suffix; NC_000962.2 and "
          ".3 do not share coordinates"),
    files=(
        DatabaseFile(
            name="{0}.fasta".format(H37RV_ACCESSION),
            url=H37RV_EFETCH,
            size_bytes=4_500_000,
            verify=VERIFY_FASTA_LENGTH,
            expect_sequence_length=H37RV_LENGTH,
            expect_header=H37RV_ACCESSION,
            note=("verified by sequence length and accession, not by checksum: "
                  "efetch's header text is not byte-stable and the sequence is "
                  "the thing that has to be right"),
        ),
    ),
))


# ---------------------------------------------------------------------------
# Kraken2 (design §8) - described, refused, never fetched
# ---------------------------------------------------------------------------

KRAKEN2_PANGENOME = _register(Database(
    name="kraken2-mycobacterial-pangenome",
    title="Kraken2 mycobacterial pangenome index (operator-built)",
    what=("The only Kraken2 index whose output Mjolnir will treat as a "
          "contamination screen. " + KRAKEN2_UNINFORMATIVE_TEXT),
    family="kraken2",
    version="operator-built",
    version_date="",
    provider="built locally with kraken2-build",
    homepage="https://github.com/DerrickWood/kraken2/wiki/Manual",
    citation=("Wood DE, Lu J, Langmead B. Improved metagenomic analysis with "
              "Kraken 2. Genome Biol 2019;20:257."),
    licence=OPERATOR_SUPPLIED,
    fetchable=False,
    auto=False,
    required_for=("the optional read-composition screen",),
    note=("Mjolnir does not download a Kraken2 index: a standard or capped index "
          "is tens of gigabytes and, measured on real Illumina data, classifies "
          "only 7.31% of M. tuberculosis reads correctly. Point ${0} or "
          "--kraken2-db at a mycobacterial pangenome index and declare it by "
          "writing {{\"mycobacterial_pangenome\": true}} into mjolnir_index.json "
          "beside its hash.k2d; without that declaration the screen is reported "
          "as uninformative rather than as clean.".format(KRAKEN2_DB_ENV_VAR)),
))


# ---------------------------------------------------------------------------
# Groups and lookup
# ---------------------------------------------------------------------------

#: Names accepted anywhere a database name is. ``default`` is what a bare
#: ``mjolnir db fetch`` installs: everything the MTBC and NTM paths need.
DB_GROUPS: Dict[str, Tuple[str, ...]] = {
    "default": tuple(name for name, spec in DATABASES.items() if spec.auto and spec.fetchable),
    "all": tuple(name for name, spec in DATABASES.items() if spec.fetchable),
    "catalogues": ("who-catalogue-v2", "mtbseq-resistance", "tbdb"),
    "references": ("h37rv", "mtbseq-ntm-references"),
    "mtbc": ("who-catalogue-v2", "mtbseq-resistance", "tbdb", "h37rv"),
    "ntm": ("mtbseq-ntm-references", "h37rv"),
}

#: The catalogue name used in a :class:`~mjolnir.records.CatalogueCall` mapped to
#: the database that supplies it, so the report can print the version and
#: checksum of the catalogue that produced each call (design §5.5).
CATALOGUE_DATABASES: Dict[str, str] = {
    CATALOGUE_WHO: "who-catalogue-v2",
    CATALOGUE_MTBSEQ: "mtbseq-resistance",
    CATALOGUE_TBDB: "tbdb",
}


def spec_for(name: str) -> Database:
    """The registry entry called *name*, or a MjolnirError listing what exists."""
    if name in DATABASES:
        return DATABASES[name]
    raise MjolnirError(
        "unknown database {0!r}. Known databases: {1}. Groups: {2}".format(
            name, ", ".join(sorted(DATABASES)), ", ".join(sorted(DB_GROUPS))))


def resolve_names(names) -> List[str]:
    """Expand group names and comma-separated lists into database names.

    An empty request resolves to the ``default`` group rather than to nothing,
    because ``mjolnir db fetch`` with no arguments meaning "fetch nothing" is a
    silence that reads as success.
    """
    if not names:
        return list(DB_GROUPS["default"])
    out: List[str] = []
    for raw in names:
        for token in str(raw).split(","):
            token = token.strip()
            if not token:
                continue
            if token in DB_GROUPS:
                candidates = list(DB_GROUPS[token])
            else:
                candidates = [spec_for(token).name]
            for candidate in candidates:
                if candidate not in out:
                    out.append(candidate)
    return out


def catalogue_database(catalogue: str) -> Database:
    """The registry entry behind a catalogue name used in resistance calls."""
    name = CATALOGUE_DATABASES.get(catalogue)
    if name is None:
        raise MjolnirError(
            "no database registered for catalogue {0!r}; known catalogues: "
            "{1}".format(catalogue, ", ".join(sorted(CATALOGUE_DATABASES))))
    return DATABASES[name]


# ---------------------------------------------------------------------------
# Families and succession
# ---------------------------------------------------------------------------

def families() -> List[str]:
    return sorted({spec.family for spec in DATABASES.values()})


def members_of(family: str) -> List[Database]:
    """Every entry in *family*, oldest edition first."""
    found = [spec for spec in DATABASES.values() if spec.family == family]
    if not found:
        raise MjolnirError(
            "no database family {0!r}; known families: {1}".format(
                family, ", ".join(families())))
    return sorted(found, key=lambda spec: (spec.edition, spec.name))


def latest_in_family(family: str, fetchable_only: bool = True) -> Database:
    """The newest entry in *family* that nothing has superseded.

    Call sites ask for ``latest_in_family("who-catalogue")`` rather than naming
    v2, so that publishing a 3rd edition is one registration here and no change
    anywhere else. If two editions were ever both current — which would be a
    registration mistake, not a real state — the higher edition wins and the
    lower one is expected to carry ``superseded_by``.
    """
    candidates = [spec for spec in members_of(family)
                  if not spec.superseded_by and (spec.fetchable or not fetchable_only)]
    if not candidates:
        raise MjolnirError(
            "every entry in database family {0!r} is marked superseded; one of "
            "them should be the current edition".format(family))
    return candidates[-1]


def superseded() -> List[Database]:
    """Entries a newer edition has replaced, for the doctor to warn about."""
    return [spec for spec in DATABASES.values() if spec.superseded_by]


# ---------------------------------------------------------------------------
# Redistribution
# ---------------------------------------------------------------------------

def redistributable(names=None) -> List[Database]:
    """The entries an MIT distribution may legally carry."""
    wanted = resolve_names(names) if names else list(DATABASES)
    return [DATABASES[name] for name in wanted if DATABASES[name].redistributable]


def must_fetch(names=None) -> List[Database]:
    """The entries that have to be obtained on the installing machine."""
    wanted = resolve_names(names) if names else list(DATABASES)
    return [DATABASES[name] for name in wanted if DATABASES[name].must_fetch]


def check_redistribution(names) -> None:
    """Refuse to bundle databases whose licence does not permit it.

    Anything that builds a distributable artefact — a wheel with data in it, a
    conda package, a tarball for an offline site — calls this first. It raises
    rather than filtering, because quietly dropping a database from a bundle
    produces an installation that is missing a catalogue and does not know it.
    """
    blocked = [spec for spec in (spec_for(name) for name in resolve_names(names))
               if spec.must_fetch]
    if not blocked:
        return
    lines = ["these databases may not be redistributed with Mjolnir "
             "({0}); fetch them at install time instead:".format(PROJECT_LICENCE)]
    for spec in blocked:
        lines.append("  {0} - {1}: {2}".format(
            spec.name, spec.licence.spdx, spec.licence.note or spec.licence.name))
    lines.append("  mjolnir db fetch {0}".format(" ".join(spec.name for spec in blocked)))
    raise MjolnirError("\n".join(lines))


def attributions(names=None) -> List[str]:
    """Every attribution string a run using these databases owes.

    The report prints these verbatim. ODC-By is an attribution licence: using
    the WHO data and not naming it is a licence breach, not a style lapse.
    """
    seen: List[str] = []
    for name in (resolve_names(names) if names else sorted(DATABASES)):
        text = DATABASES[name].licence.attribution
        if text and text not in seen:
            seen.append(text)
    return seen


# ---------------------------------------------------------------------------
# Where things live
# ---------------------------------------------------------------------------

#: Written into the database root by ``mjolnir db fetch``; read by every run to
#: report what it used. Named for the tool so that pointing ${MJOLNIR_DB} at a
#: shared directory full of other things stays unambiguous.
MANIFEST_NAME = "mjolnir_databases.json"

#: How the doctor and every "not found" message tell an operator what to do.
FETCH_HINT = "mjolnir db fetch {0}"

#: What to say when a database root has never been populated.
NO_DATABASES_TEXT = (
    "no Mjolnir databases are installed. Run 'mjolnir db fetch' to obtain them, "
    "or set ${0} to a directory that already holds them.".format(DB_ENV_VAR)
)


def fetch_hint(name: str) -> str:
    return FETCH_HINT.format(name)


def expected_paths(db_root: PathLike, name: str) -> List[Path]:
    """Every required file of *name*, where it should be under *db_root*."""
    spec = spec_for(name)
    directory = spec.directory(db_root)
    return [directory / item.name for item in spec.required_files()]
