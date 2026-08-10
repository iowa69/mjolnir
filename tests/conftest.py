"""Fixtures shared across the suite, and the synthetic WHO workbook it is built on.

Two things live here rather than in the test modules.

The first is the import path. Mjolnir is not installed in the environment these
tests run in, so ``src`` goes on ``sys.path`` before anything imports the
package. That is done here rather than in each module because a test file that
forgets it fails with an import error, which reads like a broken test rather
than like a missing install.

The second is the workbook builder. Half of the catalogue traps in the design
(§5.2) are properties of the *file*: the header on row 3, the decoy position
column, the ampersand-joined MNV names, the per-(drug, variant) grading. Testing
them needs a file with exactly those properties, and building one in each test
would let the fixtures drift apart until two tests were asserting against two
different notions of what the published workbook looks like. So there is one
builder, its default rows carry every trap at once, and each test varies the one
thing it is about.

Nothing here downloads anything or shells out. The workbook is written by
openpyxl into ``tmp_path`` and is a few kilobytes.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

SRC = Path(__file__).resolve().parent.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import pytest  # noqa: E402  - after the path fix, deliberately

from mjolnir import config  # noqa: E402
from mjolnir.records import (CALL_NO_CALL, CatalogueCall, Variant)  # noqa: E402

# ---------------------------------------------------------------------------
# The synthetic WHO workbook
# ---------------------------------------------------------------------------

#: The master sheet's columns, in the published order, including the duplicated
#: ``CHANGES vs ver1`` pair that a dict-based header index silently collapses.
WHO_MASTER_COLUMNS: List[str] = [
    "drug", "gene", "mutation", "variant", "tier", "effect", "genomic position",
    "FINAL CONFIDENCE GRADING", "INITIAL CONFIDENCE GRADING", "Comment",
    "Additional grading criteria applied", "CHANGES vs ver1", "CHANGES vs ver1",
]

#: The decoy: 38,884 of the published file's 48,152 rows carry this literal
#: string where a naive loader expects a coordinate.
DECOY = config.WHO_COORDINATE_DECOY

#: The en-dash spelling the published PDF prints. Matching on it silently
#: matches nothing in the workbook, which is why it is here as a fixture value
#: and never as an expected value.
GRADE_2_ENDASH = config.WHO_GRADE_2.replace(" - ", " – ")


def who_master_rows(include_streptomycin: bool = True) -> List[List[Any]]:
    """Master-sheet rows carrying every trap the design records as verified.

    ``inhA_c.-154G>A`` appears twice, graded differently for isoniazid and for
    ethionamide, so a loader that deduplicates by variant name loses one of them.
    ``rpoB_p.Ser450Leu`` is graded with the PDF's en-dash. ``eis_c.-14C>T``
    carries a plausible number in the decoy position column and appears nowhere
    in the coordinates sheet, so any coordinate it ends up with came from the
    column that must never be read as one.
    """
    rows: List[List[Any]] = [
        ["Isoniazid", "inhA", "c.-154G>A", "inhA_c.-154G>A", 1,
         "upstream_gene_variant", DECOY, config.WHO_GRADE_1, "", "", ""],
        ["Ethionamide", "inhA", "c.-154G>A", "inhA_c.-154G>A", 1,
         "upstream_gene_variant", DECOY, config.WHO_GRADE_2, "", "", ""],
        ["Rifampicin", "rpoB", "p.Ser450Leu", "rpoB_p.Ser450Leu", 1,
         "missense_variant", DECOY, GRADE_2_ENDASH, "", "High-level resistance", ""],
        ["Rifampicin", "rpoB", "p.Leu452Pro", "rpoB_p.Leu452Pro", 1,
         "missense_variant", DECOY, config.WHO_GRADE_3, "", "", ""],
        ["Amikacin", "eis", "c.-14C>T", "eis_c.-14C>T", 1,
         "upstream_gene_variant", "2715342", config.WHO_GRADE_1, "", "", ""],
        ["Bedaquiline", "Rv0678", "LoF", "Rv0678_LoF", 1,
         "LoF", DECOY, config.WHO_GRADE_2, "", "", ""],
    ]
    if include_streptomycin:
        rows.append(
            ["Streptomycin", "rpsL", "p.Lys43Arg", "rpsL_p.Lys43Arg", 1,
             "missense_variant", DECOY, config.WHO_GRADE_1, "", "", ""])
    return rows


def who_coordinate_rows() -> List[List[str]]:
    """Coordinate rows, including one MNV whose graded variants are ``&``-joined.

    One genomic change at 761155 decomposes into two graded variants. The
    ampersand is how the published VCF and coordinate table write that, and
    indexing the joined string produces a name no master row carries.
    """
    return [
        ["inhA_c.-154G>A", "NC_000962.3", "1673432", "G", "A"],
        ["rpsL_p.Lys43Arg", "NC_000962.3", "781687", "A", "G"],
        ["rpoB_p.Ser450Leu&rpoB_p.Leu452Pro", "NC_000962.3", "761155", "C", "T"],
    ]


def write_who_workbook(path: Path, *,
                       rows: Optional[Sequence[Sequence[Any]]] = None,
                       coordinate_rows: Optional[Sequence[Sequence[str]]] = None,
                       header_row: int = 3,
                       master_sheet: str = config.WHO_XLSX_MASTER_SHEET,
                       coordinates_sheet: str = config.WHO_XLSX_COORDINATES_SHEET,
                       ) -> Path:
    """Write a workbook shaped like the published one.

    ``header_row`` exists so a test can build the file a naive reader assumes —
    header on row 1 — and watch the loader refuse it.
    """
    from openpyxl import Workbook

    workbook = Workbook()
    master = workbook.active
    master.title = master_sheet
    if header_row == 3:
        master.append(["WHO catalogue of mutations in Mycobacterium tuberculosis "
                       "complex and their association with drug resistance"])
        master.append(["Second edition, 2023 - merged banner, two rows deep"])
    elif header_row != 1:  # pragma: no cover - defensive, no test needs it
        raise ValueError("the fixture builds the header on row 1 or row 3 only")
    master.append(list(WHO_MASTER_COLUMNS))
    for row in (who_master_rows() if rows is None else rows):
        master.append(list(row))

    coordinates = workbook.create_sheet(coordinates_sheet)
    coordinates.append(["variant", "chromosome", "position",
                        "reference_nucleotide", "alternative_nucleotide"])
    for row in (who_coordinate_rows() if coordinate_rows is None else coordinate_rows):
        coordinates.append(list(row))

    path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(str(path))
    return path


@pytest.fixture
def make_who_workbook(tmp_path):
    """Factory for a workbook: ``make_who_workbook(name=..., rows=..., ...)``."""

    def _make(name: str = "WHO-UCN-TB-2023.7-eng.xlsx", **kwargs: Any) -> Path:
        return write_who_workbook(tmp_path / name, **kwargs)

    return _make


@pytest.fixture
def who_workbook(make_who_workbook) -> Path:
    """The default workbook: published layout, every trap present."""
    return make_who_workbook()


@pytest.fixture
def who_catalogue(who_workbook):
    """The default workbook, loaded.

    ``strict=False`` because a seven-row fixture cannot match the published
    48,152; the count mismatch lands in ``discrepancies``, which several tests
    then read. Structure is still enforced — ``strict`` governs counts only.
    """
    from mjolnir.resistance import catalogues

    return catalogues.load_who(who_workbook, strict=False)


# ---------------------------------------------------------------------------
# Variant and catalogue-call builders
# ---------------------------------------------------------------------------

def make_variant(gene: str = "rpoB", hgvs: str = "p.Ser450Leu", *,
                 chrom: str = "NC_000962.3", pos: int = 761155,
                 ref: str = "C", alt: str = "T", effect: str = "missense_variant",
                 calls: Optional[Sequence[CatalogueCall]] = None,
                 **kwargs: Any) -> Variant:
    """A :class:`~mjolnir.records.Variant` with sensible H37Rv-shaped defaults."""
    variant = Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, gene=gene,
                      hgvs=hgvs, effect=effect, **kwargs)
    if calls:
        variant.catalogue_calls.extend(calls)
    return variant


def who_call(drug: str, grade: str, variant_key: str, *, comment: str = "",
             matched_by: str = "coordinate") -> CatalogueCall:
    """A WHO row for one (drug, variant), with the call its grade implies."""
    return CatalogueCall(
        catalogue=config.CATALOGUE_WHO, drug=drug, grade=grade,
        call=config.call_for_grade(grade), variant_key=variant_key,
        comment=comment, matched_by=matched_by,
        catalogue_version="WHO-UCN-TB-2023.7")


def other_call(catalogue: str, drug: str, call: str, variant_key: str,
               *, grade: str = "") -> CatalogueCall:
    """A non-WHO catalogue's statement. MTBseq has no grade, by construction."""
    return CatalogueCall(catalogue=catalogue, drug=drug, call=call,
                         grade=grade, variant_key=variant_key, matched_by="hgvs")


@pytest.fixture
def variant_factory():
    return make_variant


@pytest.fixture
def call_factory() -> Dict[str, Any]:
    return {"who": who_call, "other": other_call, "none": CALL_NO_CALL}
