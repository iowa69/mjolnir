"""The seven WHO catalogue traps of design §5.2, each pinned by a test.

Every one of these is a *verified* failure mode rather than a hypothetical, and
every one of them produces a report that looks correct. A loader that reads the
header from row 1 gets garbage column names and then finds no drugs; a loader
that reads the ``.txt`` master file loses streptomycin entirely and mis-grades
four genes; a loader that takes the ``genomic position`` column at face value
matches nothing; a loader that deduplicates by variant name reports isoniazid's
grade for ethionamide. None of those failures raises on its own, which is why
they are tested here rather than trusted to code review.
"""

from __future__ import annotations

import pytest

from conftest import (DECOY, GRADE_2_ENDASH, who_master_rows, write_who_workbook)
from mjolnir import config
from mjolnir.resistance import catalogues
from mjolnir.utils import MjolnirError


# --------------------------------------------------------------- trap 1: row 3

def test_header_is_read_from_row_three(who_catalogue):
    """The header is on row 3; rows 1-2 are a merged banner.

    If the loader had taken row 1 as the header the columns would be the banner
    text and no drug would have been recognised at all.
    """
    assert "Isoniazid" in who_catalogue.drugs
    assert "Rifampicin" in who_catalogue.drugs
    assert who_catalogue.entries, "no rows survived the header parse"
    assert config.WHO_XLSX_HEADER_ROW == 3


def test_header_on_row_one_is_refused_and_says_where_the_header_lives(make_who_workbook):
    workbook = make_who_workbook(name="WHO-UCN-TB-2023.7-eng.xlsx", header_row=1)
    with pytest.raises(MjolnirError) as excinfo:
        catalogues.load_who(workbook, strict=False)
    message = str(excinfo.value)
    assert "row 3" in message
    assert "banner" in message
    # It has to name what it could not find, or the reader looks in Excel for a
    # column that is there and concludes the tool is broken.
    assert "'drug'" in message


def test_first_two_rows_are_reported_back_when_row_three_is_not_a_header(make_who_workbook):
    workbook = make_who_workbook(header_row=1)
    with pytest.raises(MjolnirError) as excinfo:
        catalogues.load_who(workbook, strict=False)
    assert "Row 3 read as" in str(excinfo.value)


# ------------------------------------------------- trap 2: the .txt master file

@pytest.mark.parametrize("name", [
    "WHO-UCN-TB-2023.7-eng_catalogue_master_file.txt",
    "WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt",
    "catalogue_master_file.tsv",
    "master.txt",
])
def test_text_master_file_is_refused(tmp_path, name):
    """The plain-text master is not the workbook and there is no flag that says it is."""
    path = tmp_path / name
    path.write_text("drug\tgene\tvariant\n")
    with pytest.raises(MjolnirError) as excinfo:
        catalogues.load_who(path)
    message = str(excinfo.value)
    assert "Streptomycin" in message, "the refusal must name the drug that goes missing"
    assert "40,178" in message and "48,152" in message
    assert ".xlsx" in message


def test_text_master_refusal_happens_before_the_file_is_opened(tmp_path):
    """A file that does not exist still gets the .txt refusal, not a not-found error."""
    with pytest.raises(MjolnirError) as excinfo:
        catalogues.refuse_who_text_master(
            tmp_path / "WHO-UCN-TB-2023.7-eng_catalogue_master_file.txt")
    assert "refusing to read the WHO catalogue from a plain-text master file" \
        in str(excinfo.value)


def test_xlsx_is_not_refused(who_workbook):
    assert catalogues.refuse_who_text_master(who_workbook) is None


def test_absent_streptomycin_is_the_signature_of_the_text_file(make_who_workbook):
    """A workbook with 14 drugs and no streptomycin is refused under strict loading.

    Not a row count — a signature. A run built on the ``.txt`` reports no
    streptomycin result at all, which must never look like a run that found
    streptomycin susceptible.
    """
    workbook = make_who_workbook(rows=who_master_rows(include_streptomycin=False))
    with pytest.raises(MjolnirError) as excinfo:
        catalogues.load_who(workbook, strict=True)
    assert catalogues.WHO_TXT_SIGNATURE_DRUG in str(excinfo.value)
    assert catalogues.WHO_TXT_SIGNATURE_DRUG == "Streptomycin"


def test_strict_false_still_records_the_missing_drug(make_who_workbook):
    """The newer-edition override demotes the check; it does not delete the reason."""
    workbook = make_who_workbook(rows=who_master_rows(include_streptomycin=False))
    catalogue = catalogues.load_who(workbook, strict=False)
    assert any("Streptomycin" in note for note in catalogue.discrepancies)


# ------------------------------------------- trap 3: the decoy position column

def test_decoy_position_column_never_becomes_a_coordinate(who_catalogue):
    """``eis_c.-14C>T`` carries a number in ``genomic position`` and no coordinate row.

    Its coordinates must be empty. If the master sheet's column were being read,
    this entry would have acquired ``NC_000962.3:2715342`` — a coordinate that
    exists nowhere in WHO's coordinate table.
    """
    entries = who_catalogue.lookup_key("eis_c.-14C>T")
    assert entries, "the eis row did not load at all"
    for entry in entries:
        assert entry.coordinates == ()
    assert not who_catalogue.lookup_coordinate(("NC_000962.3", 2715342, "C", "T"))


def test_coordinates_come_from_the_coordinates_sheet(who_catalogue):
    entries = who_catalogue.lookup_key("inhA_c.-154G>A")
    assert entries
    for entry in entries:
        assert entry.coordinates == (("NC_000962.3", 1673432, "G", "A"),)
    assert who_catalogue.lookup_coordinate(("NC_000962.3", 1673432, "G", "A"))


def test_no_coordinate_key_carries_the_decoy_string(who_catalogue):
    for entry in who_catalogue.entries:
        for chrom, pos, ref, alt in entry.coordinates:
            assert DECOY not in "{0}{1}{2}".format(chrom, ref, alt)
            assert isinstance(pos, int)


def test_a_decoy_string_in_the_coordinate_table_is_fatal(tmp_path):
    """The decoy must not be parsed as a position even when it appears in the table."""
    path = tmp_path / "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
    path.write_text(
        "variant\tchromosome\tposition\treference_nucleotide\talternative_nucleotide\n"
        "inhA_c.-154G>A\tNC_000962.3\t{0}\tG\tA\n".format(DECOY))
    with pytest.raises(MjolnirError) as excinfo:
        catalogues.read_who_coordinates_file(path)
    assert "non-numeric position" in str(excinfo.value)
    assert config.WHO_XLSX_COORDINATES_SHEET in str(excinfo.value)


# -------------------------------------------------------- trap 4: grade strings

def test_canonical_grades_use_a_spaced_ascii_hyphen():
    for grade in config.WHO_GRADES:
        assert "–" not in grade, "en-dash in a canonical grade string"
        assert "—" not in grade
    assert config.WHO_GRADE_2 == "2) Assoc w R - Interim"
    assert config.WHO_GRADE_4 == "4) Not assoc w R - Interim"


def test_the_pdf_en_dash_form_is_not_one_of_the_five_grades():
    """Hard-coding the PDF's spelling matches nothing — which is the trap."""
    assert GRADE_2_ENDASH not in config.WHO_GRADES
    assert config.WHO_GRADE_TO_CALL.get(GRADE_2_ENDASH) is None


def test_en_dash_grade_normalises_onto_the_ascii_form():
    assert config.normalise_grade(GRADE_2_ENDASH) == config.WHO_GRADE_2
    assert config.call_for_grade(GRADE_2_ENDASH) == config.WHO_GRADE_TO_CALL[
        config.WHO_GRADE_2]


def test_a_row_spelled_with_an_en_dash_is_stored_as_the_ascii_form(who_catalogue):
    entry = who_catalogue.lookup_key("rpoB_p.Ser450Leu")[0]
    assert entry.grade == config.WHO_GRADE_2
    assert entry.call == "R-interim"


def test_an_unrecognised_grade_produces_no_call_rather_than_a_guess(make_who_workbook):
    rows = who_master_rows()
    rows.append(["Linezolid", "rrl", "n.2814G>T", "rrl_n.2814G>T", 2,
                 "non_coding_transcript_variant", DECOY, "probably resistant",
                 "", "", ""])
    catalogue = catalogues.load_who(make_who_workbook(rows=rows), strict=False)
    entry = catalogue.lookup_key("rrl_n.2814G>T")[0]
    assert entry.call == "no-call"
    assert any("probably resistant" in note for note in catalogue.discrepancies)


# ------------------------------------------------- trap 5: grading per (drug, variant)

def test_the_same_variant_carries_two_grades_for_two_drugs(who_catalogue):
    """``inhA_c.-154G>A`` is Group 1 for isoniazid and Group 2 for ethionamide.

    Deduplicating by variant name would have to pick one, and would be wrong for
    whichever drug it discarded.
    """
    grades = dict((e.drug, e.grade) for e in who_catalogue.lookup_key("inhA_c.-154G>A"))
    assert grades["Isoniazid"] == config.WHO_GRADE_1
    assert grades["Ethionamide"] == config.WHO_GRADE_2


def test_multi_graded_variants_surfaces_the_split(who_catalogue):
    multi = who_catalogue.multi_graded_variants()
    assert "inhA_c.-154G>A" in multi
    assert set(multi["inhA_c.-154G>A"]) == {"Isoniazid", "Ethionamide"}


def test_the_two_drug_rows_produce_different_calls(who_catalogue):
    calls = dict((e.drug, e.call) for e in who_catalogue.lookup_key("inhA_c.-154G>A"))
    assert calls["Isoniazid"] == "R"
    assert calls["Ethionamide"] == "R-interim"


# ------------------------------------------------------- trap 6: MNV `&` splitting

def test_mnv_graded_variants_split_on_the_ampersand(who_catalogue):
    """One genomic change, two graded variants, one coordinate each.

    Without the split the coordinate would be indexed under the joined string
    ``rpoB_p.Ser450Leu&rpoB_p.Leu452Pro``, which no master row carries, and both
    variants would be reported as ungraded.
    """
    coordinate = ("NC_000962.3", 761155, "C", "T")
    for key in ("rpoB_p.Ser450Leu", "rpoB_p.Leu452Pro"):
        entries = who_catalogue.lookup_key(key)
        assert entries, "{0} did not load".format(key)
        assert coordinate in entries[0].coordinates
    graded = set(e.variant_key for e in who_catalogue.lookup_coordinate(coordinate))
    assert graded == {"rpoB_p.Ser450Leu", "rpoB_p.Leu452Pro"}


def test_the_joined_name_is_never_itself_a_key(who_catalogue):
    assert not who_catalogue.lookup_key("rpoB_p.Ser450Leu&rpoB_p.Leu452Pro")
    assert config.WHO_MNV_SEPARATOR == "&"


def test_coordinate_file_splits_mnvs_too(tmp_path):
    path = tmp_path / "coords.txt"
    path.write_text(
        "variant\tchromosome\tposition\treference_nucleotide\talternative_nucleotide\n"
        "rpoB_p.Ser450Leu&rpoB_p.Leu452Pro\tNC_000962.3\t761155\tC\tT\n")
    coordinates = catalogues.read_who_coordinates_file(path)
    assert coordinates["rpoB_p.Ser450Leu"] == [("NC_000962.3", 761155, "C", "T")]
    assert coordinates["rpoB_p.Leu452Pro"] == [("NC_000962.3", 761155, "C", "T")]


def test_vcf_coordinates_split_on_the_ampersand_as_well(tmp_path):
    path = tmp_path / "Genomic_coordinates.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "NC_000962.3\t761155\t.\tC\tT\t.\t.\t"
        "graded_variant=rpoB_p.Ser450Leu&rpoB_p.Leu452Pro\n")
    coordinates = catalogues.read_who_coordinates_file(path)
    assert set(coordinates) == {"rpoB_p.Ser450Leu", "rpoB_p.Leu452Pro"}


# ------------------------------------------- trap 7: pooled names have no coordinates

def test_pooled_lof_names_are_matched_by_rule_not_by_coordinate(who_catalogue):
    entry = who_catalogue.lookup_key("Rv0678_LoF")[0]
    assert entry.rule_only is True
    assert entry.coordinates == ()


# ---------------------------------------------------------- structural failures

def test_a_workbook_without_the_master_sheet_is_not_the_catalogue(tmp_path):
    from openpyxl import Workbook

    path = tmp_path / "WHO-UCN-TB-2023.7-eng.xlsx"
    workbook = Workbook()
    workbook.active.title = "Sheet1"
    workbook.save(str(path))
    with pytest.raises(MjolnirError) as excinfo:
        catalogues.load_who(path, strict=False)
    assert config.WHO_XLSX_MASTER_SHEET in str(excinfo.value)


def test_duplicate_column_names_are_reported_not_silently_collapsed(who_catalogue):
    """The published master sheet carries ``CHANGES vs ver1`` twice."""
    assert any("CHANGES vs ver1" in note for note in who_catalogue.discrepancies)


def test_a_missing_workbook_says_how_to_fetch_it(tmp_path):
    with pytest.raises(MjolnirError) as excinfo:
        catalogues.load_who(tmp_path / "absent.xlsx")
    assert "mjolnir db fetch" in str(excinfo.value)


def test_catalogue_carries_its_version_and_checksum(who_catalogue):
    """Two installations that disagree must be tellable apart from the report."""
    assert who_catalogue.version == "WHO-UCN-TB-2023.7"
    assert len(who_catalogue.checksum) == 64
    version = who_catalogue.database_version()
    assert version.licence == catalogues.LICENCE_WHO
    assert version.checksum == who_catalogue.checksum
