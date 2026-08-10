"""The MTBseq list's drug column, and the rows a narrow vocabulary silently loses.

Measured against the shipped file (1,696 rows), three antibiotic labels resolved
to nothing and their rows never reached a drug column:

    22  fluoroquinolones (FQ)
    22  para-aminosalicylic acid (PAS)
     8  cycloserine (CS)

A determinant that reaches no drug column is indistinguishable in the report
from one that was never found, so losing them is not a cosmetic gap — it is the
tool reporting absence of evidence as evidence of absence, which is the single
thing this project is written to refuse.

Two of those labels are drugs the WHO catalogue v2 does not grade. The third is
a *class*, and expanding it to its members is right for the gyrA mutations
behind those rows but must be recorded, because a class-level source rendered
silently as two agent-level calls claims precision MTBseq never offered.

These tests build small synthetic tables rather than reading the shipped file:
they must pass with no database installed and no network.
"""

from __future__ import annotations

import pytest

from mjolnir import config
from mjolnir.resistance import catalogues
from mjolnir.utils import smart_open

MTBSEQ_HEADER = [
    "#Variant position genome start", "Variant position genome stop", "Var. type",
    "Number", "WT base", "Var. base", "Region", "Gene ID", "Gene Name", "Gene start",
    "Gene stop", "Gene length", "Dir.", "WT AA", "Codon nr.", "Codon nr. E. coli",
    "Var. AA", "AA change", "Codon change", "Variant position gene start",
    "Variant position gene stop", "Antibiotic", "Reference PMID",
    "High Confidence SNP", "Comment",
]


def _row(antibiotic, gene="gyrA", aa_change="Ala90Val", region="CDS",
         pos="7570", direction="+", comment="-", pmid="-", high_conf="-"):
    return [pos, pos, "SNP", "1", "C", "T", region, "Rv0006", gene, "5240", "7767",
            "2517", direction, "Ala", "90", "90", "Val", aa_change, "gcg/gtg", "269",
            "269", antibiotic, pmid, high_conf, comment]


@pytest.fixture
def mtbseq_list(tmp_path):
    """Write a small MTBseq-format table and return its path."""
    def _write(rows, encoding="latin-1"):
        path = tmp_path / "MTB_Resistance_Mediating.txt"
        with open(str(path), "w", encoding=encoding) as handle:
            handle.write("\t".join(MTBSEQ_HEADER) + "\n")
            for row in rows:
                handle.write("\t".join(row) + "\n")
        return path
    return _write


# ---------------------------------------------------------------------------
# The drugs a narrow vocabulary drops
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label,expected", [
    ("para-aminosalicylic acid (PAS)", "Para-aminosalicylic acid"),
    ("cycloserine (CS)", "Cycloserine"),
])
def test_drugs_who_does_not_grade_still_reach_a_drug_column(mtbseq_list, label,
                                                            expected):
    """PAS and cycloserine are real drugs, graded by MTBseq and not by WHO.

    They are not in the WHO v2 fifteen, which is exactly why they must survive:
    the only catalogue that says anything about them is this one.
    """
    catalogue = catalogues.load_mtbseq(mtbseq_list([_row(label)]))
    assert [e.drug for e in catalogue.entries] == [expected]


def test_a_drug_class_expands_to_its_members(mtbseq_list):
    """``fluoroquinolones (FQ)`` covers the fluoroquinolones the catalogues grade."""
    catalogue = catalogues.load_mtbseq(mtbseq_list([_row("fluoroquinolones (FQ)")]))
    assert sorted(e.drug for e in catalogue.entries) == ["Levofloxacin", "Moxifloxacin"]


def test_the_class_expansion_is_recorded_on_every_entry_it_produced(mtbseq_list):
    """The report must be able to say MTBseq graded a class, not these two agents.

    Without this the expansion is indistinguishable from MTBseq having named
    levofloxacin and moxifloxacin itself.
    """
    catalogue = catalogues.load_mtbseq(mtbseq_list([_row("fluoroquinolones (FQ)")]))
    assert catalogue.entries
    for entry in catalogue.entries:
        assert "fluoroquinolones" in entry.comment
        assert "expanded to" in entry.comment


def test_an_expansion_note_does_not_displace_the_row_s_own_comment(mtbseq_list):
    """MTBseq's comments carry MIC statements; the note is appended, not substituted."""
    catalogue = catalogues.load_mtbseq(
        mtbseq_list([_row("fluoroquinolones (FQ)", comment="moderate MIC increase")]))
    for entry in catalogue.entries:
        assert "moderate MIC increase" in entry.comment
        assert "expanded to" in entry.comment


def test_a_genuinely_unknown_label_is_skipped_rather_than_invented(mtbseq_list):
    """The vocabulary is widened, not abandoned.

    ``normalise_drug`` echoes an unknown name back, so accepting anything would
    turn a stray token into a phantom drug column. An unrecognised label must
    still produce no entry.
    """
    catalogue = catalogues.load_mtbseq(mtbseq_list([_row("not-a-drug (XYZ)")]))
    assert catalogue.entries == []


def test_phylogenetic_rows_never_reach_a_drug_call(mtbseq_list):
    """238 of the shipped rows are phylogenetic markers, not resistance."""
    catalogue = catalogues.load_mtbseq(
        mtbseq_list([_row("phylogenetic SNP"), _row("cycloserine (CS)")]))
    assert [e.drug for e in catalogue.entries] == ["Cycloserine"]


# ---------------------------------------------------------------------------
# Encoding
# ---------------------------------------------------------------------------

def test_the_list_is_read_as_latin_1(mtbseq_list):
    """The shipped file is latin-1 and raises as UTF-8 at byte 0x98, offset 7324.

    Reading it as UTF-8 with replacement does not crash — it corrupts the
    Comment column, which is where the MIC statements live, so the failure is
    silent and lands in the report.
    """
    comment = "MIC increase \x98 as reported"
    catalogue = catalogues.load_mtbseq(
        mtbseq_list([_row("cycloserine (CS)", comment=comment)]))
    assert catalogue.entries
    assert "�" not in catalogue.entries[0].comment
    assert "MIC increase" in catalogue.entries[0].comment


def test_smart_open_honours_an_explicit_encoding(tmp_path):
    """The byte that breaks UTF-8 decoding survives when the encoding is named."""
    path = tmp_path / "latin.txt"
    path.write_bytes("caf\x98\n".encode("latin-1"))
    with smart_open(path, "rt", encoding="latin-1") as handle:
        assert "�" not in handle.read()


# ---------------------------------------------------------------------------
# The vocabulary itself
# ---------------------------------------------------------------------------

def test_the_class_table_only_expands_to_drugs_the_tool_knows():
    """An expansion naming a drug outside the vocabulary would produce a phantom column."""
    for members in config.DRUG_CLASSES.values():
        for member in members:
            assert config.normalise_drug(member) == member
            assert member.lower() in config.DRUG_ALIASES


def test_parse_returns_the_class_it_expanded(mtbseq_list):
    """The third return value is what lets the loader record the expansion."""
    drugs, unresolved, expanded = catalogues._parse_mtbseq_drugs("fluoroquinolones (FQ)")
    assert sorted(drugs) == ["Levofloxacin", "Moxifloxacin"]
    assert unresolved == []
    assert expanded == ["fluoroquinolones"]


def test_parse_reports_an_unknown_label_as_unresolved():
    """Unresolved is not the same as absent: the loader counts these and logs them."""
    drugs, unresolved, expanded = catalogues._parse_mtbseq_drugs("moxifloxacin (MXF) wibble (WIB)")
    assert drugs == ["Moxifloxacin"]
    assert unresolved == ["wibble"]
    assert expanded == []
