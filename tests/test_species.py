"""Species identification must not come from a taxonomic read classifier (§6).

In current NCBI taxonomy the MTBC members are not at species rank at all:
*Mycobacterium tuberculosis* variant *bovis* (taxid 1765) has rank ``no rank``
under species *M. tuberculosis*, the members being later heterotypic synonyms at
99.21-99.92% ANI. A Kraken2 row saying "M. bovis 3.2%" is therefore not a species
identification, and the design's requirement is stronger than "do not print one":
it must not be reachable.

These tests pin the two halves of that. The call site somebody would naturally
write — classifier output into the typing layer — raises with the taxonomy
explanation, and every legitimate use of a classifier label collapses an MTBC
member to the complex on the way out.
"""

from __future__ import annotations

import pytest

from mjolnir import config
from mjolnir.contamination import purity
from mjolnir.typing import species
from mjolnir.utils import MjolnirError


# ------------------------------------------------- the refusal, at the call site

def test_species_from_a_classifier_raises_rather_than_answering():
    with pytest.raises(MjolnirError) as excinfo:
        species.species_from_classifier("Mycobacterium bovis", 3.2)
    message = str(excinfo.value)
    assert "refused" in message
    assert "taxid 1765" in message
    assert "identify_species" in message, "the refusal must name the supported route"


def test_the_refusal_explains_the_taxonomy_rather_than_asserting_a_policy():
    with pytest.raises(MjolnirError) as excinfo:
        species.species_from_classifier()
    message = str(excinfo.value)
    assert "heterotypic synonyms" in message
    assert "99.21-99.92% ANI" in message


# --------------------------------------------- what a classifier label may become

@pytest.mark.parametrize("label", [
    "Mycobacterium bovis",
    "Mycobacterium bovis BCG",
    "Mycobacterium tuberculosis",
    "Mycobacterium tuberculosis variant bovis",
    "Mycobacterium africanum",
    "Mycobacterium orygis",
    "Mycobacterium canettii",
])
def test_no_mtbc_member_survives_demotion_under_its_own_name(label):
    demoted = species.demote_classifier_label(label)
    assert demoted != label
    assert "complex" in demoted.lower()
    assert "not resolvable by a read classifier" in demoted


def test_a_non_mtbc_label_passes_through_untouched():
    """The refusal is specific, not a general distrust of taxon names."""
    assert species.demote_classifier_label("Cutibacterium acnes") == "Cutibacterium acnes"


def test_mac_members_are_demoted_too_because_ani_cannot_separate_them():
    demoted = species.demote_classifier_label("Mycobacterium chimaera")
    assert "avium complex" in demoted


def test_complex_placement_walks_down_to_the_binomial():
    assert species.complex_for(
        "Mycobacterium avium subsp. hominissuis TH135") == config.COMPLEX_MAC
    assert species.complex_for("Mycobacteroides abscessus subsp. massiliense") == \
        config.COMPLEX_ABSCESSUS
    assert species.complex_for("Escherichia coli") == ""


@pytest.mark.parametrize("name,expected", [
    ("Mycobacterium tuberculosis variant bovis", True),
    ("Mycobacterium bovis BCG", True),
    ("Mycobacterium avium", False),
    ("", False),
])
def test_is_mtbc_member(name, expected):
    assert species.is_mtbc_member(name) is expected


# ------------------------------------- the same refusal on the contamination side

@pytest.mark.parametrize("label", [
    "Mycobacterium tuberculosis",
    "Mycobacterium tuberculosis variant bovis BCG",
    "Mycobacterium orygis",
    "Mycobacterium tuberculosis complex",
])
def test_a_taxon_label_leaves_the_screen_as_the_complex(label):
    assert purity.taxon_label_for_report(label) == purity.MTBC_CLASSIFIER_LABEL


def test_a_kraken_report_cannot_emit_an_mtbc_member_name():
    """Parsed at the boundary, so no caller can print a member name downstream."""
    report = (
        "88.10\t8810\t0\tG\t1763\t  Mycobacterium\n"
        "80.00\t8000\t8000\tS\t1773\t    Mycobacterium tuberculosis\n"
        "3.20\t320\t320\tS\t1765\t    Mycobacterium tuberculosis variant bovis\n"
        "0.40\t40\t40\tS\t1717\t  Corynebacterium diphtheriae\n"
    )
    rows = purity.parse_kraken2_report(report)
    labels = [row["label"] for row in rows]
    assert "Mycobacterium tuberculosis" not in labels
    assert "Mycobacterium tuberculosis variant bovis" not in labels
    assert purity.MTBC_CLASSIFIER_LABEL in labels
    assert "Corynebacterium diphtheriae" in labels
    assert [r["collapsed_to_complex"] for r in rows] == [False, True, True, False]


def test_the_species_method_refusal_is_a_registered_threshold():
    """It is printed in the report, so it is sourced like every other number."""
    assert config.source_for("species_method_refusal") == config.SRC_DESIGN
    assert "never printed as one" in config.SPECIES_METHOD_REFUSAL
