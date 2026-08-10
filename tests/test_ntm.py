"""NTM resistance (§5.6): three genes, named citations, and a refusal to guess.

None of the three MTBC catalogues covers non-tuberculous mycobacteria, so these
calls are implemented from the primary literature — *erm(41)* sequevar typing,
*rrl* 2058/2059 and *rrs* 1408 — and every one of them names its reference. The
clinically dangerous move here is not a wrong call but a confident silence: a
species-drug pair with no published evidence base must be reported as having
none, and an *erm(41)* that was never typed must not read as an isolate without
inducible resistance.
"""

from __future__ import annotations

import pytest

from mjolnir import config
from mjolnir.records import CALL_NO_CALL, CALL_R, CALL_S, Variant
from mjolnir.resistance import ntm


def _call(species, drug, variants=(), **kwargs):
    return ntm.call_drug(species, drug, variants, **kwargs)


# --------------------------------------------------- no evidence base is a finding

def test_a_pair_with_no_evidence_base_says_so_and_is_never_susceptible():
    call, hits, checks = _call("Mycobacterium chimaera", "Linezolid")
    assert ntm.is_no_evidence_base(call)
    assert call.call == CALL_NO_CALL
    assert call.call != CALL_S
    assert "not a prediction of susceptibility" in call.note
    assert checks[0].measured is False


def test_an_organism_outside_the_table_is_answered_the_same_way():
    call, _, _ = _call("Mycobacterium gordonae", "Amikacin")
    assert ntm.is_no_evidence_base(call)
    assert "Mycobacterium gordonae" in call.note


def test_supported_pairs_are_listed_explicitly():
    assert set(ntm.supported_drugs("Mycobacteroides abscessus")) == {
        "Clarithromycin", "Amikacin"}
    assert ntm.evidence_for("Mycobacterium chimaera", "Clarithromycin") is not None
    assert ntm.evidence_for("Mycobacterium chimaera", "Rifampicin") is None


# --------------------------------------------------------- erm(41) sequevar typing

def test_an_untyped_erm41_is_not_an_absence_of_inducible_resistance():
    call, hits, checks = _call("Mycobacteroides abscessus", "Clarithromycin")
    named = dict((c.name, c) for c in checks)
    assert named["erm(41) sequevar typed"].measured is False
    assert "neither detected nor excluded" in named["erm(41) sequevar typed"].reading
    assert call.call == CALL_NO_CALL
    assert call.call != CALL_S
    assert hits == []


def test_t28_gives_inducible_macrolide_resistance():
    observation = ntm.Erm41Observation(sequevar_base="T", present=True)
    assert observation.state == ntm.STATE_INDUCIBLE
    call, hits, _ = _call("Mycobacteroides abscessus", "Clarithromycin",
                          erm41=observation)
    assert call.call == CALL_R
    assert any("erm(41)" in hit.determinant.gene for hit in hits)
    assert hits[0].determinant.citation


def test_c28_is_the_susceptible_sequevar_but_does_not_clear_the_drug():
    """The T28C polymorphism removes one mechanism; acquired *rrl* is the other."""
    observation = ntm.Erm41Observation(sequevar_base="C", present=True)
    assert observation.state == ntm.STATE_SEQUEVAR_SUSCEPTIBLE
    call, hits, _ = _call("Mycobacteroides abscessus", "Clarithromycin",
                          erm41=observation)
    assert call.call != CALL_R
    assert call.call != CALL_S, "rrl was not assessed, so nothing here is susceptible"
    assert any("rrl was not assessed" in caveat for caveat in call.caveats)


def test_a_truncated_erm41_is_the_massiliense_state():
    observation = ntm.Erm41Observation(present=True, truncated=True, deletion_bp=274)
    assert observation.state == ntm.STATE_TRUNCATED


def test_an_unassessed_erm41_defaults_to_nothing_seen():
    """A default of "present and full length" would be an invention."""
    assert ntm.Erm41Observation().state == ntm.STATE_NOT_ASSESSED
    assert ntm.Erm41Observation().assessed is False


def test_a_base_that_is_neither_sequevar_is_not_rounded_to_the_nearest_one():
    assert ntm.Erm41Observation(sequevar_base="A",
                                present=True).state == ntm.STATE_NOT_ASSESSED


def test_the_sequevar_base_must_be_a_nucleotide():
    from mjolnir.utils import MjolnirError

    with pytest.raises(MjolnirError):
        ntm.Erm41Observation(sequevar_base="T28C")


# ------------------------------------------------------------- rrl and rrs

def _rrna(gene, pos, ref="A", alt="G"):
    return Variant(chrom=gene, pos=pos, ref=ref, alt=alt, gene=gene,
                   hgvs="n.{0}{1}>{2}".format(pos, ref, alt))


@pytest.mark.parametrize("position", list(config.RRL_MACROLIDE_POSITIONS))
def test_rrl_2058_and_2059_confer_macrolide_resistance(position):
    call, hits, _ = _call("Mycobacterium chimaera", "Clarithromycin",
                          [_rrna("rrl", position)])
    assert call.call == CALL_R
    assert hits[0].determinant.positions == config.RRL_MACROLIDE_POSITIONS
    assert hits[0].determinant.numbering == ntm.ECOLI_NUMBERING


def test_rrs_1408_confers_amikacin_resistance():
    call, hits, _ = _call("Mycobacteroides abscessus", "Amikacin",
                          [_rrna("rrs", 1408)])
    assert call.call == CALL_R
    assert 1408 in hits[0].determinant.positions


def test_an_examined_gene_with_nothing_in_it_is_still_not_susceptible():
    call, hits, checks = _call("Mycobacterium chimaera", "Amikacin",
                               callable_genes=["rrs"])
    named = dict((c.name, c) for c in checks)
    assert named["rrs assessed for Amikacin"].status == "pass"
    assert call.call == CALL_NO_CALL
    assert call.call != CALL_S
    assert any("no resistance determinant detected" in caveat for caveat in call.caveats)


def test_an_unexamined_gene_is_recorded_as_unexamined():
    call, _, checks = _call("Mycobacterium chimaera", "Amikacin")
    named = dict((c.name, c) for c in checks)
    assert named["rrs assessed for Amikacin"].measured is False


def test_a_neighbouring_rrs_position_is_uncertain_rather_than_resistant():
    """The evidence base for 1406 and 1409 is thinner than for 1408."""
    call, hits, _ = _call("Mycobacteroides abscessus", "Amikacin",
                          [_rrna("rrs", 1406)])
    assert call.call != CALL_R
    assert call.call in ("Uncertain", CALL_NO_CALL)


# ------------------------------------------------------ the ONT validation gap

def test_ont_ntm_calls_carry_the_not_validated_caveat():
    """No published R10.4.1-era ONT validation of NTM genotypic DST exists."""
    call, _, _ = _call("Mycobacterium chimaera", "Clarithromycin",
                       [_rrna("rrl", 2058)], platform="ont")
    assert call.call == CALL_R, "the caveat qualifies the call, it does not remove it"
    assert any("ONT" in caveat and "validation" in caveat for caveat in call.caveats)
    assert any("homopolymer" in caveat for caveat in call.caveats)


def test_illumina_ntm_calls_carry_no_platform_caveat():
    call, _, _ = _call("Mycobacterium chimaera", "Clarithromycin",
                       [_rrna("rrl", 2058)], platform="illumina")
    assert not any("ONT" in caveat for caveat in call.caveats)


def test_the_single_rrn_operon_note_travels_with_every_rrna_call():
    call, _, _ = _call("Mycobacterium chimaera", "Clarithromycin",
                       [_rrna("rrl", 2058)])
    assert any("rrn operon" in caveat for caveat in call.caveats)


# ------------------------------------------------------------ the whole assessment

def test_the_assessment_names_every_pair_it_could_not_answer():
    assessment = ntm.call_ntm_resistance(
        "Mycobacterium chimaera", [_rrna("rrl", 2058)],
        drugs=["Clarithromycin", "Linezolid"])
    assert assessment.drug("Clarithromycin").call == CALL_R
    assert [entry["drug"] for entry in assessment.no_evidence_base] == ["Linezolid"]
    assert assessment.citations


def test_every_determinant_carries_a_citation_and_declares_its_verification():
    for row in ntm.evidence_rows():
        assert row["citation"], row


def test_unverified_determinants_are_listed_rather_than_marked_settled():
    for determinant in ntm.unverified_determinants():
        assert determinant.verified is False
