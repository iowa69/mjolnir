"""WHO's additional grading rules (§5.4), which the table alone cannot supply.

The catalogue is three components, not one table, so a lookup-only implementation
is wrong in a specific and clinically loaded direction: it under-calls. The four
borderline *rpoB* mutations, the RRDR rule, the loss-of-function rules and the
epistasis suppressions are all things the table does not say and the catalogue
does.

The epistasis tests are the ones to read twice. *mmpL5* loss of function
abrogates *Rv0678* variants for bedaquiline and clofazimine, and the design
requires the suppression to be **recorded** rather than applied silently — a
suppressed determinant that vanishes from the output is indistinguishable from a
determinant that was never there.
"""

from __future__ import annotations

import pytest

from conftest import make_variant
from mjolnir import config
from mjolnir.resistance import rules


# ----------------------------------------------- the four borderline rpoB calls

@pytest.mark.parametrize("mutation", list(config.RPOB_BORDERLINE))
def test_borderline_rpob_mutations_are_group_one_by_rule(mutation):
    variant = make_variant(gene="rpoB", hgvs="p." + mutation)
    hits = rules.rule_hits(variant)
    assert len(hits) == 1
    hit = hits[0]
    assert hit.rule == rules.RULE_RPOB_BORDERLINE
    assert hit.drug == "Rifampicin"
    assert hit.grade == config.WHO_GRADE_1
    assert hit.call == "R"


def test_there_are_exactly_four_borderline_mutations():
    assert config.RPOB_BORDERLINE == (
        "Leu430Pro", "His445Asn", "His445Ser", "Ile491Phe")


def test_ile491phe_is_outside_the_rrdr_and_still_group_one():
    """The rule exists precisely because an RRDR-bounded implementation misses it."""
    variant = make_variant(gene="rpoB", hgvs="p.Ile491Phe")
    assert rules.is_in_rpob_rrdr(variant) is False
    assert rules.rule_hits(variant)[0].grade == config.WHO_GRADE_1


def test_borderline_rule_does_not_fire_on_another_gene():
    assert rules.rpob_borderline_hit(make_variant(gene="katG", hgvs="p.Ile491Phe")) == ""


def test_borderline_takes_precedence_over_the_rrdr_rule():
    """Leu430Pro is inside the RRDR; Group 1 by rule beats Group 2 by region."""
    variant = make_variant(gene="rpoB", hgvs="p.Leu430Pro")
    assert rules.is_in_rpob_rrdr(variant) is True
    grades = [h.grade for h in rules.rule_hits(variant)]
    assert grades == [config.WHO_GRADE_1]


# ------------------------------------------------------------------- the RRDR

@pytest.mark.parametrize("codon,inside", [
    (425, False), (426, True), (435, True), (452, True), (453, False)])
def test_rrdr_bounds_are_inclusive(codon, inside):
    variant = make_variant(gene="rpoB", hgvs="p.Asp{0}Val".format(codon))
    assert rules.is_in_rpob_rrdr(variant) is inside


def test_non_synonymous_rrdr_change_is_group_two():
    hits = rules.rule_hits(make_variant(gene="rpoB", hgvs="p.Asp435Val"))
    assert [h.rule for h in hits] == [rules.RULE_RPOB_RRDR]
    assert hits[0].grade == config.WHO_GRADE_2


def test_a_silent_rrdr_change_does_not_grade_group_two():
    variant = make_variant(gene="rpoB", hgvs="p.Asp435=", effect="synonymous_variant")
    assert rules.RULE_RPOB_RRDR not in [h.rule for h in rules.rule_hits(variant)]


def test_rules_stay_silent_where_the_table_already_grades_the_variant():
    """Additional grading rules are additional, not overriding."""
    variant = make_variant(gene="rpoB", hgvs="p.Asp435Val")
    assert rules.rule_hits(variant, who_graded_drugs=["Rifampicin"]) == []


# ------------------------------------------------------------ loss of function

@pytest.mark.parametrize("gene,drugs", [
    ("katG", ("Isoniazid",)),
    ("pncA", ("Pyrazinamide",)),
    ("Rv0678", ("Bedaquiline", "Clofazimine")),
    ("pepQ", ("Bedaquiline", "Clofazimine")),
    ("ethA", ("Ethionamide",)),
])
def test_loss_of_function_grades_group_two_for_the_named_drugs(gene, drugs):
    variant = make_variant(gene=gene, hgvs="p.Trp90*", effect="stop_gained")
    hits = [h for h in rules.rule_hits(variant) if h.rule == rules.RULE_LOF]
    assert set(h.drug for h in hits) == set(drugs)
    assert set(h.grade for h in hits) == {config.WHO_GRADE_2}


def test_an_inframe_deletion_is_not_a_loss_of_function():
    """Length alone is not evidence: a one-codon deletion is usually not LoF."""
    variant = make_variant(gene="pncA", hgvs="p.Val130del", ref="CGTG", alt="C",
                           effect="inframe_deletion", variant_type="deletion")
    assert rules.is_loss_of_function(variant) is False
    assert rules.rule_hits(variant) == []


def test_pooled_lof_names_are_recognised():
    assert rules.is_loss_of_function(make_variant(gene="katG", hgvs="katG_LoF"))
    assert rules.is_loss_of_function(make_variant(gene="pncA", hgvs="pncA_deletion"))


# ------------------------------------------------------------- silent variants

def test_a_novel_silent_variant_is_group_four_for_the_genes_drugs():
    variant = make_variant(gene="katG", hgvs="p.Arg463=", effect="synonymous_variant")
    hits = rules.rule_hits(variant, gene_drugs=["Isoniazid"])
    assert [h.rule for h in hits] == [rules.RULE_SILENT]
    assert hits[0].grade == config.SILENT_VARIANT_GRADE == config.WHO_GRADE_4


def test_a_silent_variant_in_an_unassociated_gene_produces_nothing():
    """Mjolnir does not guess which drug a gene belongs to."""
    variant = make_variant(gene="Rv1258c", hgvs="p.Ala100=", effect="synonymous_variant")
    assert rules.rule_hits(variant) == []


def test_silence_is_a_positive_finding_not_a_missing_annotation():
    unannotated = make_variant(gene="katG", hgvs="c.1388G>T", effect="")
    assert rules.is_synonymous(unannotated) is False


# --------------------------------------------------------- the Comment column

def test_level_of_resistance_comes_from_the_comment_not_the_grade():
    assert rules.level_from_comment("High-level resistance") == rules.LEVEL_HIGH
    assert rules.level_from_comment("Low-level resistance") == rules.LEVEL_LOW
    assert rules.level_from_comment("") == ""


def test_an_unstated_level_is_not_softened_to_low():
    assert rules.level_from_comment("Associated with resistance in vitro") == ""


def test_cross_resistance_needs_a_cue_not_just_a_drug_name():
    assert rules.cross_resistance_from_comment(
        "graded separately for pretomanid", "Delamanid") == []
    assert rules.cross_resistance_from_comment(
        "cross-resistance with pretomanid is expected", "Delamanid") == ["Pretomanid"]


# ------------------------------------------------------------------ epistasis

def _bdq_variants(mmpl5_fraction=None, rv0678_fraction=None):
    """An *mmpL5* knockout beside an *Rv0678* mutation, the classic pattern."""
    mmpl5 = make_variant(gene="mmpL5", hgvs="p.Ala100fs", pos=778990, ref="G",
                         alt="GA", effect="frameshift_variant",
                         variant_type="insertion", allele_fraction=mmpl5_fraction)
    rv0678 = make_variant(gene="Rv0678", hgvs="p.Ser53Leu", pos=779000,
                          effect="missense_variant", allele_fraction=rv0678_fraction)
    return mmpl5, rv0678


def test_mmpl5_lof_suppresses_rv0678_for_bedaquiline_and_clofazimine():
    suppressions = rules.epistasis_suppressions(list(_bdq_variants()))
    assert set(s.drug for s in suppressions) == {"Bedaquiline", "Clofazimine"}
    for suppression in suppressions:
        assert suppression.suppressor_gene == "mmpL5"
        assert suppression.suppressed_gene == "Rv0678"
        assert suppression.suppresses("Rv0678_p.Ser53Leu")
        assert suppression.confident is True
        assert "mmpL5" in suppression.why


def test_the_suppression_is_a_record_with_both_sides_named():
    """The abrogated variant has to survive into the output to be reportable."""
    suppression = rules.epistasis_suppressions(list(_bdq_variants()))[0]
    assert suppression.suppressor_variants == ["mmpL5_p.Ala100fs"]
    assert suppression.suppressed_variants == ["Rv0678_p.Ser53Leu"]
    assert suppression.source
    assert suppression.to_dict()["rule"].startswith("epistasis:mmpL5")


def test_no_suppression_without_the_mmpl5_knockout():
    _, rv0678 = _bdq_variants()
    assert rules.epistasis_suppressions([rv0678]) == []


def test_an_mmpl5_missense_does_not_abrogate_anything():
    """The rule is about loss of function, not about any variant in the gene."""
    mmpl5 = make_variant(gene="mmpL5", hgvs="p.Ile948Val", effect="missense_variant")
    _, rv0678 = _bdq_variants()
    assert rules.epistasis_suppressions([mmpl5, rv0678]) == []


def test_a_split_subpopulation_is_reported_but_not_applied():
    """A knockout in part of the population cannot be assumed to abrogate the rest."""
    variants = list(_bdq_variants(mmpl5_fraction=0.30, rv0678_fraction=0.98))
    suppressions = rules.epistasis_suppressions(variants)
    assert suppressions
    for suppression in suppressions:
        assert suppression.confident is False
        assert suppression.caveat == rules.MIXED_SUBPOPULATION_CAVEAT


def test_eis_coding_lof_abrogates_the_eis_promoter_for_amikacin_and_kanamycin():
    coding = make_variant(gene="eis", hgvs="p.Gln50*", effect="stop_gained")
    promoter = make_variant(gene="eis", hgvs="c.-14C>T",
                            effect="upstream_gene_variant")
    suppressions = rules.epistasis_suppressions([coding, promoter])
    assert set(s.drug for s in suppressions) == {"Amikacin", "Kanamycin"}
    for suppression in suppressions:
        assert suppression.suppressed_variants == ["eis_c.-14C>T"]


def test_an_eis_promoter_variant_alone_suppresses_nothing():
    promoter = make_variant(gene="eis", hgvs="c.-14C>T",
                            effect="upstream_gene_variant")
    assert rules.epistasis_suppressions([promoter]) == []


def test_epistasis_can_be_restricted_to_a_drug_panel():
    suppressions = rules.epistasis_suppressions(list(_bdq_variants()),
                                                drugs=["Bedaquiline"])
    assert [s.drug for s in suppressions] == ["Bedaquiline"]
