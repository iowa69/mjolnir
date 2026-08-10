"""A SNP distance is a fraction, and the suite refuses to lose the denominator.

Twelve differences over 4.1 Mb of shared callable sequence and twelve over
400 kb are not the same statement, but a matrix that prints only the numerator
invites a reader to treat them as if they were — and the published 5-SNP and
12-SNP thresholds were derived over whole genomes, so a thin denominator makes
them inapplicable rather than merely uncertain.

The other half of §9 is masking. Distances counted without a mask include the
~6% of H37Rv that is repetitive, low-complexity or error-prone, which is exactly
where spurious differences accumulate; so an unmasked comparison is possible only
as an explicit decision that the report then prints.
"""

from __future__ import annotations

import pytest

from mjolnir.cohort import distance, joint
from mjolnir.records import Variant
from mjolnir.utils import MjolnirError

CHROM = "NC_000962.3"
LENGTH = 100_000


def _variant(pos, alt, allele_fraction=1.0):
    return Variant(chrom=CHROM, pos=pos, ref="A", alt=alt,
                   allele_fraction=allele_fraction,
                   is_major=allele_fraction >= 0.9)


def _sample(name, positions, callable_bp=LENGTH):
    regions = None
    if callable_bp:
        regions = joint.Regions(name).add(CHROM, 0, callable_bp)
    return joint.SampleVariants(
        sample_id=name,
        variants=[_variant(pos, alt) for pos, alt in positions],
        callable_regions=regions, reference=CHROM)


def _table(*samples):
    return joint.build_joint_table(list(samples), reference=CHROM)


UNMASKED = None  # built per test, since Mask.absent demands a stated reason


def _absent_mask():
    return distance.Mask.absent(
        reason="the test reference has no repeat mask; distances are provisional")


# ------------------------------------------------- no denominator, no distance

def test_a_pair_without_callable_regions_has_no_distance():
    table = _table(_sample("A", [(1000, "T"), (5000, "G")]),
                   _sample("B", [(1000, "T")], callable_bp=0))
    pair = distance.pairwise_distance(table, "A", "B", _absent_mask())
    assert pair.snps is None, "a pair with no denominator is not a pair at distance 0"
    assert pair.shared_callable_sites is None
    assert pair.snps_per_mb is None
    assert "no denominator" in pair.note or "denominator" in pair.note
    assert "B" in pair.note


def test_the_missing_denominator_note_says_how_to_get_one():
    table = _table(_sample("A", [(1000, "T")]),
                   _sample("B", [(1000, "T")], callable_bp=0))
    pair = distance.pairwise_distance(table, "A", "B", _absent_mask())
    assert "depth engine" in pair.note


def test_a_computed_distance_always_carries_its_denominator():
    table = _table(_sample("A", [(1000, "T"), (5000, "G")]),
                   _sample("B", [(1000, "T")]))
    pair = distance.pairwise_distance(table, "A", "B", _absent_mask())
    assert pair.snps == 1
    assert pair.shared_callable_sites == LENGTH
    assert pair.snps_per_mb == pytest.approx(10.0)


def test_the_denominator_is_the_intersection_not_the_reference_length():
    table = _table(_sample("A", [(1000, "T")]),
                   _sample("B", [], callable_bp=40_000))
    pair = distance.pairwise_distance(table, "A", "B", _absent_mask())
    assert pair.shared_callable_sites == 40_000


def test_a_thin_denominator_is_excluded_from_clustering_not_silently_used():
    table = _table(_sample("A", [(1000, "T")]), _sample("B", []))
    matrix = distance.distance_matrix(table, _absent_mask())
    assert matrix.comparable_pairs() == [], "100 kb is far below the floor"
    assert len(matrix.thin_pairs()) == 1
    assert matrix.min_shared_callable_sites == 3_000_000


def test_the_human_line_never_prints_a_naked_number():
    table = _table(_sample("A", [(1000, "T")]), _sample("B", []))
    matrix = distance.distance_matrix(table, _absent_mask())
    line = matrix.describe_pair("A", "B")
    assert "shared callable" in line
    assert distance.format_distance(matrix.pair("A", "B")) == \
        "1 SNPs / 100,000 bp shared callable"


def test_an_uncomputed_pair_prints_as_not_computed():
    table = _table(_sample("A", [(1000, "T")]),
                   _sample("B", [], callable_bp=0))
    matrix = distance.distance_matrix(table, _absent_mask())
    assert distance.format_distance(matrix.pair("A", "B")).startswith("not computed")
    assert len(matrix.uncomputed_pairs()) == 1


def test_the_matrix_check_names_the_thinnest_pair(tmp_path):
    table = _table(_sample("A", [(1000, "T")]), _sample("B", []))
    matrix = distance.distance_matrix(table, _absent_mask())
    named = dict((c.name, c) for c in matrix.checks)
    assert named["min_shared_callable_sites"].status == "warn"
    assert "not comparable to the published 5-SNP and 12-SNP thresholds" in \
        named["min_shared_callable_sites"].reading


def test_a_self_comparison_is_not_a_measurement():
    table = _table(_sample("A", [(1000, "T")]), _sample("B", []))
    with pytest.raises(MjolnirError) as excinfo:
        distance.pairwise_distance(table, "A", "A", _absent_mask())
    assert "self-comparison" in str(excinfo.value)


# ------------------------------------------------ an uncovered cell is not agreement

def test_a_position_outside_a_samples_callable_region_is_unknown_not_reference():
    """A coverage gap must not become evidence of identity."""
    a = _sample("A", [(1000, "T")])
    b = joint.SampleVariants(sample_id="B", variants=[],
                             callable_regions=joint.Regions("B").add(CHROM, 0, 500),
                             reference=CHROM)
    table = joint.build_joint_table([a, b], reference=CHROM)
    site = table.site_at(CHROM, 1000)
    assert table.allele("B", site) is None
    pair = distance.pairwise_distance(table, "A", "B", _absent_mask())
    assert pair.snps == 0, "the only variable position was not comparable"
    assert pair.shared_callable_sites == 500


def test_a_minority_allele_is_ambiguous_rather_than_a_genotype():
    a = joint.SampleVariants(
        sample_id="A", variants=[_variant(1000, "T", allele_fraction=0.30)],
        callable_regions=joint.Regions("A").add(CHROM, 0, LENGTH), reference=CHROM)
    b = _sample("B", [])
    table = joint.build_joint_table([a, b], reference=CHROM)
    site = table.site_at(CHROM, 1000)
    assert table.allele("A", site) == joint.AMBIGUOUS
    pair = distance.pairwise_distance(table, "A", "B", _absent_mask())
    assert pair.snps == 0
    assert pair.shared_callable_sites == LENGTH - 1, "the ambiguous site left the denominator"


# ----------------------------------------------------------------- the mask

def test_counting_without_a_mask_requires_an_explicit_decision():
    table = _table(_sample("A", [(1000, "T")]), _sample("B", []))
    with pytest.raises(MjolnirError) as excinfo:
        distance.distance_matrix(table, None)
    message = str(excinfo.value)
    assert "264,525" in message
    assert "Mask.absent" in message


def test_an_unmasked_comparison_must_state_why():
    with pytest.raises(MjolnirError):
        distance.Mask.absent(reason="")


def test_a_mask_that_claims_to_be_applied_needs_intervals():
    with pytest.raises(MjolnirError) as excinfo:
        distance.Mask(name="empty", regions=None, applied=True)
    assert "Mask.absent" in str(excinfo.value)


def test_an_unmasked_matrix_carries_the_sentence_that_says_so():
    table = _table(_sample("A", [(1000, "T")]), _sample("B", []))
    matrix = distance.distance_matrix(table, _absent_mask())
    assert distance.MASK_ABSENT_TEXT in matrix.caveats
    named = dict((c.name, c) for c in matrix.checks)
    assert named["mask_applied"].status == "fail"


def test_a_masked_position_is_removed_from_the_count(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("{0}\t900\t1100\trepeat\n".format(CHROM))
    mask = distance.load_mask(bed, name="test-mask")
    table = _table(_sample("A", [(1000, "T"), (5000, "G")]), _sample("B", []))
    pair = distance.pairwise_distance(table, "A", "B", mask)
    assert pair.snps == 1, "the masked difference at 1000 was not counted"
    assert pair.masked_sites == 200
    assert "sha256" in mask.describe()


def test_a_mask_naming_the_wrong_contig_is_fatal(tmp_path):
    """tbdb calls the H37Rv contig ``Chromosome``; most pipelines call it NC_000962.3."""
    bed = tmp_path / "mask.bed"
    bed.write_text("Chromosome\t900\t1100\n")
    mask = distance.load_mask(bed, name="tbdb-mask")
    table = _table(_sample("A", [(1000, "T")]), _sample("B", []))
    with pytest.raises(MjolnirError) as excinfo:
        distance.distance_matrix(table, mask)
    assert "would exclude nothing" in str(excinfo.value)


def test_an_empty_mask_file_is_not_the_same_as_no_mask(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("# nothing here\n")
    with pytest.raises(MjolnirError) as excinfo:
        distance.load_mask(bed)
    assert "not the same as no masking" in str(excinfo.value)


def test_the_proximity_rule_drops_clustered_differences():
    """No other variable position within 12 bases, in either genome."""
    table = _table(_sample("A", [(1000, "T"), (1005, "G"), (50_000, "C")]),
                   _sample("B", []))
    pair = distance.pairwise_distance(table, "A", "B", _absent_mask())
    assert pair.snps == 1
    assert "proximity rule" in pair.note


# ------------------------------------------------- small cohorts are not errors

def test_a_two_sample_cohort_is_an_ordinary_input():
    """MTBseq's TBjoin aborts here; the 2022 chimaera run died at exactly this point."""
    table = _table(_sample("A", [(1000, "T")]), _sample("B", []))
    named = dict((c.name, c) for c in table.checks)
    assert named["cohort_size"].status == "pass"
    assert "exactly one pairwise distance" in named["cohort_size"].reading


def test_a_single_sample_cohort_warns_rather_than_failing():
    table = _table(_sample("A", [(1000, "T")]))
    named = dict((c.name, c) for c in table.checks)
    assert named["cohort_size"].status == "warn"
    assert "absence of comparison rather than a finding" in named["cohort_size"].reading
    matrix = distance.distance_matrix(table, _absent_mask())
    assert matrix.pairs == []


def test_samples_called_against_different_references_cannot_be_joined():
    a = _sample("A", [(1000, "T")])
    b = _sample("B", [])
    b.reference = "M._chimaera_DSM44623"
    with pytest.raises(MjolnirError) as excinfo:
        joint.build_joint_table([a, b])
    assert "one coordinate system" in str(excinfo.value)
