"""The engines' pure layer: pileup counting, VCF parsing and the §7 thresholds.

Nothing here maps a read or calls a variant. What is tested is the arithmetic
that turns a pileup line and a VCF record into the evidence every later stage
reasons over, and the platform thresholds applied to it — 3 supporting reads on
Illumina, 5 on ONT, 90% for a major variant.

Two of these tests are about reproducing a legacy convention rather than a
correct one. MTBseq's frequency denominator includes N and GAP counts and its
tie-break puts GAP above every base, so ``--compat mtbseq`` has to do the same or
the two tools cannot be compared on the same reads. Mjolnir's own convention
counts ACGT only and reports an exact tie as no genotype.
"""

from __future__ import annotations

import pytest

from mjolnir import config
from mjolnir.engines import call as call_engine
from mjolnir.engines import pileup
from mjolnir.records import PLATFORM_FASTA, Variant
from mjolnir.utils import MjolnirError


def _site(counts, **kwargs):
    return pileup.PileupSite(chrom="NC_000962.3", pos=761155, ref_base="C",
                             counts=dict(counts), **kwargs)


# ------------------------------------------------------- an uncovered site is None

def test_an_uncovered_site_has_no_fractions_at_all():
    """Zero would say the allele was looked for and not found."""
    site = pileup.uncovered_site("NC_000962.3", 761155, "C")
    assert site.covered is False
    assert site.fraction("T") is None
    assert site.unambiguous_fraction() is None
    assert site.is_unambiguous() is None


def test_a_covered_site_with_no_alt_reads_reports_zero_not_none():
    site = _site({"C": 40})
    assert site.fraction("T") == 0.0
    assert site.fraction("C") == 1.0


# ---------------------------------------------------------- the two conventions

def test_the_mtbseq_denominator_includes_n_and_gap():
    """Design §9b: a conventional ACGT depth gives a different fraction."""
    site = _site({"C": 60, "T": 20, "N": 10, "*": 10})
    assert site.denominator(pileup.CONVENTION_ACGT) == 80
    assert site.denominator(pileup.CONVENTION_MTBSEQ) == 100
    assert site.fraction("T", pileup.CONVENTION_ACGT) == pytest.approx(0.25)
    assert site.fraction("T", pileup.CONVENTION_MTBSEQ) == pytest.approx(0.20)


def test_gap_wins_every_tie_under_the_mtbseq_convention():
    site = _site({"A": 10, "*": 10})
    assert site.major_allele(pileup.CONVENTION_MTBSEQ) == pileup.GAP


def test_an_exact_tie_is_no_genotype_under_mjolnirs_own_convention():
    """Two alleles at 50% is a real observation; picking one invents a call."""
    site = _site({"A": 10, "T": 10})
    assert site.major_allele() is None
    assert site.unambiguous_fraction() is None


def test_an_unknown_convention_raises():
    with pytest.raises(MjolnirError):
        _site({"A": 10}).denominator("whatever")


def test_the_unambiguous_fraction_is_surfaced_rather_than_applied():
    """MTBseq discards a position below 95% and says nothing about it."""
    site = _site({"C": 70, "T": 30})
    assert site.unambiguous_fraction() == pytest.approx(0.70)
    assert site.is_unambiguous(threshold_percent=config.MTBSEQ_UNAMBIG) is False


# ---------------------------------------------------- genotyping at the platform floor

@pytest.mark.parametrize("platform,reads,supported", [
    ("illumina", 2, False),
    ("illumina", 3, True),
    ("ont", 4, False),
    ("ont", 5, True),
])
def test_the_read_floor_is_the_platforms_own(platform, reads, supported):
    site = _site({"T": reads, "C": 0})
    genotype = pileup.genotype_site(site, platform=platform,
                                    min_reads=config.min_reads_for(platform))
    assert genotype.supported is supported
    if not supported:
        assert "below the {0} threshold".format(platform) in genotype.reason
        assert "Colpus" in genotype.reason, "the reason names the paper"


def test_a_site_below_the_floor_is_unsupported_not_a_reference_call():
    site = _site({"T": 2})
    genotype = pileup.genotype_site(site, platform="illumina", min_reads=3)
    assert genotype.supported is False
    assert genotype.allele == "T", "what was seen is still reported"


def test_ont_takes_the_highest_depth_allele_with_no_fraction_cut():
    site = _site({"T": 12, "C": 8})
    genotype = pileup.genotype_site(site, platform="ont", min_reads=5,
                                    min_fraction=0.9)
    assert genotype.supported is True
    assert genotype.allele == "T"
    assert "highest-depth allele" in genotype.method


def test_illumina_applies_the_fraction_as_well_as_the_reads():
    site = _site({"T": 12, "C": 8})
    genotype = pileup.genotype_site(site, platform="illumina", min_reads=3,
                                    min_fraction=0.9)
    assert genotype.supported is False
    assert "below the" in genotype.reason


def test_both_alleles_of_a_mixed_site_are_visible():
    """A barcode site at 60/40 is what a mixed infection is made of."""
    site = _site({"T": 60, "C": 40})
    present = pileup.alleles_present(site)
    assert [base for base, _f, _c in present] == ["T", "C"]


# ------------------------------------------------------------- pileup parsing

def test_a_pileup_line_is_parsed_into_bases_indels_and_depth():
    line = "NC_000962.3\t761155\tC\t9\t.,.,.,TTT\tIIIIIIIII"
    site = pileup.parse_pileup_line(line)
    assert site is not None
    assert site.chrom == "NC_000962.3" and site.pos == 761155
    assert site.count("C") == 6, "'.' and ',' are reference matches on either strand"
    assert site.count("T") == 3
    assert site.raw_depth == 9
    assert site.fraction("T") == pytest.approx(1 / 3.0)


def test_read_start_and_read_end_markers_are_not_counted_as_bases():
    """``^]`` carries a mapping quality after it and ``$`` follows its base."""
    site = pileup.parse_pileup_line("NC_000962.3\t761155\tC\t3\t^].,.$\tIII")
    assert site.count("C") == 3
    assert site.count("N") == 0


def test_an_insertion_and_a_deletion_are_counted_separately():
    line = "NC_000962.3\t761155\tC\t4\t.+2AC.-1G..\tIIII"
    site = pileup.parse_pileup_line(line)
    assert site.insertion_count("AC") == 1
    assert site.deletion_count("G") == 1


def test_a_deleted_position_is_a_gap_not_a_missing_read():
    site = pileup.parse_pileup_line("NC_000962.3\t761155\tC\t4\t..**\tIIII")
    assert site.count(pileup.GAP) == 2
    assert site.denominator(pileup.CONVENTION_ACGT) == 2
    assert site.denominator(pileup.CONVENTION_MTBSEQ) == 4


# ---------------------------------------------------------------- VCF parsing

VCF = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
    "NC_000962.3\t761155\t.\tC\tT\t228\tPASS\tDP=64\tGT:AD:DP\t1:2,62:64",
    "NC_000962.3\t2289252\t.\tGC\tG\t180\tPASS\tDP=40\tGT:AD:DP\t1:38,2:40",
])


def test_a_vcf_record_becomes_a_variant_with_its_support():
    variants = call_engine.parse_vcf_lines(VCF.splitlines(), source_caller="bcftools")
    assert len(variants) == 2
    first = variants[0]
    assert first.coordinate_key == ("NC_000962.3", 761155, "C", "T")
    assert first.alt_reads == 62
    assert first.depth == 64
    assert first.allele_fraction == pytest.approx(62 / 64.0)
    assert first.source_caller == "bcftools"


def test_an_indel_is_typed_as_one():
    variants = call_engine.parse_vcf_lines(VCF.splitlines())
    assert variants[1].is_indel is True
    assert call_engine.classify_variant_type("GC", "G") == "deletion"
    assert call_engine.classify_variant_type("G", "GC") == "insertion"
    assert call_engine.classify_variant_type("C", "T") == "snp"


# ------------------------------------------------- the platform filters (design §7)

def _variant(alt_reads, depth, fraction):
    return Variant(chrom="NC_000962.3", pos=761155, ref="C", alt="T",
                   alt_reads=alt_reads, depth=depth, allele_fraction=fraction)


def test_a_variant_below_the_illumina_read_floor_is_marked_not_deleted():
    variants = call_engine.apply_platform_filters([_variant(2, 40, 0.05)], "illumina")
    assert len(variants) == 1, "the annex must be able to show what was rejected"
    assert any("3" in f for f in variants[0].filters)


def test_the_same_two_reads_are_below_the_ont_floor_too():
    variants = call_engine.apply_platform_filters([_variant(4, 40, 0.10)], "ont")
    assert any("5" in f for f in variants[0].filters)


def test_four_reads_pass_on_illumina_and_fail_on_ont():
    illumina = call_engine.apply_platform_filters([_variant(4, 40, 0.10)], "illumina")
    ont = call_engine.apply_platform_filters([_variant(4, 40, 0.10)], "ont")
    assert not any("min_reads" in f or "reads" in f for f in illumina[0].filters)
    assert ont[0].filters


def test_the_major_variant_threshold_is_ninety_percent():
    variants = call_engine.apply_platform_filters(
        [_variant(62, 64, 0.97), _variant(20, 64, 0.31)], "illumina")
    assert variants[0].is_major is True
    assert variants[1].is_major is False


def test_an_ont_minor_variant_is_caveated_rather_than_filtered():
    """ONT under-detects minor variants; a fraction there cannot be trusted either way."""
    variants = call_engine.apply_platform_filters([_variant(9, 60, 0.15)], "ont")
    assert variants[0].note == config.ONT_MINOR_VARIANT_CAVEAT
    assert not any("minor" in f for f in variants[0].filters)


def test_an_assembly_gets_no_read_thresholds_and_says_why():
    variant = Variant(chrom="NC_000962.3", pos=761155, ref="C", alt="T")
    variants = call_engine.apply_platform_filters([variant], PLATFORM_FASTA)
    assert variants[0].filters == []
    assert "no read support" in variants[0].note


def test_each_applied_threshold_names_its_source():
    sources = call_engine.threshold_sources("ont")
    assert "Colpus" in sources["min_reads"]
    assert sources["major_variant_fraction"]
    assert sources["min_mapping_quality"]
