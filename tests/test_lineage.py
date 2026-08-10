"""Lineage from the SNP barcode, with the support printed beside the label (§6).

Inside the MTBC, ANI resolves nothing — the members sit at 99.21-99.92% of one
another — so the member is called from lineage-defining SNPs, and the number of
sites behind a call is part of the call. ``La1.2.BCG``, ``La2`` and ``La3`` are
supported by five sites each in tbdb's scheme, the thinnest evidence in the file,
and BCG carries intrinsic pyrazinamide resistance: a call made or lost by a
coverage gap over five positions therefore changes treatment.

The barcode is genotyped from a pileup at the published per-platform read minima,
so the same ONT-versus-Illumina asymmetry that governs variant calling governs
which barcode sites are callable at all.
"""

from __future__ import annotations

import pytest

from mjolnir import config
from mjolnir.typing import lineage
from mjolnir.utils import MjolnirError

CHROM = "Chromosome"

#: A miniature scheme: a lineage, one of its sublineages, and BCG on the five
#: defining sites tbdb actually gives it.
SCHEME = (
    [(CHROM, 1000 + i, "lineage4", "T", "Euro-American") for i in range(10)]
    + [(CHROM, 2000 + i, "lineage4.9", "G", "Euro-American (H37Rv-like)")
       for i in range(8)]
    + [(CHROM, 3000 + i, "La1.2.BCG", "C", "M. bovis BCG") for i in range(5)]
)


@pytest.fixture
def barcode(tmp_path):
    path = tmp_path / "barcode.bed"
    path.write_text("".join(
        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tNone\tNone\n".format(
            chrom, pos - 1, pos, taxon, allele, family)
        for chrom, pos, taxon, allele, family in SCHEME))
    return lineage.load_barcode(path)


def pileup_for(sites, supported, depth=30):
    """Counts putting the derived allele at every site of the named taxa."""
    counts = {}
    for site in sites:
        if site.taxon in supported:
            base = site.allele
        else:
            base = "A" if site.allele != "A" else "C"
        counts[(site.chrom, site.pos)] = {base: depth}
    return counts


# ------------------------------------------------------------ loading the scheme

def test_the_scheme_loads_with_its_taxa_and_families(barcode):
    assert len(barcode) == 23
    assert set(s.taxon for s in barcode) == {"lineage4", "lineage4.9", "La1.2.BCG"}
    assert "23 defining sites" in lineage.scheme_description(barcode)


def test_a_barcode_line_with_a_nonsense_allele_is_fatal(tmp_path):
    path = tmp_path / "barcode.bed"
    path.write_text("{0}\t999\t1000\tlineage4\tZ\n".format(CHROM))
    with pytest.raises(MjolnirError) as excinfo:
        lineage.load_barcode(path)
    assert "IUPAC" in str(excinfo.value)


def test_a_missing_barcode_says_how_to_fetch_it(tmp_path):
    with pytest.raises(MjolnirError) as excinfo:
        lineage.load_barcode(tmp_path / "barcode.bed")
    assert "fetch" in str(excinfo.value)


# ------------------------------------------------------------- the call and its support

def test_a_lineage_call_carries_the_sites_behind_it(barcode):
    call = lineage.call_lineage(barcode, pileup_for(barcode, {"lineage4", "lineage4.9"}),
                                "illumina")
    assert call.lineage == "lineage4"
    assert call.sublineage == "lineage4.9"
    assert call.barcode_sites_supporting == 18
    assert call.barcode_sites_callable == 18
    assert call.support_fraction == 1.0
    assert call.scheme, "the scheme and its version travel with the call"


def test_nothing_supported_is_not_determined_rather_than_lineage_four(barcode):
    """The commonest lineage must never be the fallback for an unresolved call."""
    call = lineage.call_lineage(barcode, pileup_for(barcode, set()), "illumina")
    assert call.lineage == ""
    assert call.display == "not determined"
    assert call.confidence == "none"
    assert any("not the same as lineage 4" in caveat for caveat in call.caveats)


def test_the_thresholds_the_call_had_to_clear_are_stated(barcode):
    call = lineage.call_lineage(barcode, pileup_for(barcode, set()), "illumina")
    assert any("{0:.0%}".format(config.MIN_BARCODE_CALLABLE_FRACTION) in caveat
               for caveat in call.caveats)


# ------------------------------------------------------------------------ BCG

def test_bcg_is_called_and_flagged(barcode):
    call = lineage.call_lineage(barcode, pileup_for(barcode, {"La1.2.BCG"}), "illumina")
    assert call.is_bcg is True
    assert call.is_animal is True
    assert call.animal_variant == "BCG"
    assert call.sublineage == "La1.2.BCG"


def test_bcg_carries_the_intrinsic_pyrazinamide_note(barcode):
    """A property of the organism, not a catalogue lookup."""
    call = lineage.call_lineage(barcode, pileup_for(barcode, {"La1.2.BCG"}), "illumina")
    assert config.BCG_PZA_NOTE in call.caveats
    named = dict((c.name, c) for c in lineage.lineage_checks(call))
    assert "bcg_intrinsic_pyrazinamide_resistance" in named


def test_a_five_site_call_is_never_reported_as_high_confidence(barcode):
    """At five defining sites a coverage gap or a contaminant moves the call."""
    call = lineage.call_lineage(barcode, pileup_for(barcode, {"La1.2.BCG"}), "illumina")
    assert call.confidence != "high"
    assert any("5 site" in caveat for caveat in call.caveats)


def test_the_animal_lineage_uses_the_la_nomenclature(barcode):
    call = lineage.call_lineage(barcode, pileup_for(barcode, {"La1.2.BCG"}), "illumina")
    assert call.lineage == "La1"


def test_a_scheme_missing_animal_taxa_says_which_it_cannot_call(barcode):
    call = lineage.call_lineage(barcode, pileup_for(barcode, {"lineage4"}), "illumina")
    assert any("no defining sites for" in caveat for caveat in call.caveats)


# ------------------------------------- the platform decides which sites are callable

def test_four_reads_genotype_a_site_on_illumina_but_not_on_ont(barcode):
    """3 reads on Illumina, 5 on ONT: the published clinical-DST minima."""
    counts = pileup_for(barcode, {"lineage4", "lineage4.9"}, depth=4)
    illumina = lineage.call_lineage(barcode, counts, "illumina")
    ont = lineage.call_lineage(barcode, counts, "ont")
    assert illumina.lineage == "lineage4"
    assert illumina.barcode_sites_callable == 18
    assert ont.lineage == ""
    assert ont.barcode_sites_callable == 0


def test_five_reads_genotype_a_site_on_both(barcode):
    counts = pileup_for(barcode, {"lineage4", "lineage4.9"}, depth=5)
    assert lineage.call_lineage(barcode, counts, "ont").lineage == "lineage4"


def test_a_site_with_no_reads_is_uncallable_rather_than_disagreeing(barcode):
    """Absence of coverage is not evidence against the lineage."""
    counts = pileup_for(barcode, {"lineage4", "lineage4.9"})
    for site in barcode[:4]:
        counts[(site.chrom, site.pos)] = {}
    genotypes = lineage.genotype_barcode_sites(barcode, counts, "illumina")
    blank = [g for g in genotypes if g.site.pos < 1004]
    assert len(blank) == 4
    assert all(g.is_callable is False for g in blank)
    assert all(g.supports is False for g in blank)
    assert sum(1 for g in genotypes if g.is_callable) == 19


def test_a_partly_covered_taxon_drops_out_of_the_support_count(barcode):
    """Six of ten sites is below the 80% callable floor, so lineage4 stops voting."""
    counts = pileup_for(barcode, {"lineage4", "lineage4.9"})
    for site in barcode[:4]:
        counts[(site.chrom, site.pos)] = {}
    call = lineage.call_lineage(barcode, counts, "illumina")
    assert call.barcode_sites_callable == 8, "only the sublineage's sites still count"
    assert call.barcode_sites_supporting == 8
    assert call.support_fraction == 1.0
    assert call.sublineage == "lineage4.9"


# -------------------------------------------------------------- mixed lineages

def test_two_supported_lineages_are_reported_as_mixed_not_resolved(barcode):
    counts = pileup_for(barcode, {"lineage4", "lineage4.9", "La1.2.BCG"})
    call = lineage.call_lineage(barcode, counts, "illumina")
    assert call.mixed_lineages, "both taxa were supported and both must be named"
    assert call.confidence in ("low", "none", "moderate")


def test_lineage_is_not_applicable_to_an_ntm(tmp_path):
    call = lineage.lineage_not_applicable("Mycobacterium chimaera")
    assert call.lineage == ""
    assert call.display == "not determined"
    assert call.caveats
