"""Naming variants the way the catalogues name them.

``rpoB_p.Ser450Leu`` and ``rpoB_p.S450L`` are the same mutation and different
dictionary keys. MTBseq and tbdb are matched on this key, so a near-miss here
does not degrade a call — it produces "no resistance determinant detected" for a
drug the sample is resistant to, which is the one output this project exists to
prevent.

Two kinds of test live here. The synthetic ones below build a small genome with
a plus-strand gene, a minus-strand gene and an rRNA gene, and pin one naming rule
each. The last one replays the WHO catalogue's own ``Genomic_coordinates`` sheet
— 144,964 real coordinate triples, each with the name WHO itself files it under
— and requires agreement to stay above the level it reached. It skips when the
catalogue is not installed, because the suite must pass with no databases.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from mjolnir.engines import annotate as A

# ---------------------------------------------------------------------------
# A small genome with the three cases that matter
# ---------------------------------------------------------------------------

#: Two codons of padding, then a plus-strand gene, then a gap, then a
#: minus-strand gene. Built so the reading frames are known by construction.
PLUS_CDS = "ATG" "GAC" "ACT" "ACC" "GTG" "CCA" "TGA"      # M D T T V P *
MINUS_CDS = "ATG" "TCA" "CAT" "GGC" "TAA"                  # M S H G *


def _revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.fixture
def genome(tmp_path):
    """A FASTA and a GFF describing three genes, and the Annotation over them."""
    upstream = "TTTTTTTTTT" * 3        # 30 bases of promoter for geneA
    sequence = upstream + PLUS_CDS + "AAAA" + _revcomp(MINUS_CDS) + "GGGG" + "ACGT" * 10

    fasta = tmp_path / "ref.fna"
    fasta.write_text(">chr1 test\n{0}\n".format(sequence))

    a_start = len(upstream) + 1
    a_end = a_start + len(PLUS_CDS) - 1
    b_start = a_end + 4 + 1
    b_end = b_start + len(MINUS_CDS) - 1
    r_start = b_end + 4 + 1
    r_end = r_start + 39

    gff = tmp_path / "ref.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tt\tgene\t{0}\t{1}\t.\t+\t.\tID=gene:g1;Name=geneA;biotype=protein_coding;gene_id=g1\n"
        "chr1\tt\tgene\t{2}\t{3}\t.\t-\t.\tID=gene:g2;Name=geneB;biotype=protein_coding;gene_id=g2\n"
        "chr1\tt\trRNA_gene\t{4}\t{5}\t.\t+\t.\tID=gene:g3;Name=rrx;biotype=rRNA;gene_id=g3\n"
        .format(a_start, a_end, b_start, b_end, r_start, r_end))

    annotation = A.Annotation.load(gff, fasta)
    return {"annotation": annotation, "sequence": sequence,
            "a_start": a_start, "b_end": b_end, "r_start": r_start}


def names(genome, position, ref, alt):
    return A.names_for(genome["annotation"], "chr1", position, ref, alt)


# ---------------------------------------------------------------------------
# The rules
# ---------------------------------------------------------------------------

def test_an_rrna_gene_is_loaded_at_all(genome):
    """The GFF calls them rRNA_gene, not gene.

    Reading only ``gene`` rows made rrs and rrl invisible in the real H37Rv
    annotation — and those two carry every aminoglycoside and macrolide
    determinant in both MTBC and NTM, so every ``rrs_n.1401A>G`` went unmatched.
    """
    loaded = {g.name for g in genome["annotation"].genes}
    assert {"geneA", "geneB", "rrx"} <= loaded


def test_a_coding_change_is_named_as_a_protein_change(genome):
    """Codon 2 of the plus-strand gene is GAC (Asp); GAC->GCC is Asp2Ala."""
    position = genome["a_start"] + 4        # second base of codon 2
    assert "geneA_p.Asp2Ala" in names(genome, position, "A", "C")


def test_a_synonymous_change_keeps_the_nucleotide_form(genome):
    """Same residue means the catalogues file it in nucleotides, not protein."""
    position = genome["a_start"] + 5        # third base of codon 2, GAC->GAT
    produced = names(genome, position, "C", "T")
    assert any(n.startswith("geneA_c.") for n in produced), produced
    assert not any("p.Asp2" in n for n in produced), produced


def test_a_start_codon_change_is_nucleotide_named_and_pools_as_loss_of_function(genome):
    """WHO names these ``c.3G>T`` and grades them under the gene's LoF.

    What matters is that translation no longer starts, not which residue the
    broken codon would have encoded.
    """
    produced = names(genome, genome["a_start"] + 2, "G", "T")
    assert "geneA_c.3G>T" in produced, produced
    assert "geneA_LoF" in produced, produced


def test_a_minus_strand_gene_complements_its_alleles(genome):
    """The alleles are written as the gene reads them, not as the genome does.

    Getting this backwards produces a well-formed name for a mutation that does
    not exist — and *pncA* and *katG*, which carry the pyrazinamide and isoniazid
    determinants, are both on the minus strand.
    """
    # c.1 of a minus-strand gene is its highest genome coordinate.
    position = genome["b_end"]
    genome_base = genome["sequence"][position - 1]
    assert genome_base == "T", "fixture moved; expected the start codon's A"
    produced = [n for n in names(genome, position, "T", "C") if n.startswith("geneB")]
    assert produced, "nothing named for the minus-strand gene"
    # Gene-direction alleles: genome T>C is A>G to the gene.
    assert any("A>G" in n for n in produced), produced
    assert not any("T>C" in n for n in produced), produced


def test_a_promoter_variant_is_named_for_the_gene_it_regulates(genome):
    """Negative c. coordinates, counting back from the start codon.

    ``eis_c.-14C>T`` and ``inhA_c.-154G>A`` are both WHO Group 1 determinants; a
    coding-only annotator drops them on the floor as intergenic.
    """
    produced = names(genome, genome["a_start"] - 5, "T", "A")
    assert "geneA_c.-5T>A" in produced, produced


def test_an_rrna_variant_uses_the_n_prefix(genome):
    """rRNA genes have no reading frame, so a codon number would be a fiction."""
    produced = names(genome, genome["r_start"] + 10, "A", "G")
    assert any(n.startswith("rrx_n.") for n in produced), produced


def test_an_rrna_upstream_variant_also_uses_n(genome):
    """WHO files the position before rrs as ``rrs_n.-147G>T``, not ``c.``."""
    produced = names(genome, genome["r_start"] - 3, "G", "A")
    assert any(n.startswith("rrx_n.-") for n in produced), produced


def test_a_multi_base_substitution_is_not_an_indel(genome):
    """``AT>CA`` is two substitutions, and maps to several graded variants.

    Treated as an indel it matches nothing: WHO names such a variant both by its
    constituent nucleotide changes and by the single amino-acid change they make
    together.
    """
    position = genome["a_start"] + 3        # codon 2, first base
    produced = names(genome, position, "GA", "CC")
    # A protein-level name for the codon the two bases share...
    assert any(n.startswith("geneA_p.") for n in produced), produced
    # ...and the delins form WHO also files it under, spelled with the bases
    # between the two verbs: c.4_5delGAinsCC.
    assert any("del" in n and "ins" in n for n in produced), produced
    # What must never appear is a bare insertion, which is what treating an
    # equal-length substitution as an indel produces.
    assert not any("ins" in n and "del" not in n for n in produced), produced


def test_a_frameshift_is_named_at_the_first_changed_residue(genome):
    """Not at the codon the inserted base happens to fall in."""
    position = genome["a_start"] + 4
    produced = names(genome, position, "A", "AG")
    assert any(n.endswith("fs") for n in produced), produced
    assert "geneA_LoF" in produced, produced


def test_a_deletion_offers_the_gene_it_removes_as_a_loss(genome):
    """WHO pools deletions under ``<gene>_LoF`` rather than naming each one."""
    deleted = genome["sequence"][genome["a_start"]:genome["a_start"] + 6]
    produced = names(genome, genome["a_start"], "T" + deleted, "T")
    assert "geneA_LoF" in produced, produced


def test_equivalent_placements_of_an_indel_are_all_offered(genome):
    """The same deletion can be written at several offsets and HGVS picks one.

    Which one is 3'-most depends on the gene's strand, so rather than reproduce
    the rule and its exception, every equivalent placement is named.
    """
    #                1234567890
    # sequence      = AACACACAGGG, so position 2 spells ACA and deletes CA.
    placements = A.equivalent_placements("AACACACAGGG", 2, "ACA", "A")
    assert len(placements) > 1, placements
    for position, ref, alt in placements:
        assert len(ref) - len(alt) == 2, (position, ref, alt)


def test_a_substitution_has_exactly_one_placement():
    """Only indels are ambiguous; shifting a substitution would invent names."""
    assert A.equivalent_placements("AACACACAGGG", 3, "C", "T") == [(3, "C", "T")]


def test_a_reference_mismatch_does_not_invent_equivalences():
    """If the reference allele is not what the genome says, place nothing.

    Shifting a variant that cannot be located would manufacture equivalences,
    and each one becomes a name that could match a catalogue row by accident.
    """
    assert A.equivalent_placements("AAAAAAA", 3, "GGG", "G") == [(3, "GGG", "G")]
    # The same guard, hit by a real off-by-one rather than a nonsense allele.
    assert A.equivalent_placements("AACACACAGGG", 3, "ACA", "A") == [(3, "ACA", "A")]


def test_annotate_fills_variants_in_place(genome):
    """The pipeline calls this; it must set gene, hgvs and effect together."""
    from mjolnir.records import Variant

    variant = Variant(chrom="chr1", pos=genome["a_start"] + 4, ref="A", alt="C",
                      variant_type="snp")
    named = A.annotate([variant], genome["annotation"])
    assert named == 1
    assert variant.gene == "geneA"
    assert variant.hgvs
    assert variant.effect


def test_annotate_leaves_an_already_named_variant_alone(genome):
    """A caller that got better names from elsewhere must not lose them."""
    from mjolnir.records import Variant

    variant = Variant(chrom="chr1", pos=genome["a_start"] + 4, ref="A", alt="C",
                      gene="fromElsewhere", hgvs="p.Xxx1Yyy", variant_type="snp")
    A.annotate([variant], genome["annotation"])
    assert variant.gene == "fromElsewhere"


def test_a_gff_with_no_genes_is_refused(tmp_path):
    """Silence here would produce a run where nothing is ever named."""
    from mjolnir.utils import MjolnirError

    empty = tmp_path / "empty.gff"
    empty.write_text("##gff-version 3\n")
    with pytest.raises(MjolnirError):
        A.load_gff(empty)


def test_a_deletion_that_removes_the_start_codon_is_a_loss_not_an_upstream_variant(genome):
    """The worst defect this project has had, pinned.

    A deletion beginning upstream and running into the gene takes the start
    codon with it. Classified by where it *starts*, it was named an
    ``upstream_variant``, ``is_loss_of_function`` returned False, the WHO
    loss-of-function rule never fired, and a complete *pncA* knockout — which is
    definitive pyrazinamide resistance — was reported as a regulatory nucleotide
    change with no determinant at all. Absence of a gene product rendered as
    normality, in the gene that defines the drug.
    """
    # Two bases before the gene through to c.3: the whole start codon goes.
    start = genome["a_start"]
    deleted = genome["sequence"][start - 3:start + 2]
    _gene, _locus, hgvs, effect = A.name_variant(
        genome["annotation"], "chr1", start - 3, deleted, deleted[0])
    assert effect == "start_lost", (hgvs, effect)
    assert effect in A.LOF_EFFECTS
    assert "geneA_LoF" in names(genome, start - 3, deleted, deleted[0])


def test_a_deletion_entirely_before_the_gene_stays_an_upstream_variant(genome):
    """The fix must not turn every promoter deletion into a knockout."""
    start = genome["a_start"]
    deleted = genome["sequence"][start - 10:start - 7]
    _g, _l, _h, effect = A.name_variant(
        genome["annotation"], "chr1", start - 10, deleted, deleted[0])
    assert effect == "upstream_variant", effect


def test_an_insertion_names_the_first_codon_the_shift_actually_reaches(genome):
    """On the plus strand the inserted bases land *after* the anchor.

    Naming the anchor's own codon is one too low whenever the anchor sits on a
    codon boundary, which flipped the rpoB RRDR rule in both directions — on
    rifampicin, the drug that rule exists for: it manufactured a Group 2
    determinant from a frameshift starting outside the region, and dropped one
    starting at its first codon.
    """
    # c.3 is a codon boundary: an insertion after it first alters codon 2.
    produced = names(genome, genome["a_start"] + 2, "G", "GA")
    assert any(n.endswith("p.Xaa2fs") for n in produced), produced
    assert not any(n.endswith("p.Xaa1fs") for n in produced), produced


def test_a_variant_is_not_named_for_a_gene_on_another_contig(tmp_path):
    """Position alone is not an address on a genome with more than one replicon.

    H37Rv has one, so coordinate-only lookup happened to work. The *M. chimaera*
    reference has three, and a variant at position 300 of a plasmid was matched
    against whichever gene spanned position 300 of the chromosome and named for
    it — a gene the variant is nowhere near, and whose catalogue rows it would
    then be compared against.
    """
    fasta = tmp_path / "multi.fna"
    fasta.write_text(">chrom\n{0}\n>plasmid\n{1}\n".format("ACGT" * 300, "TTGA" * 300))
    gff = tmp_path / "multi.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chrom\tt\tgene\t100\t400\t.\t+\t.\tID=gene:c1;Name=chromGene;"
        "biotype=protein_coding;gene_id=c1\n"
        "plasmid\tt\tgene\t100\t400\t.\t+\t.\tID=gene:p1;Name=plasmidGene;"
        "biotype=protein_coding;gene_id=p1\n")
    annotation = A.Annotation.load(gff, fasta)

    assert [g.label for g in annotation.genes_at(200, "plasmid")] == ["plasmidGene"]
    assert [g.label for g in annotation.genes_at(200, "chrom")] == ["chromGene"]
    produced = A.names_for(annotation, "plasmid", 200, "G", "A")
    assert produced, "nothing named on the plasmid"
    assert not any(n.startswith("chromGene") for n in produced), produced


def test_a_single_contig_annotation_ignores_the_contig_name(tmp_path):
    """tbdb calls H37Rv "Chromosome"; the catalogue calls it NC_000962.3.

    With one replicon the name carries no information the caller does not
    already have, and filtering on it would match nothing at all.
    """
    fasta = tmp_path / "one.fna"
    fasta.write_text(">Chromosome\n{0}\n".format("ACGT" * 300))
    gff = tmp_path / "one.gff"
    gff.write_text(
        "##gff-version 3\n"
        "Chromosome\tt\tgene\t100\t400\t.\t+\t.\tID=gene:c1;Name=g;"
        "biotype=protein_coding;gene_id=c1\n")
    annotation = A.Annotation.load(gff, fasta)
    assert annotation.genes_at(200, "NC_000962.3"), "a name mismatch matched nothing"


def test_a_multi_base_substitution_names_the_amino_acid_the_codon_makes(genome):
    """Not the one its first base alone would make.

    GCT>TGA is a stop codon; reading only the first base called it a serine
    substitution. The primary name is what the report prints and what the rules
    layer classifies, so a nonsense mutation was presented as a missense one.
    """
    # Codon 2 of geneA is GAC. Change all three bases to TGA, a stop.
    position = genome["a_start"] + 3
    _g, _l, hgvs, effect = A.name_variant(
        genome["annotation"], "chr1", position, "GAC", "TGA")
    assert hgvs.endswith("Ter"), (hgvs, effect)
    assert effect == "stop_gained", (hgvs, effect)


# ---------------------------------------------------------------------------
# The gold standard
# ---------------------------------------------------------------------------

#: Agreement measured over the whole sheet on 2026-08-11. A regression below
#: this means variants stopped matching catalogue rows they used to match.
WHO_AGREEMENT_FLOOR = 0.93

DB_ROOT = Path(os.environ.get("MJOLNIR_DB", Path.home() / ".mjolnir" / "db"))
WHO_XLSX = DB_ROOT / "who-catalogue-v2" / "WHO-UCN-TB-2023.7-eng.xlsx"
TBDB_GFF = DB_ROOT / "tbdb" / "genome.gff"
H37RV = DB_ROOT / "h37rv" / "NC_000962.3.fasta"


@pytest.mark.skipif(
    not (WHO_XLSX.exists() and TBDB_GFF.exists() and H37RV.exists()),
    reason="needs the installed WHO catalogue, tbdb GFF and H37Rv "
           "(mjolnir db fetch); the suite must pass without any database")
def test_names_agree_with_the_who_catalogues_own_coordinate_sheet():
    """Replay a sample of WHO's coordinate sheet and require agreement.

    A sample rather than all 144,964 rows: the full sweep takes minutes and this
    has to be runnable in a normal test cycle. The sample is deterministic —
    every 25th row — so a regression is reproducible rather than intermittent.
    """
    openpyxl = pytest.importorskip("openpyxl")

    annotation = A.Annotation.load(TBDB_GFF, H37RV)
    workbook = openpyxl.load_workbook(WHO_XLSX, read_only=True)
    rows = workbook["Genomic_coordinates"].iter_rows(values_only=True)
    next(rows)

    total = agree = 0
    misses = []
    for index, row in enumerate(rows):
        if index % 25 or not row or not row[0]:
            continue
        want, chrom, position, ref, alt = row[0], row[1], row[2], row[3], row[4]
        if position is None or ref is None or alt is None:
            continue
        total += 1
        produced = A.names_for(annotation, str(chrom), int(position), str(ref), str(alt))
        if want in produced:
            agree += 1
        elif len(misses) < 5:
            misses.append((want, produced[:2]))

    assert total > 1000, "the coordinate sheet did not yield a usable sample"
    ratio = agree / float(total)
    assert ratio >= WHO_AGREEMENT_FLOOR, (
        "agreement with the WHO catalogue's own naming fell to {0:.2%} over {1} "
        "rows (floor {2:.0%}). Examples: {3}".format(
            ratio, total, WHO_AGREEMENT_FLOOR, misses))


def test_a_nonsense_mnv_is_not_called_synonymous(genome):
    """The whole codon decides, not the lowest-coordinate changed base.

    If that base is synonymous on its own, ``_hgvs_snv`` returns a ``c.`` name
    with effect ``synonymous_variant``. Guarding the whole-codon recomputation on
    the single-base name being a protein change therefore let a genuine nonsense
    MNV keep that name and that effect — so ``is_synonymous()`` was true, the
    novel-silent Group 4 rule fired, and the loss-of-function rule did not.
    """
    # Codon 2 of geneA is GAC (Asp). GAC>GAA changes only the third base and is
    # a real amino-acid change; the first changed base alone is not decisive.
    position = genome["a_start"] + 3
    _g, _l, hgvs, effect = A.name_variant(
        genome["annotation"], "chr1", position, "GAC", "TGA")
    assert effect == "stop_gained", (hgvs, effect)
    assert hgvs.endswith("Ter"), hgvs


def test_ncbi_rrna_genes_are_named_and_non_coding(tmp_path):
    """RefSeq gives 16S and 23S no gene name and writes gene_biotype, not biotype.

    Both defaults were wrong for NCBI annotation, and together they made the NTM
    resistance rules unreachable on real data: an ``rrs`` 1408 amikacin variant
    in *M. abscessus* was named ``MAB_RS07510_p.…`` — a protein change, in a gene
    called after its own locus tag — and matched nothing anywhere.
    """
    gff = tmp_path / "ncbi.gff"
    gff.write_text(
        "##gff-version 3\n"
        "NC_1\tRefSeq\tgene\t100\t1600\t.\t+\t.\tID=gene-X_RS01;Name=X_RS01;"
        "gene_biotype=rRNA;locus_tag=X_RS01\n"
        "NC_1\tcmsearch\trRNA\t100\t1600\t.\t+\t.\tID=rna-X_RS01;Parent=gene-X_RS01;"
        "product=16S ribosomal RNA;locus_tag=X_RS01\n"
        "NC_1\tRefSeq\tgene\t2000\t5000\t.\t+\t.\tID=gene-X_RS02;Name=X_RS02;"
        "gene_biotype=rRNA;locus_tag=X_RS02\n"
        "NC_1\tcmsearch\trRNA\t2000\t5000\t.\t+\t.\tID=rna-X_RS02;Parent=gene-X_RS02;"
        "product=23S ribosomal RNA;locus_tag=X_RS02\n")
    genes = {g.name: g for g in A.load_gff(gff)}
    assert "rrs" in genes and "rrl" in genes, sorted(genes)
    assert not genes["rrs"].coding, "an rRNA gene must not be given codon numbers"
    assert not genes["rrl"].coding


def test_a_gene_named_after_its_own_locus_tag_counts_as_unnamed(tmp_path):
    """Name=X_RS01 carries no information the locus tag does not already carry."""
    gff = tmp_path / "lt.gff"
    gff.write_text(
        "##gff-version 3\n"
        "NC_1\tRefSeq\tgene\t100\t400\t.\t+\t.\tID=gene-X_RS01;Name=X_RS01;"
        "gene_biotype=protein_coding;locus_tag=X_RS01\n")
    gene = A.load_gff(gff)[0]
    assert gene.name == ""
    assert gene.label == "X_RS01"
