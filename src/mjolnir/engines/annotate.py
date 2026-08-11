"""Gene names and HGVS names for called variants.

Without this module a variant is a coordinate and nothing else, and two of the
three catalogues are unreachable: WHO is matched on ``(chrom, pos, ref, alt)``
and works, but MTBseq and tbdb are keyed on ``<gene>_<hgvs>`` and match nothing
at all. Measured on a real *M. bovis* isolate before this existed: 3,013
variants, 0 gene names, 110 WHO rows matched, 0 from either other catalogue. The
"consensus across three catalogues" was WHO alone, and the NTM ``erm(41)`` /
``rrl`` / ``rrs`` rules — all gene-keyed — could never fire.

**The naming has to match the catalogues exactly or it is worse than useless.**
``rpoB_p.Ser450Leu`` and ``rpoB_p.S450L`` are the same mutation and different
dictionary keys, and a near-miss silently produces "no determinant detected" for
a drug the sample is resistant to. So this module is written against a gold
standard rather than against a specification: the WHO catalogue ships a
``Genomic_coordinates`` sheet mapping tens of thousands of ``(position,
reference, alternative)`` triples to the graded-variant name WHO itself uses,
and :mod:`tests.test_annotate` replays that sheet through this code and requires
agreement. Every rule below exists because that comparison demanded it.

What the names look like, and why each form exists:

``rpoB_p.Ser450Leu``   a coding change, three-letter amino acids
``rpoB_c.1349C>T``     a synonymous coding change, stated in nucleotides
``eis_c.-14C>T``       a promoter change, numbered backwards from the start codon
``rrs_n.1401A>G``      a non-coding RNA gene, which has no reading frame
``pncA_p.Thr47fs``     a frameshift
``Rv0678_c.193_194del`` a deletion that does not shift the frame

The awkward part is the minus strand, and *M. tuberculosis* has plenty of it:
``pncA`` and ``katG`` both run right to left. For those genes the c. coordinate
counts in the opposite direction from the genome coordinate **and** the alleles
are complemented, so a genome ``C>T`` inside ``pncA`` is written ``G>A``. Getting
that backwards produces a well-formed name for a mutation that does not exist.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..config import PROMOTER_UPSTREAM_BP
from ..records import Variant
from ..utils import (AA_1_TO_3, LOG, MjolnirError, PathLike, require_file,
                     revcomp, smart_open, translate)

# ---------------------------------------------------------------------------
# Gene models
# ---------------------------------------------------------------------------

#: Biotypes whose variants are named with ``n.`` rather than ``c.``/``p.``. An
#: rRNA gene has no reading frame, so a codon number would be a fiction — and
#: ``rrs`` and ``rrl`` are exactly where the aminoglycoside and macrolide
#: determinants live, in both MTBC and NTM.
NON_CODING_BIOTYPES = frozenset((
    "rrna", "trna", "ncrna", "misc_rna", "tmrna", "rrna_gene", "trna_gene",
))


@dataclass(frozen=True)
class Gene:
    """One gene, in the orientation the catalogues name its variants in."""

    name: str
    locus_tag: str
    start: int          #: 1-based, inclusive, genome orientation
    end: int            #: 1-based, inclusive, genome orientation
    strand: str         #: "+" or "-"
    biotype: str = "protein_coding"
    contig: str = ""

    @property
    def coding(self) -> bool:
        return self.biotype.lower() not in NON_CODING_BIOTYPES

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def label(self) -> str:
        """The name catalogues use: the gene symbol, or the locus tag without one."""
        return self.name or self.locus_tag

    def contains(self, position: int) -> bool:
        return self.start <= position <= self.end

    def promoter_contains(self, position: int, upstream: int) -> bool:
        """Whether *position* sits in the *upstream* window before the start codon.

        "Before" is strand-relative: for a minus-strand gene the promoter is at
        higher genome coordinates than the gene, not lower.
        """
        if self.strand == "+":
            return self.start - upstream <= position < self.start
        return self.end < position <= self.end + upstream

    def coding_offset(self, position: int) -> int:
        """The c. coordinate of *position*: 1-based in the gene's own direction.

        Negative inside the promoter, counting back from the first base of the
        start codon, and with no zero — HGVS goes ``c.-1`` then ``c.1``.
        """
        if self.strand == "+":
            return position - self.start + 1 if position >= self.start \
                else position - self.start
        return self.end - position + 1 if position <= self.end \
            else self.end - position


def _attributes(field_text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for chunk in field_text.split(";"):
        if "=" not in chunk:
            continue
        key, _, value = chunk.partition("=")
        out[key.strip().lower()] = value.strip()
    return out


def load_gff(path: PathLike, *, contig: str = "") -> List[Gene]:
    """Gene models from a GFF3.

    Only ``gene`` rows are read. CDS rows would give the same span for the
    bacterial genes here and would double-count, and Mjolnir does not model
    introns because these organisms do not have them.
    """
    resolved = require_file(path, "the GFF annotation")
    genes: List[Gene] = []
    with smart_open(resolved, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            kind = fields[2].lower()
            # rRNA_gene, tRNA_gene, ncRNA_gene as well as plain gene. Reading
            # only "gene" made rrs and rrl invisible, and those two carry the
            # aminoglycoside and macrolide determinants in both MTBC and NTM -
            # so every rrs_n.1401A>G in the catalogue went unmatched.
            if len(fields) < 9 or not (kind in ("gene", "pseudogene")
                                       or kind.endswith("_gene")):
                continue
            attrs = _attributes(fields[8])
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            locus = attrs.get("gene_id") or attrs.get("locus_tag") or attrs.get("id", "")
            locus = locus.split(":")[-1]
            genes.append(Gene(
                name=attrs.get("name", ""),
                locus_tag=locus,
                start=start, end=end,
                strand="-" if fields[6] == "-" else "+",
                biotype=attrs.get("biotype", "protein_coding"),
                contig=contig or fields[0],
            ))
    if not genes:
        raise MjolnirError(
            "{0} contained no gene records; expected a GFF3 with 'gene' rows "
            "(tbdb ships one as genome.gff)".format(resolved))
    genes.sort(key=lambda g: (g.start, g.end))
    LOG.debug("loaded %d gene models from %s", len(genes), resolved)
    return genes


def load_fasta(path: PathLike) -> Dict[str, str]:
    """Contig name to sequence. Needed for codons, and for nothing else."""
    resolved = require_file(path, "the reference FASTA")
    sequences: Dict[str, str] = {}
    name = ""
    chunks: List[str] = []
    with smart_open(resolved, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(chunks)
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name:
        sequences[name] = "".join(chunks)
    return sequences


# ---------------------------------------------------------------------------
# The index
# ---------------------------------------------------------------------------

class Annotation:
    """Gene models plus reference sequence: everything naming a variant needs.

    Lookup is a linear scan over genes whose span could contain the position,
    found by bisection. *M. tuberculosis* has about 4,000 genes and a run
    annotates a few thousand variants, so an interval tree would be optimising
    something that is already imperceptible.
    """

    def __init__(self, genes: Sequence[Gene], sequences: Optional[Dict[str, str]] = None,
                 *, promoter_bp: int = PROMOTER_UPSTREAM_BP) -> None:
        self.genes = list(genes)
        self.sequences = dict(sequences or {})
        self.promoter_bp = int(promoter_bp)
        self._starts = [g.start for g in self.genes]

    @classmethod
    def load(cls, gff: PathLike, fasta: Optional[PathLike] = None,
             *, contig: str = "", promoter_bp: int = PROMOTER_UPSTREAM_BP) -> "Annotation":
        sequences = load_fasta(fasta) if fasta is not None else {}
        return cls(load_gff(gff, contig=contig), sequences, promoter_bp=promoter_bp)

    def sequence_for(self, contig: str) -> str:
        """The reference sequence, tolerating a contig-name mismatch of one.

        tbdb calls H37Rv ``Chromosome`` and NCBI calls it ``NC_000962.3``. With a
        single-contig reference the name carries no information the caller does
        not already have, so a lone sequence is returned whatever it is called —
        and a genuine multi-contig mismatch still returns nothing rather than the
        wrong contig.
        """
        if contig in self.sequences:
            return self.sequences[contig]
        if len(self.sequences) == 1:
            return next(iter(self.sequences.values()))
        return ""

    def genes_at(self, position: int) -> List[Gene]:
        """Genes containing *position*, innermost first.

        Overlapping genes are real in these genomes, and a variant inside two of
        them is named for the smaller: that is the one whose reading frame the
        catalogues quote.
        """
        hits = [g for g in self.genes if g.contains(position)]
        hits.sort(key=lambda g: g.length)
        return hits

    def promoters_at(self, position: int) -> List[Gene]:
        """Genes whose promoter window contains *position*, nearest start first."""
        hits = [g for g in self.genes
                if g.promoter_contains(position, self.promoter_bp)]
        hits.sort(key=lambda g: abs(g.coding_offset(position)))
        return hits


# ---------------------------------------------------------------------------
# Naming
# ---------------------------------------------------------------------------

def _aa3(residue: str) -> str:
    return AA_1_TO_3.get(residue.upper(), "Xaa")


def _codon_of(gene: Gene, sequence: str, coding_position: int) -> Tuple[int, int, str]:
    """``(codon number, offset in codon, codon bases)`` for a c. position.

    The codon is read in the gene's own orientation, so for a minus-strand gene
    it is the reverse complement of what the genome shows.
    """
    codon_number = (coding_position - 1) // 3 + 1
    offset = (coding_position - 1) % 3
    if gene.strand == "+":
        first = gene.start + (codon_number - 1) * 3
        codon = sequence[first - 1:first + 2]
    else:
        last = gene.end - (codon_number - 1) * 3
        codon = revcomp(sequence[last - 3:last])
    return codon_number, offset, codon.upper()


def _substitute(codon: str, offset: int, base: str) -> str:
    return codon[:offset] + base.upper() + codon[offset + 1:]


def _hgvs_snv(gene: Gene, sequence: str, position: int,
              ref: str, alt: str) -> Tuple[str, str]:
    """``(hgvs, effect)`` for a single-base substitution inside or before a gene."""
    coding_position = gene.coding_offset(position)

    # Alleles are quoted in the gene's direction, so a minus-strand gene
    # complements them. A genome C>T inside pncA is written G>A.
    if gene.strand == "-":
        ref, alt = revcomp(ref), revcomp(alt)

    if coding_position < 0:
        # rrs and rrl keep the n. prefix upstream as well as inside: WHO files
        # the position 147 bases before rrs as rrs_n.-147G>T, not c.
        prefix = "c." if gene.coding else "n."
        return "{0}{1}{2}>{3}".format(prefix, coding_position, ref, alt), \
            "upstream_variant"

    if not gene.coding:
        return "n.{0}{1}>{2}".format(coding_position, ref, alt), "non_coding_variant"

    if not sequence:
        # No reference sequence: the nucleotide form is still correct and still
        # joins against the catalogues for synonymous rows. Inventing a protein
        # change without the codon would not be.
        return "c.{0}{1}>{2}".format(coding_position, ref, alt), "coding_variant"

    codon_number, offset, codon = _codon_of(gene, sequence, coding_position)
    if len(codon) != 3:
        return "c.{0}{1}>{2}".format(coding_position, ref, alt), "coding_variant"
    mutant = _substitute(codon, offset, alt)
    from_aa = translate(codon, start_is_met=False)
    to_aa = translate(mutant, start_is_met=False)
    if not from_aa or not to_aa:
        return "c.{0}{1}>{2}".format(coding_position, ref, alt), "coding_variant"
    if codon_number == 1:
        # WHO names start-codon changes in nucleotides (Rv0010c_c.3G>T) and pools
        # them under the gene's LoF, because what matters is that translation no
        # longer starts, not which residue the broken codon would encode.
        return "c.{0}{1}>{2}".format(coding_position, ref, alt), "start_lost"
    if from_aa == to_aa:
        return "c.{0}{1}>{2}".format(coding_position, ref, alt), "synonymous_variant"
    effect = "stop_gained" if to_aa == "*" else "missense_variant"
    return "p.{0}{1}{2}".format(_aa3(from_aa), codon_number, _aa3(to_aa)), effect


def _oriented(gene: Gene, bases: str) -> str:
    """Bases as the gene reads them: complemented and reversed on the minus strand."""
    return revcomp(bases) if gene.strand == "-" else bases


def _hgvs_indel(gene: Gene, position: int, ref: str, alt: str) -> Tuple[str, str]:
    """``(hgvs, effect)`` for an indel, in the gene's own direction.

    The names carry their bases - ``c.-85delG`` rather than ``c.-85del`` - and
    insertions carry the flanking range, because that is the form the catalogues
    file them under and a shorter name matches nothing.
    """
    shift = len(alt) - len(ref)
    first_changed = position + 1
    size = abs(shift)
    coding_first = gene.coding_offset(first_changed)

    if shift < 0:
        deleted = _oriented(gene, ref[1:1 + size])
        last = gene.coding_offset(first_changed + size - 1)
        low, high = min(coding_first, last), max(coding_first, last)
        span = "{0}".format(low) if size == 1 else "{0}_{1}".format(low, high)

        if not gene.coding:
            return "n.{0}del{1}".format(span, deleted), "non_coding_variant"

        # Classify by what the deletion actually removes, not by where it
        # starts. A deletion that begins upstream and runs into the gene takes
        # the start codon with it, and calling that an upstream variant is how a
        # complete pncA knockout - definitive pyrazinamide resistance - was
        # reported as a regulatory nucleotide change with no determinant at all.
        if high < 1:
            return "c.{0}del{1}".format(span, deleted), "upstream_variant"
        if low < 1:
            return "c.{0}del{1}".format(span, deleted), "start_lost"
        if shift % 3 != 0:
            codon_number = (low - 1) // 3 + 1
            return "p.Xaa{0}fs".format(codon_number), "frameshift_variant"
        return "c.{0}del{1}".format(span, deleted), "inframe_deletion"

    inserted = _oriented(gene, alt[1:1 + size])
    anchor = gene.coding_offset(position)
    left, right = (anchor, anchor + 1) if gene.strand == "+" else (anchor - 1, anchor)
    if not gene.coding:
        return "n.{0}_{1}ins{2}".format(left, right, inserted), "non_coding_variant"
    # The insertion lands after the anchor on the plus strand and before it on
    # the minus strand, so only the minus strand's first changed base is the
    # anchor itself. Using the anchor's own codon on the plus strand names a
    # codon the shift provably does not touch whenever the anchor sits on a
    # codon boundary - which flipped the rpoB RRDR rule in both directions, on
    # the one drug that rule exists for.
    first_coding = anchor + 1 if gene.strand == "+" else anchor
    if first_coding < 1:
        return "c.{0}_{1}ins{2}".format(left, right, inserted), "upstream_variant"
    if shift % 3 != 0:
        codon_number = (first_coding - 1) // 3 + 1
        return "p.Xaa{0}fs".format(codon_number), "frameshift_variant"
    return "c.{0}_{1}ins{2}".format(left, right, inserted), "inframe_insertion"


def _protein_indel(gene: Gene, sequence: str, position: int,
                   ref: str, alt: str) -> str:
    """WHO's protein-level name for an in-frame indel, or "".

    Built by mutating the coding sequence, translating it, and diffing against
    the reference protein — not by translating the inserted bases on their own.
    An insertion at a codon boundary translates cleanly; one in the middle of a
    codon does not, and *dnaA* c.65 is mid-codon. Translating ``CGA`` there gives
    Arg, while what the genome actually produces is a duplicated Asp, which is
    what WHO files it as.

    The diff is then shifted as far C-terminal as it will go, because HGVS names
    the 3'-most equivalent and equivalence is at the protein level: two codons
    encoding the same residue are interchangeable however different their bases.
    """
    shift = len(alt) - len(ref)
    if not sequence or not gene.coding or shift == 0 or shift % 3 != 0:
        return ""
    reference = _protein(gene, sequence)
    if not reference:
        return ""

    coding = sequence[gene.start - 1:gene.end]
    offset = position - gene.start
    if offset < 0 or offset + len(ref) > len(coding):
        return ""
    if coding[offset:offset + len(ref)].upper() != ref.upper():
        return ""
    mutated_dna = coding[:offset] + alt + coding[offset + len(ref):]
    if gene.strand == "-":
        mutated_dna = revcomp(mutated_dna)
    mutant = translate(mutated_dna, start_is_met=False)
    if not mutant:
        return ""

    # Trim the identical head and tail; what is left is the event.
    head = 0
    limit = min(len(reference), len(mutant))
    while head < limit and reference[head] == mutant[head]:
        head += 1
    tail = 0
    while (tail < limit - head
           and reference[len(reference) - 1 - tail] == mutant[len(mutant) - 1 - tail]):
        tail += 1

    removed = reference[head:len(reference) - tail]
    added = mutant[head:len(mutant) - tail]

    if removed and not added:
        first = head + 1
        last = head + len(removed)
        while (last < len(reference)
               and reference[first - 1] == reference[last]):
            first += 1
            last += 1
        if len(removed) == 1:
            return "p.{0}{1}del".format(_aa3(reference[first - 1]), first)
        return "p.{0}{1}_{2}{3}del".format(
            _aa3(reference[first - 1]), first,
            _aa3(reference[last - 1]), last)

    if added and not removed:
        size = len(added)
        before = head
        while (before < len(reference)
               and reference[before:before + size] == added):
            before += 1
        start_res = before - size + 1
        if start_res >= 1 and reference[start_res - 1:before] == added:
            if size == 1:
                return "p.{0}{1}dup".format(_aa3(added[0]), before)
            return "p.{0}{1}_{2}{3}dup".format(
                _aa3(added[0]), start_res, _aa3(added[-1]), before)
        if before < 1 or before >= len(reference):
            return ""
        return "p.{0}{1}_{2}{3}ins{4}".format(
            _aa3(reference[before - 1]), before,
            _aa3(reference[before]), before + 1,
            "".join(_aa3(r) for r in added))
    return ""


def _frameshift_name(gene: Gene, sequence: str, position: int,
                     ref: str, alt: str) -> str:
    """``p.<Aa><n>fs`` at the first residue the frameshift actually changes.

    Not at the codon the inserted or deleted base happens to fall in. WHO names
    ``dnaA_p.Glu486fs`` where the nucleotide sits in codon 505: a frameshift can
    be introduced in a run of bases that leaves the next several residues intact,
    and the name marks where the protein stops matching, which is what a reader
    needs and what the catalogue is keyed on.
    """
    shift = len(alt) - len(ref)
    if not sequence or not gene.coding or shift == 0:
        return ""
    reference = _protein(gene, sequence)
    if not reference:
        return ""
    coding = sequence[gene.start - 1:gene.end]
    offset = position - gene.start
    if offset < 0 or offset + len(ref) > len(coding):
        return ""
    if coding[offset:offset + len(ref)].upper() != ref.upper():
        return ""
    mutated = coding[:offset] + alt + coding[offset + len(ref):]
    if gene.strand == "-":
        mutated = revcomp(mutated)
    mutant = translate(mutated, start_is_met=False)
    if not mutant:
        return ""
    # An in-frame indel is not a frameshift - unless what it inserts carries a
    # stop codon, which truncates the product just as a frameshift does. WHO
    # files that as fs too, and a 24-base insertion containing TAG in dnaA is
    # exactly such a case.
    truncated = len(mutant) < len(reference)
    if shift % 3 == 0 and not truncated:
        return ""
    limit = min(len(reference), len(mutant))
    index = 0
    while index < limit and reference[index] == mutant[index]:
        index += 1
    if index >= len(reference):
        return ""
    return "p.{0}{1}fs".format(_aa3(reference[index]), index + 1)


def _protein(gene: Gene, sequence: str) -> str:
    """The gene's translated product, cached on first use."""
    cached = _PROTEIN_CACHE.get((gene.locus_tag, gene.start, gene.end))
    if cached is not None:
        return cached
    coding = sequence[gene.start - 1:gene.end]
    if gene.strand == "-":
        coding = revcomp(coding)
    product = translate(coding, start_is_met=False)
    _PROTEIN_CACHE[(gene.locus_tag, gene.start, gene.end)] = product
    return product


_PROTEIN_CACHE: Dict[Tuple[str, int, int], str] = {}


#: Effects that pool under WHO's ``<gene>_LoF`` graded variant. WHO does not
#: grade every frameshift individually; it grades the loss of the gene, and the
#: Genomic_coordinates sheet maps each specific frameshift and large deletion to
#: that one name. A tool emitting only the precise coordinate name matches no
#: catalogue row at all for the mutations that matter most in pncA and katG.
LOF_EFFECTS = frozenset((
    "frameshift_variant", "stop_gained", "start_lost", "large_deletion",
    "stop_lost",
))

#: A deletion removing at least this much coding sequence is a candidate loss of
#: that gene. SOURCE: Mjolnir policy, calibrated against WHO's own coordinate
#: sheet - WHO files deletions under ``<gene>_LoF`` even when they remove only
#: part of the gene and even when they also remove a neighbour, so a fractional
#: threshold missed 16,622 of its rows. This produces a *candidate* name for
#: matching, never a printed call.
LOF_DELETION_BP = 1


def _decompose(position: int, ref: str, alt: str) -> List[Tuple[int, str, str]]:
    """An equal-length substitution as its individual changed bases.

    ``AT>CA`` at position 8 is two substitutions, not an indel. WHO names such a
    variant both by its constituent nucleotide changes and, when they fall in one
    codon, by the single amino-acid change they produce together - so a caller
    that treats it as an indel matches neither.
    """
    return [(position + index, ref[index], alt[index])
            for index in range(len(ref)) if ref[index] != alt[index]]


def equivalent_placements(sequence: str, position: int, ref: str, alt: str,
                          window: int = 24) -> List[Tuple[int, str, str]]:
    """Every genomic placement of an indel that yields the identical sequence.

    ``ACT`` deleted from ``CACTACT`` can be written at two offsets and both
    describe the same genome. HGVS resolves this by naming the 3'-most, but
    "3'-most" is relative to the *gene*, so a minus-strand gene picks the other
    one - and a catalogue lookup that offers only Mjolnir's choice misses every
    row filed under the other. Rather than reproduce the shifting rule and its
    strand exception, this enumerates the equivalent placements and lets the
    caller offer a name for each.

    Returns ``[(position, ref, alt), ...]``, the input included, or just the
    input when the variant is a substitution or sits at a contig edge.
    """
    if len(ref) == len(alt) or not sequence:
        return [(position, ref, alt)]
    start = max(1, position - window)
    stop = min(len(sequence), position + len(ref) + window)
    if stop <= start:
        return [(position, ref, alt)]
    original = sequence[start - 1:stop]
    offset = position - start
    if sequence[position - 1:position - 1 + len(ref)].upper() != ref.upper():
        # The reference allele does not match the genome; shifting a variant we
        # cannot place would invent equivalences.
        return [(position, ref, alt)]
    altered = original[:offset] + alt + original[offset + len(ref):]

    out: List[Tuple[int, str, str]] = []
    shift = len(alt) - len(ref)
    for candidate in range(start, stop - abs(shift)):
        index = candidate - start
        if shift < 0:
            new_ref = original[index:index + 1 - shift]
            new_alt = original[index:index + 1]
        else:
            new_ref = original[index:index + 1]
            new_alt = new_ref + alt[1:]
        if len(new_ref) < 1 or len(new_alt) < 1:
            continue
        rebuilt = original[:index] + new_alt + original[index + len(new_ref):]
        if rebuilt == altered:
            out.append((candidate, new_ref, new_alt))
    return out or [(position, ref, alt)]


def names_for(annotation: Annotation, contig: str, position: int,
              ref: str, alt: str) -> List[str]:
    """Every ``<gene>_<hgvs>`` name this variant legitimately answers to.

    A catalogue lookup should match a variant against any name the catalogues
    could have filed it under, not only the one Mjolnir would print. One
    multi-nucleotide variant maps to several graded variants in WHO's own
    coordinate sheet - the ampersand-joined names in its VCF - and a frameshift
    maps to both its own coordinate name and the gene's pooled ``_LoF``.
    """
    names: List[str] = []

    def _add(gene_label: str, hgvs: str) -> None:
        if gene_label and hgvs:
            candidate = "{0}_{1}".format(gene_label, hgvs)
            if candidate not in names:
                names.append(candidate)

    sequence = annotation.sequence_for(contig)
    primary_gene, _locus, primary_hgvs, effect = name_variant(
        annotation, contig, position, ref, alt)
    _add(primary_gene, primary_hgvs)

    # A position can sit inside one gene and upstream of another - 1471699 is
    # inside the ncRNA mcr3 and 147 bases before rrs, and WHO files it as
    # rrs_n.-147G>T. Both names are offered rather than one chosen.
    for gene in annotation.promoters_at(position):
        if gene.label == primary_gene:
            continue
        if len(ref) == len(alt):
            # Equal lengths are a substitution, not an indel. Routing an MNV
            # through the indel namer produced "n.-41_-40ins" with nothing
            # inserted - a name for an event that did not happen.
            changed = _decompose(position, ref, alt)
            if not changed:
                continue
            pos, ref_base, alt_base = changed[0]
            hgvs, _effect = _hgvs_snv(gene, sequence, pos, ref_base, alt_base)
        else:
            hgvs, _effect = _hgvs_indel(gene, position, ref, alt)
        _add(gene.label, hgvs)

    # In-frame indels are filed at protein level by WHO, and the same indel can
    # be written at several genomic offsets - so every equivalent placement is
    # named, at both nucleotide and protein level.
    if len(ref) != len(alt):
        for placed_pos, placed_ref, placed_alt in equivalent_placements(
                sequence, position, ref, alt):
            for gene in (annotation.genes_at(placed_pos)
                         or annotation.promoters_at(placed_pos)):
                nucleotide, _effect = _hgvs_indel(
                    gene, placed_pos, placed_ref, placed_alt)
                _add(gene.label, nucleotide)
                _add(gene.label, _protein_indel(
                    gene, sequence, placed_pos, placed_ref, placed_alt))
                frameshift = _frameshift_name(
                    gene, sequence, placed_pos, placed_ref, placed_alt)
                if frameshift:
                    _add(gene.label, frameshift)
                    _add(gene.label, "LoF")

    # Equal-length multi-base substitution: name every changed base, and every
    # codon those bases land in.
    if len(ref) == len(alt) > 1:
        codons: Dict[Tuple[str, int], Gene] = {}
        for pos, ref_base, alt_base in _decompose(position, ref, alt):
            genes = annotation.genes_at(pos) or annotation.promoters_at(pos)
            if not genes:
                continue
            gene = genes[0]
            single_hgvs, single_effect = _hgvs_snv(
                gene, sequence, pos, ref_base, alt_base)
            _add(gene.label, single_hgvs)
            # The same base can be upstream of a different gene, which is the
            # one WHO files it under (fgd1_c.-39T>A, not Rv0406c_c.-39A>T).
            for upstream in annotation.promoters_at(pos):
                if upstream.label == gene.label:
                    continue
                other, _e = _hgvs_snv(upstream, sequence, pos, ref_base, alt_base)
                _add(upstream.label, other)
            if gene.coding and sequence and single_effect != "upstream_variant":
                coding_position = gene.coding_offset(pos)
                if coding_position > 0:
                    codons[(gene.label, (coding_position - 1) // 3 + 1)] = gene
        for (label, codon_number), gene in codons.items():
            protein = _codon_change(gene, sequence, codon_number, position, ref, alt)
            if protein:
                _add(label, protein)
                if protein.endswith("Ter") or codon_number == 1:
                    _add(label, "LoF")
        # WHO files a run of substituted bases as a delins as well.
        first_gene = (annotation.genes_at(position) or [None])[0]
        if first_gene is not None:
            low = first_gene.coding_offset(position)
            high = first_gene.coding_offset(position + len(ref) - 1)
            prefix = "c." if first_gene.coding else "n."
            _add(first_gene.label, "{0}{1}_{2}del{3}ins{4}".format(
                prefix, min(low, high), max(low, high),
                _oriented(first_gene, ref), _oriented(first_gene, alt)))

    # Loss of function pools under the gene's own name.
    if effect in LOF_EFFECTS and primary_gene:
        _add(primary_gene, "LoF")
    if primary_hgvs.startswith("p.") and primary_hgvs.endswith("Ter") and primary_gene:
        _add(primary_gene, "LoF")
    # Loss of function is pooled per gene, and adjacent genes overlap in these
    # genomes: a change in the stop codon of mmpS5 is filed by WHO as
    # mmpL5_LoF, because the pair is what the efflux phenotype needs.
    disruptive = (effect in LOF_EFFECTS
                  or (primary_hgvs.startswith("p.") and primary_hgvs.endswith("Ter"))
                  # Ter143Cys is a stop *lost*: translation runs past the end of
                  # the gene, which WHO pools under the operon's LoF just as it
                  # pools a premature stop.
                  or re.match(r"^p\.Ter\d+", primary_hgvs) is not None)
    if disruptive:
        for gene in (annotation.genes_at(position) + annotation.promoters_at(position)):
            if gene.coding:
                _add(gene.label, "LoF")
    for gene in _genes_lost(annotation, position, ref, alt):
        _add(gene.label, "LoF")
    return names


def _codon_change(gene: Gene, sequence: str, codon_number: int,
                  position: int, ref: str, alt: str) -> str:
    """The protein name for a codon carrying several substituted bases at once."""
    if not sequence or not gene.coding:
        return ""
    _n, _offset, codon = _codon_of(gene, sequence, (codon_number - 1) * 3 + 1)
    if len(codon) != 3:
        return ""
    mutant = list(codon)
    for pos, _ref_base, alt_base in _decompose(position, ref, alt):
        if not gene.contains(pos):
            continue
        coding_position = gene.coding_offset(pos)
        if coding_position <= 0 or (coding_position - 1) // 3 + 1 != codon_number:
            continue
        base = revcomp(alt_base) if gene.strand == "-" else alt_base
        mutant[(coding_position - 1) % 3] = base.upper()
    from_aa = translate(codon, start_is_met=False)
    to_aa = translate("".join(mutant), start_is_met=False)
    if not from_aa or not to_aa or from_aa == to_aa:
        return ""
    return "p.{0}{1}{2}".format(_aa3(from_aa), codon_number, _aa3(to_aa))


def _genes_lost(annotation: Annotation, position: int,
                ref: str, alt: str) -> List[Gene]:
    """Coding genes a deletion removes enough of to be called a loss."""
    if len(alt) >= len(ref):
        return []
    first = position + 1
    last = position + (len(ref) - len(alt))
    lost: List[Gene] = []
    for gene in annotation.genes:
        if not gene.coding or gene.end < first or gene.start > last:
            continue
        overlap = min(gene.end, last) - max(gene.start, first) + 1
        if overlap >= LOF_DELETION_BP:
            lost.append(gene)
    return lost


def name_variant(annotation: Annotation, contig: str, position: int,
                 ref: str, alt: str) -> Tuple[str, str, str, str]:
    """``(gene, locus_tag, hgvs, effect)`` for one variant, or empties.

    A variant inside a gene is named for that gene. One that is not, but sits in
    a promoter window, is named for the gene it regulates — which is where
    ``eis_c.-14C>T`` and ``inhA_c.-154G>A`` come from, both Group 1 determinants
    that a coding-only annotator would drop on the floor.
    """
    sequence = annotation.sequence_for(contig)
    genes = annotation.genes_at(position)
    if not genes:
        genes = annotation.promoters_at(position)
    if not genes:
        return "", "", "", "intergenic_variant"
    gene = genes[0]

    if len(ref) == 1 and len(alt) == 1:
        hgvs, effect = _hgvs_snv(gene, sequence, position, ref, alt)
    elif len(ref) == len(alt):
        # A substitution of several bases at once, not an indel. Named for the
        # first base that actually changed; names_for() supplies the rest.
        changed = _decompose(position, ref, alt)
        if not changed:
            return gene.label, gene.locus_tag, "", "no_change"
        pos, ref_base, alt_base = changed[0]
        hgvs, effect = _hgvs_snv(gene, sequence, pos, ref_base, alt_base)
    else:
        hgvs, effect = _hgvs_indel(gene, position, ref, alt)
    return gene.label, gene.locus_tag, hgvs, effect


def annotate(variants: Iterable[Variant], annotation: Annotation) -> int:
    """Attach gene, locus tag, HGVS and effect to each variant, in place.

    Returns how many were named. A variant that is already named is left alone,
    so a caller that got names from somewhere better does not lose them.
    """
    named = 0
    for variant in variants:
        if variant.gene and variant.hgvs:
            named += 1
            continue
        gene, locus, hgvs, effect = name_variant(
            annotation, variant.chrom, variant.pos, variant.ref, variant.alt)
        if not hgvs:
            variant.effect = variant.effect or effect
            continue
        variant.gene = gene
        variant.locus_tag = variant.locus_tag or locus
        variant.hgvs = hgvs
        variant.effect = effect
        named += 1
    return named


__all__ = [
    "Annotation", "Gene", "annotate", "load_fasta", "load_gff", "name_variant",
    "names_for", "LOF_EFFECTS", "NON_CODING_BIOTYPES",
]
