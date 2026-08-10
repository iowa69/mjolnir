"""Turning three catalogues' spellings of the same mutation into one key.

Every observed variant leaves this module with **two** identities, because the
design needs both and neither can do the other's job.

The **coordinate key** — ``(NC_000962.3, POS, REF, ALT)``, left-aligned and
parsimonious — is WHO's own documented matching protocol. WHO grades a variant
name, then publishes a separate table saying which genomic changes produce that
name; matching on the name alone would miss every alternative codon spelling
they enumerate (``dnaA_p.Gly6Ser`` has six distinct REF/ALT pairs at position
16). So the coordinate is the primary path against WHO.

The **HGVS key** — ``<gene>_<hgvs>`` with three-letter amino acids — is the
cross-catalogue join key, because MTBseq and tbdb do not share WHO's coordinate
table and never will. It is the only string all three sources can be made to
agree on.

The third thing here is the alias table, and it exists to stop Mjolnir
manufacturing doubt. ``rpoB_p.Ser450Leu`` and ``rpoB_p.Ser531Leu`` are the same
mutation written in *M. tuberculosis* and *Escherichia coli* codon numbering;
``rrs n.1401A>G`` and ``rrs`` position 1408 are the same substitution in
mycobacterial and *E. coli* 16S numbering. A tool that reports those as a
catalogue disagreement is telling a clinician the sources conflict when they do
not. :func:`classify_difference` therefore returns ``nomenclature`` rather than
``biological`` for exactly these cases, and every alias carries the source that
establishes it.

Left alignment needs the reference sequence. When one is not supplied the
alleles are still made parsimonious, but an indel inside a homopolymer or a
tandem repeat cannot be shifted to its canonical position — so
:class:`NormalisedVariant` records ``left_aligned=False`` instead of pretending,
and a coordinate match that fails for that reason is visible rather than silent.

Thresholds live in ``config.py``. What lives here is catalogue *nomenclature*
data — offsets between two published numbering schemes — and each entry carries
its own source string and ``verified`` flag in the same spirit, so a number in
this file is no more anonymous than one in the registry.
"""

from __future__ import annotations

import re
from collections.abc import Mapping as _Mapping
from dataclasses import dataclass, replace
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from ..config import (
    H37RV_ACCESSION,
    SRC_PRAMMANANAN_1998,
    SRC_WHO_V2,
)
from ..records import (
    DISAGREEMENT_BIOLOGICAL,
    DISAGREEMENT_NOMENCLATURE,
    DISAGREEMENT_NONE,
    VARIANT_DEL,
    VARIANT_INDEL,
    VARIANT_INS,
    VARIANT_MNV,
    VARIANT_SNP,
    Variant,
)
from ..utils import AA_1_TO_3, AA_3_TO_1, LOG, MjolnirError

#: ``(chromosome, position, reference allele, alternative allele)`` — WHO's
#: matching key, and the tuple :attr:`records.Variant.coordinate_key` produces.
CoordinateKey = Tuple[str, int, str, str]

#: A reference sequence, in any of the three shapes a caller is likely to have:
#: a ``{chrom: sequence}`` mapping, a pysam ``FastaFile`` (or anything else with
#: a 0-based half-open ``fetch``), or a plain callable taking 1-based inclusive
#: coordinates.
Reference = Any

# ---------------------------------------------------------------------------
# Graded variant names that are rules rather than coordinates
# ---------------------------------------------------------------------------
#
# SOURCE: WHO-UCN-TB-2023.7 Catalogue_master_file `variant` column. Pooled
# graded variants are named `<gene>_LoF` and `<gene>_deletion`. The design says
# no coordinates exist for these; that is true of `<gene>_deletion` but *not*
# universally of LoF — `fgd1_LoF` carries 1,487 coordinate rows in
# Genomic_coordinates_7May2024. So Mjolnir marks these names as rule-matched and
# still indexes whatever coordinates the file supplies, rather than assuming the
# coordinate table is empty for them.

RULE_LOF = "LoF"
RULE_DELETION = "deletion"

#: Rule-shaped variant names, from all three sources. tbdb writes ``frameshift``
#: and ``large_deletion`` in the same column where it otherwise writes HGVS.
_RULE_NAMES = ("lof", "deletion", "frameshift", "large_deletion", "any_indel")
_RULE_PATTERNS = (
    re.compile(r"^any_missense_codon_\d+$", re.IGNORECASE),
    re.compile(r"^any_indel_(codon|nucleotide)_\d+$", re.IGNORECASE),
)

_HGVS_PREFIXES = ("p.", "c.", "n.", "r.", "g.", "m.")

_AA3 = "|".join(sorted(AA_1_TO_3.values()))
_AA1 = "".join(sorted(k for k in AA_1_TO_3 if k.isalpha()))

#: An amino acid immediately before a codon number: ``Ser450``, ``S450``.
_AA_BEFORE_NUMBER = re.compile(r"(?P<aa>" + _AA3 + r"|[" + _AA1 + r"*])(?=-?\d)")
#: An amino acid immediately after a codon number: the ``Leu`` of ``Ser450Leu``.
_AA_AFTER_NUMBER = re.compile(r"(?<=\d)(?P<aa>" + _AA3 + r"|[" + _AA1 + r"*])(?![a-z])")

#: A bare protein change with no ``p.`` prefix, as legacy MTBseq rows write it.
_BARE_PROTEIN = re.compile(
    r"^(?:" + _AA3 + r"|[" + _AA1 + r"*])-?\d+(?:" + _AA3 + r"|[" + _AA1 + r"*]|fs|\*|=)?$"
)
#: A bare nucleotide change: MTBseq writes ``-14c>t`` and ``1401a>g``.
_BARE_NUCLEOTIDE = re.compile(r"^-?\d+[ACGTUacgtu]+>[ACGTUacgtu]+$")

_SUBSTITUTION_BASES = re.compile(r"(?P<ref>[ACGTUNacgtun]+)>(?P<alt>[ACGTUNacgtun]+)")
_INDEL_BASES = re.compile(r"(?P<op>ins|del|dup|inv)(?P<bases>[ACGTUNacgtun]+)")

_PROTEIN_SPLIT = re.compile(r"^(?P<prefix>p\.)(?P<ref>" + _AA3 + r")(?P<num>\d+)(?P<rest>.*)$")
_NUCLEOTIDE_SPLIT = re.compile(r"^(?P<prefix>[cnrgm]\.)(?P<num>-?\d+)(?P<rest>.*)$")


# ---------------------------------------------------------------------------
# Coordinates
# ---------------------------------------------------------------------------

def variant_class(ref: str, alt: str) -> str:
    """Classify an allele pair into the :mod:`records` variant vocabulary.

    Length equality decides substitution versus indel, and length one decides
    SNP versus MNV. WHO enumerates multi-nucleotide spellings of a single codon
    change as first-class rows, so ``GGT>TCC`` must stay one MNV rather than
    being split into three SNPs that the catalogue does not contain.
    """
    ref_up = (ref or "").upper()
    alt_up = (alt or "").upper()
    if len(ref_up) == len(alt_up):
        return VARIANT_SNP if len(ref_up) == 1 else VARIANT_MNV
    if len(ref_up) < len(alt_up):
        return VARIANT_INS if ref_up and alt_up.startswith(ref_up) else VARIANT_INDEL
    return VARIANT_DEL if alt_up and ref_up.startswith(alt_up) else VARIANT_INDEL


def trim_alleles(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Make an allele pair parsimonious without touching the reference.

    Shared bases are removed from the right first and then from the left, and
    never to the point of emptying either allele — VCF requires both to keep at
    least one base. This is the half of normalisation that needs no reference
    genome; :func:`left_align` is the half that does.
    """
    ref_up = (ref or "").upper()
    alt_up = (alt or "").upper()
    pos = int(pos)
    while len(ref_up) > 1 and len(alt_up) > 1 and ref_up[-1] == alt_up[-1]:
        ref_up, alt_up = ref_up[:-1], alt_up[:-1]
    while len(ref_up) > 1 and len(alt_up) > 1 and ref_up[0] == alt_up[0]:
        ref_up, alt_up = ref_up[1:], alt_up[1:]
        pos += 1
    return pos, ref_up, alt_up


def _fetcher(reference: Reference) -> Optional[Callable[[str, int, int], str]]:
    """Adapt whatever reference the caller has into 1-based inclusive fetching.

    Three shapes are accepted because three are in use: a ``{chrom: sequence}``
    dict in tests, a ``pysam.FastaFile`` in a real run, and a bare callable when
    another module has already wrapped one. Anything else raises rather than
    being ignored, since a silently unused reference would leave indels
    unaligned while the record still claimed they were.
    """
    if reference is None:
        return None
    if callable(reference) and not hasattr(reference, "fetch"):
        return reference  # type: ignore[return-value]
    fetch = getattr(reference, "fetch", None)
    if fetch is not None:
        def _from_pysam(chrom: str, start: int, end: int) -> str:
            # pysam's fetch is 0-based half-open; ours is 1-based inclusive.
            return str(fetch(chrom, start - 1, end)).upper()
        return _from_pysam
    if isinstance(reference, _Mapping):
        def _from_mapping(chrom: str, start: int, end: int) -> str:
            seq = reference.get(chrom)
            if seq is None:
                raise MjolnirError(
                    "reference has no sequence named {0!r}; it carries {1}".format(
                        chrom, ", ".join(sorted(reference)) or "nothing")
                )
            return str(seq[start - 1:end]).upper()
        return _from_mapping
    raise MjolnirError(
        "cannot read a reference sequence from a {0}; pass a {{chrom: sequence}} "
        "mapping, a pysam FastaFile, or a callable(chrom, start, end) using "
        "1-based inclusive coordinates".format(type(reference).__name__)
    )


def left_align(chrom: str, pos: int, ref: str, alt: str,
               reference: Reference = None) -> Tuple[int, str, str, bool]:
    """Shift an allele pair to its leftmost equivalent representation.

    This is vt's algorithm: trim from the right while the last bases agree,
    extend one base to the left whenever an allele would otherwise be emptied,
    then trim from the left. Two callers writing the same deletion as
    ``AGG>A`` at 100 and ``GGA>GA`` at 101 produce the same key afterwards, and
    without that a catalogue lookup on an indel is a coin toss.

    Returns ``(pos, ref, alt, left_aligned)``. The flag is False when no
    reference was supplied and the alleles still share a terminal base, because
    that is precisely the case where the answer might have moved and did not.
    """
    fetch = _fetcher(reference)
    pos = int(pos)
    ref_up = (ref or "").upper()
    alt_up = (alt or "").upper()
    if not ref_up or not alt_up:
        raise MjolnirError(
            "variant at {0}:{1} has an empty allele ({2!r}>{3!r}); VCF alleles "
            "must both carry at least one base".format(chrom, pos, ref, alt)
        )

    left_aligned = True
    guard = 0
    while ref_up[-1] == alt_up[-1]:
        if len(ref_up) > 1 and len(alt_up) > 1:
            ref_up, alt_up = ref_up[:-1], alt_up[:-1]
        elif fetch is not None and pos > 1:
            base = fetch(chrom, pos - 1, pos - 1)
            if not base:
                left_aligned = False
                break
            pos -= 1
            ref_up, alt_up = base + ref_up, base + alt_up
        else:
            # A shared terminal base with nothing to extend into: either this is
            # the start of the contig, or we have no reference. Say which.
            left_aligned = pos <= 1
            break
        guard += 1
        if guard > 1000:
            raise MjolnirError(
                "left-alignment of the variant at {0}:{1} did not converge after "
                "1000 shifts; the reference sequence does not match these "
                "alleles".format(chrom, pos)
            )
    while len(ref_up) > 1 and len(alt_up) > 1 and ref_up[0] == alt_up[0]:
        ref_up, alt_up = ref_up[1:], alt_up[1:]
        pos += 1
    return pos, ref_up, alt_up, left_aligned


def normalise_coordinate(chrom: str, pos: int, ref: str, alt: str,
                         reference: Reference = None) -> CoordinateKey:
    """The coordinate key WHO's protocol matches on."""
    new_pos, new_ref, new_alt, _ = left_align(chrom, pos, ref, alt, reference)
    return (chrom or H37RV_ACCESSION, int(new_pos), new_ref, new_alt)


def coordinate_string(key: CoordinateKey) -> str:
    """``NC_000962.3:761155C>T`` — the annex spelling of a coordinate key."""
    return "{0}:{1}{2}>{3}".format(key[0], key[1], key[2], key[3])


def parse_coordinate_string(text: str) -> CoordinateKey:
    """Inverse of :func:`coordinate_string`, for reading back an artefact."""
    match = re.match(r"^(?P<chrom>[^:]+):(?P<pos>\d+)(?P<ref>[A-Za-z]+)>(?P<alt>[A-Za-z]+)$",
                     str(text or "").strip())
    if not match:
        raise MjolnirError(
            "{0!r} is not a coordinate key; expected CHROM:POS<REF>><ALT>, "
            "for example NC_000962.3:761155C>T".format(text)
        )
    return (match.group("chrom"), int(match.group("pos")),
            match.group("ref").upper(), match.group("alt").upper())


# ---------------------------------------------------------------------------
# HGVS
# ---------------------------------------------------------------------------

def is_rule_variant(hgvs: str) -> bool:
    """Whether this graded-variant name is matched by rule, not by coordinate.

    ``katG_LoF``, ``pncA_deletion`` and tbdb's ``frameshift`` describe a class of
    change rather than one allele. They are still real catalogue entries; they
    just cannot be looked up in a coordinate table, and a loader that treated
    them as ordinary HGVS would silently drop the pooled loss-of-function
    grades that §5.4 depends on.
    """
    text = str(hgvs or "").strip()
    if not text:
        return False
    tail = text.rsplit("_", 1)[-1].lower()
    if tail in _RULE_NAMES or text.lower() in _RULE_NAMES:
        return True
    return any(pattern.match(text) for pattern in _RULE_PATTERNS)


def three_letter(hgvs: str) -> str:
    """Rewrite every amino acid in an HGVS string in three-letter form.

    The design fixes three-letter as the cross-catalogue join key. tbdb and WHO
    already use it; MTBseq's older rows and several downstream tools write
    ``S450L``, and joining those two spellings by string equality yields nothing.
    """
    text = str(hgvs or "")
    if not text:
        return ""

    def _expand(match: "re.Match") -> str:
        aa = match.group("aa")
        if aa in AA_3_TO_1:
            return aa
        return AA_1_TO_3.get(aa.upper(), aa)

    return _AA_AFTER_NUMBER.sub(_expand, _AA_BEFORE_NUMBER.sub(_expand, text))


def one_letter(hgvs: str) -> str:
    """The one-letter spelling, kept only so an alias lookup can find it."""
    text = str(hgvs or "")
    if not text:
        return ""

    def _shrink(match: "re.Match") -> str:
        aa = match.group("aa")
        return AA_3_TO_1.get(aa, aa)

    return _AA_AFTER_NUMBER.sub(_shrink, _AA_BEFORE_NUMBER.sub(_shrink, text))


def normalise_hgvs(hgvs: str, default_prefix: str = "c.") -> str:
    """One canonical spelling of a variant name, whichever source wrote it.

    Bases are upper-cased (MTBseq writes ``-14c>t``), amino acids are expanded to
    three letters, and a bare change is given the prefix its context implies —
    ``default_prefix`` is ``n.`` for an rRNA gene and ``c.`` for a coding one, a
    distinction the caller knows and this function cannot.

    Anything unrecognised is returned stripped but otherwise untouched. Guessing
    at a name Mjolnir does not understand would produce a key that joins against
    the wrong catalogue row, which is worse than a key that joins against none.
    """
    text = " ".join(str(hgvs or "").split())
    if not text:
        return ""
    if is_rule_variant(text):
        return text

    prefix = ""
    body = text
    for candidate in _HGVS_PREFIXES:
        if text.lower().startswith(candidate):
            prefix, body = text[:len(candidate)].lower(), text[len(candidate):]
            break
    if not prefix:
        # A bare change is ambiguous in one direction only: ``A1401G`` could be
        # Ala1401Gly or a nucleotide substitution, and the amino-acid reading
        # wins here because that is the form legacy MTBseq rows use. None of the
        # three catalogues actually writes an unprefixed nucleotide change in
        # that shape — MTBseq writes ``1401a>g`` — so the ambiguity is recorded
        # rather than defended against, and callers with a coding context should
        # pass prefixed HGVS.
        if _BARE_PROTEIN.match(body):
            prefix = "p."
        elif _BARE_NUCLEOTIDE.match(body):
            prefix = default_prefix
        else:
            return text

    if prefix == "p.":
        return "p." + three_letter(body)

    body = _SUBSTITUTION_BASES.sub(
        lambda m: "{0}>{1}".format(m.group("ref").upper(), m.group("alt").upper()), body)
    body = _INDEL_BASES.sub(
        lambda m: "{0}{1}".format(m.group("op"), m.group("bases").upper()), body)
    return prefix + body


def hgvs_key(gene: str, hgvs: str) -> str:
    """``<gene>_<hgvs>``, or "" when either half is missing.

    Empty is a real answer — an intergenic position outside every catalogued
    region has no such key — and it must not be faked with a coordinate string,
    or it would join against nothing while looking like a missing grade.
    """
    gene_name = str(gene or "").strip()
    variant = str(hgvs or "").strip()
    if not gene_name or not variant:
        return ""
    return "{0}_{1}".format(gene_name, variant)


def _looks_like_variant_part(text: str) -> bool:
    lowered = text.lower()
    if lowered in _RULE_NAMES:
        return True
    if any(lowered.startswith(prefix) for prefix in _HGVS_PREFIXES):
        return True
    return any(pattern.match(text) for pattern in _RULE_PATTERNS)


def split_key(key: str) -> Tuple[str, str]:
    """Split ``rpoB_p.Ser450Leu`` into ``("rpoB", "p.Ser450Leu")``.

    The split point is the *first* underscore whose right-hand side looks like a
    variant name, not the first or the last underscore in the string. Both naive
    choices break on real keys: gene names contain underscores
    (``mmpS5_mmpL5_c.-74G>A``, which ``split("_", 1)`` would amputate) and so do
    rule names (``katG_any_missense_codon_450``, which ``rpartition("_")`` would
    cut at the codon number).
    """
    text = str(key or "").strip()
    if not text or "_" not in text:
        return ("", text)
    for index, char in enumerate(text):
        if char != "_" or index == 0:
            continue
        if _looks_like_variant_part(text[index + 1:]):
            return (text[:index], text[index + 1:])
    head, _, tail = text.rpartition("_")
    return (head, tail)


# ---------------------------------------------------------------------------
# Legacy numbering aliases
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NumberingAlias:
    """An offset between two published numbering schemes for one gene.

    ``offset`` is added to the canonical (M. tuberculosis) number to obtain the
    legacy one. The validity range is not decoration: these offsets are
    alignment artefacts, constant only over the stretch of the gene where the
    two organisms' sequences align without an indel, so applying one genome-wide
    would invent equivalences outside the region anyone has checked.
    """

    gene: str
    kind: str  # "protein" (codon numbers) or "nucleotide" (c./n. positions)
    offset: int
    legacy_scheme: str
    first: int
    last: int
    source: str
    note: str = ""
    verified: bool = True

    def to_legacy(self, number: int) -> Optional[int]:
        if self.first <= number <= self.last:
            return number + self.offset
        return None

    def to_canonical(self, number: int) -> Optional[int]:
        candidate = number - self.offset
        if self.first <= candidate <= self.last:
            return candidate
        return None

    def describe(self) -> str:
        mark = "" if self.verified else " [correspondence unverified]"
        return (
            "{0} {1} numbering: {2} = M. tuberculosis number {3:+d} over {4}-{5} "
            "({6}){7}".format(self.gene, self.kind, self.legacy_scheme, self.offset,
                              self.first, self.last, self.source, mark)
        )


#: SOURCE for the rpoB offset: the standard M. tuberculosis / E. coli RpoB
#: correspondence used throughout the TB literature and in WHO-UCN-TB-2023.7,
#: which grades rpoB in M. tuberculosis numbering while most pre-2015 papers and
#: several commercial assays report E. coli numbering. Four independent
#: checkpoints agree on +81 across the RRDR: Ser450Leu = Ser531Leu,
#: His445Asp = His526Asp, Asp435Val = Asp516Val, Leu430Pro = Leu511Pro. The
#: range is deliberately confined to the RRDR and its margins, because the
#: offset is an alignment artefact and is not constant over the whole protein.
_RPOB_ECOLI = NumberingAlias(
    gene="rpoB", kind="protein", offset=81, legacy_scheme="E. coli RpoB",
    first=400, last=510, source=SRC_WHO_V2,
    note="verified at codons 430, 435, 445 and 450; outside the RRDR margins the "
         "offset is not established and no alias is produced")

#: SOURCE: Prammananan et al. 1998 and the E. coli 16S numbering convention used
#: throughout the aminoglycoside-resistance literature. MTBseq's own file
#: carries both schemes side by side ("Codon nr." and "Codon nr. E. coli"), and
#: its rows give the two anchors this offset is taken from: mycobacterial rrs
#: 1401 = E. coli 1408 and rrs 1484 = E. coli 1491. Both give +7. The range is
#: confined to the 3' end of the molecule where those anchors sit; the offset is
#: not assumed to hold across the whole 16S.
_RRS_ECOLI = NumberingAlias(
    gene="rrs", kind="nucleotide", offset=7, legacy_scheme="E. coli 16S",
    first=1350, last=1537, source=SRC_PRAMMANANAN_1998,
    note="anchored on rrs 1401 = E. coli 1408 (the amikacin position) and rrs "
         "1484 = E. coli 1491, both present in MTBseq's dual-numbered table",
    verified=True)

#: No alias is defined for gyrA. The gyrA numbering discrepancy is real and is
#: discussed in the literature, but no source on this machine establishes the
#: offset, and an invented one would silently merge two different codons into a
#: single call. An unmatched key is recoverable; a wrong match is not.
NUMBERING_ALIASES: Dict[str, Tuple[NumberingAlias, ...]] = {
    "rpob": (_RPOB_ECOLI,),
    "rrs": (_RRS_ECOLI,),
}


def numbering_aliases(gene: str) -> Tuple[NumberingAlias, ...]:
    """Documented legacy numbering schemes for a gene, empty when none."""
    return NUMBERING_ALIASES.get(str(gene or "").strip().lower(), ())


def _renumber(hgvs: str, alias: NumberingAlias, forward: bool) -> str:
    """Rewrite the one number in an HGVS string under a numbering alias."""
    if alias.kind == "protein":
        match = _PROTEIN_SPLIT.match(hgvs)
        if not match:
            return ""
        number = int(match.group("num"))
        moved = alias.to_legacy(number) if forward else alias.to_canonical(number)
        if moved is None:
            return ""
        return "{0}{1}{2}{3}".format(match.group("prefix"), match.group("ref"),
                                     moved, match.group("rest"))
    match = _NUCLEOTIDE_SPLIT.match(hgvs)
    if not match:
        return ""
    number = int(match.group("num"))
    if number < 0:
        # Promoter positions count back from the start codon and are not part of
        # the mature-transcript numbering the alias describes.
        return ""
    moved = alias.to_legacy(number) if forward else alias.to_canonical(number)
    if moved is None:
        return ""
    return "{0}{1}{2}".format(match.group("prefix"), moved, match.group("rest"))


def alias_hgvs(gene: str, hgvs: str) -> Tuple[str, ...]:
    """Every other spelling of this variant name Mjolnir recognises.

    Two kinds of alias, and they are different in kind. The one-letter form is a
    pure notation difference with no possibility of error. The legacy numbering
    forms are a claim about two organisms' sequences, so each one is bounded by
    the range its source establishes and produces nothing outside it.
    """
    canonical = normalise_hgvs(hgvs)
    if not canonical:
        return ()
    out: List[str] = []

    short = one_letter(canonical)
    if short and short != canonical:
        out.append(short)

    for alias in numbering_aliases(gene):
        for forward in (True, False):
            moved = _renumber(canonical, alias, forward)
            if moved and moved != canonical and moved not in out:
                out.append(moved)
    return tuple(out)


def alias_keys(key: str) -> Tuple[str, ...]:
    """Alternative ``<gene>_<hgvs>`` keys for a join key."""
    gene, hgvs = split_key(key)
    if not gene or not hgvs:
        return ()
    return tuple(hgvs_key(gene, alt) for alt in alias_hgvs(gene, hgvs) if alt)


def classify_difference(key_a: str, key_b: str) -> str:
    """Why two catalogues named different variants for the same drug.

    Returns one of the :mod:`records` disagreement kinds. ``nomenclature`` is
    the whole point of this module: ``rpoB_p.Ser450Leu`` against
    ``rpoB_p.Ser531Leu`` is one mutation written twice, and reporting it as a
    biological disagreement would tell a clinician the evidence conflicts when
    it agrees exactly.
    """
    left = normalise_hgvs_key(key_a)
    right = normalise_hgvs_key(key_b)
    if not left or not right:
        return DISAGREEMENT_NONE if left == right else DISAGREEMENT_BIOLOGICAL
    if left == right:
        return DISAGREEMENT_NONE
    if right in alias_keys(left) or left in alias_keys(right):
        return DISAGREEMENT_NOMENCLATURE
    return DISAGREEMENT_BIOLOGICAL


def normalise_hgvs_key(key: str, default_prefix: str = "c.") -> str:
    """Canonicalise a whole ``<gene>_<hgvs>`` join key."""
    gene, hgvs = split_key(key)
    if not gene:
        return normalise_hgvs(hgvs, default_prefix)
    return hgvs_key(gene, normalise_hgvs(hgvs, default_prefix))


# ---------------------------------------------------------------------------
# The two keys together
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NormalisedVariant:
    """One observed variant with both of the identities the design requires.

    ``left_aligned`` False is the honest outcome when no reference sequence was
    available to shift an indel with. The record keeps the flag rather than
    dropping it, so a catalogue lookup that missed can be explained instead of
    appearing as an absence of any catalogued mutation — which is the one
    sentence this project is most careful about.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = ""
    hgvs: str = ""
    #: The spelling the caller supplied, when normalisation changed it.
    original_hgvs: str = ""
    left_aligned: bool = True
    variant_type: str = VARIANT_SNP
    aliases: Tuple[str, ...] = ()

    @property
    def coordinate_key(self) -> CoordinateKey:
        return (self.chrom, self.pos, self.ref, self.alt)

    @property
    def hgvs_key(self) -> str:
        return hgvs_key(self.gene, self.hgvs)

    @property
    def keys(self) -> Tuple[str, ...]:
        """The HGVS key and all of its aliases, for a cross-catalogue lookup."""
        primary = self.hgvs_key
        if not primary:
            return ()
        return (primary,) + tuple(k for k in self.aliases if k != primary)

    @property
    def note(self) -> str:
        if self.left_aligned:
            return ""
        return (
            "indel could not be left-aligned: no reference sequence was supplied, "
            "so a catalogue lookup on its coordinate may miss an entry written at "
            "the leftmost equivalent position"
        )


def normalise(chrom: str, pos: int, ref: str, alt: str, gene: str = "",
              hgvs: str = "", reference: Reference = None,
              default_prefix: str = "c.") -> NormalisedVariant:
    """Both keys for one observed variant."""
    new_pos, new_ref, new_alt, aligned = left_align(
        chrom or H37RV_ACCESSION, pos, ref, alt, reference)
    canonical = normalise_hgvs(hgvs, default_prefix)
    original = str(hgvs or "").strip()
    return NormalisedVariant(
        chrom=chrom or H37RV_ACCESSION,
        pos=new_pos,
        ref=new_ref,
        alt=new_alt,
        gene=str(gene or "").strip(),
        hgvs=canonical,
        original_hgvs="" if canonical == original else original,
        left_aligned=aligned,
        variant_type=variant_class(new_ref, new_alt),
        aliases=alias_keys(hgvs_key(gene, canonical)),
    )


def normalise_variant(variant: Variant, reference: Reference = None,
                      default_prefix: str = "c.") -> Variant:
    """A copy of a :class:`records.Variant` carrying normalised keys.

    A copy rather than a mutation, because the caller's record may already be
    inside a :class:`records.SampleResult` and a variant whose coordinates
    changed under it is a debugging problem nobody enjoys. The original HGVS
    spelling is preserved in ``hgvs_alias`` so the annex can show what the
    caller wrote as well as what Mjolnir matched on.

    ``variant_type`` is recomputed from the alleles, but a caller who declared
    ``lof`` keeps it: loss of function is a rule-derived class (§5.4), not
    something the allele lengths can tell you.
    """
    result = normalise(variant.chrom, variant.pos, variant.ref, variant.alt,
                       gene=variant.gene, hgvs=variant.hgvs, reference=reference,
                       default_prefix=default_prefix)
    filters = list(variant.filters)
    if not result.left_aligned and "not-left-aligned" not in filters:
        filters.append("not-left-aligned")
    keep_type = variant.variant_type if variant.variant_type in (
        VARIANT_INDEL, VARIANT_INS, VARIANT_DEL) and result.variant_type in (
        VARIANT_INDEL, VARIANT_INS, VARIANT_DEL) else result.variant_type
    if variant.variant_type == "lof":
        keep_type = variant.variant_type
    note = variant.note
    if result.note and result.note not in note:
        note = "; ".join(part for part in (note, result.note) if part)
    return replace(
        variant,
        pos=result.pos,
        ref=result.ref,
        alt=result.alt,
        hgvs=result.hgvs,
        hgvs_alias=variant.hgvs_alias or result.original_hgvs,
        variant_type=keep_type,
        filters=filters,
        note=note,
    )


def normalise_variants(variants: Sequence[Variant], reference: Reference = None,
                       default_prefix: str = "c.") -> List[Variant]:
    """:func:`normalise_variant` over a list, keeping order."""
    out: List[Variant] = []
    unaligned = 0
    for variant in variants:
        normalised = normalise_variant(variant, reference, default_prefix)
        if "not-left-aligned" in normalised.filters:
            unaligned += 1
        out.append(normalised)
    if unaligned:
        LOG.warning(
            "%d indel(s) could not be left-aligned because no reference sequence "
            "was supplied; their coordinate lookups may miss catalogue entries",
            unaligned)
    return out


__all__ = [
    "CoordinateKey", "NormalisedVariant", "NumberingAlias", "NUMBERING_ALIASES",
    "RULE_DELETION", "RULE_LOF",
    "alias_hgvs", "alias_keys", "classify_difference", "coordinate_string",
    "hgvs_key", "is_rule_variant", "left_align", "normalise",
    "normalise_coordinate", "normalise_hgvs", "normalise_hgvs_key",
    "normalise_variant", "normalise_variants", "numbering_aliases", "one_letter",
    "parse_coordinate_string", "split_key", "three_letter", "trim_alleles",
    "variant_class",
]
