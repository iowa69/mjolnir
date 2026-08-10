"""The grading rules that live outside the catalogue table.

The WHO catalogue is three components, not one spreadsheet: a table of graded
variant-drug rows, a set of *additional grading rules* that apply to variants the
table does not list, and a Comment column that carries clinical meaning the grade
does not. A tool that implements only the lookup is wrong in a way that is
invisible from its output — it will report "no resistance determinant detected"
for a frameshift in *pncA* and for a novel indel inside the rifampicin
resistance-determining region, which are two of the least ambiguous resistance
findings in tuberculosis.

So this module implements the rules (design §5.4):

* a novel **silent** variant grades Group 4 — not "no call", and not Group 3
  either, because a synonymous change is positively expected not to matter;
* any non-synonymous change or indel inside the *rpoB* RRDR, codons 426-452,
  grades Group 2 for rifampicin;
* loss of function grades Group 2 for the drugs whose mechanism the gene carries
  — *katG*/isoniazid, *pncA*/pyrazinamide, *Rv0678* and *pepQ*/bedaquiline and
  clofazimine, the nitroreductase set/delamanid and pretomanid, *ethA*/ethionamide;
* the four borderline *rpoB* mutations are Group 1 by rule, because they sit at
  the edge of phenotypic detectability and are the classic source of a
  "phenotypically susceptible, genotypically resistant" argument at the bench.

Two design decisions here are worth stating because they are the ones a reviewer
will want to argue with.

**A rule never overrides a graded row.** ``rule_hits`` takes the drugs WHO has
already graded *this* variant for and stays silent on them. These are *additional*
grading rules: they exist because the table is finite, not because it is wrong.
Where the table speaks, the table is the answer, and a rule that quietly promoted
a Group 3 row to Group 2 would put WHO's authority behind a number WHO did not
publish.

**Epistasis suppresses, it does not delete.** :func:`epistasis_suppressions`
returns records rather than editing anything. *mmpL5* loss of function abrogates
the *Rv0678* efflux phenotype, and a coding loss of function in *eis* abrogates
*eis* promoter mutations — but a report that simply stopped mentioning the
*Rv0678* variant would be indistinguishable from a report that never found it.
The suppression is a first-class object with the rule, the abrogating variant,
the abrogated variants and the reason, and ``consensus.py`` is required to print
it.

One honest limit is encoded rather than assumed. Epistasis is a statement about
one genome. If the abrogating loss of function is a major variant and the
abrogated mutation is a minor one, or the reverse, the two may sit in different
subpopulations and the suppression may not hold — so the record is emitted with
``confident=False`` and the consensus engine reports it without applying it.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ..config import (
    CATALOGUE_WHO,
    COMMENT_HIGH_LEVEL,
    CROSS_RESISTANCE_PAIRS,
    DRUG_CODES,
    EPISTASIS_RULES,
    LOF_EFFECTS,
    LOF_GROUP2_GENES,
    RPOB_BORDERLINE,
    RPOB_RRDR_CODONS,
    SILENT_VARIANT_GRADE,
    SRC_WHO_V2,
    WHO_GRADE_1,
    WHO_GRADE_2,
    call_for_grade,
    is_major_variant,
    normalise_drug,
    source_for,
)
from ..records import CatalogueCall, Variant
from ..utils import AA_1_TO_3, AA_3_TO_1, natural_key, to_jsonable

# ---------------------------------------------------------------------------
# Rule identities. Spelled once so the report, the tests and the CatalogueCall
# evidence field cannot drift apart.
# ---------------------------------------------------------------------------

RULE_SILENT = "who-additional:novel-silent-group4"
RULE_RPOB_RRDR = "who-additional:rpoB-RRDR-group2"
RULE_RPOB_BORDERLINE = "who-additional:rpoB-borderline-group1"
RULE_LOF = "who-additional:loss-of-function-group2"

#: All rule identities, in the order they are evaluated.
RULES: Tuple[str, ...] = (RULE_RPOB_BORDERLINE, RULE_RPOB_RRDR, RULE_LOF, RULE_SILENT)

#: The drug the RRDR rule is about. Named rather than inlined so that a reader
#: who does not know that RRDR means rifampicin does not have to guess.
RRDR_DRUG = normalise_drug("Rifampicin")

#: SOURCE: WHO-UCN-TB-2023.7, Comment column. ``config.COMMENT_HIGH_LEVEL`` fixes
#: the high-level phrase; the low-level and additive phrasings below are the
#: counterparts that appear in the same column. They are spelled here because
#: config.py does not carry them, and they are matched case-insensitively as
#: substrings so that a trailing qualifier in the cell does not defeat the match.
#:
#: UNVERIFIED: the xlsx is not on this machine (no downloads permitted during
#: this build), so these two strings are from the design and the published
#: catalogue text rather than from a cell that was opened and read. The
#: high-level phrase, which comes from config, is the one the report leans on.
COMMENT_LOW_LEVEL = "Low-level resistance"
COMMENT_ADDITIVE = "additive"

LEVEL_HIGH = "high-level"
LEVEL_LOW = "low-level"
LEVEL_ADDITIVE = "additive low-level"

#: Cues that a Comment is talking about cross-resistance rather than merely
#: naming another drug in passing. Requiring a cue is what stops "not associated
#: with resistance to bedaquiline" from being read as clofazimine cross-resistance.
CROSS_RESISTANCE_CUES: Tuple[str, ...] = (
    "cross-resistance", "cross resistance", "crossresistance",
    "cross-resistant", "cross resistant",
)

#: Effects that mean the variant sits upstream of the coding sequence. The
#: *eis* epistasis rule is specifically about promoter mutations, which act by
#: over-expressing the protein, so telling promoter from coding is load-bearing
#: rather than cosmetic here.
PROMOTER_EFFECTS: Tuple[str, ...] = (
    "upstream", "promoter", "5_prime_utr", "five_prime_utr", "regulatory",
)

#: Pooled graded-variant names WHO uses where no coordinate exists: these are
#: matched by rule, not by coordinate (design §5.2, last bullet).
POOLED_LOF_NAMES: Tuple[str, ...] = (
    "lof", "loss_of_function", "loss of function", "deletion", "any_lof",
    "any_indel", "feature_ablation",
)


def _canonical_gene(gene: str) -> str:
    return " ".join(str(gene or "").split()).lower()


def _gene_is(gene: str, wanted: str) -> bool:
    """Gene identity, case-insensitively. *Rv0678* and *rv0678* are one gene."""
    return _canonical_gene(gene) == _canonical_gene(wanted)


#: gene -> the drugs whose loss-of-function rule that gene triggers. Inverted
#: from ``config.LOF_GROUP2_GENES``, which is keyed by drug because that is how
#: the catalogue documents it, while every caller here has a gene in hand.
LOF_GENE_DRUGS: Dict[str, Tuple[str, ...]] = {}
for _drug, _genes in LOF_GROUP2_GENES.items():
    for _gene in _genes:
        _key = _canonical_gene(_gene)
        LOF_GENE_DRUGS[_key] = tuple(
            sorted(set(LOF_GENE_DRUGS.get(_key, ()) + (normalise_drug(_drug),)))
        )


def lof_drugs_for_gene(gene: str) -> Tuple[str, ...]:
    """Drugs for which loss of function in *gene* grades Group 2, in name order."""
    return LOF_GENE_DRUGS.get(_canonical_gene(gene), ())


# ---------------------------------------------------------------------------
# HGVS
#
# Every rule below is a question about a codon, a consequence or a region, and
# the only thing that reliably carries all three across the three catalogues is
# the HGVS name. So the parsing lives here, once, rather than as a regex
# scattered through each rule.
# ---------------------------------------------------------------------------

_PROTEIN_HGVS = re.compile(
    r"^p\.\(?"
    r"(?P<ref>Ter|[A-Z][a-z]{2}|[A-Z*])"
    r"(?P<codon>\d+)"
    r"(?P<rest>[^)]*)"
    r"\)?$"
)
_BARE_PROTEIN = re.compile(r"^(?P<ref>[A-Z*])(?P<codon>\d+)(?P<rest>[A-Z*=]?)$")
_NUCLEOTIDE_HGVS = re.compile(r"^(?P<prefix>[cngr])\.(?P<pos>[-+*]?\d+)")
_LEADING_AA = re.compile(r"^(Ter|[A-Z][a-z]{2}|[A-Z*])")


@dataclass(frozen=True)
class ProteinChange:
    """One amino-acid substitution as three comparable pieces.

    ``ref`` and ``alt`` are one-letter codes because comparing three-letter and
    one-letter spellings is the single most common way two catalogues appear to
    disagree when they do not. ``rest`` keeps whatever followed the codon —
    ``fs``, ``del``, ``=``, ``?`` — because those suffixes are the evidence for
    frameshift, deletion, synonymy and start loss respectively.
    """

    ref: str
    codon: int
    alt: str
    rest: str = ""

    @property
    def is_synonymous(self) -> bool:
        if self.rest.startswith("="):
            return True
        return bool(self.alt) and self.alt == self.ref

    @property
    def is_frameshift(self) -> bool:
        return "fs" in self.rest.lower()

    @property
    def is_stop_gained(self) -> bool:
        return self.alt == "*" and self.ref != "*"

    @property
    def is_start_lost(self) -> bool:
        if self.codon != 1:
            return False
        if self.rest.startswith("?"):
            return True
        return bool(self.alt) and self.ref == "M" and self.alt != "M"

    def three_letter(self) -> str:
        """``Leu430Pro`` — the spelling ``config.RPOB_BORDERLINE`` uses."""
        ref = AA_1_TO_3.get(self.ref, self.ref)
        alt = AA_1_TO_3.get(self.alt, self.alt)
        return "{0}{1}{2}".format(ref, self.codon, alt)


def _one_letter(token: str) -> str:
    """One-letter amino acid for a 1- or 3-letter token, or "" if it is neither."""
    if not token:
        return ""
    if token in ("*", "Ter", "ter", "TER"):
        return "*"
    if len(token) == 1 and token.upper() in AA_1_TO_3:
        return token.upper()
    return AA_3_TO_1.get(token.capitalize(), "")


def parse_protein_change(hgvs: str) -> Optional[ProteinChange]:
    """Parse a ``p.`` HGVS name, or return None when it is not one.

    Accepts the bare legacy form too (``S450L``): MTBseq's resistance list uses
    it, and refusing to read it would turn a matching problem into a silent
    absence of evidence.
    """
    text = " ".join(str(hgvs or "").split())
    if not text:
        return None
    match = _PROTEIN_HGVS.match(text) or _BARE_PROTEIN.match(text)
    if match is None:
        return None
    ref = _one_letter(match.group("ref"))
    if not ref:
        return None
    rest = match.group("rest") or ""
    alt = ""
    lead = _LEADING_AA.match(rest)
    if lead is not None:
        alt = _one_letter(lead.group(0))
        if alt:
            rest = rest[lead.end():]
    return ProteinChange(ref=ref, codon=int(match.group("codon")), alt=alt, rest=rest)


def nucleotide_position(hgvs: str) -> Optional[int]:
    """The first coordinate in a ``c.``/``n.``/``g.`` HGVS name, signed.

    Negative means upstream of the start codon, which is how the catalogues
    spell a promoter mutation — ``eis_c.-14C>T`` and ``inhA_c.-154G>A`` are both
    promoter variants and both must be recognised as such.
    """
    text = " ".join(str(hgvs or "").split())
    match = _NUCLEOTIDE_HGVS.match(text)
    if match is None:
        return None
    raw = match.group("pos")
    if raw.startswith("*"):
        # 3' UTR coordinate; positive but not a coding position.
        return None
    return int(raw.replace("+", ""))


def codon_of(variant: Variant) -> Optional[int]:
    """The codon a variant falls in, or None when it cannot be established.

    A protein name answers directly. A coding nucleotide name answers by
    division. A promoter or non-coding name has no codon at all, and returns
    None rather than 0 — codon 0 does not exist, and a rule that treated it as
    one would place every promoter variant at the start of the gene.
    """
    change = parse_protein_change(variant.hgvs)
    if change is not None:
        return change.codon
    pos = nucleotide_position(variant.hgvs)
    if pos is None or pos < 1:
        return None
    return (pos + 2) // 3


def _normalised_effect(effect: str) -> str:
    text = " ".join(str(effect or "").split()).lower().replace(" ", "_")
    if text.endswith("_variant"):
        text = text[: -len("_variant")]
    return text


def is_promoter_variant(variant: Variant) -> bool:
    """Whether the variant sits upstream of the coding sequence."""
    effect = _normalised_effect(variant.effect)
    for token in PROMOTER_EFFECTS:
        if token in effect:
            return True
    pos = nucleotide_position(variant.hgvs)
    return pos is not None and pos < 0


def is_coding_variant(variant: Variant) -> bool:
    """Whether the variant changes the coding sequence rather than its promoter."""
    if is_promoter_variant(variant):
        return False
    if parse_protein_change(variant.hgvs) is not None:
        return True
    if _canonical_gene(variant.hgvs) in POOLED_LOF_NAMES:
        return True
    pos = nucleotide_position(variant.hgvs)
    if pos is not None and pos >= 1:
        return True
    effect = _normalised_effect(variant.effect)
    return any(token in effect for token in
               ("missense", "synonymous", "frameshift", "stop", "start",
                "inframe", "coding"))


def is_synonymous(variant: Variant) -> bool:
    """Whether the variant leaves the protein sequence unchanged.

    Silent is a positive finding, not a missing annotation: it is the trigger
    for the Group 4 rule. So this answers True only on evidence — an effect
    string that says synonymous, or a protein name whose reference and
    alternative residues are the same — and never by default.
    """
    effect = _normalised_effect(variant.effect)
    if effect in ("synonymous", "silent", "synonymous_codon"):
        return True
    if effect.startswith("synonymous"):
        return True
    change = parse_protein_change(variant.hgvs)
    if change is not None and change.is_synonymous:
        return True
    return False


def is_loss_of_function(variant: Variant) -> bool:
    """Whether the variant destroys the gene product.

    Deliberately conservative about indels. A frameshift is loss of function; an
    in-frame deletion of one codon usually is not, and treating every indel as
    LoF would hand out Group 2 calls for *pncA* on the strength of a three-base
    deletion nobody characterised. So this needs an effect annotation, an HGVS
    suffix, or one of WHO's pooled ``<gene>_LoF`` / ``<gene>_deletion`` names —
    never the length difference alone.
    """
    if variant.variant_type == "lof":
        return True
    effect = _normalised_effect(variant.effect)
    if effect:
        for known in LOF_EFFECTS:
            token = _normalised_effect(known)
            if effect == token or effect.startswith(token) or token in effect:
                return True
    hgvs = " ".join(str(variant.hgvs or "").split())
    lowered = hgvs.lower()
    if _canonical_gene(hgvs) in POOLED_LOF_NAMES:
        return True
    if lowered.endswith("_lof") or lowered.endswith("_deletion"):
        return True
    change = parse_protein_change(hgvs)
    if change is not None:
        if change.is_frameshift or change.is_stop_gained or change.is_start_lost:
            return True
    if "fs" in lowered and "p." in lowered:
        return True
    return False


def is_in_rpob_rrdr(variant: Variant) -> bool:
    """Whether the variant falls inside the rifampicin resistance-determining region.

    SOURCE for the codon bounds: ``config.RPOB_RRDR_CODONS``, from
    WHO-UCN-TB-2023.7. The range is inclusive at both ends.
    """
    if not _gene_is(variant.gene, "rpoB"):
        return False
    codon = codon_of(variant)
    if codon is None:
        return False
    first, last = RPOB_RRDR_CODONS
    return first <= codon <= last


def rpob_borderline_hit(variant: Variant) -> str:
    """The borderline *rpoB* mutation this variant is, or "".

    SOURCE: ``config.RPOB_BORDERLINE`` — Leu430Pro, His445Asn, His445Ser and
    Ile491Phe, Group 1 by rule. Ile491Phe sits outside the RRDR, which is exactly
    why it needs a rule of its own: an RRDR-bounded implementation misses it.
    """
    if not _gene_is(variant.gene, "rpoB"):
        return ""
    change = parse_protein_change(variant.hgvs)
    if change is None:
        return ""
    spelled = change.three_letter().lower()
    for known in RPOB_BORDERLINE:
        if spelled == known.lower():
            return known
    return ""


# ---------------------------------------------------------------------------
# Rule hits
# ---------------------------------------------------------------------------

@dataclass
class RuleHit:
    """One grading rule firing on one variant for one drug.

    Carries the same provenance a table row does — the grade, the resulting
    call, the source document and a sentence saying why — because a rule-derived
    Group 2 and a table-derived Group 2 end up in the same column of the same
    report, and the reader is entitled to know which is which.
    """

    rule: str
    drug: str
    gene: str
    grade: str
    call: str
    why: str
    source: str
    variant_key: str = ""
    coordinate_key: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable(asdict(self))

    def to_catalogue_call(self, catalogue_version: str = "",
                          catalogue_checksum: str = "") -> CatalogueCall:
        """As a :class:`CatalogueCall` attributed to the WHO catalogue.

        Attributed to WHO, not to a fourth pseudo-catalogue, because these are
        WHO's own additional grading rules. ``matched_by="rule"`` is how the
        annex distinguishes them from a coordinate match against the table.
        """
        return CatalogueCall(
            catalogue=CATALOGUE_WHO,
            drug=self.drug,
            grade=self.grade,
            comment=self.why,
            source=self.source,
            call=self.call,
            variant_key=self.variant_key,
            catalogue_version=catalogue_version,
            catalogue_checksum=catalogue_checksum,
            matched_by="rule",
            evidence=self.rule,
        )


def _keys(variant: Variant) -> Tuple[str, str]:
    coordinate = "{0}:{1}{2}>{3}".format(
        variant.chrom, variant.pos, (variant.ref or "").upper(), (variant.alt or "").upper())
    return variant.hgvs_key or variant.display, coordinate


def rule_hits(variant: Variant, *,
              gene_drugs: Sequence[str] = (),
              who_graded_drugs: Sequence[str] = ()) -> List[RuleHit]:
    """Every additional grading rule that fires on one variant.

    *who_graded_drugs* is the set of drugs for which the WHO table already grades
    this exact variant. Rules stay silent on those: they are additional grading
    rules, applied where the table is finite, and a rule that overrode a
    published row would put WHO's authority behind a grade WHO did not give.

    *gene_drugs* is the set of drugs the variant's gene is catalogued for, and it
    is used only by the silent rule. It has to be supplied by the caller —
    ``catalogues.py`` derives it from the xlsx it actually loaded — because
    Mjolnir will not guess which drug a gene belongs to. With no gene-drug
    association supplied, a silent variant produces no rule hit, which is the
    correct answer for a synonymous change in a gene nobody has associated with a
    drug.
    """
    hits: List[RuleHit] = []
    gene = str(variant.gene or "")
    already = set(normalise_drug(d) for d in who_graded_drugs)
    variant_key, coordinate_key = _keys(variant)

    def _add(rule: str, drug: str, grade: str, why: str, source: str) -> None:
        canonical = normalise_drug(drug)
        if canonical in already:
            return
        hits.append(RuleHit(
            rule=rule, drug=canonical, gene=gene, grade=grade,
            call=call_for_grade(grade), why=why, source=source,
            variant_key=variant_key, coordinate_key=coordinate_key))

    borderline = rpob_borderline_hit(variant)
    if borderline:
        _add(RULE_RPOB_BORDERLINE, RRDR_DRUG, WHO_GRADE_1,
             "rpoB {0} is one of the four borderline mutations graded Group 1 by "
             "rule rather than by table lookup; they sit at the edge of "
             "phenotypic detectability and are the usual source of a "
             "phenotypically-susceptible, genotypically-resistant argument"
             .format(borderline),
             source_for("rpob_borderline"))
    elif is_in_rpob_rrdr(variant) and not is_synonymous(variant):
        codon = codon_of(variant)
        _add(RULE_RPOB_RRDR, RRDR_DRUG, WHO_GRADE_2,
             "codon {0} lies inside the rpoB rifampicin resistance-determining "
             "region (codons {1}-{2}); any non-synonymous change or indel there "
             "grades Group 2 by rule".format(codon, *RPOB_RRDR_CODONS),
             source_for("rpob_rrdr_codons"))

    if is_loss_of_function(variant):
        for drug in lof_drugs_for_gene(gene):
            _add(RULE_LOF, drug, WHO_GRADE_2,
                 "loss of function in {0} grades Group 2 for {1} by rule: the "
                 "gene product is required for the drug's activity or for "
                 "keeping its efflux route closed".format(gene, drug),
                 source_for("lof_group2_genes"))

    if is_synonymous(variant):
        for drug in sorted(set(normalise_drug(d) for d in gene_drugs if d),
                           key=natural_key):
            _add(RULE_SILENT, drug, SILENT_VARIANT_GRADE,
                 "a novel synonymous variant grades Group 4 by rule: it leaves "
                 "the protein sequence unchanged, so the catalogue's default "
                 "expectation is that it is not associated with resistance",
                 source_for("silent_variant_grade"))

    return hits


def rule_catalogue_calls(variant: Variant, *,
                         gene_drugs: Sequence[str] = (),
                         who_graded_drugs: Sequence[str] = (),
                         catalogue_version: str = "",
                         catalogue_checksum: str = "") -> List[CatalogueCall]:
    """:func:`rule_hits`, already converted for the consensus engine."""
    return [hit.to_catalogue_call(catalogue_version, catalogue_checksum)
            for hit in rule_hits(variant, gene_drugs=gene_drugs,
                                 who_graded_drugs=who_graded_drugs)]


def annotate_variant(variant: Variant, *,
                     gene_drugs: Sequence[str] = (),
                     catalogue_version: str = "",
                     catalogue_checksum: str = "") -> List[CatalogueCall]:
    """Append the rule-derived calls to *variant* and return the new ones.

    The drugs WHO has already graded are read off the variant's existing WHO
    calls, so the caller cannot forget to pass them and accidentally let a rule
    shadow a published row. Mutates ``variant.catalogue_calls`` — the one place
    in this module that changes anything — because the pipeline wants the rule
    grades to travel with the variant into the annex table.
    """
    graded = [c.drug for c in variant.catalogue_calls
              if c.catalogue == CATALOGUE_WHO and c.matched_by != "rule"]
    added = rule_catalogue_calls(
        variant, gene_drugs=gene_drugs, who_graded_drugs=graded,
        catalogue_version=catalogue_version, catalogue_checksum=catalogue_checksum)
    variant.catalogue_calls.extend(added)
    return added


# ---------------------------------------------------------------------------
# The Comment column
#
# Level of resistance and cross-resistance are read from the Comment column, not
# from the grade (design §5.4). Group 1 says a mutation is associated with
# resistance; it does not say whether the resistance is high-level, which is the
# part that changes whether a drug can still be used at an increased dose.
# ---------------------------------------------------------------------------

def level_from_comment(comment: str) -> str:
    """Level of resistance a Comment states, or "" when it states none.

    Returns "" rather than "low-level" for a silent comment. A comment that does
    not name a level has not measured one, and defaulting to the milder reading
    would be exactly the "unmeasured described as fine" failure the house rules
    forbid.
    """
    text = " ".join(str(comment or "").split()).lower()
    if not text:
        return ""
    if COMMENT_ADDITIVE.lower() in text:
        return LEVEL_ADDITIVE
    if COMMENT_HIGH_LEVEL.lower() in text:
        return LEVEL_HIGH
    if COMMENT_LOW_LEVEL.lower() in text:
        return LEVEL_LOW
    return ""


#: Whole-word patterns for every drug name and code Mjolnir knows, so that a
#: Comment mentioning "CFZ" or "clofazimine" is found and one mentioning
#: "clofazimine-sparing" is still found while "PA" inside a word is not.
_DRUG_PATTERNS: List[Tuple[str, Any]] = []
for _name, _code in sorted(DRUG_CODES.items()):
    _alternatives = "|".join(re.escape(t) for t in (_name, _code))
    _DRUG_PATTERNS.append((_name, re.compile(r"\b(?:{0})\b".format(_alternatives),
                                             re.IGNORECASE)))


def drugs_named_in(comment: str) -> List[str]:
    """Canonical names of every drug a Comment mentions, in catalogue order."""
    text = str(comment or "")
    if not text:
        return []
    return [name for name, pattern in _DRUG_PATTERNS if pattern.search(text)]


def documented_cross_resistance(drug: str) -> Tuple[str, ...]:
    """Partners of *drug* in ``config.CROSS_RESISTANCE_PAIRS``."""
    canonical = normalise_drug(drug)
    partners: List[str] = []
    for a, b in CROSS_RESISTANCE_PAIRS:
        if normalise_drug(a) == canonical:
            partners.append(normalise_drug(b))
        elif normalise_drug(b) == canonical:
            partners.append(normalise_drug(a))
    return tuple(sorted(set(partners)))


def cross_resistance_from_comment(comment: str, drug: str) -> List[str]:
    """Drugs a Comment says this variant also confers resistance to.

    Requires a cross-resistance cue as well as a drug name. Comments name other
    drugs for all sorts of reasons — a variant graded differently for a second
    drug, a note about which agent the phenotype was measured against — and
    reading every co-occurrence as cross-resistance would manufacture findings.

    When the cue is present but no partner drug is named, the documented pairs in
    ``config.CROSS_RESISTANCE_PAIRS`` supply the partner, since those are the
    pairs WHO documents: delamanid with pretomanid, bedaquiline with clofazimine.
    """
    text = " ".join(str(comment or "").split())
    if not text:
        return []
    lowered = text.lower()
    if not any(cue in lowered for cue in CROSS_RESISTANCE_CUES):
        return []
    canonical = normalise_drug(drug)
    named = [d for d in drugs_named_in(text) if d != canonical]
    if not named:
        named = list(documented_cross_resistance(canonical))
    return sorted(set(named), key=natural_key)


# ---------------------------------------------------------------------------
# Epistasis
# ---------------------------------------------------------------------------

@dataclass
class Suppression:
    """One epistasis rule, the variants it involves, and whether it was applied.

    A record rather than a mutation. The design requires the suppression to be
    *stated* in the report — "an Rv0678 mutation is present but mmpL5 is
    knocked out, so no bedaquiline resistance is predicted" is a different and
    more useful sentence than either "resistant" or "no determinant detected",
    and it can only be written if the abrogated variant survives into the output.
    """

    rule: str
    drug: str
    suppressor_gene: str
    suppressor_variants: List[str] = field(default_factory=list)
    suppressed_gene: str = ""
    suppressed_variants: List[str] = field(default_factory=list)
    why: str = ""
    source: str = ""
    #: False when the two variants may sit in different subpopulations, in which
    #: case the consensus engine reports the rule without applying it.
    confident: bool = True
    caveat: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return to_jsonable(asdict(self))

    def suppresses(self, variant_key: str) -> bool:
        return variant_key in self.suppressed_variants


#: SOURCE: Mjolnir policy, following design §5.4. Epistasis is a statement about
#: one genome: an abrogating loss of function that is present in only part of the
#: population cannot be assumed to abrogate a mutation carried by the rest. The
#: major/minor split used to detect that is ``config.MAJOR_VARIANT_FRACTION``, so
#: no new number is introduced here.
MIXED_SUBPOPULATION_CAVEAT = (
    "the abrogating loss of function and the mutation it would abrogate are not "
    "both major variants, so they may sit in different subpopulations; the "
    "epistasis rule is reported but not applied"
)


def _rule_id(rule: Dict[str, Any]) -> str:
    region = rule.get("suppressed_region")
    target = "{0}{1}".format(rule.get("suppressed_gene", ""),
                             " " + str(region) if region else "")
    return "epistasis:{0}-{1}-abrogates-{2}".format(
        rule.get("suppressor_gene", ""), rule.get("suppressor_effect", "lof"),
        target.replace(" ", "-"))


def _matches_suppressor(variant: Variant, rule: Dict[str, Any]) -> bool:
    if not _gene_is(variant.gene, rule.get("suppressor_gene", "")):
        return False
    if rule.get("suppressor_effect", "lof") != "lof":
        return False
    if not is_loss_of_function(variant):
        return False
    # The eis rule is specifically about a *coding* loss of function: a promoter
    # variant cannot knock the protein out, it over-expresses it, which is the
    # opposite direction and is the thing being suppressed.
    return is_coding_variant(variant)


def _matches_suppressed(variant: Variant, rule: Dict[str, Any]) -> bool:
    if not _gene_is(variant.gene, rule.get("suppressed_gene", "")):
        return False
    region = rule.get("suppressed_region")
    if region == "promoter":
        return is_promoter_variant(variant)
    if region is None:
        # A loss of function in the suppressor gene does not suppress itself.
        return not _matches_suppressor(variant, rule)
    return True


def epistasis_suppressions(variants: Sequence[Variant],
                           drugs: Optional[Iterable[str]] = None) -> List[Suppression]:
    """Every epistasis rule that fires over a sample's variants.

    One :class:`Suppression` per (rule, drug), so a report can print the
    bedaquiline and the clofazimine consequences of the same *mmpL5* knockout as
    the two separate clinical statements they are.

    *drugs* restricts the output to a drug panel; None means every drug the rules
    name. Nothing is mutated and nothing is dropped: a rule that fires with
    ``confident=False`` is still returned, because "we found the pattern and
    chose not to act on it" is information the reader needs.
    """
    wanted = None if drugs is None else set(normalise_drug(d) for d in drugs)
    found: List[Suppression] = []

    for rule in EPISTASIS_RULES:
        suppressors = [v for v in variants if _matches_suppressor(v, rule)]
        if not suppressors:
            continue
        suppressed = [v for v in variants if _matches_suppressed(v, rule)]
        if not suppressed:
            continue

        # Both sides major, or the rule is reported rather than applied.
        confident = True
        for left in suppressors:
            for right in suppressed:
                left_major = is_major_variant(left.allele_fraction)
                right_major = is_major_variant(right.allele_fraction)
                if left_major is None or right_major is None:
                    continue
                if left_major != right_major:
                    confident = False

        for drug in rule.get("drugs", ()):
            canonical = normalise_drug(drug)
            if wanted is not None and canonical not in wanted:
                continue
            found.append(Suppression(
                rule=_rule_id(rule),
                drug=canonical,
                suppressor_gene=str(rule.get("suppressor_gene", "")),
                suppressor_variants=sorted(
                    (_keys(v)[0] for v in suppressors), key=natural_key),
                suppressed_gene=str(rule.get("suppressed_gene", "")),
                suppressed_variants=sorted(
                    (_keys(v)[0] for v in suppressed), key=natural_key),
                why=str(rule.get("why", "")),
                source=source_for("epistasis_rules") or SRC_WHO_V2,
                confident=confident,
                caveat="" if confident else MIXED_SUBPOPULATION_CAVEAT,
            ))

    return found
