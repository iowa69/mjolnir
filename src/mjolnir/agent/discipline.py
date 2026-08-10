"""The rules a model answer has to obey, checked in code, one sentence at a time.

The verdicts in a Mjolnir report are rules; the prose is the model. This module
is the boundary between the two. It reads an answer and decides whether it may
be printed, and when it may not, the report prints the rule-derived summary and
names the reason. Nothing here asks the model to behave — the system prompt does
that, and a system prompt is a request, not a guarantee.

What is rejected:

- a number the input did not contain, including one computed from two numbers
  that it did — the discipline is to cite, not to calculate;
- anything in ``not_measured`` described as absent, normal, clean or fine;
- "susceptible", where the rule said *no resistance determinant detected* —
  absence of a determinant is absence of evidence, and the two sentences send
  patients to different regimens;
- a claim that reverses a drug call or the sample-validity verdict;
- an organism named that the evidence never contained, which is the §6 rule
  against printing a species that the ANI and barcode evidence did not support;
- for a gate, anything outside the closed option set.

**Everything works a sentence at a time, and concessive clauses are exempt.**
That is the expensive part, inherited whole from tesseract-ai. Its first
discipline pass searched entire answers for words, and failed in both directions
at once. It missed almost everything — matching the literal check name beside a
literal "no", so "contamination" was caught and "contaminant" was not, and every
paraphrase ("came back clean", "none detected", "the isolate appears pure")
walked through. And it rejected good work: the verdict rule fired on the word
"acceptable" anywhere in the text, discarding 86 of 96 valid readings shaped
"While the genome size and contiguity are acceptable, the low read depth means
consensus errors may be present" — a correct, nuanced reading that *agrees* with
the verdict it was accused of contradicting. (That count is as reported for the
run that found the bug; ``tesseract_ai/discipline.py`` records 72 of 87 for the
corpus subset it measured. The two figures are the same failure counted over
different sets, and neither is a threshold — nothing here is derived from
them.)

So: split into sentences, decide what each sentence is about, let a concession
stand, and know each check by its synonyms rather than by its identifier. A
guard that cannot tell a concessive clause from a verdict is not measuring
discipline, it is measuring vocabulary.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from mjolnir import config
from mjolnir.records import (CALL_NO_CALL, CALL_UNCERTAIN, MIXTURE_NOT_ASSESSED,
                             VALIDITY_INVALID, VALIDITY_NOT_ASSESSED,
                             VALIDITY_SUSPECT, VALIDITY_VALID, CohortResult,
                             SampleResult)

#: Numbers an answer may use without its input containing them: small integers
#: doing duty as ordinals or as counts of things the answer itself lists ("two
#: of the three catalogues").
ALLOWED_BARE = frozenset(("0", "1", "2", "3", "4", "5"))

_NUMBER = re.compile(r"\d[\d,]*\.?\d*")
_SENTENCE = re.compile(r"[^.!?;\n]+[.!?;\n]?")

# ---------------------------------------------------------------------------
# What a check is called when a human writes about it
# ---------------------------------------------------------------------------

#: Topic -> the words prose actually uses for it. Checks are matched to topics
#: by substring on their name, so this survives whichever identifiers the QC,
#: contamination and typing modules end up choosing: a check called
#: ``mean_depth``, ``depth_at_target`` or ``median_depth_1x`` all land on
#: ``depth`` and all of them are recognised in a sentence that says "coverage".
CHECK_TOPICS: Dict[str, Tuple[str, ...]] = {
    "depth": ("depth", "coverage", "sequenced to", "reads per site"),
    "breadth": ("breadth", "genome covered", "reference covered", "callable "
                "fraction", "covered at"),
    "mapped": ("mapped", "mapping rate", "alignment rate", "reads aligned",
               "on-target"),
    "evenness": ("evenness", "even coverage", "uniformity", "uneven"),
    "gc": ("gc content", "gc"),
    "unambiguous": ("unambiguous", "ambiguous base", "mixed base",
                    "majority rule"),
    "het": ("heterozygos", "hsnp", "het snp", "minor allele", "mixed infection",
            "mixture", "mixed"),
    "f2": ("f2", "lineage-defining", "lineage defining"),
    "f47": ("f47", "lineage-defining", "lineage defining"),
    "contamination": ("contaminat", "purity", "pure", "second organism",
                      "secondary organism", "foreign read", "off-target read",
                      "other species", "co-infection", "coinfection",
                      "another organism"),
    "kraken": ("kraken", "read classifier", "taxonomic classif", "screen"),
    "target": ("non-target", "non target", "off-target"),
    "barcode": ("barcode", "lineage support", "defining snp", "defining site"),
    "species": ("species", "ani", "identification", "identified"),
    "lineage": ("lineage", "sublineage", "strain family"),
    "minor": ("minor variant", "heteroresistance", "subpopulation",
              "minority population"),
    "duplicate": ("duplicate", "pcr duplicate"),
    "quality": ("base quality", "read quality", "phred", "q30"),
    "length": ("read length",),
    "callable": ("callable", "target region", "gene coverage", "evaluable"),
    "cluster": ("cluster", "transmission", "snp distance", "linkage"),
    "mask": ("mask", "repetitive region", "pe/ppe"),
    "shared": ("shared callable", "denominator", "sites compared"),
    "catalogue": ("catalogue", "catalog", "who grade", "grading"),
    "validity": ("validity", "valid sample", "usable sample", "fit for"),
}

#: Check names whose topic is not visible in the name. Substring match on the
#: left, topics on the right.
_NAME_HINTS: Dict[str, Tuple[str, ...]] = {
    "purity": ("contamination",),
    # A non-target read fraction is the contamination measurement; prose calls
    # it "a second organism" and never calls it by its identifier.
    "non target": ("contamination", "target"),
    "foreign": ("contamination",),
    "kraken": ("kraken", "contamination"),
    "screen": ("kraken", "contamination"),
    "f2": ("f2", "het"),
    "f47": ("f47", "het"),
    "ani": ("species",),
    "mixture": ("het",),
    "hsnp": ("het",),
    "unambig": ("unambiguous",),
    "even": ("evenness",),
    "barcode": ("barcode", "lineage"),
    "distance": ("cluster",),
}

#: Words that assert a thing is absent, settled or fine. ``none`` and
#: ``nothing`` are the ones a naive version misses entirely.
_SETTLED = re.compile(
    r"\b(no|none|neither|nothing|nil|zero|without|free of|absent|clean|pure|"
    r"normal|fine|negative|unremarkable|sound|satisfactory|good|excellent|"
    r"not present|not detected|no evidence|no sign|ruled out|excluded)\b",
    re.IGNORECASE)

#: Phrasings that correctly report an absence of *measurement*. A sentence
#: carrying one of these is hedged and may say whatever else it likes:
#: "contamination was not screened, so purity is unknown" mentions purity and
#: says unknown in the same breath, and must pass.
_HEDGED = re.compile(
    r"\b(not screened|never screened|not measured|unmeasured|not computed|"
    r"not calculated|not counted|not assessed|not available|not established|"
    r"not evaluated|not evaluable|no reference|unknown|undetermined|"
    r"cannot be|could not be|was not run|did not run|not informative|"
    r"uninformative|unverified|not verified|absence of evidence|"
    r"rather than absent|not the same as|remains? (unknown|unverified)|"
    r"no measurement|nothing was measured)\b", re.IGNORECASE)

#: A sentence that starts by conceding, or turns on a contrastive conjunction,
#: is not making a global claim whatever else it contains. This single
#: exemption is what saved 86 of 96 readings in the tesseract-ai corpus.
_CONCESSIVE = re.compile(
    r"^\s*(while|although|though|whilst|despite|notwithstanding|even though|"
    r"aside from|other than|apart from|beyond)\b|"
    r"\b(but|however|nevertheless|nonetheless|yet|except|other than|apart from|"
    r"beyond that)\b", re.IGNORECASE)

_NEGATED = re.compile(
    r"\b(no|not|never|cannot|can't|without|neither|nor|absence of|lacks?|"
    r"does not|do not|is not|are not|was not|were not|rather than|"
    r"must not|should not|may not)\b", re.IGNORECASE)

#: "the sample is fine" as an actual predicate of an actual subject, not the
#: words "sample" and "acceptable" somewhere in the same sentence. Scoping to
#: the construction is the second correction the corpus forced: a version that
#: looked for a subject and a positive adjective anywhere still fired on
#: "insufficient to support **reliable** consensus calling", which is a negative
#: claim containing a positive word, and that is how people write.
_GLOBAL_SUBJECT = r"(sample|isolate|result|run|data|sequencing|analysis|report)"
_GLOBAL_POSITIVE = (r"(fit for (purpose|use|this)|suitable|good enough|"
                    r"acceptable|adequate|usable|reliable|excellent|valid|"
                    r"of high quality|sound|clean|fine|trustworthy)")
_GLOBAL_CLAIM = re.compile(
    r"\b(the|this|it)\s+" + _GLOBAL_SUBJECT + r"?\s*"
    r"\b(is|was|are|were|appears? to be|seems? to be|remains?|looks?)\s+"
    r"(?!not\b|un|in|insufficient|inadequate|unsuitable|unfit|unreliable|"
    r"unacceptable|invalid|too\b|barely\b|hardly\b|only\b|marginal)"
    r"(a |an |very |fully |entirely )?" + _GLOBAL_POSITIVE, re.IGNORECASE)
_NO_PROBLEMS = re.compile(
    r"\bno\s+(major\s+|significant\s+|apparent\s+)?"
    r"(issues|problems|concerns|defects|caveats)\b", re.IGNORECASE)

#: A susceptibility claim, as a predicate. "susceptibility testing is required"
#: does not match, and must not: it is the sentence the report wants.
_SUSCEPTIBLE = re.compile(
    r"\b(susceptible|sensitive)\s+to\b|"
    r"\b(is|are|was|were|remains?|appears? to be|likely|probably)\s+"
    r"(fully\s+|still\s+|therefore\s+)?(susceptible|sensitive)\b|"
    r"\b(retains?|retain)\s+(full\s+)?(susceptibility|activity)\b",
    re.IGNORECASE)

#: A positive assertion of resistance, for the drugs the rules did not call R.
_RESISTANT_CLAIM = re.compile(
    r"\b(resistant to|resistance to|confers resistance|predicted resistan|"
    r"is resistant|are resistant|resistance is predicted|"
    r"resistance determinant (was |were )?(detected|found|present))\b",
    re.IGNORECASE)

#: Sentences that clear a drug the rules called resistant.
_CLEARS_DRUG = re.compile(
    r"\b(remains?|stays?|is still|are still)\s+(active|effective|usable|an "
    r"option)\b|\bcan (still )?be used\b|\bno resistance\b|"
    r"\bnot resistant\b|\bwithout resistance\b|\bactivity is (retained|preserved)\b",
    re.IGNORECASE)

#: `M. bovis`, `Mycobacterium chimaera`, `Mycobacteroides abscessus`.
_ORGANISM = re.compile(
    r"\b(M|Mycobacterium|Mycobacteroides|Mycolicibacterium)\.?\s+"
    r"([a-z][a-z\-]{3,})\b")

#: Short drug aliases that are also ordinary words or fragments; matching them
#: as prose tokens produces false positives, and the full name still matches.
_AMBIGUOUS_DRUG_ALIASES = frozenset(("cap", "str", "pa", "sm", "eth", "km"))


# ---------------------------------------------------------------------------
# Text mechanics
# ---------------------------------------------------------------------------

def sentences(text: str) -> List[str]:
    return [s.strip() for s in _SENTENCE.findall(text or "") if s.strip()]


def numbers_in(text: str) -> Set[str]:
    """Every number in a string, normalised so 5,015,531 and 5015531 match."""
    return set(m.group(0).replace(",", "").rstrip(".")
               for m in _NUMBER.finditer(text or ""))


def _is_rounding_of(token: str, supported: str) -> bool:
    """Whether *token* is *supported* written to fewer decimal places.

    This used to be a string-prefix test, which is the same thing for 26.2
    against 26.24 and catastrophically not the same thing for 180 against 18:
    given an input of 18x, "depth is 180x" was accepted, as was "4000 variants"
    from an input containing 40. Rounding is a numeric relation, so it is tested
    numerically, and an integer may never gain digits.
    """
    try:
        answer_value, source_value = float(token), float(supported)
    except (TypeError, ValueError):
        return False
    if answer_value == source_value:
        return True
    decimals = len(token.partition(".")[2])
    if len(supported.partition(".")[2]) <= decimals:
        return False
    return round(source_value, decimals) == answer_value


def unsupported_numbers(answer: str, context: str) -> Set[str]:
    """Numbers the answer asserts that its input never gave it.

    A percentage computed from two supplied numbers is still a number that was
    not supplied. The observation carries the percentages the report speaks in
    precisely so that citing one is possible without calculating it.
    """
    supported = numbers_in(context)
    out: Set[str] = set()
    for token in numbers_in(answer):
        if token in supported or token in ALLOWED_BARE:
            continue
        if any(_is_rounding_of(token, value) for value in supported):
            continue
        out.add(token)
    return out


def checkable_prose(answer: str) -> str:
    """The parts of an answer that make claims, without the JSON scaffolding.

    The number rule must not see fields the system itself demands: a
    ``confidence`` the model was told to emit is not a claim about the evidence,
    and reading it as one rejected more than half of otherwise-good decisions on
    the first corpus run.
    """
    try:
        parsed = json.loads(answer)
    except (ValueError, TypeError):
        return answer or ""
    if not isinstance(parsed, dict):
        return answer or ""
    parts: List[str] = []
    for key in ("headline", "body", "reason", "reading", "note", "summary"):
        value = parsed.get(key)
        if isinstance(value, str):
            parts.append(value)
    return "\n".join(parts) if parts else (answer or "")


def topics_for(check_name: str) -> Tuple[str, ...]:
    """Which topics a check name belongs to, by substring rather than by equality."""
    lowered = str(check_name).replace("_", " ").lower()
    topics: List[str] = []
    for topic in CHECK_TOPICS:
        if topic.replace("_", " ") in lowered:
            topics.append(topic)
    for fragment, extra in _NAME_HINTS.items():
        if fragment in lowered:
            topics.extend(t for t in extra if t not in topics)
    return tuple(topics)


def _contains_term(sentence: str, term: str) -> bool:
    """Substring for long terms, word-boundary for short ones.

    ``gc`` and ``f2`` inside another word are not mentions of GC content or of
    the F2 lineage set, and a plain ``in`` test says they are.
    """
    if len(term) <= 4:
        return re.search(r"\b" + re.escape(term) + r"\b", sentence,
                         re.IGNORECASE) is not None
    return term in sentence.lower()


def mentions(sentence: str, check_name: str) -> bool:
    """Whether a sentence is talking about this check, by name or by paraphrase."""
    lowered = sentence.lower()
    name = str(check_name).lower()
    if name and (name in lowered or name.replace("_", " ") in lowered):
        return True
    for topic in topics_for(check_name):
        for term in CHECK_TOPICS.get(topic, ()):
            if _contains_term(sentence, term):
                return True
    return False


def drug_terms(drug: str) -> Tuple[str, ...]:
    """The names prose uses for one drug: its own, its code, its aliases."""
    canonical = config.normalise_drug(drug) or drug
    terms = set([str(drug).lower(), canonical.lower()])
    code = config.DRUG_CODES.get(canonical, "")
    if code:
        terms.add(code.lower())
    for alias, target in config.DRUG_ALIASES.items():
        if target == canonical and len(alias) >= 3 \
                and alias not in _AMBIGUOUS_DRUG_ALIASES:
            terms.add(alias)
    return tuple(sorted(t for t in terms if t))


def mentions_drug(sentence: str, drug: str) -> bool:
    return any(_contains_term(sentence, term) for term in drug_terms(drug))


# ---------------------------------------------------------------------------
# The rules
# ---------------------------------------------------------------------------

def treats_as_settled(prose: str, check_name: str) -> Optional[str]:
    """The sentence that calls an unmeasured check absent or fine, if there is one.

    Sentence-scoped and hedge-aware, which is the whole point. "The Kraken2
    screen was not informative, so purity is unknown" mentions purity and says
    unknown; it passes. "The contamination screen came back clean" is the
    failure, and a word-level guard lets it through because it is looking for
    the literal word "no".
    """
    for sentence in sentences(prose):
        if not mentions(sentence, check_name):
            continue
        if _HEDGED.search(sentence):
            continue
        if _SETTLED.search(sentence):
            return sentence
    return None


def claims_susceptibility(prose: str, drugs: Iterable[str]) -> Optional[str]:
    """A sentence calling something susceptible where the rules did not.

    Both forms are caught: naming a drug the rules gave no determinant for, and
    the unqualified "the isolate is susceptible" with no drug attached. Negated
    and hedged sentences are exempt, because "this is not a prediction of
    susceptibility" is the sentence the report most wants written.
    """
    names = list(drugs)
    for sentence in sentences(prose):
        match = _SUSCEPTIBLE.search(sentence)
        if not match:
            continue
        if _NEGATED.search(sentence) or _HEDGED.search(sentence):
            continue
        if not names or any(mentions_drug(sentence, drug) for drug in names):
            return sentence
    return None


def contradicts_drug_calls(prose: str, resistant: Iterable[str],
                           not_resistant: Iterable[str]) -> Optional[Tuple[str, str]]:
    """A sentence that reverses a drug call, with the reason it was rejected.

    Qualifying a call is allowed and expected — "the rifampicin call rests on a
    single indel on ONT, where 16.6% of indel calls were uncorroborated" agrees
    with the call while weakening it. Reversing it is not.
    """
    resistant_names = list(resistant)
    other_names = list(not_resistant)
    for sentence in sentences(prose):
        if _CONCESSIVE.search(sentence):
            continue
        for drug in resistant_names:
            if not mentions_drug(sentence, drug):
                continue
            if _CLEARS_DRUG.search(sentence) or (
                    _SUSCEPTIBLE.search(sentence) and not _NEGATED.search(sentence)):
                return sentence, ("clears {0}, which the rules called "
                                  "resistant".format(drug))
        for drug in other_names:
            if not mentions_drug(sentence, drug):
                continue
            if _RESISTANT_CLAIM.search(sentence) and not _NEGATED.search(sentence):
                return sentence, ("asserts resistance to {0}, which the rules "
                                  "did not call resistant".format(drug))
    return None


def contradicts_validity(prose: str, verdict: str) -> Optional[str]:
    """A sentence calling the whole sample fine when the verdict says otherwise.

    Concessive sentences are skipped, and for ``suspect`` — which already
    concedes that something is wrong — only an unqualified superlative counts.
    """
    if verdict not in (VALIDITY_INVALID, VALIDITY_SUSPECT, VALIDITY_NOT_ASSESSED):
        return None
    for sentence in sentences(prose):
        if _CONCESSIVE.search(sentence):
            continue
        if _NO_PROBLEMS.search(sentence):
            return sentence
        if verdict == VALIDITY_SUSPECT:
            if re.search(r"\b(the|this)\s+" + _GLOBAL_SUBJECT +
                         r"\s+(is|was)\s+(excellent|fully (valid|suitable)|"
                         r"of high quality|clean|pure)\b", sentence, re.IGNORECASE):
                return sentence
            continue
        if verdict == VALIDITY_NOT_ASSESSED:
            # Nothing was measured, so any global endorsement is unmeasured.
            if _GLOBAL_CLAIM.search(sentence) and not _HEDGED.search(sentence):
                return sentence
            continue
        if _GLOBAL_CLAIM.search(sentence):
            return sentence
    return None


def invented_organisms(prose: str, context: str) -> Set[str]:
    """Organism names the answer states that its evidence never contained.

    §6 of the design in one rule. MTBC members sit at 99.21-99.92% ANI and are
    heterotypic synonyms of *M. tuberculosis*, so a species name that the ANI
    and barcode evidence did not produce is an invention however plausible it
    reads — and "M. bovis" printed from nowhere is exactly the failure the
    design forbids.
    """
    supported = set(m.group(2).lower() for m in _ORGANISM.finditer(context or ""))
    return set(m.group(2).lower() for m in _ORGANISM.finditer(prose or "")
               if m.group(2).lower() not in supported)


def chosen_option(answer: str, options: Sequence[str]) -> Optional[str]:
    """The option the answer chose, or None if it chose outside the closed set."""
    try:
        parsed = json.loads(answer)
    except (ValueError, TypeError):
        return None
    if not isinstance(parsed, dict):
        return None
    choice = str(parsed.get("choice", "")).strip()
    return choice if choice in list(options) else None


# ---------------------------------------------------------------------------
# The facts an answer is checked against
# ---------------------------------------------------------------------------

@dataclass
class RuleFacts:
    """The rule-derived state of one result, in the shape the rules need it.

    Built from a :class:`~mjolnir.records.SampleResult` rather than passed
    around by hand, so that a drug call added anywhere in the pipeline is
    policed here without anyone remembering to wire it up.
    """

    unmeasured: List[str] = field(default_factory=list)
    resistant_drugs: List[str] = field(default_factory=list)
    no_determinant_drugs: List[str] = field(default_factory=list)
    uncertain_drugs: List[str] = field(default_factory=list)
    unevaluable_drugs: List[str] = field(default_factory=list)
    validity: str = VALIDITY_NOT_ASSESSED
    mixture_class: str = MIXTURE_NOT_ASSESSED
    screen_informative: bool = False
    resolved_to_species: bool = False
    platform: str = ""

    @property
    def not_resistant_drugs(self) -> List[str]:
        """Everything the rules did not call resistant, for the reverse check."""
        return sorted(set(self.no_determinant_drugs) | set(self.uncertain_drugs)
                      | set(self.unevaluable_drugs))


def facts_from_sample(result: SampleResult) -> RuleFacts:
    return RuleFacts(
        unmeasured=list(result.unmeasured()),
        resistant_drugs=[d.drug for d in result.resistant_drugs()],
        no_determinant_drugs=[d.drug for d in result.drugs
                              if d.call == CALL_NO_CALL
                              and d.target_covered is not False],
        uncertain_drugs=[d.drug for d in result.drugs if d.call == CALL_UNCERTAIN],
        unevaluable_drugs=[d.drug for d in result.drugs
                           if d.target_covered is False],
        validity=result.contamination.verdict,
        mixture_class=result.contamination.mixture_class,
        screen_informative=result.contamination.screen_informative,
        resolved_to_species=result.species.resolved_to_species,
        platform=result.platform,
    )


def facts_from_cohort(cohort: CohortResult) -> RuleFacts:
    """A cohort has no drug calls; it has checks, and a mask that may be absent."""
    return RuleFacts(
        unmeasured=[c.name for c in cohort.checks if not c.measured],
        validity=VALIDITY_NOT_ASSESSED if not cohort.samples else VALIDITY_VALID,
    )


@dataclass
class Verdict:
    """Whether an answer may be printed, and if not, why — in the report's words.

    ``reason`` is written to be shown: the design requires the report to say
    that a model answer was discarded and to name the rule it broke, so the
    string has to make sense to a clinician reading a footnote.
    """

    ok: bool = True
    reason: str = ""
    rule: str = ""
    sentence: str = ""

    def __bool__(self) -> bool:
        return self.ok


def review(answer: str, context: str, facts: Optional[RuleFacts] = None,
           options: Sequence[str] = (), expected_choice: str = "") -> Verdict:
    """Judge one answer against the evidence it was given.

    Order matters only in what gets reported first; every rule is independent.
    Numbers come first because a fabricated number is the failure a reader is
    least able to catch.
    """
    if not (answer or "").strip():
        return Verdict(False, "the model returned nothing", "empty")

    facts = facts or RuleFacts()
    prose = checkable_prose(answer)

    invented = unsupported_numbers(prose, context)
    if invented:
        return Verdict(False,
                       "states numbers that are not in the evidence: {0}".format(
                           ", ".join(sorted(invented)[:5])),
                       "unsupported-number")

    for name in facts.unmeasured:
        sentence = treats_as_settled(prose, name)
        if sentence:
            return Verdict(False,
                           "describes {0!r}, which was not measured, as absent "
                           "or normal".format(name),
                           "unmeasured-as-settled", sentence)

    sentence = claims_susceptibility(
        prose, facts.no_determinant_drugs + facts.unevaluable_drugs)
    if sentence:
        return Verdict(False,
                       "calls the isolate susceptible where the rules found no "
                       "resistance determinant; absence of a determinant is not "
                       "susceptibility", "susceptible-for-no-determinant", sentence)

    contradiction = contradicts_drug_calls(prose, facts.resistant_drugs,
                                           facts.not_resistant_drugs)
    if contradiction:
        sentence, why = contradiction
        return Verdict(False, "contradicts a rule-derived drug call: {0}".format(why),
                       "contradicts-drug-call", sentence)

    sentence = contradicts_validity(prose, facts.validity)
    if sentence:
        return Verdict(False,
                       "declares the sample sound where the rule-derived "
                       "validity verdict is {0!r}".format(facts.validity),
                       "contradicts-validity", sentence)

    organisms = invented_organisms(prose, context)
    if organisms:
        return Verdict(False,
                       "names an organism the evidence does not contain: {0}".format(
                           ", ".join(sorted(organisms)[:3])),
                       "invented-organism")

    if options:
        choice = chosen_option(answer, options)
        if choice is None:
            return Verdict(False, "did not choose an option from the closed set",
                           "outside-option-set")
        if expected_choice and choice != expected_choice:
            return Verdict(False,
                           "chose {0!r} where the rules require {1!r}".format(
                               choice, expected_choice),
                           "overrides-rule")
    return Verdict(True)
