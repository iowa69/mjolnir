"""Exactly what the model is allowed to see, and the rule-derived text it writes over.

Three things are built here, and they are one idea.

**The observation.** Finished checks with their thresholds and the paper each
threshold came from; gene names and HGVS variant names; per-drug calls with
their WHO grade and the catalogues that agreed; coverage, purity and lineage
metrics; the caveats that already apply. Tables and names — never bases. The
rule that the model never sees raw sequence is enforced here in code and not
asked for in a prompt: :func:`assert_no_sequence` walks every string in the
structure and raises :class:`SequenceLeak` on a gene-length nucleotide run, so a
tool that starts putting a consensus into a note breaks the build rather than
quietly widening what the model reads.

**The rule-derived summary.** Every verdict in a Mjolnir report is computed in
Python before a model is contacted, and :func:`rule_summary` renders those
verdicts as prose. It is handed to the model as the statement to write over, and
it is what replaces the model's answer when the discipline rules discard it. The
same text in both roles is what makes the fallback honest: nothing is lost when
the model is absent except fluency.

**The playbook.** Which drugs matter for this organism, which caveats always
apply, and what a clinician has to be told. That is organism knowledge, not
program logic, so it lives in ``playbooks/*.yaml`` and is validated at load time
— a playbook whose gate default is not among its own options must fail before a
run starts, not three hours in.

The observation is size-capped, and a capped observation says so. Silent
truncation would let the model reason from a table it believes is complete,
which is the same failure as an unmeasured metric printed as a pass.

YAML is read with PyYAML when it is installed. It is not a Mjolnir dependency —
``pyproject.toml`` carries pandas, numpy and openpyxl, and a prose layer must
not add a parser to a resistance-calling tool — so a small strict reader for the
subset the shipped playbooks use stands in. That reader refuses anything it does
not fully understand instead of guessing, because a half-parsed playbook is a
silently wrong set of clinical caveats.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from mjolnir import config
from mjolnir.records import (CALL_NO_CALL, CALL_UNCERTAIN, MIXTURE_NOT_ASSESSED,
                             NO_DETERMINANT_TEXT, PLATFORM_FASTA, PLATFORM_ONT,
                             STATUS_FAIL, STATUS_WARN, VALIDITY_NOT_ASSESSED,
                             VALIDITY_VALID, Check, CohortResult, SampleResult,
                             SpeciesCall)
from mjolnir.utils import (MjolnirError, looks_like_sequence, natural_key,
                           percentage, plural, round_or_none, to_jsonable)

# ---------------------------------------------------------------------------
# Limits
#
# Presentation limits, not scientific thresholds: nothing in a report is derived
# from them, so they are not in the config.py registry. Each says why it holds
# the value it does, and every one of them is announced when it bites.
# ---------------------------------------------------------------------------

#: SOURCE: agent policy. A 32k-token context (agent/client.py) is roughly 100 kB
#: of English; a quarter of it leaves the model room to write and keeps the
#: prompt cache warm across the samples of a cohort.
MAX_OBSERVATION_BYTES = 24000

#: SOURCE: agent policy. The annex prints every variant; the reading needs the
#: graded ones. Catalogue-relevant variants are kept first and the count of what
#: was dropped is stated in the observation itself.
MAX_VARIANT_ROWS = 60

#: Keys that could only ever hold sequence. The nucleotide-run test below is the
#: real guard; this catches an empty or short field named as sequence, which is
#: a tool about to grow into a leak. ``total_reads`` and ``mapped_reads`` are
#: counts and are deliberately not in this set.
FORBIDDEN_KEYS = frozenset((
    "sequence", "seq", "sequences", "fasta", "fastq", "bases",
    "raw_sequence", "consensus_sequence", "read_sequence", "contig_sequence",
))


class SequenceLeak(MjolnirError):
    """An observation would have carried raw sequence to the model.

    A subclass of :class:`~mjolnir.utils.MjolnirError` so the CLI prints it
    without a traceback, but it is a programming error rather than an operator
    error: the fix is in whichever module put bases in a field, not in anything
    the user typed.
    """


def assert_no_sequence(obj: Any, path: str = "observation") -> None:
    """Raise if any string anywhere in *obj* is a gene-length nucleotide run.

    Walks keys as well as values, and names the path it found the leak at, since
    "the observation contains sequence" is not actionable and
    ``observation.variants[3].note`` is.
    """
    if isinstance(obj, dict):
        for key, value in obj.items():
            if str(key).strip().lower() in FORBIDDEN_KEYS:
                raise SequenceLeak(
                    "{0}.{1} is a sequence field; the model is given tables, "
                    "gene names and HGVS names, never bases".format(path, key))
            assert_no_sequence(value, "{0}.{1}".format(path, key))
    elif isinstance(obj, (list, tuple)):
        for index, value in enumerate(obj):
            assert_no_sequence(value, "{0}[{1}]".format(path, index))
    elif isinstance(obj, str) and looks_like_sequence(obj):
        raise SequenceLeak(
            "{0} carries a run of {1} nucleotides; the model is given tables, "
            "gene names and HGVS names, never bases".format(path, len(obj)))


# ---------------------------------------------------------------------------
# The observation
# ---------------------------------------------------------------------------

@dataclass
class Observation:
    """One structure, containing everything the model may read and nothing else.

    ``truncated`` and ``note`` exist so a shrunken observation announces itself.
    A model told it has the whole variant table when it has half of it will
    write "no other catalogued variants were found", which is a claim nobody
    measured.
    """

    kind: str = "sample"
    subject: str = ""
    data: Dict[str, Any] = field(default_factory=dict)
    truncated: bool = False
    note: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {"kind": self.kind, "subject": self.subject, "data": self.data,
                "truncated": self.truncated, "note": self.note}

    def json(self, indent: Optional[int] = None) -> str:
        return json.dumps(self.to_dict(), sort_keys=True, indent=indent,
                          default=str)

    def nbytes(self) -> int:
        return len(self.json().encode("utf-8"))

    def capped(self, limit: int = MAX_OBSERVATION_BYTES) -> "Observation":
        """A copy inside *limit* bytes, saying what it had to drop.

        Shrinks the largest list or dict first — a per-drug table, a variant
        list, a barcode-support table — halving it until the whole fits, and
        records the counts. Dicts count as well as lists: a per-catalogue map
        can carry thousands of keys and no list at all.
        """
        assert_no_sequence(self.data)
        out = Observation(kind=self.kind, subject=self.subject,
                          data=json.loads(json.dumps(self.data, default=str)),
                          truncated=self.truncated, note=self.note)
        if out.nbytes() <= limit:
            return out

        dropped: List[str] = []
        while out.nbytes() > limit:
            biggest = _largest_container(out.data)
            if biggest is None or len(biggest) <= 1:
                break
            before = len(biggest)
            keep = max(1, before // 2)
            if isinstance(biggest, list):
                del biggest[keep:]
            else:
                for key in list(biggest)[keep:]:
                    del biggest[key]
            dropped.append("{0}->{1}".format(before, keep))
        out.truncated = True
        out.note = (out.note + " " if out.note else "") + (
            "observation truncated to fit {0} bytes; tables were halved ({1}). "
            "Absence of a row here is truncation, not a measurement.".format(
                limit, ", ".join(dropped) or "nothing shrinkable"))
        return out


def _largest_container(obj: Any) -> Optional[Any]:
    best: Optional[Any] = None
    best_len = 1
    stack: List[Any] = [obj]
    while stack:
        current = stack.pop()
        if isinstance(current, dict):
            if len(current) > best_len:
                best, best_len = current, len(current)
            stack.extend(current.values())
        elif isinstance(current, list):
            if len(current) > best_len:
                best, best_len = current, len(current)
            stack.extend(current)
    return best


# ---------------------------------------------------------------------------
# Playbooks
# ---------------------------------------------------------------------------

PLAYBOOK_DIR = Path(__file__).resolve().parent / "playbooks"

_PLAYBOOK_KEYS = frozenset((
    "name", "organism", "describe", "audience", "applies_to", "drugs",
    "always_caveats", "must_state", "must_not_say", "headline_focus", "gates",
))
_DRUG_KEYS = frozenset(("drug", "why", "priority"))
_GATE_KEYS = frozenset(("id", "question", "options", "default", "describe"))


@dataclass
class DrugFocus:
    """One drug this organism's reader is expected to look for, and why."""

    drug: str
    why: str = ""
    priority: int = 5


@dataclass
class Gate:
    """A question the rules genuinely cannot settle, with a closed answer set.

    The closed set is the leash: the model picks among alternatives someone has
    thought about, and a declared default is taken whenever it does not. The
    default is stored beside the options rather than assumed to be the first one
    — those differ, and taking ``options[0]`` while recording "took the declared
    default" is a lie in exactly the cases where it matters.
    """

    id: str
    question: str
    options: List[str] = field(default_factory=list)
    default: str = ""
    describe: str = ""


@dataclass
class Playbook:
    """The organism-specific reading, as data.

    ``must_state`` and ``must_not_say`` are the clinical content of the report
    that does not come from a threshold: that absence of a determinant is not
    susceptibility, that a BCG call carries intrinsic pyrazinamide resistance,
    that an ONT run under-detects minor variants. They are handed to the model
    as instructions *and* used by the report when the model is absent.
    """

    name: str = ""
    organism: str = ""
    describe: str = ""
    audience: str = ""
    applies_to_complex: List[str] = field(default_factory=list)
    applies_to_species: List[str] = field(default_factory=list)
    drugs: List[DrugFocus] = field(default_factory=list)
    always_caveats: List[str] = field(default_factory=list)
    must_state: List[str] = field(default_factory=list)
    must_not_say: List[str] = field(default_factory=list)
    headline_focus: str = ""
    gates: List[Gate] = field(default_factory=list)
    path: str = ""

    def drug_names(self) -> List[str]:
        return [d.drug for d in sorted(self.drugs, key=lambda d: (d.priority, d.drug))]

    def gate(self, gate_id: str) -> Optional[Gate]:
        for gate in self.gates:
            if gate.id == gate_id:
                return gate
        return None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "organism": self.organism,
            "describe": self.describe,
            "audience": self.audience,
            "drugs_that_matter": [
                {"drug": d.drug, "why": d.why}
                for d in sorted(self.drugs, key=lambda d: (d.priority, d.drug))],
            "caveats_that_always_apply": list(self.always_caveats),
            "must_state": list(self.must_state),
            "must_not_say": list(self.must_not_say),
            "headline_focus": self.headline_focus,
        }


def _require_keys(mapping: Dict[str, Any], allowed: frozenset, what: str,
                  where: str) -> None:
    unknown = sorted(set(mapping) - set(allowed))
    if unknown:
        raise MjolnirError(
            "{0}: unknown {1} field(s) {2}; known fields are {3}".format(
                where, what, ", ".join(unknown), ", ".join(sorted(allowed))))


def _as_list(value: Any, where: str, what: str) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    raise MjolnirError("{0}: {1} must be a list, got {2}".format(
        where, what, type(value).__name__))


def playbook_from_dict(raw: Dict[str, Any], where: str = "playbook") -> Playbook:
    """Validate a parsed playbook completely, or refuse it completely."""
    if not isinstance(raw, dict):
        raise MjolnirError("{0}: a playbook must be a mapping".format(where))
    _require_keys(raw, _PLAYBOOK_KEYS, "playbook", where)

    organism = str(raw.get("organism") or "").strip()
    if not organism:
        raise MjolnirError("{0}: `organism` is required".format(where))

    applies = raw.get("applies_to") or {}
    if not isinstance(applies, dict):
        raise MjolnirError("{0}: `applies_to` must be a mapping with `complex` "
                           "and/or `species`".format(where))

    drugs: List[DrugFocus] = []
    for entry in _as_list(raw.get("drugs"), where, "`drugs`"):
        if not isinstance(entry, dict):
            raise MjolnirError("{0}: each `drugs` entry must be a mapping with "
                               "`drug` and `why`".format(where))
        _require_keys(entry, _DRUG_KEYS, "drug", where)
        drug = str(entry.get("drug") or "").strip()
        if not drug:
            raise MjolnirError("{0}: a `drugs` entry has no `drug`".format(where))
        drugs.append(DrugFocus(drug=drug, why=str(entry.get("why") or "").strip(),
                               priority=int(entry.get("priority") or 5)))
    if not drugs:
        raise MjolnirError(
            "{0}: a playbook with no drugs cannot say what matters for this "
            "organism".format(where))

    gates: List[Gate] = []
    for entry in _as_list(raw.get("gates"), where, "`gates`"):
        if not isinstance(entry, dict):
            raise MjolnirError("{0}: each `gates` entry must be a mapping".format(where))
        _require_keys(entry, _GATE_KEYS, "gate", where)
        gate = Gate(id=str(entry.get("id") or "").strip(),
                    question=str(entry.get("question") or "").strip(),
                    options=[str(o) for o in _as_list(entry.get("options"), where,
                                                      "gate `options`")],
                    default=str(entry.get("default") or "").strip(),
                    describe=str(entry.get("describe") or "").strip())
        if not gate.id or not gate.question:
            raise MjolnirError("{0}: a gate needs both `id` and `question`".format(where))
        if len(gate.options) < 2:
            raise MjolnirError(
                "{0}: gate {1!r} needs at least two options; a gate with one "
                "answer is a rule and belongs in Python".format(where, gate.id))
        if gate.default not in gate.options:
            raise MjolnirError(
                "{0}: gate {1!r} declares default {2!r}, which is not among its "
                "options {3}; the default is what a missing model falls back "
                "to, so it cannot be unreachable".format(
                    where, gate.id, gate.default, gate.options))
        gates.append(gate)
    if len(set(g.id for g in gates)) != len(gates):
        raise MjolnirError("{0}: duplicate gate id".format(where))

    return Playbook(
        name=str(raw.get("name") or organism).strip(),
        organism=organism,
        describe=str(raw.get("describe") or "").strip(),
        audience=str(raw.get("audience") or "").strip(),
        applies_to_complex=[str(x).strip() for x in _as_list(
            applies.get("complex"), where, "`applies_to.complex`")],
        applies_to_species=[str(x).strip() for x in _as_list(
            applies.get("species"), where, "`applies_to.species`")],
        drugs=drugs,
        # Folded YAML blocks keep a trailing newline; these strings are printed
        # inline in a prompt and in the report, so they are stripped once here
        # rather than at every use.
        always_caveats=[str(x).strip() for x in _as_list(
            raw.get("always_caveats"), where, "`always_caveats`")],
        must_state=[str(x).strip() for x in _as_list(raw.get("must_state"), where,
                                                     "`must_state`")],
        must_not_say=[str(x).strip() for x in _as_list(raw.get("must_not_say"),
                                                       where, "`must_not_say`")],
        headline_focus=str(raw.get("headline_focus") or "").strip(),
        gates=gates,
        path=where,
    )


def load_playbook(path: Any) -> Playbook:
    """Parse and fully validate a playbook file."""
    resolved = Path(path)
    if not resolved.exists():
        raise MjolnirError(
            "no playbook at {0}; the shipped playbooks are {1}".format(
                resolved, ", ".join(available_playbooks()) or "missing"))
    return playbook_from_dict(load_yaml(resolved.read_text(encoding="utf-8"),
                                        where=str(resolved)),
                              where=str(resolved))


def available_playbooks(directory: Optional[Any] = None) -> List[str]:
    folder = Path(directory or PLAYBOOK_DIR)
    if not folder.exists():
        return []
    return sorted(p.stem for p in folder.glob("*.yaml"))


def playbook_named(name: str, directory: Optional[Any] = None) -> Playbook:
    folder = Path(directory or PLAYBOOK_DIR)
    return load_playbook(folder / "{0}.yaml".format(name))


def playbook_for(species: Optional[SpeciesCall],
                 directory: Optional[Any] = None) -> Playbook:
    """The playbook for this species call: MTBC, otherwise the NTM one.

    Defaulting to NTM rather than to MTBC is deliberate. The MTBC playbook talks
    about WHO grades and lineage barcodes, none of which apply to *M. chimaera*,
    and an unresolved organism is far more likely to be an NTM — MTBC members
    are called from lineage SNPs and would have been resolved. The NTM playbook
    is also the one that says "no evidence base exists for this species-drug
    pair", which is the right thing to say about an organism nobody identified.
    """
    folder = Path(directory or PLAYBOOK_DIR)
    if species is not None:
        complex_name = (species.complex or "").strip().upper()
        name = (species.name or "").strip().lower()
        if species.is_mtbc or complex_name in ("MTBC", "M. TUBERCULOSIS COMPLEX"):
            return playbook_named("mtbc", folder)
        mtbc_playbook = folder / "mtbc.yaml"
        if mtbc_playbook.exists():
            candidate = load_playbook(mtbc_playbook)
            for prefix in candidate.applies_to_species:
                if name.startswith(prefix.strip().lower()) and prefix.strip():
                    return candidate
    return playbook_named("ntm", folder)


# ---------------------------------------------------------------------------
# YAML: PyYAML when present, a strict subset reader when not
# ---------------------------------------------------------------------------

def load_yaml(text: str, where: str = "playbook") -> Dict[str, Any]:
    """Parse playbook YAML, with PyYAML if it is installed.

    The built-in reader handles the subset the shipped playbooks use: nested
    mappings, lists of scalars, lists of single-level mappings, ``|`` literal
    and ``>`` folded blocks, ``[a, b]`` flow lists, quoted strings, booleans,
    integers and floats. Anything else raises rather than being guessed at, and
    the message says to install PyYAML.
    """
    try:
        import yaml  # noqa: WPS433 - optional, see the module docstring
    except ImportError:
        parsed = _mini_yaml(text, where)
    else:
        try:
            parsed = yaml.safe_load(text)
        except yaml.YAMLError as exc:
            raise MjolnirError("{0}: could not parse YAML: {1}".format(where, exc))
    if parsed is None:
        raise MjolnirError("{0} is empty".format(where))
    if not isinstance(parsed, dict):
        raise MjolnirError("{0}: a playbook must be a mapping".format(where))
    return parsed


def _strip_comment(line: str) -> str:
    out: List[str] = []
    quote = ""
    for index, char in enumerate(line):
        if quote:
            out.append(char)
            if char == quote:
                quote = ""
            continue
        if char in "\"'":
            quote = char
            out.append(char)
            continue
        if char == "#" and (index == 0 or line[index - 1] in " \t"):
            break
        out.append(char)
    return "".join(out).rstrip()


def _split_flow(text: str) -> List[str]:
    parts: List[str] = []
    current: List[str] = []
    quote = ""
    for char in text:
        if quote:
            current.append(char)
            if char == quote:
                quote = ""
            continue
        if char in "\"'":
            quote = char
            current.append(char)
            continue
        if char == ",":
            parts.append("".join(current))
            current = []
            continue
        current.append(char)
    if current:
        parts.append("".join(current))
    return [p.strip() for p in parts if p.strip()]


def _scalar(token: str) -> Any:
    text = token.strip()
    if not text:
        return ""
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
        return text[1:-1]
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        return [] if not inner else [_scalar(part) for part in _split_flow(inner)]
    lowered = text.lower()
    if lowered in ("true", "yes"):
        return True
    if lowered in ("false", "no"):
        return False
    if lowered in ("null", "~"):
        return None
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        pass
    return text


def _scan(text: str, where: str) -> List[Tuple[int, str, Optional[str]]]:
    """Significant lines as (indent, content, block-scalar-or-None)."""
    raw = text.splitlines()
    items: List[Tuple[int, str, Optional[str]]] = []
    index = 0
    while index < len(raw):
        line = raw[index]
        if not line.strip() or line.lstrip().startswith("#"):
            index += 1
            continue
        indent = len(line) - len(line.lstrip(" "))
        if "\t" in line[:indent]:
            raise MjolnirError(
                "{0} line {1}: tab used for indentation; YAML forbids it".format(
                    where, index + 1))
        content = _strip_comment(line.strip())
        block: Optional[str] = None
        if content.endswith(("|", ">")) and (
                content.endswith(": |") or content.endswith(": >")
                or content in ("- |", "- >")):
            style = content[-1]
            content = content[:-1].rstrip()
            index += 1
            body: List[str] = []
            while index < len(raw):
                nxt = raw[index]
                if not nxt.strip():
                    body.append("")
                    index += 1
                    continue
                if len(nxt) - len(nxt.lstrip(" ")) <= indent:
                    break
                body.append(nxt)
                index += 1
            while body and not body[-1].strip():
                body.pop()
            base = min((len(b) - len(b.lstrip(" ")) for b in body if b.strip()),
                       default=0)
            body = [b[base:] if b.strip() else "" for b in body]
            if style == ">":
                paragraphs: List[str] = []
                current: List[str] = []
                for entry in body:
                    if entry.strip():
                        current.append(entry.strip())
                    elif current:
                        paragraphs.append(" ".join(current))
                        current = []
                if current:
                    paragraphs.append(" ".join(current))
                block = "\n\n".join(paragraphs)
            else:
                block = "\n".join(body)
            # YAML's default "clip" chomping keeps exactly one trailing
            # newline. Matching it means this reader and PyYAML produce
            # byte-identical playbooks, which is what the parity test asserts.
            block += "\n"
            items.append((indent, content, block))
            continue
        items.append((indent, content, None))
        index += 1
    return items


def _mini_yaml(text: str, where: str) -> Any:
    items = _scan(text, where)
    if not items:
        return None
    value, position = _parse_node(items, 0, items[0][0], where)
    if position != len(items):
        raise MjolnirError(
            "{0}: could not parse the whole file (stopped at {1!r}); install "
            "PyYAML if the playbook uses YAML beyond the shipped "
            "subset".format(where, items[position][1][:60]))
    return value


def _parse_node(items: Sequence[Tuple[int, str, Optional[str]]], position: int,
                indent: int, where: str) -> Tuple[Any, int]:
    text = items[position][1]
    if text == "-" or text.startswith("- "):
        return _parse_sequence(items, position, indent, where)
    return _parse_mapping(items, position, indent, where)


def _parse_mapping(items: Sequence[Tuple[int, str, Optional[str]]], position: int,
                   indent: int, where: str) -> Tuple[Dict[str, Any], int]:
    out: Dict[str, Any] = {}
    while position < len(items):
        line_indent, content, block = items[position]
        if line_indent < indent:
            break
        if line_indent > indent:
            raise MjolnirError("{0}: unexpected indentation at {1!r}".format(
                where, content[:60]))
        if content == "-" or content.startswith("- "):
            break
        if ":" not in content:
            raise MjolnirError(
                "{0}: expected `key: value`, got {1!r}".format(where, content[:60]))
        key, _, rest = content.partition(":")
        key = key.strip()
        rest = rest.strip()
        position += 1
        if block is not None:
            out[key] = block
            continue
        if rest:
            out[key] = _scalar(rest)
            continue
        if position < len(items) and items[position][0] > line_indent:
            out[key], position = _parse_node(items, position, items[position][0], where)
        elif position < len(items) and items[position][0] == line_indent \
                and (items[position][1] == "-" or items[position][1].startswith("- ")):
            out[key], position = _parse_sequence(items, position, line_indent, where)
        else:
            out[key] = None
    return out, position


def _parse_sequence(items: Sequence[Tuple[int, str, Optional[str]]], position: int,
                    indent: int, where: str) -> Tuple[List[Any], int]:
    out: List[Any] = []
    while position < len(items):
        line_indent, content, block = items[position]
        if line_indent < indent:
            break
        if line_indent > indent:
            raise MjolnirError("{0}: unexpected indentation at {1!r}".format(
                where, content[:60]))
        if not (content == "-" or content.startswith("- ")):
            break
        rest = content[1:].strip()
        position += 1
        if block is not None:
            out.append(block)
            continue
        if not rest:
            if position < len(items) and items[position][0] > line_indent:
                value, position = _parse_node(items, position, items[position][0], where)
                out.append(value)
            else:
                out.append(None)
            continue
        if ":" in rest and rest[0] not in "\"'[":
            key, _, value_text = rest.partition(":")
            if not value_text.strip():
                raise MjolnirError(
                    "{0}: a list item whose first key has a nested value is "
                    "outside the built-in parser's subset ({1!r}); install "
                    "PyYAML".format(where, rest[:60]))
            entry: Dict[str, Any] = {key.strip(): _scalar(value_text)}
            if position < len(items) and items[position][0] > line_indent:
                more, position = _parse_mapping(items, position, items[position][0],
                                                where)
                entry.update(more)
            out.append(entry)
            continue
        out.append(_scalar(rest))
    return out, position


# ---------------------------------------------------------------------------
# Rule-derived prose — the anchor, and the substitute
# ---------------------------------------------------------------------------

def _check_dict(check: Check) -> Dict[str, Any]:
    return {
        "name": check.name,
        "category": check.category,
        "measured": check.measured,
        "status": check.status,
        "value": round_or_none(check.value, 4) if isinstance(check.value, float)
                 else check.value,
        "threshold": check.threshold,
        "comparison": check.comparison,
        "unit": check.unit,
        "source": check.source,
        "reading": check.reading,
    }


def _threshold_block(names: Sequence[str]) -> List[Dict[str, Any]]:
    """Registered thresholds with their sources, for the model to quote.

    Every number a Mjolnir report prints names where it came from, and the model
    is given the sources so it can attribute rather than assert. ``verified``
    False travels with the number: a citation written from memory and not
    checked against the primary document is worse than none, because it looks
    settled.
    """
    out: List[Dict[str, Any]] = []
    for name in names:
        entry = config.THRESHOLDS.get(name)
        if entry is None:
            continue
        out.append({"name": entry.name, "value": entry.value, "unit": entry.unit,
                    "source": entry.source, "citation_verified": entry.verified})
    return out


def _sample_threshold_names(result: SampleResult) -> List[str]:
    names = ["min_depth", "degraded_depth_floor", "min_breadth",
             "min_mapped_fraction", "major_variant_fraction",
             "min_unambiguous_fraction", "het_snp_fraction_warn",
             "het_snp_fraction_mixed", "f2_mixture_threshold",
             "f47_mixture_threshold", "min_barcode_support_fraction",
             "min_barcode_callable_fraction", "ani_species_floor"]
    if result.platform == PLATFORM_ONT:
        names += ["min_reads_ont", "ont_fbic_discordance_fraction",
                  "ont_indel_uncorroborated_fraction"]
    elif result.platform == PLATFORM_FASTA:
        names += ["fasta_capability_loss"]
    else:
        names += ["min_reads_illumina"]
    return names


def rule_summary(result: SampleResult) -> Tuple[str, str]:
    """The headline and body the rules alone support, as prose.

    This is the text the model is asked to improve on and the text that replaces
    its answer when the discipline rules discard it, so it has to be complete on
    its own: what the organism is, what was found, what was not measured, and
    what the platform could not have seen. Every number in it comes from a
    field of *result*, which is also why the model may safely quote it — the
    number check downstream is satisfied by the same observation this text
    travels in.
    """
    species = result.species
    lineage = result.lineage
    identity = species.display
    if not species.resolved_to_species and species.complex:
        identity = "{0} (not resolved below complex)".format(species.complex)
    head = ["{0}: {1}".format(result.sample_id, identity)]
    if lineage.lineage or lineage.sublineage:
        head.append("lineage {0}".format(lineage.display))
    if lineage.is_bcg:
        head.append("BCG")
    elif lineage.animal_variant:
        head.append(lineage.animal_variant)

    resistant = sorted(result.resistant_drugs(), key=lambda d: natural_key(d.drug))
    evaluated = [d for d in result.drugs if d.target_covered is not False]
    if resistant:
        head.append("resistance determinants detected for {0}".format(
            ", ".join(d.drug for d in resistant)))
    elif evaluated:
        head.append("{0} across the {1} drug{2} evaluated".format(
            NO_DETERMINANT_TEXT, len(evaluated), "" if len(evaluated) == 1 else "s"))
    else:
        head.append("no drug was evaluable in this sample")
    if result.contamination.verdict != VALIDITY_VALID:
        head.append("sample validity: {0}".format(result.contamination.verdict))
    headline = "; ".join(head) + "."

    body: List[str] = []
    if not species.resolved_to_species:
        body.append(
            "The species call did not resolve to species level ({0}); "
            "{1}.".format(species.method or "method not recorded",
                          "; ".join(species.caveats) if species.caveats
                          else "the evidence supports the complex only"))
    if lineage.barcode_sites_callable:
        body.append(
            "Lineage rests on {0} of {1} callable barcode sites ({2} in the "
            "scheme); support and callability are reported because a label "
            "alone hides how thin the evidence is.".format(
                lineage.barcode_sites_supporting, lineage.barcode_sites_callable,
                lineage.barcode_sites_total))

    if resistant:
        for call in resistant:
            detail = [call.label]
            if call.who_grade:
                detail.append("WHO grade {0}".format(call.who_grade))
            elif not call.who_graded:
                detail.append("not graded by the WHO catalogue")
            if call.supporting_variants:
                detail.append("from {0}".format(", ".join(call.supporting_variants)))
            body.append("{0}: {1}.".format(call.drug, "; ".join(detail)))
    uncertain = [d for d in result.drugs if d.call == CALL_UNCERTAIN]
    if uncertain:
        body.append(
            "Variants of uncertain significance were found for {0}; an "
            "uncertain grade is neither a resistance call nor a clearance."
            .format(", ".join(sorted(d.drug for d in uncertain))))
    no_call = [d for d in result.drugs
               if d.call == CALL_NO_CALL and d.target_covered is not False]
    if no_call:
        # The canonical wording, verbatim: the report, the TSV and this sentence
        # must say the same thing, and "no determinant" must never soften into
        # "susceptible" between them.
        body.append(
            "For {0}: {1}. That is an absence of a determinant and not a "
            "prediction of susceptibility; phenotypic testing remains the "
            "arbiter.".format(", ".join(sorted(d.drug for d in no_call)),
                              NO_DETERMINANT_TEXT))
    unevaluable = [d for d in result.drugs if d.target_covered is False]
    if unevaluable:
        body.append(
            "The target regions for {0} were not callable, so those drugs were "
            "not evaluated at all.".format(", ".join(sorted(d.drug for d in unevaluable))))
    suppressed = [d for d in result.drugs if d.suppressed_by]
    for call in suppressed:
        body.append("{0}: a call was suppressed by {1}.".format(call.drug,
                                                               call.suppressed_by))
    disagreeing = result.disagreements()
    if disagreeing:
        body.append(
            "Catalogues disagree for {0}; the annex prints all three side by "
            "side.".format(", ".join(sorted(d.drug for d in disagreeing))))

    failing = [c for c in result.all_checks() if c.status == STATUS_FAIL]
    warning = [c for c in result.all_checks()
               if c.status == STATUS_WARN and c.measured]
    if failing:
        body.append("Failed checks: {0}.".format(
            "; ".join(_render_check(c) for c in failing)))
    if warning:
        body.append("Checks in warning: {0}.".format(
            "; ".join(_render_check(c) for c in warning)))
    unmeasured = result.unmeasured()
    if unmeasured:
        body.append(
            "Not measured in this run: {0}. Nothing in that list is being "
            "reported as absent or normal.".format(", ".join(unmeasured)))

    contamination = result.contamination
    if contamination.verdict == VALIDITY_NOT_ASSESSED:
        body.append("Sample validity was not assessed.")
    else:
        body.append("Sample validity is {0}{1}.".format(
            contamination.verdict,
            ": " + contamination.verdict_reason if contamination.verdict_reason else ""))
    if not contamination.screen_informative and contamination.screen_note:
        body.append(contamination.screen_note)
    if contamination.mixture_class == MIXTURE_NOT_ASSESSED:
        body.append("Mixture status was not assessed.")

    for caveat in result.caveats:
        body.append(caveat if caveat.endswith(".") else caveat + ".")
    for warning_text in result.warnings:
        body.append("Warning: {0}.".format(warning_text.rstrip(".")))

    return headline, " ".join(body)


def _render_check(check: Check) -> str:
    if not check.measured:
        return "{0} (not measured)".format(check.name)
    value = check.value
    if isinstance(value, float):
        value = round(value, 4)
    if check.threshold is None:
        return "{0} = {1}{2}".format(check.name, value, check.unit)
    return "{0} = {1}{2} against {3} {4}".format(
        check.name, value, check.unit, check.comparison or "vs", check.threshold)


def rule_summary_cohort(cohort: CohortResult) -> Tuple[str, str]:
    """The same, for a cohort: clusters, the threshold, and its denominator."""
    clustered = [c for c in cohort.clusters if c.size > 1]
    head = ["{0} samples compared".format(len(cohort.samples))]
    if cohort.threshold is not None:
        head.append("clustering at {0} SNPs".format(cohort.threshold))
    head.append("{0} of more than one sample".format(plural(len(clustered), "cluster")))
    headline = "; ".join(head) + "."

    body: List[str] = []
    if cohort.threshold_basis:
        body.append("Threshold basis: {0}.".format(cohort.threshold_basis.rstrip(".")))
    if cohort.mask_name:
        body.append(
            "Distances were counted after masking with {0}{1}. Masking is not a "
            "solved constant and the mask used is named here for that "
            "reason.".format(
                cohort.mask_name,
                " ({0} sites, {1}% of the reference)".format(
                    cohort.masked_sites, percentage(cohort.masked_fraction))
                if cohort.masked_sites else ""))
    else:
        body.append("No mask was applied, so these distances include repetitive "
                    "and error-prone regions and are not comparable with masked "
                    "distances.")
    for cluster in clustered:
        body.append("{0}: {1} at a maximum distance of {2} SNPs over at least "
                    "{3} shared callable sites.".format(
                        cluster.cluster_id, ", ".join(cluster.members),
                        cluster.max_distance, cluster.min_shared_callable_sites))
    if not clustered:
        body.append("No pair fell within the threshold, so no cluster is claimed. "
                    "A distance above the threshold is not evidence against "
                    "epidemiological linkage on its own.")
    uncompared = [p for p in cohort.pairs if p.snps is None]
    if uncompared:
        body.append("{0} could not be compared and {1} absent from the matrix "
                    "rather than being scored as zero.".format(
                        plural(len(uncompared), "pair"),
                        "is" if len(uncompared) == 1 else "are"))
    for caveat in cohort.caveats:
        body.append(caveat if caveat.endswith(".") else caveat + ".")
    return headline, " ".join(body)


# ---------------------------------------------------------------------------
# Building the observations
# ---------------------------------------------------------------------------

def _drug_dict(call: Any) -> Dict[str, Any]:
    return {
        "drug": call.drug,
        "call": call.call,
        "label": call.label,
        "confidence": call.confidence,
        "who_graded": call.who_graded,
        "who_grade": call.who_grade,
        "level": call.level,
        "cross_resistance": list(call.cross_resistance),
        "disagreement": call.disagreement,
        "disagreement_kind": call.disagreement_kind,
        "suppressed_by": call.suppressed_by,
        "target_covered": call.target_covered,
        "supporting_variants": list(call.supporting_variants),
        "per_catalogue": [{"catalogue": c.catalogue, "call": c.call,
                           "grade": c.grade, "matched_by": c.matched_by}
                          for c in call.catalogue_calls],
        "caveats": list(call.caveats),
    }


def _variant_dict(variant: Any) -> Dict[str, Any]:
    return {
        "name": variant.display,
        "gene": variant.gene,
        "hgvs": variant.hgvs,
        "type": variant.variant_type,
        "effect": variant.effect,
        "depth": variant.depth,
        "allele_fraction": round_or_none(variant.allele_fraction, 4),
        "allele_percent": percentage(variant.allele_fraction, 1),
        "is_major": variant.is_major,
        "masked": variant.masked,
        "filters": list(variant.filters),
        "graded_by": sorted(set(c.catalogue for c in variant.catalogue_calls)),
    }


def build_sample_observation(result: SampleResult,
                             playbook: Optional[Playbook] = None,
                             run_config: Any = None,
                             limit: int = MAX_OBSERVATION_BYTES) -> Observation:
    """Everything the model may read about one sample, and nothing else.

    Percentages travel beside their fractions on purpose. The discipline layer
    rejects any number the input did not contain, and a reader writing "96%
    breadth" from an input that only said 0.96 would be computing rather than
    citing — correct arithmetic, but the same habit that produces a computed
    number nobody can trace.
    """
    qc = result.qc
    contamination = result.contamination
    if playbook is None:
        playbook = playbook_for(result.species)

    variants = sorted(result.variants,
                      key=lambda v: (not bool(v.catalogue_calls),
                                     natural_key(v.display)))
    variant_note = ""
    if len(variants) > MAX_VARIANT_ROWS:
        variant_note = ("{0} of {1} variants are listed, catalogue-graded ones "
                        "first; the rest are in the annex and their absence "
                        "here is truncation, not a measurement.".format(
                            MAX_VARIANT_ROWS, len(variants)))
        variants = variants[:MAX_VARIANT_ROWS]

    headline, body = rule_summary(result)
    data: Dict[str, Any] = {
        "sample": result.sample_id,
        "platform": result.platform,
        "profile": result.profile,
        "reference": result.reference,
        "mjolnir_version": result.mjolnir_version,
        "rule_derived_summary": {"headline": headline, "body": body},
        "species": {
            "display": result.species.display,
            "name": result.species.name,
            "complex": result.species.complex,
            "resolved_to_species": result.species.resolved_to_species,
            "method": result.species.method,
            "confidence": result.species.confidence,
            "ani_percent": round_or_none(result.species.ani, 3),
            "aligned_fraction": round_or_none(result.species.aligned_fraction, 4),
            "caveats": list(result.species.caveats),
        },
        "lineage": {
            "display": result.lineage.display,
            "lineage": result.lineage.lineage,
            "sublineage": result.lineage.sublineage,
            "scheme": result.lineage.scheme,
            "is_bcg": result.lineage.is_bcg,
            "animal_variant": result.lineage.animal_variant,
            "sites_supporting": result.lineage.barcode_sites_supporting,
            "sites_callable": result.lineage.barcode_sites_callable,
            "sites_in_scheme": result.lineage.barcode_sites_total,
            "support_fraction": round_or_none(result.lineage.support_fraction, 4),
            "mixed_lineages": list(result.lineage.mixed_lineages),
            "confidence": result.lineage.confidence,
            "caveats": list(result.lineage.caveats),
        },
        "drugs": [_drug_dict(d) for d in sorted(result.drugs,
                                                key=lambda d: natural_key(d.drug))],
        "variants": [_variant_dict(v) for v in variants],
        "variants_note": variant_note,
        "variant_count": len(result.variants),
        "qc": {
            "mean_depth": round_or_none(qc.mean_depth, 2),
            "median_depth": round_or_none(qc.median_depth, 2),
            "breadth_1x": round_or_none(qc.breadth_1x, 4),
            "breadth_10x": round_or_none(qc.breadth_10x, 4),
            "breadth_min_depth": round_or_none(qc.breadth_min_depth, 4),
            "breadth_min_depth_percent": percentage(qc.breadth_min_depth, 2),
            "coverage_evenness": round_or_none(qc.coverage_evenness, 4),
            "evenness_definition": qc.evenness_definition,
            "mapped_fraction": round_or_none(qc.mapped_fraction, 4),
            "mapped_percent": percentage(qc.mapped_fraction, 2),
            "gc_content": round_or_none(qc.gc_content, 4),
            "unambiguous_fraction": round_or_none(qc.unambiguous_fraction, 4),
            "total_reads": qc.total_reads,
            "mapped_reads": qc.mapped_reads,
            "mean_read_length": round_or_none(qc.mean_read_length, 1),
            "mean_base_quality": round_or_none(qc.mean_base_quality, 1),
            "reference_length": qc.reference_length,
            "checks": [_check_dict(c) for c in qc.checks],
        },
        "contamination": {
            "verdict": contamination.verdict,
            "verdict_reason": contamination.verdict_reason,
            "mixture_class": contamination.mixture_class,
            "f2": round_or_none(contamination.f2, 4),
            "f47": round_or_none(contamination.f47, 4),
            "lineage_het_sites": contamination.lineage_het_sites,
            "lineage_sites_examined": contamination.lineage_sites_examined,
            "het_snp_fraction": round_or_none(contamination.het_snp_fraction, 5),
            "het_snp_count": contamination.het_snp_count,
            "unambiguous_fraction": round_or_none(contamination.unambiguous_fraction, 4),
            "non_target_fraction": round_or_none(contamination.non_target_fraction, 4),
            "non_target_resolution": contamination.non_target_resolution,
            "screen_informative": contamination.screen_informative,
            "screen_method": contamination.screen_method,
            "screen_note": contamination.screen_note,
            "caveats": list(contamination.caveats),
            "checks": [_check_dict(c) for c in contamination.checks],
        },
        "checks": [_check_dict(c) for c in result.checks],
        "not_measured": result.unmeasured(),
        "caveats": list(result.caveats),
        "warnings": list(result.warnings),
        "status": result.status or result.overall_status(),
        "thresholds": _threshold_block(_sample_threshold_names(result)),
        "platform_caveats": list(config.platform_caveats(result.platform)),
        "playbook": playbook.to_dict(),
        "database_versions": [{"name": d.name, "version": d.version,
                               "checksum": d.checksum[:16]}
                              for d in result.database_versions],
    }
    if run_config is not None:
        overridden = getattr(run_config, "overridden_thresholds", None)
        if callable(overridden):
            changed = overridden()
            if changed:
                data["operator_overrides"] = to_jsonable(changed)

    observation = Observation(kind="sample", subject=result.sample_id,
                              data=to_jsonable(data))
    assert_no_sequence(observation.data)
    return observation.capped(limit)


def build_cohort_observation(cohort: CohortResult,
                             playbook: Optional[Playbook] = None,
                             limit: int = MAX_OBSERVATION_BYTES) -> Observation:
    """What the model may read about a cohort: distances with denominators."""
    headline, body = rule_summary_cohort(cohort)
    data: Dict[str, Any] = {
        "samples": list(cohort.samples),
        "reference": cohort.reference,
        "rule_derived_summary": {"headline": headline, "body": body},
        "threshold": cohort.threshold,
        "threshold_basis": cohort.threshold_basis,
        "mask": {"name": cohort.mask_name, "masked_sites": cohort.masked_sites,
                 "masked_fraction": round_or_none(cohort.masked_fraction, 4),
                 "masked_percent": percentage(cohort.masked_fraction, 2)},
        "joint_sites": cohort.joint_sites,
        "pairs": [{"a": p.sample_a, "b": p.sample_b, "snps": p.snps,
                   "shared_callable_sites": p.shared_callable_sites,
                   "snps_per_mb": round_or_none(p.snps_per_mb, 3),
                   "note": p.note}
                  for p in cohort.pairs],
        "clusters": [{"id": c.cluster_id, "members": list(c.members),
                      "max_distance": c.max_distance,
                      "min_shared_callable_sites": c.min_shared_callable_sites,
                      "note": c.note}
                     for c in cohort.clusters],
        "checks": [_check_dict(c) for c in cohort.checks],
        "caveats": list(cohort.caveats),
        "warnings": list(cohort.warnings),
        "thresholds": _threshold_block(
            ["cluster_snp_strict", "cluster_snp_relaxed", "chimaera_local_distance",
             "default_cluster_distance", "snp_proximity_window",
             "min_shared_callable_sites", "masked_loci_h37rv",
             "masked_fraction_h37rv"]),
    }
    if playbook is not None:
        data["playbook"] = playbook.to_dict()
    observation = Observation(kind="cohort",
                              subject=", ".join(cohort.samples[:8]),
                              data=to_jsonable(data))
    assert_no_sequence(observation.data)
    return observation.capped(limit)


# ---------------------------------------------------------------------------
# Prompts
# ---------------------------------------------------------------------------

_BASE_RULES = """\
Rules you must follow. An answer that breaks one of them is discarded by a \
checker and replaced with the rule-derived summary, so breaking one costs the \
reader your whole answer.
- Use only the evidence given below. Never state a number that is not in it, \
and do not compute new numbers from the ones that are: cite them.
- Never describe anything in `not_measured` as absent, normal, clean or fine. \
Say it was not measured.
- Never write that an organism or a sample is susceptible to a drug. Where no \
determinant was found, the finding is "no resistance determinant detected", \
which is an absence of evidence and not evidence of susceptibility.
- Never contradict a call in `drugs` or the verdict in `contamination`. You may \
qualify one; you may not reverse it.
- Attribute thresholds to the source given with them.
- Write for the audience named in the playbook. Be brief: a headline of one \
sentence and a body of three to six.

Reply with JSON only: {"headline": "<one sentence>", "body": "<three to six \
sentences>"}"""


def system_prompt(playbook: Optional[Playbook] = None) -> str:
    """The instructions, with the playbook's organism-specific ones folded in."""
    lines = ["You are writing the interpretation for a mycobacterial genome "
             "report. A deterministic pipeline has already reached every "
             "verdict; your job is to read them back to a clinician, not to "
             "re-decide them."]
    if playbook is not None:
        if playbook.audience:
            lines.append("Audience: {0}".format(playbook.audience))
        if playbook.describe:
            lines.append(playbook.describe)
        for statement in playbook.must_state:
            lines.append("- You must state: {0}".format(statement))
        for statement in playbook.must_not_say:
            lines.append("- You must not say: {0}".format(statement))
    lines.append(_BASE_RULES)
    return "\n".join(lines)


def reading_prompt(observation: Observation) -> str:
    """The user turn: the evidence, and the summary to write over."""
    return ("Evidence for {0} (JSON):\n{1}\n\n"
            "The rule-derived summary below is what the report prints if your "
            "answer is discarded. Improve its readability without changing any "
            "of its claims.\n{2}\n").format(
                observation.subject, observation.json(indent=1),
                json.dumps(observation.data.get("rule_derived_summary", {}), indent=1))


def gate_prompt(observation: Observation, gate: Gate) -> str:
    """The user turn for a closed-set decision."""
    return ("Evidence for {0} (JSON):\n{1}\n\n"
            "Decision point: {2}\n{3}\n\nOptions: {4}\n"
            "Reply with JSON only: {{\"choice\": \"<one of the options>\", "
            "\"reason\": \"<why, citing the evidence>\"}}").format(
                observation.subject, observation.json(indent=1), gate.id,
                gate.question, json.dumps(gate.options))


def context_text(observation: Observation) -> str:
    """The string the discipline rules check the answer's numbers against.

    The whole observation, flattened. Anything the model may legitimately quote
    is in here by construction, which is what makes "states a number absent from
    its input" a decidable question rather than a judgement call.
    """
    return observation.json(indent=1)
