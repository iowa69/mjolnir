"""The local-LLM interpretation layer: a client, an observation, and a checker.

The order of those three is the design. A deterministic pipeline reaches every
verdict first; :mod:`~mjolnir.agent.observation` builds the finished checks,
thresholds, sources, gene names and drug calls the model is allowed to read, and
the rule-derived summary it is asked to write over;
:mod:`~mjolnir.agent.client` speaks to whatever local server the operator runs;
:mod:`~mjolnir.agent.discipline` reads the answer a sentence at a time and
decides whether it may be printed.

Two outcomes are normal and neither is an error. If no host is reachable the
report is rule-only and says so. If the answer breaks a discipline rule it is
discarded, the rule-derived summary is printed instead, and
``Interpretation.discarded_reason`` names the rule — the report never silently
loses the model and never silently keeps a bad answer.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple

from mjolnir.agent.client import (LLMClient, LLMResult, ScriptedClient,
                                  client_from_config)
from mjolnir.agent.discipline import (RuleFacts, Verdict, facts_from_cohort,
                                      facts_from_sample, review)
from mjolnir.agent.observation import (Gate, Observation, Playbook, SequenceLeak,
                                       assert_no_sequence,
                                       available_playbooks,
                                       build_cohort_observation,
                                       build_sample_observation, context_text,
                                       gate_prompt, load_playbook,
                                       playbook_for, playbook_named,
                                       reading_prompt, rule_summary,
                                       rule_summary_cohort, system_prompt)
from mjolnir.records import CohortResult, Interpretation, SampleResult
from mjolnir.utils import LOG

__all__ = [
    "LLMClient", "LLMResult", "ScriptedClient", "client_from_config",
    "Observation", "Playbook", "Gate", "SequenceLeak", "assert_no_sequence",
    "build_sample_observation", "build_cohort_observation", "context_text",
    "system_prompt", "reading_prompt", "gate_prompt", "rule_summary",
    "rule_summary_cohort", "load_playbook", "playbook_for", "playbook_named",
    "available_playbooks", "RuleFacts", "Verdict", "review",
    "facts_from_sample", "facts_from_cohort", "interpret_sample",
    "interpret_cohort", "decide_gate",
]

#: How many times a rejected answer is handed back with the rule it broke. One
#: retry, because a model that invented a number tends to invent a different one
#: rather than to stop, and the rule-derived summary is a perfectly good report.
MAX_ATTEMPTS = 2


def _rule_only(headline: str, body: str, reason: str, playbook: str,
               client: Optional[LLMClient] = None) -> Interpretation:
    return Interpretation(
        headline=headline, body=body, rule_only=True, discarded_reason=reason,
        model=client.model if client is not None else "",
        host=client.host if client is not None else "",
        playbook=playbook)


def _interpret(observation: Observation, facts: RuleFacts, playbook: Playbook,
               client: Optional[LLMClient], headline: str, body: str
               ) -> Interpretation:
    """Shared body of the sample and cohort paths.

    Never raises on the model's account. A sequence leak raises — that is a bug
    in whichever module filled the observation, and it must stop the run — but a
    host that is down, an empty answer or a rejected one all degrade to the
    rule-derived summary with the reason recorded.
    """
    if client is None:
        return _rule_only(headline, body,
                          "no model was used; the report is rule-derived",
                          playbook.name)

    context = context_text(observation)
    system = system_prompt(playbook)
    prompt = reading_prompt(observation)
    last_reason = ""
    for attempt in range(MAX_ATTEMPTS):
        result = client.complete(prompt, system=system, json_mode=True)
        if not result.ok:
            return _rule_only(headline, body, result.error, playbook.name, client)
        verdict = review(result.text, context, facts)
        if verdict.ok:
            parsed = result.as_json() or {}
            model_headline = str(parsed.get("headline") or "").strip()
            model_body = str(parsed.get("body") or "").strip()
            if not model_headline and not model_body:
                last_reason = ("the model did not return the requested "
                               "headline/body JSON")
            else:
                return Interpretation(
                    headline=model_headline or headline,
                    body=model_body or body,
                    rule_only=False, discarded_reason="",
                    model=result.model, host=result.host, playbook=playbook.name)
        else:
            last_reason = verdict.reason
            LOG.debug("model answer rejected (%s): %s", verdict.rule,
                      verdict.sentence or verdict.reason)
        if attempt + 1 < MAX_ATTEMPTS:
            prompt = ("{0}\n\nYour previous answer was rejected because it "
                      "{1}. Write it again without that.".format(
                          reading_prompt(observation), last_reason))
    return _rule_only(headline, body,
                      "the model's answer was discarded: {0}".format(last_reason),
                      playbook.name, client)


def interpret_sample(result: SampleResult, client: Optional[LLMClient] = None,
                     playbook: Optional[Playbook] = None,
                     run_config: Any = None) -> Interpretation:
    """The prose for one sample: the model's, or the rules', and which one it is.

    The caller assigns the return value to ``SampleResult.interpretation`` and
    the report prints ``rule_only`` and ``discarded_reason`` on the page. Both
    fields are the point — a reader must be able to tell reasoning from
    reassurance without leaving the document.
    """
    if playbook is None:
        playbook = playbook_for(result.species)
    observation = build_sample_observation(result, playbook=playbook,
                                           run_config=run_config)
    headline, body = rule_summary(result)
    return _interpret(observation, facts_from_sample(result), playbook, client,
                      headline, body)


def interpret_cohort(cohort: CohortResult, client: Optional[LLMClient] = None,
                     playbook: Optional[Playbook] = None) -> Interpretation:
    """The prose for a cohort: clusters, the threshold, and its denominator.

    The playbook is optional and defaults to a neutral one rather than to MTBC.
    A cohort may mix organisms, and quietly handing the model the tuberculosis
    reading for a set of *M. chimaera* isolates would put WHO grades and lineage
    barcodes in front of it that do not apply to a single sample in the run.
    """
    if playbook is None:
        playbook = Playbook(name="cohort", organism="as submitted",
                            audience="whoever ordered the comparison",
                            describe="Pairwise SNP distances and clusters. Say "
                                     "which mask and which threshold produced "
                                     "them, and give the shared callable-site "
                                     "denominator beside every distance.")
    observation = build_cohort_observation(cohort, playbook=playbook)
    headline, body = rule_summary_cohort(cohort)
    return _interpret(observation, facts_from_cohort(cohort), playbook, client,
                      headline, body)


def decide_gate(gate: Gate, observation: Observation,
                client: Optional[LLMClient] = None,
                facts: Optional[RuleFacts] = None) -> Tuple[str, str, str]:
    """Choose one option from a gate's closed set: ``(choice, source, reason)``.

    ``source`` is ``"llm"`` or ``"default"`` and is never fudged. Taking the
    declared default while recording it as a model decision would let a run look
    as though it reasoned when it did not, and the gate defaults exist precisely
    for the runs where nothing reasoned.
    """
    if not gate.options:
        raise ValueError("gate {0!r} has no options".format(gate.id))
    if client is None:
        return gate.default, "default", "no model was used; took the declared default"

    prompt = gate_prompt(observation, gate)
    context = context_text(observation)
    last = ""
    for _attempt in range(MAX_ATTEMPTS):
        result = client.complete(prompt, system=system_prompt(None), json_mode=True)
        if not result.ok:
            return gate.default, "default", result.error
        verdict = review(result.text, context, facts, options=gate.options)
        if verdict.ok:
            parsed = result.as_json() or {}
            return (str(parsed.get("choice")), "llm",
                    str(parsed.get("reason") or "")[:300])
        last = verdict.reason
    return gate.default, "default", "the model's answer was discarded: {0}".format(last)
