"""What the model is allowed to see, and the leak that must be impossible.

House rule 3 says the model never sees raw sequence, and the design requires it
to be enforced in code rather than in a prompt: an observation containing a long
nucleotide run raises before it can reach the model. The check has to walk keys
as well as values, because the field that leaks is usually one somebody added
downstream — ``variants[3].note`` holding a consensus stretch — and not a field
anybody planned.

The second thing tested here is that a shrunken observation announces itself. A
model told it has the whole variant table when it has half of it writes "no other
catalogued variants were found", which is a claim nobody measured.
"""

from __future__ import annotations

import pytest

from mjolnir.agent import observation as obs
from mjolnir.records import (CALL_NO_CALL, CALL_R, ContaminationResult, DrugCall,
                             QCMetrics, SampleResult, SpeciesCall)


@pytest.fixture
def result():
    sample = SampleResult(sample_id="226-18", platform="illumina",
                          reference="NC_000962.3",
                          species=SpeciesCall(complex="MTBC", confidence="moderate"),
                          qc=QCMetrics(mean_depth=64.0, breadth_min_depth=0.99),
                          contamination=ContaminationResult(verdict="valid"))
    sample.drugs = [DrugCall(drug="Rifampicin", call=CALL_R, confidence="high"),
                    DrugCall(drug="Isoniazid", call=CALL_NO_CALL)]
    return sample


# ------------------------------------------------------------------- the leak

def test_a_gene_length_nucleotide_run_raises():
    with pytest.raises(obs.SequenceLeak) as excinfo:
        obs.assert_no_sequence({"note": "ACGT" * 60})
    assert "nucleotides" in str(excinfo.value)


def test_the_leak_names_the_path_it_was_found_at():
    payload = {"variants": [{"hgvs": "p.Ser450Leu"}, {"note": "GATTACA" * 40}]}
    with pytest.raises(obs.SequenceLeak) as excinfo:
        obs.assert_no_sequence(payload)
    assert "variants[1].note" in str(excinfo.value)


@pytest.mark.parametrize("key", sorted(obs.FORBIDDEN_KEYS))
def test_a_field_named_as_sequence_raises_even_when_empty(key):
    """A short sequence field is a leak that has not grown yet."""
    with pytest.raises(obs.SequenceLeak):
        obs.assert_no_sequence({key: ""})


def test_hgvs_and_gene_names_are_not_mistaken_for_sequence():
    obs.assert_no_sequence({
        "gene": "rpoB", "hgvs": "p.Ser450Leu", "variant": "inhA_c.-154G>A",
        "coordinate": "NC_000962.3:761155C>T", "lineage": "lineage4.9",
        "total_reads": 2_400_000, "mapped_reads": 2_350_000,
    })


def test_the_check_runs_before_an_observation_is_capped():
    leaky = obs.Observation(kind="sample", subject="226-18",
                            data={"consensus": "ACGTACGTAC" * 30})
    with pytest.raises(obs.SequenceLeak):
        leaky.capped()


# --------------------------------------------------------- truncation announces itself

def test_a_truncated_observation_says_absence_is_truncation(result):
    observation = obs.build_sample_observation(result)
    small = obs.Observation(kind=observation.kind, subject=observation.subject,
                            data={"rows": [{"n": i, "pad": "x" * 50}
                                           for i in range(400)]}).capped(limit=2000)
    assert small.truncated is True
    assert "not a measurement" in small.note
    assert small.nbytes() <= 4000


def test_an_observation_inside_the_limit_is_not_marked_truncated(result):
    observation = obs.build_sample_observation(result)
    assert observation.truncated is False


# ------------------------------------------------------ what the observation carries

def test_the_observation_carries_the_rule_derived_summary(result):
    observation = obs.build_sample_observation(result)
    summary = observation.data["rule_derived_summary"]
    assert summary["headline"], "the model writes over a verdict, not instead of one"


def test_the_observation_repeats_the_no_determinant_wording(result):
    observation = obs.build_sample_observation(result)
    assert "no resistance determinant detected" in obs.context_text(observation).lower()
    # The word "susceptible" may appear only where the playbook forbids saying
    # it — never as the description of a drug call.
    for entry in observation.data["drugs"]:
        rendered = " ".join(str(v) for v in entry.values()).lower()
        assert "susceptible" not in rendered, entry


def test_percentages_travel_beside_their_fractions(result):
    """So a reader can cite a percentage instead of computing one."""
    text = obs.context_text(obs.build_sample_observation(result))
    assert "0.99" in text


def test_the_context_is_what_the_discipline_rules_check_numbers_against(result):
    observation = obs.build_sample_observation(result)
    from mjolnir.agent import discipline

    context = obs.context_text(observation)
    facts = discipline.facts_from_sample(result)
    assert discipline.review("Mean depth was 64.0x.", context, facts)
    assert not discipline.review("Mean depth was 12.7x.", context, facts)


# ------------------------------------------------------------------- playbooks

def test_both_shipped_playbooks_load():
    names = obs.available_playbooks()
    assert set(names) >= {"mtbc", "ntm"}
    for name in ("mtbc", "ntm"):
        playbook = obs.playbook_named(name)
        assert playbook.organism
        assert playbook.must_not_say, "each playbook forbids something explicitly"


def test_an_mtbc_sample_gets_the_mtbc_playbook():
    playbook = obs.playbook_for(SpeciesCall(complex="MTBC"))
    assert "mtbc" in playbook.name.lower() or "tuberculosis" in playbook.organism.lower()


def test_an_ntm_sample_does_not_get_the_mtbc_playbook():
    playbook = obs.playbook_for(SpeciesCall(name="Mycobacterium chimaera",
                                            complex="MAC",
                                            resolved_to_species=True))
    assert "mtbc" not in playbook.name.lower()


def test_a_gate_declares_a_closed_option_set_and_a_default():
    for name in ("mtbc", "ntm"):
        for gate in obs.playbook_named(name).gates:
            assert gate.options, "a gate with no options is not a closed set"
            assert gate.default in gate.options, \
                "the default is what a run without a model takes"


def test_a_playbook_with_a_gate_default_outside_its_options_is_refused():
    raw = {
        "organism": "test",
        "gates": [{"id": "g1", "question": "?", "options": ["a", "b"],
                   "default": "c"}],
    }
    from mjolnir.utils import MjolnirError

    with pytest.raises(MjolnirError):
        obs.playbook_from_dict(raw)
