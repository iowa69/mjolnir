"""Resistance calling: three catalogues, one normal form, one consensus.

The subpackage is deliberately layered, and the layering is what keeps the
consensus rule auditable:

* :mod:`~mjolnir.resistance.normalise` turns an observed variant into the two
  keys the design requires — a left-aligned, parsimonious coordinate key for
  WHO's own matching protocol, and a three-letter HGVS key for the
  cross-catalogue join — plus the legacy-numbering alias table that stops a
  codon-numbering difference being reported as a biological disagreement.
* :mod:`~mjolnir.resistance.catalogues` loads WHO v2, MTBseq's flat list and
  tbdb into :class:`~mjolnir.resistance.catalogues.CatalogueEntry` rows, keyed
  on the (drug, variant) **pair**, and refuses inputs that are not what they
  claim to be.
* ``rules``, ``consensus`` and ``ntm`` sit above both and are imported directly
  by their consumers.

Only the two bottom layers are re-exported here. The higher modules import this
one, so pulling them into the package's ``__init__`` would make every import of
``mjolnir.resistance.normalise`` drag the whole consensus engine in behind it —
and would turn a syntax error in one of them into an unexplainable failure in
the other. ``import mjolnir.resistance.consensus`` still works exactly as
expected.

The bare ``normalise()`` function is also left out on purpose: re-exporting it
here would shadow the submodule of the same name, so ``from mjolnir.resistance
import normalise`` would hand back a function where every reader expects a
module. Its typed wrappers — :func:`normalise_variant`,
:func:`normalise_coordinate` — are exported instead.
"""

from __future__ import annotations

from .catalogues import (
    Catalogue,
    CatalogueEntry,
    calls_for_variant,
    database_versions,
    load_catalogues,
    load_mtbseq,
    load_tbdb,
    load_who,
    read_who_coordinates_file,
    refuse_who_text_master,
)
from .normalise import (
    CoordinateKey,
    NormalisedVariant,
    NumberingAlias,
    alias_keys,
    classify_difference,
    coordinate_string,
    hgvs_key,
    is_rule_variant,
    normalise_coordinate,
    normalise_hgvs,
    normalise_variant,
    normalise_variants,
    split_key,
    three_letter,
)

__all__ = [
    "Catalogue",
    "CatalogueEntry",
    "CoordinateKey",
    "NormalisedVariant",
    "NumberingAlias",
    "alias_keys",
    "calls_for_variant",
    "classify_difference",
    "coordinate_string",
    "database_versions",
    "hgvs_key",
    "is_rule_variant",
    "load_catalogues",
    "load_mtbseq",
    "load_tbdb",
    "load_who",
    "normalise_coordinate",
    "normalise_hgvs",
    "normalise_variant",
    "normalise_variants",
    "read_who_coordinates_file",
    "refuse_who_text_master",
    "split_key",
    "three_letter",
]
