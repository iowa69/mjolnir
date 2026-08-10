"""Species and lineage typing (design §6).

Two modules with one rule between them: each answers only the question its
evidence can carry. ``species.py`` names a species from ANI and stops at the
complex for MTBC and, absent a marker panel, for MAC. ``lineage.py`` names the
MTBC member from lineage-defining SNPs read out of a pileup, and reports the
barcode support behind every label.

Neither will accept a taxonomic read classifier as evidence of species. The MTBC
members are not at species rank in current NCBI taxonomy, so a classifier row
naming one of them is not an identification;
:func:`~mjolnir.typing.species.species_from_classifier` raises rather than
returning, and :func:`~mjolnir.typing.species.demote_classifier_label` is what
classifier output may legitimately become.
"""

from __future__ import annotations

from .lineage import (
    BarcodeSite,
    PileupCounts,
    SiteGenotype,
    TaxonSupport,
    barcode_path,
    barcode_positions,
    barcode_taxa,
    call_lineage,
    describe_call,
    family_for,
    genotype_barcode_sites,
    lineage_checks,
    lineage_not_applicable,
    load_barcode,
    major_allele,
    scheme_animal_coverage,
    scheme_description,
    summarise_taxa,
    taxon_ancestors,
    taxon_root,
)
from .species import (
    AniMatch,
    MarkerResult,
    MarkerSnp,
    ReferenceGenome,
    ani_matches,
    complex_for,
    demote_classifier_label,
    describe_reference_set,
    genotype_markers,
    identify_species,
    is_mtbc_member,
    load_mac_markers,
    load_reference_set,
    run_mash,
    run_skani,
    species_checks,
    species_from_classifier,
)

__all__ = [
    "AniMatch",
    "BarcodeSite",
    "MarkerResult",
    "MarkerSnp",
    "PileupCounts",
    "ReferenceGenome",
    "SiteGenotype",
    "TaxonSupport",
    "ani_matches",
    "barcode_path",
    "barcode_positions",
    "barcode_taxa",
    "call_lineage",
    "complex_for",
    "demote_classifier_label",
    "describe_call",
    "describe_reference_set",
    "family_for",
    "genotype_barcode_sites",
    "genotype_markers",
    "identify_species",
    "is_mtbc_member",
    "lineage_checks",
    "lineage_not_applicable",
    "load_barcode",
    "load_mac_markers",
    "load_reference_set",
    "major_allele",
    "run_mash",
    "run_skani",
    "scheme_animal_coverage",
    "scheme_description",
    "species_checks",
    "species_from_classifier",
    "summarise_taxa",
    "taxon_ancestors",
    "taxon_root",
]
