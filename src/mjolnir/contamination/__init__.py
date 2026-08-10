"""Contamination: what can honestly be measured, and what is refused outright.

Design §8. The refusals are as much of this subpackage's job as the
measurements, so they are enforced by the shape of the API rather than
documented and hoped for:

- No MTBC member is ever reported from a taxonomic classifier. Every taxon label
  leaves through :func:`~mjolnir.contamination.purity.taxon_label_for_report`,
  which collapses the complex.
- A Kraken2 run against a standard, PlusPF or capped index is not a
  contamination screen for mycobacteria.
  :class:`~mjolnir.contamination.purity.TaxonomicScreen` carries
  ``informative``, its rows are unreachable when it is False, and its check is
  an unmeasured one — so the only thing such a screen can contribute is the
  sentence saying it could not be run meaningfully.
- Kraken2's ``--confidence`` default of 0.0 is refused in ``config.py``.
- CheckM, CheckM2 and ConFindr raise from
  :func:`~mjolnir.contamination.purity.assert_mixture_method_supported` rather
  than returning a plausible-looking number.

What is measured: F2/F47 across lineage-defining positions, the genome-wide
heterozygous-SNP fraction under the MixInfect filters as a two-tier class, the
unambiguous-base fraction MTBseq would have discarded, and the read-level
composition signals — mapped fraction, breadth, evenness, GC, and the non-target
fraction by ANI when a reference set exists.

The headline is a sample-validity verdict per intended use, never a purity
percentage: 99.84% pure still produced 13 false-positive SNPs, and resistance
calling tolerates that where an outbreak SNP distance does not.
"""

from __future__ import annotations

from .heterozygosity import (
    F2_TOP_SITES,
    F47_TOP_SITES,
    F_STATISTIC_DEFINITION,
    HeterozygosityResult,
    LineageMixtureStats,
    LineageSite,
    MIXINFECT_FILTER_TEXT,
    ONT_SINGLE_STRAIN_REFUSAL,
    SiteObservation,
    UNAMBIGUOUS_SURFACED_TEXT,
    assess_heterozygosity,
    callable_sites,
    classify_mixture,
    f_statistic,
    genome_wide_heterozygosity,
    heterozygosity_checks,
    is_heterozygous,
    lineage_mixture_statistics,
    unambiguous_fraction,
)
from .purity import (
    INTENDED_USES,
    MAX_NON_TARGET_FRACTION,
    MTBC_CLASSIFIER_LABEL,
    REFUSED_MIXTURE_METHODS,
    SCREEN_INFORMATIVE,
    SCREEN_UNINFORMATIVE,
    SCREEN_UNINFORMATIVE_HEADLINE,
    USE_RESISTANCE,
    USE_TRANSMISSION,
    AniAssignment,
    NonTargetResult,
    PurityPanel,
    TaxonomicScreen,
    ValidityVerdict,
    assert_mixture_method_supported,
    assess_contamination,
    assess_non_target,
    evaluate_kraken2_screen,
    is_mtbc_member_name,
    measure_purity,
    no_screen,
    parse_kraken2_report,
    sample_validity,
    taxon_label_for_report,
    worst_validity,
)

__all__ = [
    # entry point
    "assess_contamination",
    # heterozygosity
    "SiteObservation", "LineageSite", "HeterozygosityResult", "LineageMixtureStats",
    "assess_heterozygosity", "callable_sites", "classify_mixture", "f_statistic",
    "genome_wide_heterozygosity", "heterozygosity_checks", "is_heterozygous",
    "lineage_mixture_statistics", "unambiguous_fraction",
    "F2_TOP_SITES", "F47_TOP_SITES", "F_STATISTIC_DEFINITION",
    "MIXINFECT_FILTER_TEXT", "UNAMBIGUOUS_SURFACED_TEXT", "ONT_SINGLE_STRAIN_REFUSAL",
    # purity, screening and refusals
    "PurityPanel", "measure_purity",
    "AniAssignment", "NonTargetResult", "assess_non_target", "MAX_NON_TARGET_FRACTION",
    "TaxonomicScreen", "evaluate_kraken2_screen", "no_screen", "parse_kraken2_report",
    "SCREEN_INFORMATIVE", "SCREEN_UNINFORMATIVE", "SCREEN_UNINFORMATIVE_HEADLINE",
    "taxon_label_for_report", "is_mtbc_member_name", "MTBC_CLASSIFIER_LABEL",
    "assert_mixture_method_supported", "REFUSED_MIXTURE_METHODS",
    # verdict
    "ValidityVerdict", "sample_validity", "worst_validity",
    "INTENDED_USES", "USE_RESISTANCE", "USE_TRANSMISSION",
]
