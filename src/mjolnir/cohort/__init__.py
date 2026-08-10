"""Cohort mode: joint variant table, masked SNP distances, clusters.

Three stages, in this order, and nothing here runs an external tool — it is
arithmetic over records the single-sample pipeline already produced::

    from mjolnir.cohort import (SampleVariants, Regions, build_joint_table,
                                load_mask, distance_matrix, cluster_samples)

    table = build_joint_table([
        SampleVariants(sample_id="A", variants=variants_a,
                       callable_regions=callable_a, reference=ref),
        SampleVariants(sample_id="B", variants=variants_b,
                       callable_regions=callable_b, reference=ref),
    ])
    mask = load_mask(db_dir / "tbdb" / "mask.bed", name="tbdb mask.bed 2025-08")
    matrix = distance_matrix(table, mask)
    clusters = cluster_samples(table.samples, matrix.pairs, threshold=12)

The output maps onto :class:`~mjolnir.records.CohortResult` field by field:
``pairs`` from ``matrix.pairs``, ``clusters`` from ``clusters.clusters``,
``threshold`` and ``threshold_basis`` from the assignment, ``mask_name`` /
``masked_sites`` / ``masked_fraction`` from the :class:`~.distance.Mask`,
``joint_sites`` from ``table.site_count``, and ``checks`` / ``caveats`` from the
concatenation of all three stages' own.

Two properties hold across the whole subpackage. A distance never leaves it
without its shared-callable denominator attached — every function returns
:class:`~mjolnir.records.PairwiseDistance` records rather than integers. And a
cohort of one or two samples is an ordinary input that produces a table and a
stated consequence, not the abort MTBseq's ``TBjoin`` produced on the local
*M. chimaera* run.
"""

from __future__ import annotations

from .cluster import (CLUSTER_ID_PREFIX, ClusterAssignment, cluster_samples,
                      clusters_at)
from .distance import (MASK_ABSENT_TEXT, DistanceMatrix, Mask, distance_matrix,
                       format_distance, iter_comparable, load_mask,
                       masked_fraction, pairs_for_cohort, pairwise_distance)
from .joint import (AMBIGUOUS, UNKNOWN_SYMBOL, JointSite, JointTable, Regions,
                    SampleVariants, build_joint_table, callable_summary,
                    cohort_size_check, merged_positions)

__all__ = [
    "AMBIGUOUS",
    "CLUSTER_ID_PREFIX",
    "ClusterAssignment",
    "DistanceMatrix",
    "JointSite",
    "JointTable",
    "MASK_ABSENT_TEXT",
    "Mask",
    "Regions",
    "SampleVariants",
    "UNKNOWN_SYMBOL",
    "build_joint_table",
    "callable_summary",
    "cluster_samples",
    "clusters_at",
    "cohort_size_check",
    "distance_matrix",
    "format_distance",
    "iter_comparable",
    "load_mask",
    "masked_fraction",
    "merged_positions",
    "pairs_for_cohort",
    "pairwise_distance",
]
