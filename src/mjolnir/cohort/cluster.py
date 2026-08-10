"""Single-linkage clustering at a stated threshold, with the basis attached.

A cluster is a claim about transmission, and the threshold is most of the claim.
So the number is a parameter and never a constant: the TB conventions are 5 SNPs
for recent transmission and 12 for a wider epidemiological link, while the prior
*M. chimaera* run on this machine used 6 — and a report that prints "cluster"
without printing which of those was applied has withheld the part a reader needs.
:func:`cluster_samples` therefore carries ``threshold_basis`` through to the
output, filled from :func:`mjolnir.config.cluster_threshold_basis`, which says in
words where the number came from and says plainly when it came from nowhere.

Two exclusions are structural rather than incidental.

**A pair with no distance does not link.** ``snps is None`` means the pair had no
shared callable denominator, which is an absence of comparison. It is not a
distance of zero and it is not a large distance, so it forms no edge, and the
count of such pairs is reported beside the clusters.

**A pair over too little shared sequence does not link either.** 12 differences
over 400 kb is not the 12 the published thresholds are about. Pairs below
:data:`~mjolnir.config.MIN_SHARED_CALLABLE_SITES` are excluded from the graph and
listed, so a cluster is never held together by an edge whose denominator would
have embarrassed it.

Single linkage is what the SNP-threshold literature uses, and it chains: two
isolates 20 SNPs apart sit in one cluster if something lies 10 SNPs from each. It
is used here for comparability, and the chaining is measured rather than hidden —
every cluster reports the largest distance inside it, and a cluster whose
internal maximum exceeds the threshold says so in its note.

Pure functions over pair records. Nothing here reads a file or runs a tool.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ..config import (CLUSTER_SNP_RELAXED, CLUSTER_SNP_STRICT,
                      CLUSTER_THRESHOLD_BASIS, DEFAULT_CLUSTER_DISTANCE,
                      MIN_SHARED_CALLABLE_SITES, cluster_threshold_basis,
                      source_for)
from ..records import (Check, Cluster, PairwiseDistance, STATUS_PASS,
                       STATUS_WARN, pair_key)
from ..utils import LOG, MjolnirError, natural_key

#: Prefix for generated cluster ids. Short, because it is printed in every row
#: of the cohort table.
CLUSTER_ID_PREFIX = "C"


@dataclass
class ClusterAssignment:
    """Clusters, the threshold that made them, and everything left out.

    The excluded lists are part of the result, not diagnostics: a reader who is
    told "3 clusters" and not told that 40 pairs had no denominator has been
    given a number that means less than it appears to.
    """

    clusters: List[Cluster] = field(default_factory=list)
    singletons: List[str] = field(default_factory=list)
    threshold: int = DEFAULT_CLUSTER_DISTANCE
    threshold_basis: str = ""
    min_shared_callable_sites: int = MIN_SHARED_CALLABLE_SITES
    linkage: str = "single"
    #: Pairs that could have linked but were not allowed to, and why.
    excluded_uncomputed: List[PairwiseDistance] = field(default_factory=list)
    excluded_thin: List[PairwiseDistance] = field(default_factory=list)
    checks: List[Check] = field(default_factory=list)
    caveats: List[str] = field(default_factory=list)

    def cluster_of(self, sample: str) -> Optional[Cluster]:
        for cluster in self.clusters:
            if sample in cluster.members:
                return cluster
        return None

    def cluster_id_of(self, sample: str) -> str:
        """The cluster id for a sample, or "" — never a fabricated singleton id."""
        found = self.cluster_of(sample)
        return found.cluster_id if found else ""

    @property
    def clustered_samples(self) -> List[str]:
        return sorted((s for c in self.clusters for s in c.members), key=natural_key)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "threshold": self.threshold,
            "threshold_basis": self.threshold_basis,
            "linkage": self.linkage,
            "min_shared_callable_sites": self.min_shared_callable_sites,
            "n_clusters": len(self.clusters),
            "clusters": [c.to_dict() for c in self.clusters],
            "singletons": list(self.singletons),
            "excluded_pairs": {
                "no_denominator": [
                    {"sample_a": p.sample_a, "sample_b": p.sample_b, "note": p.note}
                    for p in self.excluded_uncomputed],
                "insufficient_shared_callable": [
                    {"sample_a": p.sample_a, "sample_b": p.sample_b,
                     "snps": p.snps,
                     "shared_callable_sites": p.shared_callable_sites}
                    for p in self.excluded_thin],
            },
            "checks": [c.to_dict() for c in self.checks],
            "caveats": list(self.caveats),
        }


def _components(samples: Sequence[str],
                edges: Iterable[Tuple[str, str]]) -> List[List[str]]:
    """Single-linkage connected components, by union-find.

    Deterministic: members come back in natural order and the components in the
    order of their first member, so two runs of the same cohort produce the same
    cluster ids and a diff of two reports is about the biology.
    """
    parent: Dict[str, str] = dict((s, s) for s in samples)

    def find(node: str) -> str:
        root = node
        while parent[root] != root:
            root = parent[root]
        while parent[node] != root:  # path compression
            parent[node], node = root, parent[node]
        return root

    for a, b in edges:
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[root_b] = root_a

    groups: Dict[str, List[str]] = {}
    for sample in samples:
        groups.setdefault(find(sample), []).append(sample)
    ordered = [sorted(members, key=natural_key) for members in groups.values()]
    ordered.sort(key=lambda members: (-len(members), natural_key(members[0])))
    return ordered


def cluster_samples(samples: Sequence[str], pairs: Iterable[PairwiseDistance],
                    threshold: int = DEFAULT_CLUSTER_DISTANCE,
                    threshold_basis: str = "",
                    min_shared_callable_sites: int = MIN_SHARED_CALLABLE_SITES,
                    include_singletons: bool = False) -> ClusterAssignment:
    """Group samples whose pairwise distance is at or below *threshold*.

    *pairs* are :class:`~mjolnir.records.PairwiseDistance` records — the whole
    record, not a matrix of integers — because the denominator decides whether a
    pair may form an edge at all, and a function that took bare distances could
    not apply that rule.

    ``include_singletons`` False, the default, returns only groups of two or
    more and lists the rest in :attr:`ClusterAssignment.singletons`. A singleton
    is not a cluster, and giving it a cluster id makes it look like a finding.
    """
    if threshold is None:
        raise MjolnirError(
            "clustering needs an explicit SNP threshold; the TB conventions are {0} "
            "(recent transmission) and {1} (epidemiological link), and the prior "
            "local M. chimaera run used 6".format(CLUSTER_SNP_STRICT, CLUSTER_SNP_RELAXED))
    threshold = int(threshold)
    if threshold < 0:
        raise MjolnirError("clustering threshold must not be negative (got {0})".format(threshold))

    names = list(samples)
    duplicates = sorted(set(n for n in names if names.count(n) > 1))
    if duplicates:
        raise MjolnirError(
            "sample(s) {0} appear more than once in the cohort".format(", ".join(duplicates)))
    known = set(names)

    assignment = ClusterAssignment(
        threshold=threshold,
        threshold_basis=threshold_basis or cluster_threshold_basis(threshold),
        min_shared_callable_sites=int(min_shared_callable_sites))

    edges: List[Tuple[str, str]] = []
    used: Dict[Tuple[str, str], PairwiseDistance] = {}
    considered = 0
    for entry in pairs:
        if entry.sample_a not in known or entry.sample_b not in known:
            # A pair naming a sample outside the cohort is a wiring mistake, not
            # data: it would silently join two cohorts.
            raise MjolnirError(
                "pair {0}/{1} names a sample that is not in this cohort".format(
                    entry.sample_a, entry.sample_b))
        considered += 1
        if entry.snps is None:
            assignment.excluded_uncomputed.append(entry)
            continue
        if (entry.shared_callable_sites or 0) < assignment.min_shared_callable_sites:
            assignment.excluded_thin.append(entry)
            continue
        used[pair_key(entry.sample_a, entry.sample_b)] = entry
        if entry.snps <= threshold:
            edges.append((entry.sample_a, entry.sample_b))

    groups = _components(names, edges)
    counter = 0
    for members in groups:
        if len(members) < 2 and not include_singletons:
            assignment.singletons.extend(members)
            continue
        counter += 1
        cluster = Cluster(
            cluster_id="{0}{1}".format(CLUSTER_ID_PREFIX, counter),
            members=members, threshold=threshold)
        internal = [used[pair_key(a, b)]
                    for i, a in enumerate(members) for b in members[i + 1:]
                    if pair_key(a, b) in used]
        distances = [p.snps for p in internal if p.snps is not None]
        denominators = [p.shared_callable_sites for p in internal
                        if p.shared_callable_sites is not None]
        cluster.max_distance = max(distances) if distances else None
        cluster.min_shared_callable_sites = min(denominators) if denominators else None
        cluster.note = _cluster_note(cluster, threshold, len(internal), len(members))
        assignment.clusters.append(cluster)

    assignment.singletons.sort(key=natural_key)
    assignment.checks.extend(_assignment_checks(assignment, names, considered))
    assignment.caveats.extend(_assignment_caveats(assignment))
    LOG.debug("clustering at %d SNPs: %d clusters, %d singletons, %d pairs excluded",
              threshold, len(assignment.clusters), len(assignment.singletons),
              len(assignment.excluded_uncomputed) + len(assignment.excluded_thin))
    return assignment


def _cluster_note(cluster: Cluster, threshold: int, n_internal: int,
                  n_members: int) -> str:
    """What a reader has to know about this particular group.

    Two things: whether single linkage chained it — an internal maximum above
    the threshold means no single pair of extremes ever met the criterion — and
    whether every internal pair was actually measured, since a group held
    together by a spanning path rather than a complete graph is a weaker claim.
    """
    parts: List[str] = []
    expected = n_members * (n_members - 1) // 2
    if cluster.max_distance is not None and cluster.max_distance > threshold:
        parts.append(
            "single-linkage chaining: the largest distance inside this cluster is {0} "
            "SNPs, above the {1}-SNP threshold. The members are connected through "
            "intermediates, not all directly within the threshold.".format(
                cluster.max_distance, threshold))
    if n_internal < expected:
        parts.append(
            "{0} of {1} internal pairs were comparable; the rest had no usable "
            "shared callable denominator.".format(n_internal, expected))
    if cluster.min_shared_callable_sites is not None:
        parts.append("smallest shared callable denominator inside the cluster: "
                     "{0:,} bp.".format(cluster.min_shared_callable_sites))
    return " ".join(parts)


def _assignment_checks(assignment: ClusterAssignment, samples: Sequence[str],
                       considered: int) -> List[Check]:
    checks: List[Check] = []

    # A threshold with a published basis passes; one an operator invented is a
    # warning, because the number is most of the claim a cluster makes.
    published = assignment.threshold in CLUSTER_THRESHOLD_BASIS
    checks.append(Check(
        name="cluster_threshold", value=assignment.threshold, unit="SNPs",
        status=STATUS_PASS if published else STATUS_WARN, category="cohort",
        source=source_for("default_cluster_distance"),
        reading="clusters were formed by single linkage at {0} SNPs. Basis: {1}".format(
            assignment.threshold, assignment.threshold_basis)))

    if considered == 0:
        checks.append(Check.not_measured(
            "clustering",
            "no pairwise distances were available, so no sample was placed in or "
            "excluded from a cluster. {0} sample(s) are unclustered because nothing "
            "was compared, not because they are unrelated.".format(len(samples)),
            category="cohort"))
        return checks

    excluded = len(assignment.excluded_uncomputed) + len(assignment.excluded_thin)
    checks.append(Check(
        name="pairs_used_for_clustering", value=considered - excluded,
        threshold=considered, comparison="==", unit="pairs", category="cohort",
        status=STATUS_PASS if not excluded else STATUS_WARN,
        source=source_for("min_shared_callable_sites"),
        reading="{0} of {1} pairs formed the clustering graph. {2} had no shared "
                "callable denominator and {3} were compared over less than {4:,} bp "
                "of shared callable sequence; neither kind can support a link at a "
                "published SNP threshold, so neither was allowed to form one.".format(
                    considered - excluded, considered,
                    len(assignment.excluded_uncomputed), len(assignment.excluded_thin),
                    assignment.min_shared_callable_sites)))

    chained = [c for c in assignment.clusters
               if c.max_distance is not None and c.max_distance > assignment.threshold]
    checks.append(Check(
        name="single_linkage_chaining", value=len(chained), threshold=0,
        comparison="<=", unit="clusters", category="cohort",
        status=STATUS_PASS if not chained else STATUS_WARN,
        source=source_for("cluster_snp_relaxed"),
        reading=("no cluster is held together by chaining: every cluster's internal "
                 "maximum is within the threshold."
                 if not chained else
                 "{0} cluster(s) are chained — {1} — where the largest internal "
                 "distance exceeds the {2}-SNP threshold and the members are linked "
                 "through intermediates.".format(
                     len(chained), ", ".join(c.cluster_id for c in chained),
                     assignment.threshold))))

    return checks


def _assignment_caveats(assignment: ClusterAssignment) -> List[str]:
    caveats: List[str] = []
    if assignment.excluded_uncomputed:
        caveats.append(
            "{0} pair(s) had no shared callable denominator and could form no link; "
            "an absent comparison is not evidence that two samples are unrelated"
            .format(len(assignment.excluded_uncomputed)))
    if assignment.excluded_thin:
        caveats.append(
            "{0} pair(s) were compared over less than {1:,} bp of shared callable "
            "sequence and were excluded from clustering; their SNP counts are not "
            "comparable to the published thresholds".format(
                len(assignment.excluded_thin), assignment.min_shared_callable_sites))
    if any(c.max_distance is not None and c.max_distance > assignment.threshold
           for c in assignment.clusters):
        caveats.append(
            "single linkage chains: at least one cluster contains a pair further "
            "apart than the {0}-SNP threshold, connected through intermediates"
            .format(assignment.threshold))
    return caveats


def clusters_at(samples: Sequence[str], pairs: Iterable[PairwiseDistance],
                thresholds: Sequence[int] = (CLUSTER_SNP_STRICT, CLUSTER_SNP_RELAXED),
                min_shared_callable_sites: int = MIN_SHARED_CALLABLE_SITES,
                ) -> Dict[int, ClusterAssignment]:
    """The same cohort clustered at several thresholds.

    Both TB conventions at once is the honest presentation when the epidemiology
    is not settled: 5 SNPs and 12 SNPs answer different questions, and showing
    which links survive the stricter one is more informative than defending a
    single cut. The pairs are computed once and reused, since the threshold only
    decides which edges are drawn.
    """
    materialised = list(pairs)
    return dict(
        (int(value), cluster_samples(samples, materialised, threshold=int(value),
                                     min_shared_callable_sites=min_shared_callable_sites))
        for value in thresholds)
