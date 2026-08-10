"""Render a Mjolnir cohort report from a synthetic outbreak, with no external tools.

The shape is taken from a real multi-site *M. chimaera* investigation: clinical
isolates from several hospitals, plus water and heater-cooler-unit isolates from
the same sites. The question such a cohort is run to answer is whether the
environmental isolates and the clinical ones are the same organism — so the
example deliberately includes one isolate that is *not* part of the cluster, and
one pair that could not be compared at all.

Two things this example exists to demonstrate:

Every distance carries its denominator. ``PairwiseDistance`` cannot be built
without ``shared_callable_sites``, because 6 differences over 4.9 Mb of shared
callable sequence and 6 over 900 kb are not the same statement, and a report
that prints only the left-hand number invites the reader to treat them alike.

A pair that could not be compared is absent from the matrix rather than scored
as zero — which would place two unrelatable isolates at the centre of the
outbreak.

    PYTHONPATH=src python examples/demo_cohort.py OUTDIR
"""

import os
import sys

from mjolnir.records import (Check, Cluster, CohortResult, DatabaseVersion,
                             PairwiseDistance)
from mjolnir.report import write_cohort_pdf, write_cohort_html

REFERENCE = "NZ_CP015272.1 (M. chimaera ZUERICH-1)"
THRESHOLD = 6

#: Site of origin, and whether the isolate came from a patient or from plumbing.
SAMPLES = [
    ("OSR-30-20", "clinical"),
    ("OSR-31-20", "clinical"),
    ("PAVIA-590-19", "clinical"),
    ("ACQUA-12-20", "water"),
    ("SCAMB-04-20", "heater-cooler"),
    ("VENETO-771-19", "clinical"),
]

#: The outbreak clade: five isolates within a few SNPs of each other, spanning
#: three hospitals and both environmental sources. VENETO-771-19 sits far
#: outside it, which is what a genuine unrelated isolate looks like.
_DISTANCES = {
    ("OSR-30-20", "OSR-31-20"): (1, 4_912_004),
    ("OSR-30-20", "PAVIA-590-19"): (4, 4_887_551),
    ("OSR-30-20", "ACQUA-12-20"): (3, 4_901_233),
    ("OSR-30-20", "SCAMB-04-20"): (2, 4_895_120),
    ("OSR-31-20", "PAVIA-590-19"): (5, 4_880_902),
    ("OSR-31-20", "ACQUA-12-20"): (4, 4_893_774),
    ("OSR-31-20", "SCAMB-04-20"): (3, 4_888_610),
    ("PAVIA-590-19", "ACQUA-12-20"): (6, 4_871_449),
    ("PAVIA-590-19", "SCAMB-04-20"): (5, 4_866_003),
    # Compared over a much smaller shared span. The distance is inside the
    # threshold, but the denominator is what tells the reader why to hesitate.
    ("ACQUA-12-20", "SCAMB-04-20"): (2, 902_774),
    ("OSR-30-20", "VENETO-771-19"): (1_482, 4_842_119),
    ("OSR-31-20", "VENETO-771-19"): (1_479, 4_838_667),
    ("PAVIA-590-19", "VENETO-771-19"): (1_501, 4_830_115),
    ("ACQUA-12-20", "VENETO-771-19"): (1_488, 4_845_002),
    # Not comparable: this pair shared too little callable sequence to score.
    ("SCAMB-04-20", "VENETO-771-19"): (None, 0),
}

MASKED_SITES = 291_044


def pairs():
    out = []
    for (a, b), (snps, shared) in _DISTANCES.items():
        note = ""
        if snps is None:
            note = ("shared callable sequence below the reporting floor; "
                    "absent from the matrix rather than scored as zero")
        elif shared < 3_000_000:
            note = ("compared over {0:,} shared callable sites, far below the "
                    "cohort median: read this distance with that in mind".format(shared))
        out.append(PairwiseDistance(sample_a=a, sample_b=b, snps=snps,
                                    shared_callable_sites=shared,
                                    masked_sites=MASKED_SITES, note=note))
    return out


CLUSTERS = [
    Cluster(cluster_id="C1",
            members=["OSR-30-20", "OSR-31-20", "PAVIA-590-19",
                     "ACQUA-12-20", "SCAMB-04-20"],
            threshold=THRESHOLD, max_distance=6,
            min_shared_callable_sites=902_774,
            note=("spans three hospitals and includes both a water isolate and a "
                  "heater-cooler isolate. The smallest shared span behind any "
                  "member pair is 902,774 sites, which is the weakest link in "
                  "this cluster and the one to check before drawing conclusions")),
]

CHECKS = [
    Check(name="cohort size", value=len(SAMPLES), threshold=2, source="design §9",
          status="pass", comparison=">=", measured=True,
          reading="Large enough for a joint variant table."),
    Check(name="masked fraction of the reference", value=0.0594, threshold=None,
          source="tbdb mask.bed; the mask is a named, versioned input, not a constant",
          status="pass", measured=True,
          reading="5.94% of the reference was masked before any distance was computed."),
    Check(name="pairs compared", value=14, threshold=15, source="design §9",
          status="warn", comparison=">=", measured=True,
          reading="One of 15 pairs shared too little callable sequence to compare."),
]


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(out, exist_ok=True)
    cohort = CohortResult(
        samples=[name for name, _source in SAMPLES],
        pairs=pairs(),
        clusters=CLUSTERS,
        threshold=THRESHOLD,
        threshold_basis=("6 SNPs, the threshold used in the original investigation "
                         "of this collection. The 5- and 12-SNP conventions come "
                         "from M. tuberculosis transmission studies and are not "
                         "established for M. chimaera"),
        mask_name="tbdb mask.bed (merged Marin + Modlin, 2025-08)",
        masked_sites=MASKED_SITES,
        masked_fraction=0.0594,
        joint_sites=4_912_004,
        reference=REFERENCE,
        checks=CHECKS,
        caveats=["one pair could not be compared and is absent from the matrix",
                 "a SNP distance is not a transmission event"],
        tool_versions={"bwa-mem2": "2.2.1", "samtools": "1.19", "bcftools": "1.19"},
        database_versions=[
            DatabaseVersion(name="M. chimaera reference", version="NZ_CP015272.1",
                            checksum="sha256:9ac41b…", licence="NCBI-public",
                            citation="GenBank NZ_CP015272.1"),
        ],
        mjolnir_version="0.1.0")

    print("PDF :", write_cohort_pdf(os.path.join(out, "DEMO-COHORT.pdf"), cohort))
    print("HTML:", write_cohort_html(os.path.join(out, "DEMO-COHORT.html"), cohort))


if __name__ == "__main__":
    main()
