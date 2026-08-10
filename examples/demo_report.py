"""Render a Mjolnir PDF from a synthetic MDR-TB sample, to prove the deliverable works."""
import sys
# Run from a source checkout without installing:
#   PYTHONPATH=src python examples/demo_report.py OUTDIR

from mjolnir.records import (SampleResult, DrugCall, SpeciesCall, LineageCall,
                             QCMetrics, Check, ContaminationResult, Variant,
                             CatalogueCall, DatabaseVersion)
from mjolnir.config import call_for_grade
from mjolnir.resistance.catalogues import (CATALOGUE_MTBSEQ, CATALOGUE_TBDB,
                                            CATALOGUE_WHO)
from mjolnir.report import write_pdf, write_html

def V(gene, hgvs, pos, ref, alt, af, depth, calls):
    return Variant(chrom="NC_000962.3", pos=pos, ref=ref, alt=alt, gene=gene,
                   hgvs=hgvs, depth=depth, allele_fraction=af, is_major=af >= 0.9,
                   source_caller="bcftools", alt_reads=int(depth * af),
                   ref_reads=int(depth * (1 - af)), qual=222.0,
                   variant_type="snp", effect="missense_variant",
                   catalogue_calls=calls)

def CC(cat, drug, grade, comment="", call=None):
    """One catalogue's own answer for one drug.

    ``call`` is derived from the WHO grade when there is one, because a grade of
    "1) Assoc w R" beside a call of "no-call" is the kind of internal
    inconsistency this report is built to make impossible. MTBseq is flat and
    ungraded, so its rows assert resistance directly.
    """
    if call is None:
        call = call_for_grade(grade) if grade.startswith(("1)", "2)", "3)", "4)", "5)")) \
            else ("R" if cat in (CATALOGUE_MTBSEQ, CATALOGUE_TBDB) else "no-call")
    return CatalogueCall(catalogue=cat, drug=drug, grade=grade, comment=comment,
                         source=cat, call=call)

variants = [
    V("rpoB", "p.Ser450Leu", 761155, "C", "T", 0.99, 88, [
        CC(CATALOGUE_WHO, "Rifampicin", "1) Assoc w R", "High-level resistance"),
        CC(CATALOGUE_MTBSEQ, "Rifampicin", ""), CC(CATALOGUE_TBDB, "Rifampicin", "high")]),
    V("katG", "p.Ser315Thr", 2155168, "G", "C", 0.98, 76, [
        CC(CATALOGUE_WHO, "Isoniazid", "1) Assoc w R", "High-level resistance"),
        CC(CATALOGUE_MTBSEQ, "Isoniazid", ""), CC(CATALOGUE_TBDB, "Isoniazid", "high")]),
    V("gyrA", "p.Ala90Val", 7570, "C", "T", 0.41, 61, [
        CC(CATALOGUE_WHO, "Levofloxacin", "1) Assoc w R"),
        CC(CATALOGUE_WHO, "Moxifloxacin", "2) Assoc w R - Interim")]),
    V("pncA", "c.-11A>C", 2289251, "T", "G", 0.97, 54, [
        CC(CATALOGUE_WHO, "Pyrazinamide", "2) Assoc w R - Interim")]),
    V("Rv0678", "p.Asp47fs", 779382, "GA", "G", 0.95, 47, [
        CC(CATALOGUE_MTBSEQ, "Bedaquiline", ""), CC(CATALOGUE_TBDB, "Bedaquiline", "moderate")]),
    V("rrs", "n.1401A>G", 1473246, "A", "G", 0.12, 66, [
        CC(CATALOGUE_WHO, "Amikacin", "1) Assoc w R")]),
]

def D(drug, call, conf, who_grade="", who=True, dis=False, caveats=(), lvl="",
      sup="", note="", determinants=(), cats=()):
    """One drug's finished call.

    ``determinants`` and ``cats`` are populated together on purpose: a report
    showing "R" beside "none matched" would be internally inconsistent, and the
    point of the example is to exercise the layout a real run produces.
    """
    return DrugCall(drug=drug, call=call, confidence=conf, who_graded=who,
                    who_grade=who_grade, disagreement=dis, level=lvl,
                    caveats=list(caveats), suppressed_by=sup, note=note,
                    target_covered=True, supporting_variants=list(determinants),
                    catalogue_calls=list(cats))

drugs = [
    D("Rifampicin", "R", "high", "1) Assoc w R", lvl="High-level resistance",
      determinants=["rpoB_p.Ser450Leu"],
      cats=[CC(CATALOGUE_WHO, "Rifampicin", "1) Assoc w R", "High-level resistance"),
            CC(CATALOGUE_MTBSEQ, "Rifampicin", ""), CC(CATALOGUE_TBDB, "Rifampicin", "high")]),
    D("Isoniazid", "R", "high", "1) Assoc w R", lvl="High-level resistance",
      determinants=["katG_p.Ser315Thr"],
      cats=[CC(CATALOGUE_WHO, "Isoniazid", "1) Assoc w R", "High-level resistance"),
            CC(CATALOGUE_MTBSEQ, "Isoniazid", ""), CC(CATALOGUE_TBDB, "Isoniazid", "high")]),
    D("Levofloxacin", "R", "moderate", "1) Assoc w R",
      determinants=["gyrA_p.Ala90Val"],
      cats=[CC(CATALOGUE_WHO, "Levofloxacin", "1) Assoc w R")],
      caveats=["allele fraction 0.41 - minority subpopulation, not a fixed variant"]),
    D("Moxifloxacin", "R-interim", "moderate", "2) Assoc w R - Interim",
      determinants=["gyrA_p.Ala90Val"],
      cats=[CC(CATALOGUE_WHO, "Moxifloxacin", "2) Assoc w R - Interim")],
      caveats=["allele fraction 0.41 - minority subpopulation"]),
    D("Pyrazinamide", "R-interim", "moderate", "2) Assoc w R - Interim",
      determinants=["pncA_c.-11A>C"],
      cats=[CC(CATALOGUE_WHO, "Pyrazinamide", "2) Assoc w R - Interim")]),
    D("Bedaquiline", "R-outside-WHO", "low", "", who=False, dis=True,
      determinants=["Rv0678_p.Asp47fs"],
      cats=[CC(CATALOGUE_MTBSEQ, "Bedaquiline", ""), CC(CATALOGUE_TBDB, "Bedaquiline", "moderate")],
      caveats=["called by MTBseq and tbdb; not graded by the WHO catalogue v2",
               "loss of function in Rv0678"]),
    D("Amikacin", "no-call", "low", "1) Assoc w R",
      determinants=["rrs_n.1401A>G"],
      cats=[CC(CATALOGUE_WHO, "Amikacin", "1) Assoc w R")],
      caveats=["rrs n.1401A>G present at allele fraction 0.12 - below the 0.90 "
               "major-variant threshold; reported as a minority variant, not a call"]),
    D("Ethambutol", "no-call", "high"),
    D("Linezolid", "no-call", "high"),
    D("Clofazimine", "no-call", "low",
      caveats=["Rv0678 loss of function present; WHO grades this for clofazimine "
               "only in combination with an intact mmpL5"]),
    D("Delamanid", "no-call", "high"),
    D("Streptomycin", "no-call", "high"),
    D("Ethionamide", "no-call", "high"),
    D("Capreomycin", "no-call", "high"),
    D("Kanamycin", "no-call", "high"),
]

qc = QCMetrics(mean_depth=68.4, median_depth=71.0, breadth_1x=0.988,
               breadth_10x=0.974, breadth_min_depth=0.961, coverage_evenness=0.91,
               evenness_definition="fraction of sites within 0.5-2x of median depth",
               mapped_fraction=0.947, gc_content=0.654, unambiguous_fraction=0.968,
               total_reads=2_845_112, mapped_reads=2_694_321,
               duplicate_fraction=0.038, mean_read_length=148.2,
               mean_base_quality=35.1, reference="NC_000962.3",
               reference_length=4_411_532)

cont = ContaminationResult(
    f2=0.021, f47=0.014, lineage_het_sites=3, lineage_sites_examined=1111,
    het_snp_fraction=0.019, het_snp_count=57, snp_sites_examined=3012,
    mixture_class="single-strain",
    unambiguous_fraction=0.968, non_target_fraction=None,
    non_target_resolution="not measured - no mycobacterial pangenome index present",
    mapped_fraction=0.947, coverage_breadth=0.974, coverage_evenness=0.91,
    gc_content=0.654, verdict="valid",
    verdict_reason=("F2 0.021 and a genome-wide heterozygous-SNP fraction of 0.019 are "
                    "below the mixed-infection tier but above zero. Resistance calling "
                    "tolerates this; a 5-SNP transmission threshold does not."),
    screen_informative=False, screen_method="kraken2",
    screen_note=("no mycobacterial pangenome index was supplied. Measured Kraken2 "
                 "sensitivity for M. tuberculosis reads against a standard index is "
                 "0.0731 on real Illumina data, so no composition screen was run and "
                 "none is reported."),
    caveats=["read-composition screen not performed - see screen note"])

checks = [
    Check(name="mean depth", value=68.4, threshold=25.0, source="Hall et al. eLife 2024",
          status="pass", comparison=">=", unit="x", measured=True,
          reading="Depth is well above the 25x floor below which precision and recall degrade."),
    Check(name="coverage breadth at 10x", value=0.974, threshold=0.95,
          source="design section 8", status="pass", comparison=">=", measured=True),
    Check(name="mapped fraction", value=0.947, threshold=0.90, source="design section 8",
          status="pass", comparison=">=", measured=True),
    Check(name="heterozygous SNP fraction", value=0.019, threshold=0.05,
          source="Sobkowiak et al. BMC Genomics 2018 (MixInfect)", status="warn",
          comparison="<=", measured=True,
          reading="Below the mixed-infection tier, but not zero."),
    Check(name="read composition screen", value=None, threshold=None,
          source="design section 8", status="warn", measured=False,
          reading="Not measured: no mycobacterial pangenome index present."),
]

result = SampleResult(
    sample_id="DEMO-MDR-001",
    platform="illumina",
    inputs=["DEMO-MDR-001_S1_R1_001.fastq.gz", "DEMO-MDR-001_S1_R2_001.fastq.gz"],
    reference="NC_000962.3 (M. tuberculosis H37Rv)",
    species=SpeciesCall(name="Mycobacterium tuberculosis", complex="MTBC",
                        method="skani ANI + lineage barcode", ani=99.91,
                        confidence="high", resolved_to_species=True,
                        reference="H37Rv NC_000962.3", aligned_fraction=0.981),
    lineage=LineageCall(lineage="lineage2", sublineage="lineage2.2.1",
                        barcode_sites_supporting=41, barcode_sites_total=1111,
                        barcode_sites_callable=1094, is_bcg=False, is_animal=False,
                        scheme="tbdb barcode.bed", method="direct pileup",
                        confidence="high", support=[
                            {"position": 1834177, "lineage": "lineage2", "expected": "G",
                             "observed": "G", "depth": 71, "allele_fraction": 1.0,
                             "supports": True, "note": "East-Asian (Beijing), RD105"},
                            {"position": 2505085, "lineage": "lineage2.2.1", "expected": "A",
                             "observed": "A", "depth": 66, "allele_fraction": 1.0,
                             "supports": True, "note": "Beijing RD207"},
                        ]),
    variants=variants, drugs=drugs, qc=qc, contamination=cont, checks=checks,
    caveats=["no read-composition screen was run: no mycobacterial pangenome index present"],
    tool_versions={"bwa-mem2": "2.2.1", "samtools": "1.19", "bcftools": "1.19"},
    database_versions=[
        DatabaseVersion(name="WHO catalogue", version="v2 (2023) WHO-UCN-TB-2023.7",
                        checksum="sha256:3f9c1a…", licence="ODC-By v1.0",
                        citation="WHO. Catalogue of mutations in M. tuberculosis complex, 2nd ed. 2023."),
        DatabaseVersion(name="tbdb", version="7fe4364e", checksum="sha256:b21e77…",
                        licence="see repo", citation="Phelan JE et al. Genome Med 2019;11:41."),
        DatabaseVersion(name="MTBseq ResSeq", version="v1.1.0", checksum="sha256:dcc3fa…",
                        licence="GPL-3.0-or-later", citation="Kohl TA et al. PeerJ 2018;6:e5895."),
    ],
    interpretation=None, mjolnir_version="0.1.0", profile="clinical",
    status="complete")

import os
out = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(out, exist_ok=True)
p = write_pdf(out + "/DEMO-MDR-001.pdf", result)
h = write_html(out + "/DEMO-MDR-001.html", result)
print("PDF :", p)
print("HTML:", h)
