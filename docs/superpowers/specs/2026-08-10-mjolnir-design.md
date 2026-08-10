# Mjolnir — design

**Mycobacterial Junction and Omics Locus Nucleotide Identification for Resistance**

Resistance, lineage, species and contamination calling for the *M. tuberculosis*
complex and for non-tuberculous mycobacteria, from Illumina reads, ONT reads or
assemblies — one command, one report.

Status: design. Date: 2026-08-10. Licence: MIT. Repo: `github.com/iowa69/mjolnir`.

---

## 1. Why this exists

The two tools being consolidated are [MTBseq](https://github.com/ngs-fzb/MTBseq_source)
and [NTMseq](https://github.com/ngs-fzb/NTMseq), both from FZ Borstel.

MTBseq is a real variant-calling pipeline, built around H37Rv. It ships four
references — `M._tuberculosis_H37Rv_2015-11-13`, `M._abscessus_CIP-104536T`,
`M._chimaera_DSM44623`, `M._fortuitum_CT6` — so NTM is nominally supported by
swapping `--ref`. But every interpretation layer is MTB-only:
`var/res/MTB_Resistance_Mediating.txt`, `var/res/MTB_Extended_Resistance_Mediating.txt`,
`var/cat/MTB_Gene_Categories.txt`, `var/res/MTB_Base_Calibration_List.vcf`.

The consequence is visible in a real run on this machine
(`M_Chimaera_TN/samp/MTBseq_2022-11-11_iowa.log`, MTBseq 1.0.3, 2022-11-11):

```
--ref       M._chimaera_DSM44623_2016-01-28.fasta
--resilist  NONE
--categories NONE
--basecalib NONE
...
<ERROR> No joint variant file try_joint_cf5_cr5_fr75_ph4_samples2.tab to amend!
```

A *M. chimaera* run gets no resistance, no gene categories, no base
recalibration, no lineage — and then fails at `TBamend`. NTMseq does not fill
that gap: it is a collection of bash starter scripts driving other people's
tools (shovill, SRST2, mashtree, NTM-Profiler, Platon, abricate) plus
`kraken_parse_results.v2.0.0.pl` for contamination. Its only shipped database is
208 mycobacterial plasmids.

So today, answering "what is this isolate, is it resistant, is it contaminated,
is it part of an outbreak" for a mycobacterium means running two dissimilar
tools, neither of which answers all of it, and hand-assembling the result.

**Mjolnir is one tool where NTM is first-class rather than a reference swap**,
where resistance is called from graded catalogues rather than a flat list, where
contamination is measured by methods that actually work on mycobacteria, and
where the output is a report a clinician can read.

## 2. Scope

In scope for v1:

- Inputs: Illumina paired-end FASTQ, ONT FASTQ, assembled FASTA.
- Organisms: MTBC (including animal lineages and BCG) and NTM.
- Outputs: species, lineage/sublineage, per-drug resistance with confidence,
  contamination and sample-validity verdict, cohort SNP distances and clusters,
  a clinician-first PDF, and machine-readable TSV/JSON.
- A local-LLM interpretation agent that writes prose over rule-derived verdicts.

Out of scope for v1: de novo assembly (accept assemblies, do not produce them),
metagenomic/direct-from-sputum input, spoligotype prediction, phenotypic DST
integration, any cloud service.

## 3. Non-negotiable house rules

Inherited from `hydra` and `tesseract-ai`, because they are what make the output
trustworthy:

1. **Every number carries its provenance.** A threshold in the report names its
   source. Numbers in the README are generated from a measured file, never typed.
2. **The verdict is a rule; the prose is the model.** Pass/warn/fail is computed
   in Python from a stated threshold before the LLM is called. The model receives
   finished checks and writes the reading. A model answer that violates the
   discipline rules is discarded and replaced by the rule-derived summary, and
   the report says so.
3. **The model never sees raw sequence.** Enforced in code: an observation
   containing a long nucleotide run raises before it can reach the model.
4. **It runs without the model.** If no LLM host is reachable, every gate takes
   its declared default and the report states that the interpretation is
   rule-only.
5. **Nothing unmeasured is described as fine.** Absence of a call is reported as
   absence, never as susceptibility.

## 4. Architecture

```
mjolnir/
  cli.py            command surface
  config.py         thresholds, each with a documented source
  records.py        the data model every stage reads and writes
  seqio.py          input detection: paired FASTQ, ONT FASTQ, FASTA
  shell.py          external tool runner, version capture
  doctor.py         reports the whole environment before running anything
  pipeline.py       single-sample orchestration
  engines/
    map.py          bwa-mem2 / minimap2 by platform
    call.py         variant calling by platform (see §7)
    pileup.py       direct pileup for barcode sites and AF at catalogue positions
    depth.py        coverage, breadth, evenness
  typing/
    species.py      mycobacterial species ID (§6)
    lineage.py      MTBC lineage/sublineage barcode; NTM subspecies
  resistance/
    catalogues.py   WHO v2, MTBseq ResSeq, tbdb loaders -> one normal form
    normalise.py    variant normalisation to (CHROM,POS,REF,ALT) + HGVS
    rules.py        WHO additional grading rules, epistasis, LoF
    consensus.py    the three-catalogue consensus engine (§5)
    ntm.py          NTM resistance: erm(41), rrl, rrs
  contamination/
    purity.py       read-level composition, honestly bounded (§8)
    heterozygosity.py  F2/F47 and genome-wide hSNP fraction
  cohort/
    joint.py        joint variant table across samples
    distance.py     masked pairwise SNP distance
    cluster.py      threshold clustering
  agent/
    client.py       ollama native + OpenAI /v1, auto-detected
    observation.py  what the model is allowed to see
    discipline.py   rejects violating answers
    playbooks/      mtbc.yaml, ntm.yaml
  report/
    pdf.py          the deliverable
    html.py         same content, self-contained HTML
    tables.py       TSV/JSON artefacts
  db/
    fetch.py        download + verify catalogues and references
    registry.py     what each database is, its version, its licence
```

Data flow, single sample:

```
input -> detect platform -> QC/trim -> map -> call variants -> pileup at
catalogue + barcode positions -> species -> lineage -> resistance consensus ->
contamination -> rule-derived verdicts -> agent prose -> PDF
```

Cohort adds: joint variant table -> masked distance matrix -> clusters -> cohort PDF.

## 5. Resistance: the consensus engine

Three catalogues, called independently, then reconciled.

### 5.1 The sources

| Source | What it is | Licence | Notes |
|---|---|---|---|
| **WHO catalogue v2** | `WHO-UCN-TB-2023.7-eng.xlsx`, GTB-tbsequencing/mutation-catalogue-2023, `Final Result Files/` | **ODC-By v1.0** — redistributable with attribution | 48,152 variant–drug rows, 30,699 unique variants, 65 genes, 15 drugs |
| **MTBseq ResSeq** | `var/res/MTB_Resistance_Mediating.txt` (275 KB) + `MTB_Extended_Resistance_Mediating.txt` | GPL-3.0 (MTBseq) | flat list, no confidence grading |
| **tbdb** | `mutations.csv` (4.7 MB) from jodyphelan/tbdb | check repo LICENCE at fetch time | TB-Profiler's library; also supplies `barcode.bed` and `mask.bed` |

WHO v2 is the **latest edition as of 2026-08-10**. A 3rd edition was called for
on 2024-08-26 (deadline 2024-10-15) but is unreleased; `db fetch` must not
hard-code the assumption that v2 is final.

### 5.2 Traps that must be handled in code

Each of these is a verified failure mode, not a hypothetical:

- The `.xlsx` header is on **row 3** (rows 1–2 are a merged banner). Naive
  `read_excel` yields garbage columns.
- The repo's `.txt` master file is **not** equivalent to the `.xlsx`: 40,178 rows
  vs 48,152, 14 drugs vs 15 (Streptomycin absent entirely), rpoB/rpoC/rpsL/gid
  affected. **Mjolnir reads the xlsx.**
- The `genomic position` column is a decoy — 38,884 of 48,152 rows contain the
  literal string `(see "Genomic_coordinates" sheet)`.
- Grade strings are numeric-prefixed with a spaced ASCII hyphen:
  `1) Assoc w R`, `2) Assoc w R - Interim`, `3) Uncertain significance`,
  `4) Not assoc w R - Interim`, `5) Not assoc w R`. The PDF prints an en-dash
  form; hard-coding the PDF form fails to match.
- Grading is per **(drug, variant)**. `inhA_c.-154G>A` is Group 1 for isoniazid
  and Group 2 for ethionamide. Deduplicating by variant name is wrong.
- MNVs are decomposed; one genomic variant maps to several graded variants,
  joined by `&` in the VCF INFO. Must split.
- No coordinates exist for `<gene>_deletion` and pooled LoF graded-variants.
  These are matched by rule, not by coordinate.

### 5.3 Matching

WHO's own documented protocol is **coordinate-based**, not name-based: exact
match on `(chromosome, position, reference_nucleotide, alternative_nucleotide)`
against NC_000962.3 via the `Genomic_coordinates` sheet, then look up the graded
variant in `Catalogue_master_file`. Mjolnir follows that protocol as the primary
path, and uses normalised HGVS only as the cross-catalogue join key, since
MTBseq and tbdb do not share WHO's coordinate table.

`normalise.py` therefore produces, for every observed variant, both:
- a coordinate key `(NC_000962.3, POS, REF, ALT)` — left-aligned, parsimonious;
- an HGVS key `<gene>_<hgvs>` with 3-letter amino acids.

Cross-catalogue comparison happens on the HGVS key with a documented alias table
for legacy numbering; disagreements caused by numbering alone are reported as
**nomenclature**, not as biological disagreement.

### 5.4 Rules beyond the table

The catalogue is three components, not one table. Table lookup alone is
incorrect. Mjolnir implements:

- Any novel **silent** variant → Group 4.
- RIF: any non-synonymous mutation or indel in the rpoB RRDR (codons 426–452) → Group 2.
- INH: LoF in *katG* → Group 2. PZA: LoF in *pncA* → Group 2.
- BDQ/CFZ: LoF in *Rv0678*, *pepQ* → Group 2. DLM/PMD: LoF in the nitroreductase set → Group 2.
- **Epistasis**: *mmpL5* LoF abrogates *Rv0678* effects (BDQ/CFZ); *eis* coding
  LoF abrogates *eis* promoter mutations (AMK/KAN). Encoded in the `Comment`
  column and applied as a suppression step, with the suppression stated in the report.
- The four borderline *rpoB* mutations (Leu430Pro, His445Asn, His445Ser,
  Ile491Phe) are Group 1 by rule.
- Level of resistance and cross-resistance come from the `Comment` column, not
  the grade — `High-level resistance`, additive low-level mutations, DLM–PMD
  cross-resistance.

### 5.5 The consensus rule

Per drug, per sample:

1. Each catalogue independently yields a call in {R, R-interim, Uncertain, S-interim, S, no-call}.
   WHO grades map directly. MTBseq's flat list maps to R or no-call — it has no
   grading, and this asymmetry is stated in the report rather than hidden.
   tbdb maps via its own confidence field.
2. **WHO is the anchor.** Where WHO grades the variant, the WHO grade is the
   Mjolnir call. This is defensible: it is the only source with a published,
   systematically-derived grading.
3. Where WHO does **not** grade it and another catalogue calls R, the drug is
   reported as **`R (outside WHO catalogue)`** — surfaced, never silently dropped,
   and never presented as equivalent to a WHO Group 1 call.
4. Where catalogues conflict, the report shows all three side by side in the
   annex, and the front page carries a disagreement flag for that drug.
5. Absence of any catalogued mutation is reported as **"no resistance determinant
   detected"**, never as "susceptible".

Known failure modes of this rule, stated in the docs: it inherits WHO's blind
spots by construction; it cannot adjudicate a true novel mechanism; and a
catalogue-version mismatch between installations changes calls, so the report
prints the version and checksum of every catalogue used.

### 5.6 NTM resistance

Not covered by any of the three catalogues. Implemented directly from the
literature, per organism, with explicit gene targets:

- **Macrolides (*M. abscessus*)**: `erm(41)` sequevar typing — the T28C
  polymorphism separates inducible-resistant (T28) from susceptible (C28);
  truncated `erm(41)` in *M. abscessus* subsp. *massiliense*. Plus acquired
  `rrl` 2058/2059 mutations for constitutive resistance.
- **Amikacin**: `rrs` 1408 and neighbouring positions.
- **MAC / *M. chimaera***: `rrl` 2058/2059 (macrolide), `rrs` 1408 (amikacin).

Every NTM call names its supporting reference in the report. Where no evidence
base exists for a species–drug pair, the tool says so instead of guessing.

## 6. Species and lineage

**Species ID must not come from a taxonomic read classifier.** In current NCBI
taxonomy the MTBC members are not at species rank at all — `Mycobacterium
tuberculosis variant bovis` (taxid 1765) has rank `no rank` under species
*M. tuberculosis*, because MTBC members are later heterotypic synonyms of
*M. tuberculosis* (ANI 99.21–99.92%). A Kraken2 row saying "M. bovis 3.2%" is
not a species identification and Mjolnir must never print one.

So:

- **Genus/species level**: ANI-based (skani/mash) against a curated mycobacterial
  reference set, with an explicit ANI floor for a species claim and a documented
  "cannot resolve below complex" outcome for MAC and MTBC.
- **Within MTBC**: lineage-defining SNP barcodes only. `barcode.bed` from tbdb
  supplies the scheme. Animal lineages and BCG are called from their defining
  SNPs, with the caveat that *M. bovis* is defined by very few phylogenetic SNPs
  (23 in SNP-IT) and is therefore highly sensitive to coverage gaps and
  contamination — reported as a confidence caveat on the call.
- **BCG matters clinically** (intrinsic pyrazinamide resistance) and is called
  and flagged explicitly.
- **Within MAC**: *M. chimaera* vs *M. intracellulare* vs *M. avium* resolved by
  ANI plus marker SNPs, since this is exactly the distinction the outbreak data
  on hand requires.

Barcode genotyping is done from a **direct pileup**, not from the variant caller,
following pathogen-profiler's approach — and on ONT the highest-depth allele is
taken at each barcode site.

## 7. Platform handling

| | Illumina | ONT | FASTA |
|---|---|---|---|
| Map | bwa-mem2 (bwa fallback) | minimap2 `-x map-ont` | — |
| Call | bcftools/freebayes | **Clair3** preferred | direct comparison |
| Min reads for a variant | **≥3** | **≥5** | n/a |
| Major variant | ≥90% read support | ≥90% read support | n/a |
| Min depth | 25× target | **25× floor**, 10× degraded | n/a |
| Minor/heteroresistance | reported | **reported as under-detected** | **not available** |

Grounding: R10.4.1 + Dorado `sup` is the minimum credible ONT configuration
(`fast` is not acceptable); Clair3/DeepVariant lead on ONT bacterial data while
BCFtools is specifically weak on indels; precision and recall degrade notably
below 25×; the ≥5 reads (ONT) / ≥3 reads (Illumina) and ≥90% major-variant
thresholds are the published clinical-DST thresholds from a 508-isolate
ONT-vs-Illumina MTBC study.

Three consequences the report must state, per platform:

1. **ONT under-detects minor resistance variants.** In the 508-isolate study, 26
   of 27 Illumina-only minor SNPs were visible in the ONT pileup but not called.
   An ONT report must not present absence of a minor variant as absence of a
   subpopulation.
2. **`fbiC` tandem-repeat deletions cause spurious delamanid resistance on ONT** —
   47.2% of all discordant drug classifications in that study. Mjolnir suppresses
   this specific call on ONT and says why.
3. **ONT indel calls are ~16.6% uncorroborated by Illumina.** Indel-driven
   resistance calls on ONT (notably LoF rules) carry a platform caveat.

FASTA input loses allele fractions entirely: no heteroresistance, no mixed
infection detection, no contamination heterozygosity metric. The report states
this as a capability loss, not as a clean result.

## 8. Contamination — what can honestly be measured

This is where the user's warning bites, and the research is unambiguous.

**What Mjolnir will not do:**

- It will not report MTBC members from a taxonomic classifier (see §6).
- It will not present Kraken2 output with a standard/PlusPF index as a
  contamination screen for mycobacteria. Measured Kraken2 sensitivity for
  *M. tuberculosis* reads with the standard database is **0.0731 on real Illumina
  data** — ~93% of true target reads are unclassified or misassigned. A
  mycobacterial pangenome database reaches ~0.97. If a capped or standard index
  is all that is available, the report says the screen is uninformative rather
  than reporting a clean result.
- It will not leave `--confidence` at Kraken2's 0.0 default and print a tail of
  low-abundance NTM as "co-infections".
- It will not use CheckM/CheckM2 contamination as a same-species mixture
  detector — that metric is a multi-copy marker-gene statistic.
- It will not use ConFindr, which has no mycobacterial scheme and is validated
  only on rMLST databases for three genera.

**What it will measure:**

1. **Heterozygosity at lineage-defining positions (F2/F47)** — minor-allele
   frequency across lineage-defining SNP sets, the operationally validated
   approach for separating mixed infection from cross-contamination via batch
   patterns.
2. **Genome-wide heterozygous-SNP fraction**, with MixInfect-style filters
   (Q ≥ 20, DP ≥ 10), reported as a two-tier classification rather than a single
   cutoff.
3. **Mapped fraction, coverage breadth, coverage evenness, GC** — the same
   signals MTBseq's `TBstats` emits, which are its only purity-adjacent output.
4. **Unambiguous-base fraction**, MTBseq's de-facto heterozygosity filter, but
   **reported rather than silently discarded**. MTBseq's 75% majority rule throws
   minority signal away; Mjolnir surfaces it.
5. **Non-target read fraction** by ANI-based assignment, only when a suitable
   mycobacterial reference set is present, and labelled by what it can resolve.

**The headline the report gives is a sample-validity verdict**, not a purity
percentage — because a sample that is 99.84% *M. tuberculosis* still produced 13
false-positive SNPs across 12 genes, and 5% *M. avium* contamination produced
3,325 false-positive variant SNPs. Any gate at 1% or 5% is a coarse instrument
and the report says so.

## 9. Cohort mode

- Joint variant table across samples (MTBseq's `TBjoin` equivalent, without its
  crash on small cohorts and non-MTB references).
- **Masked** pairwise SNP distance. Masking is mandatory: the 508-isolate study
  masked 264,525 loci (~6% of H37Rv) covering repetitive, low-complexity and
  error-prone regions, and counted only SNPs with no other SNP within 12 bases.
  `mask.bed` from tbdb supplies the MTB mask; NTM references get their own
  repeat mask computed at database build time.
- Clustering at a stated threshold, with the threshold's basis printed. Defaults
  follow the 5-SNP / 12-SNP TB conventions; the prior chimaera run on this
  machine used `--distance 6`, so the value is a flag, not a constant.
- **Shared-callable-sites denominator beside every distance.** As in
  `tesseract-ai`'s cgMLST output: 12 differences over 4.1 Mb callable and over
  400 kb callable are not the same statement.

## 10. The agent

Pattern taken directly from `tesseract-ai`:

- `client.py` speaks ollama's native protocol and OpenAI `/v1/chat/completions`,
  detects which the configured host offers, and harvests whichever field the
  server fills. `MJOLNIR_LLM_HOST` points at ollama, vLLM, SGLang or llama.cpp
  and nothing else changes.
- `observation.py` builds what the model may see: finished checks, thresholds and
  their sources, gene names, drug calls, metrics. It **raises** if any field
  contains a long nucleotide run.
- `discipline.py` enforces, per sentence: no number absent from the input;
  nothing unmeasured called fine; no contradiction of a rule-derived verdict; no
  "susceptible" where the rule said "no determinant detected"; for a gate, a
  choice from the closed set. A violating answer is discarded, the rule-derived
  summary replaces it, and the output names the reason.
- Playbooks (`mtbc.yaml`, `ntm.yaml`) carry the organism-specific reading:
  which drugs matter, which caveats always apply, what a clinician needs told.

## 11. The PDF

Clinician-first, research annexes behind. `--profile research` reorders it.

- **Page 1** — sample identity; species with confidence; a drug-by-drug S/R
  table with WHO grade, catalogue agreement, and platform caveats; the
  sample-validity verdict; and the single most important sentence the agent
  produced.
- **Page 2** — lineage/sublineage with the barcode evidence, BCG/animal-lineage
  flags, and the contamination panel.
- **Annexes** — every variant with its coordinate, HGVS, depth, allele fraction
  and per-catalogue grade; catalogue disagreements in full; QC metrics against
  their thresholds with sources; cohort distances and clusters; methods, tool
  versions, database versions and checksums.

Aggressive graphics, in the hydra idiom: a drug-grid heatmap, a coverage/depth
strip, an allele-fraction plot at catalogue positions, a cluster dendrogram.
Self-contained, no network at view time. Same content available as HTML.

## 12. Databases and licences

Every database is fetched by `mjolnir db`, never vendored blindly, and the
registry records each one's version, checksum, licence and citation.

| Database | Licence | Redistributable |
|---|---|---|
| WHO catalogue v2 data (xlsx/VCF) | ODC-By v1.0 | yes, with attribution |
| WHO 2nd-edition PDF | CC BY-NC-SA 3.0 IGO | no |
| tbdb (`mutations.csv`, `barcode.bed`, `mask.bed`) | verify at fetch time | resolved before shipping |
| MTBseq ResSeq lists | GPL-3.0 | verify before vendoring |
| MTBseq NTM references | as MTBseq | — |
| H37Rv NC_000962.3 | public | yes |

Mjolnir itself is MIT. If any source's licence proves incompatible with
redistribution, it is fetched at install time rather than vendored, and that is
recorded in the registry.

## 13. Testing

- Unit tests for every catalogue trap in §5.2 — header row, txt/xlsx divergence,
  decoy coordinate column, grade-string form, per-drug grading, MNV `&` splitting.
- Unit tests for the consensus rule, including the WHO-ungraded-but-called-R path
  and the epistasis suppression path.
- Unit tests for discipline: a model answer inventing a number is rejected; a
  concessive clause is not falsely rejected.
- Golden-file test for the PDF's data layer.
- CI: pyflakes, pytest on 3.9/3.11/3.13, CLI smoke tests, conda build.
- **Validation on real data** (deferred to a free machine): 159 Illumina NTM
  samples across 11 *M. chimaera* sites plus 3 novel-species sets, with the known
  outbreak structure as the clustering test; the prior MTBseq chimaera run as the
  comparison baseline; MTB validation samples to be obtained for the TB path.

## 14. Open questions

- Which mycobacterial reference set ships as the ANI database, and its size
  budget. Needs a decision on breadth vs download weight.
- Whether tbdb's licence permits vendoring `mask.bed` and `barcode.bed`, or
  whether they must be fetched.
- NTM repeat-mask generation at database build time — method not yet fixed.
- MTB validation data source, since the drives on hand carry NTM only.
