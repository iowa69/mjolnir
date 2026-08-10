# Validation plan

What Mjolnir will be measured against, how, and what each measurement can and
cannot establish. Nothing in this document is a result. Results land in
`analysis/metrics.json` and are rendered into the README by `tools/sync_readme.py`,
which the test suite re-runs with `--check` so a hand-edited figure fails the build.

Deferred to a free machine — the runs below are CPU-heavy and the box was under
load when the code was written.

## 1. The data on hand

**150 Illumina paired-end *M. chimaera* isolates across 11 sites**, plus 9
isolates in three "novel species" sets. All on `/media/user/DATA`.

| Site | Samples | What it is |
|---|---|---|
| CHIMAERA-OSR | 43 | clinical |
| CHIMAERA-VENETO | 42 | clinical |
| CHIMAERA-PAVIA | 28 | clinical |
| CHIMAERA-ACQUA | 13 | **water** — environmental |
| CHIMAERA-SAMBROGIO | 8 | clinical |
| CHIMAERA-SCAMBIATORI | 6 | **heater-cooler units** — environmental |
| CHIMAERA-R-EMILIA | 4 | clinical |
| CHIMAERA-POLICLINICO | 3 | clinical |
| CHIMAERA-CESENA / MONZINO / SDONATO | 1 each | clinical |
| NUOVA-SPECIE-NEGRAR | 6 | candidate novel species |
| NUOVA-SPECIE-ZURIGO | 2 | candidate novel species |
| NUOVA-SPECIE-MODENA | 1 | candidate novel species |

Also available: `/media/user/WD_BLACK/M_Chimaera_TN/samp/`, a 2022 MTBseq run on
the same organism — its log, its `chim_ref.fasta` (*M. chimaera* ZUERICH-1,
NZ_CP015272.1), its skesa assemblies and its trees. That is the baseline to
reconcile against.

This collection has a property that makes it unusually good for validating a
clustering module: it is a **multi-site investigation with environmental
isolates included**. Water and heater-cooler isolates should cluster with the
clinical isolates they seeded if the outbreak hypothesis holds, and the tool's
job is to reproduce that structure without being told it.

## 2. What each run establishes

### 2.1 Species — the MAC discrimination test

*M. chimaera* against *M. intracellulare* against *M. avium* is a hard call and
the reason §6 of the spec exists. 159 isolates go through `mjolnir run`; the
species call and its ANI are recorded.

- **Measures**: whether the ANI floor and the marker-SNP layer resolve within
  MAC, and how often the honest answer is *cannot resolve below complex*.
- **Cannot establish**: correctness against a gold standard, because there isn't
  one for these isolates beyond the original study's own calls. Disagreements
  are reported as disagreements, not as errors on either side.
- The three NUOVA-SPECIE sets are the negative control: a tool that confidently
  names a species for a candidate *novel* species is wrong in a way that matters.

### 2.2 Clustering — the outbreak structure test

Cohort mode over all 150 chimaera isolates: joint variants against the
*M. chimaera* reference, masked pairwise distances, clustering.

- **Measures**: whether site structure and the environmental-to-clinical links
  emerge; sensitivity of cluster membership to the threshold (5, 6, 12 SNPs —
  the 2022 run used 6); and how much the shared-callable-sites denominator
  varies across pairs, which is the thing a bare distance hides.
- **Cannot establish**: transmission. A distance is not a chain of infection,
  and the report will not say it is.

### 2.3 Contamination — the honest-metric test

Every isolate gets the full contamination panel.

- **Measures**: the distribution of heterozygous-SNP fraction and
  unambiguous-base fraction across a real collection; how many samples the
  sample-validity verdict would gate for outbreak work but pass for resistance
  work — the intended-use distinction in §8 is the claim being tested.
- **Deliberate negative control**: run the Kraken2 standard index
  (`k2_standard_20250714.tar.gz`, already on WD_BLACK) and confirm Mjolnir
  reports `uninformative` rather than a purity figure. A tool that prints a
  clean percentage there has failed this test.

### 2.4 MTBseq reconciliation

Re-run the 2022 chimaera set under `--mtbseq-compat` and compare position by
position against the preserved MTBseq output.

- **Measures**: whether the compat path reproduces MTBseq's allele frequencies
  given its denominator (N and GAP included) and its tie-break (GAP wins), and
  whether the differences outside compat mode are explained entirely by the MAPQ
  floor and ACGT-only depth.
- **Also tests**: that `mjolnir cohort` completes on the 1- and 2-sample cohorts
  where `TBamend` failed in 2022.

### 2.5 The TB path

The drives on hand carry NTM only, so MTB validation needs data that is not
here. Options, in preference order:

1. A published benchmark set with phenotypic DST, so resistance calls can be
   scored against measured phenotype rather than against another tool.
2. Isolates with known lineage, to check the barcode path and the animal/BCG
   calls.
3. Failing both, cross-tool comparison against TB-Profiler — which measures
   agreement, not correctness, and will be labelled as such.

**Open**: which set, and whether it can be obtained without the compute budget
of a full download. To be decided before the TB path claims any accuracy figure.

## 3. Order of operations

1. `mjolnir doctor` on the target machine — establish what is actually installed
   before anything long runs.
2. Two isolates end to end, one clinical and one environmental, PDF included.
   Read the PDFs by eye. This is the point at which design errors are cheap.
3. The Kraken2 negative control (§2.3) — fast, and it fails loudly if the
   refusal logic is wrong.
4. Species over all 159.
5. Cohort over the 150 chimaera isolates.
6. MTBseq reconciliation.
7. TB path, once data is chosen.

Steps 4–6 are the expensive ones and should run unattended.

## 4. What gets published

Only measured numbers, each generated from `analysis/metrics.json`. Where a
measurement establishes agreement rather than correctness, the README says
agreement. Where a claim could not be tested for want of data — the TB
resistance path, ONT for NTM — the README says untested rather than leaving the
reader to assume.
