<p align="center">
  <img src="docs/logo.svg" alt="Mjolnir" width="110" height="110">
</p>

<h1 align="center">Mjolnir</h1>

<p align="center"><em>Mycobacterial Junction and Omics Locus Nucleotide Identification for Resistance</em></p>

<p align="center">
  <strong>Resistance, lineage, species and contamination for the M. tuberculosis
  complex <em>and</em> for non-tuberculous mycobacteria —<br>
  from Illumina reads, nanopore reads or an assembly, in one command,
  ending in a report a clinician can read.</strong>
</p>

<p align="center">
  <a href="#what-it-does">What it does</a> ·
  <a href="#install">Install</a> ·
  <a href="#use-it">Use it</a> ·
  <a href="#what-it-refuses-to-say">What it refuses to say</a> ·
  <a href="#how-resistance-is-called">Resistance</a> ·
  <a href="#contamination">Contamination</a> ·
  <a href="#comparability-with-mtbseq">MTBseq</a>
</p>

---

## Why this exists

Two tools cover this ground, and neither answers the whole question.

[MTBseq](https://github.com/ngs-fzb/MTBseq_source) is a real variant-calling
pipeline built around H37Rv. It ships NTM references — *M. abscessus*,
*M. chimaera*, *M. fortuitum* — so non-tuberculous mycobacteria look supported:
you swap `--ref`. But every interpretation layer above the caller is MTB-only.
The resistance list, the gene categories and the base-calibration set are all
named `MTB_*` and none of them applies.

Here is what that costs, from a real run on the machine this was written on:

```
--ref        M._chimaera_DSM44623_2016-01-28.fasta
--resilist   NONE
--categories NONE
--basecalib  NONE
...
<ERROR> No joint variant file try_joint_cf5_cr5_fr75_ph4_samples2.tab to amend!
```

A *M. chimaera* isolate gets no resistance, no gene categories, no
recalibration and no lineage — and then the joint step fails outright.

[NTMseq](https://github.com/ngs-fzb/NTMseq) does not close that gap. It is a set
of bash starter scripts driving other people's tools, and it screens
contamination with Kraken2 — whose measured sensitivity for *M. tuberculosis*
reads on a standard index is **0.0731 on real Illumina data**. Roughly 93% of
true target reads never get classified at all.

So today, answering *what is this isolate, is it resistant, is it contaminated,
is it part of an outbreak* means running two dissimilar tools, getting a partial
answer from each, and assembling the rest by hand.

**Mjolnir is one tool where NTM is first-class rather than a reference swap.**

## What it does

| | |
|---|---|
| **Species** | ANI-based, with an explicit *cannot resolve below complex* outcome instead of a false-precision species name |
| **Lineage** | MTBC lineage and sublineage from the tbdb SNP barcode — 1,111 SNPs, 126 taxa — including animal lineages La1/La2/La3 and BCG |
| **Resistance** | Consensus across the WHO catalogue v2, MTBseq's ResSeq list and tbdb, with per-drug agreement shown, not hidden |
| **NTM resistance** | `erm(41)` sequevar typing, `rrl` 2058/2059, `rrs` 1408 — with the citation behind every call |
| **Contamination** | Heterozygosity at lineage-defining sites, not a taxonomic classifier read-out |
| **Cohort** | Masked SNP distances with their shared-callable-sites denominator, and threshold clustering |
| **Report** | A clinician-first PDF: drugs on page one, research annexes behind |
| **Interpretation** | A local model that writes the prose over rule-derived verdicts, and never sees a nucleotide |

## Install

Not on bioconda yet — the recipe is written and builds, but the submission has
not been made, so the conda one-liner does not work and is not printed here.
From a checkout:

```bash
git clone https://github.com/iowa69/mjolnir && cd mjolnir
conda env create -f environment.yml     # the external tools
conda activate mjolnir
pip install .

mjolnir doctor            # says exactly what is present and what is missing
mjolnir db list           # every database, its licence and its citation
mjolnir db fetch          # obtain them
```

The conda recipe is in `conda-recipe/` and `conda build conda-recipe` succeeds;
once it is on bioconda this section becomes one line.

`doctor` never dies at the first missing tool. It reports the whole environment
up front, separating what is required from what is optional, so you find out
once rather than one subprocess at a time.

## Use it

```bash
# one isolate, Illumina
mjolnir run -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o out/

# one isolate, nanopore
mjolnir run --ont reads.fastq.gz -o out/

# an assembly
mjolnir run --fasta assembly.fasta -o out/

# a cohort: joint variants, masked distances, clusters, cohort PDF
mjolnir cohort reads/ -o out/ --distance 12

# regenerate the report from a finished run
mjolnir report out/run.json --profile research
```

## What it refuses to say

Most of the engineering here is in what the tool declines to claim. Each of
these is a documented failure mode, not a hypothetical.

**It will not name an MTBC member from a taxonomic classifier.** In current NCBI
taxonomy the MTBC members are not at species rank at all — taxid 1765,
*Mycobacterium tuberculosis* variant *bovis*, has rank `no rank` beneath species
*M. tuberculosis*, because these are later heterotypic synonyms separated by
ANI 99.21–99.92%. A Kraken2 row reading `M. bovis 3.2%` is not a species
identification, and Mjolnir cannot be made to print one. Within MTBC, identity
comes from lineage-defining SNPs or it does not come at all.

**It will not report "susceptible".** Absence of a catalogued determinant is
reported as *no resistance determinant detected*, which is a different
statement, and the report keeps the two visually distinct everywhere they
appear. Phenotypic susceptibility is not something a genome can assert.

**It will not present a standard Kraken2 index as a contamination screen.** If
the available index is standard, PlusPF or size-capped, the contamination result
carries the status `uninformative` and the report says the screen could not be
run meaningfully — rather than reporting a clean bill of health obtained from a
tool that missed 93% of the target reads.

**It will not give a SNP distance without its denominator.** Twelve differences
over 4.1 Mb of shared callable sequence and twelve over 400 kb are not the same
statement. The API makes the bare number impossible to obtain.

**It will not describe anything unmeasured as fine.** If a capability was lost —
a missing tool, an absent database, FASTA input with no allele fractions — that
loss reaches the report as a loss.

## How resistance is called

Three catalogues, called independently, then reconciled.

| Source | What it is | Licence |
|---|---|---|
| **WHO catalogue v2** | 48,152 variant–drug rows, 30,699 variants, 65 genes, 15 drugs | ODC-By v1.0 — redistributable with attribution |
| **MTBseq ResSeq** | flat list, no confidence grading | GPL-3.0 |
| **tbdb** | TB-Profiler's library; also the barcode and the mask | verified at fetch time |

**WHO is the anchor**, because it is the only source with a published,
systematically derived grading. Where WHO grades a variant, that grade is the
call. Where WHO does not grade it and another catalogue calls resistance, the
drug is reported as **`R (outside WHO catalogue)`** — surfaced, never silently
dropped, and never presented as equivalent to a Group 1 call. Where the three
disagree, page one carries a disagreement flag and the annex shows all three
side by side.

Table lookup alone would be wrong, so the WHO *rules* are implemented too: the
RRDR rule for rifampicin, loss-of-function rules for `katG`, `pncA`, `Rv0678`,
`pepQ` and the nitroreductases, the four borderline `rpoB` mutations that are
Group 1 by decree, and the epistasis suppressions — `mmpL5` loss abrogating
`Rv0678`, `eis` coding loss abrogating the `eis` promoter. A suppression is
recorded and stated, never applied invisibly.

The catalogue carries traps that turn into wrong clinical answers if missed, and
each has a test: the real header is on row 3 of the spreadsheet; the repo's
`.txt` twin is missing streptomycin and four genes and is refused; the
`genomic position` column is a decoy in 38,884 of 48,152 rows; the grade strings
use a spaced ASCII hyphen and not the en-dash the PDF prints; and grading is per
*(drug, variant)*, so `inhA_c.-154G>A` is Group 1 for isoniazid and Group 2 for
ethionamide and must never be deduplicated to one row.

## Contamination

The honest measurements, not the convenient ones.

A sample that was 99.84% *M. tuberculosis* still produced 13 false-positive SNPs
across 12 genes; 5% *M. avium* produced 3,325. Any gate at 1% or 5% is a coarse
instrument, so the headline is a **sample-validity verdict judged against
intended use** — resistance calling tolerates contamination that outbreak SNP
distances do not — rather than a purity percentage.

What is actually measured: minor-allele frequency at lineage-defining positions
(F2/F47), a genome-wide heterozygous-SNP fraction under MixInfect filters
(Q ≥ 20, DP ≥ 10) reported as two tiers rather than one cutoff, mapped fraction,
coverage breadth and evenness, and the unambiguous-base fraction — which MTBseq
computes and then throws away under a 75% majority rule, and which Mjolnir
surfaces instead.

CheckM contamination is a multi-copy marker-gene statistic and is not used as a
same-species mixture detector. ConFindr has no mycobacterial scheme and is not
used at all.

## Platforms

| | Illumina | ONT | FASTA |
|---|---|---|---|
| Min reads for a variant | ≥ 3 | ≥ 5 | — |
| Major variant | ≥ 90% support | ≥ 90% support | — |
| Depth | 25× target | 25× floor, 10× degraded | — |
| Minor variants | reported | **reported as under-detected** | **unavailable** |

Nanopore gets three caveats printed on its reports, because they are measured:
ONT under-detects minor resistance variants — 26 of 27 Illumina-only minor SNPs
were visible in the ONT pileup but uncalled; `fbiC` tandem-repeat deletions
drive spurious delamanid resistance, 47.2% of all discordant drug
classifications in the reference study, so Mjolnir suppresses that call and says
why; and ~16.6% of ONT major indels are uncorroborated by Illumina, so
indel-driven loss-of-function calls carry a platform caveat.

Those thresholds come from a bioRxiv preprint and a benchmark built on a
pseudo-real truthset. They are the best available and Mjolnir uses them — and
cites them as what they are rather than as settled standards. There is no
published R10.4.1-era ONT validation of NTM species ID or NTM genotypic DST at
all, so ONT NTM resistance calls say so on their face.

## Comparability with MTBseq

Numbers will differ, and the reasons are specific. MTBseq calls variants in pure
Perl — a majority-allele caller over a position table, not GATK and not
bcftools. Its frequency denominator includes N and GAP counts, and GAP wins
every tie. No mapping-quality filter is applied at any stage, and `-A` re-admits
anomalous read pairs. `samtools mpileup` caps depth at 250 and MTBseq never
passes `-d`, so deep samples are downsampled before frequencies are computed.

`mjolnir run --mtbseq-compat` reproduces those thresholds *and* that denominator
and tie-break, so the two tools can be reconciled directly. Outside that flag
Mjolnir uses conventional ACGT depth and a mapping-quality floor, and the report
states which convention produced the numbers.

## The model

Interpretation is written by a local model. It changes no verdict.

Every check is computed in Python from a threshold with a stated source, and
pass/warn/fail is fixed before the model is called. The model receives finished
checks and writes the reading. What it may see is built by one function that
**raises** if any field contains a gene-length run of nucleotides — the rule
that it never sees raw sequence is enforced in code, not by convention.

Answers are validated a sentence at a time: a number absent from the input, an
unmeasured thing called fine, a contradiction of a rule-derived verdict, or
"susceptible" where the rule said "no determinant detected" all cause the answer
to be discarded, the rule-derived summary to replace it, and the report to name
the reason.

Point `MJOLNIR_LLM_HOST` at ollama, vLLM, SGLang or llama.cpp — it speaks
ollama's native protocol and OpenAI `/v1/chat/completions` and detects which one
it is talking to. **If no model is reachable the tool still runs**: every gate
takes its declared default and the report says the interpretation is rule-only.

## Licence

MIT. Every database is fetched by `mjolnir db`, and the registry records each
one's version, checksum, licence and citation. Anything whose licence does not
permit redistribution is fetched at install time rather than vendored.

## Status

Version 0.1.0. 42 modules, 613 unit tests, pyflakes clean, a wheel that installs
and runs from a clean environment.

**It has now been run on real isolates**, and it got the answers right on the
three whose identity is not in question — *M. bovis* BCG, *M. bovis* and H37Rv —
including the `pncA_p.His57Asp` pyrazinamide call that is the *M. bovis*
hallmark, and it separated *M. chimaera* from *M. intracellulare* at 99.40%
against 97.32% ANI. Doing so found five defects that 607 unit tests had not,
four of them now fixed, and one missing feature that is not:

> **No module attaches gene names to called variants.** WHO matches on genomic
> coordinates and works; MTBseq and tbdb match on `<gene>_<hgvs>` and therefore
> match nothing. **The three-catalogue consensus is WHO-only in practice, and
> the NTM `erm(41)`/`rrl`/`rrs` rules cannot fire at all.** The report shows
> those catalogue columns as `--`, "not consulted", rather than as `ND`.

Read [docs/VALIDATION-RESULTS.md](docs/VALIDATION-RESULTS.md) before trusting a
consensus call. What is still untested is listed there too: no phenotypic DST
truth, no ONT data, and only four of the 159 chimaera isolates.

Run the examples to see the output shape without installing a single database:

```bash
PYTHONPATH=src python examples/demo_report.py out/     # a clinical PDF + HTML
PYTHONPATH=src python examples/demo_cohort.py out/     # an outbreak cohort
```

<!-- METRICS:BEGIN -->

## Validation

No measurements yet. This block is generated from `analysis/metrics.json` by
`tools/sync_readme.py`, and the test suite re-runs it with `--check`, so a
figure typed in by hand fails the build rather than reaching a clinician.

Until those runs happen, nothing here claims accuracy — including the
tuberculosis resistance path, which has no phenotypic DST data on hand to be
measured against.

<!-- METRICS:END -->
