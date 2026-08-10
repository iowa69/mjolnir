# Design-Grounding Document: Consolidated MTBseq + NTMseq + LLM Interpretation + PDF Report

**Status:** synthesis of 7 adversarially-verified research briefs (referenced below as **B1** MTBseq internals, **B2** NTMseq/NTM-Profiler, **B3** WHO catalogue v2, **B4** cross-catalogue prior art, **B5** contamination, **B6** ONT, **B7** lineage/clustering).
**Rule applied throughout:** claims are marked `[V]` verified in ≥1 brief against a primary source, `[C]` **contradiction between briefs — resolve before spec freeze**, `[U]` uncertain/unverified.
**Pin targets:** MTBseq_source master `7d543f604c0df7fc67bfda9ff2913021ecd3e835` (NOT the `v1.1.0` tag `14d6b61617371d716ff655440a8608221b26125a`; `$VERSION` string is stale) `[V, B1]`; NTMseq main `8d69fb96` (2026-08-07) `[V, B4]`; TBProfiler `47e6c639`; pathogen-profiler `0bd6ddc5`; tbdb `7fe4364e`; NTM-Profiler `24b830b9` `[V, B4]`.

---

## 1. What MTBseq does, precisely, and what of it must be preserved

### 1.1 Architecture (all `[V]`, B1)

A single Perl program `MTBseq` (35,661 B) + 11 `lib/*.pm` modules (`TBtools.pm` alone is 3,094 lines / 155,403 B). No config file, no workflow engine, no sample sheet.

- **11 `--step` values**, validated at `MTBseq#L184-L195`, dispatched by Perl `goto` at `L336-L346`, with labels falling through in fixed linear order: `TBfull`(L349) → `TBbwa`(354) → `TBrefine`(409) → `TBpile`(450) → `TBlist`(491) → `TBvariants`(532) → `TBstats`(573) → `TBstrains`(613) → `TBjoin`(636) → `TBamend`(701) → `TBgroups`(743).
- **Every step ends with `if($continue == 0 && $step ne 'TBfull') { exit 1; }`** (L405, 446, 487, 528, 569, 609, 632, 697, 739, 780). A single-step run exits **status 1 even on success**. Any wrapper treating non-zero exit as failure will misreport.
- 10 fixed output dirs created in CWD (`MTBseq#L36-L45`, mkdir at L367-376): `Bam/ GATK_Bam/ Mpileup/ Position_Tables/ Called/ Statistics/ Joint/ Amend/ Classification/ Groups/`.
- FASTQ discovery regex `MTBseq#L359`: `/^\w.*R\d+\.f(ast)?q\.gz/`. But `TBbwa.pm#L38` then `die`s unless the remainder matches `(R1|R2).f(ast)?q.gz`. The resume check at `L385` strips only `_R\d\.fastq.gz`, so `.fq.gz` inputs are re-mapped every run. **Two mutually inconsistent filename validators + a third for resume.**

### 1.2 The data model — this is the thing to preserve

**MTBseq's canonical intermediate is a per-position count table, not a VCF.** `TBlist` writes one row per reference base with 21 columns (`TBtools.pm#L551-L553`, verbatim):

```
#Pos Insindex RefBase As Cs Gs Ts Ns GAPs as cs gs ts ns gaps Aqual_20 Cqual_20 G_qual_20 Tqual_20 Nqual_20 GAPqual_20
```

(note the inconsistent `G_qual_20`). Positions absent from the mpileup are emitted all-zero (`L593-L599`). Parallelised with MCE, `chunk_size => 1`, ordered gather via a `.tmp` file (`L303`).

**Everything downstream — variants, statistics, lineage, joint tables, clustering — is derived from this table by pure Perl.** If a reimplementation does not keep an equivalent per-position, per-strand, per-base count table with Q20 sub-counts (Parquet/Arrow/numpy), its outputs are not comparable to MTBseq's. `[V, B1]`

### 1.3 Variant calling — pure Perl majority-allele caller (`TBtools::call_variants`, L2043-L2406) `[V, B1]`

- `$cov = adenosin + cytosin + guanosin + thymin + nucleosin + gaps` (L2087) — **Ns and GAPs are in the frequency denominator.**
- `$maximum = max(...)`; ties resolved by six sequential `if($x == $maximum)` blocks in order **A < C < G < T < N < GAP**, later assignment overwriting (L2171-L2225).
- `type` = `SNP` unless allele eq ref; `Del` when winner is GAP; `Ins` when `insertion_index != 0` (freq recomputed against the *previous* position's coverage, L2260-L2277); `Unc`/allele `U` when cov == 0.
- **Unambiguous call (L2281-L2285) requires ALL of:** `cov_f >= --mincovf` AND `cov_r >= --mincovr` AND `freq >= --minfreq` AND `qual20 >= --minphred20` AND allele matches `/[ACGTacgt]/` or eq `'GAP'`.

**Exact defaults (`MTBseq#L198-L218`):** `--minbqual 13, --mincovf 4, --mincovr 4, --minphred20 4, --minfreq 75, --unambig 95, --window 12, --distance 12, --threads 1`. `--samples/--project/--resilist/--intregions/--categories/--basecalib` default to the literal string `'NONE'`. `[V, B1]`

**Must preserve:** the **per-strand minimum** (`mincovf` AND `mincovr` on the winning allele). Almost no modern caller enforces this; dropping it changes resistance calls at low coverage and in homopolymers `[V, B1]`.

**CLI-name trap:** the flag is `--minphred20` (`Getopt` spec `'minphred20:i'` at `MTBseq#L101`); MANUAL.md documents `--minphred` `[C, B4 vs B5]` — B5 and the MTBseq MANUAL both say `--minphred`; B4 verified the source says `--minphred20`. **Source wins.**

### 1.4 Read processing stack `[V, B1]`

`TBbwa.pm#L70-L129`, verbatim: `bwa index` → `samtools faidx` → `bwa mem -t N -R '@RG\tID:..\tSM:..\tPL:Illumina\tLB:..'` → `samtools view -b` → `samtools sort` → `samtools index` → **`samtools rmdup`** (no `-s`/`-S`) → index → move to `Bam/`.

- `samtools rmdup` is declared obsolete in the samtools 1.6 man page and **does nothing for single-end libraries** as invoked, although the MANUAL claims single-end support.
- **There is no trimming, adapter removal or read QC anywhere.** `grep -rinE 'kraken|mash|trimmomatic|fastqc|cutadapt|fastp|trim'` over `MTBseq`, `lib/`, `*.md` → 0 hits.

`TBrefine.pm`: GATK 3.8 `RealignerTargetCreator` (`--downsample_to_coverage 10000`) + `IndelRealigner` (`--defaultBaseQualities 12 --noOriginalAlignmentTags`) + optional BQSR against `var/res/MTB_Base_Calibration_List.vcf` (374 records, VCFv4.1, custom INFO `REG/SVTYPE/GENE/AMINO/RES`). `TBrefine.pm#L106` contains a verbatim stray `-T` before `--analysis_type PrintReads`.

`TBpile.pm#L83-84`: `samtools mpileup -B -A -x -f REF BAM`. **No MAPQ filter (`-q` absent → min MAPQ 0).** `-B` clears BAQ, `-x` clears smart-overlaps, `-A` re-includes orphans. mpileup's default `rflag_filter` still removes unmapped/secondary/QC-fail/duplicate-flagged reads.

**CORRECTED depth cap:** samtools 1.6 `bam_plcmd.c` L557-565: `if (max_depth * sm->n < 8000) max_depth = 8000 / sm->n;`. With one BAM, **the effective cap is 8000×, not 250×** `[V, B1 — explicit correction of a prior brief]`.

### 1.5 Resistance annotation `[V, B1, B4]`

`var/res/MTB_Resistance_Mediating.txt`: 25 columns, **1,696 data rows** (1,691 SNP / 3 `del` / 2 `ins`), latin-1 encoded (`utf-8` raises `UnicodeDecodeError` at byte `0x98`, offset 7324; md5 `dcc3fa0512ba757cdf0b7ea39f49af27`).

`parse_variant_infos` (`TBtools.pm#L753-L855`) consumes only 0-indexed cols **0,1,2,3,5,7,8,12,21**. It upper-cases the Var. base and applies `tr/ATGC/TACG/` **unless Dir. eq '+'** — a per-base complement (not reverse-complement) because the file stores gene-orientation alleles. **1,195 of 1,696 rows are minus-strand**, so naive loading mismatches the majority of the file.

- Rows whose Antibiotic matches `/phylo/` → `PHYLO` (238 rows); all others → `RESI` (1,458 rows).
- **20 distinct non-phylo Antibiotic label strings**, including one drug *class* (`fluoroquinolones (FQ)`, 22 rows) and three multi-drug single strings.
- Per-gene dominance: `pncA/Rv2043c` 891 rows of 1,458; `ethA` 109; `rpoB` 108; `Rv0678` 89; `katG` 72.
- **Columns 22 (Reference PMID), 23 (High Confidence SNP — 224 'yes' / 1,472 '-') and 24 (Comment) are parsed into `@line` and never read.** Every hit is reported with equal weight. Provenance is a mixed bag: `YADON` 628+547, `WHO2021` 79, `WHO2018` 42, `SONNENKALB` 35, `CRyPTIC` 27, `WALKER` 24, plus ~40 DOIs. **It is not the WHO catalogue.**

**CONFIRMED BUG — catalogue indel annotation is completely dead.** The file writes lowercase `del`/`ins`; `parse_variant_infos` tests `eq 'Del'` (L820) and `eq 'Ins'` (L836). All 5 indel rows are silently dropped: `fbiC` 1305494-1305556 (delamanid), `embA` 4243206 (ethambutol), `fabG1` 1673424 (isoniazid), plus 2 phylo insertions. Even with the case fixed, the insertion path stores a 3-level key that `print_variants`' 2-level lookup (`L1209`) can never match, and indexes `$insertion_string[$i]` from 1 over a 0-based split `[V, B1 and independently B4]`.

`var/res/MTB_Extended_Resistance_Mediating.txt`: 10 columns, 43 lines = 1 header + **42 data rows**; parser uses only cols 1,2,6 and tags every position in `[start,stop]` with an antibiotic string → the `InterestingRegion` output column. `tlyA` has a CDS row only (no upstream); two rows are multi-gene windows (`aftA,embC` and `embA,embB`); `ahpC` and `mshA` appear here but have no SNP rows in the main list. Parser uses `split(/\t+/)`, which collapses consecutive tabs and shifts columns on any empty field `[V, B1]`.

`var/cat/MTB_Gene_Categories.txt`: 3,996 lines, 7 starting `#` (skipped). **Effective counts: nonessential 2,878 / essential 760 / repetitive 273 / epitope 78** (total 3,989). Only `repetitive` has any effect — it excludes positions from the phylogeny filter (`TBtools.pm#L1861`). Three genes (`Rv1818c`, `Rv2654c`, `Rv3873`) are commented out of `repetitive` and are therefore **not** excluded `[V, B1; B7 independently confirms 78 epitope, correcting an earlier 76]`.

**Amino-acid substitution strings** are recomputed in-code (`TBtools.pm#L2327-L2398`) by naive codon-rank arithmetic `int(($gene_pos+2)/3)`, formatted `S450L (TCG/TTG)`. No frameshift/indel handling; non-coding features (rrs, rrl, tRNA) get no substitution string; overlapping ORFs silently keep only the last gene (`parse_annotation` comment at `L736`: *"reminder: overlapping ORFs on same strand are ignored… This should be changed."*) `[V, B1]`.

### 1.6 Reference and annotation `[V, B1]`

Default `M._tuberculosis_H37Rv_2015-11-13`, contig name `M.tuberculosis_H37Rv`, 4,411,532 bp, md5 `57b12ff2773c5fd3a0f879972f176e50` — **byte-identical in sequence to NCBI NC_000962.3** (independently recomputed). The `.dict` contains a developer's home path.

Annotation is a **custom 14-column TSV** (`*_genes.txt`), not GFF/GenBank: `# ID name start stop frame product description function cogcats status_region status_function type region_number function_number`. H37Rv file: 3,906 CDS / 45 tRNA / 3 rRNA = 3,954 rows; **1,968 rows encode reverse strand as start > stop**.

**Three NTM references ship with full FASTA + `_genes.txt`:** `M._abscessus_CIP-104536T_2014-02-03`, `M._chimaera_DSM44623_2016-01-28`, `M._fortuitum_CT6_2016-01-08` `[V, B1, B2]`. **But** resistance/intregions/categories/basecalib are auto-loaded only when the reference filename *string* equals `'M._tuberculosis_H37Rv_2015-11-13.fasta'` (`MTBseq#L226-L231`). NTM provenance/accessions of the three references are `[U]` — headers carry only strain names.

**Single-contig by construction.** `parse_reference`/`parse_fasta` (L206-L241) concatenate all non-`>` lines into one coordinate space; `parse_mpile` assigns `my $chromosome = $fields[0];` at L313 and never uses it. **A multi-contig draft NTM assembly is accepted without error and silently merged into one fake replicon with wrong coordinates and wrong gene assignments.** This is the single biggest trap on the NTM side `[V, B1]`.

### 1.7 Lineage (`TBstrains`) `[V, B1, B4, B7]`

96 positions **hard-coded as Perl literals** in `TBtools::parse_classification` (L1050-L1153): 25 Homolka 2012 + 62 Coll 2014 + 9 Merker 2015 Beijing. There is no external barcode file.

`Classification/Strain_Classification.tab` has 16 columns and reports **three schemes, two of them independent Coll implementations**: `Date, SampleID, LibraryID, FullID, Homolka species, Homolka lineage, Homolka group, Quality, Coll lineage (branch), Coll lineage_name (branch), Coll quality (branch), Coll lineage (easy), Coll lineage_name (easy), Coll quality (easy), Beijing lineage (easy), Beijing quality (easy)`. **Every cell is prefixed with a literal apostrophe** (Excel text-forcing, `TBstrains.pm#L157-187`).

Quality flags (`TBstrains.pm#L132-L155`), repeated per scheme:
```perl
$quality_X = "ugly" unless($allel1 =~ /[AGCTagct]/);
$quality_X = "bad"  unless(($freq1 >= 75) && ($count1 >= 10)); # hard coded
```
Both use `unless`, so they fire only on failure and nothing restores `good` — **`good` genuinely means every barcode position passed** (this refutes an earlier "last position wins" claim). What *is* order-dependent (Perl hash order) is only the bad-vs-ugly label. The regex is unanchored, so `GAP` matches `/[AGCTagct]/` and is never flagged `ugly` `[V, B1]`.

**Cannot call:** L8, L9, M. orygis, BCG, or any La* taxon (`grep -ci 'orygis'` = 0, `grep -ci 'BCG'` = 0) `[V, B7]`.

**Documentation bug worth not propagating:** MTBseq's inline comments transpose L1 and L3 — pos 615938 commented `# L1 East-African-Indian`, pos 3273107 commented `# L3 Indo-Oceanic`. Functional output is correct via `translate_coll2homolka` (`1=>'EAI'`, `3=>'Delhi-CAS'`). Any migration script scraping comments for display strings emits clinically wrong names `[V, B7]`.

### 1.8 Cohort layer: TBjoin / TBamend / TBgroups `[V, B1]`

- **TBjoin** requires `--samples` (SampleID<TAB>LibID) and `--project`. It pre-allocates a hash with 3 keys per genome position (~13.2M keys) before doing anything (`TBjoin.pm#L49-L53`). It builds a union-of-variant-positions scaffold, then re-parses each sample's *full* position table and re-runs `call_variants` with `all_vars=1` restricted to those positions. `TBtools.pm#L1369`: `$allel1 = lc($allel1) unless(exists $strains->{$pos.$index.$id})` — **alleles not in that sample's own variant file are lower-cased.**
- **TBamend** phylogeny filter (`TBtools.pm#L1861`): `($pure_SNP == 1) && ($perc_unambigous >= $unambigous) && ($category ne "repetitive") && ($resistance_gene ne "yes")` — where `resistance_gene` is 'yes' whenever the **gene ID** appears anywhere in the resilist. Swapping in a larger catalogue silently removes more genome from the phylogeny.
- **`--window` filter** (`filter_wlength`, L1885-L1981): flags a position if the *same* strain has another SNP within ±12 bp, **but then drops the position for the entire cohort if ANY strain is flagged** (L1956-L1958) — contradicting the code's own comment.
- **TBgroups distance** = count of positions where `uc(allele)` strings differ. **No missing-data handling**; `-` placeholders and re-derived lower-cased alleles participate. Clustering is agglomerative single-linkage at `<= --distance` (12), scanning for the global minimum with strict `<` over an **unordered Perl hash → nondeterministic tie-breaking**; repeat runs can produce different group memberships.

**Filename-encoded parameters are the ONLY provenance/idempotency mechanism** (`_cf<mincovf>_cr<mincovr>_fr<minfreq>_ph<minphred20>_outmode0<snp><lowfreq>`), and that is exactly what causes:

**CONFIRMED BUG — `MTBseq#L657` and `L660` substitute `$micovf` into the `cr` field** while `TBvariants.pm#L57` writes `_cr$micovr`. Any run with `mincovf != mincovr` aborts at TBjoin with "No files to create joint variant tables!" `[V, B1]`.

**`--all_vars` is dead.** `TBvariants.pm#L29` and `TBstats.pm#L32` shadow it with a local `0`; `TBjoin.pm#L30-31` and `TBstrains.pm#L30` hard-code it to `1`. The outmode string's first digit is always `0`; MANUAL's `outmode100` example describes an unreachable state `[V, B1]`.

**TBstats trap:** `Mapping_and_Variant_Statistics.tab` has 24 columns; `% Mapped Reads` IS ×100, but `% (Any) Total Bases` and `% (Unambiguous) Total Bases` are computed **without ×100** (`TBtools.pm#L1246, L1267`) despite the `%` headers — they are fractions `[V, B1]`.

### 1.9 Environment `[V, B1]`

README/MANUAL: Perl v5.22.1; **Oracle/OpenJDK Java 8 only** — `MTBseq#L167-L173`: `unless ($java_version == 1.8) { die "...Need exatly java 1.8..." }` (typo verbatim); bwa 0.7.17, GATK 3.8 (**user-supplied, not bundled**), picard 2.17.0, samtools 1.6; MCE 1.833, Statistics::Basic 1.6611. Bundled `opt/samtools` links `libncurses.so.5`/`libtinfo.so.5` (absent on modern distros). README: bioconda v1.0.4 broken because picard ≥3 requires newer Java. MANUAL: 20–25 GB RAM, Ubuntu 16.04 LTS, max 8 threads (documented but **not enforced** anywhere in code).

**Maintenance status:** last release v1.1.0 2023-08-02; last commit 7d543f60 2023-08-23. **It predates the WHO v2 catalogue entirely.** 10 open issues corroborate the structural problems (#16 rmdup, #18/#27 GATK/Java, #97 "How to update to the newest WHO mutation catalogue?", #100 same mutation multiple times, #102 TBjoin speed) `[V, B1, B4]`.

### 1.10 Regression fixture `[V, B1]`

`MTBseq_source/test/` ships `Test-20_MTBSeq_nextseq_151bp_R{1,2}.fastq.gz` with a documented expected result: MDR (RMP + INH) plus EMB, SM, PZA; lineage 4.1.2.1. **Caveat:** the expected resistance set depends on the broken indel path never firing, so a corrected reimplementation may legitimately report *more*.

### 1.11 What MUST be preserved

1. The 21-column position table as canonical intermediate (or accept non-comparability).
2. Per-strand minimums as a first-class filter.
3. The single parameter object `(mincovf, mincovr, minphred20, minfreq)` applied consistently across call/parse/join/amend — but expressed as a **machine-readable run manifest**, not encoded in filenames.
4. Explicit decisions (documented as divergences) on: MAPQ filtering, depth cap (8000), N/GAP in the coverage denominator, the A<C<G<T<N<GAP tie-break, duplicate handling, and how missing data contributes to SNP distance.
5. The per-sample / cohort split: resistance + lineage + QC are per-sample; clustering + phylogeny FASTA require the joint pass.
6. The bundled test dataset as an end-to-end fixture.

---

## 2. What NTMseq does, precisely, and what of it must be preserved

### 2.1 Architecture (all `[V]`, B2)

`scripts/starter_NTMseq.sh`, `VERSION="2.0.0"`, **750 lines** — a bash orchestrator that shells out to per-tool `starter_*.sh` wrappers, each activating its own conda env. Startup order: arg parse (L247) → source config (L290) → required-var check (L296) → derived vars (L318) → validation (L326) → filetype/suffix detection (L349) → Yes/No validation (L372) → pipeline (L394).

`run_module()` (L106) writes `$LOGDIR/<safe_label>.log`, appends `label\tstart\tend\truntime_min` to `$PATH_output/time.txt`, and on non-zero exit prints `WARNING: $label failed with exit code $status` + `Pipeline will continue with remaining modules.` **Module failure is non-fatal.**

**REFUTED:** there is no `failed_modules.txt`. Failure is recorded only as `STATUS: FAILED` inside that module's own log. **A consolidator cannot detect NTMseq failures from a manifest — it must parse every `logs/*.log`** `[V, B2, explicit correction]`.

Module order in code: multiqc → subsampling → fastp → kraken2 (+perl parser) → NTMprofiler → *[SUBSPECIES block commented out]* → MLST(SRST2) → assembly(shovill) → mash fastQ → mash fastA → fastANI → fastANI all-vs-all → plasmidspades(+Platon) → SRST2 customDB → platon → amrfinder → abricate.

### 2.2 Configuration `[V, B2]`

`config/NTMseq.config` (`#version=v1.0.0`) is a **sourced bash file** — no schema, no CLI other than the config path. **16** `Do_*` switches (the 17th, `Do_subspecies`, and its `validate_yes_no` at L380 are commented out). Required vars (L299-308): `PATH_scripts, PATH_fastQ, PATH_output, species, genome_size, cpu, PATH_tmp, ass`. Derived: `PATH_fastA="${PATH_output}/Assemblies/FinalAssemblies"`; run ID `set="$(date +%Y%m%d_%H%M)"`.

**Sample naming, verbatim from the config:** *"Output folders and result files are named using only the first part of the FASTQ filename (before the first `_`)."* `detect_read_suffixes()` accepts only `_R1_001/_R2_001`, `_R1/_R2`, `_1/_2`.

**This is directly incompatible with MTBseq's `[SampleID]_[LibID]_[*]_[Direction]` scheme** (MTBseq takes the *first two* underscore fields as SampleID and LibID) `[V, B2]`.

### 2.3 Species / subspecies / resistance — fully delegated `[V, B2]`

`starter_NTMprofiler.sh` (`version="2.0.0"`, 447 lines, `platform="illumina"` hardcoded at L32) runs:
```
ntm-profiler profile --read1 $fastq --read2 <derived R2> --platform illumina --dir $PATH_output --threads $cpu --csv --txt -p $SampleName
ntm-profiler collate
```
**NTMseq contains zero marker-gene or ANI-cutoff species logic of its own.** `grep` of the whole tree for `rpoB|hsp65|16S|23S|erm|rrl|rrs|2058|1408` hits only `logo/ntmseq_logo.svg`. `fastANI` writes `summary_fastANi.txt` with no threshold and no species call.

#### NTM-Profiler speciation
`ntm_profiler/cli.py#L396`: `--taxonomic_software` default **`sylph`**, choices `[sourmash, sylph]`. `SylphSketch.classify` runs `sylph profile -t N {ref_db}/* {sketch}` → `TaxonomicHit(ani=round(Adjusted_ANI,2), abundance=Taxonomic_abundance)`. `SourmashSig.classify` filters `intersect_bp<500000` and `f_match<0.1`. `get_species_prediction()` (cli.py L459-499) drops hits with `relative_abundance<1`; hits below `--min_species_relative_abundance` (**default 2.0**) get `relative_abundance=None` and go to `qc_fail_taxa`.

`[C]` **B4/B5 report that NTM-Profiler's README documents *mash* against GTDB** ("If no species is found using this method, mash is run…"), and B4 notes ntm-db ships both a sylph DB and mash sketches. **B2 verified from code that sylph is the default.** README is stale; **design from the code**.

#### ntm-db content
`db/species/accessions.csv` = 3,784 lines (**3,783 genome accessions**, **850 distinct GTDB species labels**); `db/species/sketches/` = 3,794 `.sig` files; `variables.json` = `{"db-schema-version":"4.0.0","type":"species","gtdb-version":"232"}`; sylph submodule `pathogen-profiler/ntm-sylph-db`. `taxonomy.csv` = 33 lines (32 species). Per-species genome counts: M. tuberculosis 500, M. abscessus 500, M. avium 324, M. intracellulare 131, M. marinum 114, M. gwanakae 87, M. fortuitum 73, M. smegmatis 69, M. kansasii 42.

**Exactly 7 per-species resistance folders**, of which only 5 have a `variants.csv`:

| Species | variants.csv | barcode.bed | mask.bed | declared drugs |
|---|---|---|---|---|
| M. abscessus | yes | 450 markers on CU458896 | yes | macrolides, amikacin, fluoroquinolones |
| M. avium | yes | 40 markers on CP000479.1 | yes | rifampicin, macrolides |
| M. intracellulare | yes | 162 markers on CP015278.1 | yes | macrolides, amikacin |
| M. fortuitum | yes (see `[C]`) | no | no | — |
| M. leprae | yes (24 rows) | no | no | rifampicin, dapsone, fluoroquinolones |
| M. malmoense | **no** | no | yes | `[]` |
| M. marinum | no (watchlist only) | 20 markers `ulcerans` | no | `ntm-profiler-ignore: true` |

`[C]` **Contradiction on catalogue sizes.** B2: abscessus "28 data rows", fortuitum "variants.csv = 1 data row" (erm(39) `functionally_normal`, literature `10.1093/cid/ciae421`), leprae 24 rows. B4: "line counts: abscessus 28, avium 24, intracellulare 24, leprae 25, **fortuitum 1 (HEADER ONLY — zero data rows)**", design-implication totals "abscessus 27, avium 23, intracellulare 23, leprae 24, fortuitum 0". These are off by one throughout (line-count vs data-row confusion) **and disagree flatly on whether M. fortuitum has an erm(39) rule.** *Resolve by `wc -l` + header check on a pinned commit before the spec asserts NTM drug coverage.* Total is ~97–100 variant rows across 4–5 species either way.

**No folder exists for M. kansasii, M. chelonae, M. gordonae, M. xenopi, M. simiae, M. scrofulaceum or M. ulcerans.**

#### erm(41) — the clinically decisive, inverted rule `[V, B2, B4]`
`db/Mycobacterium_abscessus/variants.csv` contains `erm(41),functionally_normal,drug_resistance,macrolides,10.1038/s41467-021-25484-9` plus **two** curated LoF alleles (`p.Trp10Arg`, `p.Arg7*`, gene_id `MAB_2297`, drug column **empty** on both LOF rows), commented *"Causes a loss of function in erm(41) and removes the inducible macrolide resistance phenotype"*.

`pathogenprofiler/mutation_db.py::get_functionally_normal_genes` (L106): a gene is `intact=False` if ANY called consequence has type in `(loss_of_function_variant, stop_gained, frameshift_variant, feature_ablation, transcript_ablation)`; **intact** genes are emitted as `Gene(type='functionally_normal')` carrying the drug. `apply_lof_annotation()` (L73) force-sets `csq.type = annotation['so_term']` when `annotation['type']=='loss_of_function'` — this is what turns the missense `p.Trp10Arg` into `loss_of_function_variant`. **`p.Trp10Arg` IS the T28C sequevar**, stated explicitly in `pathogenprofiler/rules.py#L203` doctest: `Variant(chrom='CU458896', pos=2345982, gene_id='MAB_2297', gene_name='erm(41)', nucleotide_change='c.28T>C', protein_change='p.Trp10Arg')`.

**`functionally_normal` = intact gene = inducible macrolide resistance.** A parser reading "normal" as wild-type/susceptible **inverts the clinical call** `[V, B2]`.

**Dead rules path:** pathogen-profiler ships a rules DSL (`inactivates_resistance`, `make_interaction_note`) and NTM-Profiler applies it `if 'rules' in args.conf`, but **no ntm-db `variables.json` has a `rules` key** — `grep` returns nothing `[V, B2]`.

#### Subspecies barcodes `[V, B2]`
`barcode.bed` is 5-column (chrom, start, end, label, allele). Verified label counts: **abscessus 450 markers = exactly 150 each for `subsp. massiliense`, `subsp. bolletii`, `subsp. abscessus`**; **avium 40 = 20 hominissuis + 20 paratuberculosis**; **intracellulare 162 = subsp. chimaera 20, L1 20, L2.1 20, L2.2 20, L3 20, L4 20, L4.1 20, L4.2 14, `HCU_outbreak_strain` 7**; marinum 20 all `ulcerans`.

`pathogenprofiler/barcode.py`: a taxon is skipped if `num_good_sites == 0` (a "good" site requires `target_allele_percent >= 2` AND `all_allele_count >= 5`, **hardcoded in the body despite `min_percent`/`min_allele_count` parameters**), if the IQR of target-allele-% exceeds `iqr_cutoff` (default 15), or if the median frequency is below `freq_cutoff` (default 2). **`BarcodeResult.frequency` is the MEDIAN target allele percent.**

**M. marinum's barcode is dead code:** `update_db` skips the whole folder because `variables.json` sets `"ntm-profiler-ignore": true`, so no resistance DB — and therefore no barcode — is ever built `[V, B2]`.

#### Hard gate on mixtures `[V, B2]`
`ntm_profiler/cli.py#L101-125`: if `number_of_species > 1` → `args.conf = None`; then `if (args.conf is None) or (args.species_only):` → log ("No resistance database found for X" / "Multiple species found, analysis can't continue." / "No species prediction was made") → QC only → `create_species_result()` → **`quit(0)`**.

**`SpeciesResult.result_type == 'Species'` vs `ProfileResult.result_type == 'Profile'`. `dr_variants == []` is indistinguishable from "never assessed" unless you branch on `result_type` first.**

### 2.4 Contamination module `[V, B2]`

`starter_kraken2.sh` (`version="1.2.1"`) L177: `kraken2 --threads $cpu --db $db_kraken2 --output …kraken --report …report --paired $input <R2>` — **`--minimum-base-quality` is commented out; no `--confidence` is set anywhere in the repo.** Then `ktImportTaxonomy`.

`kraken_parse_results.v2.0.0.pl` (internal `$version = "1.0.22"`, 184 lines, "Christian Utpatel, modified by Margo Diricks") emits a fixed 19-column TSV with **no thresholds and no verdict**: `SampleID, Unclassified_perc, Unclassified_reads, Total_reads, Human_perc, Bacteria_perc, Mycobacterium_perc, Species, Species_perc, Species_other_max, Species_other_max_perc, Species_Myco_max, Species_Myco_max_perc, Species_nonMyco_max, Species_nonMyco_max_perc, Genus_max, Genus_max_perc, Genus_max_nonMyco, Genus_max_nonMyco_perc`.

**NEW DEFECT:** L121 computes `Species_Myco_max` with `($name =~ /Mycobacterium/ or $name !~ /Mycolicibacterium/ or $name !~ /Mycobacteroides/ or …)` — **an OR of negations, true for essentially every name**. That column is effectively "highest-abundance species excluding Homo sapiens", **not** a mycobacteria-restricted maximum.

Config default `species="Mycobacteroides abscessus"` **does not match the GTDB label `Mycobacterium abscessus`** that NTM-Profiler reports, so `Species`/`Species_perc` read 0 for the organism of interest unless the user edits the config `[V, B2]`.

### 2.5 Other modules — verified commands and thresholds `[V, B2]`

- **fastp**: `fastp -i -I -o … --cut_tail --cut_tail_window_size 1 --thread N --json --html`, `--phred64` appended if `guess-encoding.py` says so.
- **Subsampling** (`version="2.0.0"`): seqkit to config `cov` (default 100), `seqkit sample -p $Factor -s 11`; writes `<set>_readstats.txt` (11 cols). **`Min_BQ=3` is declared at L30 but its only consumer (`seqtk fqchk -q`) is COMMENTED OUT** — dead code.
- **Assembly**: `shovill --R1 --R2 --outdir --gsize $genome_size --depth $cov --trim --assembler $ass --tmpdir --cpus`; defaults `ass="skesa"`, `genome_size="5.1M"`, `cov=100`.
- **Phylogeny**: `mashtree --mindepth 5 --kmerlength 21 --sketch-size 10000 --outmatrix <set>_distance --numcpus N`; the bootstrap branch prints *"This functionality does not work yet"* (L138).
- **AMRFinderPlus** (`version="1.1.0"`): `amrfinder -n $fasta --threads --name --plus -i 0.5 -c 0.5 -d $db` — **no `--organism`**. `[C]` **B2 explicitly corrects an earlier brief: `coverage_min=0.5` IS the AMRFinderPlus default; only `ident_min=0.5` (vs 0.9) is non-default.** The script itself carries the warning *"It is advised just to run AMRfinder with default values."*
- **abricate**: `--minid 80 --mincov 80` — **both are abricate's own defaults**; config `db_abricate="vfdb"` = **virulence factors, not AMR**.
- **MLST**: SRST2/pubMLST, M. abscessus only; README verbatim: *"only ST types and profiles submitted up to December 2024 will be reported"*. The documented `tseemann/mlst` replacement uses config keys `do_MLST`/`do_MLST_update` **that do not exist** in `NTMseq.config` and a conda env `NTMseq_mlst` **that the installer never creates**.

**Pinned versions** (`installation/Installation_NTMseq.sh`, `VERSION="1.1.0"`): fastqc=0.12.1, multiqc=1.34, seqkit=2.13.0, pigz=2.8, fastp=1.3.2, fastani=1.34, kraken2=2.17.1, krona=2.8.1, srst2=0.2.0, **ntm-profiler=0.8.1** (upstream main is 0.8.2), shovill=1.4.2, mashtree=1.4.6, ncbi-amrfinderplus=4.2.7, abricate=1.4.0, spades=4.2.0, platon=1.7. Tested on Ubuntu 20.04.6, conda 22.11.1 `[V, B2]`.

### 2.6 Shipped data `[V, B2]`

The **only** reference database inside NTMseq is `db/2023_11_03_v2_PLSDB_mycobacteriaceae_208plasmids.fasta` (208 sequences) + `Plasmid_metadata_208.xlsx`. **Eight DB paths must be user-supplied:** `db_kraken2, db_krona, db_custom_SRST2, db_platon, db_AMRfinder, fastANI_ref, db_abricate, db_MLST_SRST2`.

### 2.7 Confirmed defects to not reproduce `[V, B2]`

1. **`starter_NTMseq.sh#L694` calls `starter_platon.sh` (lowercase p); the file is `starter_Platon.sh`.** Silent non-fatal failure on case-sensitive filesystems.
2. **The FASTA branch of `starter_NTMprofiler.sh` (L138-150) is broken and unreachable**: it iterates `$PATH_input/*$filetype` but computes `SampleName` from the stale `$fastq`, and passes the FASTA to `-a` which is `--bam` in NTM-Profiler's CLI (`--fasta` is `-f`, output dir is `--dir/-d`).
3. **Preprocessing outputs are never fed forward.** `-i "$PATH_fastQ"` (raw input) is passed at L403, 420, 438, 455, 478, 522, 540, 562, 631, 673 — i.e. to multiqc, subsampling, fastp, kraken2, NTMprofiler, MLST, assembly, mash-fastQ, plasmidspades and SRST2. `FastQ_subsampled/` and `FastP/` are **written and never read**. Only shovill trims (internally, `--trim`).
4. `starter_fastANI.sh` writes a 4-column header but pastes 5 fields per row — the ANI column is unlabelled and every header after col 2 is shifted.

### 2.8 What MUST be preserved

1. **NTM-Profiler + ntm-db as the required upstream engine** for NTM species, subspecies and point-mutation resistance. Do not rebuild it.
2. The per-species `barcode.bed` semantics (median target-allele %, three filters) and the extra clinically valuable taxa: M. avium hominissuis/paratuberculosis, M. intracellulare subsp. chimaera + L1–L4.2, and the dedicated **`HCU_outbreak_strain`** (heater-cooler unit) barcode — worth a named alert.
3. The `functionally_normal` mechanism, with its inverted logic spelled out in prose.
4. The hard mixture gate (report "not attempted", never "none detected").
5. Kraken2 as the orthogonal read-level contamination check — but re-implemented with a threshold and a verdict.
6. `ntm-profiler update_db --branch --commit` pinning (**it exists**; do not build a shim) and recording `pipeline.species_db_version` / `pipeline.resistance_db_version` from the result JSON. `[C]` B2 explicitly corrects an earlier claim that `NTMprofiler/info.txt` records the DB version — **it does not**; the DB version lives solely in each `<sample>.results.json`.

### 2.9 NTMtools (archive) `[V, B2]`

Explicitly superseded ("For the latest version of NTMseq, please go to…"). Contains the **M. abscessus cgMLST manual** (SeqSphere+, **2,904 loci = 59% of the ATCC 19977 gene set**, good-target definition: same length ±3 triplets, no ambiguities, no frameshifts, ≥90% identity; QC >95% good targets; `Cluster_10/25/250`; subspecies by minimum allele distance to ATCC19977T `NC_010397.1` / JCM 15300 `NZ_AP014547.1` / BD `NZ_AP018436.1`, ties → "Undecided"). **Commercial-only today.** Also contains the M. abscessus marker-gene subspecies SRST2 module (Steindor2019/Minias2020), which is **commented out of NTMseq v2 and whose required marker FASTA is not distributed anywhere.** **NTMtools has no LICENSE file at all.**

---

## 3. Where the two overlap — the shared core

### 3.1 Genuine shared substrate

| Layer | MTBseq | NTMseq | Shared? |
|---|---|---|---|
| Input | paired FASTQ, `[Sample]_[Lib]_*_R{1,2}.fq.gz` | paired FASTQ, name-before-first-`_` | **No — incompatible naming** `[V, B1, B2]` |
| Read QC/trim | none | fastp (output discarded) | **No — both must be rebuilt** |
| Reference | H37Rv + 3 NTM FASTA+`_genes.txt` | none (delegated) | **Yes — MTBseq's `--ref` is a real bridge** |
| Mapping | bwa mem | (shovill internal only) | Partially |
| Variant evidence | 21-col position table | none | MTBseq only |
| Resistance | flat pos+allele catalogue | NTM-Profiler/ntm-db | Different mechanisms, same output need |
| Lineage/species | 96 hardcoded SNPs | sylph ANI + barcode.bed | Different mechanisms, same output need |
| Contamination | none | Kraken2+Krona, no verdict | Neither is adequate |
| Clustering | TBgroups, 12 SNPs | mashtree / fastANI all-vs-all | Different; but see below |
| Output | apostrophe-prefixed Excel TSV | per-module TSV + iTOL | Both unusable as an API |
| Config | CLI only, no config file | sourced bash, no schema | Neither |
| Provenance | filename-encoded params | `time.txt` + per-module `info.txt` | Neither |

### 3.2 The shared core that justifies consolidation

1. **Both are Mycobacteriaceae WGS pipelines over the same organism family, with a common human-facing question:** *what organism, what lineage/subspecies, what resistance, is the sample clean, is it related to anything else.* Both currently answer this as a pile of TSVs with no verdict layer.
2. **Both share people and provenance.** `kraken_parse_results.v2.0.0.pl` authorship line reads "Christian Utpatel, modified by Margo Diricks" — Utpatel is an MTBseq copyright holder and paper author `[V, B2]`.
3. **MTBseq already ships NTM references with full annotation** (abscessus CIP-104536T, chimaera DSM44623, fortuitum CT6). MTBseq's mapping/variant/TBjoin/TBgroups machinery can run against an NTM reference **today**. The missing piece is a per-NTM-species `--resilist`/`--intregions`/`--categories` bundle. **This is a cheaper work item than porting SeqSphere+ cgMLST** `[V, B2]`.
4. **Both need the same five infrastructure layers, which neither has:** a validated sample manifest; a workflow engine with per-process containers; structured per-step status; a machine-readable run manifest (tool versions, DB commit SHAs, thresholds, reference MD5); and a single result schema.
5. **The interpretation and PDF layers are entirely shared.** Resistance evidence classes, mixture gating, QC verdicts, provenance footers, citation handling and the LLM constraint envelope are organism-agnostic.
6. **SNP-distance clustering is shared and already exists on both sides.** `[C]` **B2 explicitly refutes B4's framing that the NTM stack has no SNP-distance capability:** NTM-Profiler ships an experimental `--snp_dist` / `--dist_db_name` (default `ntm-profiler-dists.db`) backed by `pathogenprofiler/snp_dist.py` (sqlite + filelock; `extract_variant_set()` via `bcftools view | bcftools query -f '%POS[\t%GT]\n'`), and `collate --distance_cutoff` (default 10) writes `<outfile>.snp-dist.json`. **NTMseq never invokes it.**

### 3.3 The proposed consolidated shape

```
sample manifest (validated schema)
  → read QC/trim (fastp + FastQC) [NEW — neither tool has a usable one]
  → contamination/purity module [NEW — see §6]
  → organism router: MTBC | NTM | mixed | non-mycobacterial
      MTBC branch:  map(H37Rv) → position table → variant call → catalogue engine (§4) → MTBC lineage (§5.1)
      NTM  branch:  NTM-Profiler (species+subspecies+resistance) [+ optional MTBseq-style mapping to a per-species reference bundle]
  → shared: cohort/clustering pass (optional) → structured Result JSON → deterministic adjudication → LLM narration (constrained) → PDF
```

Per-sample outputs (resistance, species/lineage, QC, contamination) must be producible **with no cohort**; clustering and phylogeny FASTA require the joint pass `[V, B1]`.

---

## 4. Resistance calling: three catalogues, normalisation, and a consensus rule

### 4.1 The three MTBC catalogues, precisely

#### A. WHO catalogue of mutations, 2nd edition `[V, B3, B4]`

- **Only two editions exist.** 1st ed IRIS `dc.date.issued` 2021-06-24; 2nd ed 2023-11-14 (WHO publication page: 15 November 2023). **No 3rd edition published**; a call for data was posted 26 Aug 2024, deadline 15 Oct 2024, target drugs BDQ, LZD, DLM, PMD, D-cycloserine. GTB-tbsequencing org still has exactly 1 repo, last activity 2025-08-28.
- **Correct title:** *"Catalogue of mutations in Mycobacterium tuberculosis complex **and their association with** drug resistance"*. The phrase *"associated with drug resistance"* belongs to the Lancet Microbe companion paper (Walker TM et al. 2022, PMID 35373160, DOI 10.1016/S2666-5247(21)00301-3) — **do not conflate**.
- Corrigendum 17 Jan 2024: *"Due to formatting issues single letters were removed from 148 words randomly throughout the document."*
- **The machine-readable data is NOT on IRIS.** The IRIS ORIGINAL bundle contains only the PDF. Data lives at `github.com/GTB-tbsequencing/mutation-catalogue-2023/Final Result Files/` — exactly 7 files: `Genomic_coordinates_7May2024.vcf.gz`, the "Instruction of use" PDF, `LICENSE.md`, `MCNV.pdf`, `WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt`, `WHO-UCN-TB-2023.7-eng.xlsx`, `WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt`.
- **`WHO-UCN-TB-2023.7-eng.xlsx`**: 2 sheets (`Catalogue_master_file`, `Genomic_coordinates`); **header on row 3**; **48,152 data rows × 114 columns**; coordinates sheet 144,964 data rows × 5 cols (`variant, chromosome, position, reference_nucleotide, alternative_nucleotide`). sha256 `0cd2be5e90d9d43bffb7972d377ea7e2a66e57c77a9243f4fb85c9dd51f680db`, git blob sha1 `c7cc287d5396fd5243404279d828fe4ee4c78823`, 30,943,740 bytes.
- **Exactly 5 FINAL grade strings** (ASCII hyphen, spaces, capital I): `1) Assoc w R` 253; `2) Assoc w R - Interim` 1,130; `3) Uncertain significance` 33,906; `4) Not assoc w R - Interim` 12,379; `5) Not assoc w R` 484. The PDF renders these with an **en-dash and lowercase** (`Assoc w R–interim`) — **never parse grades from the PDF**.
- **15 drugs, 65 genes, 30,699 unique variants, 48,152 (drug,variant) pairs** (unique). Tier 1 = 21,589 / tier 2 = 26,563. Row counts: INH 7286, RIF 7274, EMB 7118, SM 3061, LFX 3047, ETO 2748, MFX 2723, CAP 2654, PZA 2572, KAN 2240, AMK 2190, CFZ 1903, BDQ 1465, DLM 943, LZD 928. **No D-cycloserine, no pretomanid, no rifapentine, no prothionamide.**
- **14 `effect` values**: missense 23065, synonymous 12385, upstream_gene 5240, non_coding_transcript_exon 4028, frameshift 2100, stop_gained 586, inframe_deletion 297, inframe_insertion 168, LoF 91, stop_lost 50, initiator_codon 42, start_lost 41, feature_ablation 40, stop_retained 19.
- **Three components, not one list:** (1) the graded list; (2) *additional grading rules*; (3) *other interpretation criteria*.
  - **Silent rule:** any novel silent variant → Group 4, all genes, all drugs. (`Additional grading criteria applied == 'Silent mutation'` on 12,316 rows.)
  - **RRDR rule:** non-silent variants in rpoB codons **426–452** → Group 2 for RIF (109 rows carry `RRDR`).
  - **LoF rule 11 (verbatim):** premature stop, gene deletion (feature_ablation), frameshift **or start-loss** — **in-frame indels excluded** — in `ethA`(ETO), `gid`(STM), `katG`(INH), `pncA`(PZA), `Rv0678`(BDQ,CFZ), `pepQ`(BDQ,CFZ), `ddn`/`fbiA`/`fbiB`/`fbiC`/`fgd1`/`Rv2983`(DLM), `tlyA`(CAP) → Group 2. (764 rows.)
  - **Borderline rule:** rpoB `Leu430Pro, His445Asn, His445Ser, Ile491Phe` hard-coded to Group 1 (exactly 4 rows flagged `Borderline`).
  - **Epistasis:** *"Rv0678 mutations cannot confer resistance if genetically linked with LoF variants in mmpL5"*; *"eis promoter mutations cannot confer resistance if genetically linked with LoF variants in eis coding region"*.
  - **Level/cross-resistance:** katG → high-level INH; fabG1-inhA → low-level INH **and** INH-ETO cross-resistance; gyrA Gly88Cys/Asp94Asn/Gly/His/Tyr → high-level MFX, all other gyrA/gyrB → low-level MFX; **LFX is explicitly not stratified**. Discrete encodings exist: `BDQ-CFZ cross-resistance` 8, `FQ cross-resistance` 7, `INH-ETO cross-resistance` 3.
  - **Species overrides:** M. canettii → intrinsic PZA resistance (**no plausible pncA marker — must come from species ID**). **M. bovis needs no special case** — its intrinsic PZA resistance comes from `pncA_p.His57Asp`, already a Group 1 row.
- **Drug inheritance, with the correct asymmetry:** RIF→**rifapentine** and ETO→**prothionamide** are FULL Group 1–5 inheritance, stated verbatim in Table 1. **DLM→pretomanid is NOT** — WHO had set no critical concentration, *"Standard analysis for PMD was therefore not possible"*, and only LoF guidance in the six activation genes is offered (+6 rows commented `Confers DLM-PMD cross-resistance`).
- **Allele-frequency basis: ≥75%** (v1 used ≥90%). 75%→25% sensitivity gains, verbatim: RIF **+1.1%**, LFX and MFX **~+4.5% each**, LZD **+0.4%**, INH **+0.5%**; BDQ *"detection of heteroresistance down to at least 25% is important"*. **There is no verified AMK figure — do not use "+2.6%".**
- **rpoB position 761152 / nucleotide 1346 / codon 449 is a documented Illumina artefact site** (footnote **b**, not c). Codon 449 is inside the RRDR, so the RRDR auto-rule would call an artefact RIF-resistant. Deeplex's own mitigation is a 10% floor at that nucleotide.

#### B. tbdb / TB-Profiler `who_v2+` `[V, B4]`

- The compiled DB (`TBProfiler/db/who_v2+/`) is **10 files + `snpeff/`**: `barcode.bed, genes.bed, genome.fasta, genome.gff, mask.bed, mutations.json (27,416,890 B), rules.yml, spoligotype_list.csv, spoligotype_spacers.txt, variables.json`. **No `watchlist.csv`** — watchlist associations are compiled into `genes.bed` column 6 and looked up **by gene NAME** (`gene2drugs.get(var.gene_name, [])`, `reformat.py#L163`).
- `mutations.json`: **73 locus tags, 31,296 distinct (gene,variant) keys, 59,405 annotations, 50,670 distinct genome positions.** `variables.json`: `db-schema-version 2.1.0`, `tb-profiler-version ">=6.6.0,<7.0.0"`, 20 drugs, `version.commit 7fe4364e`, `status "modified"`.
- **The source of truth is TWO files on the `who_v2+` branch:** `mutations.csv` (58,409 rows, 100% `source='WHO catalogue v2'`) **and `additional_mutations.csv` (996 rows, 100% `source='tbdb'`, 100% `type='drug_resistance'`)**. The second file **does not exist on master** (created by commit 7fe4364e, "split mutations into two files"). **A loader reading only `mutations.csv` gets zero tbdb-curated resistance calls** — including all 9 `mmpR5` BDQ/CFZ rows and all cycloserine (15) and para-aminosalicylic_acid (44) coverage.
- **Only `type == 'drug_resistance'` produces a call** (`models.py#L125`, `L310`). Cross-tab: WHO/drug_resistance/`Assoc w R` 317; WHO/drug_resistance/`Assoc w R - Interim` 1,593; WHO/who_confidence/{Not assoc w R 566, Not assoc w R - Interim 15,464, Uncertain significance 40,469}; tbdb/drug_resistance/{'' 644, Not assoc w R - Interim 3, Uncertain significance 349}. **Totals: 2,906 actionable / 56,499 informational; 352 actionable annotations carry a WHO grade of Uncertain or Not-assoc-Interim.**
- **`source` is NOT provenance.** 10,105 rows labelled `source='WHO catalogue v2'` were never graded by WHO: **rifapentine 7,303 rows is a proven exact clone of rifampicin** (identical (Gene,Mutation) keysets, identical type+confidence pairs), **prothionamide 2,772 an exact clone of ethionamide**, and **pretomanid 30 rows** = 6 genes × 5 LoF SO terms with the DLM-PMD comment.
  - `[C]` **B4 frames all 10,105 as "fabricating WHO authority"; B3 verifies that WHO's Table 1 explicitly states RIF classifications also apply to RPT and ETO classifications also apply to PTO.** *Resolution:* the RPT and PTO grades are **legitimate by WHO's own stated inheritance rule** but are **not WHO-graded rows**; the pretomanid 30 are the genuine overreach. The design must attribute them as "WHO rule-derived (RIF→RPT)" and "TB-Profiler inference (DLM LoF→PMD)" respectively — two different labels, not one.
- **Matching** (`mutation_db.py::get_annotation`, L81-104): returns `[]` if `csq.gene_id` not in DB; then **unions** annotations for **three** HGVS keys `(gene_id, nucleotide_change)`, `(gene_id, protein_change)`, `(gene_id, sequence_hgvs)`, then each SO term in `csq.type.split('&')`, then `check_for_so_wildcard` regex `f"{t}_([pcn]).(\d+)_(\d+)"`. **All hits accumulate into a `DictSet` — not first-match-wins.**
- `genome_positions` is **codon-level and carries no ref/alt**: `Rv0667 p.Ser450Leu → [761154,761155,761156]`; reverse lookup of 761155 returns **12 entries**. **Never join on position alone.**
- **`rules.yml` in the shipped `who_v2+` DB has FOUR rules**, not six: `mmpR5_rule` (epistasis), `eis_rule` (epistasis), `canettii_pza_rule`, `rpoB1346Rule`. `rpt_rule`/`eto_rule` (CrossResistanceRule) exist **only on tbdb master** — in `who_v2+` the cross-resistance is pre-baked as cloned rows. **Implementing both double-counts.**
- **`--implement_rules` is DEAD CODE.** It is declared once (`tb-profiler:517`) and never read. The dispatch loop (`tb-profiler:155-159`, `195-199`) runs every rule unconditionally. **Rules are ON by default.** Epistasis semantics: `if source_vars_total_freq >= source_inactivation_freq_cutoff: inactivate_drug_resistance(target_vars)` — **unconditional on target frequency**; `target_escape_freq_cutoff` (rules.yml 50, code default 10) only appends a note.
- **Undocumented `SetConfidence` plugin** (`tbprofiler/rules.py`) fabricates WHO-style confidences for variants absent from the catalogue: `who_confidence` = `'Not assoc w R - Interim'` if synonymous else `'Uncertain significance'`, with `comment: 'Not found in WHO catalogue'`. Runs unconditionally (`tb-profiler:144-145`). **Filter on that comment before reporting anything as a WHO grading.**
- **fabG1 re-basing:** `Rv1483` is absent from `mutations.json`. fabG1 promoter variants are stored under `Rv1484`/inhA coordinates. **The offset is 762** (inhA CDS 1674202 − fabG1 CDS 1673440 = 762 = 777 − 15), **not 777**. `1673425 → Rv1484 c.-777C>T`, comment `Alias fabG1_c.-15C>T`. **Exactly 14 rows contain 'Alias'; `original_mutation` does NOT preserve the WHO name.**
- **Abstract-term expansion:** 601 rows on `who_v2+` have `Mutation != original_mutation` — `LoF` 545, `deletion` 46, **`RRDR non-silent` 10**. SO-term counts: feature_ablation 155, transcript_ablation 110, start_lost 109, stop_gained 109, frameshift_variant 109. *(The 530/{LoF,deletion} and 137/98/97/97/97 figures are master-branch numbers.)*
- `get_drtypes` (`reformat.py#L109-137`) is a hard-coded 8-branch function; ciprofloxacin and ofloxacin in its FLQ set are **unreachable**; the `HR-TB` branch has no fluoroquinolone guard.

#### C. MTBseq `MTB_Resistance_Mediating.txt` — see §1.5

**Position overlap:** 827 distinct MTBseq resistance genome start positions; **789 (95.4%) appear as a genome_position in `who_v2+`; 38 do not** — thyA 8, **pncA 8**, alr 7, fbiC 4, ddn 4, rrl 1, Rv0678 1, Rv2670c 1, fbiA 1, cycA 1, fbiB 1, katG 1. `[C]` **B4 explicitly corrects an earlier "46 missing, pncA 16" figure that could not reconcile with 827−789.**
Locus-tag level: 31 MTBseq resistance genes, **5 with no tbdb counterpart** (Rv1483/fabG1, Rv1704c/cycA, Rv2670c, Rvnr01/rrs, Rvnr02/rrl); **47 of tbdb's 73 locus tags have no MTBseq counterpart** (incl. EBG00000313325/rrs, EBG00000313339/rrl, Rv0676c, Rv0677c, Rv1918c/PPE35, Rv1979c, Rv2477c, Rv2680, Rv2681, Rv3236c, Rv3793, Rv3805c, Rv3806c, Rv0001, Rv1258c, Rv2245, Rv3457c). **Three-way consensus is meaningful for ~26 genes only.**

### 4.2 Normalising a mutation across all three

**Canonical join key: `(chrom='NC_000962.3', pos, ref, alt)`, left-normalised.** WHO states this explicitly: normalise with `bcftools norm`; *"matching for these genomic-variant will only be guaranteed if your own list of variants is properly normalized as well"* `[V, B3]`.

Ingestion per source:

| Source | Coordinate carrier | Transform required |
|---|---|---|
| WHO | `Genomic_coordinates` sheet **or** `Genomic_coordinates_7May2024.vcf.gz` — **provably equivalent** (identical 41,358-name set, identical 118,431 coordinate-key set, zero difference either way; TSV is long-format, VCF is collapsed with `&`-joined `graded_variant`) `[V, B3]` | Split INFO on `&`; many-to-many both ways |
| tbdb | **NOT `genome_positions`** (codon-level, allele-free). Re-derive from `(locus_tag, HGVS)` with snpEff against `Mycobacterium_tuberculosis_h37rv_tbprofiler`, or parse TB-Profiler's `results.json` | Ingest **both** `mutations.csv` and `additional_mutations.csv` |
| MTBseq | cols 1 (genome start), 6 (WT base), 7 (Var base), 13 (Dir.), 22 (Antibiotic) | **latin-1 decode**; `tr/ATGC/TACG/` on Var base when Dir != '+'; drop `/phylo/i` rows; 37 rows are MNVs (Number 2/3) keyed on a single start with a 2–3 char allele string; the 5 indel rows are dead in MTBseq itself |

**Display identity is separate from the join key.** Keep `(locus_tag, gene_symbol, HGVS)` plus **every source's verbatim string**. Ship a version-pinned, CI-asserted-bijective alias table with at minimum:
`Rv1483 ↔ fabG1 ↔ (TB-Profiler's inhA/Rv1484 re-basing, offset 762)`; `Rv3919c ↔ gid ↔ gidB` (MTBseq itself is inconsistent: 11 rows `gidB`, 2 rows `gid`); `Rv0678 ↔ mmpR5` (tbdb: 461 rows `Rv0678` from WHO, 9 rows `mmpR5` from tbdb curation; `rules.yml` refers only to `mmpR5`); `rrs ↔ Rvnr01 ↔ EBG00000313325`; `rrl ↔ Rvnr02 ↔ EBG00000313339`; `Rv1704c ↔ cycA`; `Rv2670c` `[V, B4]`.

**Always print both `inhA c.-777C>T` and `fabG1 c.-15C>T`.** Clinicians and MTBseq know only the second; TB-Profiler JSON emits only the first `[V, B4]`.

**Drug vocabulary table** must reconcile: WHO's 15 Capitalised full names; tbdb's 20 lowercase with `para-aminosalicylic_acid` (underscore); MTBseq's 20 label strings with `para-aminosalicylic acid (PAS)` (space), one class `fluoroquinolones (FQ)` (22 rows) that must fan out to **levofloxacin + moxifloxacin only**, and three multi-drug single strings (`amikacin (AMK) kanamycin (KAN) capreomycin (CPR)`, `bedaquiline (BDQ) clofazimine (CFZ)`, `isoniazid (INH) ethionamide (ETH)`); GARC's 3-letter codes `[V, B4]`.

**Read the WHO xlsx by COLUMN POSITION, never by name.** **41 of 114 header names are duplicated** between the `DATASET ALL` and `DATASET WHO` blocks (every statistic from `algorithm_pass` through `Spec_SOLO_ub`, plus `CHANGES vs ver1`). Name-keyed reads silently merge the two datasets `[V, B3]`.
`[C]` **B3 and B4 give different column indices** — B3 (1-based): `drug=1, gene=2, mutation=3, variant=4, tier=5, effect=6, genomic position=7, Additional grading criteria applied=105, FINAL CONFIDENCE GRADING=106, Comment=107, CHANGES vs ver1 (text)=108, Relaxed simulation=109`. B4 (0-based): `104 Additional grading criteria applied, 105 FINAL CONFIDENCE GRADING, 106 Comment, 111 Additional grading`. **These are consistent modulo a 1-offset; pin the base explicitly in the spec and assert header strings at load.**

### 4.3 Proposed consensus rule

**Design axiom: WHO v2 is the sole grading authority. tbdb and MTBseq can add evidence and can never upgrade a WHO grade.** This follows from B3 (WHO is the only source with a graded evidence model), B4 (tbdb's `source` column is unreliable and its `SetConfidence` plugin fabricates grades) and B1 (MTBseq's `High Confidence SNP` column exists but is never read — every hit is equal-weight).

**Stage 0 — gates.** If contamination verdict is FAIL (§6), or organism is not MTBC, or the drug's genes are not callable at the required depth → emit `NOT_ASSESSED(reason)` for that drug and stop. Never emit `S`.

**Stage 1 — per-variant, per-source verdict.** Each source returns one of:
`GRADED(group)` | `PRESENT_UNGRADED` | `NOT_IN_CATALOGUE` | `GENE_NOT_COVERED` | `NOT_REPRESENTABLE`.
- `GENE_NOT_COVERED` = the source's catalogue has no rows for that locus tag (MTBseq: 47 tbdb loci; tbdb: 5 MTBseq loci).
- `NOT_REPRESENTABLE` = the source structurally cannot express this class of evidence (MTBseq for any LoF/frameshift/whole-gene-deletion; MTBseq for its own 5 dead indel rows).
- `PRESENT_UNGRADED` = matched a WHO coordinate but has no `Catalogue_master_file` row — **10,711 of 41,358 coordinate names are in this state**, heaviest in rrl 952, rpoC 783, rpoB 781, embC 585, gyrA 539, rrs 487, embB 401, clpC1 400, embA 354, gyrB 340 `[V, B3]`.

**Stage 2 — WHO engine, three stages in order** (§4.1A): (a) exact coordinate lookup → grade; (b) additional grading rules (silent→4; RRDR 426–452 non-silent→2 for RIF; LoF in the 12 rule genes→2); (c) other interpretation criteria (level, cross-resistance, epistasis, species override).

**Stage 3 — deterministic modifiers, applied in this fixed order.**
1. **Epistasis (suppression).** Rv0678/mmpR5 BDQ+CFZ calls are inactivated if a linked mmpL5 LoF is present; eis-promoter KAN/AMK calls are inactivated if a linked eis-CDS LoF is present. **If phasing cannot be resolved, do not decide — emit `EPISTASIS_UNRESOLVED` and print both branches.**
2. **Additivity.** Multiple genetically linked low-level INH mutations → high-level INH (WHO Table 1, verbatim).
3. **Artefact guard.** Sub-75%-AF calls at NC_000962.3:761152 (rpoB codon 449) on Illumina → suppress the RRDR auto-rule, flag.
4. **Inheritance.** RIF grade → rifapentine, ETO grade → prothionamide, each labelled "WHO rule-derived". DLM LoF in the 6 activation genes → pretomanid **guidance only**, never a graded prediction.
5. **Species override.** M. canettii → intrinsic PZA resistance. (M. bovis handled by normal `pncA_p.His57Asp` lookup.)
6. **Platform guard** (§7). ONT-derived indels in the suspect-locus list → downgrade to `CONFIRMATION_REQUIRED`.

**Stage 4 — drug-level roll-up (worst-case, with provenance).**

| Report state | Trigger |
|---|---|
| `RESISTANCE_PREDICTED` | ≥1 variant at WHO Group 1 or 2 (direct or rule-derived), passing AF/depth, not suppressed by Stage 3 |
| `RESISTANCE_PREDICTED — INTERIM` | best evidence is Group 2 only |
| `UNCERTAIN` | WHO Group 3, or a `PRESENT_UNGRADED` hit, or a tbdb-only / MTBseq-only actionable hit |
| `NO_MARKER_DETECTED` | all catalogued positions for the drug's genes callable, nothing above Group 3 |
| `NOT_ASSESSED` | Stage-0 gate, or coverage failure over the drug's genes |

**Stage 5 — the agreement matrix, not a vote.** For every reported variant, print a 3×2 table: source × (verdict, coverage state). Print the coverage denominators up front: 789/827 (95.4%) of MTBseq resistance positions exist in tbdb; ~26 genes have three-way coverage; everything else is a **single-source call**. Never take a majority vote across sources with different coverage — that manufactures agreement `[V, B4]`.

**Stage 6 — mandated report categories.**
- **352 overrides:** the tbdb-curated actionable annotations that WHO grades Uncertain (349) or Not-assoc-Interim (3) get their own section — these are the highest-value disagreements and the strongest argument against treating TB-Profiler as a WHO proxy.
- **4 coordinate-less Group 1/2 deletions:** `ethA_deletion, fgd1_deletion, katG_deletion, pncA_deletion` — the **complete, exhaustively enumerated** set a coordinate-only engine would miss. Requires a coverage/depth-based gene-ablation caller.
- **`NOT_REPRESENTABLE in MTBseq`** — katG frameshifts, pncA truncations, ethA disruption, Rv0678 indels. Report as a structural limitation, **never as disagreement**.

### 4.4 Failure modes of this rule

1. **Ingest only `mutations.csv`** → zero tbdb-curated calls, silently. `[V, B4]`
2. **Trust tbdb `source`** → attribute 10,105 rows to WHO; assert WHO authority for pretomanid that does not exist. `[V, B4; nuanced by B3 — see `[C]` in §4.1B]`
3. **Implement rules.yml's 4 rules AND the RPT/PTO cross-resistance rules** → double-counted rifapentine/prothionamide. `[V, B4]`
4. **Pass TB-Profiler `who_confidence` annotations through as WHO gradings** → ungraded variants reported as WHO-assessed "Uncertain significance". Filter on `comment == 'Not found in WHO catalogue'`. `[V, B4]`
5. **Read the WHO xlsx by header name** → the ALL and WHO datasets merge across 41 duplicated columns. `[V, B3]`
6. **Use the `genomic position` column as a coordinate** → 38,884 of 48,152 rows contain the literal string `(see "Genomic_coordinates" sheet)`. The rule is exact: an integer appears **iff** effect is `upstream_gene_variant` (5,240) or `non_coding_transcript_exon_variant` (4,028). `[V, B3]`
7. **Use `WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt`** → 40,178 rows, **14 drugs (Streptomycin entirely absent), 61 genes (gid, rpoB, rpoC, rpsL missing), Rifampicin only 2,361 rows vs 7,274.** A tool parsing this TSV **will never call rifampicin resistance.** `[V, B3]`
8. **Join on genome position alone** → 761155 maps to 12 entries; S450L, S450W, S450F, S450Y fabricate agreement. `[V, B4]`
9. **Join on HGVS string alone** → 14 aliased fabG1 rows report the most important low-level INH marker as a TB-Profiler omission. **Hard-coding 777 instead of 762 shifts every re-based promoter variant by 15 bp.** `[V, B4]`
10. **MCNV/phasing.** WHO's `MCNV.pdf` puts the burden on the genotyper: *"all variants distant less than 2bp away will be grouped together, irrespective of whether they fall on the same codon. Ideally, we only need to phase variants that fall on the same codon."* Mis-phasing changes the amino-acid call itself. The coordinate table carries **both** atomic and merged representations (`dnaA_p.Asp3Ala` at POS 8 as both `A>C` and `AT>CA`). `[V, B3]`
11. **Indel coordinates are only guaranteed for observed changes.** Constant-length SNV/MNV missense and nonsense are exhaustively enumerated even if never observed, so absence is informative for SNVs; **it is NOT informative for indels.** `[V, B3]`
12. **`Comment` parsed by equality instead of substring** → MmpL5-conditional rows counted as 112 instead of **127**; low-level rows as 14 instead of **21**. `[V, B3]`
13. **MTBseq loaded as utf-8** → `UnicodeDecodeError` at byte 0x98, offset 7324. **Loaded without the Dir.-conditional complement** → 1,195 of 1,696 rows mismatch. `[V, B4]`
14. **MTBseq MNV rows.** Whether MTBseq's own per-position pileup can ever produce a key matching its 37 multi-nucleotide rows (e.g. rpoB 761139 `CA>TG` → His445Cys) is **`[U]` — the DB-loading side was traced, the calling side was not.** Test against real data before asserting.
15. **Deduplicating by variant name** destroys the drug axis: `inhA_c.-154G>A` is Group 1 for INH and Group 2 for ETO. `[V, B3]`
16. **Surfacing the `Relaxed thresholds simulation` column** (177 rows) as a call — WHO labels these *"not endorsed"*. `[V, B3]`
17. **The `Footnote` column** holds bare letters whose legends live only in graphically-laid-out per-drug PDF tables, **and the letters are reused with different meanings across drugs**. `[U]` — treat as unusable unless hand-transcribed.
18. **The catalogue predicts a binary phenotype at a critical concentration.** There are **no MIC values or breakpoints in any of the 114 columns.** Do not imply MIC-level precision. `[V, B3]`
19. **Partial application is the dominant published implementation error** — Laurent S, Phelan JE, Chindelevitch L, Walker TM, Cirillo DM, Suresh A, Rodwell TC, Miotto P, Köser CU. *Microbiol Spectr* 2025;13(5), DOI 10.1128/spectrum.02157-24, PMID 40172190. **Cite Laurent, not Walker** (Walker is 4th author). The three omitted components: RRDR non-silent, LoF in non-essential genes (worked example: `pncA Ala134fs` must read as PZA-resistant despite not being listed), and epistasis. `[V, B3]`
20. **Swapping in a larger resistance catalogue silently degrades MTBseq-style phylogeny** — TBamend excludes ALL positions inside any gene appearing anywhere in the resilist (keyed by Gene ID), so SNP distances change. `[V, B1]`

### 4.5 NTM resistance — a structurally different model

Three evidence classes must be modelled separately `[V, B2]`:
1. **Point mutations with an E. coli-equivalent column.** rrl A2057/A2058/A2059/G2069/A2082 (macrolides); rrs T1406/A1408/C1409 (amikacin). **Display via the E.coli column but never key on it** — M. intracellulare's `variants.csv` has a copy-paste defect where `n.1388T>A`, `n.1388T>C` and `n.1388T>G` **all** carry `T1406A`; deduplicating loses two alleles. Coordinates differ per species: A2058 = abscessus `n.2270`, avium `n.2279`, intracellulare `n.2266`.
2. **`functionally_normal` gene rules** (erm(41)/M. abscessus, erm(39)/M. fortuitum `[C]` — see §2.3). Inverted logic, must be spelled out in prose.
3. **Acquired genes** from AMRFinderPlus/abricate, with the `-i 0.5, no --organism` caveat bound to each row, and abricate's vfdb default flagged as virulence not AMR.

**Write a per-species `variants.csv` parser, not one schema:** abscessus and avium have 10 columns; intracellulare has 7 (no gene_id/so_term/comment); fortuitum and leprae have 5; marinum's `gene_watchlist.csv` uses a different `Gene,Info` key=value schema entirely `[V, B2]`.

**Declared-vs-actual drug mismatches:** M. avium `variables.json` declares `[rifampicin, macrolides]` while `variants.csv` contains **amikacin** rrs variants and **no rifampicin variant at all**. Since `collate` builds drug columns from the declared list, an amikacin finding may have no column and a rifampicin column will always be empty `[V, B2]`. M. malmoense's `variables.json` has `"species":"Mycobacterium_malmoense"` **with an underscore** while every other folder uses a space — species-string lookup likely misses it `[V, B2]`.

**M. intracellulare resistance is called against a M. chimaera reference** (`snpEff_db Mycobacterium_chimaera_ZUERICH-1`, barcode chrom CP015278.1). Never present that as an M. intracellulare type strain `[V, B2]`.

---

## 5. Species / lineage calling

### 5.1 MTBC lineage

**Barcode sources, ranked by machine-readability:**

| Scheme | Size | Machine-readable form | Notes |
|---|---|---|---|
| tbdb `barcode.bed` | **1,111 rows, 126 taxa, 8 columns** | yes, in-repo | **Recommended source of truth** `[V, B4, B7]` |
| Coll 2014 | 62 SNPs, 7 lineages + 55 sublineages | `fast-lineage-caller/snp_schemes/coll.tsv` (62 rows, 4 cols) | `[V, B7]` |
| Homolka 2012 + Merker 2015 | 25 + 9 | **Perl literals only** (MTBseq) | must be ported by hand `[V, B1]` |
| Napier 2020 | **90 SNPs covering 90 (sub-)lineage groups** (85 M.tb sublineages + 2 M. africanum + 3 animal species), 27 new sublineages | **none — no barcode file in `GaryNapier/tb-lineages`** | `[V, B7 — explicitly corrects a "9 lineages"/"30 new" claim]` |
| Lipworth 2019 / SNP-IT | **13,610 SNPs, 26 taxa** | `snp_schemes/lipworth.tsv` | includes bcg 146, microti 127, pinipedii 296, dassie 319, mungi 595, suricattae 505, canetti 6745, **ghana 185** `[V, B7]` |
| Shitikov & Bespiatykh 2023 | **213 SNPs, 169 lineages + 9 animal species, 5 levels** | **YES — `pip download tblg` → `tblg/data/levels.tsv`, 213 rows, 174 labels** | `[V, B7 — refutes an earlier "not machine-readable" pitfall]` |

**tbdb `barcode.bed` columns 6/7/8 encode lineage → common name → spoligotype family → RD deletion**, so no hand-built lookup table is needed. Verified strings include `lineage2|East-Asian|Beijing|RD105`; `lineage2.2.1.1|East-Asian (Beijing)|Beijing-RD150|RD105;RD207;RD181;RD150`; `lineage4.3|Euro-American (LAM)|mainly-LAM|None`; `lineage5|West-Africa 1|AFRI_2;AFRI_3|RD711`; `La1|M.bovis`; `La1.2.BCG|M.bovis|BCG`; `La2|M.caprae`; `La3|M.orygis`; `M.canetti` `[V, B7]`.

**Support depth is wildly uneven — this must be surfaced per call:**
- **1 SNP:** `lineage3.1`, `lineage6.1`.
- **4 SNPs:** La1, La1.1, La1.2, La1.3, La1.6, La1.7, La1.8, lineage2.2.1.2, lineage6.2.1.
- **5 SNPs:** **La1.2.BCG**, La1.4, La1.5, La1.7.1, La1.7.X-unk4/5, La1.8.1/.2, La1.8.X-unk6, **La2 (M. caprae)**, **La3 (M. orygis)**.
- All other taxa: 10 (the per-taxon cap). `[V, B7]`

**Calling algorithm** (`pathogenprofiler/barcode.py`): median target-allele % per taxon; skip if zero "good" sites (`target_allele_percent >= 2` AND `all_allele_count >= 5`, hardcoded), skip if IQR > 15, skip if median < `freq_cutoff` (2). **Natively returns multiple taxa with frequencies — use this for mixed-infection reporting.**
**Implementation defect to fix on reuse:** `barcode_rows_quantile` uses `sorted(x)[int(len(x)*q)]` — nearest-rank, no interpolation. With n=10 the "IQR" is percentile[2]..percentile[7]; **with the many 1/4/5-SNP taxa the IQR filter is statistically meaningless and may never fire** `[V, B7]`.

**Animal lineages and BCG — require TWO independent signals before printing "BCG":**
(a) La* barcode taxa (5 SNPs for La1.2.BCG); (b) the Lipworth/SNP-IT scheme (bcg 146 SNPs); (c) **RD coverage read directly from the BAM** — RD1 for BCG; RD105/RD207/RD181/RD150/RD142 for Beijing subtypes; RD239/RD750/RD702/RD711/RD724/RD726/RD182 all appear in `barcode.bed` col 8.
**RD1 absence is not perfectly BCG-specific:** Mahairas et al. *J Bacteriol* 1996;178(5):1274-82 (PMID 8631702) established RD1 (9.5 kb) is deleted from all BCG substrains and conserved in virulent M. bovis and M. tuberculosis — but *Microbiol Resour Announc* 14(4), **2025**, DOI 10.1128/mra.00837-24 describes clinical M. bovis strain 3488 (from a cat, GenBank CP139557) with a **14.2 kb deletion encompassing RD1** `[V, B7]`.
**BCG substrains (Pasteur, Danish, Tokyo, Sofia, Russia, Moreau, Glaxo) are NOT resolved** by La1.2.BCG and are absent from `tblg levels.tsv` entirely `[V, B7]`.

**Clinical consequence of a correct M. bovis/BCG call:** intrinsic PZA resistance via `pncA` His57Asp — **98.2% (162/165) of isolates** per Dong et al., *Transbound Emerg Dis* 2024 (PMC12017247). Mechanism: Petrella et al. *PLoS One* 2011 (PMC3025910) — His57 is axially coordinating the active-site metal and *"is naturally replaced by an aspartate residue in the PZA resistant species Mycobacterium bovis"*. `[C]` **B7 explicitly corrects a prior brief that attributed the 98.2% figure to Petrella — Petrella gives no prevalence.** The therapy consequence (drop PZA, extend to ~9 months) is **`[U]`** — the cited ASM chapter (10.1128/microbiolspec.tnmi7-0021-2016) is behind a login wall; **obtain a primary guideline before printing it.**

**Nomenclature is forked, with no single consensus document.** Coll 2014 → Napier 2020 → Thawornwattana 2021 (L2: L2.2.M1–M6, L2.2.A/AA1–AA4/B–E; *Microb Genom* 7(11):000697, PMID 34787541) / Coscolla 2021 (L5, L6, **L9**; *Microb Genom* 7(2):000477) / Zwyer 2021 (La1/La2/La3; *Open Res Europe* 1:100) → Shitikov 2023 (*mSphere*, 213 SNPs).
**`[C]` Direct disagreement between tools:** tbdb says `lineage4.2.1 = Euro-American (TUR)` and `lineage4.2.2 = Euro-American (Ural)`; MTBseq's `specificator_coll_easy` comments say 4.2.1 = Ural (pos 783601) and 4.2.2.1 = TUR (pos 1455780). **Printing both vocabularies without a crosswalk reads as a contradiction** `[V, B7]`.

**Label normalisation is mandatory at load:** `coll.tsv` uses `lineage4**` and `lineage4.9**` (double asterisk = ancestral/wild-type allele) — **there is no clean `lineage4` row** — plus `lineageBOV` and `lineageBOV_AFRI`. `lipworth.tsv` uses `pinipedii` (one n) and `canetti` (one t). tbdb uses `M.canetti`. `tblg levels.tsv` uses `L4` and contains **malformed labels** (`2.2.1.2` at level 3 vs `lin2.2.1.2` at level 5, plus `AA1SA`, `Bmyc3`), only **five** animal-species labels (not the 9 claimed in the paper), and **no BCG label at all** `[V, B7]`.

**Barcode calling must be independent of the variant caller** — call from a direct pileup at barcode positions. See §7 for the ONT-specific reason.

### 5.2 SNP-distance clustering and masking

**Thresholds and their real meaning:**
- Walker et al. *Lancet Infect Dis* 2013 (PMC3556524): clock **0.5 SNPs/genome/year (95% CI 0.3–0.7)**; *"None of 69 epidemiologically linked and two (15%) of 13 possibly epidemiologically linked patients were separated by more than five SNPs"*; expectation of linkage ≤5, no linkage >12. Mean reference coverage **88.5%**, mapped with Stampy to NC_000962.**2**.
- Meehan et al. *EBioMedicine* 2018;37:410-416 (PMC6284411): **5-SNP clusters span a median 10.86 y (95% HPD 0–47.07); 12-SNP clusters a median 23.63 y (95% HPD 0–102.58).** Cohort **309 L4 + 15 L5** from Kinshasa `[V, B7 — corrects "324 L4"]`, analysed with MTBseq, resistance genes excluded per PhyResSE v27.
- Xiao et al. *Microbiol Spectr* 2024 (PMC11302064): ≤5 = definite, ≤15 = probable; **a 12-SNP threshold would have excluded epi-linked cases in 3 of 4 Taiwanese clusters.**

**→ Emit BOTH 5-SNP and 12-SNP clusterings with their time windows. Never present `group_N` as a confirmed transmission chain** `[V, B1, B7]`.

**Masks differ by ~2.7×:**

| Mask | Intervals | bp | % H37Rv | Contig name |
|---|---|---|---|---|
| Marin 2022 RLC | 773 | 177,077 | 4.014% | `NC_000962.3` |
| tbdb `mask.bed` (merged Marin+Modlin, 2025-08-05) | 2,299 | 311,910 | 7.070% | **`Chromosome`** |
| Coscolla exclusion list | 384 | 472,554 | 10.712% | `NC_000962.3` |

The Coscolla breakdown: **PE/PPE 168 regions, 283,344 bp = 6.42%** (not >10%); InsertionSeqs_And_Phages 147 regions 108,663 bp = 2.46%; "Coscolla Repetative Genes" [sic] 69 regions 80,547 bp = 1.83% `[V, B7]`. **Meehan 2018's "> 10% for the PE/PPE genes alone" is wrong — do not propagate.**

Marin et al. (*Bioinformatics* 2022;38(7):1781-1787, plus a published Correction PMC9997699) found **52/168 PE/PPE genes have perfect sequence-uniqueness**, that PE-PGRS and PPE-MPTR caused 45.4% of false-positive calls, that the top 30 regions (65 kb, 1.5% of the genome) contained 89.4% of false positives, and that an **MQ ≥ 40 filter alone gives recall 85.8% / precision 99.1%** — i.e. MQ tuning outperformed blanket region masking. **But mixing MQ tuning with a published mask makes your distances incomparable to everyone else's** `[V, B7]`.

**tbdb's mask has changed 8 times** (2022-10-04 ×3, 2023-03-27 switch to Modlin blind spots DOI 10.1099/mgen.0.000465, 2024-10-07 ×3, 2025-08-05). **A tool that silently updates its mask silently reclassifies transmission clusters.** Version it, record intervals + bp in the report `[V, B7]`.

**Do not copy TB-Profiler's `snp_dists`.** Two verified defects (`tbprofiler/snp_dists.py`, 149 lines): (1) `SELECT sample, diffs, missing FROM variants WHERE lineage=?` keyed on `result.sub_lineage` — clusters straddling a sublineage boundary are invisible; (2) `dist = self.diffs.symmetric_difference(...); dist -= self.missing; dist -= pickle.loads(m)` — subtracting **both** samples' missing sets systematically shrinks distances for low-coverage samples. Default cutoff is 20, not 12 `[V, B4, B7]`.

**Contig-name skew must be normalised at load:** tbdb BEDs use `Chromosome`; farhat-lab BEDs use `NC_000962.3`; Walker 2013 used NC_000962.2; MTBseq bundles its own FASTA with contig `M.tuberculosis_H37Rv`. **Assert the reference MD5 (`57b12ff2773c5fd3a0f879972f176e50`) before any position-keyed lookup** `[V, B1, B7]`.

### 5.3 NTM species and subspecies

Identification covers **~850 GTDB species labels across 3,783 genomes** via sylph/sourmash ANI + abundance. Resistance interpretation exists for **at most 5 species** (§2.3). **Print a species-capability matrix.** For M. kansasii, M. chelonae, M. gordonae, M. xenopi, M. simiae, M. ulcerans, M. scrofulaceum: identification only. **Say "not assessed", never imply susceptibility** `[V, B2]`.

**Nomenclature normalisation is not hypothetical:** ntm-db ships a directory literally named `Mycobacteroides_abscessus` while NTMseq's README writes `M. abscessus` and GTDB/NTM-Profiler emit `Mycobacterium abscessus`. Gupta et al. *Front Microbiol* 2018;9:67 split the genus into emended *Mycobacterium* (Tuberculosis-Simiae clade, all major human pathogens) + *Mycolicibacterium* (Fortuitum-Vaccae), *Mycolicibacter* (Terrae), *Mycolicibacillus* (Triviale), *Mycobacteroides* (Abscessus-Chelonae) `[V, B5]`.

**M. chimaera is formally a subspecies:** LPSN records `Mycobacterium intracellulare subsp. chimaera (Tortoli et al. 2004) Nouioui et al. 2018`. **MAC sub-species discrimination is a within-species problem, not a Kraken2 species-row task** `[V, B5]`. `[C]` **B5 explicitly REFUTES a fabricated quote** attributed to Khieu et al. *Pathogens* 2021;10:879 about M. chimaera being indistinguishable, and refutes an unsourced "~97.7% ANI" figure — **neither appears in that paper**. What the paper does say: Kraken2 had the highest sensitivity/specificity (100%, 98.23%) for mixed-NTM detection; misidentification of three species *"might have been due to limitations of the database…or to the presence of closely related species"*; all four tools agreed on only 120/155 (77.41%) samples.

**MTBC members are not species in NCBI taxonomy.** taxid 1765 (`Mycobacterium tuberculosis variant bovis`) has rank **"No rank"** with parent species *M. tuberculosis*. **Kraken2/Bracken species-level output can never report M. bovis as a species.** Riojas et al. *IJSEM* 2018;68(1):324-332 (PMID 29205127) formally synonymised M. africanum, M. bovis, M. caprae, M. microti, M. pinnipedii into M. tuberculosis: *"dDDH: 91.2-99.2 %, ANI: 99.21-99.92 %"*. **M. orygis is NOT in that list** `[V, B5]`. **`[U]`: only taxid 1765 was individually verified — verify the others before relying on the "no species rank" rule generally.**

---

## 6. Contamination: what can be honestly measured, and what must never be claimed

### 6.1 Three separable questions `[V, B5]`

1. **Is there non-mycobacterial DNA?** → read-fate / taxonomic classification.
2. **Is there more than one mycobacterial genome?** → heterozygosity at phylogenetically informative sites.
3. **Which organism is dominant?** → marker-SNP / k-mer scheme.
**Kraken2 answers only (1).**

### 6.2 What can be honestly measured

**Read-level (must be a closed budget summing to 100%): host / target clade / other bacteria / unclassified.**

- **Do NOT ship the Kraken2 Standard or PlusPF index as the purity screen.** Verified M. tuberculosis read-level sensitivity (*GigaScience* 2024;13:giae010, PMC10993716): Kraken standard **real Illumina 0.0731** (0.073–0.0732), simulated Illumina 0.1634, real Nanopore 0.7114, simulated Nanopore 0.7408. With a custom Mycobacterium DB: **0.9715 / 0.9731 / 0.9773 / 0.9953**. RAM: **66.8 GB** standard vs **8.2 GB** custom Myco DB.
- `--confidence` **default is 0.0** (Kraken2 MANUAL). The manual itself says *"database false positive errors occur in less than 1% of queries, and can be compensated for by use of confidence scoring thresholds."* **Set it explicitly and record it.**
- Kraken2 assigns **all k-mers the LCA of their minimizer** — a documented coarsening, verbatim: *"All k-mers are considered to have the same LCA as their minimizer's database LCA value."*
- **NTMseq sets no confidence threshold anywhere** (grep of README returns zero occurrences of "confidence") and directs users to a **capped** index; the upstream Langmead page states capped indexes are smaller *"at the expense of some sensitivity and accuracy"*. `[C]` **B5 corrects the attribution: that quote is on benlangmead.github.io/aws-indexes/k2, not in the NTMseq README.** Also: NTMseq hardcodes `k2_pluspf_16_GB_20260226.tar.gz` while the Langmead listing renders it `k2_pluspf_16gb_20260226.tar.gz`. **Do not hardcode either — resolve at runtime with a checksum.**
- **Suppress Bracken's downward redistribution for the Mycobacterium clade.** Reads at genus *Mycobacterium* or at "M. tuberculosis complex" are genuinely ambiguous. Bracken's README documents no behaviour for highly similar genomes or for species absent from the database.

**Variant-level MTBC mixture (two published parameterisations — report separately, never averaged):**
- **MixInfect / Sobkowiak et al. *BMC Genomics* 2018;19:613 (PMC6092779):** filters Q≥20, DP≥10; **pe/ppe and known antibiotic-resistance genes excluded**; only samples with >10 het sites considered; **>20 hSNPs → mixed**; **11–19 hSNPs plus >1.5% het/total-SNP proportion → mixed**. Detection limit: reliable only above **~10% minor strain**; at 0.95/0.05 mixtures, *"between 0 and 2 sites in all samples"*. Mixed-infection frequency in 1,963 Malawi genomes: ~10%.
- **F2/F47 / Wyllie et al. *J Clin Microbiol* 2018;56:e00923-18 (PMC6204665):** minor-allele frequency across the top-two vs lowest-47 lineage-defining SNP sets; thresholds >10× and >5× the development medians = **4.7% (F2)** and **0.2% (F47)**. Production proportions over n=1,794: **97% / 2.8% / 0.001%**. Batch/control patterns are **the only verified means of separating laboratory cross-contamination from true mixed infection**; raised F47 across a batch *"likely reflects process failures"*.

**→ Accept a batch/run identifier and ingest negative/no-template controls as first-class CLI inputs from day one.**

**Assembly-level (NTM branch only):** CheckM `contamination` = *"presence of multi-copy marker genes and the expected collocalization"*; `strain heterogeneity` = multi-copy marker pairs exceeding an **amino-acid identity threshold, default 90%** — report both. CheckM2 auto-selects a completeness model by taxonomic novelty via cosine similarity, but uses **a single universal gradient-boost contamination model** regardless — **its contamination number carries no taxon-specific calibration** `[V, B5]`.

**Published warn/fail gate to reuse as prior art** (Bogaerts et al. *J Clin Microbiol* 2021, PMC8316078, validated clinical MTBC WGS workflow):

| Metric | Warn | Fail |
|---|---|---|
| Contamination | 1.00% | 5.00% |
| Median coverage | 20× | 10× |
| Reads mapping to H37Rv | 95% | 90% |
| GC deviation (expected 65.5%) | 2.00% | 4.00% |
| Average read quality | Q30 | Q25 |

Their heterozygosity policy, verbatim: *"Interpretation of mixed mutations is left to the end user."* `[U]` — the exact metric-definition string for the contamination percentage was **not** verified; define your own explicitly.

**Signals already derivable from MTBseq's own data model** `[V, B1, B5]`: `% Mapped Reads` (flagstat), `(Any)` vs `(Unambiguous)` total bases and GC content and coverage, count of positions with intermediate allele frequency straight from the position table, and conflicting lineage barcode alleles. MTBseq has **zero** contamination code (`grep -rin contamin` → 0 hits across `MTBseq`, `lib/*.pm`, `*.md`).

### 6.3 What must NEVER be claimed

1. **"M. bovis: 3.2%" from a Kraken report.** taxid 1765 has no species rank; that row is a strain-level DB artefact `[V, B5]`.
2. **"M. tuberculosis: 96%" ⇒ clean.** It says nothing about whether one or two MTBC strains are present.
3. **"No contamination detected"** from a Standard/PlusPF index — 93% of true target reads land unclassified or at a higher rank; the derived "percent non-target" is a precise-looking nonsense number.
4. **"High purity ⇒ clean variants."** Goig et al. *BMC Biology* 2020;18:24 (PMC7053099): a sample **99.84% MTB** still had *"13 false positive vSNPs in 12 different genes"*; **5% simulated M. avium contamination produced 3,325 false-positive vSNPs and 51 false-negative fSNPs**, reduced by taxonomic filtering to 24 and 9. In a high-depth cohort, 7 of 63 samples changed after filtering, mean 16.9 removed FP vSNPs (range 2–42). *(Study scope: 4,194 Illumina runs from 20 studies, 1,553 MTB + 2,641 other. `[C]` B5 explicitly corrects an earlier "1,553 across 8 studies" and rejects an unverifiable per-cohort breakdown.)*
5. **"A taxonomic filter removes NTM contamination."** M. avium — the most plausible contaminant — defeats it, and *"non-MTB alignments are not only produced in these regions but across the reference genome."*
6. **"Pure" as a negative result.** The honest string is **"no mixture detected above ~10% minor-strain fraction."**
7. **"Mixed infection diagnosed."** The defensible string is **"evidence of more than one MTBC genome in this library."** The F2/F47 authors could not perform confirmatory microbiology on multiple picks.
8. **"A species hit ⇒ physical presence."** Breitwieser et al. *Genome Res* 2019;29(6):954-960 (PMID 31064768): **2,250 genomes are contaminated by human sequence**, primarily from high-copy human repeats *"not adequately represented in the current human reference genome, GRCh38"*, plus **3,437 spurious protein entries in nr/TrEMBL**. Steinegger & Salzberg, *Genome Biol* 2020;21:115: 2,161,746 / 114,035 / 14,148 contaminated sequences in RefSeq / GenBank / NR. `[U]` — the sentence "nearly all contaminants occurred in small contigs in draft genomes" is **not in the abstract and could not be verified**; the Conterminator counts came from abstract text only.
9. **"CheckM/CheckM2 contamination ⇒ no MTBC mixed infection."** No documented statement that either detects same-species contamination; **no published benchmark of either on deliberately mixed MTBC read sets.** This is an evidence gap, not reassurance.
10. **"ConFindr says clean ⇒ reassuring."** ConFindr *"has only been validated using rMLST databases"* and runs automatically only on *Escherichia, Salmonella, Listeria*. **Neither README nor docs mention Mycobacterium.** `[U]` — whether PubMLST rMLST covers Mycobacterium at all, and whether MTBC has any rMLST allelic diversity, is **unresolved**. `[C]` B5 corrects the attribution: the "53 genes", "3 Contaminating SNVs", "5 percent" quotes are on the docs site (olc-bioinformatics.github.io/ConFindr), not the README; readthedocs 404s.
11. **Sketching (sourmash/mash) for MTBC-internal discrimination.** The sourmash FAQ says at k=31 *"99% (or more) of 31-mers are genome, species, or genus specific"* — **it does not document an upper ANI bound**; the `~85%–99%` range in circulation is **unverified**. The justification for not using sketching within MTBC is the **Riojas 99.21–99.92% ANI**, not any sourmash-documented ceiling `[V, B5]`.
12. **"H37Rv generates false NTM hits."** `[U]` — no primary source documents this. The documented adjacent mechanisms are NTM reads mapping across the whole MTB reference, human repeats inside RefSeq bacterial genomes, and LCA collapse.
13. **NTM contamination thresholds copied from MTBC.** `[U]` — **there are no documented per-metric QC thresholds for the NTM branch anywhere.** NTMseq states none; ntm-db documents none; no published NTM equivalent to the Bogaerts table was found. The 65.5% GC expectation is M. tuberculosis-specific. **Derive and validate in-house or mark advisory-only.**

### 6.4 Asymmetry on human reads

Treat "human reads present" and "human reads absent" differently. minimap2 misclassified **2,476 real M. tuberculosis reads (0.01%)** as human in the GigaScience benchmark. Report the human fraction as an input-quality metric and **never let host depletion run silently before purity accounting** `[V, B5]`.

### 6.5 One place ONT genuinely wins

Colpus et al. 2026 found NTM co-infection inside MTB samples in 9 samples (7 M. avium, 1 `M. novom` [rendered truncated in the preprint — **do not expand it**], 1 M. marseillense): *"In all cases, strong support was found by mapping ONT reads, whereas, with Illumina there was strong support in five and weak support in the remaining four."* Their stack — **Hostile v1.1.0 or Deacon v0.7.0 → fastp v0.24.1 → Kraken2 v2.1.3 restricted to Mycobacteriaceae → minimap2 v2.28 against a Mycobacterium multi-FASTA** — is a validated template `[V, B6]`.

---

## 7. Illumina vs ONT vs FASTA

### 7.1 What changes per input type

| Dimension | Illumina | ONT | FASTA (assembly) |
|---|---|---|---|
| Upstream gate | none beyond QC | **basecall model from RG/header; hard-refuse `fast`** | assembler + read platform provenance |
| Aligner | bwa mem (MTBseq), bwa/minimap2/bowtie2/bwa-mem2 (TB-Profiler) | minimap2 `-x map-ont` | n/a |
| Caller | MTBseq Perl majority-allele; TB-Profiler default freebayes | **Clair3 v2 haploid** (recommended) — **TB-Profiler ships NO Clair3 subclass and silently falls back to bcftools** | n/a — variants from alignment |
| AF/depth thresholds | MTBseq 4/4/4/75%; TB-Profiler `--depth 0,10 --af 0,0.1 --strand 0,3` | asymmetric: **≥5 reads and ≥5% AF; major/minor split at 90% AF** (Colpus); or tbpore's min_depth 5 / min_frs 0.90 / min_mq 30 / min_qual 25 | binary presence |
| Minor/subclonal variants | detectable | **systematically under-detected** | **impossible** |
| Barcode calling | full AD vector retained | **pathogen-profiler collapses each site to the single max-AD allele** — minor-allele fraction discarded | presence only |
| Indels | trustworthy | **separate confidence class** | assembler-dependent |
| Gene ablation | coverage-based | coverage-based | contig-break-confounded |

### 7.2 ONT specifics `[V, B6]`

**Basecalling.** Hall MB, Wick RR, Judd LM et al. *eLife* 2024, DOI 10.7554/eLife.98300, PMID 39388235, PMC11466455. Median unfiltered read identity: **duplex sup 99.93% (Q32), duplex hac 99.79% (Q27), simplex sup 99.26% (Q21), simplex hac 98.31% (Q18), simplex fast 94.09% (Q12)**. Verbatim: *"Fast model ONT data has a lower F1 score than Illumina, only achieving parity in the best case for SNPs."*

**Callers.** *"SNP F1 scores of 99.99% are obtained from Clair3 and DeepVariant on sup-basecalled data. For indel calls, Clair3 achieves F1 scores of 99.53% and 99.20% for sup simplex and duplex… DeepVariant scores 99.61% and 99.22%."* Illumina/Snippy reference: median best SNP F1 **99.45%**, indel **95.76%**. *"FreeBayes matches Illumina for indel calls, but BCFtools shows reduced indel accuracy across all models and read types."* **Longshot calls no indels.**

**CRITICAL CAVEAT:** the truthset is **pseudo-real** — variants from a ~99.5%-ANI donor grafted onto the sample's own assembly, with **indels >50 bp and all structural variation removed**. The MTB row is `AMtb_1__202402`, ANI 99.73%, 2,102 SNPs / 95 ins / 84 del. **99.99% is not a clinical-isolate resistance-calling accuracy figure.**

**Depth.** *"Precision and recall decrease as read depth is reduced, notably below 25x."* *"Clair3 or DeepVariant on 10x ONT sup simplex data provides F1 scores consistent with, or better than, full-depth Illumina."* Recommendation: **minimum 25× for clinical/public-health**; 5× duplex sup reaches Illumina SNP parity. `[U]` — the claim that the authors "agreed 25× in peer review" could not be verified; the published text alone supports it.

**The strongest negative result, which must not be omitted:** *PLOS One* 2024, DOI 10.1371/journal.pone.0303938 — R10.4.1 + V14 + SUP + DeepVariant + TBProfiler v6.2.0 nanopore mode on 17 samples: **ONT *"consistently underreported the level of drug resistance"*, 12/17 (70.5%) wrongly pan-susceptible, lineage concordance run-dependent (1/6, 6/6, 1/6)** — because **only 5/17 (29%) libraries exceeded 10× median depth. Depth, not chemistry, is the dominant failure mode.**

**The largest comparison (unreviewed preprint):** Colpus et al., bioRxiv DOI **10.64898**/2026.04.08.717216 (the `10.64898` prefix is genuine — the new openRxiv prefix; **do not "correct" it to 10.1101**). 508 samples SA+Vietnam, over-selected for resistance, **pre-filtered to ≥50× Illumina depth**, 425 retained (≥10× depth AND ≥90% coverage; 75 excluded for insufficient ONT reads). Overall **VME 1.0% (0.6–1.5%), ME 1.7% (1.3–2.2%), unclassified 6.9% (6.3–7.5%)** across 15 drugs, against CLSI tolerances VME ≤1.5% / ME ≤3%. **Samples with an unknown/failed classification on either platform are EXCLUDED from VME/ME.** Reference standard is Illumina, not phenotypic DST.

**ONT-specific failure modes to encode:**
- **Minor-variant under-detection.** *"we looked at the 27 minor SNPs resolved by Illumina only. In 26/27 cases there is evidence of the SNP in the ONT pileup but it is either below the arbitrary five read minimum threshold for calling (n=8) or is not called by Clair3 (n=18) despite there being sufficient reads."* Mitigation stated: ~50× depth. Rv0678 short indels account for 23 of the 29 platform-exclusive minor indels `[C]` — **B6 explicitly corrects a prior "23 minor indels missed by ONT": the 23 are drawn from a pool that is 2 ONT-exclusive + 27 Illumina-exclusive.**
- **A documented ONT FALSE POSITIVE.** *"About half, 60/127 (47.2%), are due to a detected fbiC deletion triggering a resistance call in ONT for delamanid, though we suggest that the nature of the deletion would not in fact give rise to resistance."* **Delamanid is the sole ME outlier:** *"All drugs, except delamanid, had ME rates ≤1.5%."* Corroboration: **TB-Profiler's own test suite subtracts an `('fbiC','c.2565_*117del')` call before asserting expected DR variants** `[V, B6]`.
- **Genome-wide indel asymmetry.** ONT finds 97.6% of Illumina's major indels; **Illumina finds only 83.4% of ONT's** — i.e. **16.6% of ONT major indels are short-read-unconfirmable.**
- **katG homopolymer deletions.** Hall et al. *Lancet Microbe* 2022/2023 (PMC9892011, DOI 10.1016/S2666-5247(22)00301-9): all **four** ONT-vs-Illumina resistance-genotype discordances across 66,537 resistance-conferring positions were **katG 1 bp deletions, three of them consecutive positions within a katG homopolymer**. `[C]` **B6 explicitly refutes a prior claim that no mycobacterial ONT homopolymer evidence exists.** This is R9.4.1-era (upper bound), but a katG/homopolymer flag is evidence-backed.
- **Barcode collapse.** `pathogenprofiler/bam.py` L302-303: `if platform in ("nanopore","pacbio"): caller="bcftools"` — **before caller dispatch, so it also overrides an explicit `--barcode_caller mpileup`**. L333-338: on ONT, `idx = ad.index(max(ad)); d[genotypes[idx]] = ad[idx]` vs the Illumina branch storing AD for every allele. **Per-site minor-allele fraction at barcode positions is discarded on ONT.** `[U]` — code-reading inference, not benchmarked.
- **Lineage/mixture.** 406/425 (95.5%) lineage-concordant in Results text vs **407/425 (95.8%) in the abstract — the preprint is internally inconsistent by one sample.** All 19 discordances were mixtures; heterogeneity-based mixture detection disagreed both ways (21 both, 4 Illumina-only, 18 ONT-only, all ONT-only from South Africa).
- **Masking.** Colpus built a **264,525-locus (~6% of H37Rv)** ONT/Illumina-comparability mask, added 9 positions post-hoc, and applied a no-other-SNP-within-12-bases density filter → mean cross-platform SNP distance **0.13**, 376/382 replicate pairs ≤1 SNP (vs Hall 2023's **0.75 mean**, SD 1.33). **The authors state they could not independently validate the mask.**

**Tooling state.** Clair3 v2.0.2 (25 Jun 2026) — **PyTorch migration in v2.0.0; v1 TensorFlow models (including ONT Rerio TF models) are incompatible.** Haploid recipe, verbatim: `--platform=ont --model_path=... --no_phasing_for_fa --include_all_ctgs --haploid_precise (or --haploid_sensitive) --enable_variant_calling_at_sequence_head_and_tail`. Bacteria model `r1041_e82_400bps_sup_v430_bacteria_finetuned` (*"Fine-tuned on 12 bacterial genomes"* per README — **note the upstream README says 12 while the linked eLife study analysed 14 samples**). Move-table models `*_with_mv` require `--enable_dwell_time` and Dorado `--emit-moves`.

Dorado: **current release v2.1.1 (30 Jul 2026)**, not v2.0.0. **v2.0.0 shipped a new DNA HAC v6.0 model but NO new DNA SUP model** (the only new SUP models are RNA). Independent benchmark (Wick, 11 Jun 2026, 5 bacterial species): hac@v6.0.0 Q18.1 / 11 median assembly errors vs sup@v5.2.0 Q20.6 / 4 errors; recommendation *"stick with sup@v5.2.0 if they can."* This is **in direct tension with ONT's own v2.0.0 note** claiming HAC v6.0 gives *"polished assemblies comparable to SUP v5.2"* — treat the vendor claim as marketing until reproduced.

**Dorado `smallvar`** (new in v2.0.0, replacing `dorado variant`; bug fixes through v2.1.1 for missing alt alleles, incorrect genotypes, overlapping variants) is documented for **HAC v6.0/v5.2 (not SUP)**, its haploid path is a human-oriented `--hemizygous-regions` flag, and **no bacterial validation exists.** Watch item, not a dependency.

**tbpore** (`.config.yaml`, HEAD fb94398): **SNP-only by design** — `bcftools mpileup -x -I -Q 13 -h100 -M10000` + `bcftools call --ploidy 1 -V indels -m`; filters `min_depth 5, min_frs 0.90, min_mq 30, min_qual 25, min_vdb 0.00001`; `minimap2 -a -L --sam-hit-only --secondary=no -x map-ont`.

**TB-Profiler's own nanopore regression test overrides its defaults hard:** `--caller bcftools --af '0.5,0.7' --depth '0,5'` — i.e. a majority-allele-only policy — versus the Illumina defaults `af 0,0.1 / depth 0,10`.

**MTBseq and NTMseq have ZERO nanopore support.** `grep -ril 'nanopore|minion|map-ont|long.read|oxford'` over both trees returns **zero hits**. The ONT path is new build, not port. **NTM-Profiler exposes `--platform nanopore` but its filtering defaults are unchanged from Illumina** (`--depth 0,10 --af 0,0.1 --strand 0,3`, `--caller freebayes`) and it routes into the same collapse code path `[V, B6]`.

**The WHO catalogue itself is Illumina-derived**, stated verbatim by Colpus: *"The dataset used to derive this catalogue entirely used samples sequenced using Illumina."* This is a report-face disclosure, not just a caveat `[V, B6]`.

**Direct-from-sputum** (Saleeb et al., bioRxiv 10.1101/2025.09.23.678181): *"Optimal results (>90% genome covered, mean coverage >45x and >70% genome covered >20x) were obtained from 33.8% of cases… A further 12.6% of samples yielded suboptimal results (15.5%-90.92% at >10x)."* `[C]` **B6 explicitly refutes two prior claims: there is NO ">25% of a cluster's marker SNPs" rescue threshold in the text (what is reported is 30%/40%/87.5%/91.3% recovery in four samples), and the claimed journal publication `Microb Genom 10.1099/mgen.0.001709` does NOT appear on the bioRxiv record — do not cite it.**

**ONT mycobacterial species ID IS validated** (Baker CS, Colpus M et al., bioRxiv 10.64898/2026.02.04.703726, v2 2026-05-08, from primary MGIT culture): species concordance **98.3% (95.8–99.5%)**, all differences from potential mixed infections; SNP agreement mean 0.3, median 0. `[C]` **an unresolved discrepancy: a search snippet reports 95% concordance and mean 1.0 SNPs, while the bioRxiv API abstract for both v1 and v2 says 98.3% and 0.3. Resolve against the PDF before citing.** **NTM genotypic DST on ONT remains unvalidated.**

### 7.3 What the report must say differently

**Illumina.** State depth, per-strand support, and AF for every call. State the mask and its bp. State that heteroresistance below 75% AF is extrapolation beyond the catalogue's evidence base, with WHO's own 25% gain figures (RIF +1.1%, INH +0.5%, LZD +0.4%, LFX/MFX ~+4.5% each; BDQ qualitative). Apply the rpoB 761152 artefact guard.

**ONT.** Everything above, plus a mandatory provenance block (flow cell, kit, basecaller + exact model string, caller + model) and:
- Three-tier depth semantics: **<10× refuse to report resistance; 10–24× "reduced sensitivity, minor variants not assessable"; ≥25× full report; ≥50× before claiming any minor/subclonal detection.**
- **"Minor/subclonal resistance variants were not assessable at this depth on this platform"** — never "absent".
- **"Mixed infection not excluded"** whenever platform is ONT and depth <50×, because the barcode AD vector is collapsed.
- Every ONT indel in pncA / katG / Rv0678 / ethA / gid / fbiC flagged **"ONT indel — orthogonal confirmation recommended"**, printing the homopolymer run length and tandem-repeat status.
- **fbiC-derived delamanid calls suppressed to a flag by default.**
- A statement that the WHO catalogue's evidence base is entirely Illumina.
- Distances labelled as mask-conditioned and platform-conditioned; state whether the mask was the Colpus ONT-comparability mask.

**FASTA.** No allele frequencies exist → **no heteroresistance, no mixture detection, no minor variants, no LoF-vs-low-coverage disambiguation.** The report must state: *"Assembly input: all variants are reported as fixed; subclonal resistance and mixed infection cannot be assessed."* Gene-ablation calls (the 4 coordinate-less WHO deletions) must be re-derived from contig structure and flagged as assembly-dependent. NTM-Profiler applies `filter_low_coverage_genes(..., cutoff=90)` only when `input_type in ('fasta','bam')` `[V, B2]`. **NTMseq's own FASTA branch is broken and unreachable** (§2.7) — do not copy it.

---

## 8. Unresolved questions a designer must decide

**Data-model and caller**
1. Reimplement MTBseq's position table faithfully (comparable outputs, inherits N/GAP-in-denominator, A<C<G<T<N<GAP tie-break, 8000× cap) **or** adopt a modern caller and accept documented non-comparability? A hybrid (modern caller + position-table export for QC) is possible but doubles the surface.
2. MAPQ: MTBseq filters none. Adding a filter changes resistance calls in PE/PPE with no error raised; keeping secondary alignments diverges the other way. **Decide both explicitly.**
3. Duplicate handling: `samtools rmdup` (obsolete, no-op for SE) vs `samtools markdup` vs Picard. This changes allele frequencies. MTBseq's historical single-end results were **never** duplicate-filtered.
4. Drop GATK 3.8 + Java 8 (removes the biggest packaging and licensing obstacle) — but BQSR against the 374-record calibration VCF is lost. **Measure before discarding.**
5. Do MTBseq's 37 MNV catalogue rows ever match at runtime? Untraced.

**Catalogue engine**
6. Column indexing base for the WHO xlsx (§4.2 `[C]`), and whether to pin by position + header assertion or by a hand-maintained schema file.
7. Attribution wording for the 10,105 tbdb rows: "WHO rule-derived" (RPT/PTO, defensible per WHO Table 1) vs "TB-Profiler inference" (PMD). §4.1B `[C]`.
8. Does the tool ingest tbdb at all, or only WHO + MTBseq? Ingesting tbdb buys 996 curated calls (incl. all mmpR5 BDQ/CFZ, cycloserine, PAS) at the cost of unreliable provenance and the `SetConfidence` artefacts.
9. When phasing cannot resolve an epistasis pair, does the report show both branches or default to the un-suppressed (conservative-for-safety) branch? Both are defensible; pick one.
10. Whether to build the pncA/katG/ethA/Rv0678/gid LoF detector as first-class (recommended) or rely on the catalogue.
11. Whether MTBseq's 5 dead indel rows are fixed (calling variants MTBseq never calls) or faithfully reproduced. **These are different tools — say which in the footer.**
12. WHO catalogue swap-ability: a 3rd edition will land eventually. How is the loader versioned and how are historical reports re-derivable?

**Lineage / clustering**
13. Which barcode is authoritative: tbdb 1,111 (broad, includes La*/BCG) vs Shitikov 213 (169 sublineages, but only 5 animal labels, no BCG, malformed strings) vs both side-by-side with a crosswalk.
14. The TUR/Ural disagreement (§5.1 `[C]`) — pick one vocabulary or print both with an explicit note.
15. Cohort-scoped (MTBseq: recompute everything) vs pairwise-deletion distances. Materially changes the 12-SNP threshold.
16. **How do N/gap/missing contribute to SNP distance?** Neither MTBseq nor TB-Profiler documents a convention. This can move a pair across the 5- or 12-SNP boundary.
17. Which mask ships as default (4.01% / 7.07% / 10.71% / ONT-comparability 6%), and whether MQ≥40 tuning is used instead (better accuracy, worse comparability).
18. cgMLST: does the tool attempt the 2,904-locus M. abscessus scheme via chewBBACA/BIGSdb, or does it stop at SNP distances? SeqSphere+ is commercial.

**NTM**
19. §2.3 `[C]`: does M. fortuitum have an erm(39) rule? Determines whether the tool claims 4 or 5 NTM resistance species.
20. Does the tool build per-NTM-species MTBseq-style `--resilist`/`--intregions`/`--categories` bundles (real work, real payoff) or stay purely delegated to ntm-db?
21. Multi-contig NTM references: reject, or implement real contig-aware coordinates? (MTBseq silently concatenates.)
22. NTM QC thresholds: derive in-house or ship advisory-only? None exist anywhere.

**Contamination / platform**
23. Which Kraken2 DB ships (custom Myco 8.2 GB vs standard 66.8 GB), and what `--confidence` value.
24. Batch/control ingestion: mandatory input or optional? F2/F47 batch analysis is the only verified cross-contamination discriminator.
25. Are ONT and Illumina results ever placed in the same distance matrix? If so, which mask, and is a "cross-platform" flag mandatory?
26. Clair3 v2 as a hard dependency (PyTorch, model-version coupling, v1 incompatibility) vs Medaka vs Dorado smallvar behind a flag.

**Product / governance**
27. **Is this RUO or does it target diagnostic use?** NTMseq's README says *"For Research use only. Not for use in diagnostic procedures."* NTM-Profiler's README says *"in alpha testing and should not yet be considered for production use."* Neither has a methods paper (Europe PMC search for the literal string "NTM-Profiler" returns 13 hits, all uses, no methods paper; NTMseq README says *"Paper in progress…"*). **There are no published sensitivity/specificity, LoD or validation figures for sylph/sourmash speciation, barcode subspecies calls, or NTMseq end-to-end.**
28. The LLM's exact envelope: which fields it sees, which sentences it may generate, and — critically — that the epistasis/artefact/suppression branch decisions are **deterministic code**, with the LLM permitted only to narrate and to *downgrade* confidence, never upgrade, and never to convert `Unknown`/`Fail`/`CONFIRMATION_REQUIRED` into a call.
29. Whether the tool ingests legacy MTBseq/NTMseq outputs (requires a dedicated parser for apostrophe-prefixed, ragged multi-row-header Excel TSVs) or is greenfield-only.
30. The M. bovis/BCG therapy statement (§5.1 `[U]`) — obtain a primary guideline or omit it.
31. Whether the WHO catalogue's treatment of pncA in M. bovis/BCG flags H57D as phylogenetic — **unchecked**, and it gates the PZA-suppression rule.
32. What the `CHANGES vs ver1` integer column (0–6) at position 114 means — **no legend found in the xlsx or the PDF. Do not guess.**

---

## 9. Hard constraints: licensing and dependencies

### 9.1 Databases and catalogues to be redistributed

| Asset | Licence | Verified? | Consequence |
|---|---|---|---|
| **WHO catalogue data** (`Final Result Files/*.xlsx`, `.vcf.gz`, `.txt`) | **ODC-By v1.0** — verbatim: *"All data published here, in excel or VCF format, are licensed under the Open Data Commons Attribution License (ODC-By) v1.0."* | `[V, B3]` | **Redistributable, incl. commercially, with attribution. BUNDLE IT.** Include the ODC-By notice + "WHO Global TB Programme, Catalogue of mutations in Mycobacterium tuberculosis complex and their association with drug resistance, 2nd edition (2023)". **ODC-By preamble caveat:** *"this license only governs the rights over the Database, and not the contents of the Database individually"* — WHO's copyright in individual contents is **not** addressed. **Get counsel before commercial distribution.** |
| GTB-tbsequencing repo **code** (root `LICENSE`) | **MIT** (`Copyright (c) 2023 GTB-tbsequencing`) | `[V, B3]` | fine |
| **WHO 2nd-ed PDF** (text, figures, tables) | **CC BY-NC-SA 3.0 IGO** | `[V, B3]` | **DO NOT bundle PDF text or figures.** Quote sparingly with attribution; NC clause blocks commercial reuse. |
| **MTBseq_source** (code + `var/` data) | **GPLv3 "or any later version"** (`LICENSE.md` + `GPL.md`); copyright 2018 Kohl, Koch, Utpatel, De Filippo, Schleusener, Beckert, Cirillo, Niemann. *"The included third party programs are redistributed under their own license."* | `[V, B1, B2]` | GPLv3 is **viral** — vendoring `MTB_Resistance_Mediating.txt` / `MTB_Gene_Categories.txt` / the barcode literals makes the derived work GPLv3. `[U]` whether the `var/` data files are separately licensed — **assume GPLv3.** |
| **NTMseq** | **GPLv3** (full 674-line text) | `[V, B2]` | viral, same as above |
| **NTMtools** | **NO LICENSE FILE ANYWHERE** (full tree listing verified) | `[V, B2]` | **All rights reserved by default. The cgMLST manual, `Reference_genomes_for_classification.xlsx` and publication scripts CANNOT be vendored without written permission from mdiricks@fz-borstel.de.** |
| **NTM-Profiler** | **GPL-3.0-or-later** (`pyproject.toml` `license = "GPL-3.0-or-later"`) | `[V, B2]` | viral if linked |
| **pathogen-profiler** | **LGPLv3** (`LICENSE` = "GNU LESSER GENERAL PUBLIC LICENSE Version 3"); `setup.py`, no pyproject | `[V, B2]` | LGPL — linkable without viral effect if used as a library |
| **ntm-db** | **NO LICENSE FILE** (confirmed by `git ls-tree`, root contains only `.github/`, `.gitignore`, `.gitmodules`, `README.md`, `db/`) | `[V, B2]` | **Do NOT vendor. Pin-and-cite: `ntm-profiler update_db --commit <SHA>` and record `pipeline.species_db_version` / `resistance_db_version` in every report.** No releases exist. |
| **tbdb** (`mutations.csv`, `additional_mutations.csv`, `barcode.bed`, `mask.bed`, `rules.yml`) | **NOT VERIFIED IN ANY BRIEF** | `[U]` | **BLOCKER — must be checked before shipping.** The whole §4 consensus design depends on it. |
| **TBProfiler** (incl. compiled `db/who_v2+/`, `default_template.docx`) | **NOT VERIFIED IN ANY BRIEF** | `[U]` | **BLOCKER.** |
| **Oxford GARC** `parse-who-tb-catalogue-v2` (3,715-row CSV, `expert-rules.csv`) | **NOT VERIFIED** | `[U]` | check before use as a WHO-parsing shortcut |
| **fast-lineage-caller** `snp_schemes/*.tsv` | **NOT VERIFIED** | `[U]` | check |
| **`tblg`** (PyPI, Shitikov 213-SNP `levels.tsv`) | **NOT VERIFIED** | `[U]` | check |
| **farhat-lab masks** (`RLC_Regions.H37Rv.bed`, `200624_CoscollaExcludedGenes.bed`) | **NOT VERIFIED** | `[U]` | check |
| **NTMseq plasmid FASTA** (208 seqs, PLSDB-derived) | repo GPLv3, but **PLSDB's own terms not verified** | `[U]` | check upstream PLSDB terms |
| **Kraken2 index** (Langmead prebuilt, or self-built from RefSeq) | not verified | `[U]` | RefSeq content is public domain; index build terms unchecked |
| **H37Rv NC_000962.3 reference** | NCBI/GenBank — public | `[V, B1]` (sequence identity verified; terms not) | fine |

### 9.2 External tool dependencies

**Licences below are marked `[V]` only where a brief actually verified them. Everything marked `[U]` MUST be confirmed against the shipped LICENSE file before redistribution — do not take the annotation on faith.**

**MTBseq branch:** bwa 0.7.17 `[U]`; samtools 1.6 `[U]`; **GATK 3.8 `[U]` — NOT redistributable, user-supplied, not bundled in `opt/`; upstream issues #18 and #27 confirm this is a live pain point** `[V, B1]`; picard 2.17.0 `[U]` (must be <3 for Java 8); Perl ≥5.22.1; CPAN `MCE` 1.833, `Statistics::Basic` 1.6611 `[U]`. **Hard constraint: `java -version` must report exactly 1.8 or MTBseq dies** `[V, B1]`.

**NTMseq branch (pinned versions verified `[V, B2]`; licences `[U]` throughout):** fastqc 0.12.1, multiqc 1.34, seqkit 2.13.0, pigz 2.8, fastp 1.3.2, fastani 1.34, kraken2 2.17.1, krona 2.8.1, srst2 0.2.0 (+ GNU parallel 20260122), ntm-profiler 0.8.1, shovill 1.4.2, mashtree 1.4.6, ncbi-amrfinderplus 4.2.7, abricate 1.4.0, spades 4.2.0, platon 1.7. **Note: SPAdes ships a non-standard licence with usage restrictions — verify. AMRFinderPlus is NCBI (US Government work) — verify. abricate is commonly GPLv2 — verify.**

**TB-Profiler / pathogen-profiler runtime:** bwa | minimap2 | bowtie2 | bwa-mem2 (mapper choices), freebayes | bcftools | gatk | pilon | lofreq | freebayes-haplotype (caller subclasses — **no Clair3**), samclip, delly (on by default), trimmomatic (on by default), snpEff, FastK | kmc | dsk (kmer counters), bedtools | samtools (coverage). **All licences `[U]`.**

**ONT branch (all `[U]`):** Dorado (**ONT-proprietary EULA — likely the single hardest redistribution constraint on this branch; verify before containerising**), Clair3 v2.x (+ PyTorch), DeepVariant, Medaka, minimap2, rasusa, SeqKit, vcfdist. Optionally Hostile / Deacon for host depletion.

**Contamination branch (all `[U]`):** Kraken2, Bracken, Krona, sourmash, sylph, CheckM/CheckM2, ConFindr (**do not use for mycobacteria — §6.3**), Mash.

### 9.3 Licence-consolidation blockers, summarised

1. **The consolidated tool will be GPLv3** if it vendors any MTBseq or NTMseq code or data. This forecloses proprietary distribution.
2. **NTMtools is all-rights-reserved.** Its cgMLST scheme documentation, reference spreadsheet and publication scripts are off-limits without written permission.
3. **ntm-db has no licence.** Pin-and-cite; never vendor.
4. **tbdb and TBProfiler licences are unverified** and gate the entire cross-catalogue design. **Resolve first.**
5. **GATK 3.8** is non-redistributable and forces exactly Java 8. Dropping it is the single biggest containerisation win — measure the BQSR loss first.
6. **Dorado** is vendor-licensed; the ONT branch may need to require a user-supplied basecaller rather than shipping one.
7. **WHO PDF text is CC BY-NC-SA 3.0 IGO** while WHO catalogue *data* is ODC-By. Do not mix them in the report.

### 9.4 Mandatory report-face disclosures

- "Research use only. Not for use in diagnostic procedures." (NTMseq's own wording.)
- NTM-Profiler is self-declared alpha with **no methods paper and no published sensitivity/specificity** for speciation or barcoding; NTMseq is unpublished ("Paper in progress"); ntm-db is a live git clone with no releases.
- **The WHO catalogue's evidence base is entirely Illumina-derived and reflects a ≥75% allele-frequency threshold.**
- MTBseq's evidence base is **pre-WHO-v2** (frozen since 2023-08-23).
- AMRFinderPlus identity forced to 0.5; abricate defaulted to **vfdb (virulence, not AMR)**; mashtree bootstrap unimplemented; pubMLST STs frozen at December 2024 under SRST2.
- Catalogue provenance as a **checksum, not a version string** — the WHO repo silently rewrote its files twice (Feb 2024 MCNV fix, May 2024 stop-codon fix) with no filename change. **The xlsx filename cited *inside* the WHO PDF is `WHO-UCN-TB-2023.5-eng.xlsx`, which does not match the distributed `WHO-UCN-TB-2023.7-eng.xlsx` — never build a download URL from the PDF text** `[V, B3]`.
- Full provenance block: basecaller + exact model, aligner, caller + model, catalogue version + sha256, barcode file + commit SHA, mask file + interval count + bp, reference accession + MD5, clustering threshold, every filter value, and the DB commit SHAs from `pipeline.*_db_version`.

### 9.5 Operational note on IRIS

DSpace bitstream content endpoints (`server/api/core/bitstreams/<uuid>/content`) work for GET but **return HTTP 500 for HEAD**, and WebFetch is blocked on iris.who.int. Any availability check or cache layer probing with HEAD will conclude the file is gone. **Do not put IRIS on the critical path — the data you need is on GitHub** `[V, B3]`.

---

### Appendix: contradictions requiring resolution before spec freeze

| # | Contradiction | Briefs | Resolution path |
|---|---|---|---|
| 1 | M. fortuitum `variants.csv`: 1 data row (erm(39)) vs header-only | B2 vs B4 | `wc -l` + header check at a pinned ntm-db commit |
| 2 | ntm-db counts off-by-one throughout (lines vs data rows) | B2 vs B4 | same |
| 3 | NTM-Profiler speciation: sylph (code) vs mash (README) | B2 vs B4/B6 | code wins; README is stale |
| 4 | ntm-db canonical repo: `pathogen-profiler/ntm-db` vs `jodyphelan/ntmdb` | B2/B4 vs B5 | jodyphelan/ntmdb README says "Repo moved"; **B5's correction is itself stale** |
| 5 | 10,105 tbdb rows: "fabricated WHO authority" vs WHO's own stated RIF→RPT / ETO→PTO inheritance | B4 vs B3 | split labels: rule-derived (RPT/PTO) vs inference (PMD) |
| 6 | WHO xlsx column indices 1-based vs 0-based | B3 vs B4 | pin the base, assert header strings |
| 7 | MTBseq flag `--minphred` vs `--minphred20` | B5/MANUAL vs B4/source | source wins (`'minphred20:i'`) |
| 8 | tbdb `barcode.bed`: 1,111 (shipped/master/who_v2) vs 1,114 (who_v2+ HEAD) | B4 vs B7 | pin by commit + md5 (`cb70a99a42dc1a2cbc8f289cc21d0b11`) |
| 9 | L4.2.1/L4.2.2 = TUR vs Ural | tbdb vs MTBseq | pick one vocabulary, print a crosswalk |
| 10 | Colpus lineage concordance 406 (405.5%→95.5%) vs 407 (95.8%) | within B6's source | preprint is internally inconsistent; cite both |
| 11 | Baker et al. species concordance 98.3%/0.3 SNPs (API abstract) vs 95%/1.0 (search snippet) | within B6 | resolve against the PDF |
| 12 | Hall 2023 katG homopolymer evidence: "no mycobacterial ONT homopolymer evidence exists" (prior) vs verified 4 discordances | B6 internal refutation | evidence exists; flag katG |
| 13 | "NTM stack has no SNP-distance clustering" | B4 design implication vs B2 verified `--snp_dist` | capability exists, is unused |