# Validation — first results

Run 2026-08-11 on the machine Mjolnir was written on. Everything below was
measured; nothing is projected. Where a run establishes agreement rather than
correctness, it says so.

**Headline: the tool works end to end. Validating it found one missing feature
and thirteen defects that the unit suite did not — because every one of them
lived in the plumbing between tested functions, or needed real data to appear.**
All are now fixed, and each carries a regression test. Three adversarial review
passes over the code found six more, including the worst defect this project has
had; those are in §5.

---

## 1. What was validated against known truth

Three isolates whose identity is not in question, downloaded from ENA.

| Sample | Expected | Mjolnir | |
|---|---|---|---|
| `ERR3376643` *M. bovis* BCG | MTBC; animal lineage; intrinsic PZA resistance | MTBC; **La1.2**; **Pyrazinamide R** | ✓ |
| `ERR11267966` *M. bovis* | MTBC; La1.x; intrinsic PZA resistance | MTBC; **La1.6**; **Pyrazinamide R** | ✓ |
| `ERR10370893` H37Rv | MTBC; lineage 4.9; pan-susceptible | MTBC; **lineage4.9**; **no determinant detected** | ✓ |

The pyrazinamide call is `pncA_p.His57Asp`, graded `1) Assoc w R`, selected from
eleven candidate variants in graded genes. That is the textbook *M. bovis*
mutation and the consensus engine reduced to it correctly.

**What this establishes:** species-complex assignment, MTBC lineage barcoding
including the animal lineages, and WHO-catalogue resistance calling on a mutation
with an unambiguous expected answer. **What it does not establish:** accuracy
against phenotypic DST, which no sample on hand has.

### Species discrimination inside MAC

A *M. chimaera* isolate from the local outbreak collection
(`CHIMAERA-OSR/30-20`), against the ten-genome ANI reference set built for this
run:

| Reference | ANI |
|---|---|
| *M. intracellulare* subsp. *chimaera* | **99.402%** |
| *M. intracellulare* | 97.318% |
| *M. avium* subsp. *hominissuis* | 88.893% |
| *M. avium* subsp. *avium* | 88.522% |
| *M. kansasii* | 82.559% |

The ordering is right and the *chimaera*/*intracellulare* margin is two points.
The call is reported as **MAC, not resolved to species**, because §6 requires
marker SNPs inside MAC and no marker file is installed. That is the designed
answer, not a failure.

### The cohort now compares, and the denominator promise held on real data

*Superseded by the 2026-08-11 afternoon run; the refusal described below was
correct and its cause is now fixed.* Four *M. chimaera* isolates — two clinical
from one hospital, one from **water**, one from a **heater-cooler unit** — run at
a 6-SNP threshold against a computed repeat mask (277,412 bp, 4.52%):

| pair | SNPs | shared callable sites |
|---|---|---|
| 30-20 vs 31-20 (both clinical, same site) | **12** | 5,572,304 |
| 240-19 (water) vs 30-20 | 20,510 | 5,156,129 |
| 240-19 (water) vs 31-20 | 20,509 | 5,155,550 |
| 102-20 (heater-cooler) vs any | **not computed** | — |

No clusters at 6 SNPs: the two same-site clinical isolates sit at 12, outside the
threshold the original investigation used, and the water isolate is four orders
of magnitude away — a clean negative control.

**102-20 is the result that matters most.** Its callable regions were lost when
`samtools depth` was killed, so every pair involving it is reported `NA` with the
reason — *"the shared callable region of this pair is unknown, so a distance
would have no denominator"* — rather than as a distance of zero. A zero there
would have placed an unrelatable isolate at the centre of the outbreak. The
promise that a distance cannot be obtained without its denominator held on real
data, under a real failure, unprompted.

Every row also itemises what was excluded: variant positions inside the mask,
positions with no comparable allele, indels not counted as SNPs, differences
dropped by the 12 bp proximity rule, and positions outside the shared callable
region.

### The earlier cohort refusal, and why it was right

Four chimaera isolates — two clinical from one hospital, one from **water**, one
from a **heater-cooler unit** — were run as a cohort at a 6-SNP threshold, the
value the 2022 local investigation used. All four:

- were called **MAC**, and
- were mapped to the same *M. chimaera* reference by ANI proximity, so their
  coordinates are comparable, and
- produced **no pairwise distances at all**, because:

> the mask `tbdb mask.bed` names contig `Chromosome` but this cohort was called
> against `NZ_CP015278.1`, `NZ_CP015279.1`, `NZ_CP015280.1`; under those names
> the mask would exclude nothing and every distance would silently include the
> repetitive regions it exists to remove.

**This is the design working, not failing.** The alternative — applying an MTB
mask that matches no contig, excluding nothing, and printing distances that look
like masked distances — is exactly the silent wrongness §9 was written against.
No `cohort.json`, no distance matrix and no clusters were written.

At the time, this meant the outbreak-structure test could not run at all: §9 said
NTM references get their own mask computed at database build time and nothing
implemented it. `cohort/mask.py` now does, and the run recorded above is the
result.

One further observation, recorded then as unexplained: on sample `102-20`,
`samtools depth` was killed with exit `-9`. **The kill was ours.**
`iter_output`'s cleanup ran on normal completion of the read loop and killed a
process that had written all its output but not yet been reaped, turning a
successful run into a failure — and costing that sample its callable regions,
and with them every distance it took part in. Found by review pass 3 and fixed;
streaming tools now get 30 seconds to exit on their own.

### NTM resistance, on a real *M. abscessus* isolate

Until 2026-08-11 the `erm(41)` / `rrl` / `rrs` rules had never fired once: no NTM
reference carried gene models, so every determinant was keyed on a gene name
that did not exist. With a 27-species panel carrying gene models, `DRR245157`:

| | |
|---|---|
| Species | **Mycobacteroides abscessus**, resolved to species |
| Variants named | 24,300 of 24,487 |
| `erm(41)` position 28 | **C**, read from a pileup at 39× |
| `rrl` / `rrs` assessed | yes / yes |
| **Clarithromycin** | **S**, moderate — *"a statement about inducible resistance only"* |
| **Amikacin** | **no-call** — *"rrs examined and no determinant found"* |

Both calls are the right shape. C28 `erm(41)` is non-functional, so there is no
inducible macrolide resistance — but the caveat says the call is about inducible
resistance and not about acquired `rrl` mutations, which is the distinction that
decides whether a patient gets a macrolide. And amikacin is a **no-call**
although `rrs` was examined and clean, because a clean `rrs` is absence of a
known determinant and not evidence of susceptibility.

The sequevar is read from a pileup rather than from the variant list on purpose:
T28C is a polymorphism, so against a subsp. *abscessus* reference C28 is a
variant and T28 is invisible, and against a subsp. *massiliense* reference the
reverse. Only the base itself answers the question.

### The reference panel, and six taxids that were wrong

Species identification is ANI against a reference set, and that set existed only
as ten files assembled by hand. `mjolnir db panel` now fetches **27 mycobacterial
species, each with its gene models**, and the gene models are the part that
mattered: every NTM resistance rule is keyed on a gene name, so against an
unannotated reference `erm(41)`, `rrl` and `rrs` are written, tested and dead.

Each genome is checked against the organism NCBI reports for its accession
before it enters the panel, and that check immediately caught **six wrong
taxids** in the first draft of the species table:

| taxid | intended | what NCBI actually returns |
|---|---|---|
| 1698 | *M. chelonae* | ***Brevibacterium epidermidis*** |
| 1096 | *M. leprae* | ***Chlorobium phaeobacteroides*** (a green sulfur bacterium) |
| 1769 | *M. fortuitum* | *M. leprae* |
| 1770 | *M. haemophilum* | *M. avium* subsp. *paratuberculosis* |
| 1794 | *M. gordonae* | another organism |
| 1305 | *M. phlei* | another organism |

**A wrong taxid does not fail.** It returns a real, well-formed genome for a
different organism, which then sits in the panel under the name that was asked
for — and every isolate matching it is given that name. The `M. haemophilum`
entry had already shipped, and was caught only because a real *M. avium* isolate
matched it at 98.53% next to *M. avium* at 98.57%: those two species are nowhere
near that close, and the number was biologically impossible.

The taxids are corrected, but the durable fix is that they are no longer
trusted. The manifest carries the name NCBI reports, not the one that was typed.

### The deliberate negative control

No Kraken2 index was configured. The contamination screen reported
`screen_informative: false` and the report said the screen was not performed —
it did not report the sample as clean. The refusal in §8 works in production.

---

## 2. The gap that validation found, and how it was closed

When this document was first written, **no module attached gene names to called
variants**. WHO is matched on genomic coordinates and worked; MTBseq and tbdb are
keyed on `<gene>_<hgvs>` and matched nothing at all:

| Sample | Variants | With a gene name | WHO | MTBseq | tbdb |
|---|---|---|---|---|---|
| *M. bovis* `ERR11267966`, before | 3,013 | **0** | 110 | **0** | **0** |
| *M. bovis* `ERR11267966`, after | 3,013 | **2,976** | 120 | 1 | **99** |

So the "consensus across three catalogues" was WHO alone, and the NTM
`erm(41)` / `rrl` / `rrs` rules — all gene-keyed — could never fire.

`engines/annotate.py` closes it. Naming has to match the catalogues *exactly* —
`rpoB_p.Ser450Leu` and `rpoB_p.S450L` are the same mutation and different
dictionary keys — so it is written against a gold standard rather than a
specification: the WHO catalogue ships a `Genomic_coordinates` sheet mapping
144,964 coordinate triples to the name WHO itself uses. Agreement went
**17.95% → 79.54% → 90.55% → 95.19%** as each class of mismatch was read and
fixed, and the suite now replays a deterministic sample of that sheet and fails
below 93%.

What the gold standard taught, in order:

- The GFF calls rRNA genes `rRNA_gene`, not `gene`, so **`rrs` and `rrl` were
  invisible** — the two genes carrying every aminoglycoside and macrolide
  determinant in both MTBC and NTM.
- A multi-base substitution is not an indel, and maps to several graded variants.
- WHO pools loss of function under `<gene>_LoF` rather than naming each
  frameshift, so a precise coordinate name matches nothing for exactly the
  `pncA` and `katG` mutations that matter most.
- Start-codon changes are named in nucleotides, not as `p.Met1Ile`.
- The same indel can be written at several offsets, and equivalence is at the
  **protein** level: `dnaA` codons 111 and 112 are `ACT` and `ACC`, different DNA
  and both threonine, so WHO files the deletion as `Thr112del`.

**Still open:** MAC marker SNPs, so *M. chimaera* can be resolved below the
complex; and no NTM reference has a GFF, so NTM variants remain unnamed and the
NTM resistance rules still do not fire on real NTM data. That is now a
data-availability problem rather than a missing capability.

---

## 3. Defects found by running the tool, and fixed

Each of these passed the unit suite.

**The anchor catalogue could not be found.** `load_catalogues()` hard-coded both
the directory and the filename and looked in `who/`, while `db fetch` installs
to `who-catalogue-v2/`. With no `db_dir` it did not resolve a root at all and
named the path as `None`. Every unit test passes a path in explicitly. The
directory now comes from the registry the fetcher installs against.

**The species call was measuring the wrong thing.** `mash dist a b c query`
names only `a` as the reference and treats every later path, the query included,
as a query against it. Mjolnir passed the reference set as a list, so the
candidate table was the *reference genomes' distances to each other*. An
*M. chimaera* isolate and an *M. bovis* isolate both reported an identical
`99.1268% against H37Rv` — the distance from H37Rv to *M. bovis*, in which
neither sample took any part. **The chimaera isolate was called MTBC.** The
output was well-formed, sorted and plausible throughout; only the sample was
missing. The references are now sketched into one file and the query compared
against that, and the new tests assert the shape of the invocation rather than
the shape of the result.

**The PDF failed to render on a real isolate.** Enough graded variants matched
per drug to build a table row 4,076 points tall — taller than the page — and
reportlab raised `LayoutError`. The HTML and TSVs, which do not paginate, were
written and looked complete. Page one now names the determinants that drive the
call and counts the rest.

**Cohort mode crashed before joining a sample.** `callable_regions` takes a
`name` argument that collided with `_stage`'s own parameter. The parameters
before `func` are now positional-only.

**A fetched reference could not be used.** `db fetch` installs the reference and
`run` refused it for having no `.fai`, treating a one-second `samtools faidx` as
if it were the multi-minute aligner index. The `.fai` is now built on demand;
the aligner index still waits for `--build-index`.

---

## 4. What the reference data cost

`mjolnir db fetch` pulled 61.8 MB — WHO catalogue v2 (48,152 graded rows over 15
drugs and 65 genes, exactly the published count), tbdb (49,330 rows), MTBseq's
list (1,547 rows) and its three NTM references, and H37Rv.

The ANI reference set is **not** shipped and had to be assembled by hand: ten
mycobacterial genomes, 54 MB, listed in `<db>/ani/references.tsv`. Design §14
left its composition open; this run answers that question with a working set
covering MTBC, MAC, *M. abscessus*, *M. kansasii*, *M. marinum* and
*M. fortuitum*. It should become a `mjolnir db fetch` target.

---

## 5. Still not done

- **MAC marker SNPs**, so *chimaera* can be resolved below the complex.
- **Phenotypic DST truth** for a resistance-accuracy claim. The drives on hand
  carry no MTB at all and the ENA samples have no DST attached.
- **The full 159-isolate chimaera collection.** Four samples were run here; the
  outbreak-structure test in `VALIDATION.md` §2.2 needs the whole set. The two
  blockers that prevented it — annotation and the NTM mask — are now closed.
- **Why `samtools depth` was killed (exit -9) on two samples.** Memory was not
  exhausted when checked afterwards. It costs that sample its coverage metrics
  and its cohort comparability, so it needs reproducing rather than guessing.
- **ONT.** No nanopore mycobacterial data was run. Every ONT threshold in the
  tool is still untested.

---

## 5. What three adversarial review passes found

Each pass ran several independent reviewers over a defined area, and every
finding was then attacked by a separate verifier whose default was that the
finding was wrong. Only survivors are listed.

### Pass 1 — the new science code

**A deletion removing a gene's start codon was an "upstream variant".** The
worst defect this project has had. `_hgvs_indel` classified a deletion by where
it *started*, so one beginning upstream and running into the gene was named
`c.-2_3delGTATG` with effect `upstream_variant` — even when it removed the whole
start codon, even when it removed most of the CDS. `is_loss_of_function()` then
returned False, the WHO loss-of-function rule never fired, and **a complete
`pncA` knockout — definitive pyrazinamide resistance — produced no determinant
at all.** A real WHO row, the 950 bp deletion at 2288689 that WHO grades
`pncA_p.Met1?`, came out as a regulatory nucleotide change. Absence of the gene
product rendered as normality, in the gene that defines the drug.

**A frameshift insertion named a codon one too low at codon boundaries.** On the
plus strand the inserted bases land *after* the anchor, so naming the anchor's
own codon is wrong whenever the anchor sits on a boundary. That flipped the rpoB
RRDR rule in both directions — on rifampicin, the drug the rule exists for —
manufacturing a Group 2 determinant from a frameshift beginning outside the
region and dropping one beginning at its first codon.

**Gene lookup ignored the contig.** `genes_at()` filtered on coordinate alone
and never read the contig it had already parsed. H37Rv has one replicon so this
happened to work; the *M. chimaera* reference has three, and a variant at
position 300 of a plasmid was matched against whichever gene spanned position
300 of the chromosome.

**A multi-base substitution took its name from its first base.** `GCT>TGA` is a
stop codon; read one base at a time it was called a serine substitution — a
nonsense mutation presented as missense.

**Three separate "it did not happen, so it was fine" defects.** Drug rows printed
`ND` — no determinant detected — for a sample whose variant caller had died;
the annotation check reported pass whatever it produced, including nothing; and
the mask refusal, which exists because distances through an over-large mask
would not mean anything, was caught and downgraded to a stderr warning.

**The reference note never reached a report.** `result.caveats.extend(...)` ran
before the note explaining that a sample was mapped to the nearest genome by ANI
had been created — a note that changes what every coordinate in that report
means.
