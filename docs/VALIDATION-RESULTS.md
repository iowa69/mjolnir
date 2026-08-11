# Validation — first results

Run 2026-08-11 on the machine Mjolnir was written on. Everything below was
measured; nothing is projected. Where a run establishes agreement rather than
correctness, it says so.

**Headline: the tool works end to end, and validation found one missing feature
and five defects that 607 unit tests did not.** Four of the defects are fixed;
the missing feature is characterised below and is not fixed.

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

It also means the outbreak-structure test in `VALIDATION.md` §2.2 **has not been
run**, and cannot be until an NTM repeat mask exists. §9 says NTM references get
their own mask computed at database build time; that is unimplemented. Until it
is, cohort mode works only for MTBC.

One further observation, recorded rather than explained: on sample `102-20`,
`samtools depth` was killed with exit `-9` and coverage was reported as not
measured. Memory was not exhausted when checked afterwards, so the cause is
unknown and should be reproduced before anything is concluded from it.

### The deliberate negative control

No Kraken2 index was configured. The contamination screen reported
`screen_informative: false` and the report said the screen was not performed —
it did not report the sample as clean. The refusal in §8 works in production.

---

## 2. The missing feature — read this before trusting a consensus

**No module in Mjolnir attaches gene names to called variants.** Measured:

| Sample | Variants called | With a gene name | WHO rows matched | MTBseq | tbdb |
|---|---|---|---|---|---|
| *M. bovis* `ERR11267966` | 3,013 | **0** | 110 | **0** | **0** |
| *M. chimaera* `OSR-30-20` | 5,050 | **0** | 0 | **0** | **0** |

WHO is matched on genomic coordinates, so it works. MTBseq and tbdb are matched
on the `<gene>_<hgvs>` key, so with no annotation they cannot match anything.

Three consequences, all live:

1. **The three-catalogue consensus is WHO-only in practice.** The
   `R (outside WHO catalogue)` result — the whole point of consulting the other
   two — cannot currently arise from real data.
2. **The NTM resistance rules never fire.** `erm(41)` sequevar typing, `rrl`
   2058/2059 and `rrs` 1408 are all gene-keyed. On the chimaera isolate they
   were reported as not established, which is honest but means NTM resistance is
   presently unimplemented in effect.
3. **A chimaera isolate matches no catalogue at all**, so its resistance section
   is empty rather than wrong.

The report no longer disguises this: catalogue columns that could not be
consulted render as `--` with a legend entry saying so, instead of `ND`, which
would have read as "that catalogue looked and found nothing".

**This is a missing feature, not a bug to patch quickly.** A correct annotator
has to get codon numbering, strand, promoter offsets and indel representation
right, and getting them wrong produces confidently wrong clinical calls. It is
the top item of remaining work.

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

- **Variant gene annotation** (§2 above). Everything else is downstream of it.
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
