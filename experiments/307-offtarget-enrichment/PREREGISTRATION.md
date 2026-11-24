# 307 — PREREGISTRATION: ONE aggregate candidate-repression statistic, the CONSTRUCT as the unit, against a selection-matched joint null

Written and hashed **before any per-construct candidate expression value was computed**. The hash of
this file is `work/PREREGISTRATION.sha256`; it is appended to `docs/CONFIRMATION-LEDGER.md` as
**Claim C8** before any outcome exists. Every number in `RESULTS.md` is produced by code written
after this file was hashed.

## 0. What has been read at the moment of hashing, and what has not

**Read (sequence and metadata only, all of it 301's or 305's saved products, no h5ad opened):**
301's `promoters_cage_tss_hg19.fa.gz`, `vcc2025_spacers.csv`, `horlbeck_supp3.pkl` (the whole
hCRISPRi-v2 library table), 298's `seed298.csv`, and the *metadata* JSONs 305 wrote beside its sums
(`sums305_<split>.json`: construct names, target labels, cell and UMI counts, lane composition).
The three sequence-side products already computed and frozen are `work/seq307.py`,
`work/pool307.py`, `work/probe307.py` and their outputs (`constructs307.csv`, `genes307.csv`,
`poolhits307.pkl`, `seq307_meta.json`, `pool307_meta.json`, `probe307.json`).

**Not read:** no element of any `X/data` array; **no value of `sums305_<split>.npz`'s `sums` array
has been used in any calculation** (its `keys`, `shape`, `dtype` and the `genes` vector were
inspected to confirm 18,080 columns and that all 16,936 corpus genes are measured — no count was
read). No candidate effect, no control expression value, no basal decile has been computed.

**This is an extension on already-used data, not a reserve and not a blinded replication.** The three
2025 h5ads were read by 293 and 305; 305's integer sufficient sums are reused here rather than
re-streamed, so **no new access of any kind is declared or consumed**. The candidate list is 301's,
frozen from sequence alone.

## 1. What a positive would and would not establish

A positive result would say: *genes selected by a 12-nt-seed + NGG promoter match are, in aggregate
and with the construct as the unit, more repressed in the cells carrying that construct than
sequence-, promoter- and expression-matched genes that the same rule did not select.* It would be
**observational corroboration of a sequence-linked association**. It would **not** establish
dCas9-KRAB off-target binding: expression here is observational, the 34 constructs differ in
on-target knockdown efficacy (`kd_lfc` −3.04 … +0.03), in cell number (95 … 3,415), in population
composition and in handling, and a candidate gene's repression is equally compatible with
intended-target downstream biology, with a construct-specific population effect, and with chance
that the null happens not to model. A negative result means "no aggregate expression corroboration
under this test", never "off-targets do not occur" or "off-targeting does not explain failing
constructs". Nothing measured here can distinguish a binding event from a correlate of one.

## 2. The experimental unit and the exact population

**The unit is the targeting construct.** Programme rule 1: every interval and every p below is
labelled with it.

**Primary population — the 34 gene-targeting constructs of the 2025 library** (13 failing + 21
working by 293's `|kd_lfc| < 1.0` rule, frozen in 298's `seed298.csv`), each measured only in its own
split (train 21 / validation 7 / test 6 by the frozen mapping). These 34 constructs cover 17 paired
intended targets.

**Verified before hashing (KA-2, already run, PASS):** re-running 301's rule from its saved hits
reproduces its published counts exactly — failing **4/13**, working **4/21**, non-targeting **4/31**
constructs with at least one L=12 candidate. The 8 targeting constructs carrying an L=12 candidate
are `CAST_P1`, `GRB10_P2`, `KDM2B_ENST00000377071.4`, `KLF10_P1`, `DPF2_P1`, `IGF1R_P1`, `MED15_P1`,
`TMSB4X_P2` — **the coordinator's count of 8 is confirmed**, and their 11 construct–gene pairs are
exactly 305's primary family.

**Separately declared population — the 31 non-targeting constructs.** Reported by the identical
procedure as a declared secondary with no verdict. Following 305's coordinator review they are
treated as *candidate-bearing observations selected by the same rule*, **not** as negative controls.

## 3. The candidate set of every construct (301's rule, unchanged)

301's frozen selection rule, reproduced verbatim from its saved promoter corpus: exact match of the
PAM-proximal L nt of a protospacer (either of the construct's two spacers), PAM-anchored, canonical
NGG immediately 3′, inside a CAGE-TSS promoter window, offset from TSS in the **CRISPRi active
window −50 … +300**, own-target symbol excluded, and 301's **10 kb distal filter** (a candidate whose
promoter window lies within 10 kb of any window of the construct's own target is the construct's own
locus under a neighbour's symbol, not an off-target).

**Candidate set `C_c` := the hit set at `L*(c)`, the longest L ∈ {12, 11, 10, 9} yielding ≥ 1 hit.**
This is 301's own candidate definition (`seed_scan.py`'s `cand`; `seed298_plan.md` step 5) applied
after the active-window and distal filters, and it is what makes the statistic computable for
**every** construct: measured before hashing, all 34 targeting constructs have a non-empty `C_c`
(L\* = 12 for 8, 11 for 16, 10 for 7, 9 for 3; k = |C_c| runs 1 … 18, median 1). No construct is
dropped for having no strong candidate; a construct whose best evidence is a 9-mer contributes an
observation that will sit at its own null unless it does not. One of the 31 NTC constructs has no
hit at any L; **declared rule: a construct with k = 0 contributes z = 0**, i.e. exactly the null
expectation, to both the observed statistic and to every null draw.

The L\* rule is frozen here and may not be changed after any expression value is seen. Fixed-L
readings (L = 12 only, etc.) are declared as secondaries in §8, not as alternative primaries.

## 4. The response formula (305's, unchanged)

`e_c(g) = prof(S_c)[g] − prof(S_ctrl(split(c)))[g]` with `prof(S)[g] = log1p(50000 · S[g] / Σ_j S[j])`,
the denominator always the full 18,080-gene library total. `S_c` = the exact integer count sums over
**both** halves of construct `c`'s cells; `S_ctrl` = the exact integer sums over **all 31
non-targeting constructs** of that split (305's pooled A+B reference). For a non-targeting construct
under test (the §2 secondary population only), the reference is **leave-that-construct-out**, 305's
rule. Units: natural-log CP50k. All of this is read from 305's saved
`work/sums305_<split>.npz`; **the expression data are not re-streamed** — the sums contain every
quantity this design needs, because every statistic below is a function of per-construct and
per-control integer column sums.

`basal(g) = prof(S_ctrl(split))[g]`, the control-pool expression used as a matching covariate.

## 5. THE JOINT NULL

The null resamples **candidate sets, not genes independently**. For construct `c` at level
`L = L*(c)`, each real candidate `G ∈ C_c` is replaced by a gene drawn from a pool matched to `G`.

### 5.1 What each gene's eligibility and matching covariates are

Measured before hashing, for all 16,936 genes with a CAGE promoter window (all 16,936 are among the
18,080 measured genes):

| covariate | definition | axis it serves |
|---|---|---|
| `pi_off_L(g)` | **the exact marginal probability that a protospacer drawn uniformly from the hCRISPRi-v2 library (209,070 spacers — the library the 2025 guides come from; all 68 real spacers found verbatim in it) hits `g` at seed length `L`, OFF TARGET, under 301's identical rule** (own-target symbol excluded, 10 kb distal filter applied) | (a) the sequence rule: seed length, PAM class, and the active-window offset constraint, weighted by the library's own empirical seed composition |
| `n_sites`, `mean_offset` | count and mean TSS-offset of the gene's active-window PAM-anchored 12-mer anchors | (a) where in the promoter window a match can fall |
| `gc` | GC fraction of the active window | (b) promoter GC |
| `n_windows`, `tss_density_10kb` | number of CAGE promoter windows (P1/P2/P1P2), and number of **other** genes with a TSS within 10 kb | (b) promoter architecture / TSS density |
| `map_mult` | mappability proxy: mean number of **distinct other genes** sharing each of the gene's active-window 12-mers inside the promoter corpus (1.0 = every 12-mer unique) | (c) mappability |
| `lowcomp` | fraction of the gene's active-window 12-mers that are ≤ 2 distinct bases or contain a ≥ 6-nt homopolymer | (c) repeat/low-complexity proxy |
| `basal(g)` | control-pool `log1p`-CP50k in the construct's own split | (d) basal expression |

`pi_off` is the *measured* propensity, not an approximation of it: it is computed by rescanning all
209,070 hCRISPRi-v2 protospacers through 301's rule. Measured calibration of that pool against the
real library, before hashing: per-spacer P(≥ 1 off-target hit) at L = 12 is **0.132 for the 68 real
2025 targeting spacers vs 0.159 for the pool** (0.412/0.389 at L=11, 0.750/0.694 at L=10,
0.956/0.920 at L=9); all 68 real spacers carry Horlbeck "off-target stringency" **0**, so the
secondary spacer-resampling null of §8 restricts the pool to stringency 0 (202,064 spacers). A
pseudo-construct built from a random pool spacer **pair** has L\* = 12/11/10/9 with frequency
0.297/0.330/0.281/0.087 versus the real 34 constructs' 8/16/7/3 = 0.235/0.471/0.206/0.088.

### 5.2 The eligible universe (the candidate-specific exclusion 305's null omitted)

A gene `g` is eligible as a pseudo-candidate for construct `c` iff **all** of:

1. `pi_off_{L*(c)}(g) > 0` — the rule could have selected it at that seed length;
2. **`g` is not a gene that ANY construct in the 2025 library intends to target** (300 such genes,
   read from the three splits' `target_of_construct` metadata) — this is the exclusion 305's all-gene
   null omitted;
3. `g` is not `c`'s own intended target;
4. `g` is not in `c`'s hit set at **any** L ∈ {12, 11, 10, 9} — a gene the rule did touch is not a
   clean non-candidate;
5. `g` is measured in `var/_index` (all 16,936 are).

### 5.3 Matching and drawing

For candidate `G` of construct `c` at level `L`:

- **hard strata (declared):** the eligible universe is cut into 10 × 5 = **50 cells** by
  `basal` decile (quantiles over the eligible universe of `c`'s split) × `pi_off_L` quintile
  (quantiles over eligible genes with `pi_off_L > 0`). `G`'s pool starts as its own cell.
- **soft matching inside the cell:** the pool is the **M = 200** genes of that cell nearest to `G` by
  Euclidean distance in the z-scored space `[gc, log1p(tss_density_10kb), map_mult, lowcomp,
  mean_offset/100, log1p(n_windows), log1p(n_sites)]` (z-scored over the eligible universe). If the
  cell holds fewer than 50 eligible genes, the whole cell is used; if fewer than **20**, the pool is
  widened to the same `pi_off` quintile ∩ the adjacent `basal` deciles and the widening is recorded
  per candidate in `work/enrichment307.csv`.
- **a null draw** replaces every `G ∈ C_c` by one gene drawn uniformly from `G`'s own pool, without
  replacement across the k slots of `c`.

Matching quality is reported as standardized mean differences (KA-6).

## 6. THE ONE PRIMARY STATISTIC

For candidate `G` of construct `c` with pool `P_G`:

```
q(G)  = ( #{g in P_G : e_c(g) < e_c(G)} + 0.5 · #{g in P_G : e_c(g) = e_c(G)} ) / |P_G|
r_c   = mean over G in C_c of q(G)                      (k_c = |C_c|; r_c = 0.5 under the null)
z_c   = sqrt(12 · k_c) · (0.5 − r_c)                    (POSITIVE = candidates more repressed)
```

**PRIMARY STATISTIC**

```
S = (1 / 34) · Σ over the 34 targeting constructs of z_c
```

`S` is a **per-targeting-construct standardized candidate-repression enrichment**. It is a
within-construct rank statistic: a shift common to all of a construct's genes (a shared control
pool, a library-size effect, a lane effect on that construct) cancels, so `S` is insensitive to the
between-split control sharing 305 disclosed. It is computable for every construct, and a construct
with `k_c = 0` contributes `z_c = 0` identically to the observed value and to every null draw.

**Primary p-value — the joint null.** `B = 100,000` joint draws; in draw `b` every construct's whole
candidate set is redrawn as in §5.3 and `S^(b)` recomputed by the identical code path.

```
p_joint = (1 + #{ b : S^(b) >= S_obs }) / (B + 1)          one-sided, repression direction
```

**Interval.** A **construct-level bootstrap**: 20,000 resamples of the 34 targeting constructs with
replacement, percentile 95% CI on `S`. It is labelled *n = 34 targeting constructs* and it measures
construct-to-construct variability, **not** the null; the verdict rests on `p_joint`.

**Interpretable effect size (reported with the primary, no separate verdict).**
`Δ = (1/34) Σ_c [ mean_{G∈C_c} e_c(G) − mean_{G∈C_c} median(e_c over P_G) ]` in log1p-CP50k units,
with its own construct-level bootstrap CI.

## 7. PASS RULE, VERDICT CELLS AND ATTAINABLE RESOLUTION — fixed before computing anything

- **Number of verdict cells: ONE.** No multiplicity correction is required; every other quantity in
  §8 is reported without a verdict and cannot change this one.
- **α = 0.05, one-sided** (candidates *more* repressed than matched non-candidates).
- **C8 PASSES iff both hold:**
  1. `p_joint ≤ 0.05`; **and**
  2. **leave-one-construct-out robustness:** recomputing `S` and `p_joint` with each of the 34
     constructs deleted in turn (denominator 33), the **maximum** of the 34 leave-one-out p-values is
     `≤ 0.05`. No single construct may carry the result. (This condition only makes the rule
     stricter and is stated before any outcome; it may not be dropped afterwards.)
- **C8 FAILS** if condition 1 or condition 2 does not hold.
- **Attainable-resolution floor.** The joint null's smallest attainable p is `1/(B+1) = 1.0e-5`,
  three orders of magnitude below α. **Unlike 305's same-split rank test (floor 11/158 ≈ 0.0696
  after Holm), this design CAN attain significance at its stated α**, so `FAIL` here is informative
  and `UNRESOLVABLE-AT-THIS-RESOLUTION` is reserved for the single case in which a known-answer row
  aborts the analysis. The *power* of the design is a separate matter and is not knowable before the
  outcome; the reported CI bounds what the data can say either way, and a FAIL will be reported as a
  bounded estimate with that interval, not as evidence of absence.
- **Sign convention frozen:** `S > 0` = candidates repressed relative to matched non-candidates.

## 8. DECLARED SECONDARIES — reported, no verdict, cannot change §7

1. **The failing-vs-working construct contrast (301's original T7 question):** `z_c` for the 13
   failing vs the 21 working constructs — difference in means with a construct-level bootstrap CI and
   a Mann–Whitney p. Unit: the targeting construct.
2. **The 31 non-targeting constructs** through the identical pipeline (`S_NTC`, its joint-null p,
   per-construct table). Labelled as candidate-bearing observations selected by the same rule, not
   negative controls; the two nominal 305 hits (`DCTN2`, `PRKCQ`) are reported inside it.
3. **`GRB10_P2 → SNX8`** with its partner construct `GRB10_P1`: both effects, both matched
   percentiles, and `GRB10_P1`'s own candidate (`USPL1`, L\* = 11).
4. **Per-construct cell counts, control UMIs and the measurement floor**: how many candidates sit
   below 305's `basal ≥ 0.5` CP50k floor, per construct; and `S` recomputed on candidates above the
   floor only.
5. **A second, independent null — real-spacer resampling.** Instead of matched gene draws, a
   pseudo-construct is a random **pair of stringency-0 hCRISPRi-v2 protospacers** rescanned through
   301's exact rule; draws are kept by rejection only if they yield `L* = L*(c)` and exactly `k_c`
   candidates. This preserves the sequence rule *jointly* (co-occurrence within a set) but does not
   match expression; it is the complement of the primary null, not a replacement.
6. **Level and population restrictions:** `S` on the 8 L\*=12 constructs alone; `S` computed at
   fixed L = 12, 11, 10, 9.
7. **Is the 2025 library unusually promiscuous?** The 34 constructs' L\* distribution against the
   pool's pseudo-construct distribution (§5.1), χ².
8. **KA-4's cross-construct matrix** as a descriptive scale for `z_c`.

## 9. KNOWN-ANSWER ROWS — run before any candidate outcome; a failure aborts the conclusions

| id | row | pass rule | status |
|---|---|---|---|
| KA-1 | 305's 15 published pooled effects and both half effects recomputed from `sums305_<split>.npz` by code written here | max abs difference ≤ 1e-9 vs `results305.json` | to run |
| KA-2 | 301's selection reproduced from its saved hits with 301's exact filters | failing 4/13, working 4/21, NTC 4/31; 8 L=12 targeting constructs; the 11 pairs are 305's | **RUN, PASS** (`seq307.log`) |
| KA-3 | 305's finding that the 31 NTC constructs' sums are identical across the three splits | exact equality of all 31 constructs' per-half integer sums across train/validation/test | to run |
| KA-4 | **no-signal row**: `S` recomputed with each construct scored against **another** construct's candidate set (200 random derangements of the 34 constructs, pools rebuilt for the scoring construct) | \|mean S_perm\| ≤ 0.05 and the 2.5–97.5 % interval brackets 0 | to run |
| KA-5 | **null calibration**: 200 null draws each treated as "observed" and given a p against the other draws | KS test of uniformity p ≥ 0.05, and the fraction ≤ 0.05 inside [0.02, 0.10] | to run |
| KA-6 | **matching balance**: standardized mean difference between real candidates and their matched pools on all seven covariates of §5.1 plus `basal` | \|SMD\| ≤ 0.25 on every axis | to run |
| KA-7 | 293's on-target `kd_lfc` for all 34 constructs recomputed from 305's sums | 34/34 within 1e-9 of `seed298.csv` | to run |

KA-2 is already run and recorded above because it is pure sequence work and contains no expression
value. KA-1, KA-3, KA-6 and KA-7 touch expression but only reproduce numbers **already published**
by 293/305 or describe covariates, never a new candidate outcome; they run before the primary.

## 10. Seeds, determinism, resources, and what is not done

`SEED307 = 20260907`; `numpy.random.default_rng`. `B = 100,000` joint null draws, 20,000 construct
bootstrap resamples, 200 derangements, 200 calibration draws — all fixed here. Everything is read
from 305's saved `.npz` sums and 301's saved sequence products; **no h5ad is opened, nothing is
re-streamed, no download, no reserve access, no new seed search.** Peak RSS budget < 6 GiB (the
largest object is a 34 × 18,080 float64 effect matrix, ≈ 5 MB). Nothing is written outside
`experiments/307-offtarget-enrichment/` except the C8 entry in `docs/CONFIRMATION-LEDGER.md`.
No `vcc`, no `cell-eval2`, no build, no submission, no model, no production change.

## 11. Deviations

Every deviation from this document will be listed in `RESULTS.md` §Deviations, each marked as fixed
**before** or **after** any outcome was seen.
