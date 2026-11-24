# 305 — PREREGISTRATION: expression evidence for 301's frozen T7 seed candidates

Written and hashed **before any candidate-gene expression value was read**. Everything read up to this
point is metadata: `obs/{guide_id,target_gene,batch}`, `var/_index`, `X.attrs['shape']`, file sizes and
sha256s, plus 301's saved seed hits and 298's `seed298.csv`. **No element of any `X/data` array has been
read.** The hash of this file is `PREREGISTRATION.sha256`; every number in `RESULTS.md` is produced by
code written after it.

## 0. What this is, and what it is not

**This is an EXTENSION ON ALREADY-USED DATA, not an untouched confirmation.** The three 2025 VCC h5ads
were read by 293 (the construct-ceiling measurement) and their per-construct knockdown scalars are
already published in `experiments/293-challenge-construct-ceiling/work/run293_*.json` and
`construct293.csv`. The failing/working labels tested here were **derived from those same files** by 293
(`|kd_lfc| < 1.0`, 13 of 34 arms) and frozen by 298. The candidate list was frozen by 301 from sequence
alone. Nothing here is a reserve; nothing here is a blinded replication. The only thing that is genuinely
new is the **candidate genes' own columns**, which no run has looked at.

**It is observational corroboration, not causal proof.** A repressed candidate gene is consistent with a
seed-driven off-target, and also with intended-target downstream biology, with a construct-specific
population/handling effect, and with chance. Nothing measured here can distinguish a dCas9-KRAB binding
event from a correlate of one. A negative result means "no expression corroboration under this test",
never "off-targets do not explain failing constructs".

**The seed scan itself is inconclusive and stays that way.** 301's primary reading is 4/13 failing vs
4/21 working constructs with a candidate (background 4/31 non-targeting), Fisher p = 0.679. That is not
evidence of equivalence and not evidence of absence. This run does not re-open, re-tune or re-run the
scan: seed length (12 nt), PAM rule (NGG, PAM-anchored), active window (TSS −50 … +300), distal filter
(≥ 10 kb from any window of the arm's own gene), candidate gene set, the failing definition
(`|kd_lfc| < 1.0` on 293's full-depth on-target lfc) and the effect threshold below are **frozen** and may
not be changed after any expression value is seen.

## 1. The frozen candidate table

Rebuilt by `work/candidates305.py` from 301's `seed_hits301_cage.csv` with 301's exact filters;
duplicate spacers/hits collapse to one construct–gene pair. The rebuild **reproduces 301's published
counts exactly: failing 4/13, working 4/21, non-targeting 4/31**, 24 self-locus hits dropped over 3
construct–gene pairs. `work/candidates305.csv` (sha256 in `work/input_manifest.json`) is the frozen
table. **All 15 evaluated pairs are frozen, not only the failing arm's:**

| # | group | construct | intended target | split | candidate gene | offsets from TSS |
|---:|---|---|---|---|---|---|
| 1 | failing | `CAST_P1` | CAST | train | UTP18 | +259 |
| 2 | failing | `GRB10_P2` | GRB10 | test | NEUROG3 | +274 |
| 3 | failing | `GRB10_P2` | GRB10 | test | SNX8 | +291, +299 |
| 4 | failing | `KDM2B_ENST00000377071.4` | KDM2B | train | WRN | −24 |
| 5 | failing | `KLF10_P1` | KLF10 | train | ST3GAL3 | −27 |
| 6 | working | `DPF2_P1` | DPF2 | test | CEP120 | +247 |
| 7 | working | `IGF1R_P1` | IGF1R | validation | IFITM5 | +114 |
| 8 | working | `MED15_P1` | MED15 | test | AFF4 | +106 |
| 9 | working | `MED15_P1` | MED15 | test | TMCC1 | +92 |
| 10 | working | `MED15_P1` | MED15 | test | ZSWIM5 | +208 |
| 11 | working | `TMSB4X_P2` | TMSB4X | train | GNG10 | +141 |
| 12 | NTC | `non-targeting_00010\|non-targeting_02551` | — | all | DDB1 | +278 |
| 13 | NTC | `non-targeting_00018\|non-targeting_00127` | — | all | DCTN2 | +177 |
| 14 | NTC | `non-targeting_00095\|non-targeting_03143` | — | all | CLNS1A | +295 |
| 15 | NTC | `non-targeting_00121\|non-targeting_00339` | — | all | PRKCQ | +132 |

**Primary family = pairs 1–11** (11 pairs, 9 constructs, 5 failing + 6 working). **Negative-comparator
family = pairs 12–15** (the non-targeting guides' candidates), evaluated by the identical procedure and
Holm-corrected inside their own family; they are comparators, not primary evidence.

All 15 candidate genes are present in `var/_index` (18,080 genes, identical gene list and digest
`da2b3839…` in all three splits). **No construct anywhere in the 2025 library targets any of the 15
candidate genes**, so the own-targeting exclusion of §6 removes nothing; the rule is kept in force
regardless.

## 2. Sources, fixed by identity

The three 2025 files (`data/external/vcc2025/{train,validation,test}/adata_*.h5ad`), read-only, with
bytes and sha256 recorded in `work/input_manifest.json`; 221,273 / 98,927 / 170,846 cells × 18,080 genes;
189 / 85 / 136 constructs (158 / 54 / 105 gene-targeting + 31 non-targeting in each); 48 lanes each.
These are the files 293 used; they are **not** a new reserve and this run declares no reserve.

**Measured fact recorded before the test:** the 31 non-targeting constructs, their 38,176 cells and
their per-lane composition are **byte-identical across the three splits** — the control pool is the same
cells in every file. Controls are therefore *shared between splits* (never between the two halves of one
comparison; see §5), and this is a limitation of the design, stated up front.

## 3. Split / construct mapping

Each arm is measured **only in its own split** (verified: 34/34 arms occur in exactly one split, and
their full cell counts match `seed298.csv` to the cell). `construct_of()` is 293's, verbatim:
`GENE_P1_A|GENE_P1_B → GENE_P1`; a non-targeting pair keeps its full id. The mapping is frozen in
`work/mapping305.json`.

## 4. The response formula (frozen)

`lib305.effect()`, which is 285's `lib285.prof`/`delta` on one row pair — the same arithmetic 293 used
for its on-target `kd_lfc`:

```
prof(S)[g] = log1p( 50000 * S[g] / sum_j S[j] )          # sum over ALL 18,080 genes
e(c, g)    = prof(S_c)[g] - prof(S_ctrl)[g]              # the DIRECT, control-referenced effect
```

`S_c` = exact integer count sums over the construct's cells; `S_ctrl` = exact integer count sums over the
referenced non-targeting cells. **The normalising denominator is always the full library total**, whatever
subset of columns is afterwards inspected. Units: natural-log CP50k; e = −0.10 is a ≈ 9.5 % reduction for
a well-expressed gene.

**Residual (reported separately, never as the repression criterion).** With 159's panel centring,
`r(c, g) = e(c, g) − mean_{c' in the split's targeting constructs} e(c', g)`. It is declared in advance
that within one split the rank of `r` over constructs is *identical* to the rank of `e` (the panel mean
is a common shift), so the residual is reported for its **sign and magnitude**, to keep "below the panel
mean" and "below controls" visibly separate; the empirical rank is reported once and labelled as applying
to both.

## 5. Halves and controls

**Construct halves.** Deterministic, no RNG: within each batch, the construct's cells are taken in
ascending file-row order and assigned alternately to half 0 and half 1 (`lib305.halves`). This is a pure
function of the file's row order and is batch-balanced by construction. Cell and UMI counts and the
per-lane composition of each half are reported.

**Control halves.** 293's partition, verbatim (`run293.py:92`): the split's 31 non-targeting constructs
are sorted, and `ntc_par[c] = i % 2` assigns them to control half A (16 constructs, 21,581 cells) and
control half B (15 constructs, 16,595 cells). The two control halves share **no cell and no guide
construct**. Half 0 of a construct is referenced to control half A, half 1 to control half B; the pooled
effect is referenced to A + B.

**Non-targeting constructs under test** (pairs 12–15) are themselves inside one control half. For those,
and only those, the containing half is taken **leave-that-construct-out** (its own exact sums subtracted),
so a tested NTC construct never appears in its own reference. This asymmetry is declared here and
reported.

**No cell bootstrap will be labelled biological replication.** The 2025 files carry no biological
replicate of a construct: the two halves are two halves of one construct's cells in one experiment, and
every interval computed from them is a cell-sampling interval. No biological replicate is invented.

## 6. Statistics

For each pair (construct `c`, candidate gene `G`, split `s`):

- **pooled direct effect** `e_pool = e(c, G)` on all of `c`'s cells vs control A+B;
- **half effects** `e_0 = e(c half 0, G)` vs control A, `e_1 = e(c half 1, G)` vs control B;
- **empirical null**: the pooled direct effect of gene `G` under every **other gene-targeting construct
  of the same split**, excluding (i) `c` itself and (ii) any construct whose `target_gene` is `G`.
  Non-targeting constructs are not in the null. Minimum null size **N_null ≥ 30**; a pair whose null is
  smaller is **NOT TESTABLE** and no other split is substituted for it;
- **empirical rank p-value**, one-sided (more negative is more extreme), with the +1 correction:
  `p_raw = (1 + #{null <= e_pool}) / (1 + N_null)`;
- **Holm** over all 11 primary pairs (and, separately, over the 4 comparator pairs);
- **robust outlier size** `z_rob = (e_pool − median(null)) / (1.4826 · MAD(null))`, reported always;
- **secondary null (declared here, reported as secondary)**: the same statistic against the 317
  gene-targeting constructs of all three splits, each referenced to its own split's control pool. This
  exists only because the within-split null's resolution floor makes the primary Holm criterion
  unattainable (§8); the primary reading remains the within-split one.

## 7. Pass rules (frozen)

A primary pair is **CORROBORATED** iff all four hold:

1. **C1 direction** — `e_pool < 0`;
2. **C2 magnitude** — `e_pool <= −0.10` (≈ 9.5 % repression; one tenth of the 1.0-lfc failing threshold);
3. **C3 both halves** — `e_0 < 0` **and** `e_1 < 0`, each against its own disjoint control half;
4. **C4 unusualness** — `p_raw <= 0.01` in the within-split empirical null of §6.

The Holm-adjusted p is reported for every pair alongside C4 and is **not** part of the conjunction, for
the reason given in §8. Effect magnitude is reported for every pair whether or not any criterion fires.

**Testability.** A pair is NOT TESTABLE if any of: the candidate gene is absent from `var/_index`
(none are); the pooled control expression of `G` is below **CP50k ≥ 0.5** (≈ 10 CPM — the standard
expressed-gene floor; a gene below it has too few control counts for the formula to mean anything); the
construct has fewer than **30 cells in either half** (all 34 arms have ≥ 95 cells, so none are expected
to fail); or `N_null < 30`. Untestable pairs are reported as such and are excluded from the Holm family
(with the family size recorded both ways).

## 8. Declared power limit (before any expression value)

The empirical rank p-value has a hard floor of `1/(N_null + 1)`: **0.0063** in train (158 targeting
constructs), **0.0185** in validation (54), **0.0095** in test (105). With 11 primary pairs, the smallest
attainable Holm-adjusted p is `11 × 0.0063 = 0.070 > 0.05`. **Holm significance at family-wise α = 0.05 is
therefore unattainable in the primary null by construction, for any effect however large**, and a report
of "0 pairs Holm-significant" would carry no information. This is why C4 is stated as a raw-rank
criterion, why the Holm p is reported rather than thresholded, and why the all-splits null (N ≈ 315,
floor 0.0032, smallest attainable Holm 0.035) is declared in advance as the secondary reading. All of
this is fixed before any expression value is read.

The design is also underpowered at the level it is being asked about: 11 pairs over 9 constructs, and
uncertainty clusters at the **intended target (17 paired targets)**, not at 34 independent constructs.
No population off-target rate will be estimated from it.

## 9. Specificity and negative comparators

- **Paired-construct comparator (item 7)**: the same candidate gene `G` under the **other construct of the
  same intended target** (e.g. UTP18 under `CAST_P2` for the `CAST_P1` pair), with its own effect, halves
  and rank. A candidate repressed under *both* constructs of a target is more consistent with
  intended-target downstream biology than with a construct-specific event; a candidate repressed under
  only the seed-carrying construct is *consistent with* — never proof of — a construct-specific
  association. Both readings are reported; neither establishes a mechanism.
- **Non-targeting comparators**: pairs 12–15.
- **Control-vs-control comparator**: `e(control half A, G) − ` referenced to control half B, i.e. the
  effect of one control half against the other, for all 15 candidate genes. Its distribution is the
  scale on which "unusually negative" must be read.

## 10. Known-answer rows — run BEFORE any candidate conclusion; a failure aborts the conclusions

| id | row | pass rule |
|---|---|---|
| KA-1 | 34 construct identities and the 13/21 failing/working labels, read back from `seed298.csv` and cross-checked against `mapping305.json` cell counts | 34/34 identities, 13/21 labels, 34/34 cell counts equal |
| KA-2 | gene order and missingness: `var/_index` identical (same order, same sha256) across the three splits; every candidate gene present | identical digest; 15/15 present |
| KA-3 | 301's contrast reproduced from its saved hits with its exact filters | failing 4/13, working 4/21, NTC 4/31 |
| KA-4 | synthetic known fold change: a construct's stored sums are multiplied by a known factor `f` on one gene column and the formula re-run; the recovered effect must equal the analytic value | \|recovered − analytic\| ≤ 1e-9 |
| KA-5 | chunked streamed sums vs a direct small non-reserve slice: for one construct with < 1,500 cells, its rows are re-read directly by index (no chunk loop) and summed | exact equality of all 18,080 columns |
| KA-6 | same-NTC-halves null: control half A vs control half B on all 18,080 genes; the effect distribution must be centred on ~0 | \|median\| ≤ 0.01 and the 15 candidate genes' control-vs-control effects reported |
| KA-7 | selected own-target effects through 293's **original** cell selection and regime: the on-target `kd_lfc` of all 34 arms recomputed from this run's full-depth sums must reproduce `run293_<split>.json`'s `kd_full[c]['lfc']` | 34/34 within **1e-9 absolute** (dtype-aware: float64 accumulation in a different group order; both are full-depth, un-thinned, pooled-NTC-referenced, so the configurations are exactly comparable) |

**KA-7's comparable configuration is stated explicitly**: 293's `kd_full` uses *all* of a construct's
cells (not its 1,600-cell extraction cap) at *native depth* (`thin = 1.0`) against `prof(cA_full +
cB_full)`, i.e. the pooled NTC set — identical to this run's pooled reference. 293's *thinned*
(`thin = 0.368`) and *rung* (`kd_n`, 25–200 cells/side) numbers are **different configurations and are
not used as known answers**; no full-depth re-extraction here is required to equal a differently thinned
historical number.

## 11. Seeds, determinism, resources

Half assignment: deterministic, no RNG. Bootstrap/CI: 285's `lib285` with `SEED = 20260906`,
`NBOOT = 20000`. Synthetic-KA seed `SEED305 = 20260907`. Streaming: `CHUNK = 2048` rows,
one split at a time, **serially**, with the repo memory gate checked before each stream. Budget: peak RSS
< 6 GiB, 2 h wall clock. Nothing is written outside `experiments/305-t7-expression-test/`. No `vcc`, no
`cell-eval2`, no build, no submission, no download, no new seed search, no reserve access.

## 12. Deviations

Every deviation from this document will be listed in `RESULTS.md` §Deviations, each marked as fixed
**before** or **after** any outcome was seen.
