# 307 — RESULTS: the aggregate candidate-repression test, construct as the unit

Preregistration `PREREGISTRATION.md`, sha256
`81f7cda4533bc522995d2168e05542c673edacf161f5368080f0223ac801d9e7`, hashed **2026-09-07 16:00:01Z**
and appended to `docs/CONFIRMATION-LEDGER.md` as Claim C8 before any per-construct candidate
expression value was computed. No reserve declared or spent; no h5ad opened; 305's saved integer
sufficient sums reused throughout. Peak RSS 1.17 GiB.

---

## Abstract

**What was measured.** One aggregate statistic, with the **targeting construct** as the experimental
unit, over **all 34** gene-targeting constructs of the 2025 VCC library: are the genes that 301's
seed-match rule selected for a construct more repressed in that construct's cells than
sequence-, promoter- and expression-matched genes the same rule did not select? The joint null
resamples **candidate sets**, matched on the measured probability that a real hCRISPRi-v2
protospacer hits each gene off-target under 301's exact rule, on promoter GC, TSS count and density,
on a mappability proxy, and on basal control expression — with every gene the 2025 library targets,
the construct's own target and every gene its own spacers touch excluded.

**The number.** `S = +0.786` standardized units per targeting construct, construct-level bootstrap
95% CI **[+0.387, +1.166]** (n = 34 constructs); joint-null p = **1.0e-5**, the null's floor;
`Δ = −0.0790 log1p-CP50k` per construct **[−0.1366, −0.0310]**, i.e. candidate genes sit about
**7.6 % below** their matched non-candidates. The same observation survives two structurally
different nulls (a real-spacer resampling null, p = 0.0015; a wrong-recipient construct permutation,
p = 2.0e-4) and a partner-construct specificity test: the identical candidate genes read at the
wrong-recipient level (+0.041 [−0.295, +0.369]) in the cells of the **partner construct that knocks
down the same intended gene with different protospacers**.

**The verdict, nonetheless, is not PASS.** Two of the seven preregistered known-answer rows failed —
KA-4's no-signal row measured a small but real residual bias of **+0.068** (rule: ≤ 0.050) and KA-6's
matching balance missed on promoter-window count, SMD **+0.291** (rule: ≤ 0.250). The frozen rule
says a known-answer failure aborts the conclusions, and §7 maps that case to
`UNRESOLVABLE-AT-THIS-RESOLUTION`. The standard is not weakened after the fact, so the result is
reported as a **bounded estimate with an interval, not a verdict**. The bias is ~9 % of the
observed effect and the bias-immune permutation reading is unchanged, so this is a formal abort on
a marginal calibration miss, not a retraction of the estimate.

**What a positive would and would not establish — stated here, not only in the limitations.** Even a
clean PASS would say only that genes selected by a 12-nt-seed + NGG promoter match are, in
aggregate, more repressed in the cells carrying that construct than matched genes the rule did not
select. It would **not** establish dCas9-KRAB off-target binding. Expression here is observational;
the 34 constructs differ in on-target knockdown (`kd_lfc` −3.04 … +0.03), in cell number (95 …
3,415), in population composition and in handling; and the statistic is computed on one experiment's
cells with no biological replicate of any construct. The partner-construct test excludes
intended-target downstream biology and generic candidate-gene properties as explanations, and
narrows what is left to **something carried by the protospacer sequence itself** — which includes
seed-mediated binding but also any construct-level artefact correlated with spacer identity. It
cannot separate those. A negative would have meant "no aggregate expression corroboration under this
test", never "off-targets do not occur".

---

## THE LINE

**C8 (aggregate candidate-repression enrichment, construct as unit):
UNRESOLVABLE-AT-THIS-RESOLUTION** — the pass rule's two conditions were both met (`p_joint` =
9.9999e-06; max leave-one-construct-out p = 1.99998e-05) but 2 of 7 preregistered known-answer rows
failed (KA-4 residual bias +0.0683 vs a ≤ 0.05 rule; KA-6 promoter-window SMD +0.2905 vs a ≤ 0.25
rule), and §9's frozen abort clause governs. **Bounded estimate: S = +0.786 [+0.387, +1.166]
standardized units per targeting construct (n = 34), Δ = −0.0790 [−0.1366, −0.0310] log1p-CP50k;
null: selection-matched candidate-set resampling, 100,000 joint draws, null mean +0.00005, sd
0.17173, 95 % [−0.334, +0.335]; attainable floor p = 1.0e-5.**

---

## 1. Population and unit (verified, not assumed)

| | |
|---|---|
| **Experimental unit** | the gene-targeting construct (programme rule 1) |
| **Primary population** | **all 34** targeting constructs — 13 failing + 21 working by 293's `\|kd_lfc\| < 1.0`, covering **17 paired intended targets**, each measured only in its own split (21 train / 7 validation / 6 test) |
| **Candidate-bearing at L = 12** | **8 constructs** — `CAST_P1`, `GRB10_P2`, `KDM2B_ENST00000377071.4`, `KLF10_P1`, `DPF2_P1`, `IGF1R_P1`, `MED15_P1`, `TMSB4X_P2`, carrying exactly 305's 11 primary pairs. **The coordinator's count of 8 is confirmed** (305's §1/§8 said 9) |
| **Candidate set used** | 301's own `L*` rule — the hit set at the longest seed length in {12, 11, 10, 9} that yields a hit in the CRISPRi active window after the 10 kb distal filter. L\* = 12 for 8, 11 for 16, 10 for 7, 9 for 3; **all 34 have a non-empty set**, so no construct is dropped for lacking a strong candidate |
| **Candidate genes** | 83 construct–gene pairs over the 34 constructs (k = 1 for 19 constructs, median 1, max 18 for `STAG2_P1`) |
| **Secondary population** | the 31 non-targeting constructs, run identically; one has no hit at any L and contributes z = 0 by the declared rule |

## 2. Known-answer rows

| id | row | result | rule | |
|---|---|---|---|---|
| KA-1 | 305's 23 published rows (pooled effect, both half effects, control CP50k) recomputed from `sums305_<split>.npz` by code written here | max \|diff\| **0.000e+00** on all four columns | ≤ 1e-9 | **PASS** |
| KA-2 | 301's selection rule reproduced from its saved hits | failing **4/13**, working **4/21**, NTC **4/31**; 8 L=12 targeting constructs; the 11 pairs are 305's | exact | **PASS** |
| KA-3 | the 31 NTC constructs' sums across the three splits | identical names; max \|sum diff\| **0**; max \|cell-count diff\| **0**; 38,176 control cells | exact | **PASS** |
| KA-4 | **no-signal row** — each construct scored against *another* construct's candidate set, 200 within-split derangements | mean **+0.0683**, 95 % [−0.2406, +0.3769] | \|mean\| ≤ 0.05 **and** interval brackets 0 | **FAIL** (interval condition met; mean missed by 0.018) |
| KA-5 | null calibration — 200 null draws each p-valued against the rest | KS uniformity p = **0.3682**, fraction ≤ 0.05 = **0.0600** | KS p ≥ 0.05, fraction ∈ [0.02, 0.10] | **PASS** |
| KA-6 | matching balance, real candidates vs their matched pools | `pi_off` **+0.001**, basal **−0.027**, GC **−0.129**, TSS density **+0.199**, mappability **+0.242**, low-complexity **+0.213**, mean offset **+0.082**, n_sites **+0.114**, **n_windows +0.291** | \|SMD\| ≤ 0.25 on every axis | **FAIL** (one axis of nine) |
| KA-7 | 293's on-target `kd_lfc` for all 34 constructs recomputed from 305's sums | max \|diff\| **2.220e-16** (worst `DAXX_P2`) | ≤ 1e-9 | **PASS** |

KA-1 comes back **bit-identical** from the saved integer sums. That is not a stronger check than
the coordinator's independent audit of 305 (2.22e-16): that audit deliberately avoided importing
305's helpers, whereas KA-1 re-implements the same `prof` arithmetic in the same summation order, so
exact equality is what it should produce and a non-zero difference would have meant a bug here.
KA-7's 2.22e-16 is the honest cross-implementation figure — it compares against 293's differently
ordered float64 accumulation. KA-3 confirms 305's observation and
tightens what it implies: **the 38,176 control cells are the same cells in all three files**, so the
control pool is shared *between* splits (never between the two halves of one comparison), every
construct's effect in every split is referenced to one and the same control transcriptome, and the
three splits are therefore not three independent experiments for any control-referenced statistic.
The primary statistic is a **within-construct rank**, which is why it is largely insensitive to
this: a shift common to all of a construct's genes cancels.

**Why the two failures matter, and how much.** KA-4 says the machinery carries a small
construct-independent bias: any construct's candidate genes read slightly repressed in *other*
constructs' cells. Extended post-hoc to 5,000 derangements it is **+0.0717 ± 0.0022** — real, and
**9 % of the observed +0.786**. KA-6 says the matching leaves genes with more promoter windows
over-represented among real candidates. Both are marginal misses of thresholds fixed in advance;
neither is repaired by loosening the threshold, and neither was loosened.

## 3. The null and its attainable resolution

Genes are not resampled independently; **candidate sets are**. For construct `c` at level `L = L*(c)`,
each real candidate is replaced by a gene drawn from its own matched pool:

- **hard strata** — 10 basal-expression deciles × 5 `pi_off_L` quintiles (50 cells);
- **`pi_off_L`** is *measured*, not modelled: all 209,070 hCRISPRi-v2 protospacers (the library the
  2025 guides come from; **68/68 real spacers found verbatim in it**) were rescanned through 301's
  exact rule, own-target excluded and the 10 kb distal filter applied, giving the empirical
  probability that a real library spacer hits each gene off-target at each seed length. Calibration:
  per-spacer P(≥ 1 off-target hit) at L = 12 is **0.132 for the 68 real 2025 spacers vs 0.159 for the
  pool** (0.412/0.389 at L = 11, 0.750/0.694 at L = 10, 0.956/0.920 at L = 9);
- **soft matching** inside the cell — the 200 nearest genes on GC, TSS density, mappability proxy,
  low-complexity proxy, mean TSS offset, promoter-window count and anchor count;
- **exclusions (the ones 305's null omitted)** — every gene any 2025 construct intends to target
  (300 genes), the construct's own target, and every gene its own spacers hit at any seed length.
  Eligible universe: 14,132 genes at L = 12, 16,263 / 16,614 / 16,636 at L = 11 / 10 / 9. Mean pool
  size 196.7 genes (minimum 134); **no pool needed the declared widening rule**.

**Resolution.** `B = 100,000` joint draws ⇒ smallest attainable p = **1/(B+1) = 1.0e-5**, three
orders of magnitude below α = 0.05. This is the design difference from 305, whose same-split rank
test had a Holm floor of 11/158 ≈ 0.0696 and therefore could not reject *for any effect however
large*. Here a FAIL would have been informative. Null distribution of `S`: mean +0.00005, sd 0.17173,
95 % [−0.33393, +0.33489]; KA-5 confirms the machinery is calibrated (KS p = 0.368).

## 4. The ONE primary statistic

`q(G)` = mid-rank percentile of `e_c(G)` within `G`'s matched pool (lower = more repressed);
`r_c` = mean over the construct's candidates; `z_c = sqrt(12·k_c)·(0.5 − r_c)`;
**`S` = mean of `z_c` over the 34 targeting constructs.**

| quantity | value | unit |
|---|---|---|
| `S` observed | **+0.78567** | per targeting construct |
| construct-level bootstrap 95 % CI (20,000 resamples) | **[+0.38669, +1.16604]** | n = 34 targeting constructs |
| joint-null p (one-sided, repression) | **9.9999e-06** (= the floor) | 100,000 joint draws |
| leave-one-construct-out max p | **1.99998e-05** (deleting `GRB10_P2`) | 34 deletions |
| `Δ` effect size | **−0.079012** log1p-CP50k, 95 % **[−0.136647, −0.030962]** | per targeting construct |
| candidates / constructs | 83 genes / 34 constructs | |

Both pass conditions of §7 are satisfied. **The verdict is nevertheless
`UNRESOLVABLE-AT-THIS-RESOLUTION`** because §9's abort clause fired. `Δ = −0.079` on the
log1p-CP50k scale is roughly a **7.6 % reduction** relative to matched genes — an aggregate over
83 candidates, not a per-gene claim.

## 5. Declared secondaries (no verdict)

| # | secondary | result |
|---|---|---|
| **S1** | **failing vs working constructs — 301's original T7 question, at the aggregate level** | failing (n = 13) `z` = **+0.7383**, working (n = 21) **+0.8150**; difference **−0.077 [−0.869, +0.699]**, Mann–Whitney p = 0.595. **No contrast.** The aggregate candidate repression is present in *both* classes; T7's premise — that off-targeting is what distinguishes failing constructs — remains unsupported, now at the aggregate level as well as 301's 4/13-vs-4/21 count level |
| **S2** | the 31 non-targeting constructs, identical pipeline | `S_NTC` = **+0.2857 [−0.0691, +0.6358]**, joint-null p = **0.0537**. Nominally positive, not separated from the null. `PRKCQ` (z = +1.732, p = 1e-5) and `DCTN2` (z = +1.697, p = 0.0102) are the two 305 observations; `CDCP1` (z = +1.680) is a third. As 305's review said, these are candidate-bearing observations selected by the same rule, **not negative controls** — and their aggregate sits at about a third of the targeting constructs' |
| **S3** | `GRB10_P2 → SNX8` with its partner | `e` = **−0.492396**, matched-pool percentile **0.0000** (the most repressed of its 200-gene matched pool); the partner `GRB10_P1`, same intended target, different protospacers, reads **+0.045325** on the same gene. `GRB10_P2` overall z = **+2.418** (p = 1.1e-4); its partner scored on the same set: **−1.922** |
| **S4** | measurement floor | **33 of 83** candidates sit below 305's RAW CP50k ≥ 0.5 floor. `S` on above-floor candidates only = **+0.6907** — the result is not a low-expression artefact. Per-construct cells 95 … 3,415 (`DPF2_P1` is the 95-cell construct 305's review flagged) |
| **S5** | **second, independent null — real-spacer resampling.** A pseudo-construct is a random pair of stringency-0 hCRISPRi-v2 protospacers rescanned through 301's exact rule, kept only if it yields `L* = L*(c)` and exactly `k_c` candidates after `c`'s exclusions; scoring is matching-free (percentile in the whole eligible universe) | **all 34 constructs resolvable**. `W_obs` = **+0.74470** vs null mean −0.01033, sd 0.17883, 95 % [−0.360, +0.340]; **p = 0.00149** (672 joint draws, floor 1.49e-3). This null preserves the sequence rule *jointly* and does **not** match expression — the complement of the primary null, and it agrees |
| **S6** | by seed length | L\* = 12 (n = 8) `S` = **+0.633**, p = 0.0366 · L\* = 11 (n = 16) **+1.014**, p = 3.0e-5 · L\* = 10 (n = 7) **+0.496**, p = 0.094 · L\* = 9 (n = 3) **+0.652**, p = 0.129. **The effect is not concentrated at L = 12.** A strictly seed-length-graded mechanism predicts the opposite ordering; this is the single most awkward fact for a simple seed-binding reading, and it is reported as such |
| **S7** | is the 2025 library unusually promiscuous? | observed L\* counts [8, 16, 7, 3] vs the hCRISPRi-v2 pseudo-construct expectation [10.1, 11.3, 9.6, 3.0]; χ² = 3.13, **p = 0.371**. **No** — the 2025 constructs carry a typical number of seed matches for their library |
| **S8** | KA-4's cross-construct matrix as a scale for `z_c` | mean +0.068, sd 0.157 over 200 derangements; extended in §6 |

## 6. Post-hoc, explicitly labelled, carrying no verdict

None of these can change §4's verdict; they exist because the mandate asks what the data can and
cannot resolve.

- **P1 — wrong-recipient permutation, 5,000 derangements** (a Monte-Carlo extension of the
  preregistered KA-4 statistic, not a new statistic). Null mean **+0.07168** (SE 0.00222), sd
  0.15665, 95 % [−0.2367, +0.3840]. `S_obs` = +0.78567 ⇒ **permutation p = 2.0e-4** (0 of 5,000
  derangements reached it; floor 2.0e-4). **This reading is immune to the matched-pool
  construction**, because any matching bias enters the observed and the permuted statistic alike.
- **P2 — repaired matching.** Adding promoter-window count to the hard strata fixes the axis KA-6
  failed on (SMD +0.2905 → **0.0000**) and changes nothing: `S` = **+0.78813**, joint-null
  p = 1.0e-5, residual KA-4 bias **+0.081**. The KA-6 failure is therefore not the cause of the
  estimate; the residual bias is not removed by better covariate balance either.
- **P3 — partner-construct specificity** (a generalisation of declared S3 to all 34 constructs).
  Each construct's candidate genes re-scored in the cells of its **partner construct — same intended
  target, different protospacers**: `z` = **+0.04106 [−0.29505, +0.36933]**, versus **+0.78567** in
  their own construct's cells. Own − partner = **+0.74461 [+0.22730, +1.24514]**, paired Wilcoxon
  p = **0.00459**, n = 34 paired constructs. **This is the strongest single result here.** The
  partner knocks down the same gene, in the same split, referenced to the same controls; the
  candidate genes are the same genes; only the protospacer sequence differs. It excludes
  intended-target downstream biology and generic candidate-gene properties, and leaves the
  association attached to the spacer sequence.
- **P4 — sign count.** 24 of 34 targeting constructs have `z > 0` (exact binomial p = 0.0122);
  median `z` = +1.221. Top: `HIRA_P2` +2.429, `GRB10_P2` +2.418, `SOX4_P2` +2.262, `SMARCB1_P1`
  +2.168. Bottom: `TMSB4X_P2` −1.732, `MED15_P2` −1.386, `SOX4_P1` −1.212. (For a k = 1 construct
  `z` is bounded in [−1.732, +1.732].)
- **P6 — 305's two cell halves.** The construct's cells are split deterministically (batch-balanced,
  no RNG) and referenced to **disjoint control constructs** — half 0 to control half A (16 NTC
  constructs, 21,581 cells), half 1 to control half B (15 NTC constructs, 16,595 cells). Half 0
  gives `S` = **+0.75776 [+0.38759, +1.09989]**, half 1 **+0.69101 [+0.28470, +1.07550]**, against
  the pooled +0.78567. Per construct the two halves agree at Pearson **+0.718** (p = 1.8e-6),
  Spearman +0.726, same sign in **27/34**, both positive in 23/34. This is a **cell-sampling**
  split, not a biological replicate — 305's preregistration says so and this one repeats it — so it
  shows the aggregate is a stable feature of the measurement rather than an accident of which cells
  and which control constructs were used. It cannot turn one experiment into two.
- **P5 — bias-corrected bounded estimate.** `S` − (P1's wrong-recipient level) = **+0.71399**,
  construct-bootstrap 95 % [+0.31501, +1.09436], n = 34 targeting constructs.

## 7. What this does not establish, and what would

**Not established.**
1. **Off-target binding.** Nothing here is an occupancy measurement. The association is with the
   protospacer *sequence*, not with a demonstrated dCas9-KRAB binding event.
2. **A mechanism graded by seed length.** S6 shows the L\* = 11 constructs carrying the largest
   aggregate, L\* = 12 a smaller one, L\* = 10 and 9 unresolved. A seed-binding model predicts the
   opposite; either the 12-nt threshold is not the operative one, or part of the aggregate is
   something else that travels with spacer identity.
3. **A failing-construct explanation.** S1 finds no failing-vs-working contrast (−0.077
   [−0.869, +0.699]). T7 stays where 301 left it: **off-targeting is neither proved as the
   failing-construct mechanism nor excluded.**
4. **Independent replication.** Every number rests on one experiment's cells, the three splits share
   one control pool (KA-3), and no construct has a biological replicate. The construct-level
   bootstrap resamples constructs, not experiments.
5. **A population off-target rate.** 83 candidates over 34 constructs and 17 intended targets cannot
   support one.

**What would resolve it, and what it would cost.**
- **A sequence-mutant control (the decisive cheap experiment).** Re-clone each of ~8 candidate-bearing
  constructs with 2–3 silent seed mutations that abolish the predicted off-target match while
  preserving the on-target site, and read the candidate genes. If repression tracks the seed and not
  the construct, the mechanism is established; if it survives the mutation, the association is a
  construct artefact. ~24 constructs × ~500 cells, one 10x lane, order **$3–6 k** plus cloning.
  This is the only design here that separates hypothesis 1 from hypothesis 3.
- **An occupancy assay** (dCas9 ChIP-seq / CUT&RUN on 3–4 of these constructs) tests binding
  directly rather than its transcriptional shadow; order **$5–10 k**, and it answers "is it bound",
  not "is it repressed".
- **An independent screen carrying the same library.** 285/293 already establish that no source on
  disk shares the 2026 constructs, and 293 measured the 2025 library's own two-construct structure;
  a published CRISPRi screen using hCRISPRi-v2 in another line, with ≥ 400 cells per construct, would
  let the identical statistic be recomputed with **cell line as a second unit** — free if such a
  screen exists (306's remit), otherwise not commissionable at this budget.
- **More constructs, not more cells.** With the construct as the unit the interval is set by n = 34.
  Doubling cells per construct does not narrow it; 34 → 100 candidate-bearing constructs would.

## 8. Deviations from the preregistration

| # | deviation | before or after any outcome |
|---|---|---|
| 1 | KA-1's control-expression comparison initially used `log1p` CP50k where 305's `ctrl_cp50k` column is **raw** CP50k, and did not apply 305's leave-that-construct-out reference for tested NTC rows. Fixed; the row then reproduced at 0.000e+00. The pooled and half-effect comparisons — KA-1's substance — reproduced at 0.000e+00 both times | **before** (no candidate statistic had been computed) |
| 2 | KA-4, KA-5 and KA-6 were computed in the same script invocation as the primary rather than in a prior one. The code for all of them was written before the script was first run, so no threshold or definition was chosen after seeing a number; but the outputs were seen together | **stated for the record**; no rule was changed afterwards |
| 3 | The KA-4 derangement count was extended from the frozen 200 to 5,000 for the P1 reading. The frozen 200-draw row is reported unchanged in §2 and is the one that governs the verdict; P1 is labelled post-hoc | **after** (Monte-Carlo resolution only, same statistic) |
| 4 | P2, P3, P4, P5 are post-hoc analyses added after the outcome; all are labelled and none carries a verdict | **after** |
| 5 | `seq307.py`'s `pi12…pi9` columns in `genes307.csv` are superseded: they were computed before the Excel-escaping bug in Horlbeck's gene labels was found (`="SYMBOL"`), which stopped the own-target exclusion firing and inflated the per-spacer L = 12 hit rate from 0.107 to 0.597. **`pool307.py`'s `pi_off` is the propensity used everywhere**; the bug was found and fixed by the pre-hash calibration check, and no analysis used the superseded columns | **before** (pure sequence work, pre-hash) |

## 9. Files

| file | what |
|---|---|
| `PREREGISTRATION.md`, `work/PREREGISTRATION.sha256` | the frozen design, hashed 16:00:01Z |
| `work/seq307.py`, `constructs307.csv`, `genes307.csv`, `seq307_meta.json` | 301's rule reproduced; `L*(c)` and candidate sets for all 65 constructs; per-gene covariates |
| `work/pool307.py`, `poolhits307.pkl`, `pool307_meta.json` | the 209,070-spacer rescan, `pi_off`, the pseudo-construct distribution |
| `work/probe307.py`, `probe307.json` | the real-vs-pool spacer calibration that found the label bug |
| `work/lib307.py`, `ka307.py`, `ka307.json` | 305's sums reader and KA-1 / KA-3 / KA-7 |
| `work/enrich307.py`, `enrich_pools307.py`, `results307.json`, **`enrichment307.csv`** | the primary statistic, the joint null, KA-4 / KA-5 / KA-6, S1–S4, S6, S7; per-construct table |
| `work/sens307.py`, `sens307.json` | S5, P1, P2 |
| `work/extra307.py`, `extra307.json` | P3, P4, P5 |
| `work/halves307.py`, `halves307.json` | P6 |

`work/enrichment307.csv` carries, per construct: unit, group, intended target, split, `kd_lfc`,
cells, `L*`, candidate count and names, mean matched-pool size, pools widened, null draws, `r_c`,
`z_c`, `Δ`, the null distribution's mean / sd / 2.5 % / 97.5 %, the per-construct p, the number of
candidates below the measurement floor and the minimum candidate CP50k.

No model, no submission, no production change. Nothing was written outside this directory except the
C8 entry in `docs/CONFIRMATION-LEDGER.md`.
