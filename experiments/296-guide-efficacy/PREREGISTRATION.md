# 296 — PREREGISTRATION: is F14's "construct term" guide EFFICACY?

Written and hashed (`work/PREREGISTRATION.sha256`) **before any stratified number was computed**. Only
metadata was read first: file sizes, `obs` column names and label counts, the 280 bank's guide/target index,
the CD4 supplementary table's column names, 268's `halves/*.npz` key list, and the already-published numbers
in `E285/RESULTS.md`, `E293/RESULTS.md`, `E285/work/*.log`. No expression value, no cosine, no knockdown
value was read before this file was hashed.

Mandate: `experiments/296-guide-efficacy/MANDATE.md`. Analysis only; no build, no submission, nothing written
outside `experiments/296-guide-efficacy/`.

---

## 0. The question and the two outcomes

285 (F14) measured the **guide construct** as the dominant irreducible term of a perturbation response
(ρ_guide 0.09 K562-GWPS / 0.17 Nadig Jurkat / 0.41 Nadig HepG2 / 0.02–0.06 CD4 / 0.76–0.81 Broad-16 Cas9).
293 measured the same term in Arc's own 2025 library and found it **bimodal in the weaker construct's
on-target knockdown**: +0.687 where both constructs reach |lfc| ≥ 1.0, −0.003 where the weaker does not.

* **If** the public libraries split the same way → F14's headline becomes *"guide **efficacy** is the
  dominant irreducible term"*: a validated construct transfers, a failing one is a reproducible response
  to the vector.
* **If** they do not → the construct term is construct-specific biology and **F14 stands as written**.

Both outcomes are publishable. A well-supported negative is the result.

---

## 1. Known-answer rows (must reproduce before any new number is reported)

| id | row | source | rule |
|---|---|---|---|
| **KA-1** | 285's K562 guide row at n = 25, 748 pairs: (b) +0.0006, (a) +0.0063, `pds` (b) 0.5093 / (a) 0.6330 | `E285/work/k562_guide.{npz,json}` | medians recomputed from the SAVED per-target arrays must equal the published values to **0.0e+00** (they are the same arrays); the `pds` values re-read from `k562_guide.json` |
| **KA-2** | 285 §6b's three knockdown terciles (T1 +0.0009 / T2 +0.0003 / T3 +0.0021 guide; +0.0072 / +0.0064 / +0.0089 cells; ρ 0.124 / 0.041 / 0.236; weaker-KD medians +0.015 / −0.119 / −0.300; n 180 / 179 / 180) | `E285/work/k562_kd.{npz,json}` | the tercile rule re-applied to the saved `kd`, `guide`, `cells` arrays must reproduce every cell at **0.0e+00** |
| **KA-3** | 293's ALL / both-≥1.0 / weaker-<1.0 rows at n = 200, challenge depth: ρ_guide +0.075 / +0.687 / −0.003; (a) +0.312 / +0.150 / +0.317; (b) +0.023 / +0.103 / −0.001; k = 16 / 5 / 11 | `E293/work/construct293.csv` | the stratification rule re-implemented **here** (not by importing `strat293.py`) must reproduce every cell to **≤ 1e-9** |
| **KA-4** | union identity of halves: Broad-16 `A + B == full` exactly and `v1 + v2 ≤ full` on all 12 lines | `E268/work/halves/*.npz` | max|A+B−full| = 0 and max(v1+v2−full) ≤ 0 |
| **KA-5** | the CD4 bank's guide rows join 1:1 to the authors' published per-guide knockdown table | bank `obs280`, `guide_kd_efficiency.suppl_table.csv` | all bank targeting guides present in the table; my bank-computed on-target |lfc| must correlate with the authors' own `t_statistic` at Spearman > +0.3 (a cross-estimator check, not a threshold on the science) |
| **KA-6** | K562-GWPS on-target medians: per-construct median −0.320, weaker-construct median −0.119, both negative on 79.2 % of 539 pairs | `E285/work/k562_kd.{npz,log}` | recomputed here from the raw file at n-independent full depth must land within ±0.02 of the saved values (the construct cap differs, see §2.1) |

**Abort rule.** If KA-1, KA-2, KA-3 or KA-4 does not reproduce, the run stops and the failure is the result.
KA-5 and KA-6 are cross-estimator checks: a miss is reported and the affected library is labelled, not
silently dropped.

---

## 2. The frozen measurement

### 2.0 The response object (unchanged from 285/293 — nothing mathematical is new)

`prof(S) = log1p(5e4 · S / ΣS)`; `delta = prof(side) − prof(control)`; residual = delta minus the **matched**
panel mean on the same side; per-target cosine on the **scored** columns (every construct's own target gene
column removed from every scored statistic); `pds` = 218's `_per_target`, imported verbatim through
`E285/work/lib285.py` (273 KA-2: equals `cell-eval2` 0.16.0's `pds_cosine` to 5.6e-17). Seed 20260906,
8 draws, paired target bootstrap 20,000 draws — all `lib285.py` constants.

### 2.1 The five libraries, their construct pairs and their budgets

| library | modality | construct unit | pairs | rungs / budgets | control |
|---|---|---|---|---|---|
| **K562-GWPS** `data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad` | CRISPRi | `obs/guide_id` inside one `obs/perturbation` | 775 two-construct perturbations (748 usable at n = 25) | n = 25 and 50 cells/side | `perturbation == "control"` cells |
| **Nadig Jurkat** `data/external/nadig2025/jurkat_raw_singlecell.h5ad` | CRISPRi | `obs/sgID_AB` inside one `obs/gene` | ~124 at n = 25 | n = 25 (and 50 if ≥ 30 targets qualify) | `gene` starting `non-target` |
| **Nadig HepG2** `…/hepg2_raw_singlecell.h5ad` | CRISPRi | as above | ~70 at n = 25 | n = 25 (and 50 if ≥ 30 qualify) | as above |
| **Broad-16** `E268/work/halves/<line>.npz` | **Cas9 KO** | the `v1` / `v2` guide sums | ~97 targets × 12 lines | T = 2×10⁵ and 5×10⁵ UMIs/side | the 10 OR knockouts summed (268's `CONTROLS`) |
| **CD4 GWCD4i** `E280/work/{bank280_X.npy, obs280.npz}` | CRISPRi | guide slot 0 vs slot 1 inside one (target, donor, condition) | 2,969 (cond, donor, target) cells with two guides | T = 1×10⁵ and 2×10⁵ UMIs/side | the (cond, donor) NTC pool |

Memory: every raw pass preallocates its output (no `concatenate` doubling), caps cells per construct at
**300** and control cells at **25,000** (K562) / **20,000** (Nadig); one raw pass at a time; the K562 and
Nadig extractions are the only passes above 4 GiB and they run serially, gated on `memory_pressure ≥ 40 %`.
The cap is the only design difference from 285's K562 pass and is the reason KA-6 is a ±0.02 band, not 0.

### 2.2 The five statistics, per target (computed ONCE per library on the full cohort, then stratified)

Exactly 293's construction, so that the strata share one centring and one distractor set:

* **(a)** the **larger** construct's two disjoint cell halves (or, where only sums exist, its binomial count
  split) — the matched cell-sampling reference.
* **(a2)** the **smaller** construct's own two halves — its own reliability.
* **(b)** larger construct's side A vs the smaller construct's cells — the construct term.
* **(x)** the cross-target null: side A of target *t* vs the smaller construct of the **next** target in the
  cohort (a different gene).
* **ρ_guide = median(b) / median(a)**; **ρ_dis = median(b) / √(median(a)·median(a2))**.
* the scorer's **`pds`** for (a), (b), (x), every cohort target passed as a distractor.

Bootstraps: `lib285.boot_median` / `boot_median_diff` (paired over targets where the strata are the same
targets). Between two **disjoint** strata a paired bootstrap does not exist; the difference
ρ_guide(both) − ρ_guide(weaker fails) is bootstrapped **two-sample**, resampling targets independently
within each stratum, 20,000 draws — labelled *stratum bootstrap* everywhere it appears.

### 2.3 The knockdown measure and the cuts (frozen)

**Measure.** Each construct's on-target log delta at **its own gene's column**, at **full depth** (every cell
of that construct), against the library's control profile, on the `log1p(5e4·p)` scale — 293's definition.
That column is **excluded from every scored statistic**, so no stratification here is circular.
`|lfc|` = its absolute value. (A positive on-target delta = no knockdown at all.)

**Cuts, applied to the WEAKER construct of each pair** (`kd_weak = min(|lfc_big|, |lfc_small|)`):

1. **C1 — 293's absolute cut**: strata `both constructs |lfc| ≥ 1.0` vs `weaker construct |lfc| < 1.0`.
2. **C2 — library-relative cut**: the same, with the cut set to **the library's median single-construct
   |lfc|** (median over every individual construct measured in that library).
3. **C3 — expression-free robustness cut** (declared here, before any value is read, because the
   `log1p(5e4·p)` scale compresses low-expressed genes and the public libraries are not expression-selected
   the way the 2026 panel is): knockdown **fraction** `1 − p_construct / p_control` at the gene's own column,
   cut at **0.5** (≥ 50 % knockdown).

C1 and C2 are the headline (the mandate's two). C3 is reported beside them; if C1 and C3 disagree the
disagreement is the finding and both are reported.

### 2.4 Frozen predictions and their decision rules

| | prediction | rule (stated before the measurement) |
|---|---|---|
| **P1** | in **every CRISPRi library** (K562-GWPS, Nadig Jurkat, Nadig HepG2, CD4) ρ_guide in the both-≥cut stratum exceeds the weaker-<cut stratum, interval above 0 | PASS for a library iff, at cut C1, ρ_guide(both) − ρ_guide(weaker fails) has its **2.5th stratum-bootstrap percentile > 0**. P1 as a whole PASSES only if all four CRISPRi libraries pass; libraries with < 8 targets in either stratum are reported as **UNDERPOWERED**, not as failures, and P1 is then stated over the libraries that can be tested |
| **P2** | the both-≥1.0 ρ_guide is **≥ 0.4** in K562-GWPS (293's +0.687 minus its width) | PASS iff the point estimate ≥ 0.40 at n = 25, cut C1 |
| **P3** | **Broad-16** (Cas9, validated guides) shows **no such split** | PASS iff P1's rule FAILS on Broad-16 (interval on the difference includes 0) at T = 2×10⁵, cut C1 |
| **P4** | the weaker construct's knockdown predicts ρ_dis at **Spearman > +0.3** in each CRISPRi library | per-target ρ_dis = b_t/√(a_t·a2_t) where a_t > 0 and a2_t > 0; Spearman of `kd_weak` (as a **strength**, i.e. |lfc|) against ρ_dis; PASS per library iff > +0.30. The companion Spearman against (b) itself is reported beside it (293 reported +0.626 / +0.644) |
| **P5** | **failing constructs are reproducible and share a vector signature**: (i) constructs with |lfc| < 0.5 have same-construct split-half reliability > 0 with the bootstrap interval above 0; (ii) their cross-construct agreement (failing × failing, **different genes**) exceeds the working constructs' (x) null | (i) PASS per library iff the 2.5th percentile of the median split-half > 0; (ii) PASS iff median cos(failing_i, failing_j) − median cos(working_i, working_j) has its 2.5th stratum-bootstrap percentile > 0, reported both in residual space (strict: the library mean already removed) and in un-centred delta space (permissive: includes the mean) |

Every outcome is stated in `RESULTS.md` whether it passes or fails. No prediction is added, removed or
re-worded after a number is seen.

---

## 3. The failing-construct phenomenon (§3 of the mandate)

Per library, over constructs with **|lfc| < 0.5** ("failing") against constructs with **|lfc| ≥ 1.0**
("working"):

1. same-construct split-half reliability (its own (a) or count split) at the library's rung;
2. cos(failing construct's response, its **working partner's** response) — the (b) row of those pairs;
3. cos(failing construct's response, the **library's mean response** over all its targets), **un-centred**
   (the residual removes the mean by construction), vs the same cosine for working constructs — is the
   failing response a generic stress / vector signature?
4. cos across **failing constructs of different genes** vs cos across **working constructs of different
   genes** — the (x) null restricted by stratum — in both residual and un-centred space.

---

## 4. Library failure rates (§4 of the mandate)

Share of **single-construct** targets with |lfc| < 1.0 and < 0.5, per library, beside 293's 18.7 % / 4.3 %
for the 2025 VCC panel:

* K562-GWPS: **2026-panel targets** (`data/challenge/pert_counts.csv` ∩ GWPS) and **non-panel** separately;
* Nadig Jurkat, Nadig HepG2;
* Broad-16 (per line and pooled; Cas9, so a "failure" is a target whose mRNA does not drop);
* CD4: both my bank-computed |lfc| **and** the authors' own `signif_knockdown` share from
  `data/external/cd4_cell2026/suppl_tables_2026-08-24/guide_kd_efficiency.suppl_table.csv`.

Plus, **indicative only**: the fraction of the 300 2026-panel targets whose **K562-GWPS** construct fails
(the 2026 constructs are different vectors in a different line; this is a base rate, not a prediction about
the challenge's own truths).

---

## 5. Descriptive tail (POST-HOC, labelled as such in `RESULTS.md`)

Whether promoter choice / transcript-level constructs (293's `KDM2B`, `_P1P2` stems) explain any of the
split; the 20 most and least efficient constructs per library. Anything discovered while looking is labelled
POST-HOC and carries no frozen rule.

---

## 6. Deliverables

`PREREGISTRATION.md` (this file) + `work/PREREGISTRATION.sha256`; `RESULTS.md` with the KA rows, the
stratified table per library at both headline cuts, P1–P5 with their rules and outcomes, the
failing-construct phenomenon, the failure-rate table, the descriptive tail; `work/efficacy296.csv`
(one row per **(library, target, construct)**: cells, |lfc|, (a), (a2), (b), (x), ρ_guide, ρ_dis, `pds`);
and the one-line verdict:

> **GUIDE EFFICACY EXPLAINS THE CONSTRUCT TERM: YES / PARTLY / NO** (K562 both-≥1.0 ρ_guide … vs
> weaker-<1.0 …; Jurkat …; HepG2 …; Broad-16 …; CD4 …); **FAILING CONSTRUCTS REPRODUCIBLE: YES / NO**
> (split-half …; shared vector signature: YES / NO); **FAILURE RATES: …**

Verdict rule, frozen: **YES** if P1 passes in every testable CRISPRi library and P2 passes; **PARTLY** if
P1 passes in at least half of the testable CRISPRi libraries; **NO** otherwise.
</content>
</invoke>
