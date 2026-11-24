# 282 — PREREGISTRATION: how conserved is a perturbation response across activation states of the SAME donors' CD4 T cells?

Written 2026-09-06, experiment 282, **before any across-state cosine, any state-dependence index, any family
contrast or any decomposition number was computed**. Everything read up to this point is *design metadata only*:
`obs280.npz`'s `.obs` columns (donor / condition / guide / target / cell counts / `keep_for_DE`), the HGNC
gene-group table, the *names* of the columns of `conservation276.csv`, `cons276.npz` and `znf279.csv`, and the
counts printed by `work/meta282.py` (target availability, cells per side, guide-set stability, family sizes).
No expression value, no residual row and no cosine has been opened.

Nothing below may be changed after this file's sha256 is recorded (`work/prereg282.sha256`); anything added
later is labelled **POST-HOC** in `RESULTS.md`.

---

## 0. The question and the frozen scope

> Within one individual, how much of a perturbation response is preserved when the recipient cell changes
> **state** (Rest → Stim8hr → Stim48hr of the same donor's CD4 T cells), measured against the **across-donor
> reproducibility within a state** as the ceiling?

Every cross-line comparison in the programme confounds cell state with genotype, culture history and platform.
GWCD4i measures the same 4 donors' primary CD4 T cells in three activation states, one lab, one chemistry
(10x GEM-X Flex), ~500 cells per (target, state) pooled over donors. 280 established the across-donor ceiling
(+0.054 at ~250 cells/side, F10). This experiment measures the across-**state** number on the *same rows*, the
same geometry and the same instrument, so the two are directly comparable.

**Data.** 280's extracted bank only: `experiments/280-gwcd4i-donor-replicates/work/bank280_X.npy` (17,495 ×
18,129 f32, memmapped), `bank280_rows.npy`, `obs280.npz`, `axis280.npz`, `de_stats_panel280.csv`.
**No download, no new extraction, no build, no `vcc`, no `cell-eval2`, nothing written outside
`experiments/282-activation-state-dependence/`.** `extra280_X.npy` carries 500 random NON-panel targets but
**Stim48hr only** (measured: `extra280.py:9` filters `cond == "Stim48hr"`), so it cannot support any
across-state statistic; it is therefore **out of scope** for 282 and the analysis is the panel targets only.
(Recorded here so the mandate's "if `extra280` carries non-panel targets, those too" is answered with a fact.)

## 1. Objects — the residual effect rows

Exactly 280's arithmetic (`lib280.py`, `q1_280.py:38-53`), reused by import, not re-implemented:

* `p = log1p(5e4 · S / ΣS)` on the raw pseudobulk counts `S` of the selected rows;
* control `b` = the same transform of the **summed non-targeting rows of the same donors in the same state**
  (`Bank.ctrl_profile(cond, donors)`) — *each state is centred on its own controls*;
* `delta = p − b`, on the 17,771 mapped columns of the official 18,533-gene axis;
* residual `r = delta − mean(delta over the analysis target set, on the same side)` — *each (donor pair, state)
  matrix is centred on its own across-target mean*; the 300 panel positions are then zeroed and cosines are
  taken on the **17,471 scored columns**.

**Unit of a side.** A *donor pair* `P` ∈ {D1D2, D3D4, D1D3, D2D4, D1D4, D2D3} (D1 `CE0008162`, D2 `CE0010866`,
D3 `CE0008678`, D4 `CE0006864`), in one state `s` ∈ {Rest, Stim8hr, Stim48hr}. Each side is ~2 donors,
~215 cells (measured, `meta282.json`). The three *pairings* are the three disjoint splits
(D1D2|D3D4, D1D3|D2D4, D1D4|D2D3).

**Analysis target set.** The panel targets present with usable guides in **all 4 donors × all 3 states**, so
that every state's centring is over the same rows and every cosine is over the same targets.

**Guide filter — PRIMARY: the state-common set.** A `(target, donor, guide)` row is used only if it passes the
authors' `keep_for_DE` in **all three states**. Measured cost (metadata only): **241 targets** vs 242 under
280's per-state filter, and 490–511 pooled cells per (target, state) vs 536–587. Reason, declared in advance:
under 280's per-state filter the `keep_for_DE` guide set is identical across the three states for only
**154 / 242** targets, so an across-state cosine would partly measure a **guide-construct** swap — which
280/F1 show is a real unit of variance and which would bias state dependence **upward** (toward the
hypothesis). The state-common set holds the construct fixed on both sides of every comparison.
**SENSITIVITY arms, both reported: (a) 280's per-state `keep_for_DE` exactly (242 targets); (b) all guides
(273 targets).**

## 2. Statistics, frozen

Write `r(P, s, t)` for the residual row of donor pair `P`, state `s`, target `t`; `P̄` for the complement of
`P` within its pairing; `cos(·,·)` the cosine on the 17,471 scored columns.

* **(S1) across-state cosine, donor held fixed** — `cos(r(P, s_i, t), r(P, s_j, t))` for the three state pairs
  (Rest·Stim8hr, Rest·Stim48hr, Stim8hr·Stim48hr), averaged over the **6 donor pairs**. Reported per state pair
  and pooled over the three. Both sides come from the **same two people**, so donor variation is held fixed.
* **(S2) across-donor cosine within a state — THE CEILING and the known answer** —
  `cos(r(P, s, t), r(P̄, s, t))` over the 3 pairings, per state, and pooled over states. This is 280's §2.1
  statistic; its reproduction on 280's own filter and target sets is **KA-1**.
* **(S12) across-state AND across-donor** — `cos(r(P, s_i, t), r(P̄, s_j, t))`, both changed. The floor of the
  2×2, and the comparator that turns S1 and S2 into two attributable terms:
  `S2 − S12` = the value of sharing the **state**, `S1 − S12` = the value of sharing the **donors**.
* **(W) counting floor** — `cos` of two disjoint binomial half-splits of the *same* `(P, s, t)` counts
  (seeded, `lib280.thin` arithmetic). Half depth, so it is context, not a matched comparator; labelled as such.
* **(S3) state-dependence index** — `SD(t) = 1 − S1(t) / S2(t)` on targets with `S2(t) > 0`, where `S1(t)` and
  `S2(t)` are that target's pooled S1 and S2. **Headline = median SD with a paired target bootstrap**
  (B = 10,000, seed 20280906, percentile [2.5, 97.5], the resampled unit is the TARGET — `lib280.boot_median`).
  Declared robustness: the pooled ratio `SD_pooled = 1 − median(S1)/median(S2)` (stable where per-target `S2`
  is near 0) and `SD` on the subset with `S2(t) ≥ 0.05`.
* **(S4) family / class contrast** — `SD` by HGNC gene group (`data/external/xcell/hgnc_protein_coding.txt`,
  the table 279 used) for every group with ≥ 8 analysis targets, and by tercile of cross-line conservation.

**Cross-line conservation predictor `C_cross`.** 279's `work/znf279.csv` is the **held, NON-panel** gene set
(measured: **0 of the 300 panel targets appear in it**), so the panel's conservation comes from the same frozen
276 statistic on the panel: `experiments/276-conservation-map/work/`.
**PRIMARY `C_noCD4`** = 276's `C_raw` recomputed from `cons276.npz`'s per-pair cosines with **every pair
involving `cd4_rest` / `cd4_stim8hr` / `cd4_stim48hr` removed** (69 of 315 pairs), because those three sources
are derived from the *same GWCD4i cells* this experiment measures and would make the test partly circular;
296 targets keep ≥ 3 pairs. **SECONDARY `C_raw`** = 276's published column, unchanged (the F8 object).

**FROZEN PREDICTION (S4).** *Targets with low cross-line conservation have high within-individual state
dependence.* **Pass rule: Spearman(SD, `C_noCD4`) < 0 with a paired target bootstrap 95 % interval entirely
below 0, AND the same sign under the 280-filter sensitivity arm.** The secondary `C_raw` and the tercile
medians of SD are reported alongside; they do not change the verdict.

**Nulls and controls (all preregistered).**
* **N1 target-shuffled state pairing** — `cos(r(P, s_i, t), r(P, s_j, σ(t)))` with `σ` a seeded derangement of
  the analysis targets: the S1 construction with the target correspondence destroyed. Expectation ≈ 0.
* **N2 = S12** — the "different state, different donors" cell; `S1 > S12` is what shows the across-state
  agreement is not an artefact of the target set, and `S2 > S12` likewise.
* **N3 depth matching** — every `(P, s, t)` count vector thinned to a per-target budget equal to the **minimum
  total over the 18 (pair, state) cells** (measured cost: 0.441 of the median, so this arm is ~44 % depth) so
  that S1, S2 and S12 are read at *identical* counting depth. Answers "is S1 < S2 just a depth difference?".
* **N4 pairing symmetry** — S1 is reported per donor pair; a state effect must be present in most of the 6.

## 3. The shared vs state-specific decomposition (mandate item 3)

Per target `t`, using a pairing's two complementary donor pairs as **independent replicates** `r₁ = r(P, s, t)`,
`r₂ = r(P̄, s, t)` (averaged over the 3 pairings). All inner products are **cross-replicate**, so they are
unbiased for the reproducible part and carry no noise term:

| quantity | definition |
|---|---|
| total observed energy `T` | `mean_s ½(‖r₁ₛ‖² + ‖r₂ₛ‖²)` |
| reproducible energy `R` | `mean_s ⟨r₁ₛ, r₂ₛ⟩` |
| **shared across states** `Sh` | `⟨m₁, m₂⟩`, `mᵢ = mean_s rᵢₛ` |
| **state-specific & reproducible** `Sp` | `R − Sh` |
| **noise** `Nz` | `T − R` |

Reported as fractions of `T`, pooled (energy-weighted, `Σ_t` component / `Σ_t T`) and per target (median), by
HGNC family and by `C_noCD4` tercile. **Scope note, frozen:** the residual centring removes each state's
across-target mean, so a *global* activation programme shared by every target is already gone; `Sp` is the
**target-specific** part of state dependence, and the number is a lower bound on total state dependence.

**Effective rank of the reproducible state-specific part.** Build `Dᵢ[(t,s), :] = rᵢ,ₛ,ₜ − mᵢ,ₜ`, form the
cross-replicate Gram `G = ½(D₁D₂ᵀ + D₂D₁ᵀ)` in the (target, state) space (≤ 726 × 726; its non-zero spectrum
equals the gene-space one), and report the **participation ratio** `(Σλ₊)² / Σλ₊²` over the positive
eigenvalues, plus the same for the shared part `Mᵢ[t,:] = mᵢ,ₜ` and for each pairwise state contrast
`dᵢ(t) = rᵢ,ₛⱼ,ₜ − rᵢ,ₛᵢ,ₜ`. A participation ratio near 1 means one shared axis; near `n` means
target-specific.

## 4. Known-answer rows (must reproduce before any 282 number is read)

| row | requirement |
|---|---|
| **KA-1 (the ceiling reproduces 280)** | running 280's `q1_280.py` §3.2 code path on 280's own filter and per-condition target sets must reproduce R280 §2.1 **exactly**: Stim48hr pairings +0.0421 / +0.0538 / +0.0599 and POOLED +0.0539 [+0.0427, +0.0620]; Rest POOLED +0.0478; Stim8hr POOLED +0.0397. Bar: max\|Δ\| ≤ 5e-5 on the medians (they are printed to 4 d.p. in R280). |
| **KA-2 (extraction integrity)** | 280's KA-0 on a 4,000-row sample: every bank row's gene sum equals `obs.total_counts`, max\|Δ\| = 0. |
| **KA-3 (conservation join)** | every analysis target carries a `conservation276.csv` row and the re-derived `C_raw` from `cons276.npz`'s pair matrix equals the published column, max\|Δ\| ≤ 1e-6 (regression check on the `C_noCD4` machinery). |
| **KA-4 (null)** | N1's mean across-state cosine within ±0.02 of 0. |
| **KA-5 (design)** | pooled-4-donor cells per (target, state) on the analysis set reconciles with 280 §1.4 / KA-5 (280: 500 / 476 / 481 all-guides, 537 / 536 / 587 measured here on 280's filter over the common set). |

**Abort rule:** if KA-1 does not reproduce, the experiment stops and reports the failure; no 282 number is
quoted from a bank that does not reproduce 280's ceiling.

## 5. What would falsify the headline

* `S1 ≈ S2` (SD interval covering 0) ⇒ *a perturbation response in one person is as reproducible across
  activation states as it is across people* — a real, publishable negative for the organising question.
* `S1 < S2` but also `S1 ≈ S12` and `S2 ≈ S12` ⇒ the ordering is a depth or centring artefact, not state.
* `SD` uncorrelated with `C_noCD4` (S4's interval covering 0) ⇒ the cross-line conservation axis (F8/F9) does
  **not** reduce to within-individual state dependence; the two are different properties.

## 6. Reporting

`RESULTS.md` carries: the KA table; S1 / S2 / S12 / W with intervals and the nulls; the SD table with the
family and tercile contrasts and the S4 pass rule; the decomposition and effective ranks; a **post-hoc,
labelled** descriptive section (the 20 most and 20 least state-dependent targets, and whether they are
immune-programme genes — labelled by HGNC gene group plus a frozen keyword match on the HGNC `name` field:
`interleukin|chemokine|cytokine|interferon|T cell receptor|CD\d|tumor necrosis|JAK|STAT|NFKB|MAP kinase|
cyclin|proliferat`); and the single line
**STATE DEPENDENCE WITHIN INDIVIDUALS: median SD … [lo, hi]; PREDICTED BY CROSS-LINE CONSERVATION: YES / NO
(ρ …)**. Per-target table `work/state282.csv`. No model, no submission, no priced arm.

Seed 20280906 throughout (280's seed). `uv run python`, `PYTHONDONTWRITEBYTECODE=1`, peak RSS reported.
