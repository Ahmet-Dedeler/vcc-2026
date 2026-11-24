# 518 — The bulk mse gap: where our error lives, and the next transferable levers

2026-09-28, research worker (Claude Opus 5.5). Analysis only: no build, no `vcc`, no `cell-eval2` CLI (the H1 legs use
096's `msekit`, which calls cell_eval2's own primitives), no git, nothing written outside this directory. The one
builder patch here (`patch_builder518.py`) is **not applied**; it was checked with `--check` on a scratch copy.

## Answer

1. **At the cliff, raw mse reads the aggregate cosine: raw ≈ 1.032 − 2X/sumD, and 2X/sumD = 0.406 for 446bj2.** The
   re-level pins each context's spread energy B at the claim (the truth's own jackknife), so the credit cancels our
   energy, and only X (energy-weighted direction, including how the energy is split across targets) moves the member.
   **+1 % relative X = −0.0041 raw mse = +0.0007 overall.** At our constants the rivals' raw mse implies +13.8 %
   (cqawesome), +13.3 % (Testh) and +33.6 % (Illumina) more aggregate cosine than ours.
2. **The magnitude pot, at our directions and pinned energy:** a per-target magnitude oracle raises the aggregate cosine
   by **+25.1 %** on the cancer PLOLO folds (+23.2 / +25.1 / +27.1 %), **+52 %** on kolf2 and **+26.5 %** on H1
   (noise-corrected, exact legs).
   - Board value at k 1 (1.31): **−0.10 (−0.13) raw mse, +0.018 (+0.023) overall.**
   - It is pds-free by construction, and on its own it is larger than the whole mse gap to cqawesome or Testh. That fits
     Testh's profile: worse pds, much better mse, and a `learnedgate_tgtw` model string whose step moved its mse
     0.810 → 0.698 (read as per-target weights, FINDINGS §1).
3. **The direction pot is far larger but information-limited.** On H1 the direction oracle at our per-target energies
   takes the aggregate cosine from 0.113 to 0.818 (×7.3), against ×1.27 for the magnitude oracle. The mean-response
   component is small: the truth mean is 1.1–1.7 % of A/B/C's signal (read off the metric's own baseline anchors).
   Sub-gate genes (< 10 CPM) carry 11 % of H1's noise-corrected signal. We predict nothing there.
4. **One new transferable lever, and it overlaps the cap in flight.** A **strength key** scales each target's bulk row by
   a robust response-strength rank. The rank is a noise-normalised count of strongly responding genes in the donors' own
   rows, plus CD4 T cells, VIPerturb-K562 and the KOLF2 atlas.
   - On 446bj2's key it prices **+0.0016 … +0.0024 overall** (exact PLOLO +2.40 %, 3/3 folds; kolf2 +5.6 %; H1 −0.0126 raw,
     beating 20/20 permutations).
   - That is the same pot and the same size as 495's cap 6 (+0.0019 … +0.0029). **On top of cap 6 it adds only +0.0001 …
     +0.0003.** Two independent estimators (the LODO key's tail and the donors' own response strength) find the same
     thing: strong targets are under-scaled.
5. **What remains is not reachable with any per-target information we hold.**
   - After cap 6 plus the strength key, the magnitude oracle still offers ≈ +20–23 % on the cancer folds.
   - Half of it sits in 11–19 line-idiosyncratic targets per held-out line.
   - Nothing predicts them out of fold: donor features, CD4 alone, context expression, or a correction learned across
     folds.
6. **Instrument finding: PLOLO mis-signs render and gate changes on low-expressed genes.** It renders in the challenge
   frames but scores against the held line's own frame. Lowering the 10-CPM render gate reads +1.3 % on PLOLO (3/3).
   When prediction and truth share a frame it loses −0.9 % (7 CPM) … −8 % (0 CPM) on the Orion folds, and is net
   negative on H1 (gate 2: mse −0.005 raw but pds −0.0045; gate 0: both worse).

**Verdict:** no new lever clears +0.002 on two surfaces *on top of 495*. The strength key clears it against 446bj2
(cancer at k_X 1.31, kolf2, H1). Its build spec (§6) is a **contingency** for a failed or partial 495 read, and a second,
bounded estimator for October. It is not a slot to stack on cap 6.

## 1. Known answers (every script aborts on failure)

| row | expected | measured |
|---|---|---|
| KA1 187's cliff identity at 446bj2's measured energies, X solved from the board | raw 0.6264081539046438 | \|d\| 0 (`ka518.py`) |
| KA2 540's recorded cap-6 price, re-derived from the saved inputs (analyse540's arithmetic) | r 0.026192735825910935, dmse k1.31 −0.01395937378276535, net_K +0.00288936059067574 | \|d\| 0 |
| KA3 420's saved per-target arrays → comp420's realised-16 cos_agg, 446bj2's bulk unlevelled and levelled, 4 folds | comp420 / comp420L | \|d\| ≤ 3e-16 |
| KA4 H1: msekit on 435's harness stats + 323's reference | raw mse 0.9253584049461503 | \|d\| 0 |
| KA-dump the fold dumps (unkeyed PP/stack compositions + LODO key vector), rows rebuilt without the operator | arms540's V1_0.5_3 and V1_0.5_6 cos_agg / pds / pds_lam0.5 | \|d\| 0, 4/4 folds (`dump518.py`, and on every `score518.Fold()` load) |
| KA-H1 legs: my N_p, D_p vs msekit per target; identity terms rebuild the raw | 0.9253584049461 | rel 1.4e-11 / 2.2e-11; exact |
| KA-H1 jk: per-gene delete-1 jackknife streamed from the cells, summed | the reference jk per target and control | rel 2.7e-16 |
| KA-H1 pds: my cosine discrimination on the realised rows | the harness's 0.8952380952 | 0.89524 |
| KA-abc: cap 6's A/B/C energy ratios and dK in `abc518.py` | abc540.json `V1_0.5_6` | \|d\| 0 |
| KA-analyse: skey518's exact cap-6 r per cancer fold | arms540's recorded r | \|d\| < 1e-12 |
| KA-art: the strength artifacts at θ 0 | 8c82bef3's w (cap 3) and 495's dc374281 w (cap 6) | bit-identical |

## 2. Step 1 — the decomposition

### 2.1 The board's own identity for 446bj2 (`board518.py`; per context, units of sumD = 0.2291)

| ctx | A/sumD (our energy) | mean row (A−B)/sumD | Kc/sumD | credit/sumD | ρ = B/claim | 2X/sumD | cos_agg = X/√(A·sumD) | X/A | raw |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 1.019 | 0.0079 | 0.022 | 1.011 | 0.991 | 0.406 | 0.201 | 0.199 | 0.6243 |
| B | 1.037 | 0.0086 | 0.023 | 1.028 | 0.990 | 0.406 | 0.199 | 0.196 | 0.6256 |
| C | 1.136 | 0.0094 | 0.026 | 1.127 | 0.987 | 0.406 | 0.190 | 0.179 | 0.6294 |

- **Our energy equals the truth's signal energy (A ≈ B ≈ claim ≈ sumD), and the credit pays for all of it except the
  mean row.** The energy-optimal global scale without the credit would be X/A ≈ 0.19 of today's. The cliff, not the
  signal, sets our level.
- Below the cliff, raw = 1 + [(A−B) + Kc − 2X]/sumD. So: spread energy is free, the mean row costs, Kc (the observed
  control's noise, 2.2–2.6 %) is a floor, and **the member is 1.032 − 2X/sumD**.
- **Rivals at our constants** (2X/sumD they need; ratio to ours): cqawesome 0.462 (×1.138), Testh 0.460 (×1.133),
  Illumina 0.542 (×1.336). The mse-only overall gap is +0.0098 / +0.0095 / +0.0239.
  - Caveat: this assumes their mean-row and Kc terms equal ours. Their claims may differ; 227 showed a 2.4 % claim
    difference reads as ~5 % cosine.
- **Truth mean response** from the metric's baseline (every row = the pooled perturbed profile, so B = 0 and there is no
  credit): 1 − b = 0.8–1.4 % of sumD, or 1.1–1.7 % with the jackknife's +0.32 %. Our mean row is 0.8–0.9 %.

### 2.2 The pots per surface (relative X at pinned energy; board raw = −0.406 × k × r)

| surface | aggregate cos now | magnitude oracle r | of which: zero the anti-aligned targets | direction oracle (at our energies) | noise ceiling on the oracle |
|---|---:|---:|---:|---:|---:|
| PLOLO k562_gwps (446bj2's bulk, levelled) | 0.0643 | **+23.2 %** | +1.3 % (11 % of targets m < 0) | n/a (truth noise) | ≤ 0.9 % of Σm² |
| PLOLO orion_hct116 | 0.0526 | **+25.1 %** | +3.4 % (15 %) | n/a | ≤ 1.3 % |
| PLOLO orion_hek293t | 0.0407 | **+27.1 %** | +4.1 % (15 %) | n/a | ≤ 2.1 % |
| PLOLO kolf2 | 0.0285 | **+52.2 %** | +12.4 % (31 %) | n/a | ≤ 3.0 % |
| H1 (435 harness bulk, exact legs, jackknife-corrected) | 0.113 | **+26.5 %** (raw +26.6 %) | 13 % of targets m < 0 | **0.818 (×7.3)** | 0.18 % measured |
| board (mapped; k 1 / k 1.31) | 0.20 | **−0.102 / −0.134 raw** | | the rest of 0.626 → ~0.03 | |

- Fold numbers come from 420's saved per-target arrays (`folds518.py`), the direction oracle from `h1_518.py`, and the
  post-cap-6 pot below from the fold dumps (`remain518.py`).
- **Which pot is larger?** Direction, by far, wherever it can be measured: H1's energy-weighted RMS per-target cosine is
  0.143.
  - On the PLOLO folds the direction pot is unmeasurable: the donor-pseudobulk truth is noise-dominated (cos_agg
    0.03–0.06), so only the magnitude pot (unbiased: its noise ceiling is ≤ 3 %) can be read there.
  - The magnitude pot is the pds-free one. It is the one a per-target-weight model (Testh) can win without better
    ranking.
  - It is ~2× the mse gap to cqawesome or Testh.
- **How the pots map onto the members.**
  - Magnitude reallocation is exactly pds-invariant (pds is per-target scale-free) and moves mse only.
  - Direction moves both. pds weighs every target equally; mse weighs each by its truth energy. On H1 the top decile of
    targets carries 51 % of the signal and 62 % of X.

### 2.3 Where the magnitude pot lives

- **Over-scaled targets carry 57–60 % of the allocation loss** on every fold. The attribution is additive:
  X* − X = Σ_p (s_p − λm⁺_p)²/2λ + s_p·max(−m_p, 0).
- **By our energy decile (folds):**
  - The top two deciles hold 71–75 % of our energy and 70–87 % of X.
  - They also hold 38–59 % of the loss.
  - Their median ideal/current scale is 0.72–0.80, i.e. over-scaled.
  - The ideal scale has a heavy upper tail: 10/50/90 % quantiles 0.00 / 0.62–0.71 / 1.7–2.2.
- **By TRUE energy (H1, noise-corrected):**
  - The top decile holds 51 % of the truth signal, 62 % of X, and only 34 % of our energy.
  - Its median ideal/current is **1.31, under-scaled**.
  - Spearman(ideal scale, true strength) = +0.69; Spearman(our energy, true strength) = +0.50.
- Both views say the same thing: **our per-target energies track true strength only partially and with a compressed
  range.** Truly strong targets are under-scaled; some targets we predict strong are not.
- **After cap 6 (495) and the strength key the pot is barely dented** (`remain518.py`):
  - cancer +22.6 / +24.4 / +26.4 % on 446bj2 → +23.0 / +20.8 / +23.6 % on cap 6 → +22.9 / +20.4 / +23.3 % adding the
    strength key; kolf2 +52 → +41 → +38 %.
  - **Half the remaining loss sits in 11–19 targets of ~270–294 per held-out line** (5–6 on kolf2).
- **By context:** the pot is uniform across the A/B/C frames (±2 points on every fold). On the board the three contexts
  differ only in energy: C has the most (A/sumD 1.14) and the lowest cos_agg (0.190).

### 2.4 Gene classes (H1 exact, noise-corrected; `h1_518.py`, 85 targets)

| class | genes | truth signal share | our energy share | X share | class cosine |
|---|---:|---:|---:|---:|---:|
| gated (> 10 CPM in the control) | 9,598 | 89.0 % | 99.6 % | 100.7 % | 0.121 |
| **ungated (< 10 CPM): the builder's bulk is zero there** | 8,482 | **11.0 %** | 0.4 % | −0.7 % | −0.04 |
| other H1 targets (analogue of the panel columns the mse scores) | 85 | 1.1 % | 1.1 % | 1.1 % | 0.115 |
| the 300 A/B/C panel genes (ordinary genes in H1) | 295 | 2.2 % | 3.0 % | 2.7 % | 0.120 |
| gated expression deciles, lowest → highest | 960 each | 6.2 → 11.4 % | 8.5 → 12.7 % | 4.3 → 17.8 % | **0.067 → 0.166** |

- **Sub-gate signal is target-specific**, not a shared component: the across-target mean is 2–4 % of each sub-gate class.
  Its signal-to-noise is 1.1–1.2 (≥ 10 CPM: 2.5). So it is a direction problem.
- **Mean component on H1:** the truth mean is 4.0 % of the signal (the A/B/C baseline says 1.1–1.7 %). Ours is 1.6 %, and
  **cos(our mean, truth mean) = 0.06**. The cancer folds against an uncentred truth give 0.09–0.11 (`meanfold518.py`).

## 3. Step 2 — levers (all one variable vs 446bj2's bulk; board units = overall unless marked raw)

| lever | class (why) | surfaces | price | verdict |
|---|---|---|---|---|
| **Strength key S1**: w × exp(1.25 z), z = robust strength rank over six screens | key: a per-target magnitude multiplier, same mechanism as the LODO key and the cap | exact PLOLO +2.40 % (3/3; bootstrap CIs [−0.5, +4.5] / [+0.03, +5.2] / [+0.3, +5.1] %), permuted control −0.9 / −1.5 / −1.2 %; kolf2 +5.6 % (control −3.8 %); H1 −0.0126 raw vs null +0.0031 ± 0.0074 (0/20 permutations as good) | **+0.0016 … +0.0024 vs 446bj2** (k 1 … 1.31, both ballast models); kolf2 gross +0.0038; A/B/C energy ×0.98, max row ×1.91, ballast K −0.61 | **contingency** (§6) |
| Strength key on cap 6 (S2, θ 0.5) | key | PLOLO +0.23 % over cap 6 (3/3) | **+0.0001 … +0.0003 over 495** | not a slot: same pot as the cap |
| src-only strength (the 3 production donors, October-portable) | key | exact +1.94 % (3/3), kolf2 +3.1 % | +0.0013 … +0.0020 vs 446bj2 | the October form, if ever |
| CD4 T-cell DE count alone (Marson GWCD4i DESeq2) | key | LOFO +0.92 % (3/3), kolf2 +5.3 % (analytic) | ≈ +0.0007 | inside the strength key |
| Donor-side features: knockdown depth, crude DE count, cross-donor cosine, norm dispersion, target expression in the context or held line, panel share | key | Spearman with the ideal scale ≤ 0.2, mixed signs; LOFO −0.6 … −1.4 % | ≤ 0 | **negative** |
| Correction learned across folds (other folds' residual ideal scale) | key | Orion +2.0 / +2.3 % but k562 −1.2 %, kolf2 −1.2 % (clean test); on cap 6 kolf2 −1.7 % | ≤ 0 | **negative**: the cross-fold consistency is the two Orion lines' shared lab |
| Across-target strength tilts, top-10 % ×1.25 / ×1.5, bottom-50 % shrink | key | folds mixed (cancer −3.2 … +0.2 %, kolf2 up); H1 tilt optimum at 0 | ≤ +0.0001 | negative (540's headroom result, again) |
| Remove our mean row (centring) | composition (the common component) | H1 −0.0012 raw, pds +0.003; folds (uncentred truth) −0.17 / +0.07 / −0.91 %, pds +0.008 / +0.002 / −0.005 | mse ≤ +0.001 (393's best case) ± pds ~0.0006 | small; 393 stands |
| Add a predicted mean response | composition | pds-toxic (the fold oracle mean: pds −0.14 … −0.24) | < 0 | closed (§2 "coherence is the price") |
| Lower the render gate (additive render) | shape | PLOLO (challenge frame) gate 2: +1.29 % (3/3); **same frame: Orion −3.3 / −2.7 %, H1 −0.005 raw with pds −0.004** | ≤ 0 | **negative, and a PLOLO artefact** (frame mismatch) |
| Hybrid render: production ≥ 10 CPM, fold-change reading of h below | shape | same-frame Orion +0.41 / +0.40 % (pds +0.0004 / +0.0002); H1 −0.0075 raw (+4.0 % X), pds −0.0013 | +0.0001 (Orion-based, k 0.4–0.65) … +0.0026 (H1-based, k 1) | unresolved; below the two-surface bar |
| Predict the other panel genes' columns (the mse scores them, pds drops them) | new coordinates | H1: zeroing the 295 A/B/C-panel columns costs +0.0018 raw (≈ 1 % X) | ≈ +0.0007, pds-free; needs a builder change | small |
| Expression tilt (cpm/median)^0.1 after sharpening | shape | PLOLO +0.30 % (3/3), kolf2 +1.25 % | ≈ +0.0002 | negligible |
| q 2 | shape | PLOLO +0.84 % (3/3), kolf2 −0.84 %, H1 +0.002 raw worse | — | the coordinator's q stops at 1.75 |
| ρ* (the re-level target 0.9793) | level | the EV optimum under 187's cliff-position sd is +0.2 … +0.5 % in scale above 446bj2's factors; EV raw 0.6252 vs 0.6264 | ≤ +0.0002 | closed |

**Levers in the direction pot proper (per-target cosines).** They are outside this worker's reach by construction: the
per-target direction is the composition and operator, and §2's closed list covers donor transforms, operator geometry,
a fifth source, neighbourhoods and consensus signs through the bulk. What this worker could touch are the mse-specific
coordinate classes:
- sub-gate genes (hybrid render, unresolved);
- the panel-gene columns (≈ +0.0007);
- the expression tilt (≈ +0.0002);
- the mean component (small).
None clears the bar.

**The ρ* question.** It cannot be read on the folds. PLOLO has no jackknife and no claim, so there is no cliff, and the
fold-optimal pure-energy scale (X/A) is meaningless when the credit sets the level.
- On H1, with exact claims, the optimum sits exactly on the cliff (s 0.930 vs the cliff at 0.929). The 0.9793 hedge is
  the price of the NTC-proxy claim's uncertainty on A/B/C.
- 0.9793 on the NTC proxy corresponds to 187-model ρ 0.987–0.991.

### 3.1 The strength key, in detail (`strength518.py`, `skey518.py`, `abc518.py`, `analyse518.py`, `robust518.py`, `h1skey518.py`)

**The statistic, per screen.** Only each screen's panel rows are read.
- Genes: those > 10 CPM in the screen's control, minus the panel genes.
- Per-gene noise scale: 1.4826 · MAD across the panel targets (most are weak).
- z = delta / scale; per-target noise = median |z| / 0.6745 (cells per perturbation vary 20–1,963).
- **strength = log1p #{g : |z| > 4 × noise}.**
- On CD4 it agrees with DESeq2's nDE at Spearman +0.60. A crude |delta| > 0.25 count, which ignores per-target noise,
  loses (LOFO −1.4 %).
- The key feature is the centred percentile rank averaged over the screens that carry the target.

**Out of fold** (one θ fitted on two cancer folds, scored on the third; kolf2 scored with θ from the three cancer folds):

| feature | held-out cancer (k562 / hct116 / hek293t) | kolf2 |
|---|---|---:|
| sources' robust strength | **+2.44 %** (+2.93 / +2.21 / +2.19) | +3.0 % |
| meta (sources + CD4 + VIPerturb + KOLF2) | **+2.51 %** (+1.88 / +2.87 / +2.79) | +6.1 % |
| CD4 robust / CD4 DESeq2 / VIPerturb / KOLF2 alone | +1.93 / +0.92 / +0.81 / −0.05 % | +3.8 / +5.3 / +2.0 / 0 % |

**Caveats.**
- The count threshold matters: |z| > 3 gives +0.99 %, 5 gives +3.03 % but kolf2 −1.0 %. Four was fixed before looking
  and is not tuned.
- On the k562 fold the bootstrap interval touches zero (P ≤ 0 = 0.06).
- **Mechanism:** our per-target energy comes from donor row norms, and those norms include each donor's noise (few-cell
  perturbations have large noise norms). The count separates real response breadth from noise, so it moves energy from
  noise-dominated weak rows to broad, strong responses. That is the same move the cap makes with the LODO key's upper
  tail.
- **Why it adds nothing to cap 6:** both boost the same strong, consistent targets. The combination is +2.86 % against
  +2.62 % for cap 6 alone.
- **Substitute or complement?** (`subst518.py`, exact) The LODO key itself is worth +18.1 % over no key (cancer mean).
  The strength key alone recovers only 36–61 % of that (θ 1 … 2: −11.5 … −7.1 % vs 446bj2). It overlaps the LODO key
  heavily and cannot replace it, even in October. It is a margin on top of the key (+2.4 %) and nearly nothing on top of
  cap 6 (+0.2 %).

## 4. Instrument notes for the next worker

- **PLOLO cannot price render or gate changes on low-expressed genes.** It renders in A/B/C's frames against a truth in
  the held line's frame. Genes low in A/B/C but expressed in the held line make lower gates look good (+1.3 %).
  - In a shared frame they lose: Orion −0.9 % (7 CPM) … −8 % (0 CPM); H1 net negative.
  - `h1arms518`'s variant, which also sharpens on the lowered gate, reads gate 2 at −0.019 raw with pds −0.009. The pds
    cost is again as large as the mse gain.
  - Use same-frame surfaces: the Orion folds rendered in the held control (`render518.py`), or H1. The K562-GWPS, KOLF2
    and Nadig banks were pre-filtered to expressed genes and cannot be used.
- **PLOLO's truth is bank-centred** (E062.raw_truth subtracts the bank's mean over all its perturbations), so it is blind
  to the NTC-referenced mean response. Score the mean component against the raw delta rows (`meanfold518.py`).
- **H1 has too few targets (85) to resolve per-target magnitude levers alone.** A random per-target reweighting moves raw
  mse by ±0.007–0.009 (sd). Always report the permutation null (`h1skey518.null`).
- **The fold dumps** (`work/dump518-<fold>.npz`, 84–102 MB each, local only) hold the unkeyed PP and stack compositions, the LODO key
  vector, truth and betas. `score518.Fold` rebuilds any per-target-key, gate, gene-multiplier or centring arm exactly, at
  ~1–2 GB with no operator, in ~1 s per arm. It re-asserts arms540's rows on every load.

## 5. Negatives, with mechanisms (so nobody re-runs them)

- **The remaining ~20–23 % magnitude pot after cap 6 is line-idiosyncratic.** It is concentrated in ~15 targets per line.
  No donor-side, CD4/VIPerturb/KOLF2, context-expression or cross-fold feature predicts it out of fold. The Orion–Orion
  cross-fold agreement (Spearman +0.39) is lab, not biology: it fails on K562 and kolf2.
- **Per-target tilts by our own energy, pre- or post-key, and top/bottom-decile rescaling** re-derive 540's headroom
  result. Our energy carries no residual information about the ideal scale (H1 Spearman +0.12; folds mixed).
- **Lowering the gate** is a PLOLO artefact (§4).
- **The mean component:** ours is nearly orthogonal to the truth's (0.06–0.11). Removing it buys ≤ 0.009 raw with an
  unsigned pds effect. Adding one is pds-toxic.

## 6. Build spec (CONTINGENCY — for a failed or partial 495 read, or as October's bounded second estimator)

- **Arm S1 = 446bj2 with its stack dumped under `VCC_PERT_SCALE=strength518_S1_lo0.5_hi3`**, re-levelled by ntclevel to
  ρ 0.9793.
  - Artifact: `work/pert_scale_518_S1_strengthmeta_th1.25_lo0.5_hi3.npz`, sha256
    `1443c1a41b285872769f20a99ee958d06266b0171e99328de22f13a4cdd00da6`.
  - Rule: w = clip(key/median(key[carried]), 0.5, 3) × exp(1.25 z), where key is 8c82bef3's own vector and z is from
    `work/strength518.npz`.
  - w ranges 0.30 … 5.22 (multiplier 0.61 … 1.74); A/B/C energy ×0.984 / 0.985 / 0.976; max per-row ratio ×1.91; ballast
    mean K 8.83 → 8.21.
- **Registration:** `patch_builder518.py <builder> [--check]` adds two `PERT_SCALE_ARTIFACTS` entries, a pure addition
  anchored after 495's `lodo415_lo0.5_hi6`. It applies cleanly (checked on a scratch copy, 2026-09-28).
- **Expected config diff vs 446bj2:** exactly `bulk_load`, `bulk_load_sha256`, `bulk_ctx_scale` (as 495's), plus the dump
  env's `VCC_PERT_SCALE` inside the dump record.
- **Predicted board row vs 446bj2** (540's calibration: literal pds × 0.34, mse at k_X 1.31 … 1, ballast K / D models):

  | member | predicted |
  |---|---|
  | pds (raw) | +0.0002 |
  | mse (raw) | −0.0128 … −0.0098 |
  | nmae (raw) | −0.0005 … +0.0007 |
  | fid / reach / jac | pinned |
  | **overall** | **+0.0016 … +0.0024** |

- **Bins:**
  - transfers: overall ≥ +0.0010 and Δmse ≤ −0.005;
  - partial: +0.0003 … +0.0010;
  - fails: < +0.0003 or Δmse > −0.002;
  - DE leak check as 495's.
- **When to use it:**
  - If 495 transfers, do not build S1. S2 (cap 6 + θ 0.5, sha `dc0856e0…`) is worth only +0.0001 … +0.0003 and is not
    a slot.
  - If 495 fails, S1 tests the same mechanism with bounded tails (no 6× weights). A second failure closes the pot's
    reachable part.
- **For October:** the "src" form needs only the three production donors' panel rows, which are always on disk.
  - It prices +0.0013 … +0.0020 on A/B/C folds.
  - Its value there is robustness against the LODO key's heavy October tails (540o: uncapped weights 19–36). That is
    untested on 540o's October-like panels. It is the next experiment if the coordinator wants it.

## 7. Memory and deviations

- `dump518.py` ran once per fold at **5.4–5.5 GiB peak**, one at a time. Before each run free+inactive was 33.0–33.5 GiB.
  That is the mandate's PLOLO exception (≥ 33 GiB, one at a time) and above the 4 GB general cap. It was a one-time
  cost; everything later reuses the dumps at 0.4–2.7 GB.
- `h1_518.py`'s first run streamed the per-gene jackknife at **4.25 GB peak**, slightly above the 4 GB cap, without a
  fresh vm_stat check. It is cached now (`work/h1_jk_pergene.npz`); reruns stay under 1 GB.
- Every other script peaked at ≤ 3.1 GB (`/usr/bin/time -l`).

## 8. Files

- **Known answers:** `ka518.py`. **Board identity:** `board518.py`.
- **Decomposition:**
  - `folds518.py` (the pots from 420's arrays);
  - `h1_518.py` (the exact H1 identity, gene classes, mean component);
  - `remain518.py` (the pot after cap 6 and the strength key);
  - `meanfold518.py` (the mean component vs an uncentred truth).
- **Machinery:**
  - `dump518.py` (per-fold dumps) and `rows518.py`;
  - `score518.py` (exact arms, no operator);
  - `screen518.py` (the one-variable screen: gate, strength tilts, centring, expression tilt, q).
- **Features and the strength key:**
  - `feat518.py` (donor features), `cd4feat518.py`, `lofo518.py`, `strength518.py`;
  - `skey518.py` (exact), `robust518.py`, `xfold518.py`, `subst518.py` (substitute vs complement);
  - `abc518.py` (A/B/C), `analyse518.py` (pricer), `h1skey518.py` (H1 with its null).
- **H1 levers:** `h1arms518.py` (the emulated recipe), `render518.py` (same-frame render test).
- **Build spec:** `mkart518.py` (artifacts plus their known answers), `patch_builder518.py` (not applied).
- **Outputs:** `work/*.json` and `work/*.log`.
  - `*.npz` is gitignored. The small arrays the build spec needs must be committed with `git add -f`:
    `work/pert_scale_518_S1_*.npz` and `work/pert_scale_518_S2_*.npz` (25 KB each, with their `.json`),
    `work/strength518.npz` (10 KB), `work/feat518_arrays.npz` (0.3 MB).
  - Local only, never to be added (> 5 MB): the four fold dumps (84–102 MB each), `h1_518_rows.npz` (29 MB),
    `h1arms518_recipe.npz` (21 MB), `h1_jk_pergene.npz` (10 MB).
