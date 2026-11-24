# 265 — The proper test of 263's `kolf2` majority result. Per-target additivity FAILS (the pipeline is not row-local: one global preserve scalar in `Operator.apply` feeds a nonlinear render), but the failure cancels on arm TOTALS to 2.6e-4 on `kolf2` (0.04 of a null SD), so the permutation test runs on a labelled additive proxy. **On `kolf2` the majority's excess over frequency-preserving random drops survives: p = 0.011 (proxy; 0.012 with the real total; z +2.3 against a null SD of 0.0068).** Its excess over the INCUMBENT does not clear the paired target bootstrap (+0.01209, 95 % [−0.0035, +0.0278], P(≤ 0) 0.064), the τ = 0 gate sits on the line (p 0.048 / 0.051), and 264's "random drop = 56 % of the arm" rested on three seeds that are an upper-tail draw of the null (joint p 0.005). **On the three cancer folds, in 264's matched construction, the same test reads chance: frequency-preserving p 0.68 / 0.12 / 0.50** — the rule beats a uniform random drop there only through its donor frequencies (always-drop-`kolf2` beats the frequency-preserving null at p 0.008 / 0.018 / 0.013; the rule itself does not). (b) − (c) is significantly negative (3-fold −0.0107 [−0.0187, −0.0030]). 263's two positive reverse cells and its best arbiter cell all include 0. **Measurement only. Recommendation in §6: the re-fitted three-donor validation is NOT worth building as a priority; H2's majority variant does not re-open on these numbers.**

**Date 2026-09-05 (12:38–13:00 UTC).** Deterministic (seed 20260905). **No build, no `vcc`, no `cell-eval2`, no staging, no fold rebuilt, no operator loaded; nothing under `submissions/`, nothing written outside this directory** (`find -newer work/.marker-start -type f` over `experiments/ data/ scripts/ src/ docs/ submissions/` returns only sibling agent 266's own work files and bytecode, 256's `backfill-content.log` and the coordinator's ledger; `git status` on 263/264 shows only the coordinator's earlier `analyse264` re-fix; `PYTHONDONTWRITEBYTECODE=1` on every run, no `__pycache__` written by this experiment). Arithmetic on saved arrays only: **22.3 s, peak RSS 0.40 GiB.** HEAD `924c681`. `uv run python` throughout. numpy 2.5.2.

| row | required | measured | verdict |
|---|---|---|---|
| **KA-0** 264's `kolf2` arrays are 263's | per-target `pds` of `cd4svd5`, `majority:raw:all:nm`, `majority:raw:tau0:nm`, the reverse, and the majority's drop vector bit-identical; Δ`pds` +0.012089 | **0.0e+00** on all five; Δ +0.012089 (|d| 1.9e-7) | **PASS** |
| **KA-1** every saved per-target array averages to its sweep JSON total (25 + 29 × 6 rows) | 1e−12 | **≤ 2.2e-16** | **PASS** |
| **H-2** (the mandate's premise, MEASURED) the always-drop-j arms leave the targets they do not touch at the base value (row-locality of the pipeline) | 0.0 | **FAILS on 7/7 folds**: `kolf2` max |Δ`pds`| 0.0123 on the 33 untouched targets (own cosine moves ≤ 9.5e-4); matched cancer folds 0.055–0.101; letter 0.48–0.84 (§2) | **FALSE** |
| **H-3** (the mandate's step 1) an arm's per-target `pds` equals its per-target choices read off the always-drop-j arms; its total equals the sum of choices at 0.0e+00 | 0.0 | **FAILS per target on 7/7 folds** (`kolf2` `:nm` arms: 98–143 of 273 targets off, max 0.017); **on the TOTALS the deviations cancel: ≤ 2.6e-4 over the 13 `kolf2` `:nm` arms** (0.04 null SD), ≤ 1.6e-3 on the matched cancer folds (0.18–0.27 null SD), ≤ 2.1e-2 on the letter folds (1.4–3.2 null SD) | **FALSE; proxy usable on `kolf2` and the matched folds, NOT on the letter folds** |
| **KA-4** bootstrap point estimates == the sweeps' deltas | 1e−12 | majority − incumbent 4.7e-17; mean(random drop) − incumbent == 264's +0.0068 (2.8e-7); 3-fold (b) − (a) / (b) − (c) == R264 −0.00509 / −0.01073 (matched) and +0.03628 / −0.01945 (letter) at ≤ 4.2e-6; 263's `union:wht:tau0:nm` −0.0004 and `oracle_wht:all:nm` +0.0459 at 1.8e-7 | **PASS** |
| **KA-5** the null's construction: E[uniform null] == mean of the three always-drop-j totals with the untouched-row coupling removed; the MC means (100,000 draws) match the closed forms within 4 MC-se | exact / 4 se | **0.0e+00 (≤ 1e-12)**; uniform 4.1e-5, frequency-preserving 1.8e-5 vs se 2.2e-5 | **PASS** |
| **KA-6** the majority's drop counts == 263 `kolfctrl` 110 / 65 / 65; the oracle drops are extreme in every null | exact; p ≈ 1/B | **110 / 65 / 65**; oracle z +8.3 / +12.2 (`kolf2`), +9.9 … +10.6 (matched folds), p_upper 1e-5 everywhere | **PASS** |
| **KA-7** the outside re-analysis (GPT, 2026-09-05) reproduced | within MC error | majority − incumbent +0.01209 [−0.0035, +0.0278] vs GPT +0.01209 [−0.0034, +0.0280]; majority − mean(random drop) +0.00526 [−0.0115, +0.0222] vs +0.0053 [−0.0115, +0.0228] | **PASS** |
| **KA-8** the null is calibrated where it can be checked: 18 real random-drop seeds on the matched cancer folds positioned in their uniform nulls | p ~ U(0, 1) | mean p 0.484, min 0.022, max 0.990, one ≤ 0.05, one ≥ 0.95 | **PASS** |

Regression checks **25/25**; measured hypotheses H-2 / H-3 **0/14** hold (that is the finding, not a failure of the run).

---

## 0. The answer in ten lines

1. **The mandate's premise is false for the pipeline, true for the scorer.** A target's `pds` is a function of its own RENDERED row against the fixed truths (`resid218._per_target`: the fraction of other truths whose cosine is below its own, ties at half, averaged over frames A/B/C), and 262 KA-3 holds. But `plolo.Operator.apply` rescales every whitened row by ONE scalar (the un-whitened median row norm over the un-preserved whitened median, `plolo.py:754-771`), `E062.add_post_whiten_residual` by another (`run.py:581-588`), and `plolo.deliver` is nonlinear in the row amplitude — so changing 240 rows moves the un-preserved median, hence every row's amplitude, hence every row's cosine matrix row. Measured: the always-drop arms move the `pds` of the 33 targets they never touch by up to 0.0123 on `kolf2` (H-2), and no drop-assignment arm equals its per-target choices (H-3). **Strict additivity fails.**
2. **The failure is a global scalar, and it cancels on totals.** Over the 13 norm-matched `kolf2` arms, |real total − sum of per-target choices| ≤ 2.6e-4 = 0.04 of a null SD; on 264's matched cancer folds ≤ 1.6e-3 (0.18–0.27 SD); on the letter folds ≤ 2.1e-2 (1.4–3.2 SD, because `kolf2` at weight 1 dominates the whitened median). The permutation test is therefore run on the **additive proxy** (an assignment's total = the sum of its per-target `pds` read off the always-drop-j arms), labelled as such, with the p-value also given for the real total, on `kolf2` and the matched folds; on the letter folds only the bootstrap is read.
3. **`kolf2`, bootstrap (exact, real arms, 20,000 resamples of the 273 targets):** majority − incumbent **+0.01209 [−0.00346, +0.02783], P(≤ 0) 0.064**; majority − mean(3 random drops) +0.00526 [−0.01146, +0.02217], P 0.27; majority − reverse +0.0218 [−0.0029, +0.0463], P 0.041; majority − always-drop K562 / HCT116 / HEK293T +0.0152 / +0.0198 / +0.0119 (P 0.038 / 0.037 / 0.13); τ 0: +0.00567 [−0.0067, +0.0180], P 0.19; the oracle drop +0.054 [+0.041, +0.068] (the positive control). Table T1.
4. **`kolf2`, permutation (100,000 assignments per null):** uniform null mean −0.00350 vs the incumbent, SD **0.00696**, 95 % band [−0.0171, +0.0101]; frequency-preserving (110 / 65 / 65 kept) mean −0.00347, SD 0.00685, band [−0.0168, +0.0100]. The majority (proxy +0.01223, real +0.01209): **p_upper 0.0116 / 0.0111 (real total 0.0122 / 0.0117), z +2.26 / +2.29.** τ 0 (129 targets): p 0.049 / 0.048 (real 0.051 / 0.049), z +1.7. The reverse sits at p 0.82 (z −0.95), the whitened-space majority at 0.48, the three always-drop arms at 0.30–0.72. Table T2.
5. **So the "8 σ" was 2.3 σ**, against the proper sampling scale of a random drop-one control on this fold (0.0068), and **264's random control was an upper-tail draw**: the three real seeds sit at p 0.015 / 0.285 / 0.038 in the null they are drawn from (τ 0: 0.11 / 0.51 / 0.031); P(mean of three null draws ≥ the seeds' mean) = **0.0047** (τ 0: 0.038). Their drop vectors are unremarkable (counts 82 / 74 / 84 …, agreement with the majority 0.28–0.35 and with the oracle 0.34–0.36 against chance 0.333 ± 0.030), and the same null is calibrated on the 18 cancer-fold seeds (mean p 0.48). Table T3. "56 % of the arm" measured three seeds, not the control.
6. **Cancer folds, 264's matched construction (the proxy readable, error 0.2–0.3 SD):** frequency-preserving p **0.68 / 0.12 / 0.50** (real total 0.68 / 0.10 / 0.58; z −0.46 / +1.19 / 0.00); uniform p 0.21 / 0.034 / 0.16. The uniform null is beaten on HCT116 only through the donor frequencies: always-drop-`kolf2` beats the FREQUENCY-PRESERVING null at p 0.008 / 0.018 / 0.013 (z +2.4 / +2.1 / +2.2) on every fold, the majority never. Bootstrap: (b) − (a) 3-fold **−0.00509 [−0.01222, +0.00201]**, (b) − (c) **−0.01073 [−0.01871, −0.00302]** (excludes 0), (b) − always-drop-`kolf2` −0.0104 [−0.0178, −0.0031] (excludes 0), (b) − mean(random drop) +0.0076 [−0.0003, +0.0155] — the last is the frequency effect. Tables T4, T4b.
7. **Letter construction:** bootstrap (b) − (a) +0.0363 [+0.0240, +0.0487] but (b) − always-drop-`kolf2` **−0.0187 [−0.0277, −0.0095]** and (b) − (c) −0.0195 [−0.0293, −0.0097], all three excluding 0 — the arm is significantly WORSE than its own control, as 264 read; the proxy permutation is not usable here (error up to 3.2 SD) and is tabulated (T5b) only to show it. Table T5.
8. **263's two-donor cells** (paired bootstrap vs `cd4svd5`, stratified 3-fold mean): `cd4d:raw:tau0:rev:nm` +0.00178 [−0.0038, +0.0076]; `union:raw:tau0:rev` +0.00157 [−0.0011, +0.0046]; `union:wht:tau0:nm` (the best arbiter) −0.00040 [−0.0024, +0.0019]; the oracle +0.0459 [+0.0383, +0.0538]. Noise-sized, as 263 said. Table T6.
9. **The one-line answer (§6): YES on `kolf2` at p = 0.011 (both nulls; 0.012 with the real total), one fold, τ 0 on the line — and NO on every cancer fold in the only construction the test can read (p 0.68 / 0.12 / 0.50).** The per-cell p-value does not account for the screen that selected the cell (263: the majority row was the best of ~130 cells on `kolf2`; 8 fold × gate cells tested here; Bonferroni over the four folds 0.044, over the eight cells 0.09).
10. **Recommendation:** a correctly re-fitted three-donor validation (three held-out operators fitted with `kolf2` as the third donor, 244's recipe, ~30 min per fold, then 264's sweep and this script) is **not worth building as a priority**: it could at best validate an iPSC third donor on cancer truths, never the production trio; the target-level information the test isolates is at chance on 3/3 cancer folds already; and the `kolf2` reading does not clear the bootstrap against the incumbent. It is acceptable as an overnight compute-only item with the pre-specified criterion in §6, nothing more. No price, no preregistration, no submission.

---

## 1. What was computed (`boot265.py`; inputs read-only; outputs `work/boot265.{json,log}`, `work/tables265.md`)

**Inputs.** 264's `kolf2` sweep (`sweep264-kolf2.json` `f914fe98…`, `pertarget264-kolf2.npz` `66c2fa2b…`, `picks264-kolf2.npz` `68be76f5…`): per-target a23 `pds` of `cd4svd5`, `majority:raw:{all,tau0}:nm`, the reverse, the three random-drop seeds per gate, the three always-drop-j arms, the whitened-space majority and the two oracle drops, plus each arm's drop vector; 263's `kolf2` arrays (`pertarget263-kolf2.npz` `7ad456f4…`, `picks263-kolf2.npz` `9609f831…`, `sweep263-kolf2.json` `81248516…`) for KA-0; 264's three cancer folds in both constructions (`sweep264m-*` `89c83aca…` / `2ede0294…` / `625e5def…`, `pertarget264m-*` `e840aa19…` / `63854541…` / `6110fb9d…`, `picks264m-*` `348efce0…` / `54afc884…` / `7e82ccb9…`; letter `sweep264-*` `ec043232…` / `1e717742…` / `6374ea2e…`, `pertarget264-*` `7ff9707d…` / `3608e37e…` / `3f3d730a…`, picks identical to the matched pass by construction); 263's cancer-fold arrays (`pertarget263-*` `0d0cb1df…` / `2b853554…` / `677b2253…`, `sweep263-*` `2521a0a7…` / `74bc18c2…` / `7654c05a…`).

**The statistic** (`resid218._per_target`, tag a23): per frame A/B/C, `pds_t = (#{q ≠ t: S[t,q] < S[t,t]} + ½ #{q ≠ t: S[t,q] = S[t,t]}) / (n − 1)` on the (rendered prediction × truth) cosine matrix S, averaged over frames; a fold's `pds` is the mean over its n targets (273 / 269 / 290 / 294; KA-1). Per target it is discrete in steps of 1 / (3 (n − 1)) ≈ 0.0012.

**Paired target bootstrap.** For two arms evaluated for real, the per-target difference d_t (exactly 0 where neither arm touched the row) is resampled over the n targets with replacement, 20,000 times; reported: point (= the sweeps' Δ, KA-4), the 2.5 / 97.5 percentiles, the bootstrap SD and the fraction of resamples ≤ 0. The 3-fold mean is resampled within each fold (stratified). This is exact: it uses no additivity.

**Permutation over drop assignments.** On the m gated targets (240; τ 0: 129 / 176 / 177 / 168), P[t, j] = the per-target `pds` of target t in the always-drop-j `:nm` arm; the proxy total of an assignment d is (Σ_t P[t, d_t] + Σ_untouched base_t) / n. Nulls: **uniform** (d_t iid over the three donors) and **frequency-preserving** (a uniformly random permutation of the arm's own drop vector across the gated targets, keeping every donor's drop count), 100,000 draws each; p_upper = (1 + #{null ≥ obs}) / (B + 1), one-sided for a gain; reported with the proxy AND the real total as `obs`; every other arm on the same gate (reverse, random seeds, always-drop, whitened-space majority, the two oracle drops) is positioned in the same nulls. The identity E[uniform] = mean of the always-drop totals (with the untouched-row coupling removed) is checked at 1e-12 (KA-5).

## 2. H-2 / H-3 — additivity, measured (table T7)

| fold | construction | untouched targets | max |Δ`pds`| on an untouched target | max |Δ own cosine| there | max per-target |Δ| among `:nm` arms | max TOTAL error among `:nm` arms | null SD (uniform, gate all) | error / SD |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `kolf2` | production donors | 33 | 0.0123 | 9.5e-4 | 0.017 | **2.6e-4** | 0.00696 | **0.04** |
| K562 | matched | 29 | 0.0547 | 1.3e-2 | 0.132 | 1.1e-3 | 0.00622 | 0.18 |
| HCT116 | matched | 50 | 0.0565 | 1.4e-2 | 0.128 | 1.6e-3 | 0.00607 | 0.27 |
| HEK293T | matched | 54 | 0.1013 | 1.4e-2 | 0.191 | 1.3e-3 | 0.00594 | 0.22 |
| K562 | letter | 29 | 0.4789 | 6.2e-2 | 0.631 | 1.1e-2 | 0.00773 | 1.38 |
| HCT116 | letter | 50 | 0.5029 | 1.1e-1 | 0.730 | 1.1e-2 | 0.00648 | 1.66 |
| HEK293T | letter | 54 | 0.8373 | 1.3e-1 | 0.803 | 2.1e-2 | 0.00641 | 3.23 |

Why it fails, from the code: every `:nm` arm keeps the un-whitened row norms (so `Operator.apply`'s target median `before` = 4.140339 and the resolved fractions 0.336448 / 0.189583 are identical to six decimals across arms, as the sweep audits show), but a pair-mean row whitens to a different un-preserved norm than a three-mean row, the median `after` over the 273 rows moves, and the preserve scalar `before / after` applied to EVERY row moves with it; the untouched rows' amplitude then shifts against the residual add and through `deliver`'s `expm1` / mass fix, and their own cosines move by up to 1e-3 (`kolf2`) — enough to flip a few ranks per target on a rank statistic with steps of 0.0012. The deviations are symmetric (mean Δ on the untouched `kolf2` targets +0.00007 / −0.00063 / −0.00019 for the three drop arms) and cancel in the sum. **Consequence:** the mandate's "0.0e+00 against 263's `majority:raw:all:nm` row" cannot hold on this pipeline for any arm, and per-target `pds` from one arm cannot be attributed to that target's own choice alone; totals can, to the error above. The proxy is used where the error is < 0.3 SD and refused where it is not.

## 3. `kolf2` (tables T1–T3 in `work/tables265.md`)

Bootstrap (exact), Δ`pds` vs the named reference, 95 % interval, P(≤ 0):

| comparison | point | 95 % | P(≤ 0) |
|---|---:|---|---:|
| majority − incumbent | **+0.01209** | [−0.00346, +0.02783] | 0.064 |
| majority − mean(random drop, 3 seeds) | +0.00526 | [−0.01146, +0.02217] | 0.27 |
| majority − random s0 / s1 / s2 | +0.0006 / +0.0117 / +0.0035 | [−0.021, +0.022] / [−0.008, +0.031] / [−0.016, +0.024] | 0.48 / 0.12 / 0.36 |
| majority − always-drop K562 / HCT116 / HEK293T | +0.0152 / +0.0198 / +0.0119 | [−0.0015, +0.0324] / [−0.0019, +0.0410] / [−0.0093, +0.0334] | 0.038 / 0.037 / 0.13 |
| majority − reverse | +0.0218 | [−0.0029, +0.0463] | 0.041 |
| mean(random drop) − incumbent | +0.00683 | [−0.00276, +0.01659] | 0.084 |
| always-drop K562 / HCT116 / HEK293T − incumbent | −0.0032 / −0.0077 / +0.0002 | all include 0 | 0.66 / 0.85 / 0.49 |
| majority τ 0 − incumbent; − mean(random drop τ 0) | +0.00567; +0.00345 | [−0.0067, +0.0180]; [−0.0091, +0.0164] | 0.19; 0.29 |
| whitened-space majority − incumbent | −0.00307 | [−0.0169, +0.0109] | 0.66 |
| oracle drop raw / whitened − incumbent | +0.0542 / +0.0816 | [+0.041, +0.068] / [+0.068, +0.095] | 0.0000 |

Permutation (proxy; Δ vs the incumbent 0.640294):

| arm | null | null mean | null SD | 95 % band | obs proxy / real | p_upper proxy / real | z |
|---|---|---:|---:|---|---|---|---:|
| `majority:raw:all:nm` (110 / 65 / 65) | uniform | −0.00350 | 0.00696 | [−0.0171, +0.0101] | +0.01223 / +0.01209 | **0.0116 / 0.0122** | +2.26 |
| | frequency-preserving | −0.00347 | 0.00685 | [−0.0168, +0.0100] | +0.01223 / +0.01209 | **0.0111 / 0.0117** | +2.29 |
| `majority:raw:tau0:nm` (45 / 41 / 43) | uniform | −0.00331 | 0.00546 | [−0.0141, +0.0074] | +0.00573 / +0.00567 | 0.0494 / 0.0508 | +1.66 |
| | frequency-preserving | −0.00343 | 0.00547 | [−0.0142, +0.0074] | +0.00573 / +0.00567 | 0.0478 / 0.0488 | +1.68 |

Positions of the other arms in the gate-all nulls (uniform / frequency-preserving p_upper, proxy; real in T2b): reverse 0.82 / 0.83; random s0 0.015 / 0.014, s1 0.29 / 0.28, s2 0.039 / 0.037; whitened-space majority 0.48 / 0.48; always-drop K562 / HCT116 / HEK293T 0.48 / 0.72 / 0.30; oracle raw / whitened 1e-5 / 1e-5 (z +8.3 / +12.2). The majority agrees with the whitened oracle's drop on 39 % of targets (chance 33 %); its gain by dropped donor is K562 (110) +0.81, HCT116 (65) −0.03, HEK293T (65) +2.53 in summed per-target `pds` (263 `kolfctrl`'s split, reproduced).

**Reading.** The rule's target-level choices carry information on this fold: shuffling which target gets which dropped donor, with the donor frequencies fixed, beats the rule in 1.1 % of 100,000 shuffles. The proper sampling scale of a drop-one control here is 0.0068, so the majority is +2.3 SD above a random drop, not "8 σ". What the excess is worth against the incumbent is a different question, and the paired bootstrap answers it: +0.0121 with 6.4 % of resamples at or below zero — not cleared. The three real random seeds were an unlucky control (joint p 0.005), which is why 264's "56 %" and the outside re-analysis's interval against mean(random) both understate the rule; the calibration on the 18 cancer-fold seeds (mean p 0.48, one each side of 0.05 / 0.95) says the null itself is right where it can be checked.

## 4. The cancer folds (264's arrays; tables T4, T4b, T5, T5b)

**Matched construction** (`kolf2` scaled to the production donors' whitened amplitude; the proxy's error 0.18–0.27 SD):

| fold | freq-preserving p (proxy / real) | z | uniform p (proxy / real) | always-drop-`kolf2` in the freq null: p | oracle (whitened) z | (b) − (a) bootstrap | (b) − (c) bootstrap |
|---|---|---:|---|---:|---:|---|---|
| K562 | **0.679 / 0.680** | −0.46 | 0.206 / 0.208 | 0.0076 | +10.2 | −0.0039 [−0.0173, +0.0098] | −0.0136 [−0.0284, +0.0009] |
| HCT116 | **0.116 / 0.096** | +1.19 | 0.034 / 0.027 | 0.0179 | +10.6 | −0.0022 [−0.0127, +0.0080] | −0.0035 [−0.0149, +0.0073] |
| HEK293T | **0.502 / 0.577** | 0.00 | 0.155 / 0.197 | 0.0129 | +9.9 | −0.0091 [−0.0223, +0.0032] | −0.0151 [−0.0299, −0.0008] |
| 3-fold mean | | | | | | **−0.00509 [−0.01222, +0.00201]**, P(≤ 0) 0.92 | **−0.01073 [−0.01871, −0.00302]**, P 0.997 |

τ 0: frequency-preserving p 0.60 / 0.074 / 0.43 (real 0.61 / 0.072 / 0.48). Also 3-fold: (a) − (c) −0.00564 [−0.01104, −0.00042]; (b) − always-drop-`kolf2` −0.01042 [−0.01779, −0.00314]; (b) − mean(random drop) +0.00756 [−0.00030, +0.01545] (P 0.030) — a frequency effect, since the same rule is at chance once the frequencies are held. The uniform null on HCT116 (p 0.03) is the only cell under 0.05 on the cancer folds and it is the donor-frequency channel: the majority drops `kolf2` on 130 / 240 targets against 80 expected, and dropping `kolf2` is worth +0.0012 … +0.0122 per fold by itself.

**Letter construction** (equal three-donor mean; the proxy refused, error 1.4–3.2 SD; bootstrap only): (b) − (a) +0.0363 [+0.0240, +0.0487]; (b) − (c) −0.0195 [−0.0293, −0.0097]; (a) − (c) −0.0557 [−0.0705, −0.0410]; **(b) − always-drop-`kolf2` −0.0187 [−0.0277, −0.0095]**; (b) − mean(random drop) +0.0378 [+0.0258, +0.0496]; (b) τ 0 − (a) +0.019 / +0.016 / +0.007 per fold, all excluding 0, and (b) τ 0 − (c) −0.039 / −0.024 / −0.063. Every interval agrees with 264 §4: the arm's gain over (a) is the removal of a junk donor, and always removing it is better.

## 5. 263's two-donor cells on the cancer folds (table T6; paired bootstrap vs `cd4svd5`, per fold and the stratified 3-fold mean)

| cell | K562 | HCT116 | HEK293T | 3-fold mean [95 %] | P(≤ 0) |
|---|---|---|---|---|---:|
| `cd4d:raw:tau0:rev:nm` | +0.0101 [−0.0012, +0.0222] | −0.0043 [−0.0126, +0.0040] | −0.0004 [−0.0100, +0.0090] | +0.00178 [−0.00382, +0.00756] | 0.27 |
| `union:raw:tau0:rev` | +0.0001 [−0.0057, +0.0067] | +0.0008 [−0.0009, +0.0032] | +0.0038 [−0.0008, +0.0101] | +0.00157 [−0.00112, +0.00463] | 0.14 |
| `union:wht:tau0:nm` (best arbiter, 263 table A) | +0.0023 [−0.0020, +0.0084] | −0.0017 [−0.0042, +0.0001] | −0.0019 [−0.0047, +0.0006] | −0.00040 [−0.00236, +0.00192] | 0.66 |
| `basal:raw:tau0:nm` | −0.0023 | −0.0000 | −0.0016 | −0.00129 [−0.00672, +0.00399] | 0.68 |
| `random:tau0:s0:nm` | +0.0001 | +0.0021 | −0.0024 | −0.00008 [−0.00708, +0.00685] | 0.51 |
| `oracle_wht:tau0:nm` / `oracle_wht:all:nm` (ceilings) | +0.0307 / +0.0490 | +0.0184 / +0.0468 | +0.0217 / +0.0417 | +0.02361 [+0.01794, +0.02971] / +0.04586 [+0.03834, +0.05378] | 0.0000 |

Every observable cell's interval includes 0 on the 3-fold mean and on every fold; the two "positive reverse" cells are +0.0016 … +0.0018 with intervals [−0.004, +0.008] and [−0.001, +0.005]. 263 §6's "noise-sized" stands, now with the interval that says so.

## 6. The one line, and what follows

**Does the majority rule's excess over frequency-preserving random drops survive at p < 0.05 on `kolf2`? YES: p = 0.011 (0.012 with the real total), z +2.3, one fold, gate all; the τ = 0 gate is on the line (p 0.048 / 0.049 real 0.051); and the same test on the three cancer folds, in 264's matched construction, reads p 0.68 / 0.12 / 0.50 — chance.**

Since the answer on `kolf2` is yes, what a correctly re-fitted three-donor validation would need:

- **Three held-out whitening operators**, one per cancer fold, each fitted on the union of that fold's two production donors AND `kolf2` (the only third donor on disk), with 244's recipe (`scripts/build_union_whitener.py`, the `union_pruned` pruning, the held line excluded; ~30 min per fold per the mandate, so ~1.5 h of fits), recorded by path and sha, never by name (A:1497).
- **264's `sweep264.py` re-run against each re-fitted operator** (`lib264.stage3` with `L.op_path_for(held)` re-pointed; the fold caches' `pre` / `truth` unchanged; R1 / R3 must still reproduce at 0.0e+00 with the production operator before the swap), equal and matched passes (with a re-fit the two should converge), then this experiment's `boot265.py` on its arrays.
- **A pre-specified criterion, written before the fits**: (b) − (c) 3-fold bootstrap interval excluding 0 on the positive side AND the frequency-preserving permutation at p < 0.05 on ≥ 2 of 3 folds AND (b) − always-drop-`kolf2` ≥ 0 on the 3-fold mean. Anything less is 264's answer with a better operator.
- Cost: compute only (~2 h, ≤ 6 GiB per fold process for the sweeps; the operator fit's memory is 244's), no slot.

**Recommendation: not as a priority.** (i) Even a 3/3 pass would validate "majority with an iPSC third donor" on cancer truths — the production rule's third donor is a cancer line that no local fold can hold out and keep, so the transfer to the production trio would remain by analogy and the 3/3-fold condition for a build would still not be met in its own terms. (ii) The target-level information the test isolates is at chance on all three cancer folds in the construction that already removes the amplitude artefact, and the rule loses to the incumbent there with an interval excluding 0. (iii) On `kolf2` itself, the excess over the incumbent is +0.0121 with P(≤ 0) 0.064, and p = 0.011 is a per-cell value for a cell that 263 selected as the best of ~130 on that fold (Bonferroni over the eight fold × gate cells tested here alone: 0.09). (iv) 263's own break-even arithmetic (a selector must be right two times in three; the majority is right 39 % of the time against the oracle on `kolf2`, chance 33 %) says the mechanism needs an arbiter that knows the held line, which is 264 §0's conclusion unchanged. If idle compute is available overnight it is an acceptable low-priority item under the criterion above; it is not a path to a 09-06 slot. **H2's majority variant does not re-open on these numbers.** No price, no preregistration, no submission.

## 7. Honest limits

- **The permutation test is on an additive proxy.** Strict additivity fails (§2); the proxy's error on totals is measured on every real arm of each gate (13 per fold) and is 0.04 SD on `kolf2`, 0.18–0.27 SD on the matched folds, and is refused on the letter folds. Every p-value is given for the proxy and for the real total; no conclusion changes between the two. An exact permutation would need real evaluations of thousands of assignments (263's staging: ~2 s per arm at 5.4 GiB RSS) — outside this mandate's scope and its 4 GiB box.
- **The `kolf2` random seeds.** Six real seeds (three per gate) sit high in a null that is calibrated on 18 cancer-fold seeds; the joint probability of their positions is 0.005 (gate all) and 0.04 (τ 0). Nothing in their drop vectors is anomalous. The reading "an unlucky draw" is the only one consistent with the proxy check (each seed's real total is within 2e-4 of its proxy), but it is a 1-in-200 event and is recorded as such.
- **One fold carries the yes.** `kolf2` is an iPSC truth with a base at `pds` 0.640 and donor agreement centred on zero (262 §2, 263 §7); the cancer folds test a different third donor through an operator that has not seen it (264 §3, §10). Nothing here measures the production trio on a cancer line.
- **The bootstrap resamples targets, not lines.** It is the uncertainty of a fold mean under target resampling on one held line; the line-to-line variation is what the 3/3 rule is for, and READ-FIRST §3's three labels apply: this is a fitted-free re-analysis of cached rows, not an independent validation of anything.
- **Selection.** The majority cell was reported by 263 because it was the best positive among ~130 cells on `kolf2`; a per-cell p-value cannot be read as if the cell were pre-specified.
- **Discreteness.** Per-target `pds` moves in steps of ≈ 0.0012; the untouched-row deviations in §2 are rank flips of that size, which is why they cancel.

## 8. Files, memory, provenance

| file | what | sha256 (12) |
|---|---|---|
| `boot265.py` | the whole analysis: loaders, KA rows, H-2 / H-3, the paired bootstrap (20,000), the two permutation nulls (100,000), the seed diagnostics, the calibration table, the markdown tables | `8d845c9e28c3` |
| `work/boot265.json` | every number in this file (bootstrap rows, permutation entries with the positions of every other arm, additivity per arm, seeds, calibration, the KA log) | `ad3576965ac8` |
| `work/boot265.log` (= `work/boot265.nohup`) | the run's transcript | `73ca7554449e` |
| `work/tables265.md` | T1–T7 as generated | `caacba96b13d` |
| `work/.marker-start` | the write-audit marker (2026-09-05 12:42Z) | — |

Peak RSS 0.40 GiB, 22.3 s, one process. Nothing outside this directory was created or modified; no fold cache copied, no operator loaded, no builder import, no `vcc`, no `cell-eval2`, no submission, no preregistration.

---

## 9. Verdict

**The majority rule's target-level choices carry information on `kolf2` (frequency-preserving permutation p = 0.011, z +2.3 against the proper null SD 0.0068), and nothing else survives: the excess over the incumbent is +0.0121 [−0.0035, +0.0278] (P(≤ 0) 0.064), the same permutation is at chance on 3/3 cancer folds (p 0.68 / 0.12 / 0.50) where the rule loses to the incumbent with an interval excluding 0 (−0.0107 [−0.0187, −0.0030]) and to always dropping `kolf2`, 264's "56 %" control was three seeds from the null's upper tail, 263's positive reverse cells include 0, and the mandate's additivity premise is false for the pipeline (one preserve scalar couples every row) while its totals-level consequence holds to 0.04 SD on `kolf2`. A re-fitted three-donor validation is specified in §6 and not recommended as a priority; H2's majority variant stays closed. Measurement only; no build, no submission.**
