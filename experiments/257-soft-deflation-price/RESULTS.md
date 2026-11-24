# 257 — A SOFTER CD4 deflation does not beat `cd4svd5`. Every cell of the s × k grid and every single-direction cell LOSES `pds` against the incumbent; the cosine it gives back is worth less than the `pds` it gives up at the board's factors. Best cell (k = 5, s = 0.75) prices **+0.00055 overall at the conservative factors, 0.29× the +0.0019 gate**, negative at k_cos 2.15 / 1.36. **DOES NOT CLEAR — verdict line at the bottom (§10).**

**Date 2026-09-05 (01:45–02:05 UTC).** Measurement only. **No build, no `vcc`, no `cell-eval2`, nothing under
`submissions/`, no model name, no description, no preregistration.** No file written outside this directory (a
`find -newer` over `experiments/ data/ scripts/ src/ docs/ submissions/` against `work/.marker-start` returns only
the sibling 258's own cache copies). `uv run python` throughout; one PLOLO fold per process, never alongside another
≥ 15 GiB job (no builder / `vcc prep` running; `ps` checked; memory 68 % free at start). Peak RSS per script in §9.

| row | required | measured | verdict |
|---|---|---|---|
| **R1** `cd4svd5` vs `h1union` through the recipe | overall +0.0043 ± 0.0006, `pds` +0.0055 ± 0.0006, `mse` −0.0008 ± 0.0003 | **overall +0.00472, `pds` +0.00552, `mse` −0.00083** (board +0.0042634 / +0.00545 / −0.00080); DE-three −0.00035 (board −0.00034), `nmae` +0.00038 (board −0.00005); Δcos_board −1.736 % (board −1.735 %) | **PASS** |
| **R2** `k5_s1` == 246 `B2_cd4_svd5_resolved` | Δ`pds` +0.012263, agg −0.4617 %, macro −0.8087 % to the digit | **+0.012262757 / −0.461748 % / −0.808675 %, max\|Δ\| 0.0e+00** vs `price246.json` and 249's `rel_agg_plolo`; the deflated arrays bit-identical to `lib246.deflate("svd", k)` for k = 3, 5, 8, 10 on all 4 folds | **PASS** |
| **R3** `k10_s1` vs `k5_s1` | 255's increment (`pds` +0.00054, agg −1.07 %), prices negative at the point factors | **+0.000544385 / −1.0671 % / macro −0.8621 %, max\|Δ\| 0.0e+00** vs `repro255-summary.json`; price **−0.00212** point (`pds` +0.00025, `mse` −0.00148, DE-three −0.00079, `nmae` −0.00010), −0.00248 conservative | **PASS** |
| sweep known-answer rows | 1e−9 | **24/24 at 0.0e+00** (per fold: STACK, `B2_cd4_svd{3,5,10}_resolved`, `B2_cd4_randE5` fixed, 249's resolved control) | **PASS** |
| identity known-answer rows (`six246.json`, 4 folds) | 5e−6 | **9/9 at 0.0e+00** (`a15`, STACK, `k5_s1`: `nmae_wald89`, top-89 precision, delivered ratio) | **PASS** |
| **the gate** | any cell > +0.0019 over `cd4svd5` at the conservative factors, shape-only | **NONE.** Best `k5_s0.75` +0.00055 (point +0.00023; k_cos 2.15 −0.00021; 1.36 −0.00042). Every other cell ≤ +0.00048 conservative; 15 of 16 negative at the point factors | **DOES NOT CLEAR** |

---

## 0. The answer in ten lines

1. **Softening is a `pds`-for-cosine trade at roughly one-to-one in board units, and `cd4svd5` already sits at the
   break-even.** Along s at k = 5 (vs the shipped stack) `pds` grows almost linearly — +0.0043 / +0.0077 / +0.0104 /
   +0.0123 at s = 0.25 / 0.5 / 0.75 / 1 (35 / 63 / 85 / 100 % of the full gain, 3/3 folds each) — while the
   aggregate cosine is flat or positive up to s = 0.75 (+0.20 / +0.29 / +0.01 %) and **the entire −0.46 % arrives in
   the last quarter of the removal** (§4.1).
2. So the softest cells lose `pds` at 3–8 PLOLO-σ for a cosine gain worth +0.0010–0.0017 overall at k_cos 3.76:
   net −0.0026 … −0.0019 at s = 0.25. The best cell, **k = 5, s = 0.75**, trades `pds` −0.00183 raw (1/3 folds, −1.2 σ) for
   +0.48 % aggregate cosine: **+0.00023 point / +0.00055 conservative**, and it turns negative as soon as the cosine
   transfer is read at 2.15 or 1.36 instead of 3.76 (§5).
3. **k = 8 is dominated by k = 5** (s = 1: `pds` −0.00065 AND cosine −0.58 % vs `cd4svd5`, −0.0016 overall); k = 3 at
   any s loses `pds` (−0.0033 … −0.0087) for +0.3–0.8 % cosine, −0.0026 … +0.0001 (§3).
4. **A hard bound on the s axis**: everything between s = 0.75 and 1 can at most recover the +0.48 % cosine (worth
   +0.0010 at 3.76×) without losing the `pds` it gives up; even that impossible best is 0.53× the gate (§5.3). A
   finer s grid cannot clear.
5. **Per direction** (k = 5, one direction removed alone): direction 1 (8.0 % of the centred energy) carries
   **+0.0064 of the +0.0123** `pds` at +0.12 % aggregate cosine (3/3 folds); directions 2–5 alone are worth
   −0.0001 … +0.0018; the five singles sum to +0.0099, i.e. the projector is **1.24× super-additive** — no
   single direction reproduces `cd4svd5`, and every single-direction cell is −0.0059 … −0.0124 `pds` against it (§4.2).
6. The cosine cost is not carried by one direction alone (singles: −0.41 … +0.12 %), but the **exploratory**
   leave-one-out family (§6, post hoc, not part of the fixed grid) shows **removing direction 4 costs +0.38 % aggregate
   cosine on 3/3 folds (+0.40 / +0.39 / +0.37 %) and buys no `pds`** (−0.0005, −0.3 σ): the projector without it
   prices +0.00059 / +0.00074, **0.39× the gate**. Direction 1 alone carries 75 % of the `pds` (dropping it: −0.0092).
7. Energy-matched random subspaces at every cell's energy buy no `pds` (−0.0055 … +0.0003 vs the stack) and cost
   0.3–3.6 % cosine: the cells' `pds` is direction-specific at every s (§4.3). Leave-top-3-out retention of the best
   cell vs the stack is 61.9 % (`cd4svd5` 52–56 %); of its increment over `cd4svd5` it is undefined (the loss is
   diffuse, 160 % / −554 % / 138 % per fold).
8. `nmae` by the validated identity moves ≤ ±0.0004 overall on every grid cell (dir-3 −0.0007); q_w and top-89
   precision move ≤ 0.001. Delivered amplitude: the softest cells sit at +0.50–0.59 % by PLOLO's proxy (beyond the
   shape-only bound; they are negative anyway); the best cell +0.12 %.
9. The recipe's fixed constants against the captured board entries: `pds` scaled/raw 2.2138 vs 2.1925
   (`h1union`→`cd4svd5`) and 2.2385 (`h1union`→ARM6); `mse` 1.051 vs 1.0558 / 1.0512; d raw/d cos −1.996 vs the
   exact −1.9943 at `cd4svd5`'s terms; none changes a sign or a verdict (§1.3).
10. **DOES NOT CLEAR.** BR249 §5's hypothesis (b) — "a deflation that spends less cosine per unit `pds`" — is closed
    on A/B/C: the projector's cosine cost lives in the last quarter of the removal and in direction 4, and the board's
    factors (`pds` 1.22×, cosine 3.76×) price that quarter within noise of break-even. No `PREREGISTRATION-DRAFT.md`.

---

## 1. Instrument and recipe (`price257.py`; fixed before any grid number was read)

### 1.1 Instrument
162's PLOLO, three cancer folds (`k562_gwps`, `orion_hct116`, `orion_hek293t`; `kolf2` reported apart), 246's
construction imported (`lib246`: 218's residuals, the shipped stack `cd4@0.34 + union(gse295214_h1, h1_2025)@0.28`,
`compose`, `score_effects`), fold caches byte copies of 246's (`work/plolo/`, shas `c13c1cec…`, `9f3e0452…`,
`f9a3cc30…`, `d86202a9…`; `bctrl-challenge-{A,B,C}.npy` `f1f04360…`, `4db8bd10…`, `cdebb5ea…`). Every cell is the
**resolved** construction: the deflated CD4 residual re-sized by `solve_fraction` (production's all-rows rule at
`f` 0.34) and merged with the union at its shipped fraction, one add, one preserve — exactly the builder's
`apply_frozen_cd4_residual` order (249 §3 proved the two paths equal at 0.0e+00). Endpoints per cell: `pds`,
`cos_agg` (243's functional), macro cosine, q_w, per-target `pds`/cosine; 188's exact `nmae` identity (wald89 pool,
frame A, no mass fix, six246's construction), top-89 precision, delivered norm (median row norm of the frame-A render
on the operator columns).

### 1.2 The operator (`lib257.deflate_shrink`)
`R_c = R − mean_rows(R)` over the carried rows (support re-zeroed), `U_k` = top-k eigenvectors of `R_c R_cᵀ`,
`R_s = R_c − s·U_S U_Sᵀ R_c` (S = all top-k, or one direction), support re-zeroed. At s = 1, S = top-k the expression
is `lib246.deflate`'s own line; the arrays are **bit-identical** to `lib246.deflate(res, "svd", k)` for k = 3, 5, 8,
10 on every fold (asserted). Energy-matched random control: `sweep246b.rand_energy_matched` verbatim (246's `rand`
mode, seed `RAND_SEED + 5`, bisection to the cell's removed energy), resolved sizing; the k = 5 one reproduces 246's
`B2_cd4_randE5` (52/54/55/52 dims) and 249's resolved control at 0.0e+00.

### 1.3 The price (mandate's recipe; constants as fixed)
| leg | formula | constants |
|---|---|---|
| `pds` | 2.2138 × k_pds × Δ`pds`(3 cancer folds, vs `cd4svd5`) / 6 | k_pds **1.22** point, **0.85** conservative |
| `mse` | Δcos_board = k_cos × rel `cos_agg` × 0.105470; Δraw = −1.996 × Δcos_board (+ 243's amplitude leg: 197's identity at `cd4svd5`'s own `A_c`/`B_c`/`Kc_c`/`cap_c` from `amp249.json`, A and B rescaled by the delivered ratio², when the delivered move exceeds 0.2 %); overall = −1.051 × Δraw / 6 | k_cos **3.76** point (re-shape class, BR249 §2); 2.15 (macro read) and 1.36 (re-dose) as sensitivity |
| DE-three (`reach`, `fid`, `jac`) | 0.24 × Δcos_board × share, share = the three members' part of 218's DE-four rate card = 0.776 (point); the whole 0.24 (conservative); split among the three by the card | 243 §4.3's union-point propagation |
| `nmae` | 216/188's exact identity, −Δ`nmae_wald89` / 0.60 / 6 (the three cancer folds concatenated) | ×1.0 point; ×1.5 on a loss and **0 credit on a gain** conservative (BR237: a lower bound on the loss) |

Sensitivity of the fixed constants against the captured entries (`data/board-entries/`): `pds` scaled-per-raw
**2.1925** on `h1union`→`cd4svd5` and **2.2385** on `h1union`→ARM6 (recipe 2.2138); `mse` **1.0558** / **1.0512**
(recipe 1.051); 197's slope at `cd4svd5`'s terms **−1.9943**, at `h1union`'s −2.0025 (recipe −1.996). Each is within
1 % and the exact-identity `mse` form reproduces every total to 1e−5 (`price257.json` `exact_identity`).

### 1.4 Retrodiction rows
- **R2** (`k5_s1` vs STACK, three cancer folds): Δ`pds` **+0.012262757** (per fold +0.00999 / +0.00945 / +0.01736,
  3/3, 8.2 σ), aggregate **−0.461748 %**, macro **−0.808675 %**, energy removed 0.2008, CD4 fraction ×1.0602 — equal to
  `price246.json`'s `B2_cd4_svd5_resolved` and to 249's `rel_agg_plolo` at **0.0e+00**.
- **R1** (the same increment through the recipe, base `h1union`: cos 0.107332, its own terms): Δcos_board −0.001863
  (−1.736 %; board −1.735 %); Δraw `mse` = cosine +0.00372 + amplitude +0.00103 (leg applied: PLOLO's delivered ratio
  reads **−0.533 %** on the three cancer folds against the package's measured −0.409 %; six246's 4-fold read was
  −0.317 %) = +0.00475 (board +0.00454). Legs: `pds` **+0.00552** (band +0.0055 ± 0.0006; board +0.00545), `mse`
  **−0.00083** (band −0.0008 ± 0.0003; board −0.00080), DE-three −0.00035 = `reach` −0.00025 / `fid` −0.00009 /
  `jac` −0.00001 (board −0.00026 / −0.00006 / −0.00003), `nmae` +0.00038 (identity, a gain; board −0.00005).
  **Overall +0.00472** (band +0.0037 … +0.0049; board +0.0042634). Conservative factors: +0.00257. At k_cos 2.15 /
  1.36: +0.00515 / +0.00536. Check: 197's identity at k_cos 3.76 with the package's own terms returns raw `mse`
  **0.82181 = the board's 0.82181**; `h1union`'s terms reproduce its 0.8172630 to 9e−7.
- **R3** (`k10_s1` vs `k5_s1`): Δ`pds` **+0.000544385** (1/3 folds: −0.00017 / +0.00558 / −0.00378; 0.36 σ),
  aggregate **−1.0671 %**, macro −0.8621 % — 255 §2.1's increment at 0.0e+00. Price at the point factors: `pds`
  +0.00025, `mse` −0.00148 (Δcos_board −0.004232, Δraw +0.00845), DE-three −0.00079, `nmae` −0.00010 = **−0.00212**;
  conservative −0.00248 (BR249 §3 estimated ≈ −0.0014 `mse` against +0.0002 `pds`). Negative: **PASS**.

---

## 2. Known-answer rows (`sweep257.py`, one fold per process; `work/sweep257-<fold>.{json,log}`)

24/24 at 0.0e+00 over four folds: STACK vs 246 `sweep246`; `k3_s1` / `k5_s1` / `k10_s1` vs 246
`B2_cd4_svd{3,5,10}_resolved` (endpoints AND the re-solved fractions, 1e−12); the deflated arrays for k = 3, 5, 8, 10
bit-identical to `lib246.deflate("svd", k)` with identical energy audits; `randE_k5_s1_fixed` vs 246 `B2_cd4_randE5`
(dims 52/54/55/52 = 246's); `randE_k5_s1` (resolved) vs 249's `repro249-<fold>.json` control (−0.00016 / −0.00541 /
+0.00144 per fold, mean −0.001377 = 249 KA 3). Identity rows (`price257.py`, 4 folds concatenated as six246 does):
`a15` 0.98475707 / 0.58994688 / 1.0, STACK 0.98569263 / 0.60312409 / 1.00333918, `k5_s1` 0.98456099 / 0.60265859 /
1.00016210 — 9/9 equal to `six246.json` at 0.0e+00. `sweep257b.py` (exploratory) 8/8.

---

## 3. The grid — every cell as an increment over `cd4svd5` (`k5_s1`), three cancer folds

Columns: energy removed (fraction of the CD4 residual's energy, mean of 3 folds); CD4 fraction relative to
`cd4svd5`'s re-solved 0.3960; Δ`pds` raw vs `cd4svd5` with the count of positive folds and its size in PLOLO LOO-σ
(0.001494); relative `cos_agg` and macro cosine; delivered amplitude (PLOLO proxy); shape-only = within 0.5 %; the
six-member legs (point k_pds 1.22 / conservative 0.85; `mse` at k_cos 3.76 identical in both; DE-three point / cons;
`nmae` point / cons); the totals; the totals at k_cos 2.15 and 1.36 (point k_pds); `kolf2` apart.

| cell | k | s / dirs | E removed | CD4 fraction × | Δ`pds` vs `cd4svd5` (raw) | folds + | σ | rel `cos_agg` | rel macro | delivered | shape-only | `pds` 1.22× / 0.85× | `mse` | DE-three (pt / cons) | `nmae` (pt / cons) | **overall point** | **overall conservative** | k_cos 2.15 | k_cos 1.36 | `kolf2` Δ`pds` / agg |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---|---|---:|---:|---:|---:|---|
| `k3_s0.25` | 3 | 0.25 | 0.066 | 0.9630 | **−0.008742** | 0/3 | −5.9 | +0.710 % | +0.823 % | +0.526 % | NO | −0.00393 / −0.00274 | +0.00117 | +0.00052 / +0.00068 | −0.00035 / −0.00052 | **−0.00259** | **−0.00142** | −0.00324 | −0.00356 | −0.00312 / +3.16 % |
| `k3_s0.5` | 3 | 0.5 | 0.113 | 0.9772 | **−0.006023** | 0/3 | −4.0 | +0.781 % | +0.731 % | +0.225 % | yes | −0.00271 / −0.00189 | +0.00116 | +0.00058 / +0.00074 | −0.00034 / −0.00052 | **−0.00132** | **−0.00050** | −0.00203 | −0.00238 | −0.00164 / +2.57 % |
| `k3_s0.75` | 3 | 0.75 | 0.142 | 0.9845 | **−0.004310** | 0/3 | −2.9 | +0.645 % | +0.442 % | −0.017 % | yes | −0.00194 / −0.00135 | +0.00089 | +0.00048 / +0.00061 | −0.00002 / −0.00002 | **−0.00059** | **+0.00013** | −0.00117 | −0.00146 | −0.00189 / +1.74 % |
| `k3_s1` | 3 | 1 | 0.151 | 0.9871 | **−0.003269** | 0/3 | −2.2 | +0.305 % | −0.053 % | −0.097 % | yes | −0.00147 / −0.00103 | +0.00042 | +0.00023 / +0.00029 | −0.00002 / −0.00004 | **−0.00085** | **−0.00035** | −0.00113 | −0.00126 | −0.00207 / +0.71 % |
| `k5_s0.25` | 5 | 0.25 | 0.088 | 0.9659 | **−0.007985** | 0/3 | −5.3 | +0.693 % | +0.928 % | +0.543 % | NO | −0.00359 / −0.00250 | +0.00115 | +0.00051 / +0.00066 | +0.00006 / 0 | **−0.00188** | **−0.00070** | −0.00251 | −0.00282 | −0.00216 / +3.12 % |
| `k5_s0.5` | 5 | 0.5 | 0.151 | 0.9880 | **−0.004515** | 0/3 | −3.0 | +0.777 % | +0.948 % | +0.230 % | yes | −0.00203 / −0.00142 | +0.00116 | +0.00057 / +0.00074 | +0.00004 / 0 | **−0.00026** | **+0.00048** | −0.00097 | −0.00132 | −0.00044 / +2.38 % |
| **`k5_s0.75`** (best) | 5 | 0.75 | 0.188 | 0.9944 | **−0.001829** | 1/3 | −1.2 | +0.482 % | +0.616 % | +0.118 % | yes | −0.00082 / −0.00057 | +0.00067 | +0.00036 / +0.00046 | +0.00003 / 0 | **+0.00023** | **+0.00055** | −0.00021 | −0.00042 | +0.00026 / +1.28 % |
| **`k5_s1`** (= `cd4svd5`) | 5 | 1 | 0.201 | 1.0000 | 0 | — | — | 0 | 0 | 0 | yes | 0 | 0 | 0 | 0 | **0** | **0** | 0 | 0 | 0 |
| `k8_s0.25` | 8 | 0.25 | 0.112 | 0.9713 | **−0.007598** | 0/3 | −5.1 | +0.614 % | +0.904 % | +0.515 % | NO | −0.00342 / −0.00238 | +0.00103 | +0.00045 / +0.00058 | +0.00005 / 0 | **−0.00188** | **−0.00077** | −0.00244 | −0.00272 | −0.00204 / +3.05 % |
| `k8_s0.5` | 8 | 0.5 | 0.191 | 1.0000 | **−0.004154** | 1/3 | −2.8 | +0.592 % | +0.861 % | +0.211 % | yes | −0.00187 / −0.00130 | +0.00089 | +0.00044 / +0.00056 | +0.00002 / 0 | **−0.00052** | **+0.00015** | −0.00106 | −0.00132 | −0.00016 / +2.16 % |
| `k8_s0.75` | 8 | 0.75 | 0.239 | 1.0150 | **−0.001671** | 1/3 | −1.1 | +0.157 % | +0.399 % | +0.220 % | yes | −0.00075 / −0.00052 | +0.00029 | +0.00012 / +0.00015 | −0.00002 / −0.00003 | **−0.00036** | **−0.00011** | −0.00050 | −0.00057 | +0.00122 / +0.88 % |
| `k8_s1` | 8 | 1 | 0.255 | 1.0247 | **−0.000645** | 1/3 | −0.4 | −0.584 % | −0.478 % | +0.079 % | yes | −0.00029 / −0.00020 | −0.00081 | −0.00043 / −0.00056 | −0.00009 / −0.00013 | **−0.00162** | **−0.00169** | −0.00108 | −0.00082 | +0.00111 / −0.67 % |
| `k5_dir1` | 5 | dir 1 only, s = 1 | 0.080 | 0.9676 | **−0.005864** | 0/3 | −3.9 | +0.608 % | −0.021 % | +0.030 % | yes | −0.00264 / −0.00184 | +0.00084 | +0.00045 / +0.00058 | −0.00035 / −0.00052 | **−0.00170** | **−0.00094** | −0.00225 | −0.00252 | −0.00546 / +1.32 % |
| `k5_dir2` | 5 | dir 2 only, s = 1 | 0.039 | 0.9556 | **−0.012410** | 0/3 | −8.3 | +0.339 % | +0.765 % | +0.590 % | NO | −0.00559 / −0.00389 | +0.00067 | +0.00025 / +0.00032 | +0.00008 / 0 | **−0.00458** | **−0.00290** | −0.00489 | −0.00504 | −0.00482 / +3.42 % |
| `k5_dir3` | 5 | dir 3 only, s = 1 | 0.033 | 0.9508 | **−0.010964** | 0/3 | −7.3 | +0.259 % | +0.652 % | +0.205 % | yes | −0.00494 / −0.00344 | +0.00043 | +0.00019 / +0.00025 | −0.00072 / −0.00108 | **−0.00504** | **−0.00385** | −0.00527 | −0.00539 | −0.00110 / +3.30 % |
| `k5_dir4` | 5 | dir 4 only, s = 1 | 0.030 | 0.9501 | **−0.011734** | 0/3 | −7.9 | +0.091 % | +0.627 % | +0.500 % | yes (at the bound) | −0.00528 / −0.00368 | +0.00030 | +0.00007 / +0.00009 | +0.00012 / 0 | **−0.00480** | **−0.00330** | −0.00488 | −0.00492 | −0.00326 / +3.70 % |
| `k5_dir5` | 5 | dir 5 only, s = 1 | 0.021 | 0.9446 | **−0.010415** | 0/3 | −7.0 | +0.451 % | +0.820 % | +0.519 % | NO | −0.00469 / −0.00327 | +0.00080 | +0.00033 / +0.00043 | −0.00031 / −0.00046 | **−0.00386** | **−0.00250** | −0.00427 | −0.00447 | −0.00225 / +2.88 % |
| `k10_s1` (255's arm, R3) | 10 | 1 | 0.283 | 1.0442 | **+0.000544** | 1/3 | +0.4 | −1.067 % | −0.862 % | −0.070 % | yes | +0.00025 / +0.00017 | −0.00148 | −0.00079 / −0.00102 | −0.00010 / −0.00015 | **−0.00212** | **−0.00248** | −0.00115 | −0.00068 | +0.00075 / −1.49 % |

Per-fold Δ`pds` vs `cd4svd5` (k562 / hct116 / hek293t): `k5_s0.75` −0.00176 / +0.00029 / −0.00401; `k8_s0.75`
−0.00204 / +0.00378 / −0.00675; `k8_s1` −0.00065 / +0.00251 / −0.00380; `k3_s1` −0.00406 / −0.00124 / −0.00451;
`k5_dir1` −0.00779 / −0.00270 / −0.00710. Per-fold aggregate cosine vs `cd4svd5`: `k5_s0.75` +0.798 / +0.884 /
−0.235 %; `k5_s0.5` +1.444 / +1.515 / −0.629 %; `k8_s1` −0.812 / +0.280 / −1.218 % (the hek293t fold, whose donors are
K562 + HCT116, is the one where the projector GAINS cosine, 246 §6; softening gives that back). Full per-fold
tables in `work/price257.json` (`grid.<cell>.vs_cd4svd5`) and the log.

**Nothing clears at any factor.** Sixteen candidate cells: 15 negative at the point factors (only `k5_s0.75`
+0.00023 is positive); four positive at the conservative factors (`k5_s0.75` +0.00055, `k5_s0.5` +0.00048, `k8_s0.5`
+0.00015, `k3_s0.75` +0.00013), all ≤ 0.29× the gate and all negative at k_cos ≤ 2.15. Five cells (`s` = 0.25 at
every k, `k5_dir2`, `k5_dir5`) sit at +0.52–0.59 % delivered amplitude by the PLOLO proxy (`k5_dir4` at +0.500 %, on
the bound) — outside the shape-only class by the letter — and are among the most negative cells anyway.

---

## 4. Reading the grid

### 4.1 The s axis: `pds` is linear in the removal, the cosine cost is not
k = 5, every cell vs the SHIPPED STACK (`h1union`), three cancer folds, 3/3 folds positive on `pds` at every s:

| s | energy removed | CD4 fraction × | Δ`pds` vs STACK | share of the full gain | agg cosine vs STACK (per fold) | macro | delivered | random control (dims; Δ`pds`, agg) |
|---:|---:|---:|---:|---:|---|---:|---:|---|
| 0.25 | 0.088 | 1.0240 | +0.004278 | 35 % | **+0.200 %** (−0.11 / −0.10 / +0.81) | +0.095 % | +0.007 % | 24; −0.00145, −1.03 % |
| 0.5 | 0.151 | 1.0474 | +0.007748 | 63 % | **+0.293 %** (−0.21 / −0.47 / +1.55) | +0.120 % | −0.304 % | 40; −0.00347, −2.27 % |
| 0.75 | 0.188 | 1.0543 | +0.010433 | 85 % | **+0.009 %** (−0.84 / −1.09 / +1.96) | −0.204 % | −0.415 % | 50; −0.00100, −1.93 % |
| 1 (`cd4svd5`) | 0.201 | 1.0602 | +0.012263 | 100 % | **−0.462 %** (−1.63 / −1.95 / +2.20) | −0.809 % | −0.533 % | 53; −0.00138, −2.82 % |

The `pds` increments per quarter are +0.0043 / +0.0035 / +0.0027 / +0.0018 (mildly concave); the cosine is flat to
s = 0.75 and then loses 0.47 % in the last quarter. At the board's factors the last quarter is priced `pds`
+0.00082 (1.22×) against cosine −0.00067 (`mse`) −0.00036 (DE-three) = −0.00103: **`cd4svd5` is within +0.0002 of
break-even, i.e. it sits on the knee** — the +0.00023 point price of `k5_s0.75` is that difference. The same holds at
k = 3 (agg +0.21 / +0.29 / +0.16 / −0.17 %; `pds` +0.0035 / +0.0062 / +0.0080 / +0.0090) and k = 8 (+0.12 / +0.10 /
−0.32 / −1.05 %; +0.0047 / +0.0081 / +0.0106 / +0.0116). 246 §3.5's exchange rate of 0.0266 `pds` per −1 %
aggregate is therefore an AVERAGE over a removal whose first three quarters are free and whose last quarter costs
0.0039 per −1 % — still 4× the dose axis (0.00093), but no longer 29×.

### 4.2 The k axis and the directions
- k = 8 vs k = 5 at s = 1: `pds` −0.00065 (1/3), cosine −0.58 % — dominated. k = 10 (R3): +0.00054 (1/3) for
  −1.07 %. The `pds` ladder saturates at k ≈ 5 (246 §3.3) while directions 6–10 keep costing cosine.
- One direction removed alone (k = 5 basis, s = 1), vs STACK:

| direction | share of centred energy | Δ`pds` vs STACK (folds +) | agg cosine vs STACK | macro | Δ`nmae_wald89` | Δ prec-89 | vs `cd4svd5`: Δ`pds` / agg |
|---|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.0798 | **+0.006399** (3/3) | **+0.122 %** | −0.840 % | −0.00012 | −0.00001 | −0.005864 / +0.608 % |
| 2 | 0.0380 | −0.000148 (2/3) | −0.158 % | −0.070 % | −0.00167 | −0.00113 | −0.012410 / +0.339 % |
| 3 | 0.0327 | +0.001299 (2/3) | −0.232 % | −0.179 % | +0.00123 | −0.00048 | −0.010964 / +0.259 % |
| 4 | 0.0295 | +0.000529 (2/3) | **−0.406 %** | −0.208 % | −0.00181 | +0.00012 | −0.011734 / +0.091 % |
| 5 | 0.0203 | +0.001848 (3/3) | −0.036 % | −0.011 % | −0.00026 | −0.00099 | −0.010415 / +0.451 % |

  Direction 1 carries half the `pds` at no aggregate-cosine cost (its cost is on the macro cosine, i.e. on the
  small rows); the five singles sum to +0.0099 against the joint +0.0123 (super-additive, ×1.24), so no
  single-direction cell approaches `cd4svd5` and all five are −0.0059 … −0.0124 against it. The joint −0.46 % is
  not attributable to any single direction from the singles (they range −0.41 … +0.12 %); §6 measures the
  complement.

### 4.3 Controls
Energy-matched random subspaces (246's `rand`, one seed, resolved sizing) at every cell's energy: vs STACK Δ`pds`
−0.0055 … +0.0003 (the +0.0003 is the 8-dimensional control at `k5_dir4`'s 3 % energy) and aggregate cosine
−0.30 … −3.63 % (the `k5_s1` one = 246/249's −0.001377 / −2.82 %). At every
s and k the cells' `pds` gain is direction-specific and their cosine behaviour is the opposite of a norm cut.
Leave-top-3-out retention (246's statistic, per-target Δ`pds`): `k5_s0.75` vs STACK **61.9 %** (50 / 58 / 72 % per
fold), `cd4svd5` 52–56 %; of the best cell's increment over `cd4svd5`: undefined in any useful sense (160 % / −554 % /
138 %: a small diffuse loss, 3-fold mean −0.0018 with the top-3 targets removed −0.0033).

---

## 5. The best cell — `k5_s0.75` (k = 5, s = 0.75; shrink the top-5 component to a quarter instead of removing it)

### 5.1 Increment over `cd4svd5`
Δ`pds` **−0.001829** raw (1/3 folds; −0.00176 / +0.00029 / −0.00401; −1.22 σ), aggregate cosine **+0.482 %**
(+0.798 / +0.884 / −0.235 %), macro +0.616 %, delivered amplitude +0.118 % (shape-only), CD4 fraction ×0.9944, energy
removed 0.188 vs 0.201, Δq_w +0.00003, Δ`nmae_wald89` −0.00011 (identity: a gain of +0.00003 overall), Δ prec-89
+0.00060. `kolf2` apart: Δ`pds` +0.00026, aggregate +1.285 %.

### 5.2 Six members
| factors | `pds` | `mse` (Δcos_board / Δraw) | `reach` | `fid` | `jac` | `nmae` | **overall** | vs +0.0019 |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| point (k_pds 1.22, k_cos 3.76, DE share, identity ×1) | −0.00082 | +0.00067 (+0.001912 / −0.00382) | +0.00026 | +0.00009 | +0.00001 | +0.00003 | **+0.00023** | 0.12× |
| conservative (k_pds 0.85, k_cos 3.76, DE-four whole, no `nmae` credit) | −0.00057 | +0.00067 | +0.00033 | +0.00011 | +0.00001 | 0 | **+0.00055** | **0.29×** |
| k_cos 2.15 (point k_pds) | −0.00082 | +0.00038 | +0.00015 | +0.00005 | +0.00001 | +0.00003 | **−0.00021** | — |
| k_cos 1.36 (point k_pds) | −0.00082 | +0.00024 | +0.00009 | +0.00003 | 0 | +0.00003 | **−0.00042** | — |

The exact-identity `mse` form gives the same totals to 1e−5. The conservative total is larger than the point total
only because the conservative factor SHRINKS the `pds` loss (0.85×) while the cosine legs do not change — an
artefact of pricing a `pds` LOSS with a discounted transfer, so the honest reading of this cell is the point row.

### 5.3 A bound on the whole s axis between 0.75 and 1
Any s in (0.75, 1) recovers a fraction of the +0.482 % cosine and loses a fraction of the −0.00183 `pds`. The most it
can be worth is the full cosine recovery at zero `pds` loss: +0.00067 + 0.00036 = **+0.00103 overall at k_cos 3.76
(0.53× the gate), +0.00059 at 2.15, +0.00037 at 1.36**. Since `pds` is near-linear in s (§4.1) and the cosine loss is
convex in s, no interior s reaches that bound. A finer s grid is not worth a fold load.

### 5.4 Controls for this cell
Energy-matched random subspace at its energy (47 / 52 / 51 dims, E 0.1905 vs 0.1883): vs STACK Δ`pds` −0.000997
(1/3), aggregate −1.930 % (the cell: +0.010433, +0.009 %); vs `cd4svd5` −0.01326 / −1.448 %. Leave-top-3-out
retention vs STACK 61.9 %; of the increment, undefined (§4.3).

---

## 6. EXPLORATORY — post hoc, not part of the fixed grid, not used for the verdict (`sweep257b.py`, `explore257.py`)

The per-direction singles (§4.2) said direction 4 alone costs −0.41 % and buys nothing, and that the directions are
not additive. The complement was measured: the k = 5 projector with ONE direction left in (`k5_drop{j}`, s = 1),
resolved sizing; known-answer rows 8/8 at 0.0e+00.

| cell | left in | E removed | CD4 fraction × | Δ`pds` vs `cd4svd5` (folds +, σ) | rel `cos_agg` (per fold) | macro | delivered | Δ`pds` vs STACK | agg vs STACK | **overall point / conservative** | k_cos 2.15 | `kolf2` |
|---|---|---:|---:|---|---|---:|---:|---:|---:|---|---:|---|
| `k5_drop1` | 1 | 0.121 | 0.9699 | **−0.009217** (0/3, −6.2) | −0.237 % (−0.14 / +0.69 / −1.26) | +0.644 % | +0.403 % | +0.003046 | −0.711 % | **−0.00442 / −0.00331** | −0.00420 | +0.00155 / +2.37 % |
| `k5_drop2` | 2 | 0.163 | 0.9875 | −0.000569 (2/3, −0.4) | +0.114 % (+0.24 / +0.16 / −0.05) | −0.060 % | +0.022 % | +0.011694 | −0.351 % | **−0.00001 / +0.00009** | −0.00011 | −0.00004 / +0.24 % |
| `k5_drop3` | 3 | 0.168 | 0.9883 | −0.001830 (0/3, −1.2) | +0.140 % (+0.44 / +0.38 / −0.40) | +0.031 % | +0.205 % | +0.010433 | −0.329 % | **−0.00041 / −0.00018** | −0.00054 | −0.00378 / +0.33 % |
| **`k5_drop4`** | 4 | 0.171 | 0.9940 | −0.000520 (1/3, −0.3) | **+0.384 % (+0.40 / +0.39 / +0.37, 3/3)** | +0.078 % | +0.019 % | +0.011743 | −0.079 % | **+0.00059 / +0.00074** | +0.00024 | −0.00106 / −0.09 % |
| `k5_drop5` | 5 | 0.180 | 0.9911 | −0.003386 (0/3, −2.3) | −0.106 % (+0.41 / +0.26 / −0.99) | −0.178 % | −0.036 % | +0.008877 | −0.579 % | **−0.00176 / −0.00133** | −0.00167 | −0.00067 / +0.78 % |

Reading: **direction 1 carries 75 % of the `pds`** (dropping it loses −0.0092 of +0.0123), directions 2 and 4 carry
none (−0.0006 / −0.0005, inside noise), directions 3 and 5 carry −0.0018 / −0.0034. **Direction 4 carries a
fold-consistent cosine cost**: leaving it in recovers +0.38 % aggregate on all three folds (the only cell in this
experiment whose cosine increment has one sign on all three folds) at zero `pds` cost — the projector without it
prices +0.00059 point / +0.00074 conservative, **0.39× the gate**, +0.00024 at k_cos 2.15. It is the best-priced
object measured here and it is not a candidate: it is post hoc (one of five, chosen after the read), a `pds`
reading inside PLOLO's noise, and under half the gate at its most favourable factor. It is recorded so the next
mandate on this class starts from "the cost is in direction 4 and the last quarter of the removal", not from zero.
Legs for `k5_drop4` (point): `pds` −0.00023, `mse` +0.00053, DE-three +0.00028, `nmae` +0.00001; Δ`nmae_wald89`
−0.00002, Δ prec-89 −0.00033, Δq_w −0.00018.

---

## 7. The builder flag such a cell would need (named, NOT made)

The builder has only `VCC_CD4_DEFLATE_K` (a full projector, `scripts/make_submission_pm_cd4.py:1993-2065`) and, for
October, `VCC_CD4_DEFLATE_BASIS` (the pinned gene-space basis, `:2177-2246`). A shrink cell would need
**`VCC_CD4_DEFLATE_S`** — a float in (0, 1], unset/1.0 byte-inert — read at import next to `CD4_DEFLATE_K` (`:740`),
recorded as `cd4_deflate_s` in the meta config (`:5302`) and in the `effects.cd4.deflation` audit, and applied by one
line in each path:
- `deflate_cd4_residual` (`:2050`): `R = Rc - Uk @ (Uk.T @ Rc)` → `R = Rc - CD4_DEFLATE_S * (Uk @ (Uk.T @ Rc))`;
- `deflate_cd4_residual_pinned` (`:2226`): `Rk = Rc - proj @ V.T` → `Rk = Rc - CD4_DEFLATE_S * (proj @ V.T)`.
A per-direction cell (§6) would instead need `VCC_CD4_DEFLATE_DIRS` (1-based indices into the top-k, selecting
columns of `Uk` / of the pinned `V`). **Neither is warranted by this result**; nothing was edited.

---

## 8. Honest limits

- **One instrument.** Every increment is 162's PLOLO on three cancer folds; the `pds` losses of the softer cells are
  2–8 σ and unambiguous, but the best cell's increment is −1.2 σ on `pds` and its cosine recovery is fold-asymmetric
  (+0.80 / +0.88 / −0.24 %). The board is the only instrument that could price a +0.0002 … +0.0006 object, and
  §8.1b's ordering noise (sd 0.0045 overall) makes such an object unreadable there in any case.
- **The recipe is a retrodiction of one point.** k_cos 3.76 and k_pds 1.22 come from the `cd4svd5` point alone
  (BR249 §2); R1 passing shows the arithmetic is wired as the board read, not that the factors generalise. The
  sensitivity rows at 2.15 and 1.36 are the honest band, and every cell is negative there.
- **DE-three propagation** is 243's union-point 0.24 per unit board cosine, split by 218's card; it retrodicts
  `cd4svd5`'s DE-three (−0.00035 vs −0.00034) but rests on two board points.
- **Delivered amplitude is a PLOLO proxy** (median row norm of the frame-A render): it read `cd4svd5` at −0.533 % on
  the three cancer folds against the package's measured −0.409 % (six246's 4-fold −0.317 %), so its cells' values
  carry ±0.15 pp; the five cells flagged "NO" for shape-only sit at +0.52–0.59 % and would need the +0.0089 gate,
  which they miss by an order of magnitude regardless.
- **`nmae`** is the exact identity on the wald89 pool, 3 cancer folds; increments are ≤ ±0.0004 overall (dir 3
  −0.0007) — priced as a magnitude with the BR237 caveat, and never the deciding leg here.
- **§6 is post hoc** and reported as such; its +0.00074 is not a price to act on.
- **Not measured:** s > 1 (over-removal), s between 0.75 and 1 (bounded in §5.3), any change to the H1 union
  (246 §3.3: deflating it costs cosine for nothing), October's disjoint panels (250: the basis is pinned; a shrink
  would apply to the pinned path unchanged).

---

## 9. Cost and memory (one fold per process; no other heavy job running)

| script | rows | wall | peak RSS |
|---|---:|---:|---:|
| `sweep257.py --fold k562_gwps` | 39 | 39.4 s | 4.52 GiB |
| `sweep257.py --fold orion_hct116` | 39 | 50.5 s | 4.72 GiB |
| `sweep257.py --fold orion_hek293t` | 39 | 52.5 s | 4.73 GiB |
| `sweep257.py --fold kolf2` | 39 | 47.6 s | 4.72 GiB |
| `sweep257b.py` (4 folds, exploratory) | 7 each | 7.9–10.4 s | 4.44–4.65 GiB |
| `price257.py` / `explore257.py` | — | < 5 s | 0.04 GiB |

Fold caches: `legacy` status on every fold (pre-252 caches validated field by field by `fold_cache_status`).

---

## 10. Files

- `lib257.py` — `deflate_shrink`, `rand_energy_matched` (246's, verbatim), `Identity` (188/six246's `nmae` identity,
  precision, delivered norm), `resolved_parts`; imports `lib246`, pins `plolo.WORK` to `work/plolo/`.
- `sweep257.py` — the grid, one fold per process; `work/sweep257-<fold>.{json,log}`, `pertarget257-<fold>.npz`,
  `identity257-<fold>.npz`.
- `price257.py` — the recipe, R1/R2/R3, the grid price, best cell, controls; `work/price257.{json,log}`.
- `sweep257b.py`, `explore257.py` — §6 (exploratory); `work/sweep257b-*.{json,log}`, `work/explore257.json`.
- `work/plolo/` — byte copies of 246's fold caches (shas in §1.1); `work/.marker-start` — the write-audit marker.
- Nothing outside this directory was created or modified.

---

## Verdict

**DOES NOT CLEAR.** No softer deflation — s ∈ {0.25, 0.5, 0.75} at k ∈ {3, 5, 8}, k = 8 at s = 1, or any of the
top-5 directions removed alone — clears +0.0019 overall over `cd4svd5` at the conservative factors (best `k5_s0.75`
+0.00055, 0.29×; point +0.00023; negative at k_cos 2.15 and 1.36), and the whole s axis between 0.75 and 1 is bounded
at +0.0010 (0.53×). The projector's cosine cost lives in the last quarter of the removal and in direction 4; the
board's factors price that quarter at break-even, so `cd4svd5` is already the right point on this axis. No
preregistration draft, no build, no flag.
