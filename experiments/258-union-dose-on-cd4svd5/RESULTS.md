# 258 — The six-source union on top of `cd4svd5`: **no dose and no source subset clears +0.0019 on the validated `nmae` identity.** The best cell (H1 pair + `dld1_lowmoi` at g = 0.35) prices **+0.00119** at the conservative factors (0.63× the gate), on a `pds` gain of 0.8 LOO-σ that leave-top-3-out turns negative. The six-source union itself is negative at every dose (−0.0027 … −0.0051 conservative; −0.0004 … −0.0029 even at the point factors). **DOES NOT CLEAR — verdict line at the bottom (§10).**

**Date 2026-09-05 (01:55–02:10 UTC).** Deterministic. **Measurement only: no build, no `vcc`, no `cell-eval2`, no
model name, no description, nothing under `submissions/`, nothing written outside this directory** (the PLOLO fold
caches are byte copies of 249's, `work/plolo/SHA256SUMS.txt`; the repo-root `work/` that `git status` shows is
another agent's, born 2026-09-04 07:04). Builder: git HEAD's `scripts/make_submission_pm_cd4.py` (`3615d20f…`),
imported with the incumbent's env and driven through its own `apply_frozen_cd4_residual` with the module globals as
the knobs (248/249's route). Peak RSS per script in §9: 4.96–5.26 GiB per fold, 3.18 GiB for the production frame,
0.25 GiB for the price. One sibling (257) ran a 4.6 GB sweep alongside; nothing ≥ 15 GiB was live.

---

## 0. Verdict

```
R1  stack6 g=0.35 on the h1union base reproduces 229's six229.json               PASS  (0.0e+00 on every key;
    nmae_wald89 1.00291120 vs the mandate's 1.0029112, overall -0.00478294 vs -0.0047829)
R2  ARM6's six-member price on BR237's board point                                6/36 factor sets inside all
    four bands -- every one of them needs k_pds at the ARM6-REALISED value (3.79 frame A / 2.48 frames ABC),
    the nmae identity x1.5, and the amplitude leg; k_pds 0.7 misses the pds band by construction
R3  cd4svd5 vs h1union at k_cos 3.76 / k_pds 1.22                                 PASS  12/12 variants,
    +0.00459 ... +0.00494 vs +0.0043 +- 0.0008 (board +0.0042634)

THE SCREEN (31 cells on top of cd4svd5, three cancer folds):
    cells clearing +0.0019 at the CONSERVATIVE factors (k_pds 0.7, nmae x1.5)     0 / 31
    cells clearing +0.0019 at the POINT factors (k_pds ARM6-realised, nmae x1.0)  2 / 31  (stack3 g=0.35 +0.00343,
                                                                                           stack3 g=0.45 +0.00272)
    best cell: H1 pair + dld1_lowmoi at g = 0.35   conservative +0.00119 (0.63x gate), point +0.00343;
               d pds +0.00119 (3/3 folds, 0.8 LOO-sigma, fold se 0.00095), leave-top-3-out retention -66 %
    the six-source union, every dose 0.10 ... 0.45:  conservative -0.0027 ... -0.0051, point -0.0004 ... -0.0029

DOES NOT CLEAR.  The family is closed on this base: the union's magnitude cost on the significant genes
(the identity's leg, +0.0143 ... +0.0159 raw at g >= 0.28, e = dM/dU -0.02 ... -0.05) is 2.5-3x its
pds + mse gain at every dose where the gain exists, and below g = 0.28 the single shared dose under-doses the
H1 pair the incumbent already carries (cosine -0.8 ... -4.8 %).
```

**Two findings beyond the gate.** (1) **The CD4 deflation already bought a third of what the union's directions
buy on `pds`:** the same six-source union at g = 0.35 reads Δ`pds` +0.002687 (frames ABC) / +0.001755 (frame A)
on the `h1union` base and **+0.001906 / +0.001144 on `cd4svd5`** (−29 % / −35 %) at an unchanged cosine gain
(+3.673 % → +3.686 %). (2) **The union's `nmae` cost is `igvf_cardiomyocyte`'s:** the H1 pair + igvf alone charges
+0.0144 raw at g = 0.28 (e −0.06) for no `pds`; the pair + `dld1_lowmoi` charges +0.0014 with e +0.20 and is the
only subset whose `pds` is positive on 3/3 folds. `jiang_5line_noK562` and `myllia_a549` are `nmae`-cheap but
`pds`-negative. An `nmae`-neutral rendering of the six directions (BR237 §3) would price **+0.0014 conservative /
+0.0034 point** on this base (§6) — under the gate at the conservative factors and needing k_pds ≥ 2.5 to clear.

---

## 1. The instrument (`price258.py`, three phases, every number checked against an artifact on disk)

**Construction.** On each of 162's PLOLO folds the production builder's `apply_frozen_cd4_residual` is called with
the fold's whitened base: the CD4 residual (`cd4_for`, complete-mean3) is **deflated at k = 5 before sizing**
(`deflate_cd4_residual`, 249's flag), sized by `solve_cd4_fraction` at f = 0.34 (the resolved variant, the
package's object), the union is prepared by `prepare_h1_union(sources=…)` from the registry in the env-var order
and sized by the carrying-row rule at g, ONE merged add, ONE global preserve — 249 §3 proved this path equals 246's
`compose` to 0.0e+00, and 248 §4.4 proved it reproduces 229's single-union arms. Each arm is scored on two
functionals from six renders (two tags × three challenge frames):

- **(a) 218's `score_effects` functional** — `cos_agg` (energy-weighted aggregate, a23 tag), macro cosine, `pds`,
  each the mean over frames A, B, C. This is the frame 243's class rule quotes (`rel cos_agg`, three cancer folds,
  mean of per-fold ratios) and the frame 246/249's Δ`pds` +0.012263 lives in.
- **(b) 216/229's frame-A functional** — 188's exact `nmae` identity on three pools (`wald89`, `lfc89`, `h89`; 4
  folds pooled, `kolf2` included, exactly as `six229.py` does), top-89 precision, the delivered norm ratio, and
  `pds` on frame A. This is the frame `six229.json` lives in.

The two `pds` functionals differ (frame A alone vs the A/B/C mean): `six229.json`'s `d_pds_3cancer` for ARM6 is
+0.001755 where `price229.json`'s is +0.002687. Both are carried everywhere, each with its own calibrated `k_pds`.

**The production frame without a 20 GiB staging (`--prod`).** 249 staged the shipped `h1union` effect matrix
(`effects-off.npy`, `abeb0749…`) together with the CD4 residual before deflation (`cd4-off.npz`) and the H1-pair
union (`union-off.npz`); the preserve step is one scalar (`final = s·(wh + merged)`), so the whitened base is
recovered exactly (`wh = final/s − merged`) and any arm can be composed on the real 300-target panel in 44 MB. Its
delivered amplitude is read as 211/243's (`A_c`, `B_c`) terms on the effect-level render (`deliver`, minnorm, the
scorer's own `_row_weights`/`_across_pert_budget`), never the floor.

**Known-answer rows — all pass, most at 0.0e+00.**

| row | against | rows | result |
|---|---|---:|---|
| per fold: `h1union` == 246 `sweep246` STACK (6 keys); `cd4svd5` == 246 `B2_cd4_svd5_resolved` (4 keys + fraction) and 249 `repro249` (3 keys + union fraction); ARM6/ARM3/ARM4 == 229 `price229.json` (`cos_agg`, `pds`, macro, union fraction); frame-A `pds` of six arms == `six229.json` `pds_per_fold` | 246 / 249 / 229 | **33 × 4 folds** | **132/132, all 0.0e+00** |
| production frame: recovered base median == `stage-off`'s 4.229295604454299; `prepare_h1_union(PAIR)` == `union-off.npz` bit for bit; re-composed `h1union` == `effects-off.npy` (max\|Δ\| 5.6e−17); **re-composed `cd4svd5` == 249's `effects-svd5.npy` `a28d5be2…` (max\|Δ\| 2.2e−16), fraction 0.34104860341941284, preserve 0.9115882187626339, energy removed 0.19866046506100588** (255 KA 2 re-hashed the same file); re-composed ARM6 == 229's `effects-stack6g035.npy` `c5721825…` (2.2e−16), union fraction 0.22985837674034879 | 249 / 229 / 255 | **14** | **14/14** |
| cross-fold: `nmae_wald89`/`lfc89`/`h89`, top-89 precision, delivered norm ratio, frame-A `pds_3cancer` for a15, cd4 f=0.34, h1union, stack6, stack4, stack3 == `six229.json`; `cd4svd5` == `six246.json` `B2_cd4_svd5_resolved`; a15 / f=0.34 == `six216.json`; **R1 by the mandate's digits** | 229 / 246 / 216 | **46** | **46/46, 0.0e+00 (R1 4.8e−09 / 3.5e−08 to the rounded digits)** |

---

## 2. The retrodiction rows

### R1 — `stack6` g = 0.35 on the `h1union` base reproduces `six229.json`: PASS

`nmae_wald89` **1.00291120** (six229 1.0029111951943475; mandate 1.0029112), `overall_nmae_vs_incumbent`
**−0.00478294** (mandate −0.0047829), frame-A `pds_3cancer` 0.78830127 (Δ +0.001755), top-89 precision
0.60427826, delivered norm ratio 1.00056456 — every one at 0.0e+00 against the json. In 246's frame: per-fold
`cos_agg` 0.03909261 / 0.03230229 / 0.02569506 and `pds` 0.80779190 / 0.79370401 / 0.76779813 equal
`price229.json`'s ARM6 rows at 0.0e+00 (rel `cos_agg` +3.6728 % = 243's +3.673 %; Δ`pds` +0.002687).

### R2 — ARM6's six-member price vs BR237: lands only with the ARM6-realised `k_pds`, `nmae` ×1.5 and the amplitude leg

Board (BR237): overall −0.0036403; `pds` +0.00248, `mse` +0.00088, `nmae` −0.00716, DE-three +0.00016.
Recipe: `mse` from 197's identity at ARM6's own terms with k_cos 0.96 on `rel cos_agg` +3.6728 % ×
cos_h1union 0.107332; `pds` = 2.2138 × k_pds × Δ`pds`/6; DE-three = rate × Δcos; `nmae` = the identity's −0.004783.

| factor set (36 tried) | `pds` | `mse` | DE-three | `nmae` | TOTAL | bands |
|---|---:|---:|---:|---:|---:|---|
| k_pds 0.7, any DE rate, `nmae` ×1.0 / ×1.5, package amplitude | +0.00045 | +0.00088 | +0.00016 … +0.00091 | −0.00478 / −0.00717 | −0.0025 … −0.0057 | **`pds` MISS** (band +0.0025 ± 0.0008 needs k_pds ≥ 2.6) |
| **k_pds ARM6-realised (3.79 A / 2.48 ABC), `nmae` ×1.5, package or calibrated effect-level amplitude, any DE rate** | **+0.00246** | **+0.00088 / +0.00089** | +0.00091 (0.24) / +0.00052 (0.136) / +0.00016 (realised) | **−0.00717** | **−0.00293 / −0.00332 / −0.00368** | **ALL PASS (6 sets)** |
| k_pds ARM6-realised, `nmae` ×1.0 | +0.00246 | +0.00088 | as above | −0.00478 | −0.0005 … −0.0013 | total MISS (band −0.0036 ± 0.0010); the `nmae` band itself passes by 0.00008 |
| no amplitude leg (any k_pds) | | **+0.00133** | | | | **`mse` MISS** (band +0.0009 ± 0.0003) |

**Which factor set does it:** the ARM6-realised `k_pds`, the identity at ×1.5, and the amplitude leg
(ARM6's −1.24 % delivered amplitude is +0.00250 raw `mse` at zero cosine change, −0.00044 overall). The DE-three
rate is irrelevant to the pass but not to the reading: **the mandate's 0.24 recipe over-predicts ARM6's realised
DE-three (+0.000161) by 5.6×**; 243 §4.3's 0.239 is a DE-FOUR rate that carries `nmae` at 43 %, so its DE-three
share is 0.136 (3.2× over); the realised rate is **0.0425 per unit board Δcos** and is what the conservative
price below carries. The `cd4svd5` step realised 0.185 per unit Δcos (§3) — the DE-three propagation is not one
number across classes, and at ±0.0003 it never decides a cell here.

**Honesty about R2.** k_pds and k_cos in the passing sets are the values BR237 itself fixed (0.962 and the
`pds` ratio), so the `pds` and cosine legs pass by construction. The non-circular content of R2 is (i) the
identity's leg, −0.004783 predicted before the slot against −0.00716 realised (ratio 1.497, the mandate's ×1.5),
(ii) the amplitude leg, without which the `mse` band fails, and (iii) the DE-three rate. R3 is circular the same
way on `pds` and `mse` (k_cos 3.76 / k_pds 1.22 are BR249's own ratios); what it tests is the arithmetic — spans,
the identity at the arm's own terms, the amplitude leg — and the identity's `nmae` on a second class.

### R3 — `cd4svd5` vs `h1union` at k_cos 3.76 / k_pds 1.22: PASS, 12/12 variants

`pds` +0.00552 (board +0.005452), `mse` −0.00079 / −0.00080 with the package / calibrated amplitude (board
−0.000799), DE-three −0.00025 … −0.00045 (board −0.000344), `nmae` +0.00031 / +0.00047 (identity, ×1.0 / ×1.5;
board −0.000045), **TOTAL +0.00459 … +0.00494** against +0.0043 ± 0.0008 (board +0.0042634). The identity at
`cd4svd5`'s own terms reproduces its board raw `mse` 0.8218062 to 0.0e+00, and reads its measured cosine as
**0.105469** (−1.735 % vs `h1union`'s 0.107332, BR249 §2).

---

## 3. The board-calibrated factors (all computed here from the captured entries and the packages on disk)

| quantity | ARM6 (new source into empty rows, BR237) | `cd4svd5` (re-shape of occupied rows, BR249) |
|---|---:|---:|
| measured board cosine (197's identity at own terms) | 0.111123 | 0.105469 |
| PLOLO rel `cos_agg`, 3 cancer folds | +3.6728 % (1.966 / 3.412 / 5.640) | −0.4617 % |
| board Δcos → **k_cos** | +0.003791 → **0.962** | −0.001863 → **3.759** |
| PLOLO Δ`pds` frame A / frames ABC | +0.001755 / +0.002687 | +0.011028 / +0.012263 |
| board raw Δ`pds` → **k_pds** (A / ABC) | +0.006656 → **3.79 / 2.48** | +0.014920 → **1.35 / 1.22** |
| DE-three realised (overall) → per unit Δcos | +0.000161 → **0.0425** | −0.000344 → 0.185 |
| `nmae`: identity → board | −0.004783 → −0.007160, **ratio 1.497** | +0.000314 → −0.000045 (flat, as priced) |

**The delivered-amplitude instrument.** Effect-level (`A_c`, `B_c`) on the production panel against the packages
streamed by 243/249/255: ARM6 vs `h1union` **−1.234 % vs package −1.242 %**; `cd4svd5` vs `h1union` −0.428 % vs
−0.409 %; `cd4svd10` vs `cd4svd5` −0.458 % vs −0.456 % — package/effect-level 1.007 / 0.954 / 0.996, mean
**0.986**, applied to every cell. 216's median-row-norm proxy (the mandate's fallback, "246's delivered ratio")
reads ARM6 at **−1.81 %** — 46 % off — because it is a median over rows and the union moves the energy-weighted
tail; it is reported alongside but not used. Every cell here moves ≤ 1.1 %, so the amplitude leg is ≤ 0.0002
overall and never decides a cell (it is charged wherever it is, > 0.2 % or not).

---

## 4. Dose — on top of `cd4svd5`, three cancer folds

Columns: carry = union carrying rows on the hek293t fold (/294; the production panel carries 176 for six sources,
125 for three, 77 for the pair); agg / mac = rel `cos_agg` / macro cosine vs `cd4svd5`; Δ`pds` A / ABC = the two
functionals (per-fold ABC values in brackets, k562 / hct116 / hek293t); + = folds with Δ`pds` > 0; σ = \|Δ`pds`
ABC\| / 0.001494; Δ`nmae` = the exact identity's raw shift (4 folds; the 3-cancer-fold value differs by < 20 % on
every cell, same sign); e = 188's ΔM/ΔU (break-even 0.5); amp = calibrated delivered amplitude; **CONS** = k_pds
0.7 (frame A), `mse` at k_cos 0.96 with the amplitude leg, DE-three 0.0425, `nmae` ×1.5; **POINT** = k_pds 3.79
(frame A), DE-three 0.24, `nmae` ×1.0; ret3 = leave-top-3-out `pds` retention (222's definition).

### 4.1 The six-source union replacing the H1 pair (229's construction)

| g | carry | agg % | mac % | Δ`pds` A | Δ`pds` ABC (k562 / hct / hek) | + | σ | Δ`nmae` | e | amp % | `pds`(0.7) | `mse` | DE3 | `nmae`×1.5 | **CONS** | `pds`(3.79) | **POINT** | ret3 |
|---:|---:|---:|---:|---:|---|:-:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | 172 | −4.768 | −2.975 | +0.00046 | +0.00065 (−0.0036 / +0.0027 / +0.0029) | 2 | 0.4 | +0.00236 | −0.33 | −0.09 | +0.00012 | −0.00172 | −0.00021 | −0.00098 | **−0.00279** | +0.00065 | −0.00288 | −274 % |
| 0.15 | 172 | −2.700 | −1.305 | +0.00134 | +0.00186 (−0.0016 / +0.0036 / +0.0036) | 2 | 1.2 | +0.00479 | −0.04 | +0.22 | +0.00035 | −0.00088 | −0.00012 | −0.00200 | **−0.00265** | +0.00187 | −0.00100 | −42 % |
| 0.20 | 172 | −0.823 | +0.091 | +0.00165 | +0.00230 (−0.0004 / +0.0036 / +0.0037) | 2 | 1.5 | +0.00761 | −0.02 | −0.20 | +0.00043 | −0.00036 | −0.00004 | −0.00317 | **−0.00313** | +0.00231 | −0.00036 | −26 % |
| 0.28 | 172 | +1.791 | +1.844 | +0.00164 | +0.00243 (+0.0005 / +0.0028 / +0.0039) | 3 | 1.6 | +0.01436 | −0.02 | −0.71 | +0.00042 | +0.00039 | +0.00008 | −0.00598 | **−0.00509** | +0.00229 | −0.00087 | −38 % |
| 0.35 | 172 | +3.686 | +2.955 | +0.00114 | +0.00191 (+0.0007 / +0.0016 / +0.0034) | 3 | 1.3 | +0.01565 | −0.03 | −1.05 | +0.00030 | +0.00093 | +0.00016 | −0.00652 | **−0.00513** | +0.00160 | −0.00092 | −98 % |
| 0.45 | 172 | +5.835 | +4.031 | +0.00001 | +0.00073 (−0.0001 / +0.0002 / +0.0021) | 2 | 0.5 | +0.01585 | −0.05 | −0.31 | +0.00000 | +0.00195 | +0.00025 | −0.00660 | **−0.00440** | +0.00002 | −0.00101 | −519 % |

**Negative at every dose, at both factor sets.** The mechanism has two halves and the single shared dose cannot
escape either: below g = 0.28 the H1 pair the incumbent carries at 0.28 is under-dosed (cosine −0.8 … −4.8 %,
`k562_gwps` `pds` negative), and from g = 0.28 up the four magnitude-heavy sources push the identity to
+0.0144 … +0.0159 raw with e ≈ −0.02 … −0.05 — "magnitudes, not signs" (188), the same failure ARM6 showed on
the old base (+0.0172, e −0.033) and the board confirmed at 150 %. The `pds` gain peaks at g = 0.28 (+0.00243 ABC,
3/3, 1.6 σ) and is already smaller than on the `h1union` base (§0, finding 1).

### 4.2 The H1 pair + `dld1_lowmoi` (229's ARM3 / 248's DLD-1 rider, as a single union)

| g | carry | agg % | mac % | Δ`pds` A | Δ`pds` ABC (k562 / hct / hek) | + | σ | Δ`nmae` | e | amp % | `pds`(0.7) | `mse` | DE3 | `nmae`×1.5 | **CONS** | `pds`(3.79) | **POINT** | ret3 |
|---:|---:|---:|---:|---:|---|:-:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | 123 | −5.498 | −3.858 | −0.00039 | −0.00029 (−0.0033 / +0.0010 / +0.0015) | 2 | 0.2 | +0.00085 | −1.21 | +0.02 | −0.00010 | −0.00194 | −0.00024 | −0.00036 | **−0.00263** | −0.00054 | −0.00405 | — |
| 0.15 | 123 | −3.686 | −2.436 | +0.00064 | +0.00072 (−0.0010 / +0.0013 / +0.0018) | 2 | 0.5 | +0.00152 | −0.16 | +0.08 | +0.00016 | −0.00128 | −0.00016 | −0.00063 | **−0.00190** | +0.00089 | −0.00170 | −123 % |
| 0.20 | 123 | −1.994 | −1.201 | +0.00142 | +0.00135 (+0.0009 / +0.0016 / +0.0016) | 3 | 0.9 | +0.00133 | +0.11 | +0.32 | +0.00037 | −0.00060 | −0.00009 | −0.00056 | **−0.00087** | +0.00199 | +0.00054 | −15 % |
| 0.28 | 123 | +0.453 | +0.408 | +0.00176 | +0.00158 (+0.0024 / +0.0011 / +0.0012) | 3 | 1.1 | +0.00144 | +0.20 | +0.14 | +0.00045 | +0.00021 | +0.00002 | −0.00060 | **+0.00008** | +0.00246 | +0.00238 | −13 % |
| **0.35** | 123 | +2.337 | +1.485 | +0.00149 | +0.00119 (+0.0031 / +0.0002 / +0.0003) | 3 | 0.8 | +0.00057 | +0.33 | +0.32 | +0.00038 | +0.00094 | +0.00010 | −0.00024 | **+0.00119** | +0.00208 | **+0.00343** | −66 % |
| 0.45 | 123 | +4.613 | +2.571 | +0.00059 | +0.00013 (−0.0013 / +0.0004 / +0.0012) | 1 | 0.1 | +0.00265 | +0.23 | −0.37 | +0.00015 | +0.00150 | +0.00020 | −0.00110 | **+0.00075** | +0.00083 | +0.00272 | — |
| 0.60 | 123 | +7.171 | +3.452 | −0.00140 | −0.00204 (−0.0044 / −0.0007 / −0.0010) | 1 | 1.4 | +0.00437 | +0.14 | −1.99 | −0.00036 | +0.00181 | +0.00031 | −0.00182 | **−0.00006** | −0.00196 | +0.00038 | — |

### 4.3 The H1 pair alone, re-dosed (the incumbent's own union; g = 0.28 is `cd4svd5` itself)

| g | agg % | mac % | Δ`pds` A | Δ`pds` ABC | + | Δ`nmae` | e | amp % | **CONS** | **POINT** |
|---:|---:|---:|---:|---:|:-:|---:|---:|---:|---:|---:|
| 0.20 | −2.513 | −1.659 | +0.00011 | −0.00002 | 2 | +0.00018 | +0.35 | −0.05 | −0.00106 | −0.00141 |
| 0.35 | +2.065 | +1.201 | −0.00027 | −0.00035 | 1 | +0.00060 | +0.26 | −0.15 | +0.00045 | +0.00064 |
| 0.45 | +4.794 | +2.559 | −0.00094 | −0.00110 | 1 | +0.00096 | +0.21 | −1.01 | +0.00090 | +0.00092 |

Most of the pair + DLD-1 cell's cosine is the pair's own re-dose (+2.07 % of the +2.34 % at g = 0.35); what
`dld1_lowmoi` adds is +0.27 % cosine and the whole of the `pds` (+0.0015 ABC), and that `pds` is one fold's
(§6). The pair's own re-dose costs `pds` (−0.0004 … −0.0011, 1/3 folds) — 222's g = 0.28 was not under-dosed on
this channel.

---

## 5. Subsets at g = 0.20 / 0.28 / 0.35 (the four other sources, one at a time on the pair; the four without the pair; five)

| cell | g | carry | agg % | Δ`pds` A | Δ`pds` ABC (k562 / hct / hek) | + | Δ`nmae` | e | **CONS** | **POINT** | ret3 |
|---|---:|---:|---:|---:|---|:-:|---:|---:|---:|---:|---:|
| pair + `jiang_5line_noK562` | 0.20 | 81 | −2.262 | +0.00020 | +0.00023 (−0.0015 / +0.0030 / −0.0008) | 1 | +0.00029 | +0.85 | −0.00105 | −0.00124 | — |
| | 0.28 | 81 | +0.307 | −0.00026 | −0.00024 (−0.0009 / +0.0021 / −0.0019) | 1 | +0.00008 | −0.47 | −0.00000 | −0.00022 | — |
| | 0.35 | 81 | +2.398 | −0.00093 | −0.00099 (−0.0010 / +0.0011 / −0.0031) | 1 | +0.00075 | +0.14 | +0.00033 | −0.00014 | — |
| pair + `igvf_cardiomyocyte` | 0.20 | 138 | −2.083 | +0.00002 | +0.00025 (−0.0029 / +0.0004 / +0.0032) | 2 | **+0.00859** | −0.09 | −0.00439 | −0.00359 | — |
| | 0.28 | 138 | +0.204 | −0.00009 | +0.00048 (−0.0023 / −0.0001 / +0.0038) | 1 | **+0.01441** | −0.06 | −0.00583 | −0.00389 | — |
| | 0.35 | 138 | +1.913 | −0.00058 | +0.00007 (−0.0022 / −0.0009 / +0.0033) | 1 | **+0.01625** | −0.08 | −0.00632 | −0.00435 | — |
| pair + `myllia_a549` | 0.20 | 83 | −2.197 | −0.00024 | −0.00036 (−0.0006 / +0.0001 / −0.0006) | 1 | +0.00016 | +0.31 | −0.00103 | −0.00172 | — |
| | 0.28 | 83 | +0.423 | −0.00042 | −0.00042 (+0.0007 / −0.0005 / −0.0015) | 1 | −0.00005 | +0.32 | +0.00007 | −0.00032 | — |
| | 0.35 | 83 | +2.569 | −0.00079 | −0.00086 (+0.0013 / −0.0013 / −0.0026) | 1 | +0.00061 | +0.24 | +0.00045 | +0.00015 | — |
| `stack4` (the four, H1 pair DROPPED) | 0.20 | 135 | −5.038 | −0.00031 | +0.00037 (−0.0046 / +0.0024 / +0.0033) | 2 | +0.00818 | −0.12 | −0.00555 | −0.00578 | — |
| | 0.28 | 135 | −3.903 | −0.00037 | +0.00031 (−0.0046 / +0.0020 / +0.0035) | 2 | +0.01291 | −0.10 | −0.00712 | −0.00653 | — |
| | 0.35 | 135 | −3.188 | −0.00068 | −0.00004 (−0.0048 / +0.0012 / +0.0035) | 2 | +0.01647 | −0.13 | −0.00862 | −0.00775 | — |
| `stack5` (pair + dld1 + jiang + igvf) | 0.20 | 168 | −1.325 | +0.00194 | +0.00257 (−0.0009 / +0.0039 / +0.0047) | 2 | +0.00770 | −0.03 | −0.00323 | −0.00021 | −13 % |
| | 0.28 | 168 | +1.081 | +0.00198 | +0.00279 (−0.0001 / +0.0032 / +0.0052) | 2 | +0.01440 | −0.03 | −0.00516 | −0.00068 | −21 % |
| | 0.35 | 168 | +2.797 | +0.00151 | +0.00228 (−0.0001 / +0.0021 / +0.0048) | 2 | +0.01570 | −0.04 | −0.00529 | −0.00082 | −66 % |
| `stack3` (pair + dld1) — §4.2 | 0.28 / 0.35 | 123 | +0.45 / +2.34 | +0.00176 / +0.00149 | +0.00158 / +0.00119 | 3 / 3 | +0.00144 / +0.00057 | +0.20 / +0.33 | +0.00008 / **+0.00119** | +0.00238 / **+0.00343** | −13 % / −66 % |
| `stack6` (all six) — §4.1 | 0.28 / 0.35 | 172 | +1.79 / +3.69 | +0.00164 / +0.00114 | +0.00243 / +0.00191 | 3 / 3 | +0.01436 / +0.01565 | −0.02 / −0.03 | −0.00509 / −0.00513 | −0.00087 / −0.00092 | −38 % / −98 % |

**Source attribution.** The identity's cost is `igvf_cardiomyocyte`'s (+0.0144 raw alone at g = 0.28, e −0.06,
`pds` +0.0005 ABC / −0.0001 A) — the four-source and six-source unions inherit it (+0.0129 / +0.0144) and `stack5`
does too (+0.0144); `jiang` and `myllia` are `nmae`-cheap (+0.0001 / −0.0001) but `pds`-negative on 2/3 folds;
dropping the H1 pair (`stack4`) loses 3–5 % of cosine at every dose, 229 §7's finding again. The `pds` the union
buys is `dld1_lowmoi`'s on `k562_gwps` and igvf's on `orion_hek293t` (+0.0032 … +0.0038 there), and the two never
come without igvf's magnitudes or without each other's dilution. `kolf2`, reported apart: every cell with
g < 0.28 loses 8–24 % of `kolf2` cosine (the pluripotent pair under-dosed); at g ≥ 0.35 the six-source union gains
+4.1 %, the pair + dld1 +2.6 %.

---

## 6. The best cell in full — H1 pair + `dld1_lowmoi` at g = 0.35 (`VCC_UNION_SOURCES=gse295214_h1,h1_2025,dld1_lowmoi`, `VCC_H1_UNION_G=0.35`)

| quantity | value |
|---|---|
| PLOLO vs `cd4svd5` | rel `cos_agg` **+2.3366 %** (2.373 / 2.242 / 2.395), rel macro +1.485 % |
| Δ`pds` frames ABC | **+0.001193** (k562 +0.00310, hct116 +0.00017, hek293t +0.00031; 3/3 positive; **0.8 LOO-σ**; fold se 0.00095); frame A +0.001488 |
| `nmae` identity | Δ`wald89` **+0.000567** (3 cancer folds +0.000726); ΔU +141.96, ΔM +46.82, **e 0.330** (magnitudes, but 5× less so than the six-source cells); top-89 precision +0.00133; Δq_w −0.00000 |
| delivered amplitude vs `cd4svd5` | effect-level +0.323 %, calibrated +0.318 % (PLOLO median-norm +0.475 %); charged: +0.00011 overall on the `mse` leg |
| production-frame composition | union 125/300 carrying (41.7 % coverage), solved fraction 0.2351252249384753, CD4 fraction 0.34104860341941284 (untouched), preserve 0.8969284728276888 (`cd4svd5` 0.9116), deflation energy 0.198660 (untouched) |
| `kolf2` (apart) | rel `cos_agg` +2.594 %, Δ`pds` −0.00303 |
| **leave-top-k-out `pds` retention** | k=1 **37.4 %** (top-1 share 212 %), **k=3 −65.9 %** (top-3 share 549 %), k=5 −146.6 % — the k562 fold's TRIM35 +0.228, LENG1 +0.214, MAP2K2 +0.208 carry it (248 saw ZNF324 / LENG1 / TRIM35 for the same source); hct116's SH2D3A / CYTH2 / TTBK2 and hek293t's ZNF219 / LPCAT3 / CYTH2 are offset by MAPK7 −0.376 and ZNF821 −0.251 |

| factor set | Δcos | raw `mse` | `mse` | `pds` A | `pds` ABC | DE3 | `nmae` | **TOTAL A** | TOTAL ABC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **CONSERVATIVE**: k_cos 0.96, k_pds 0.7, DE3 0.0425, `nmae` ×1.5 | +0.002366 | 0.81646 | +0.00094 | +0.00038 | +0.00031 | +0.00010 | −0.00024 | **+0.00119** | +0.00111 |
| conservative, `nmae` ×1.0 | | | +0.00094 | +0.00038 | +0.00031 | +0.00010 | −0.00016 | +0.00127 | +0.00119 |
| k_cos 0.94 (band low), otherwise conservative | +0.002317 | 0.81656 | +0.00092 | +0.00038 | +0.00031 | +0.00010 | −0.00024 | +0.00117 | +0.00109 |
| k_pds 1.22 (BR249's re-shape value), DE3 0.136, `nmae` ×1.0 | | | +0.00094 | +0.00067 | +0.00054 | +0.00032 | −0.00016 | +0.00177 | +0.00164 |
| 243's class rule as 248 priced CUM1 (k_cos 0.94, k_pds 0.70, DE-FOUR 0.239) minus the identity's `nmae` | +0.002317 | 0.81656 | +0.00092 | +0.00038 | +0.00031 | +0.00055 | −0.00016 | ≈ +0.0017 | ≈ +0.0016 |
| **POINT**: k_pds 3.79 A / 2.48 ABC, DE3 0.24, `nmae` ×1.0 | +0.002366 | 0.81646 | +0.00094 | +0.00208 | +0.00109 | +0.00057 | −0.00016 | **+0.00343** | +0.00244 |
| point, DE3 at ARM6's realised 0.0425 | | | +0.00094 | +0.00208 | +0.00109 | +0.00010 | −0.00016 | +0.00296 | +0.00197 |
| k_cos 1.36 (re-dose class), otherwise point | +0.003352 | 0.81449 | +0.00128 | +0.00208 | +0.00109 | +0.00080 | −0.00016 | +0.00401 | +0.00302 |
| `nmae` at ZERO (an `nmae`-neutral rendering, BR237 §3 — not this builder), point otherwise | | | +0.00094 | +0.00208 | +0.00109 | +0.00057 | 0 | +0.00359 | +0.00260 |
| `nmae` at ZERO, conservative otherwise | | | +0.00094 | +0.00038 | +0.00031 | +0.00010 | 0 | +0.00142 | +0.00135 |

**Reading.** The whole distance between +0.00119 and +0.00343 is `k_pds` (0.7 → 3.79 on a +0.00149 frame-A
gain). That point value is ARM6's, calibrated on a six-source gain that was broad on the `h1union` base; here it
is being applied to a 0.8-σ gain 549 % of which sits on three `k562_gwps` targets and which is net negative on
the other 287. Under the conditions 248 used for the same source (price ≥ 1.5× the gate AND leave-top-3-out
retention ≥ 50 %) this cell fails both at the conservative factors and the retention condition at any factor; and
+0.0012 against the ordering channel's sd 0.0045 (230) is a quarter of the noise. **It is not a candidate.**

**Re-pricing the two old-base numbers the mandate asked for.** 229's ARM3 (+0.00067 composite on `h1union`) and
248's DLD-1 rider CUM1 (+0.00216 (i) / +0.00107 (ii) on `h1union`, retention 25.2 %) become, on `cd4svd5`,
**+0.00119 conservative / +0.00343 point / ≈ +0.0017 under 243's class rule** with retention −66 %: the cosine
increment holds (+2.27 % → +2.34 %) but the `pds` increment halves (+0.00256 → +0.00119) — the deflation already
took part of DLD-1's `pds` directions as well as the union's (§0).

**And the six-source union's own ceiling on this base** (the BR237 §3 hypothesis, the identity's leg set to
zero, g = 0.35): conservative +0.00030 + 0.00093 + 0.00016 = **+0.00139**, point +0.00160 + 0.00093 + 0.00090 =
**+0.00343**. An `nmae`-neutral rendering of the six directions would need k_pds ≥ 2.5 to clear on `cd4svd5`
even if it removed the whole −0.0065 leg at no cost — the deflation has already bought a third of what those
directions had to give.

---

## 7. What the screen settles

1. **No dose of the six-source union clears on `cd4svd5`, at either factor set, on either `pds` functional.** The
   family "a single union at a single dose in place of the H1 pair" is closed on this base as it was on the old
   one (229 §0), and for the same reason the board then confirmed (BR237): the identity's leg, now at ×1.5.
2. **The only source whose magnitudes the identity tolerates is `dld1_lowmoi`** (e +0.11 … +0.33 on the pair +
   dld1 cells against −0.02 … −0.13 for anything containing igvf). Its `pds` is real on `k562_gwps` and on three
   targets; on the other two folds it is +0.0002 / +0.0003 and the top-3-out mean is negative.
3. **The union's `pds` gain shrinks by ~30 % on the deflated base** at unchanged cosine (six sources: +0.002687 →
   +0.001906 ABC; dld1: +0.00256 → +0.00119). The deflation's five directions and the union's carrying rows are
   not independent on the `pds` channel; any re-pricing of a direction rider on `cd4svd5` should expect this.
4. **The two `pds` instruments disagree by 1.5× in `k_pds`** (ARM6: 3.79 on frame A vs 2.48 on the A/B/C mean;
   `cd4svd5`: 1.35 vs 1.22) and the class factor itself is 2–3× between the two board points — `k_pds` is not a
   constant of the instrument; it is a property of the arm class (new rows vs re-shaped rows). The conservative
   0.7 is what `h1union` realised (243 §4.3).
5. **The DE-three propagation is not 0.24 per unit Δcos for a union arm**: ARM6 realised 0.0425 (5.6× under the
   mandate's recipe, 3.2× under 243's DE-three share), `cd4svd5` 0.185. At the cosines in play (Δcos ≤ 0.0034) the
   difference is ≤ 0.0007 and never flips a cell; it is charged at the realised rate in the conservative price.
6. **The delivered-amplitude fallback the mandate named (216's median-norm ratio) mis-reads a union arm by ~50 %**
   (ARM6 −1.81 % vs −1.24 %); the effect-level (`A_c`, `B_c`) read on the reconstructed production frame reproduces
   all three package amplitudes within 0.02 pp and costs 30 s. The reconstruction itself (`wh = final/s − merged`
   off 249's staged artifacts) reproduces `a28d5be2…` and `c5721825…` to 2.2e−16 and is the cheap way to put any
   residual arm on the production panel without a 20 GiB stage.

---

## 8. The builder configuration the best cell would need — recorded, NOT recommended, NOT built

All existing flags; no new flag (the group-wise dose 229 §11 sketched is not this cell):

```
VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_CD4_DEFLATE_K=5 \
VCC_H1_UNION_G=0.35 VCC_UNION_SOURCES=gse295214_h1,h1_2025,dld1_lowmoi \
  --whitener data/donors/whitener_union_pruned.npz --whitener-id union_pruned --min-cpm 10 --tail-spread 1.0 \
  --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0 --emitter rc_callset
```

Rule-10 diff it would show against `submissions/249-cd4svd5-f034.h5ad.meta.json`: config `h1_union_g` 0.28 → 0.35
and `h1_union_sources` + `dld1_lowmoi` (its four shas from the registry), `effects.h1_union.*` (carrying 77 → 125,
solved fraction 0.19207161964444322 → 0.2351252249384753), the merged preserve pair (0.9115882187626339 →
0.8969284728276888) and `effects.final_sha256`; `effects.cd4.*` including the deflation block byte-identical. On a
non-A/B/C bundle the driver would carry `VCC_CD4_DEFLATE_BASIS=V5_abc` in place of `VCC_CD4_DEFLATE_K=5` (256's
`DEF_ENV`) and `dld1_lowmoi` is the one registered source that gains coverage in October (67 → 73, 229 §10) — the
cell is not worth a slot today, and that is the only sentence about October this experiment can support.

---

## 9. Cost and memory (every phase alone in its process; nothing ≥ 15 GiB was live; 257's 4.6 GB sweep ran alongside)

| script / phase | runtime | peak RSS |
|---|---:|---:|
| `price258.py --fold k562_gwps` (269 targets, 38 arms) | 38 s | 4.96 GiB |
| `--fold orion_hct116` (290) | 49 s | 5.25 GiB |
| `--fold orion_hek293t` (294) | 50 s | 5.26 GiB |
| `--fold kolf2` (273) | 46 s | 5.19 GiB |
| `--prod` (production-frame reconstruction, 40 compositions, 3 renders each) | 29 s | 3.18 GiB |
| `--price` | 8 s | 0.25 GiB |

The first pass (without the g = 0.45 / 0.60 rows) is kept under `work/pass1/`; every number common to both passes
is identical.

---

## 10. Honest limits

1. **R2 and R3 are consistency checks, not independent retrodictions** on their `pds` and `mse` legs (§2). The
   independent content is the identity's ×1.5, the amplitude leg, and the DE-three rate. The screen's negative
   does not rest on either row: it rests on the identity's leg (validated at scale, conservative) exceeding the
   direction legs by 2.5–3× at every six-source dose, and on a 0.8-σ, three-target `pds` gain for the best cell.
2. **The point `k_pds` (3.79 / 2.48) is a one-point class estimate** from a broad six-source gain; its transfer to a
   concentrated single-source gain is unmeasured. The verdict is taken at the conservative 0.7, as the mandate asks;
   at the point value the best cell reads +0.0034 and would still fail 248's retention condition.
3. **The `nmae` identity is a proxy** (fixed top-89 pools, four folds including `kolf2`; 237 ADDENDUM-2 §2). The
   3-cancer-fold variant is within 20 % and same-signed on every cell (+0.000726 vs +0.000567 for the best cell).
4. **The amplitude read is effect-level**, calibrated on three package pairs (0.954–1.007); the arms here move
   ≤ 1.1 %, so the leg is ≤ 0.0002 overall — the calibration cannot change a verdict.
5. **`fid` / `jac` / `reach`'s yield halves are unreadable on any harness** (rule 9); their board propagation is
   carried at ARM6's realised rate, ±0.0003 at these cosines.
6. **No package, no board member measured.** 248's post-build falsifiers (rule-10 diff, 900-block profile split,
   `n_pred`, `vcc prep --dry-run`) were not reached because nothing clears.
7. The subset grid was {0.20, 0.28, 0.35}; the six-source and pair + dld1 ladders were extended to 0.45 (and 0.60)
   after the first pass put the best cell at the grid's top — the extension only lowered it.

## 11. Files

| file | what |
|---|---|
| `price258.py` | the instrument: `--fold` (both functionals through the builder's path, 33 KA rows per fold), `--prod` (production-frame reconstruction, 14 KA rows, amplitude calibration), `--price` (46 cross-fold KA rows, R1–R3, the tables, the verdict) |
| `run258.sh` | the four folds one process at a time, gated (no ≥ 15 GiB job, ≥ 40 % free) |
| `work/fold-<held>.{json,npz,stdout}` | per-fold rows (38 arms), per-target `pds`/cosine, the identity's U/M/D arrays |
| `work/prod258.{json,stdout}` | the production-frame terms per arm, calibration, KA rows |
| `work/price258.{json,log}` | every number in this file; `work/pass1/` the first pass |
| `work/plolo/` | 249's fold caches, byte copies (`SHA256SUMS.txt`) |

---

**DOES NOT CLEAR.** No dose of the six-source union and no subset of its sources clears +0.0019 over `cd4svd5` at
the conservative factors (k_pds 0.7, the validated `nmae` identity at ×1.5); the six-source union is negative at
every dose at the point factors too; the best cell (H1 pair + `dld1_lowmoi`, g = 0.35) is +0.00119 on a 0.8-σ,
three-target `pds` gain. No `PREREGISTRATION-DRAFT.md`. No model name, no description, nothing submitted.
