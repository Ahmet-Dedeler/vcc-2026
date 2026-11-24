# 261 — Deflating the H1-union residual's target-shared directions before its carrying-row sizing does NOT transfer the CD4 mechanism. The union has no hub to remove: its centred spectrum is flat (leading direction 5.2 % of the energy vs CD4's 8.0 %), every projector costs 5× the cosine per unit energy that CD4's did and buys one-tenth the `pds` (+0.0005 … +0.0015, 0.4–1.0 LOO-σ, 3 targets per fold), and 0 of 10 cells price positive even OPTIMISTICALLY beyond +0.00014. Best cell (centring only) **−0.00024 pessimistic, −0.13× the gate**. **DOES NOT CLEAR — verdict line at the bottom (§10). The analogue is closed.**

**Date 2026-09-05 (03:11–03:35 UTC).** Deterministic. **Measurement only: no build, no `vcc`, no `cell-eval2`, no
model name, no description, nothing under `submissions/`, nothing written outside this directory** (`find -newer
work/.marker-start` over `experiments/ data/ scripts/ docs/ submissions/ src/` returns only the sibling 260's own log;
no source file changed; the builder was imported with `PYTHONDONTWRITEBYTECODE=1`). Fold caches are byte copies of
257's (= 249's = 246's; `work/plolo/SHA256SUMS.txt`, shas `c13c1cec…`, `9f3e0452…`, `f9a3cc30…`, `d86202a9…`). One
PLOLO fold per process, serial; peak RSS **4.91–5.00 GiB** per fold, 1.22 GiB for the production frame, 0.04 GiB for
the price (§9). The sibling 260 ran alongside at **5.9 GB** RSS (the mandate expected < 1 GB); the machine stayed at
73–74 % memory free throughout and nothing ≥ 15 GiB was live.

| row | required | measured | verdict |
|---|---|---|---|
| **R1** k = 0 reproduces `cd4svd5` | 0.0e+00 on every fold and on the production frame | **0.0e+00 on 4/4 folds** (endpoints, q_w, identity, precision, delivered norm, union fraction; per-target `pds` arrays bit-identical); `cd4svd5` == 257 `k5_s1` == 258 / 259 `cd4svd5 <- INCUMBENT` at 0.0 on every key (**44/44 per fold**); production frame: 249's `effects-svd5.npy` (`a28d5be2…`) re-composed at **2.2e−16** through the builder's path and **≤ 1e−9** through this experiment's compose, fractions / preserve at 0.0 (**37/37**) | **PASS** |
| **R2** the price recipe retrodicts BR249 and 257's best cell | overall +0.0043 ± 0.0008; `pds` +0.0055; `mse` −0.0008; 257's `k5_s0.75` | **+0.00472** (board +0.0042634); `pds` **+0.00552** (board +0.00545); `mse` **−0.00083** (board −0.00080); DE-three −0.00035 (board −0.00034); `nmae` +0.00038 (board −0.00005) — equal to 257's R1 at **0.0e+00**; `k5_s0.75` point **+0.00023282** / conservative **+0.00055389** = 257's at 0.0e+00; this file's `price()` == `price257.price()` at 0.0e+00 | **PASS** |
| **R3** the same `deflate` code on CD4 at k = 5 | 246 `B2_cd4_svd5_resolved`: Δ`pds` +0.012263, agg −0.4617 % | **+0.012262757 / −0.461748 % / macro −0.808675 %**, max\|Δ\| **0.0e+00** vs `price246.json` and 249's `rel_agg_plolo`; 3/3 folds, 8.2 σ; energy removed 0.2008; the deflated CD4 bit-identical to `lib246.deflate("svd", 5)` on every fold | **PASS** |
| extra KA: 258's `pair g=0.2` / `pair g=0.35` (the builder's path) == my dose-reference rows | 0.0 | **56/56 keys at 0.0e+00** over 4 folds (endpoints, q_w, precision, delivered norm, union fraction) | **PASS** |
| identity known-answer rows (`six246.json`, 4 folds) | 5e−6 | **9/9 at 0.0e+00** (`a15`, STACK, `cd4svd5`) | **PASS** |
| amplitude calibration factor (package over effect-level A_c on ARM6 / `cd4svd5` / `cd4svd10`) | 258's 0.985550749393 | **0.985550749393**, 0.0e+00 | **PASS** |
| **the gate** | PESSIMISTIC ≥ +0.0019 AND leave-top-3-out ≥ 50 % AND `pds` > 0 on 3/3 folds at ≥ 2 σ AND controls < ½ arm | **0 / 10 cells on any one condition but the fold-sign count.** Best pessimistic **−0.00024** (`u_c0`), best point −0.00013, best optimistic +0.00014; every k ≥ 2 negative at every scenario | **DOES NOT CLEAR** |

---

## 0. The answer in ten lines

1. **The union residual has no target-shared hub.** Its centred carried-row spectrum is flat: the leading direction
   carries **5.2–5.7 %** of the centred energy on every fold (production panel 5.24 %), the next four 4.6 / 4.2 / 3.8 /
   3.5 %; CD4's (257) is **8.0** / 3.8 / 3.3 / 3.0 / 2.0 %. The top five sum to 21.4 % (CD4 20.0 %) but nothing stands
   out, and 257 showed CD4's `pds` lives in its one dominant direction (75 %). The union's rank-1 mean is larger
   relative to its energy (1.75 % of raw vs CD4's 0.07 %), and removing it alone buys **+0.00053** `pds` (0.35 σ) for
   −0.17 % cosine (§3).
2. **Every projector removes transfer signal, not a hub.** Per unit of energy removed the union deflation costs
   **−10.9 … −13.1 % aggregate cosine** where CD4's cost −2.3 % (5× more), and buys **0.004–0.016 `pds`** where CD4's
   bought 0.061 (4–15× less). The exchange rate against `cd4svd5` is 0.0004–0.0013 `pds` per −1 % aggregate cosine
   against `cd4svd5`'s own 0.0266 (§4).
3. **The grid** (k ∈ {1, 2, 3, 5, 8} at g = 0.28): Δ`pds` **+0.0011 / +0.0015 / +0.0011 / +0.0011 / +0.0013** (3/3 folds
   each, but 0.7–1.0 LOO-σ; fold se 0.00044) for aggregate cosine **−0.85 / −1.25 / −1.82 / −3.00 / −3.36 %**. Priced:
   PESSIMISTIC −0.0017 / −0.0025 / −0.0040 / −0.0072 / −0.0081; POINT −0.0013 … −0.0070; OPTIMISTIC +0.00004 / −0.0001 /
   −0.0006 / −0.0015 / −0.0017. **0/10 real cells clear at any scenario** (§3, §5).
4. **Re-dosing does not rescue it.** At g = 0.35 the extra dose's +2.07 % cosine cancels k = 3's −1.82 % (net −0.17 %) at
   Δ`pds` +0.00119 (0.8 σ), but the union's `nmae` identity charges +0.00178 raw: **−0.00077 pessimistic / −0.00032
   point / −0.00006 optimistic**; k = 5 at 0.35 −0.0041; at g = 0.20 everything is −0.0059 … −0.0109 (§3).
5. **The `pds` is not direction-specific.** The split-half projector (basis fitted on the other half of the union's
   rows) buys **+0.0006 … +0.0016** — as much as or more than the in-sample projector — at a third of the cosine cost;
   the shuffled projector +0.0000 … +0.0009; the energy-matched random subspace +0.0012 at k = 1 and ≈ 0 / negative
   beyond. The "controls below half the arm" condition fails at k = 1, 2, 3 and holds at k = 5, 8 only because the random
   control turns negative (§5). This is 246 §4's control pattern inverted: for CD4 the split-half basis kept 85–90 % of
   the gain *because the subspace was a property of the source*; here the controls match the arm *because the gain is
   generic* — what any small perturbation of 74–77 rows plus a re-solved fraction does — and none of it is priceable.
6. **Leave-top-3-out retention** 8 % (k = 1), 28 % (k = 2), 26 % (k = 3), −45 % (k = 5), −38 % (k = 8); the top three
   targets carry 71–167 % of each fold's gain (`ZNF397`/`PI4KB`/`HMGN3`, `ZNF670`/`PI4KB`/`PPIF`, `SFXN4`/`PRNP`/`TIMM50`
   at k = 2, §6). Three targets per fold, the same shape 258/259 found for the union's `pds`.
7. **`kolf2` (apart) is negative on every cell**: Δ`pds` −0.0008 (k = 1) … −0.0058 (k = 8), cosine −3.0 … −16.6 % —
   the projector hurts most where the union's transfer is weakest (§3).
8. `nmae` by the validated identity: +0.0002 … +0.0003 raw at k ≤ 3 (−0.0001 overall pessimistic), +0.0015 at k = 5 / 8
   (−0.0006), +0.0018 … +0.0019 at g = 0.35 (−0.0007 … −0.0008); top-89 precision −0.0004 … −0.0020; q_w ≤ 0.0002.
   Delivered amplitude on the production panel (calibrated) −0.01 … −0.32 %; the leg fires only at k = 8 and at k = 5,
   g = 0.35 and moves the totals by < 0.0002 (§7).
9. The two code paths agree: the builder's own `apply_frozen_cd4_residual` with the pre-deflated union handed in as
   `extra` (what a `VCC_UNION_DEFLATE_K` flag would do) equals this experiment's compose at ≤ 1e−9 on all 12 cells of
   the production panel; the union's solved fraction there would rise ×1.023 (k = 1) … ×1.109 (k = 8), 77 carrying rows,
   the CD4 block byte-identical (§8).
10. **DOES NOT CLEAR.** The mechanism P249 §2 named — an un-whitened residual keeps a target-shared hub that the
    whitener strips from the donors — requires the hub to exist. CD4's residual (264–291 rows, one dominant direction)
    has one; the H1 pair's (74–77 rows, flat spectrum) does not, so its top-k directions are the on-target signal the
    incumbent is paid for. The analogue is closed for the union at every k and dose tested. No
    `PREREGISTRATION-DRAFT.md`.

---

## 1. Instrument and construction (`lib261.py`, `sweep261.py`; nothing mathematical is new)

162's PLOLO, three cancer folds (`k562_gwps` 269 targets / `orion_hct116` 290 / `orion_hek293t` 294; `kolf2` 273
apart), 246's construction imported through 257 (`lib257` → `lib246`: 218's residuals and `union_rows`, the shipped
stack, `deflate`, `compose`, `score`; `lib257.Identity` for 188's exact `nmae` identity on frame A, six246's
construction). **The base** is 257's `k5_s1`: CD4 deflated by `lib246.deflate(cd4, "svd", 5)` and re-sized by
`plolo.solve_fraction` at f = 0.34 (the RESOLVED variant, the package's object), the H1-pair union
(`gse295214_h1` + `h1_2025`, 74 / 75 / 76 / 71 carrying rows) at its shipped carrying-row fraction for g = 0.28.
**The arm** applies the *same* `lib246.deflate(res, "svd", k)` to the union object — per-gene centring over ITS
carrying rows (`target_condition_count > 0`, the rows one or both sources carry), support re-zeroed,
`G = U_c U_cᵀ`, `np.linalg.eigh`, `U_k` = top-k eigenvectors, `(I − U_k U_kᵀ) U_c`, support re-zeroed — and then
re-sizes the deflated union by 218's carrying-row rule (`compose218.solve_fraction_carrying`, the builder's
`solve_carrying_row_fraction` verbatim) at g against the same whitened base: deflate **before** sizing, 249's order.
One merged add (`merged_many`), one preserve. `mode="centre"` is the rank-1 mean alone; `mode="none"` at g = 0.28 is
the identity row R1. Sizing `fixed` keeps the shipped union fraction (the projector alone).

Controls, 246's verbatim: energy-matched random subspace (`lib257.rand_energy_matched`, 246's `rand`, seed
`RAND_SEED + 5`, bisection to the cell's removed energy), shuffled projector (`sweep246b.shuffled_projector`: top-k
fitted on a row-deranged copy, applied to the real rows, seed `SHUF_SEED + k`), split-half projector
(`sweep246c.split_half_deflate`: the right-singular basis fitted on a random half of the carrying rows applied to the
other half, both ways, seeds `RAND_SEED + 31·i + k`, i ∈ {0, 1}). Per row: `cos_agg`, macro cosine, `pds`, q_w,
per-target `pds`/cosine, `nmae_wald89`, top-89 precision, delivered norm (frame A). 39 rows per fold.

## 2. The price (`price261.py`; 257 §1.3's recipe, three coherent scenarios)

`price257.py` is imported for its constants and its `price()`; this file's `price()` is the same arithmetic with the
DE-three rate and the identity multipliers as free arguments and is asserted equal to 257's on 257's two factor sets
(0.0e+00, R2c). Legs: `pds` = 2.2138 × k_pds × Δ`pds` / 6; `mse`: Δcos_board = k_cos × rel `cos_agg` × 0.105470,
Δraw = −1.996 × Δcos_board (+ 243's amplitude leg — 197's identity at `cd4svd5`'s own `A_c`/`B_c` terms from
`amp249.json` — when the delivered amplitude moves > 0.2 %), overall = −1.051 × Δraw / 6; DE-three = rate ×
Δcos_board split `reach`/`fid`/`jac` 72.1 / 25.0 / 2.9 % by 218's card; `nmae` = −Δ`nmae_wald89` / 0.60 / 6 (three
cancer folds concatenated) × multiplier.

| scenario | k_pds {0.85, 1.22} | k_cos {1.36, 3.76} | DE-three {0.0425, 0.24} | identity |
|---|---|---|---|---|
| **PESSIMISTIC** | 0.85 on a `pds` gain, 1.22 on a loss | 3.76 on a cosine loss, 1.36 on a gain (one k_cos feeds `mse` and DE-three) | 0.24 on a loss, 0.0425 on a gain | ×1.5 on a loss, ×0 credit on a gain |
| **POINT** (= 257's point set, the one that retrodicts BR249) | 1.22 | 3.76 | 0.24 × 0.776 = 0.186 | ×1.0 |
| **OPTIMISTIC** | 1.22 on a gain, 0.85 on a loss | 1.36 on a loss, 3.76 on a gain | 0.0425 on a loss, 0.24 on a gain | ×1.0 |

The amplitude leg is read primarily from PLOLO's frame-A delivered-norm proxy (257's recipe, the one R2 validates);
the production-panel calibrated read (`prod261.py`, 258/259's reconstruction) is priced alongside (§7). CLEARED =
PESSIMISTIC ≥ +0.0019 AND leave-top-3-out ≥ 50 % AND Δ`pds` > 0 on 3/3 folds with the mean ≥ 2 × 0.001494 AND the
random and shuffled controls' Δ`pds` below half the arm's.

## 3. The k × g grid — every cell as an increment over `cd4svd5` (three cancer folds)

Columns: energy removed (fraction of the union residual's energy); the union's re-solved fraction relative to the
shipped 0.2237 / 0.2225 / 0.2252; Δ`pds` with positive-fold count and size in LOO-σ; relative aggregate and macro
cosine; delivered amplitude (PLOLO proxy / production panel calibrated); Δ`nmae_wald89` raw; the six-member legs at
POINT; the three totals; leave-top-3-out retention and the top-3 share; `kolf2` apart.

| cell | E rm | frac × | Δ`pds` | folds | σ | agg % | macro % | dlv PLOLO / prod % | Δ`nmae` | `pds` | `mse` | DE-three | `nmae` | **PESS** | **POINT** | **OPT** | ret3 | top3 | `kolf2` Δ`pds` / agg |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `u_c0` (centre only) | 0.017 | 1.0024 | +0.000530 | 3/3 | +0.4 | −0.171 | −0.178 | +0.007 / −0.012 | +0.0000 | +0.00024 | −0.00024 | −0.00013 | −0.00000 | **−0.00024** | −0.00013 | +0.00014 | −19 % | — | −0.00009 / +0.27 % |
| `u_k1` | 0.070 | 1.0201 | +0.001118 | 3/3 | +0.7 | −0.849 | −0.625 | +0.209 / −0.050 | +0.0002 | +0.00050 | −0.00111 | −0.00063 | −0.00005 | **−0.00165** | −0.00129 | +0.00004 | +8 % | 94 % | −0.00084 / −3.05 % |
| `u_k2` | 0.115 | 1.0256 | **+0.001540** | 3/3 | +1.0 | −1.246 | −0.845 | +0.131 / −0.115 | +0.0002 | +0.00069 | −0.00173 | −0.00092 | −0.00005 | **−0.00251** | −0.00201 | −0.00006 | +28 % | 71 % | −0.00184 / −6.17 % |
| `u_k3` | 0.157 | 1.0361 | +0.001149 | 3/3 | +0.8 | −1.824 | −1.185 | +0.015 / −0.164 | +0.0003 | +0.00052 | −0.00253 | −0.00135 | −0.00007 | **−0.00401** | −0.00343 | −0.00058 | +26 % | 80 % | −0.00217 / −9.20 % |
| `u_k5` | 0.229 | 1.0765 | +0.001090 | 3/3 | +0.7 | −2.999 | −1.833 | +0.266 / −0.237 | +0.0015 | +0.00049 | −0.00408 | −0.00221 | −0.00041 | **−0.00721** | −0.00621 | −0.00152 | −45 % | 167 % | −0.00412 / −12.71 % |
| `u_k8` | 0.309 | 1.1134 | +0.001299 | 3/3 | +0.9 | −3.364 | −2.215 | +0.168 / −0.319 | +0.0015 | +0.00058 | −0.00466 | −0.00248 | −0.00043 | **−0.00810** | −0.00699 | −0.00173 | −38 % | 148 % | −0.00577 / −16.61 % |
| *dose axis* | | | | | | | | | | | | | | | | | | | |
| `u_k0 g=0.20` (258 `pair g=0.2`) | 0 | 0.7143 | −0.000021 | 2/3 | −0.0 | −2.513 | −1.659 | +0.138 / −0.047 | −0.0001 | −0.00001 | −0.00348 | −0.00186 | +0.00002 | −0.00589 | −0.00533 | −0.00140 | — | — | −0.00317 / −8.92 % |
| `u_k3 g=0.20` | 0.157 | 0.7401 | +0.000622 | 2/3 | +0.4 | −3.843 | −2.588 | +0.098 / −0.148 | +0.0001 | +0.00028 | −0.00533 | −0.00284 | −0.00003 | −0.00884 | −0.00792 | −0.00191 | −76 % | 29 % | −0.00526 / −15.55 % |
| `u_k5 g=0.20` | 0.229 | 0.7689 | +0.000448 | 2/3 | +0.3 | −4.690 | −3.087 | +0.168 / −0.197 | +0.0002 | +0.00020 | −0.00650 | −0.00346 | −0.00006 | −0.01091 | −0.00982 | −0.00249 | −264 % | 84 % | −0.00700 / −18.07 % |
| `u_k0 g=0.35` (258 `pair g=0.35`) | 0 | 1.2500 | −0.000354 | 1/3 | −0.2 | +2.065 | +1.201 | +0.016 / −0.146 | +0.0011 | −0.00016 | +0.00286 | +0.00153 | −0.00029 | +0.00056 | +0.00394 | +0.00442 | — | — | +0.00213 / +7.61 % |
| `u_k3 g=0.35` | 0.157 | 1.2952 | +0.001191 | 3/3 | +0.8 | −0.171 | −0.164 | +0.111 / −0.135 | +0.0018 | +0.00054 | −0.00024 | −0.00013 | −0.00050 | **−0.00077** | −0.00032 | −0.00006 | +29 % | 76 % | +0.00004 / −3.79 % |
| `u_k5 g=0.35` | 0.229 | 1.3456 | +0.001396 | 3/3 | +0.9 | −1.626 | −0.920 | +0.191 / −0.307 | +0.0019 | +0.00063 | −0.00225 | −0.00120 | −0.00052 | −0.00414 | −0.00334 | −0.00080 | −17 % | 122 % | −0.00209 / −8.14 % |
| *projector at the SHIPPED fraction* | | | | | | | | | | | | | | | | | | | |
| `u_k3 fixed` | 0.157 | 1.0000 | +0.001098 | 3/3 | +0.7 | −2.066 | −1.345 | −0.144 / — | +0.0002 | +0.00049 | −0.00286 | −0.00153 | −0.00007 | −0.00459 | −0.00396 | −0.00074 | +21 % | 92 % | −0.00255 / −10.17 % |
| `u_k5 fixed` | 0.229 | 1.0000 | +0.001017 | 3/3 | +0.7 | −3.406 | −2.122 | +0.033 / — | +0.0004 | +0.00046 | −0.00472 | −0.00252 | −0.00010 | −0.00780 | −0.00688 | −0.00156 | −56 % | 326 % | −0.00470 / −14.01 % |

Per fold (k562_gwps / orion_hct116 / orion_hek293t): `u_k2` Δ`pds` +0.00071 / +0.00167 / +0.00224 at aggregate
cosine −0.845 / −1.845 / −1.049 %; `u_k3` +0.00069 / +0.00183 / +0.00092 at −1.764 / −2.870 / −0.837 %; `u_k5`
+0.00034 / +0.00084 / +0.00209 at −2.082 / −4.455 / −2.459 %; `u_k3 g=0.35` +0.00117 / +0.00171 / +0.00069 at +0.043 /
−1.194 / +0.640 %. The `pds` gain sits on `orion_hek293t` on most cells, the cosine loss on `orion_hct116`.

**Reading.** The projector's cost is monotone in k and the `pds` is not: k = 2 buys the most (+0.00154) and k = 5 buys
less (+0.00109) than k = 1 (+0.00112). The re-solved fraction (×1.02 … ×1.11) recovers a quarter of the cosine the
projector removes (`fixed` vs resolved at k = 3: −2.07 % → −1.82 %; k = 5: −3.41 % → −3.00 %) and moves `pds` by
< 0.0001, so the projector, not the sizing, is what costs. The dose axis is 258's family (my `u_k0` rows reproduce
258's `pair g=0.2` / `pair g=0.35` at 0.0e+00): the pair at 0.35 is a pure re-dose (+2.07 % cosine, `pds` −0.00035,
`nmae` +0.0011 raw) worth +0.00056 at the re-dose class's own k_cos 1.36 — 0.3× the gate, 258's reading, not a
candidate. Deflating on top of it (`u_k3 g=0.35`) trades that cosine back for +0.0015 `pds` and lands at −0.00077.

## 4. The spectrum, and the cosine per unit energy

| | carrying rows | rank-1 mean, % of raw energy | direction 1 / 2 / 3 / 4 / 5, % of centred | cum k = 5 | cosine per unit energy removed | `pds` per unit energy |
|---|---:|---:|---|---:|---:|---:|
| **union** k562_gwps | 74 | 1.68 | 5.23 / 4.55 / 4.45 / 3.84 / 3.61 | 21.7 % | −12.2 % (k = 1) … −10.9 % (k = 8) | 0.016 (k = 1) … 0.004 (k = 8) |
| union orion_hct116 | 75 | 1.79 | 5.40 / 4.61 / 4.27 / 3.73 / 3.52 | 21.5 % | | |
| union orion_hek293t | 76 | 1.77 | 5.35 / 4.59 / 4.24 / 3.70 / 3.49 | 21.4 % | | |
| union kolf2 | 71 | 1.80 | 5.68 / 4.97 / 4.46 / 4.08 / 3.85 | 23.0 % | | |
| union production panel | 77 | 1.68 | 5.24 / 4.59 / 4.20 / 3.83 / 3.55 | 21.4 % | | |
| **CD4** (257, 3 cancer folds) | 264–286 | 0.07 | **7.98** / 3.80 / 3.27 / 2.95 / 2.03 | 20.0 % | **−2.3 %** | **0.061** |

CD4's top direction is twice its second and 257 measured it carrying 75 % of `cd4svd5`'s `pds`; the union's leading
directions are within 15 % of each other and none of them carries anything (k = 1: +0.0011, 0.7 σ; the increments k → k+1
are +0.0004, −0.0004, −0.0001, +0.0002). Removing the same ~20 % of energy costs the union 3.0 % aggregate cosine where
it cost CD4 0.46 %, i.e. **the union's directions are the transfer signal**: 74–77 rows of direct H1 measurements of the
target entering at g = 0.28 are what the incumbent's cosine on those rows is made of, and there is no shared component
above them to strip. Exchange rate vs the STACK (246 §3.5, Δ`pds` per −1 % aggregate cosine): the whole stack with the
deflated union reads 0.0103 (k = 1) … 0.0036 (k = 8) against `cd4svd5`'s 0.0266; the increment over `cd4svd5` alone
0.0013 … 0.0004.

## 5. Controls at every k (vs `cd4svd5`; the arm's Δ`pds` and pessimistic total alongside)

| k | arm Δ`pds` / agg / PESS | random (dims, E) Δ`pds` / agg / PESS | shuffled projector (E) Δ`pds` / agg / PESS | split-half seed 0 (E) Δ`pds` / agg / PESS | split-half seed 1 Δ`pds` / agg / PESS | random & shuffled < ½ arm |
|---:|---|---|---|---|---|---|
| 1 | +0.001118 / −0.85 % / −0.00165 | (4, 0.075) **+0.001198** / −0.70 % / −0.00171 | (0.031) +0.000024 / −0.49 % / −0.00116 | (−0.017) **+0.001171** / −0.34 % / −0.00043 | +0.000586 / −0.35 % / −0.00065 | **no** |
| 2 | +0.001540 / −1.25 % / −0.00251 | (8, 0.121) +0.000012 / −0.97 % / −0.00227 | (0.046) **+0.000853** / −0.36 % / −0.00062 | (−0.011) +0.000631 / −0.40 % / −0.00075 | **+0.000861** / −0.37 % / −0.00060 | **no** |
| 3 | +0.001149 / −1.82 % / −0.00401 | (11, 0.164) +0.000045 / −1.25 % / −0.00291 | (0.057) **+0.000812** / −0.39 % / −0.00069 | (−0.001) **+0.001286** / −0.22 % / −0.00013 | +0.000812 / −0.32 % / −0.00053 | **no** |
| 5 | +0.001090 / −3.00 % / −0.00721 | (16, 0.233) −0.000926 / −2.12 % / −0.00547 | (0.084) +0.000483 / −0.77 % / −0.00172 | (0.007) **+0.001596** / −0.36 % / −0.00034 | **+0.001618** / −0.32 % / −0.00077 | yes (random negative) |
| 8 | +0.001299 / −3.36 % / −0.00810 | (23, 0.325) −0.001305 / −2.04 % / −0.00535 | (0.107) −0.000350 / −0.58 % / −0.00143 | (0.016) **+0.001461** / −0.29 % / −0.00074 | +0.001225 / −0.45 % / −0.00067 | yes (both negative) |

Three things the table says. (i) The split-half projector — a basis fitted on rows it is never applied to — buys
**as much or more `pds` than the arm at every k ≥ 3** (+0.0013 … +0.0016, 3/3 folds on five of the eight seed × k
cells) at a quarter to a sixth of the cosine cost, and removes almost no energy (−1.7 … +2.2 %, the fit half's mean
shifting the other half). (ii) The shuffled projector buys 0.05–0.8 of the arm's `pds` at a third of the cost. (iii) The
random subspace buys the arm's `pds` at k = 1 and nothing beyond. So the `pds` here is not a property of the fitted
directions: it is what any small re-centring / re-shaping of these 74–77 rows followed by a re-solved fraction does,
worth +0.0005 … +0.0016 by any route, and no route prices above **−0.00013** (the best control, `split_k3`). For CD4
(246 §3.4) the split-half basis kept 85–90 % of a +0.0123 gain because the subspace was the source's; here the
controls match a +0.0011 gain because there is no subspace to keep. `kolf2` (apart): every control negative except
`shufP_k1` / `shufP_k2` (+0.00001 / +0.00007).

## 6. Retention and the per-target share (the arm's increment over `cd4svd5`)

| cell | leave-top-3-out retention (3-fold; per fold) | top-3 share | top-3 targets per fold (Δ per-target `pds`) |
|---|---|---|---|
| `u_k1` | 8 % (−7 / 21 / 4 %) | 94 % | — |
| `u_k2` | 28 % (29 / 42 / 18 %) | 71 % | k562: `ZNF397` +0.057, `PI4KB` +0.041, `HMGN3` +0.039; hct116: `ZNF670` +0.106, `PI4KB` +0.090, `PPIF` +0.088; hek293t: `SFXN4` +0.249, `PRNP` +0.158, `TIMM50` +0.134 |
| `u_k3` | 26 % (10 / 39 / 13 %) | 80 % | — |
| `u_k5` | −45 % (−112 / −67 / −26 %) | 167 % | — |
| `u_k8` | −38 % (−61 / −56 / −29 %) | 148 % | — |
| `u_k3 g=0.35` | 29 % (37 / 35 / 2 %) | 76 % | — |
| `u_c0` | −19 % (−73 / −8437 / 35 %) | (undefined: the k562 and hct116 means are ≈ 0) | k562: `ZNF397` +0.039, `TWF2` +0.024, `HMGN3` +0.021; hct116: `ZNF670` +0.097, `POR` +0.053, `UROS` +0.045; hek293t: `SFXN4` +0.127, `TIMM50` +0.080, `PKN3` +0.061 |

No cell reaches the 50 % retention condition; at k ≥ 5 the gain outside the top three is negative. `SFXN4` on
`orion_hek293t` alone is +0.25 of per-target `pds` at k = 2 — 0.0008 of the fold's +0.0022.

## 7. `nmae`, precision, amplitude — the members that could not rescue it

- **`nmae`** (188's exact identity, three cancer folds concatenated; 4-fold in parentheses): `u_c0` +0.00001 (+0.00002),
  `u_k1` +0.00019 (+0.00024), `u_k2` +0.00019 (+0.00031), `u_k3` +0.00026 (+0.00046), `u_k5` +0.00148 (+0.00150),
  `u_k8` +0.00153 (+0.00163), `u_k3 g=0.35` +0.00178 (+0.00147), `u_k5 g=0.35` +0.00185 (+0.00167) raw — a loss
  everywhere (−0.00001 … −0.00077 overall at ×1.5), small next to `mse`. Top-89 precision −0.0004 … −0.0020; q_w
  −0.0001 … +0.0002.
- **Amplitude.** PLOLO's frame-A proxy reads +0.01 … +0.27 % (the leg fires at k = 1 and k = 5, g = 0.28); the
  production-panel calibrated read (`prod261.py`, factor 0.98555 = 258's) is **−0.012 (c0) / −0.050 / −0.115 / −0.164
  / −0.237 / −0.319 %** for c0 / k = 1 / 2 / 3 / 5 / 8 at g = 0.28 — the re-solved union fraction adds norm on 77 rows
  and the preserve takes it off everywhere — and fires at k = 8 and at k = 5, g = 0.35. Re-priced on the production
  read the pessimistic totals move by ≤ 0.00015 (`u_k1` −0.00165 → −0.00171, `u_k5` −0.00721 → −0.00736, `u_k8`
  −0.00810 → −0.00819); no sign, no ordering, no verdict changes. Every cell is shape-only (|Δ| < 0.5 %).
- **`mse` and DE-three** carry the verdict: at k_cos 3.76 the cosine loss is worth −0.0011 … −0.0047 on `mse` and
  −0.0006 … −0.0032 on `reach`/`fid`/`jac`; even at k_cos 1.36 and DE-three 0.0425 (OPTIMISTIC) the `mse` leg
  (−0.0004 … −0.0017) exceeds the `pds` leg at 1.22× (+0.0002 … +0.0007) on every k ≥ 2.

## 8. The builder flag it would need — `VCC_UNION_DEFLATE_K` (described, NOT made)

`scripts/make_submission_pm_cd4.py` (HEAD `3615d20f…`), one place:

- **Env** next to `VCC_CD4_DEFLATE_K` (line 740): `UNION_DEFLATE_K = os.environ.get("VCC_UNION_DEFLATE_K", "").strip()`,
  the same integer ≥ 1 validation, unset = `None` = byte-inert.
- **Action** in `apply_frozen_cd4_residual`, after the CD4 deflation branch and before the `residuals = [...]`
  accumulation (line ~2575): `if UNION_DEFLATE_K is not None and extra is not None: extra, union_deflation_audit =
  deflate_cd4_residual(extra, UNION_DEFLATE_K)`. `deflate_cd4_residual` is source-agnostic — it reads `mean_effect`,
  `gene_votes`, `target_condition_count`, and on the union object the count is the number of sources carrying the
  row, so `> 0` selects exactly the carrying rows `lib246._carried` selects. Union 1 then goes through
  `solve_carrying_row_fraction` at `H1_UNION_G` unchanged — deflate before sizing, 249's order. Union 2 (`extra2`)
  and CD4 untouched.
- **Meta**: `config.union_deflate_k` and `effects.h1_union.deflation.{k, n_carried, energy_removed_frac,
  top_k_energy_frac_of_centred, basis_u_k_sha256}` written only when set (the `h1_union_g` convention). A rule-10 diff
  against the incumbent would read: config +1; `effects.h1_union.residual_frame.{residual_median_row_norm_over_carrying_rows,
  g_per_unit_fraction, solved_fraction, applied_fraction}` moved, the two merged-preserve scalars, `effects.final_sha256`;
  added = the five deflation keys; `effects.cd4.*` byte-identical.

**What it would build is what was measured**: `prod261.py` pushed the pre-deflated union through the builder's own
`apply_frozen_cd4_residual` as `extra` and got this experiment's compose at ≤ 1e−9 on all 12 cells (KA-P cross-path
12/12), with the union's solved fraction on the real panel 0.19207 → 0.19648 / 0.19676 / 0.19882 / 0.20522 / 0.21292
(k = 1 / 2 / 3 / 5 / 8; ×1.023 … ×1.109), 77 carrying rows, energy removed 0.068 / 0.114 / 0.155 / 0.227 / 0.310,
preserve 0.911588 → 0.911622 at k = 8, the CD4 fraction 0.34104860341941284 unchanged. Not made: the flag has
nothing to build.

## 9. Files, memory, provenance

| file | what | peak RSS |
|---|---|---|
| `lib261.py` | the construction (257 → 246 imported; `cd4svd5_parts`, `union_parts`, 246's two projector controls verbatim, `spectrum`) | — |
| `sweep261.py --fold <held>` → `work/sweep261-<held>.{json,log}`, `pertarget261-<held>.npz`, `identity261-<held>.npz` | 39 rows per fold, 44 known-answer checks + 5 row KAs per fold, exit 1 on any failure | 5.00 / 4.91 / 4.96 / 4.92 GiB; 34–47 s |
| `prod261.py` → `work/prod261.{json,log}` | the production frame off 249's staged artifacts; KA-P1–P4 + 12 cross-path cells (37/37); amplitude terms and calibration | 1.22 GiB; 17 s |
| `price261.py` → `work/price261.{json,log,stdout}` | R1–R3, the grid, the controls, retention, the three scenarios, the verdict | 0.04 GiB |
| `work/plolo/` | byte copies of 257's fold caches, `SHA256SUMS.txt` | |

Builder imported by `prod261.py` only, read-only, `3615d20f8b535c2f…` (= 258/259's bytes). Sibling 260 at 5.9 GB
alongside; memory 73–74 % free throughout. Seeds recorded in every control's audit.

## 10. Verdict

**DOES NOT CLEAR.** 0 of 10 real cells meet the pessimistic gate (best −0.00024, the centring-only cell; the projector
cells −0.0017 … −0.0081), 0 of 10 meet leave-top-3-out ≥ 50 % (best 28 %), 0 of 10 reach 2 LOO-σ on `pds` (best
1.0 σ), and the random / shuffled / split-half controls buy the same +0.0005 … +0.0016 by any route. The H1-union
residual has no target-shared hub — leading direction 5.2 % of a flat centred spectrum against CD4's 8.0 % dominant
one — so every direction removed is on-target transfer signal, costing 5× the cosine per unit energy that CD4's
deflation cost for one-tenth the `pds`. **The CD4 deflation's analogue on the union is closed for k ∈ {1, 2, 3, 5, 8},
g ∈ {0.20, 0.28, 0.35}, resolved and fixed sizing, on the `cd4svd5` base.** No `PREREGISTRATION-DRAFT.md`, no flag,
no build, no submission.
