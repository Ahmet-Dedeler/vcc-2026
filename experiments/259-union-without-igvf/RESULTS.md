# 259 — The union WITHOUT `igvf_cardiomyocyte` on top of `cd4svd5`: **the attribution is right, and removing IGVF does not buy a cleared arm.** Dropping IGVF removes the `nmae` cost (Δ`nmae` +0.0012 vs +0.0097 … +0.0154 at the same dose, e +0.26 vs −0.01) and the second-residual construction B (pair kept at 0.28 + `dld1,jiang,myllia` at g₂) is the best rendering of these sources this project has measured — **+0.00184 PESSIMISTIC at g₂ = 0.35 (0.97× the gate), +0.00211 at the out-of-grid g₂ = 0.45 (1.11×)** — but its `pds` is 2/3 folds at 0.8–1.5 LOO-σ and leave-top-3-out is negative at every dose ≥ 0.28. Construction A (single union replacing the pair) is worse at every dose (best +0.00099). **0 of 21 real cells meet the three conditions. DOES NOT CLEAR — verdict line at the bottom (§10).**

**Date 2026-09-05 (02:45–03:10 UTC).** Deterministic. **Measurement only: no build, no `vcc`, no `cell-eval2`, no model
name, no description, nothing under `submissions/`, nothing written outside this directory** (`find -newer work/.marker-start`
over `experiments/ data/ scripts/ docs/ submissions/ src/` returns only `scripts/__pycache__/make_submission_pm_cd4.cpython-312.pyc`,
the bytecode cache Python writes when the builder is imported; no source file changed). Builder: git HEAD's
`scripts/make_submission_pm_cd4.py` (`3615d20f…`, the same bytes 258 priced on), imported with the incumbent's env and driven
through its own `apply_frozen_cd4_residual` — union 1 via `extra`, union 2 via 248's `extra2` / `audits_out["union2"]` —
with the module globals as the knobs. Fold caches are byte copies of 258's (= 249's; `work/plolo/SHA256SUMS.txt` diffed
equal). Peak RSS per process in §8: 5.11–5.42 GiB per fold, 3.41 GiB for the production frame, 0.25 GiB for the price;
one fold at a time; the sibling 260 ran alongside (0.4 GB).

---

## 0. Verdict

```
R1  258's rows reproduce from THIS runner                        PASS  4 folds x 338/338 (222 scalars at 0.0e+00 and
    (cd4svd5, stack6 g=0.35, stack3 g=0.35, pair+igvf g=0.28,           77 per-target/identity arrays BIT-IDENTICAL per fold),
     h1union, ARM6, a15)                                                98/98 cross-fold aggregates and 258's own conservative /
                                                                        point totals and retention at 0.0e+00; prod frame 21/21
R2  (a) BR237's ARM6 by member at 258 s0's factor set            PASS  pds +0.00246 (board +0.00248), mse +0.00088 (+0.00088),
        k_cos 0.96, k_pds 2.48 (ABC), identity x1.5,                    DE-three +0.00016 (+0.00016), nmae -0.00717 (-0.00716),
        package amplitude -1.24 %, DE-three 0.0425                      TOTAL -0.00367 (board -0.00364)
    (b) BR249's cd4svd5 vs h1union at k_cos 3.76 / k_pds 1.22    PASS  +0.00470 vs +0.0043 +- 0.0008 (board +0.00426);
                                                                        pds +0.00552 (+0.00545), mse -0.00079 (-0.00080)
R3  248's second residual through MY union2 path                 PASS  per fold 39 rows x 4 folds at 0.0e+00 (CUM1@0.28/0.35,
                                                                        CUM4@0.35: pds, cos_agg, macro, fractions, preserve);
                                                                        cross-fold 14/14 (rel macro +1.405 %, d pds +0.00256,
                                                                        leave-top-3-out 25.2 %, share 81.9 % = 225/248's digits);
                                                                        PRODUCTION FRAME: 248's staged effects-cum1.npy
                                                                        (6225a582) and effects-cum4.npy (2f9eba51) re-composed
                                                                        to max|d| 2.2e-16, fractions/preserve at 0.0e+00

THE SCREEN (21 real cells + 9 controls on top of cd4svd5, three cancer folds; three coherent scenarios):
    PESSIMISTIC >= +0.0019            1 / 21   B3 g2=0.45 (+0.00211; beyond the mandate's grid)
    POINT       >= +0.0019            8 / 21
    OPTIMISTIC  >= +0.0019           17 / 21
    ALL THREE conditions              0 / 21   (pessimistic gate AND leave-top-3-out >= 50 % AND pds > 0 on 3/3 folds at >= 2 sigma)

    construction A (single union, no igvf):  best A5 g=0.35   PESS +0.00099  POINT +0.00157  OPT +0.00831   pds 2/3, 0.5 s, ret3 -292 %
    construction B (pair@0.28 + union2):     B3 g2=0.35       PESS +0.00184  POINT +0.00331  OPT +0.01291   pds 2/3, 1.3 s, ret3  -45 %
                                             B3 g2=0.45       PESS +0.00211  POINT +0.00305  OPT +0.01403   pds 2/3, 0.8 s, ret3 -161 %
    (i)  igvf back (B4) at the same g2:      PESS -0.0009 / -0.0031 / -0.0055 vs B3 at 0.15 / 0.28 / 0.35 -- the attribution holds
    (ii) dld1 alone (248's rider) g2=0.35:   PESS +0.00094  POINT +0.00257  OPT +0.00707   pds 3/3 (+0.00217, 1.5 s), ret3 15 %

DOES NOT CLEAR.  Removing igvf_cardiomyocyte removes the union's nmae cost, as 258 attributed; what is left is a
cosine arm (+4.3 ... +5.7 % aggregate, worth +0.0015 ... +0.0019 on mse even at k_cos 0.94) whose pds is four targets
on two folds and negative on the third, and the three clearing conditions are pds conditions.
```

**Three findings beyond the gate.** (1) **The IGVF attribution is confirmed directly** (§4.1): adding IGVF back to the
second residual costs Δ`nmae` +0.0031 / +0.0084 / +0.0142 raw at g₂ = 0.15 / 0.28 / 0.35 for +0.0004 `pds` and +0.9 pp
cosine — −0.0009 / −0.0031 / −0.0055 pessimistic. (2) **258's "the deflation took half of DLD-1's `pds`" was two effects
in one number** (§6.3): at g = 0.35 DLD-1 reads +0.00148 as a single union on `h1union` (229 ARM3), +0.00256 as a second
residual on `h1union` (248 CUM1), +0.00119 as a single union on `cd4svd5` (258), **+0.00217 as a second residual on
`cd4svd5` (here)**. The construction is worth ×1.7–1.8; the deflated base costs 15–20 %. (3) **Inside B the `pds` is DLD-1's
three `k562_gwps` targets plus ONE `jiang` target on `orion_hct116`** (`IFNAR2`, +0.88 of per-target `pds` = the whole of
that fold's gain; §5); `jiang` carries 9 rows and `myllia` 13, and `myllia` buys +1.9 % aggregate cosine for negative
`pds`. Every B cell loses `pds` on `orion_hek293t`, which is the fold IGVF's `pds` lives on.

---

## 1. The instrument (`price259.py`, three phases; 258's code paths with the union2 route added)

**Construction.** On each of 162's PLOLO folds the production builder's `apply_frozen_cd4_residual` is called with the
fold's whitened base: the CD4 residual is deflated at k = 5 before `solve_cd4_fraction` sizes it at f = 0.34; union 1
(`prepare_h1_union(sources=…)`, registry order) is sized by the carrying-row rule at `H1_UNION_G`; union 2 (`extra2`) by
the same rule at `UNION2_G` against the same base and accumulated after union 1 (`merged = cd4·f + u1·f1 + u2·f2`,
support OR'd, counts max'd — 248's order, verified in §2 R3 against 248's own fold rows and staged arrays); ONE merged
add, ONE global preserve. Construction A sets union 1 = the five sources and no union 2; construction B keeps union 1 =
the H1 pair at 0.28 and puts the three (or four, or one) sources in union 2. Each arm is scored on the two functionals
258 carried: **(a)** 218's `score_effects` — `cos_agg`, macro cosine, `pds`, mean over frames A/B/C, the frame 243's
class rule and 246/249's Δ`pds` live in; **(b)** 216/229's frame-A functional — 188's exact `nmae` identity on three
pools (4 folds pooled, `kolf2` included, as `six229.py` does), top-89 precision, delivered norm, `pds` on frame A.
Controls for B: 218's target-shuffle of union 2 (a derangement of its carried rows, seeds 20726 / 771103 — `resid218.shuffled`
verbatim) and 225's F3 row-matched control (CD4's own information on exactly the rows union 2 carries, `compose218.restrict_rows`
verbatim), each through the same path at the same g₂. **Construction A does not admit a clean energy-matched control**: a
target-shuffle of the five-source union would also scramble the H1 pair the incumbent carries, so the B controls (same
three sources, same rows) stand in for it and this is stated rather than faked.

**The production frame (`--prod`)**, 258 §1's reconstruction: `wh = final/s − merged` off 249's staged `effects-off.npy`,
every arm re-composed on the real 300-target panel, delivered amplitude read as 211/243's (`A_c`, `B_c`) terms and
calibrated on the three packages (factor 0.986, 258's value at 0.0e+00).

**The price (`--price`)** — three COHERENT scenarios per cell, never one coefficient set (READ-FIRST v2.7 §3; 257 §0.9;
258 §0). Legs: `pds` = 2.2138 × k_pds × Δ`pds`(ABC, 3 cancer folds) / 6; `mse` = 197's identity at `cd4svd5`'s own terms
(`amp249.json`; reproduces its board raw 0.8218062 at 0.0e+00) with Δcos_board = k_cos × rel `cos_agg` × 0.105469 and
the calibrated amplitude leg **always charged**; DE-three = rate × Δcos_board, split reach / fid / jac 72.1 / 25.0 / 2.9 %
by 218's card as 257 did; `nmae` = −Δ`wald89` / 0.60 / 6 × x.

| scenario | k_pds {0.70, 2.48} | k_cos {0.94, 3.76} | DE-three {0.0425, 0.24} | identity x |
|---|---|---|---|---|
| **PESSIMISTIC** | 0.70 on a `pds` gain, 2.48 on a loss | 0.94 on a cosine gain, 3.76 on a loss (feeds `mse` AND DE-three the same way) | 0.0425 on a gain, 0.24 on a loss | ×1.5 on a loss, **×0 credit on a gain** |
| **POINT** (ARM6-realised for union legs) | 2.48 | 0.96 | 0.0425 | ×1.0 |
| **OPTIMISTIC** | 2.48 on a gain, 0.70 on a loss | 3.76 on a gain, 0.94 on a loss | 0.24 on a gain, 0.0425 on a loss | ×1.0 (a gain credited in full at ×1.0; ×1.5 credit is not a validated point) |

Every cell also carries 258's conservative set (`C258`: k_pds 0.7 on the FRAME-A `pds`, k_cos 0.96, DE-three 0.0425, ×1.5)
for continuity with 258's tables; the frame-A `pds` is reported, never priced in the scenarios. **CLEARED** = PESSIMISTIC
≥ +0.0019 AND leave-top-3-out retention ≥ 50 % AND Δ`pds` > 0 on 3/3 folds with |Δ`pds`| ≥ 2 × 0.001494 = 0.002988.

---

## 2. The known-answer and retrodiction rows

### R1 — 258's rows reproduce at 0.0e+00
Per fold, for the seven shared arms (`a15`, `h1union`, ARM6 = `stack6 g=0.35 on h1union`, `cd4svd5`, `stack6 g=0.35`,
`stack3 g=0.35`, `pair+igvf g=0.28`): every numeric key of 258's `fold-<held>.json` row (222 scalars: `cos_agg`, macro,
`pds`, `pds_A`, precision, delivered norm, fractions, preserve, energy removed, …) equal at **0.0e+00**, and all 77
per-target / identity arrays (`pds_pt`, `cos_pt`, U/M/D on three pools) **bit-identical** — 4 × 338/338 (the 338 includes
the R3 rows below). Cross-fold: 98/98 — the twelve aggregates per arm (`cos_agg_3c`, `pds_3c`, `pds_A_3c`, `nmae_wald89`
4-fold and 3-fold, `lfc89`, `h89`, precision, delivered ratio, `U_sum`, `M_sum`), the amplitude calibration factor, the
measured board cosine, and for the three screen cells 258's own `conservative.total_A` (−0.00513 / +0.00119 / −0.00583),
`point.total_A` (−0.00092 / +0.00343 / −0.00389) and leave-top-3-out retention (−98 % / −66 % / −377 %) through my `legs()`
— all at 0.0e+00 (1e−9 asked). Production frame 21/21 (§2 R3 and 258's fourteen).

### R2 — the two board points, from this runner's rows (the factor set 258 §0 found, stated)
Realised factors from these rows: ARM6 k_cos 0.962, k_pds 2.48 (ABC) / 3.79 (A), DE-three 0.0425 per unit board Δcos,
identity −0.00478 vs board −0.00716 (ratio 1.497); `cd4svd5` k_cos 3.759, k_pds 1.22 (ABC).

| retrodiction | `pds` | `mse` | DE-three | `nmae` | TOTAL | board |
|---|---:|---:|---:|---:|---:|---:|
| **(a) BR237 ARM6**: k_cos 0.96, k_pds 2.48 (ABC), DE-three 0.0425, identity ×1.5, package amplitude −1.24 % | +0.00246 | +0.00088 | +0.00016 | −0.00717 | **−0.00367** | −0.00364 (`pds` +0.00248, `mse` +0.00088, DE3 +0.00016, `nmae` −0.00716) |
| (a) with DE-three 0.24 | +0.00246 | +0.00088 | +0.00091 | −0.00717 | −0.00293 | the 0.24 recipe over-predicts DE-three 5.6× (258 R2) |
| **(b) BR249 `cd4svd5` vs `h1union`**: k_cos 3.76, k_pds 1.22, DE-three at the realised 0.185, ×1.0, package amplitude −0.41 % | +0.00552 | −0.00079 | −0.00034 | +0.00031 | **+0.00470** | +0.00426 (band +0.0043 ± 0.0008; `pds` +0.00545, `mse` −0.00080, DE3 −0.00034, `nmae` −0.00005) |

Both PASS their bands. Honesty, as 258 §2 recorded: k_cos and k_pds in (a) and (b) are the values those board points
fixed, so those legs pass by construction; the non-circular content is the identity at ×1.5, the amplitude leg and the
DE-three rate.

### R3 — 248's second residual through my union2 path
Per fold (39 rows × 4 folds, all **0.0e+00**): my `h1union + union2 dld1@0.28 / @0.35` and `+ union2 CUM4@0.35` rows equal
248's `repro248-<fold>.json` arms on `pds`, `cos_agg`, macro, `a15_sgn_w`, median row norm, union-2 solved fraction and
carrying rows, CD4 fraction, union-1 fraction and the preserve scalar; my `h1union` row equals 248's base. Cross-fold
14/14: rel macro **+1.405 %** / rel `cos_agg` +2.266 % / Δ`pds` **+0.002557** (3/3) for `dld1@0.35` (the mandate's
digits), `kolf2` −0.00254, and CUM4@0.35 +4.037 % / +5.779 % / +0.003609; leave-top-3-out of `dld1@0.35` **25.19 %**, top-3
share 81.9 % = 225's `conc225.json` (248 KA-E). **Production level:** 248's staged `effects-cum1.npy` (`6225a582…`) and
`effects-cum4.npy` (`2f9eba51…`) re-composed off 249's recovered base through `extra2` to max|Δ| 2.2e−16, union-2
fractions 0.19847054819075655 / 0.20712104380637059, carrying 67 / 139, preserves 0.8944861646401419 / 0.8696966167072415
at 0.0e+00. The mixing order is 248's, not invented. The mandate's spelling of the four-source set (`dld1,jiang,myllia,igvf`)
and 248's (`dld1,jiang,igvf,myllia`) build the same object: `union_residual` is a per-gene mean, max|Δ| 5.6e−17, votes equal
on every fold; 248's order is the one priced.

---

## 3. The two constructions on top of `cd4svd5` (three cancer folds; `kolf2` apart)

Columns: c1 / c2 = union-1 / union-2 carrying rows on the `orion_hek293t` fold (/294); agg / mac = rel `cos_agg` / macro
cosine vs `cd4svd5`; Δ`pds` ABC with the per-fold values (k562 / hct116 / hek293t), the count of positive folds and its
size in LOO-σ (0.001494); Δ`pds` A = the frame-A functional; Δ`nmae` = the exact identity (4 folds, `wald89`); e = ΔM/ΔU
(break-even 0.5); amp = calibrated delivered amplitude; **PESS / POINT / OPT** = the three coherent six-member totals;
C258 = 258's conservative set; ret3 = leave-top-3-out `pds` retention; `kolf2` apart.

### 3.1 A — the single union without IGVF replacing the pair (`VCC_UNION_SOURCES=gse295214_h1,h1_2025,dld1_lowmoi,jiang_5line_noK562,myllia_a549`, `VCC_H1_UNION_G=g`)

| g | c1 | agg % | mac % | Δ`pds` ABC (k562 / hct / hek) | + | σ | Δ`pds` A | Δ`nmae` | e | amp % | **PESS** | POINT | OPT | C258 | ret3 | `kolf2` agg / Δ`pds` |
|---:|---:|---:|---:|---|:-:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0.20 | 133 | −1.256 | −0.072 | +0.00175 (+0.0011 / +0.0036 / +0.0005) | 3 | 1.2 | +0.00165 | +0.00141 | 0.09 | +0.05 | **−0.00305** | +0.00073 | +0.00074 | −0.00064 | −41 % | −11.08 % / −0.0063 |
| 0.28 | 133 | +1.429 | +1.843 | +0.00147 (+0.0022 / +0.0027 / −0.0005) | 2 | 1.0 | +0.00158 | +0.00151 | 0.14 | −0.37 | **+0.00018** | +0.00137 | +0.00414 | +0.00022 | −81 % | −3.30 % / −0.0044 |
| **0.35** | 133 | +3.489 | +3.130 | +0.00073 (+0.0027 / +0.0013 / −0.0018) | 2 | 0.5 | +0.00100 | +0.00055 | 0.25 | **−0.93** | **+0.00099** | +0.00157 | +0.00831 | +0.00109 | −292 % | +3.17 % / −0.0031 |
| 0.45 | 133 | +5.975 | +4.454 | −0.00096 (+0.0028 / −0.0013 / −0.0044) | 1 | 0.6 | −0.00048 | +0.00150 | 0.21 | **−1.63** | +0.00024 | +0.00049 | +0.01262 | +0.00104 | — | +11.83 % / −0.0016 |

**A fails at every dose, for 258's two reasons plus one.** Below 0.28 the single shared dose under-doses the pair the
incumbent carries (cosine −1.3 %, `kolf2` −11 %); at 0.35 and 0.45 the delivered amplitude moves −0.93 % / −1.63 % —
beyond the 0.5 % shape-only bound, so those cells are amplitude-changing arms facing the +0.0089 gate, not +0.0019 — and
the `pds` peaks at g = 0.20 (+0.00175, 3/3, 1.2 σ) where the cosine is negative. What removing IGVF did buy is visible in
the identity: Δ`nmae` +0.0006 … +0.0015 with e +0.09 … +0.25, against 258's six-source +0.0144 … +0.0159 with e −0.02 … −0.05.
The five-source union's `pds` is also smaller than the three-source second residual's at every dose: on the shared rows the
union takes a per-gene MEAN, so the pair and the added sources dilute each other (§6.3).

### 3.2 B — the pair kept at 0.28 + the SECOND residual (`VCC_H1_UNION_G=0.28`, `VCC_UNION2_SOURCES=dld1_lowmoi,jiang_5line_noK562,myllia_a549`, `VCC_UNION2_G=g₂`)

| g₂ | c2 | agg % | mac % | Δ`pds` ABC (k562 / hct / hek) | + | σ | Δ`pds` A | Δ`nmae` | e | amp % | **PESS** | POINT | OPT | C258 | ret3 | `kolf2` agg / Δ`pds` |
|---:|---:|---:|---:|---|:-:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0.15 | 81 | +2.596 | +2.397 | +0.00225 (+0.0025 / +0.0039 / +0.0004) | **3** | 1.5 | +0.00203 | +0.00081 | 0.21 | −0.09 | **+0.00122** | +0.00283 | +0.00787 | +0.00119 | 16 % | +3.56 % / −0.0004 |
| 0.28 | 81 | +4.271 | +3.610 | +0.00227 (+0.0031 / +0.0039 / −0.0002) | 2 | 1.5 | +0.00237 | +0.00124 | 0.24 | +0.03 | **+0.00174** | +0.00344 | +0.01173 | +0.00180 | −13 % | +6.11 % / −0.0019 |
| **0.35** | 81 | +4.951 | +3.974 | +0.00193 (+0.0030 / +0.0034 / −0.0007) | 2 | 1.3 | +0.00226 | +0.00120 | 0.26 | −0.24 | **+0.00184** | +0.00331 | +0.01291 | +0.00196 | −45 % | +7.28 % / −0.0027 |
| 0.45 (beyond the grid) | 81 | +5.665 | +4.218 | +0.00119 (+0.0027 / +0.0023 / −0.0015) | 2 | 0.8 | +0.00164 | +0.00075 | 0.28 | −0.23 | **+0.00211** | +0.00305 | +0.01403 | +0.00228 | −161 % | +8.70 % / −0.0038 |

**B is the better construction and it is a cosine arm.** Its pessimistic total is the `mse` leg at k_cos 0.94 (+0.0015 …
+0.0019 of the +0.0017 … +0.0021) plus a 0.70-weighted `pds` leg of +0.0003 … +0.0006, minus the identity at ×1.5 (−0.0003 …
−0.0005). The `pds` is 3/3 folds only at g₂ = 0.15 and never reaches 2 σ (best +0.00227 = 1.5 σ against the 0.00299 bar);
from g₂ = 0.28 up `orion_hek293t` is negative and leave-top-3-out is negative. The amplitude stays within the shape-only
bound at every g₂ (−0.24 … +0.03 %). `kolf2` gains +3.6 … +8.7 % cosine and loses `pds`, 229 §7's inflation pattern.

---

## 4. The controls that make the reading honest

### 4.1 (i) IGVF added back — `VCC_UNION2_SOURCES=dld1_lowmoi,jiang_5line_noK562,igvf_cardiomyocyte,myllia_a549` (248's CUM4 set)

| g₂ | c2 | agg % | Δ`pds` ABC (k562 / hct / hek) | + | σ | Δ`nmae` | e | amp % | **PESS** | POINT | OPT | ret3 | IGVF's marginal over B3 at the same g₂: Δ`pds` / agg / Δ`nmae` / PESS / POINT |
|---:|---:|---:|---|:-:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0.15 | 135 | +3.378 | +0.00271 (+0.0018 / +0.0037 / +0.0026) | 3 | 1.8 | +0.00389 | 0.10 | −0.14 | **+0.00034** | +0.00269 | +0.00924 | 18 % | +0.00046 / +0.78 pp / **+0.00308** / **−0.00088** / −0.00014 |
| 0.28 | 135 | +5.215 | +0.00269 (+0.0014 / +0.0033 / +0.0034) | 3 | 1.8 | +0.00966 | 0.03 | −0.06 | **−0.00133** | +0.00183 | +0.01195 | −16 % | +0.00042 / +0.94 pp / **+0.00842** / **−0.00306** / −0.00161 |
| 0.35 | 135 | +5.824 | +0.00237 (+0.0010 / +0.0026 / +0.0035) | 3 | 1.6 | +0.01536 | −0.01 | −0.48 | **−0.00370** | +0.00004 | +0.01132 | −50 % | +0.00044 / +0.87 pp / **+0.01416** / **−0.00553** / −0.00327 |

**The attribution is right.** IGVF adds +0.0004 `pds` (it is the source that is positive on `orion_hek293t`, so the
four-source cells are 3/3 where the three-source cells are 2/3) and +0.9 pp of cosine, and charges +0.003 … +0.014 raw on
the identity with e falling to 0.03 … −0.01 — the same "magnitudes, not signs" signature as ARM6 (e −0.033), which the board
realised at 150 %. Every B4 cell prices below its B3 counterpart at PESSIMISTIC and POINT alike.

### 4.2 (ii) DLD-1 alone as the second residual — 248's rider re-priced on `cd4svd5`

| g₂ | c2 | agg % | Δ`pds` ABC (k562 / hct / hek) | + | σ | Δ`nmae` | e | amp % | **PESS** | POINT | OPT | C258 | ret3 (top-3 share) |
|---:|---:|---:|---|:-:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0.15 | 65 | +1.364 | +0.00146 (+0.0021 / +0.0014 / +0.0008) | 3 | 1.0 | +0.00086 | 0.24 | +0.22 | +0.00062 | +0.00171 | +0.00437 | +0.00062 | 37 % |
| 0.28 | 65 | +2.094 | +0.00205 (+0.0032 / +0.0020 / +0.0010) | 3 | 1.4 | +0.00139 | 0.24 | +0.37 | +0.00089 | +0.00245 | +0.00652 | +0.00096 | 23 % |
| **0.35** | 65 | +2.313 | +0.00217 (+0.0036 / +0.0020 / +0.0010) | 3 | 1.5 | +0.00137 | 0.27 | +0.13 | **+0.00094** | **+0.00257** | +0.00707 | +0.00106 | **15 %** (107 %) |

248 read this rider at +0.00216 under 243's class rule on the `h1union` base (1.14× gate, retention 25.2 %). On `cd4svd5`
it keeps its `pds` sign on 3/3 folds (+0.00217, −15 % vs 248's +0.00256), keeps its cosine (+2.31 % vs +2.27 %), and
prices +0.00094 pessimistic / +0.00257 point; retention falls to 15 % (top-3 `LENG1 / TRIM35 / MAP2K2` on k562,
`TARBP2 / SH2D3A / CYTH2` on hct116, `ZNF219 / LPCAT3 / CYTH2` on hek293t — 248's targets). Its production-frame union-2
fraction is exactly 248's 0.19847054819075655 (the same base, the same rule).

### 4.3 Attribution inside B: `jiang` alone and `myllia` alone

| union 2 | g₂ | c2 | agg % | mac % | Δ`pds` ABC (k562 / hct / hek) | + | Δ`nmae` | **PESS** | POINT | ret3 |
|---|---:|---:|---:|---:|---|:-:|---:|---:|---:|---:|
| `jiang_5line_noK562` | 0.28 | 9 | +0.529 | +0.714 | +0.00058 (−0.0001 / **+0.0027** / −0.0009) | 1 | +0.00008 | +0.00030 | +0.00070 | −162 % |
| | 0.35 | 9 | +0.638 | +0.836 | +0.00033 | 1 | +0.00008 | +0.00027 | +0.00050 | −380 % |
| `myllia_a549` | 0.28 | 13 | **+1.872** | +1.340 | −0.00032 (+0.0006 / −0.0003 / −0.0012) | 1 | −0.00021 | +0.00049 | +0.00056 | — |
| | 0.35 | 13 | +2.291 | +1.534 | −0.00043 | 1 | −0.00024 | +0.00060 | +0.00068 | — |

`jiang`'s whole `pds` is one target on one fold (`IFNAR2` on `orion_hct116`, +0.878 of per-target `pds` ≈ the fold's
+0.0027); `myllia` carries 13 rows, buys +1.9 % aggregate cosine (high-energy rows) and loses `pds` on 2/3 folds. The
three legs are close to additive on cosine (2.09 + 0.53 + 1.87 = 4.49 % vs B3's 4.27 % at 0.28).

### 4.4 Energy-matched controls for B (same rows, same rule, same g₂)

| g₂ | control | agg % | mac % | Δ`pds` ABC | + | σ | **PESS** | POINT | OPT | vs the real B3 |
|---:|---|---:|---:|---:|:-:|---:|---:|---:|---:|---|
| 0.15 | SHUF 20726 | −0.108 | +0.043 | −0.00120 | 1 | 0.8 | −0.00155 | −0.00133 | −0.00054 | real +0.00122 / +0.00283 |
| | SHUF 771103 | −0.438 | −0.510 | −0.00172 | 0 | 1.2 | −0.00266 | −0.00180 | −0.00067 | |
| | row-matched CD4 | **+2.892** | +1.466 | −0.00097 | 0 | 0.6 | +0.00024 | +0.00032 | +0.00657 | |
| 0.28 | SHUF 20726 | −0.733 | −0.595 | −0.00343 | 0 | 2.3 | −0.00504 | −0.00357 | −0.00132 | real +0.00174 / +0.00344 |
| | SHUF 771103 | −1.344 | −1.507 | −0.00446 | 0 | 3.0 | −0.00741 | −0.00480 | −0.00186 | |
| | row-matched CD4 | **+4.851** | +2.012 | −0.00254 | 0 | 1.7 | −0.00076 | −0.00063 | +0.01046 | |
| 0.35 | SHUF 20726 | −1.262 | −1.137 | −0.00501 | 0 | 3.4 | −0.00760 | −0.00512 | −0.00182 | real +0.00184 / +0.00331 |
| | SHUF 771103 | −2.021 | −2.199 | −0.00612 | 0 | 4.1 | −0.01055 | −0.00648 | −0.00244 | |
| | row-matched CD4 | **+5.701** | +2.126 | −0.00365 | 0 | 2.4 | −0.00244 | −0.00198 | +0.01147 | |

The target-shuffle is negative on cosine and `pds` at every dose and both seeds (the real B3 beats it by +0.0028 …
+0.0124 pessimistic): the gain is direction-specific. **But CD4's own information on the same rows buys MORE aggregate
cosine than the three sources** (+4.85 % vs +4.27 % at 0.28, +5.70 % vs +4.95 % at 0.35) while losing `pds` — 248 §4.1's
finding again. On the functional the class rule prices, the sources' cosine is not better than the free alternative; what is
source-specific is the `pds`, and that is the concentrated, 2/3-fold channel.

---

## 5. The best cells in full

### 5.1 B3 g₂ = 0.35 — the best cell inside the mandate's grid

| quantity | value |
|---|---|
| PLOLO vs `cd4svd5` | rel `cos_agg` **+4.951 %** (3.55 / 5.74 / 5.56), rel macro +3.974 %; Δq_w +0.00050; Δ precision-89 +0.00204 |
| Δ`pds` frames ABC | **+0.001932** (k562 **+0.00301**, hct116 **+0.00343**, hek293t **−0.00066**); 2/3 positive; **1.3 LOO-σ**; frame A +0.00226 |
| `nmae` identity | Δ`wald89` **+0.001200** (3 cancer folds +0.00170); ΔU, ΔM → **e 0.263** (signs, unlike ARM6's −0.03) |
| delivered amplitude | effect-level −0.244 %, calibrated −0.241 % (shape-only; charged −0.00008) |
| production-frame composition | union 1 77/300, fraction 0.19207161964444322 (untouched); union 2 **84/300 carrying** (25 shared with union 1, 2 outside CD4 support), solved fraction **0.19843869479105122**; CD4 fraction 0.34104860341941284 (untouched); preserve **0.8897130320431748** (`cd4svd5` 0.9115882187626339); deflation energy 0.198660 (untouched) |
| `kolf2` (apart) | rel `cos_agg` +7.28 %, Δ`pds` −0.00272 |
| **leave-top-k-out `pds` retention** | k = 1 **23 %**, **k = 3 −45 %**, k = 5 −93 %; top 3 per fold: k562 `LENG1 +0.201, TRIM35 +0.188, MAP2K2 +0.165` · hct116 **`IFNAR2 +0.878`**, `TARBP2 +0.217, SH2D3A +0.170` · hek293t `ZNF219 +0.205, ZNF22 +0.204, LPCAT3 +0.152` (the top-3 "share" is ill-defined here: the hek293t fold's total is negative) |

| scenario | k_cos | k_pds | DE-three | identity × | Δcos_board | raw `mse` | `pds` | `mse` | `nmae` | `reach` | `fid` | `jac` | (amplitude inside `mse`) | **TOTAL** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **PESSIMISTIC** | 0.94 | 0.70 | 0.0425 | 1.5 | +0.004909 | 0.81250 | +0.00050 | +0.00163 | −0.00050 | +0.00015 | +0.00005 | +0.00001 | −0.00008 | **+0.00184** |
| **POINT** | 0.96 | 2.48 | 0.0425 | 1.0 | +0.005013 | 0.81230 | +0.00176 | +0.00167 | −0.00033 | +0.00015 | +0.00005 | +0.00001 | −0.00008 | **+0.00331** |
| **OPTIMISTIC** | 3.76 | 2.48 | 0.24 | 1.0 | +0.019635 | 0.78321 | +0.00176 | +0.00677 | −0.00033 | +0.00340 | +0.00118 | +0.00014 | −0.00008 | **+0.01291** |
| 258-conservative (k_pds 0.7 on frame A) | 0.96 | 0.70 | 0.0425 | 1.5 | | | +0.00058 | +0.00167 | −0.00050 | | | | | +0.00196 |

Gate: PESSIMISTIC ≥ +0.0019 → **no** (0.97×); leave-top-3-out ≥ 50 % → **no**; `pds` 3/3 at ≥ 2 σ → **no**. **NOT CLEARED.**

### 5.2 B3 g₂ = 0.45 — the only cell over the pessimistic gate (beyond the mandate's grid)
rel `cos_agg` +5.665 % (3.82 / 6.68 / 6.50); Δ`pds` ABC **+0.001191** (+0.00273 / +0.00229 / **−0.00145**; 2/3; 0.8 σ), frame A
+0.00164; Δ`nmae` +0.000752, e 0.278; amplitude −0.234 %; union 2 fraction 0.2551354647313516, preserve 0.8790175027896644.
PESSIMISTIC **+0.00211** = `pds` +0.00031 + `mse` +0.00188 + `nmae` −0.00031 + DE-three +0.00024 (amplitude −0.00008 inside
`mse`); POINT +0.00305; OPTIMISTIC +0.01403. Retention k = 1 −36 %, **k = 3 −161 %**. It is over the gate by 0.0002 on
the cosine leg alone while its `pds` falls; both `pds` conditions fail. **NOT CLEARED.**

### 5.3 A5 g = 0.35 — the best A cell
rel `cos_agg` +3.489 % (2.05 / 3.86 / 4.55), macro +3.13 %; Δ`pds` ABC +0.000728 (+0.00272 / +0.00129 / −0.00183; 2/3; 0.5 σ);
Δ`nmae` +0.000552, e 0.249; **amplitude −0.93 % (calibrated), beyond the shape-only bound**; union 136/300 carrying,
fraction 0.23357153481174606, preserve 0.8826226999155439. PESSIMISTIC **+0.00099** (`pds` +0.00019, `mse` +0.00088 of which
amplitude −0.00032, `nmae` −0.00023, DE-three +0.00015); POINT +0.00157; OPTIMISTIC +0.00831. Retention k = 3 −292 % (the
same targets as B3 plus `IFNAR2`). Fails all three conditions.

---

## 6. What the screen settles

1. **The `nmae` attribution to `igvf_cardiomyocyte` is right, measured directly** (§4.1): the same three sources with and
   without IGVF at the same g₂ differ by +0.003 … +0.014 raw on the identity and by −0.0009 … −0.0055 pessimistic. Without
   IGVF every union cell has e ≥ +0.09 — the direction channel, not magnitudes.
2. **Removing IGVF does not buy a cleared arm.** The remaining value is cosine: B3's pessimistic total is 85–90 % `mse` leg
   at k_cos 0.94; its `pds` is +0.0019 … +0.0023 (1.3–1.5 σ) on two folds, negative on `orion_hek293t` (the fold IGVF's
   `pds` lives on), and four targets (`LENG1`, `TRIM35`, `MAP2K2` on k562; `IFNAR2` on hct116) carry more than the whole
   of it. The three clearing conditions are `pds` conditions, and it fails two of them at every dose. A cosine-first reading
   would be a different mandate — and the row-matched CD4 control shows the aggregate-cosine gain is not source-specific
   (§4.4), so it would not be a strong one.
3. **Construction beats base as the reason DLD-1's `pds` looked halved on `cd4svd5`** (258 §6 read +0.00256 → +0.00119).
   At g = 0.35, three cancer folds:

   | construction \ base | `h1union` | `cd4svd5` |
   |---|---:|---:|
   | single union with the pair (229 ARM3 / 258 `stack3`) | +0.00148 (agg +2.28 %) | +0.00119 (agg +2.34 %) |
   | second residual, pair kept at 0.28 (248 CUM1 / 259 B1) | **+0.00256** (agg +2.27 %) | **+0.00217** (agg +2.31 %) |

   The second-residual construction is worth ×1.7–1.8 on `pds` at equal cosine (on the 19 rows both carry, the single union
   takes the per-gene mean and dilutes both the pair and DLD-1; union 2 adds DLD-1 at its full dose); the deflated base
   costs 15–20 %. 248's flag is the right rendering for any source rider; 229/258's single-union re-dose is the wrong one.
4. **Inside the three-source residual the sources do different things**: DLD-1 carries the `pds` (3/3 folds alone, 15 %
   retention), `jiang` is one target on one fold, `myllia` is cosine for negative `pds`. A "union without IGVF" is DLD-1 plus
   cosine riders, and DLD-1 alone was already priced at 1.14× the gate with 25 % retention on the old base (248) and reads
   +0.00094 pessimistic / 15 % retention on this one.
5. **`kolf2` inflates every B cell** (+3.6 … +8.7 % cosine) while its `pds` falls — the pluripotent fold is not the cancer
   folds and is reported apart, as since 218.
6. The frame-A `pds` functional and the ABC functional agree in sign on every real cell and differ by ≤ 0.0005; nothing
   here turns on which one is read.

---

## 7. The builder configuration the best cells would need — recorded, NOT recommended, NOT built

All existing flags; no new flag. On top of `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_CD4_DEFLATE_K=5` (A/B/C; the driver's
`DEF_ENV` carries `VCC_CD4_DEFLATE_BASIS=V5_abc` instead on any other bundle) and 222's build flags (`--whitener
data/donors/whitener_union_pruned.npz --whitener-id union_pruned --min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread
0.25 --overdispersion 4.0 --emitter rc_callset`):

```
B3 g2=0.35:  VCC_H1_UNION_G=0.28  VCC_UNION2_SOURCES=dld1_lowmoi,jiang_5line_noK562,myllia_a549  VCC_UNION2_G=0.35
A5 g=0.35:   VCC_UNION_SOURCES=gse295214_h1,h1_2025,dld1_lowmoi,jiang_5line_noK562,myllia_a549  VCC_H1_UNION_G=0.35
```

Rule-10 diff B3 would show against `submissions/249-cd4svd5-f034.h5ad.meta.json` (248 §7's pre-declared diff): config ADDED
exactly `union2_g` 0.35, `union2_sources` (three names, four shas each from the registry), `union2_sizing_rule`; audit
CHANGED exactly `effects.final_sha256` and the merged-preserve scalars under `effects.cd4` and `effects.h1_union`
(0.9115882187626339 → 0.8897130320431748); audit ADDED only `effects.union2.*` (84 carrying, 25 shared with union 1, solved
fraction 0.19843869479105122); `effects.cd4.fraction`, the deflation block, `effects.h1_union.residual_frame` and the
operator sha byte-identical. Neither cell is worth a slot today; that is the only sentence about a build this experiment
supports. The 2026-09-06 slots stay unspent on this family.

---

## 8. Cost and memory (one fold per process, gated on no ≥ 15 GiB job and ≥ 40 % free; the sibling 260 ran alongside)

| script / phase | runtime | peak RSS |
|---|---:|---:|
| `price259.py --fold k562_gwps` (269 targets, 37 arms, 11 unions) | 40 s | 5.11 GiB |
| `--fold orion_hct116` (290) | 56 s | 5.41 GiB |
| `--fold orion_hek293t` (294) | 55 s | 5.42 GiB |
| `--fold kolf2` (273) | 47 s | 5.34 GiB |
| `--prod` (production-frame reconstruction, 36 compositions + 4 files, 3 renders each) | 34 s | 3.41 GiB |
| `--price` | 2 s | 0.25 GiB |

Wall clock from mandate (02:43 UTC) to write-up: about 30 minutes, of which 4 minutes of compute; the mandate's 3-hour box was not needed.

---

## 9. Honest limits

1. **R2 is a consistency check on its `pds` and `mse` legs** (the factors are the board points' own); the verdict does not rest
   on it. It rests on the sign structure of the cells: a 2/3-fold, 1.3-σ `pds` whose leave-top-3-out is negative, under a gate
   whose three conditions are all about `pds`.
2. **The k_pds band {0.70, 2.48} is two one-point class estimates** (h1union's realised, ARM6's realised on a broad six-source
   gain). Applied to a concentrated three-source gain the truth may sit anywhere in it; that is exactly why the scenarios span
   it. Even at 2.48 (POINT) no B cell passes the retention or 2-σ condition.
3. **The `nmae` identity is a proxy** (fixed top-89 pools, four folds with `kolf2`; 237 ADDENDUM-2 §2); its 3-cancer-fold value
   is same-signed and within 50 % on every cell (B3 g₂ = 0.35: +0.00170 vs +0.00120). At ×1.5 it costs B3 −0.0005: it does not
   decide a cell.
4. **The amplitude read is effect-level**, calibrated on three package pairs (0.954–1.007); B moves ≤ 0.24 %, A −0.93 …
   −1.63 % (which puts A's g ≥ 0.35 cells under the amplitude-changing gate regardless).
5. **`fid` / `jac` / `reach` yield halves are unreadable on any harness** (rule 9); they are carried at the two realised
   DE-three rates (0.0425 / 0.24 per unit board Δcos), ±0.0003 at PESSIMISTIC / POINT.
6. **The ordering channel's noise is sd 0.0045 overall** (230); every total here is inside it. The `mse`/`pds` channel is
   deterministic, which is the only reason a +0.002 reading is discussable at all.
7. **Construction A has no clean energy-matched control** (§1); B's controls stand in.
8. **g₂ = 0.45 was added beyond the mandate's grid** because the first reading put B's pessimistic optimum at the grid's top;
   it is the only pessimistic-gate crossing and it is a cosine crossing. No finer grid was run: between 0.35 and 0.45 the
   `pds` falls monotonically (1.3 σ → 0.8 σ) and retention worsens, so no interior dose can satisfy the two `pds` conditions.
9. **No package, no board member measured**; the post-build falsifiers (rule-10 diff, 900-block digest, `n_pred`, `vcc prep
   --dry-run`) were not reached because nothing clears.

## 10. Files

| file | what |
|---|---|
| `price259.py` | the instrument: `--fold` (37 arms through the builder's union-1 + union-2 path, 338 KA rows per fold), `--prod` (21 KA rows incl. 248's staged arrays, amplitude for every arm), `--price` (R1–R3, the tables, scenarios, retention, verdict) |
| `run259.sh` | the folds one process at a time, gated (258's gate) |
| `work/fold-<held>.{json,npz,stdout,stderr}` | per-fold rows, per-target `pds`/cosine, the identity's U/M/D arrays, KA rows |
| `work/prod259.{json,stdout}` | production-frame terms per arm, calibration, KA rows |
| `work/price259.{json,log}` | every number in this file |
| `work/plolo/` | 258's fold caches, byte copies (`SHA256SUMS.txt` equal to 258's) |
| `work/run259.log` | timings and peak RSS per fold |

```bash
cd ~/Code/virtual-cell-challenge-2026
bash experiments/259-union-without-igvf/run259.sh k562_gwps orion_hct116 orion_hek293t kolf2
uv run python experiments/259-union-without-igvf/price259.py --prod
uv run python experiments/259-union-without-igvf/price259.py --price
```

---

**DOES NOT CLEAR.** The union without `igvf_cardiomyocyte` on top of `cd4svd5` clears +0.0019 at the PESSIMISTIC scenario in
one cell only (construction B, g₂ = 0.45, +0.00211, beyond the mandate's grid, on its cosine leg), and no cell of either
construction meets the retention (≥ 50 %) or the `pds` (3/3 folds at ≥ 2 LOO-σ) condition; the best grid cell, B3 g₂ = 0.35,
is +0.00184 / +0.00331 / +0.01291 pessimistic / point / optimistic with `pds` 2/3 folds at 1.3 σ and leave-top-3-out −45 %.
The attribution to IGVF is right; removing it does not buy a cleared arm. No `PREREGISTRATION-DRAFT.md`. No model name, no
description, nothing built, nothing submitted.
