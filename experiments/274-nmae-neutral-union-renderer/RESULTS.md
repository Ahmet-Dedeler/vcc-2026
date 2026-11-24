# 274 — The union's DIRECTION cannot be rendered without its MAGNITUDE cost, because on the scored DE set the two are one object. **0 of 39 buildable renderings and 0 of 11 oracle renderings clear the pessimistic gate over `cd4svd5`; the best buildable cell (259's B3 at g₂ = 0.35 with the incumbent's top-400 magnitudes held, +0.00164 pessimistic) is BELOW its un-rendered parent (+0.00184). DOES NOT CLEAR — verdict line at the bottom (§11).**

**Date 2026-09-06 (03:30–04:35 UTC).** Deterministic (bootstrap seed 20260906, 20,000 draws; permutation seeds 20260906 / 771103; 218's
derangement seed 20726). **Measurement only: no build, no `vcc`, no `cell-eval2`, no model name, no description, nothing under
`submissions/`, nothing written outside this directory** (`find … -newer work/.marker-start` over `experiments data scripts src docs
submissions` returns only this directory and the sibling 275 agent's files). Builder: git HEAD's `scripts/make_submission_pm_cd4.py`
(`3615d20f…`, the bytes 258 and 259 priced on), imported with the incumbent's env and driven through its own `apply_frozen_cd4_residual`
(258's route; union 2 through 248's `extra2` path as 259). Fold caches are byte copies of 258's (= 249's; `work/plolo/SHA256SUMS.txt`
checked). `uv run python`, `PYTHONDONTWRITEBYTECODE=1`, one fold per process; peak RSS 5.66–6.00 GiB per fold, 5.41 GiB for the
production frame, 0.85 GiB for the price (§9). Two passes: pass 1 (87 arms) and pass 2 (91 arms, the exact-pool oracle added) agree on
all 3,144 shared numbers at 0.0e+00 (`work/pass1/`).

---

## 0. Verdict

```
KA  per fold 216/216 x 4 (k562_gwps, orion_hct116, orion_hek293t, kolf2); production frame 9/9; cross-fold 139/139
KA-1  ARM6's identity leg on the h1union base        -0.0047829354 vs the mandate's -0.0047829 (|d| 3.5e-08);
      nmae_wald89 1.0029111952 vs 1.0029112 (4.8e-09)                                                   PASS
KA-2  (a) BR237 by member at 258 s0's factor set:    pds +0.00246 (board +0.00248)  mse +0.00088 (+0.00088)
      DE-three +0.00016 (+0.00016)  nmae -0.00717 (-0.00716)  TOTAL -0.00367 (board -0.00364)          PASS
      (b) BR249 at k_cos 3.76 / k_pds 1.22:          pds +0.00552 (+0.00545)  mse -0.00079 (-0.00080)
      TOTAL +0.00470 (board +0.00426; band +0.0043 +- 0.0008)                                            PASS
      (both FITTED calibrations on their pds/mse legs; the identity's x1.5, the amplitude leg and the DE-three rate
      are the independent content -- R258 s2)
THE RENDERER'S OWN ROWS  identity check (X_arm := X_ref) reproduces cd4svd5 at 0.0e+00 on every scalar and per-target
      array (4 renderers x 4 folds); re-preserve(recover(final)) == the builder's final at <= 2.2e-16; the rank
      transplant's sorted magnitudes == the incumbent's at 0.0e+00; the changed rows == the union's carried rows
      (165 / 159 / 172 / 172) with the rounding floor at 1.1e-16

THE SCREEN (three cancer folds, kolf2 apart; three coherent scenarios; 39 buildable + 11 oracle + 28 control cells):
    buildable cells   PESSIMISTIC >= +0.0019   0 / 39      POINT >= gate  5 / 39 (all five are 259's B3 rendered,
                                                            each BELOW its un-rendered parent)   OPTIMISTIC  8 / 39
    ALL conditions (pess. gate AND pds > 0 on 3/3 with the pooled interval above 0 AND leave-top-3-out >= 50 %
                    AND amplitude within 0.5 % AND every control below half)                         0 / 39
    oracle cells (the truth's wald89 pool held)   PESSIMISTIC >= gate   0 / 11;  best  B3 g2=0.35 | wald89 exact  +0.00139
    best buildable  B3 g2=0.35 | own400 ratio   +0.00164 pess. / +0.00302 point / +0.00399 opt.;  its un-rendered
                    parent (R259's cell)         +0.00184 / +0.00331 / +0.00442
    the six-source union at g = 0.35, every renderer:  pessimistic -0.00068 (oracle exact) ... -0.00505; point <= +0.00032
    controls: 24 random-direction renderings, max d pds +0.00012, max pessimistic -0.00092 (every one below half)

DOES NOT CLEAR.  The union's direction cannot be rendered without its magnitude cost: 54 % of the six-source union's
transfer-cosine gain sits on the truth's 89 Wilcoxon-pool genes (6 % of the truth's energy, 4 % of the prediction's),
where the union raises the incumbent's magnitude by 4.4 % and moves sign accuracy by 0.2 pp.  Holding those 89 genes
per row removes the identity's cost and the direction gain together; every renderer that keeps the gain keeps the cost.
```

**Three findings beyond the gate.** (1) **The union's `pds`/cosine gain is a magnitude re-weighting of the DE pool, not sign
information** (§4): on the truth's `wald89` pool the six-source union's sign accuracy is 0.5772 vs the incumbent's 0.5749, its L1
1.044× the incumbent's, and its contribution to the inner product with the truth at EQUAL magnitude (the exact-pool oracle) is
−0.00029 — worse than the incumbent's. (2) **The rank transplant is exactly `nmae`-neutral and exactly `pds`-destroying**: with the
union's order and signs and the incumbent's magnitude spectrum per row (R2b, `all rank`), Δ`nmae` +0.00087 [−0.0003, +0.0029] and
Δ`pds` −0.00249 on 0/3 folds, `cos_agg` −0.22 %: the union's directions are worth nothing without the union's magnitudes. (3) **Even
the `nmae`-cheap union pays for rendering**: every rendering of 259's B3 cell (pair kept + `dld1, jiang, myllia` at g₂ = 0.35) sits
0.0002–0.0004 below the un-rendered cell at every scenario, because its `nmae` cost was already small (+0.0012 raw, e 0.26) and the
renderer only removes cosine.

---

## 1. The renderer (`render274.py`; frozen before the first fold was measured)

**Where it sits.** The builder composes `final = s · (wh + merged)`: `E062.add_post_whiten_residual` (`experiments/062-cd4-post-whiten/
run.py:571-594`, called from `scripts/make_submission_pm_cd4.py:2647`) forms `candidate = whitened + fraction · merged_effect` (line 580)
and then ONE global preserve scalar `s = median‖wh[:, wcols]‖ / median‖candidate[:, wcols]‖` (lines 581-588). A renderer stage lives
between those two lines: it sees the pre-preserve candidate of the ARM (`X_arm = wh + merged_arm`) and of the REFERENCE (`X_ref = wh +
merged_ref`, the incumbent = CD4 deflated k = 5 at f = 0.34 + the H1 pair at g = 0.28), and hands a rendered candidate to the same
preserve. Here `X_arm` and `X_ref` are recovered from the builder's own finals as `final / s` (`recover()`); rows the arm did not change
are snapped to the reference (the two differ there only by the rounding of two preserve scalars, 1.1e-16 relative; the changed rows are
exactly the union's carried rows — KA rows in §2), and the same preserve rule is re-applied. The identity renderer reproduces the
builder's `final` at ≤ 2.2e-16.

**The family.** On every changed row, a magnitude-control set S and a matching mode:

| set | S | buildable? |
|---|---|---|
| `wald89` | the identity's truth pool: the row's top-89 genes by \|lfc\|·√cpm of the TRUTH (188's pool) | **no — ORACLE**, the ceiling |
| `ownK` | the incumbent's own top-K entries by \|X_ref\| (K = 89 / 400 / 2000; panel columns excluded) | yes |
| `supp` | the entries the union changed (\|X_arm − X_ref\| above the rounding floor; 2,536–18,001 per row, median 9,200–11,800) | yes |
| `all` | the whole row (18,233 non-panel columns) | yes |

| mode | X_out on S | reading |
|---|---|---|
| `ratio` | X_arm[S] · (Σ\|X_ref[S]\| / Σ\|X_arm[S]\|) | the L1 total on S held at the incumbent's, the direction within S kept (R1) |
| `norm` | X_arm[S] · (‖X_ref[S]‖ / ‖X_arm[S]‖) | the L2 norm on S held (`all` + `norm` = R2a: the union's unit vector at the incumbent's norm) |
| `rank` | \|X_out[S]\| = sorted \|X_ref[S]\| placed at the ranks of \|X_arm[S]\|, signs of X_arm | the incumbent's magnitude SPECTRUM on the union's order and signs (`all` + `rank` = R2b) |
| `exact` (`wald89` only) | fixed point: the pool's L1 held at the incumbent's POST-preserve value | the oracle ceiling proper (the ratio mode holds it pre-preserve and the union's off-pool magnitudes then shrink it by the preserve, 0.9645×) |

Off-S entries stay as composed; then the preserve, then the panel columns zeroed, exactly as the builder does. **R1** = `wald89` /
`ownK` sets on the six-source union at g = 0.35; **R2** = `all`; **R3** = the six renderers on `stack3` (pair + `dld1_lowmoi`) at g =
0.35 / 0.28 and on the six-source union at g = 0.28 / 0.20; **R4** = 259's construction B (the pair kept at 0.28 as union 1 + `union2`
of `dld1, jiang, myllia` at g₂ = 0.35 / 0.45, and with `igvf` back at 0.35) rendered as R2a / R2b / `own400 ratio` / `supp rank` /
`wald89`. The renderer is context-independent by construction (one effect matrix for A/B/C), so a √cpm-weighted prediction set was not
used: it would differ per context.

**Controls.** (i) The identity check: the renderer applied with X_arm := X_ref (the pair at 0.28 composed again) must reproduce the
incumbent — it does, at 0.0e+00 on every scalar and per-target array for `wald89 ratio`, `all rank`, `all norm`, `own400 ratio`. (ii)
R0: the unions at their own magnitudes (258 / 259's rows, bit-identical). (iii) Energy-matched random directions rendered the same
way: the union's per-row contribution d = X_arm − X_ref replaced by a gene permutation of itself over its own support (`PERM`, two
seeds) or by another carried row's contribution (`SHUF`, 218's derangement), on every family — 24 rendered control cells + 4 un-rendered.

**Measurement.** 258's two functionals from six renders per arm: (a) 218's `score_effects` — `cos_agg`, macro cosine, per-target `pds`
(the scorer's rank statistic, 262/273's form, a23 render, mean over frames A/B/C); (b) 216/229's frame-A functional — 188's exact
`nmae` identity on three pools (`wald89`, `lfc89`, `h89`; 4 folds pooled as `six229.py` does, and the 3 cancer folds apart), top-89
precision, delivered norm. Plus: the effect-level delivered amplitude (243's A_c terms on the production 300-target panel through 249's
staged base, calibrated 0.986 against three packages — 258 §3; the fold frame for the oracle set, which cannot be composed without
the truth); paired target bootstraps of Δ`pds` (stratified by fold, 20,000 draws) and of the per-target identity; leave-top-k-out
retention (222's definition); Δ`pds` by 262's base-cosine group and the count of targets crossing the 0.02 own-cosine knee; `kolf2`
apart. Price: 259's three coherent scenarios with the new-source bands on the union legs (k_cos 0.94 / 0.96 / 1.00, k_pds 0.7 / 2.48,
DE-three 0.0425 / 0.24 per unit board Δcos split by 218's rate card, the identity ×1.5 on a loss and ×0 on a gain pessimistic, `mse`
from 197's identity at `cd4svd5`'s own terms with the amplitude leg always charged), and 258's conservative / point sets alongside.

---

## 2. Known-answer rows

| row | against | rows | result |
|---|---|---:|---|
| per fold: 258's rows for a15 / cd4 f=0.34 / h1union / ARM6 / cd4svd5 / stack6 g=0.35, 0.28, 0.20 / stack3 g=0.35, 0.28 (8 scalars + 5 per-target arrays each); 259's rows for B3 g₂=0.35, 0.45 and B4+igvf g₂=0.35 (7 scalars + 3 arrays); 246 / 249 / 229 rows as 258 checked them | 258 / 259 / 246 / 249 / 229 | 171 × 4 folds | **all 0.0e+00** |
| per fold, the renderer: re-preserved incumbent == builder's final; `recover()` changed rows == the union's carried rows (165 / 159 / 172 / 172), rounding floor ≤ 1.1e-16; `all rank` sorted \|X_out\| == sorted \|X_ref\| on every changed row (0.0e+00) and row-norm ratio − 1 ≤ 2.2e-16; `wald89 ratio` pool L1 ratio − 1 ≤ 4.4e-16; `wald89 exact` POST-preserve pool L1 ratio − 1 ≤ 1.2e-14 with the off-pool entries untouched (0.0e+00); untouched rows == the reference (0.0e+00) and == the arm's own rows (≤ 1.1e-16); the four IDENT cells == `cd4svd5` (rank statistic exact, scalars ≤ 1.1e-16, no row changed); re-preserve(recover(stack6)) == the builder's final ≤ 2.2e-16 | this experiment | 45 × 4 folds | **216/216 per fold** |
| production frame: recovered base median 4.229295604454299; `prepare_h1_union(PAIR)` == `union-off.npz` bit for bit; re-composed `h1union` == `abeb0749…` (5.6e−17); re-composed `cd4svd5` == 249's `effects-svd5.npy` `a28d5be2…` (2.2e−16), preserve 0.9115882187626339; re-composed ARM6 == 229's `c5721825…` (2.2e−16); the renderer's re-preserved incumbent == `effects-svd5.npy` (3.3e−16) | 249 / 229 | 9 | **9/9** |
| cross-fold: 258's aggregates for 10 arms × 11 keys and 259's for 3 arms × 5 keys (0.0e+00); `six229.json` ARM6 / h1union `nmae_wald89`; **KA-1**; the amplitude calibration and `cd4svd5`'s measured board cosine 0.105469 == 258; 258's conservative total, retention and 259's POINT total for stack6 / stack3 at g = 0.35 (≤ 1e−9) | 258 / 259 / 229 | 139 | **139/139** |
| **KA-2** (a) BR237 by member, (b) BR249 | the board | 2 | **both inside the bands** (§0) |

---

## 3. The screen — every cell vs `cd4svd5`, three cancer folds

Columns: agg / mac = relative `cos_agg` / macro cosine; Δ`pds` = 246's functional (per-fold sign count; the paired target bootstrap of
the three-fold mean); Δ`nmae` = the exact identity on `wald89` (4 folds / 3 cancer); e = ΔM/ΔU (break-even 0.5); amp = calibrated
delivered amplitude; PESS / POINT / OPT = the coherent six-member totals; C258 = 258's conservative set; ret3 = leave-top-3-out `pds`
retention (shown where \|Δ`pds`\| ≥ 0.0005). The full tables, the per-fold bootstrap intervals (T6) and **every cell's six members
priced three ways (T5, 91 × 3 rows)** are in `work/tables274.md`.

### 3.1 R0 — the unions at their own magnitudes (258 / 259's rows, reproduced at 0.0e+00)

| cell | agg % | mac % | Δ`pds` | + | bootstrap | Δ`nmae` 4f | 3c | e | amp % | **PESS** | POINT | OPT | ret3 |
|---|---:|---:|---:|:-:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| stack6 g=0.35 (ARM6's union on `cd4svd5`) | +3.69 | +2.96 | +0.00191 | 3/3 | [−0.0022, +0.0063] | +0.01565 | +0.01606 | −0.03 | −1.05 | **−0.00496** | −0.00151 | −0.00068 | −98 % |
| stack6 g=0.28 | +1.79 | +1.84 | +0.00243 | 3/3 | [−0.0011, +0.0063] | +0.01436 | +0.01664 | −0.02 | −0.71 | **−0.00490** | −0.00129 | −0.00089 | −38 % |
| stack6 g=0.20 | −0.82 | +0.09 | +0.00230 | 2/3 | [−0.0007, +0.0057] | +0.00761 | +0.00783 | −0.02 | −0.20 | **−0.00315** | −0.00040 | −0.00039 | −26 % |
| stack3 g=0.35 (pair + dld1) | +2.34 | +1.49 | +0.00119 | 3/3 | [−0.0012, +0.0036] | +0.00057 | +0.00073 | +0.33 | +0.32 | **+0.00109** | +0.00197 | +0.00250 | −66 % |
| stack3 g=0.28 | +0.45 | +0.41 | +0.00158 | 3/3 | [−0.0004, +0.0036] | +0.00144 | +0.00168 | +0.20 | +0.14 | **+0.00003** | +0.00128 | +0.00138 | −13 % |
| B3 g₂=0.35 (259's best) | +4.95 | +3.97 | +0.00193 | 2/3 | [−0.0013, +0.0054] | +0.00120 | +0.00170 | +0.26 | −0.24 | **+0.00184** | +0.00331 | +0.00442 | −45 % |
| B3 g₂=0.45 | +5.67 | +4.22 | +0.00119 | 2/3 | [−0.0025, +0.0051] | +0.00075 | +0.00112 | +0.28 | −0.23 | **+0.00211** | +0.00305 | +0.00432 | −161 % |
| B4+igvf g₂=0.35 | +5.82 | +4.20 | +0.00237 | 3/3 | [−0.0015, +0.0066] | +0.01536 | +0.01563 | −0.01 | −0.48 | **−0.00370** | +0.00004 | +0.00135 | −50 % |

The six-source union's `pds` gain has its paired-bootstrap interval across zero on every dose (P(≤ 0) 0.19 at g = 0.35); its identity
leg is 3.1 bootstrap-SDs from zero (+0.01565 [+0.0068, +0.0262], P(≤ 0) 5e-5). That is the object the renderer has to work on: a
certain cost and an uncertain gain.

### 3.2 R1 / R2 — the six-source union at g = 0.35, every set × mode; the trade-off frontier

| cell | agg % | Δ`pds` | + | bootstrap | Δ`nmae` 4f | e | amp % | **PESS** | POINT | OPT | cos kept | `pds` kept | `nmae` kept |
|---|---:|---:|:-:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **`wald89` ratio** (oracle) | +0.84 | −0.00004 | 2/3 | [−0.0041, +0.0043] | +0.00251 | +0.74 | −0.63 | **−0.00097** | −0.00061 | −0.00040 | 23 % | −2 % | 16 % |
| **`wald89` exact** (oracle, post-preserve) | +1.64 | +0.00084 | 2/3 | [−0.0032, +0.0051] | +0.00314 | — | −0.63 | **−0.00068** | +0.00032 | +0.00069 | 44 % | 44 % | 20 % |
| `wald89` rank (oracle) | +0.48 | −0.00092 | 1/3 | [−0.0048, +0.0032] | −0.00107 | +0.39 | −0.64 | −0.00087 | −0.00057 | +0.00029 | 13 % | −48 % | −7 % |
| `own89` ratio | +1.67 | +0.00138 | 2/3 | [−0.0028, +0.0058] | +0.00559 | −0.13 | **−1.61** | −0.00188 | −0.00018 | +0.00019 | 45 % | 73 % | 36 % |
| `own89` rank | +0.36 | −0.00061 | 1/3 | [−0.0047, +0.0036] | −0.00091 | +0.29 | −1.62 | −0.00097 | −0.00072 | −0.00011 | 10 % | −32 % | −6 % |
| `own400` ratio | +2.55 | +0.00167 | 3/3 | [−0.0025, +0.0061] | +0.01151 | −0.05 | **−1.52** | −0.00390 | −0.00119 | −0.00061 | 69 % | 88 % | 74 % |
| `own400` rank | +0.18 | −0.00076 | 1/3 | [−0.0048, +0.0034] | −0.00021 | +0.33 | −1.71 | −0.00120 | −0.00114 | −0.00057 | 5 % | −40 % | −1 % |
| `own2000` ratio | +3.32 | +0.00186 | 3/3 | [−0.0023, +0.0063] | +0.01522 | −0.03 | **−1.36** | −0.00505 | −0.00168 | −0.00094 | 90 % | 98 % | 97 % |
| `own2000` rank | −0.03 | −0.00108 | 0/3 | [−0.0051, +0.0031] | −0.00018 | +0.29 | −1.30 | −0.00144 | −0.00139 | −0.00065 | −1 % | −56 % | −1 % |
| `supp` ratio | +2.71 | +0.00108 | 2/3 | [−0.0028, +0.0052] | +0.01459 | −0.02 | +0.20 | −0.00468 | −0.00192 | −0.00131 | 73 % | 56 % | 93 % |
| `supp` rank | −0.20 | −0.00260 | 0/3 | [−0.0061, +0.0011] | −0.00009 | +0.28 | −0.02 | −0.00251 | −0.00244 | −0.00072 | −5 % | −137 % | −1 % |
| `all` norm (R2a) | +2.49 | +0.00191 | 3/3 | [−0.0023, +0.0062] | +0.01302 | 0.00 | +0.22 | −0.00389 | −0.00081 | −0.00025 | 68 % | 100 % | 83 % |
| **`all` rank (R2b)** | −0.22 | −0.00249 | 0/3 | [−0.0062, +0.0014] | +0.00087 | +0.23 | +0.20 | −0.00271 | −0.00254 | −0.00090 | −6 % | −131 % | 6 % |
| `all` ratio | +3.28 | +0.00188 | 3/3 | [−0.0023, +0.0062] | +0.01644 | +0.02 | +0.12 | −0.00504 | −0.00150 | −0.00076 | 89 % | 99 % | 105 % |

**Reading.** The three columns on the right are the whole result. Every renderer keeps at most as much of the `pds`/cosine gain as it
keeps of the `nmae` cost, with one exception — the oracle `exact` cell, which keeps 44 % of the gain for 20 % of the cost and is still
negative at the pessimistic factors (−0.00068) and worth +0.00032 at the point factors, because the union's direction legs at the
pessimistic factors are only +0.0019 in total before the amplitude leg and +0.0015 after it (`pds` +0.00049 at k 0.7, `mse` +0.00126
from the cosine, DE-three +0.00016, amplitude −0.00036) against a −0.00652 identity leg: even the +0.0014 "identity at zero" ceiling 258 priced (§6 there) is under the gate, and no
rendering reaches it. The `rank` modes — the R2b reading of "the union's direction with the incumbent's amplitude structure" — are
exactly `nmae`-neutral (−0.0002 … +0.0009) and lose 30–140 % of the `pds` gain on 0–1 of 3 folds: the union's order and signs
carry nothing once its magnitudes are gone. The `ownK ratio` cells, the only ones that keep more than half of the `pds` gain while
cutting the cost by more than a quarter, move the delivered amplitude by −1.3 … −1.6 % (the union's own −1.05 % plus the cut) and so
belong to the amplitude-changing class (gate +0.0089); they are −0.0002 … −0.0051 pessimistic anyway.

### 3.3 R3 — the cheap subset and the other doses (six renderers each; `work/tables274.md` T2 for every row)

| base | un-rendered PESS | best rendered cell | PESS | Δ`pds` | Δ`nmae` | agg % | amp % | note |
|---|---:|---|---:|---:|---:|---:|---:|---|
| stack3 g=0.35 | +0.00109 | `own89 ratio` | +0.00097 | +0.00101 (2/3) | +0.00020 | +1.81 | +0.28 | every renderer below the parent; `supp rank` +0.00060 with Δ`nmae` −0.00022 |
| stack3 g=0.28 | +0.00003 | `supp rank` | +0.00025 | +0.00108 (3/3) | +0.00001 | −0.06 | +0.03 | the one family where a renderer edges above its parent, by 0.0002, at 0.13× the gate |
| stack6 g=0.28 | −0.00490 | `wald89 ratio` (oracle) | −0.00067 | +0.00110 (2/3) | +0.00162 | −0.24 | −0.37 | best buildable `own89 ratio` −0.00146 (amp −1.08 %); `all rank` −0.00221 with Δ`pds` −0.00124 |
| stack6 g=0.20 | −0.00315 | `wald89 ratio` (oracle) | −0.00134 | +0.00166 (2/3) | +0.00168 | −1.89 | +0.31 | best buildable `all norm` −0.00141; `all rank` −0.00240 |

### 3.4 R4 — 259's construction B rendered (the pair kept at 0.28 + `union2`)

| cell | agg % | Δ`pds` | + | bootstrap | Δ`nmae` 4f | e | amp % | **PESS** | POINT | OPT | ret3 |
|---|---:|---:|:-:|---|---:|---:|---:|---:|---:|---:|---:|
| B3 g₂=0.35 (un-rendered) | +4.95 | +0.00193 | 2/3 | [−0.0013, +0.0054] | +0.00120 | +0.26 | −0.24 | **+0.00184** | +0.00331 | +0.00442 | −45 % |
| B3 g₂=0.35 \| `own400` ratio (**best buildable**) | +4.30 | +0.00183 | 2/3 | [−0.0014, +0.0053] | +0.00102 | +0.26 | −0.22 | **+0.00164** | +0.00302 | +0.00399 | −52 % |
| B3 g₂=0.35 \| `all` norm | +3.71 | +0.00192 | 2/3 | [−0.0013, +0.0054] | +0.00123 | +0.28 | +0.35 | +0.00155 | +0.00302 | +0.00385 | −45 % |
| B3 g₂=0.35 \| `all` rank | +3.00 | +0.00162 | 2/3 | [−0.0015, +0.0051] | +0.00061 | +0.42 | +0.28 | +0.00143 | +0.00260 | +0.00328 | −73 % |
| B3 g₂=0.35 \| `supp` rank | +2.74 | +0.00142 | 2/3 | [−0.0016, +0.0048] | −0.00036 | +1.06 | −0.04 | +0.00142 | +0.00247 | +0.00314 | −87 % |
| B3 g₂=0.35 \| `wald89` exact (oracle) | +3.25 | +0.00053 | 2/3 | [−0.0027, +0.0041] | +0.00018 | +0.42 | +0.20 | +0.00139 | +0.00179 | +0.00252 | −418 % |
| B3 g₂=0.45 \| `all` norm / `all` rank | +3.94 / +3.05 | +0.00119 / +0.00081 | 2/3 | — | +0.00132 / +0.00076 | | +0.52 / +0.41 | +0.00147 / +0.00122 | +0.00247 / +0.00189 | | |
| B4+igvf g₂=0.35 (un-rendered) | +5.82 | +0.00237 | 3/3 | [−0.0015, +0.0066] | +0.01536 | −0.01 | −0.48 | **−0.00370** | +0.00004 | +0.00135 | −50 % |
| B4+igvf g₂=0.35 \| `all` rank | +1.09 | −0.00208 | 0/3 | [−0.0057, +0.0017] | +0.00061 | +0.40 | +0.21 | −0.00166 | −0.00157 | +0.00004 | 235 % |
| B4+igvf g₂=0.35 \| `own400` ratio | +4.43 | +0.00210 | 3/3 | [−0.0018, +0.0062] | +0.01085 | −0.03 | −0.78 | −0.00253 | +0.00039 | +0.00138 | −68 % |
| B4+igvf g₂=0.35 \| `wald89` ratio (oracle) | +2.78 | +0.00055 | 2/3 | [−0.0033, +0.0047] | +0.00140 | +0.61 | −0.61 | +0.00043 | +0.00101 | +0.00163 | −495 % |
| B4+igvf g₂=0.35 \| `wald89` exact (oracle) | +3.41 | +0.00126 | 2/3 | [−0.0026, +0.0053] | +0.00281 | — | −0.61 | +0.00027 | +0.00151 | +0.00228 | −163 % |

The direct test of the mandate's premise is the `igvf` row: 258/259 attributed the union's cost to `igvf_cardiomyocyte`; with `igvf`
back in the second residual the cell pays +0.01536 raw and reads −0.00370; the best rendering (the oracle exact pool) keeps 53 % of
its `pds` and 59 % of its cosine for 18 % of the cost and reads **+0.00027** — 0.14× the gate, retention −163 %; the buildable `all
rank` neutralises the cost (+0.00061) and loses the `pds` (−0.00208, 0/3). B3 without `igvf` was `nmae`-cheap to begin with (+0.0012),
so there was almost nothing for a renderer to remove, and each rendering only costs cosine.

### 3.5 Controls

The four identity cells reproduce `cd4svd5` at 0.0e+00. The 24 random-direction cells (the union's per-row contribution permuted within
its support, or deranged across carried rows, then rendered by the same set × mode) read Δ`pds` −0.0082 … +0.0001 (max +0.00012, on
B3 g₂=0.35 `all rank`, 6 % of that cell's gain) and pessimistic −0.0096 … −0.0009: no renderer manufactures a gain from a random
direction of the union's energy, and every real cell's controls are below half of it (`controls_ok` True wherever the cell is
positive). The un-rendered random contributions read −0.0028 … −0.0065 on the six-source union: replacing the union's directions by
random ones of the same energy costs 1.5–3× what the union gains, so the union's directions are real — they are just not separable
from its magnitudes.

---

## 4. The mechanism — where the union's gain sits on the row (`mech274.py`, T8; frame A, a23 render, scored columns, means over the union's carried rows, three cancer folds)

| arm | own cos | Δ own cos | Δ ip on `wald89` (89 genes) | on `own89` | on `own400` | off `own400` | pred. energy on `wald89` | L1 on `wald89` / inc | sign acc on `wald89` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `cd4svd5` | +0.03889 | — | — | — | — | — | 0.0358 | 1.000 | 0.5749 |
| stack6 g=0.35 | +0.04034 | **+0.00146** | **+0.00079** | +0.00059 | +0.00052 | +0.00094 | 0.0409 | **1.044** | **0.5772** |
| \| `wald89` ratio (oracle) | +0.03894 | +0.00006 | −0.00079 | −0.00068 | −0.00079 | +0.00084 | 0.0333 | 0.965 | 0.5774 |
| \| `wald89` exact (oracle) | +0.03940 | +0.00052 | **−0.00029** | −0.00030 | −0.00038 | +0.00089 | 0.0356 | 0.999 | 0.5772 |
| \| `own89` ratio | +0.03968 | +0.00079 | +0.00021 | −0.00019 | −0.00022 | +0.00101 | 0.0379 | 1.030 | 0.5773 |
| \| `all` norm | +0.04038 | +0.00149 | +0.00083 | +0.00062 | +0.00055 | +0.00094 | 0.0413 | 1.017 | 0.5773 |
| \| `all` rank | +0.03868 | −0.00020 | −0.00081 | −0.00130 | −0.00134 | +0.00114 | 0.0335 | 1.014 | 0.5772 |

1. **The truth's Wilcoxon pool is where the transfer cosine lives.** The 89 `wald89` genes (of ~6,000 scored) carry 6.3 % of the truth's
   unit-row energy and 3.6 % of the incumbent's prediction energy, yet ~40–47 % of the incumbent's own cosine (ip on `wald89` 0.0127 /
   0.0101 of own cos 0.0311 / 0.0262 on HCT116 / HEK293T). This is the identity's pool by construction (the reference's largest, best-
   expressed effects), and it is also where any correctly-signed prediction earns its cosine.
2. **54 % of the union's own-cosine gain sits on those 89 genes** (+0.00079 of +0.00146), 40 % on the incumbent's own top-89, 36 % on
   its top-400. The union puts it there as MAGNITUDE: L1 on the pool 1.044× the incumbent's (per-row ratio median 0.978–0.985 after
   the ratio renderer's correction, i.e. the union raises the pool on most carried rows, up to 2.6–2.8×), while sign accuracy on the pool
   moves 0.5749 → 0.5772. At EQUAL magnitude (the exact oracle) the union's direction inside the pool contributes LESS than the
   incumbent's (Δ ip −0.00029): the union has no sign information on the DE set that the incumbent lacks; it has more of the same.
3. **That is exactly what the identity charges.** 188's split on the union: ΔU/D +0.00657 (unmatched magnitude up), ΔM/D −0.00020
   (matched magnitude flat), e −0.03 — "magnitudes, not signs" (188 §1.3, 228 §0.4), and now with the direction side measured on the
   same genes: the magnitude the identity charges as unmatched is the magnitude the cosine collects as aligned, because the pool's
   sign accuracy is 0.57 and \|p − r\| on a correctly-signed gene grows once \|p\| passes \|r\| while p·r keeps growing. The two members
   read the same entries with opposite signs; a renderer can only choose the trade.
4. **The off-pool gain is real but small and `pds`-blind.** +0.00084 … +0.00094 of own cosine survives every pool-holding renderer
   (the "off `own400`" column), 58–64 % of the union's cosine gain in cosine units — but `pds` is a threshold at 0.02 own cosine on
   the marginal band (262 H3) and the union's `pds` is the anti-aligned group's (Δ +0.0235 on 46–52 targets per fold vs −0.0029 on
   the 155–209 marginal targets; T7). The renderers keep 50–90 % of the anti-aligned gain (+0.0115 … +0.0214) and lose more on the
   marginal band (−0.0039 … −0.0065): targets crossing the knee upward fall from 13 to 8–10 while the downward crossings stay at 7–10.
   The oracle exact cell keeps +0.0214 on the anti-aligned group and reads Δ`pds` +0.00084 with the interval [−0.0032, +0.0051].
5. **The rank transplant destroys the direction because the union's ORDER is not its direction.** Rank agreement between \|X_arm\| and
   \|X_ref\| on the whole row is ρ 0.96; the transplant changes each entry by the local spectrum step and every marginal target's cosine
   by noise of that size — 8 up / 9 down across the knee, Δ`pds` −0.0025, on 0/3 folds; its `nmae` neutrality (+0.0009) is the
   spectrum's, not a property of the union's signs.

---

## 5. The best cells in full

### 5.1 Best buildable: B3 g₂=0.35 \| `own400` ratio (259's construction B, the incumbent's top-400 magnitudes held per row)

| quantity | value |
|---|---|
| PLOLO vs `cd4svd5` | rel `cos_agg` **+4.297 %** (2.90 / 5.12 / 4.88), rel macro +3.611 % |
| Δ`pds` ABC | **+0.001833** (k562 +0.00292, hct116 +0.00332, hek293t **−0.00074**; 2/3; pooled bootstrap [−0.00138, +0.00533], P(≤ 0) 0.140; per fold [−0.0013, +0.0073] / [−0.0034, +0.0117] / [−0.0058, +0.0041]); frame A +0.00218 |
| `nmae` identity | Δ`wald89` **+0.001015** (3 cancer +0.001464; 4f bootstrap [−0.0005, +0.0033]); `lfc89` −0.00016, `h89` +0.00233; ΔU +233, ΔM +61, e 0.262; Δ precision89 +0.00204; Δq_w +0.00049 |
| delivered amplitude | production effect-level −0.225 %, calibrated **−0.221 %** (within 0.5 %); preserve 0.8964 vs 0.9123 |
| renderer | 76 / 81 / 81 changed rows (union 2's carried rows); per-row scale median 0.991 [0.83, 1.01]; rank agreement on S 0.76–0.82 |
| 262 groups | anti-aligned +0.0131, marginal +0.0000, saturated −0.0008; crossing 0.02: up 10 / down 4 |
| `kolf2` (apart) | rel `cos_agg` +6.53 %, Δ`pds` −0.00271 |
| retention | k=1 18.7 %, **k=3 −52 %**, k=5 −103 %; top-3: k562 LENG1 +0.20, TRIM35 +0.19, MAP2K2 +0.17; hct116 **IFNAR2 +0.88**, TARBP2 +0.22, SH2D3A +0.17 — 259's four targets again |
| controls | `own400 ratio` on the permuted union 2: Δ`pds` −0.00015, pessimistic −0.00110 (below half) |

| scenario | k_cos | k_pds | DE3 | `nmae` × | Δcos | raw `mse` | `pds` | `mse` | `nmae` | `reach` | `fid` | `jac` | (amp) | **TOTAL** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PESSIMISTIC | 0.94 | 0.70 | 0.0425 | 1.5 | +0.004260 | 0.81376 | +0.00047 | +0.00141 | −0.00042 | +0.00013 | +0.00005 | +0.00001 | −0.00008 | **+0.00164** |
| POINT | 0.96 | 2.48 | 0.0425 | 1.0 | +0.004350 | 0.81358 | +0.00168 | +0.00144 | −0.00028 | +0.00013 | +0.00005 | +0.00001 | −0.00008 | +0.00302 |
| OPTIMISTIC | 1.00 | 2.48 | 0.24 | 1.0 | +0.004532 | 0.81322 | +0.00168 | +0.00151 | −0.00028 | +0.00078 | +0.00027 | +0.00003 | −0.00008 | +0.00399 |

258's conservative (frame A) +0.00177, point (frame A, k_pds 3.79) +0.00526. Conditions: pessimistic ≥ +0.0019 **NO** (0.86×);
`pds` 3/3 with the interval above 0 **NO**; retention ≥ 50 % **NO**; amplitude within 0.5 % yes; controls below half yes. This is
259's cell (+0.00184 / +0.00331 / +0.00442, retention −45 %, `pds` 2/3) minus 0.0002–0.0004 at every scenario: the renderer removed
+0.00018 of a +0.00120 `nmae` cost and 0.65 pp of cosine.

### 5.2 Best oracle: B3 g₂=0.35 \| `wald89` exact — +0.00139 pessimistic, 0.73× the gate, `pds` +0.00053 [−0.0027, +0.0041], retention −418 %, and not buildable (the truth's DE set). The six-source oracle exact cell: −0.00068 / +0.00032 / +0.00069 (§3.2); its per-fold Δ`pds` −0.0007 / +0.0013 / +0.0019.

---

## 6. The production change it would need — named, NOT made

One stage, one place: `scripts/make_submission_pm_cd4.py::apply_frozen_cd4_residual`, between the construction of the merged
`added` residual (lines 2638–2645) and the call `EXP062.add_post_whiten_residual(whitened, added, applied, operator_columns)` (line 2647) —
equivalently inside `E062.add_post_whiten_residual` between line 580 (`candidate[available] += fraction * mean_effect`) and line 581
(the preserve). The stage would take the candidate of the arm AND the candidate of the reference composition (the pair at 0.28 — so
`prepare_h1_union` would run twice, once for `VCC_UNION_SOURCES` and once for the pair), apply `render274.render(X_arm, X_ref, set,
mode, K)` on the union's carried rows, and hand the result to the unchanged preserve; a flag `VCC_UNION_RENDER=<set>:<mode>[:K]`,
unset = byte-inert (the KA row a build would need: flag unset reproduces `aba27758…`). Its audit block would carry the changed-row count,
the per-row scale quantiles and the rank agreement, as `work/fold-*.json` `renderer` does. **Nothing here is a candidate, so nothing was
written into the builder.**

---

## 7. What the screen settles

1. **The union's direction is not separable from its magnitude on the DE set.** 54 % of its cosine gain and all of its identity cost
   sit on the same 89 genes per row; at equal magnitude on those genes the union's direction is worth −0.00029 of own cosine
   (less than the incumbent's), and its sign accuracy there is the incumbent's + 0.2 pp. The premise of BR237 §3 / R258 §6 — an
   `nmae`-neutral rendering that keeps the +0.0034 direction legs — is false as stated: the neutral renderings keep ≤ 44 % of the
   cosine (oracle) and lose the `pds` (buildable), and the "identity at zero" ceiling (+0.0014 pessimistic) was under the gate anyway.
2. **The R2 reading ("the union's unit vector with the incumbent's amplitude structure") is measured and dead**: `all norm` keeps the
   whole `pds` gain and 83 % of the cost (it only fixes the preserve); `all rank` neutralises the cost and loses 131 % of the gain.
3. **The `ownK ratio` renderers are the only buildable ones on the frontier's good side and they are amplitude-changing** (−1.3 … −1.6 %
   delivered), which puts them under the +0.0089 gate; at 36–74 % of the cost kept they are −0.0019 … −0.0051 pessimistic regardless.
4. **On the `nmae`-cheap unions (`stack3`, B3) a renderer only costs**: their identity leg was +0.0006 … +0.0014 raw to begin with
   (e +0.20 … +0.33), so there was ≤ +0.0003 overall to recover and every rendering spends more cosine than that. The best buildable
   cell in this experiment is 259's cell with a renderer subtracted from it.
5. **The `igvf` attribution and its remedy are both measured**: with `igvf` back, the renderers can neutralise its +0.0154 raw cost
   (`all rank` +0.0006; the oracle exact pool +0.0028) but not keep its `pds` (−0.0021 buildable; +0.0013 [−0.0026, +0.0053] oracle);
   the oracle reads +0.00027, 0.14× the gate. The sources' magnitudes on the DE set are their information.
6. **The controls behave**: identity at 0.0e+00, random directions rendered the same way at −0.0009 … −0.0096 pessimistic — the
   renderer is not a hidden source of `pds`, and the union's real directions are worth 1.5–3× a random direction of the same energy.

This closes the last measured positive direction on the `cd4svd5` base: 258 closed the union at its own magnitudes (31 cells), 259
closed it without `igvf` in both constructions (21 cells), 261 closed its deflation (10 cells), and this experiment closes its
rendering (39 buildable + 11 oracle cells, one production-shaped stage). What remains of BR237's +0.0034 is the arithmetic that an
`nmae`-free direction gain would be worth at the pessimistic factors — +0.0015 at most — which is under the gate before any renderer
touches it. **A six-source union earns a slot only if a source with sign information on the recipient's DE set exists; none of the
five on disk has it (sign accuracy on the pool 0.5772 vs 0.5749), and that is a source question, not a rendering question.**

---

## 8. Honest limits

1. **Per-target `pds` on PLOLO is noise-dominated** (267 / 273: 150–180-cell references, split-half residual cosine ≈ +0.03 … +0.05):
   the six-source union's own Δ`pds` +0.00191 has its paired-bootstrap interval across zero (P(≤ 0) 0.19), so every "fraction of the
   `pds` gain kept" in §3.2 is a ratio of two uncertain numbers; the cosine fractions (deterministic aggregates) and the identity's
   leg (3.1 bootstrap-SDs from zero) carry the closure. The board realised the union's `pds` at +0.00666 raw on `h1union` (BR237),
   so the gain is real; the question was its separability, and the mechanism (§4) is measured on the deterministic quantities.
2. **The identity is a proxy** (fixed top-89 pools, four folds including `kolf2`; validated at scale as a conservative bound, ×1.5,
   BR237); the 3-cancer-fold variant is within 20 % and same-signed on every cell here (T1 / T2 columns).
3. **The oracle cells use the truth's DE set** and are ceilings only; the exact-pool cell also holds the pool L1 at the post-preserve
   level, which no buildable renderer can compute (the preserve depends on the whole matrix). Their amplitude is read on the fold
   frame (243's A_c on the fold's targets), not the production frame.
4. **The class factors** (k_cos 0.94–1.00, k_pds 0.7 / 2.48, DE-three 0.0425 / 0.24) are the two board points' realised values on
   arms of this class; a rendered union is a new sub-class whose transfer has never been measured on the board. The verdict does not
   depend on them: at k_pds 2.48 and k_cos 1.00 (optimistic) the best six-source rendering is +0.00069 and the best buildable cell of
   any family +0.00442 is 259's cell itself un-rendered (+0.00442 vs +0.00399 rendered).
5. **Not tried**: a per-context renderer inside the emitter (`profile_from_effect`), a renderer that adds back the off-pool direction
   at a larger dose, and a learned magnitude map on the union's contribution only. The first cannot see the truth either and acts on
   the same entries; the second is 258's dose axis (the off-pool cosine is 58–64 % of the union's, its `pds` is the anti-aligned group's,
   and the dose ladder was negative at every g); the third is 188's family restricted to the union — 188 measured that our predicted
   magnitude carries no information about the true magnitude (+9.5 % of median \|r\| over a 36× range of \|h\|), and the union's pool
   magnitudes here are sign-uninformative (item 2 of §4). None was priced; each is named so it is not re-proposed as new.
6. **Retention** is reported by 222's definition and is meaningless where Δ`pds` ≈ 0 (the ratio explodes); it is shown only for
   \|Δ`pds`\| ≥ 0.0005 in the tables and was never the deciding condition (the pessimistic gate is failed by every cell first).
7. **No package, no board member measured**; nothing clears, so 248/249's post-build falsifiers were not reached.

---

## 9. Cost and memory (every phase alone in its process; nothing ≥ 15 GiB was live on the machine)

| script / phase | runtime | peak RSS |
|---|---:|---:|
| `render274.py --fold k562_gwps` (269 targets, 91 arms) | 102 s | 5.66 GiB |
| `--fold orion_hct116` (290) | 120 s | 5.81 GiB |
| `--fold orion_hek293t` (294) | 121 s | 6.00 GiB |
| `--fold kolf2` (273) | 112 s | 5.87 GiB |
| `--prod` (production-frame reconstruction, 80 compositions, 3 renders each) | 63 s | 5.41 GiB |
| `--price` (bootstraps included) | 14 s | 0.85 GiB |
| `mech274.py --fold <held>` × 3 | 7 s each | 5.18–5.43 GiB |

Wall clock 65 min from launch to this file (two passes; pass 1 under `work/pass1/`, identical on every shared number).

## 10. Files

| file | what |
|---|---|
| `render274.py` | the instrument: `--fold` (compose through the builder, recover the pre-preserve compositions, render, both functionals, 216 KA rows per fold), `--prod` (production-frame amplitudes of every buildable cell, 9 KA rows), `--price` (139 cross-fold KA rows, KA-1/KA-2, the screen, bootstraps, retention, 262 groups, the verdict), `--list` |
| `mech274.py` | the mechanism (§4, T8): inner product / energy / L1 / sign accuracy on the truth pool and the incumbent's own sets, per target |
| `tables274.py` → `work/tables274.md` | T0–T8: every cell, six members priced three ways (T5), per-fold bootstraps (T6), 262 groups (T7), the mechanism (T8) |
| `run274.sh` | the folds one process at a time, gated as 258 |
| `work/fold-<held>.{json,npz,stdout,stderr}` | per-fold rows (91 arms), renderer diagnostics, per-target `pds` / cosine / U-M-D / pool L1 / changed-row masks, 262's base cosine aligned |
| `work/prod274.{json,stdout}` | production-frame terms, renderer diagnostics on the 300-target panel, calibration |
| `work/price274.{json,log}` | every number in this file; `work/pass1/` the first pass (87 arms) |
| `work/mech274-<fold>.json`, `work/mech274-summary.json` | the mechanism arrays and their carried-row means |
| `work/plolo/` | 258's fold caches, byte copies (`SHA256SUMS.txt`) |

---

## 11. Verdict

**DOES NOT CLEAR.** No rendering of the six-source union's directions on top of `cd4svd5` — pool-restricted (R1), whole-row (R2), on
the cheap subset or the other doses (R3), on the second-residual construction (R4), buildable or oracle — clears +0.0019 at the
coherent pessimistic factors; 0/39 buildable and 0/11 oracle cells, 0 cells on all conditions. The union's direction cannot be
rendered without its magnitude cost: on the scored DE set they are the same entries. No `PREREGISTRATION-DRAFT.md`. No model name, no
description, nothing built, nothing submitted.
