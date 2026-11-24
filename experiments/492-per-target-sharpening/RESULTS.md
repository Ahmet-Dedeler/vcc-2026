# 492 — Per-target sharpening exponent by the 410 key tercile: GO as a candidate (preregistered rule passes), worth ~+0.001 … +0.003

2026-09-26 23:15Z, third session. Analysis only (no build, no builder edit, no `vcc`, no `cell-eval2`). Preregistered in
`PREREGISTRATION.md` (commit `9086df2a`, before any arm number). Scripts: `arms492.py` (PLOLO, one fold per process),
`abc492.py` (A/B/C energies, derived from 473's `abc473.py`), `analyse492.py` (473's pricer with 492's arms), `k492.py`
(the ballast's nmae cost, calibrated on the board's q step).

## Answer
- **P2 prices +0.0034 overall vs 446b** (k 1, both arms re-levelled). P2 sharpens the low-key third of targets at
  q 2.5, the middle third at 2.0, and keeps the high-key third at 1.75.
  - **3/3 cancer folds positive**, each fold's target-bootstrap interval above zero. kolf2 +0.0011.
  - Pooled relative cos +0.95 % [+0.67 %, +1.27 %].
  - A/B/C render energy ×0.99 in every context (no re-level needed beyond ~1 %), and no per-target ratio above ×1.00.
- **The gain is targeted.** The same q values dealt to random targets (the control) keep the pds gain but lose direction:
  pooled cos −0.56 % [−2.3 %, +1.0 %], kolf2 −0.0019, overall +0.0023 (below P2 on 3/3 folds).
- **Why:** low-key targets are the ones whose donors disagree. Their rows are the noisiest, and sharpening extracts
  their peaks. High-key rows are already saturated at q 1.75 (420's arrays: +0.002 pds from q 1.75 to 2 on the high
  tercile, vs +0.006 … +0.008 on the low).
- **Board-calibrated expectation: +0.0013 … +0.0032 overall.**
  - The upper end assumes k 1. The lower end uses 417a's sharpening transfer (k 0.65 on mse, k_pds 0.42).
  - Both ends include the ballast's nmae cost, **−0.0002** (below).
  - Global q 2 (G2) prices +0.0029 here but kolf2 −0.0003. It moves A/B/C render energy to ×0.89, so it needs a large
    re-level. P2 is the better way to spend the same sharpening.
- **Hand-off:** a builder knob (the ballast bulk sharpened with a per-target q from the pinned K3 artifact's tercile) and
  a build are the coordinator's. The October day-of path computes the key per panel, so the tercile is panel-relative
  and ports as is.

## Known answers (all exact, before any arm number)
| row | expected | measured |
|---|---|---|
| per-target path at q 1.75 for every row vs `comp420` `stack\|K3+sh1.75` (designed, realised 16) | 1e-12 | \|d\| 0, 4/4 folds |
| per-target path at q 2 for every row vs `comp420` `stack\|K3+sh2` | 1e-12 | \|d\| 0, 4/4 folds |
| KA-PP (compose_arm == compose_pp), K3 weights == 418's W3 | bit-identical | PASS 4/4 |
| A/B/C: 473's KA-a / KA-c (input shas; emulated ballast K median 6, max 16) | as 473 | PASS (KA-b is 473's documented ~2 % level offset; only ratios are used) |
| emulated K medians for 433k / 446b vs the build metas (4 / 6) | — | 3 / 6 (the calibration below uses mean K, same route for every arm) |

## PLOLO (realised, KMAX 16) and the price (446b-anchored cliff model, re-levelled, k 1; 95 % target CI per fold)
| arm | q (low, mid, high key tercile) | k562 | hct116 | hek293t | cancer mean (k 1 / k 1.407) | kolf2 | pooled rel cos |
|---|---|---|---|---|---|---|---|
| P1 | 2.0, 1.75, 1.75 | +0.0009 [+0.0002, +0.0017] | +0.0012 [+0.0005, +0.0021] | +0.0007 [−0.0002, +0.0016] | +0.0009 / +0.0010 | +0.0007 | +0.22 % |
| **P2** | **2.5, 2.0, 1.75** | **+0.0024 [+0.0004, +0.0045]** | **+0.0042 [+0.0021, +0.0065]** | **+0.0035 [+0.0009, +0.0063]** | **+0.0034 / +0.0036** | **+0.0011** | **+0.95 % [+0.67, +1.27]** |
| P3 | 3.0, 2.25, 1.75 | +0.0030 [−0.0001, +0.0063] | +0.0061 [+0.0031, +0.0095] | +0.0044 [+0.0005, +0.0084] | +0.0045 / +0.0048 | +0.0008 | +1.15 % |
| G2 | 2.0 everywhere | +0.0018 | +0.0031 | +0.0039 | +0.0029 / +0.0032 | −0.0003 | +1.06 % [+0.17, +1.97] |
| CTRL | P2's q, permuted | +0.0018 | +0.0016 | +0.0034 | +0.0023 / +0.0021 | −0.0019 | −0.56 % [−2.3, +1.0] |

- P2's board-emulated pds: +0.0044 / +0.0094 / +0.0082 (cancer folds).
- P2's raw mse at k 1, re-levelled: −0.0037 (cancer mean).
- A/B/C render (realised, arm / base), P2: energy ×0.990 / ×0.990 / ×0.989, max per-target ×1.00, mass share unchanged.
  The ballast-capped share is 45 / 43 / 46 % vs 41 %.

## The ballast's nmae cost (`k492.py`), which the PLOLO instrument does not see
- The board's own q step (433k → 446b, q 1.5 → 1.75) cost nmae +0.001346 raw, as the deep ballast cells entered the
  per-cell means. It moved the emulated mean ballast count K from 7.28 to 8.65.
- Scaled by each arm's mean ΔK:
  - P2 +0.72 → nmae **+0.0007 raw, −0.0002 overall**;
  - P3 −0.0004; G2 −0.0002; P1 0.
- P2's extra sharpening lands mostly on low-key rows. Their bulk is small, so their ballast need is small: K 2.4 → 3.0 on
  the low tercile and 8.7 → 10.3 on the middle, with the high tercile unchanged at 14.8.

## GO rule (preregistered, for P2) → all four pass
| item | rule | P2 | pass |
|---|---|---|---|
| 1 | cancer mean ≥ +0.0015 (k 1), ≥ 2/3 folds positive, kolf2 ≥ 0 | +0.0034, 3/3, kolf2 +0.0011 | yes |
| 2 | target-bootstrap lower bound of Δpds_bd > 0 | per-fold lower −0.0008 / +0.0041 / +0.0012, mean +0.0015 | yes |
| 3 | CTRL below P2 on ≥ 2/3 cancer folds | 3/3 (hek293t by only 0.0001) | yes |
| 4 | A/B/C: max per-target ratio ≤ 1.5, capped share ≤ base + 10 points | ×1.00; +2 … +5 points | yes |

## Caveats
- The transfer class is sharpening's, unmeasured for a per-target form. The board has measured global q steps:
  417a q 1.25 at k_sh 0.42–0.65; 433k and 446b at the calibrated rate.
- P3 prices higher on cancer folds but lower on kolf2. The q curve saturates, so P2, the preregistered primary, is
  the arm.
- The tercile uses the K3 key, whose floor holds 82 of 300 A/B/C targets. The low tercile is essentially the floor.
- Item 3's margin on hek293t is 0.0001: the targeting's value is carried by direction (cos), not pds.

## Knob definition (answer to the coordinator, 23:30Z; what the instrument actually did)
- **Key array:** the per-target weight AFTER the clip.
  - On A/B/C: `w` of the pinned `pert_scale_410_lo0.3_hi3.npz` (`wshare_lo0.3_hi3`, sha `d7f76a22…`), clipped to
    [0.3, 3], with non-carried targets at 1.0.
  - On PLOLO: the fold's own `weights(key410, carried410, 0.3, 3.0)`. Not the raw key before the clip.
- **Cut points:** `edges = np.quantile(w, [1/3, 2/3])` (numpy's default linear interpolation), then
  `tercile = np.digitize(w, edges)`: low if w < q1, mid if q1 ≤ w < q2, high if w ≥ q2.
- **On A/B/C this gives exactly 100 / 100 / 100** (q1 0.5234, q2 1.4605, no ties at either cut).
  - All 82 floor targets (w 0.3) are low, together with the 18 lowest above the floor.
  - Rank-based thirds give the identical split here.
  - PLOLO folds: 90/89/90, 97/96/97, 98/98/98, 91/91/91.
- **October guard (not exercised on A/B/C):** if ≥ 1/3 of a panel sits at the floor, q1 equals the floor value and
  `digitize` moves every floor target to MID, leaving LOW empty.
  - The knob should assign w ≤ floor to LOW explicitly: `low = (w < q1) | (w <= 0.3 + 1e-12)`.
  - That is what "sharpen the least reliable rows hardest" means, and it is identical on A/B/C.
- The q values per tercile are (2.5, 2.0, 1.75). Each row keeps its own gated non-panel norm after sharpening, exactly as
  the global knob does, and the panel columns are untouched.

### Knob definition in LOAD mode (the coordinator's r-tercile, checked 04:45Z on 446bj2's dump; the knob's first KA)
- In LOAD mode the arm's own VCC_PERT_SCALE stays 410's key, so the key tercile is read from the loaded dump itself:
  r_t = ‖bulk_t‖ / ‖carriers_t‖ over non-panel columns (= c · w_t; the global c cancels in the terciles).
- **On `dump_h1x2_cd4x15_lodo05.npz` (sha `00ccfcba…`) the r-terciles reproduce the LODO-w terciles exactly.** Counts
  100 / 100 / 100 with an identical assignment, guarded and unguarded.
  - The floor guard "r within 1e-9 relative of min(r) → LOW" selects the same 97 floor targets as w ≤ 0.5.
  - r / w is constant to 1.6e-14 relative (c = 0.5771677758632). No carrier row is zero (min ‖carriers_t‖ 2.47).
  - The nearest target to a cut point is 0.06 % (relative) away, so rounding cannot move a target across an edge.
- `stercile` = ‖carriers_t‖ over non-panel columns: 100 / 100 / 100 (unchanged from 492d).
- **October caveat (493g).** Non-carried targets (fewer than 2 donors) have w = 1 exactly, i.e. r = c up to ~1e-14
  rounding. When they are numerous (efp: 140 of ~300 carried) the tie block straddles a cut point.
  - In w-space the whole block lands in one tercile, and on efp hct116 / hek293t the MID tercile is empty.
  - In r-space the rounding noise would split the block arbitrarily.
  - A key tercile used in October needs an explicit rule for the tie block (being measured in 493h), or the strength
    tercile, which has no ties.

### The tie rule for non-carried targets (493h, 04:50Z): what the key tercile does with the w = 1 block
Measured on 438's October-like panels under the October LODO key (`493/comp493h.py`; known answer = 493g's P2g row,
|d| 0 on 9/9). Only efp has the tie: 153–156 of ~300 targets are non-carried on two of its folds. efr and real carry all
300, so every rule coincides there. Direction vs the q 1.75 base on efp (3 folds):

| rule (floor → LOW throughout) | efp rel cos (folds +) | ΔK | on A/B/C |
|---|---|---|---|
| plain (cut points over all targets) | +0.48 % (2/3), MID empty on 2 folds | +0.83 | as priced (100 / 100 / 100) |
| **tieMID**: plain cut points, non-carried → MID | +0.63 % (3/3) | +1.60 | **identical to the priced arm** (the 6 non-carried are MID already) |
| tieLOW: non-carried → LOW | +0.69 % (2/3) | +2.88 (most ballast) | — |
| **carried-only**: cut points over carried targets, non-carried → MID | **+0.84 % (3/3)** | +1.61 | moves 4 targets (98 / 104 / 98) |
| S_P2 (strength tercile, for reference) | **+1.32 % (3/3)** | +1.34 | 100 / 100 / 100 |

- **479 on A/B/C: use the plain rule as priced.** Its 6 non-carried targets are MID by value and far from the cut
  points, so the tie rule is moot there.
- **October:**
  - Strength (`stercile`) is the default candidate: best or tied on every panel, and tie-free.
  - If the key tercile is used, apply carried-only cut points with non-carried → MID.
  - Read `carried` from the dump's own key artifact (its env records `VCC_PERT_SCALE`; the builder's registry resolves
    the npz and its `carried` array). Do not infer it from r: r-space has a second tie block at the cap (w = 3, 49
    targets on A/B/C) that no rule can tell apart from the w = 1 block without the artifact.

## Addendum 492b (23:15Z, EXPLORATORY: written after the preregistered read, not a selection rule)
`arms492.py <fold> --b`, `analyse492.py --b` (`work/arms492b-*`, `work/analyse492b.json`). Same instrument; all KAs exact
again. Price vs 446b (k 1, re-levelled; "cal." = 417a's k 0.65 / k_pds 0.42):

| arm | q (low, mid, high) | cancer mean k 1 (folds +) | kolf2 | cal. | A/B/C energy | pooled rel cos |
|---|---|---|---|---|---|---|
| L25 | 2.5, 1.75, 1.75 | +0.0020 (3/3) | +0.0013 | +0.0009 | ×1.000 | +0.41 % [+0.24, +0.59] |
| **L30** | **3.0, 1.75, 1.75** | **+0.0023 (3/3)** | **+0.0017** | **+0.0010** | **×1.000** | +0.40 % [+0.18, +0.64] |
| L40 | 4.0, 1.75, 1.75 | +0.0020 (3/3) | +0.0018 | +0.0009 | ×1.000 | +0.27 % |
| L30M20 | 3.0, 2.0, 1.75 | +0.0036 (3/3) | +0.0014 | +0.0017 | ×0.990 | +0.93 % [+0.63, +1.29] |
| P2 (preregistered) | 2.5, 2.0, 1.75 | +0.0034 (3/3) | +0.0011 | +0.0016 | ×0.990 | +0.95 % |
| CTRL L30's q permuted | — | **−0.0011 (1/3)** | −0.0035 | −0.0010 | ×0.89 | −3.5 % |

- **Sharpening only the low-key third (L30) is the clean version.** It gains +0.0023 with no energy change and no ballast
  change (low-key rows are low-amplitude), and it has the best kolf2 of the group.
  - Its control loses outright: dealing q 3.0 to random targets hits high-amplitude rows. That is sharper evidence for
    the targeting than P2's control.
- **The middle tercile at q 2.0 adds ~+0.0013** (P2, L30M20 vs L25, L30). Its cost is the ballast, −0.0002 overall
  (§ above), and a ×0.99 energy that the re-level refunds.
- The build stays **P2** as preregistered. If a minimal-risk variant is preferred (no energy or ballast change,
  kolf2-best), L30 is the alternative, at about 2/3 of P2's cancer-fold gain.

## Addendum 492c (23:25Z, EXPLORATORY): the key floor × per-target sharpening
`arms492.py <fold> --c`, `abc492.py`, `analyse492.py --c` (`work/arms492c-*`, `work/analyse492c.json`).
- New known answer: the builder's `apply_pert_scale` on the stack dump's unkeyed rows with the pinned K3 w reproduces the
  dump's keyed bulk exactly (|d| 0). Re-clipping the artifact's key at 0.3 reproduces its w.
- Other floors re-clip the same key vector. The q tercile stays the K3 one.

| arm | floor, q (low, mid, high) | cancer mean k 1 (folds +) | kolf2 | cal. | Δpds board-emul. / literal | Δraw mse | rel cos |
|---|---|---|---|---|---|---|---|
| F5 | 0.5, 1.75 everywhere | +0.0015 (3/3) | −0.0005 | +0.0007 | +0.0031 / **+0.0000** | −0.0020 | +0.51 % |
| F4L30 | 0.4, (3.0, 1.75, 1.75) | +0.0033 (3/3) | +0.0018 | +0.0015 | +0.0076 / +0.0063 | −0.0030 | +0.76 % |
| F5L30 | 0.5, (3.0, 1.75, 1.75) | +0.0039 (3/3) | +0.0020 | +0.0018 | +0.0087 / +0.0063 | −0.0039 | +1.00 % |
| **F5P2** | **0.5, (2.5, 2.0, 1.75)** | **+0.0051 (3/3)** | +0.0012 | **+0.0024** | +0.0109 / +0.0076 | −0.0062 | +1.58 % |
| F2L30 | 0.2, (3.0, 1.75, 1.75) | +0.0005 (2/3) | +0.0015 | +0.0002 | +0.0015 / +0.0062 | +0.0002 | −0.06 % |

- **The floor stacks with the per-target sharpening almost additively.** F5P2 +0.0051 against P2's +0.0034 plus F5's
  +0.0015.
- A lower floor (0.2) loses what L30 gains.
- The floor raises the floor rows' render energy ×2.7 by design. Context energy moves ≤ 0.5 %.
- **Caveat, as 420 found:** the floor's pds value is the board-emulated β term alone (literal pds +0.0000).
  - That term is a mechanism fitted to one board point (410b − 395).
  - So the floor's measured value is its mse (−0.002 raw) plus a one-point emulation.
- **Sequencing suggestion (the coordinator's call):** the floor is a second, separate variable. Put it on the board only
  after P2/L30 have read, so the two are not confounded.

## Addendum 492d (23:35Z, EXPLORATORY): the same q assigned by response strength instead of the key
`arms492.py <fold> --d`, `analyse492.py --d`. Strength = the unkeyed stack row's non-panel norm. Its rank correlation
with the K3 key is only 0.15–0.23 on the folds, so the two measures pick different rows.

| arm | assignment, q (low, mid, high) | cancer mean k 1 (folds +) | kolf2 | pooled rel cos | A/B/C energy |
|---|---|---|---|---|---|
| S_L30 | strength tercile, (3.0, 1.75, 1.75) | +0.0026 (3/3) | +0.0033 | +0.83 % [+0.15, +1.54] | ×0.95 |
| **S_P2** | strength tercile, (2.5, 2.0, 1.75) | **+0.0039 (3/3)** | **+0.0032** | **+1.76 % [+1.14, +2.43]** | ×0.93 |
| P2 (preregistered) | key tercile, (2.5, 2.0, 1.75) | +0.0034 (3/3) | +0.0011 | +0.95 % | ×0.99 |
| Sinv_L30 (q 3.0 on the STRONGEST third) | — | cos −10 % on every fold | −0.021 cos | — | — |

- Weak-response rows benefit from harder sharpening as the unreliable ones do. The strength-based assignment has the
  better kolf2 and direction. It moves more A/B/C energy (×0.93; the re-level refunds it inside the price).
- Sharpening the strongest rows harder is destructive, which confirms the mechanism: extra sharpening denoises weak or
  unreliable rows and distorts strong ones.
- This is one more exploratory fork after 492b/c. The preregistered build stays P2 (key).
  - A strength-based or combined assignment is a later variant, to be preregistered on its own if P2/L30 read positive.

### Knob definition, strength assignment (`stercile`; answer to the coordinator, 23:40Z; what 492d/493 computed)
- **strength_t = the L2 norm of target t's UNKEYED stack row, over ALL non-panel columns.** It is not gated, and it is
  computed once on the shared matrix, not per context.
  - On A/B/C: `np.linalg.norm(dump["carriers"][:, nonpanel], axis=1)` on `dump_h1x2_cd4x15.npz` (sha `bb19f0e5…`), whose
    `carriers` are the stack composition before `apply_pert_scale`. The dump's `bulk` is `apply_pert_scale(carriers, K3 w)`
    exactly (492c's known answer), so strength_t = ‖bulk_t‖ / (w_t · c) equally. c is the one global constant, which
    does not change ranks.
  - On PLOLO (`arms492.py --d`, `comp493.py`): the same norm of the unkeyed composition rows (`fST`, `FC.final("h2c1")`),
    over the fold's non-panel columns, whose panel coordinates are zero.
- Cut points as for the key: `np.quantile(strength, [1/3, 2/3])` + `np.digitize`. Low (weakest third) gets q 2.5, mid
  2.0, high 1.75. Strength is continuous with no ties (A/B/C gives 100/100/100), so no floor guard is needed.
- Rows keep their own gated non-panel norm after sharpening, and panel columns are untouched, as for `tercile:`.

## Board-calibrated predictions (`calib492.py`, 23:50Z): the instrument fitted on the board's own levelled q steps
Calibration rows are the same stack, key and DE channel, re-levelled, from 420's PLOLO rows and A/B/C energies:
| board step | board Δpds / Δraw mse | PLOLO cancer Δpds (literal / board-emul.) | k_pds (lit / bd) | k_X (the levelled cliff price's X factor that reproduces the board's Δmse) |
|---|---|---|---|---|
| 433h → 433k (q 1.25 → 1.5) | +0.00462 / −0.01790 | +0.01335 / +0.01273 | 0.346 / 0.363 | 0.12 (the 446b-anchored energy route fits this step poorly) |
| 433k → 446b (q 1.5 → 1.75) | +0.00317 / −0.00764 | +0.00950 / +0.00882 | 0.334 / 0.359 | **0.61** (≈ 417a's k_sh 0.65) |

- **pds transfer is stable at k_pds ≈ 0.36 across both board steps.** kolf2 reads 0.75–0.84× of the board's Δpds.
- For mse, the step adjacent to 446b gives k_X 0.61. 0.12 is kept as a pessimistic bound.

**Predicted board deltas vs 446b** (k_pds 0.361; k_X 0.61 / pessimistic 0.12; before the ballast's nmae cost, −0.0002 for
P2 / S_P2 per `k492.py`):
| arm | Δpds raw | Δraw mse | Δoverall (k_X 0.61 / 0.12) | net of ballast |
|---|---:|---:|---|---|
| **P2 / P2g** (key tercile) | +0.0027 | −0.0022 | **+0.0014 / +0.0011** | **+0.0012 / +0.0009** |
| **S_P2** (strength tercile) | +0.0026 | −0.0042 | **+0.0017 / +0.0011** | **+0.0015 / +0.0009** |
| L30 | +0.0020 | −0.0009 | +0.0009 / +0.0008 | +0.0009 / +0.0008 |
| F5P2 (floor 0.5 + P2) | +0.0039 | −0.0037 | +0.0021 / +0.0016 | +0.0019 / +0.0014 |
| G2 (global q 2) | +0.0022 | −0.0023 | +0.0012 / +0.0009 | +0.0010 / +0.0007 |
| CTRL (P2's q permuted) | +0.0027 | +0.0017 | +0.0007 / +0.0009 | — |

For bins: P2g's Δpds of ~+0.0027 raw is readable, because bulk members reproduce on the board to 16 significant figures.
The predicted Δoverall is about +0.001. The targeting's share, P2 minus its permuted control, is about +0.0007.

## Addendum 492e (00:00Z, EXPLORATORY): robust to the key — the same arms under the LODO-V1 key (`--e`)
Here the bulk is keyed by LODO-V1 (415's key on PP, clip 0.3–3; the October switch if 446bj2 reads `key_ok`), and the
tercile comes from V1's w with the floor guard. Known answers: V1 weights == 418's saved `V1__PP`, and
`stack|V1+sh1.75` == comp420, both |d| 0 on 4/4 folds. V1's floor holds 22–25 % of targets.

| fold | P2g (V1 tercile) rel cos / Δpds_bd | S_P2 | CTRL (P2g permuted) |
|---|---|---|---|
| k562_gwps | +1.48 % / +0.0060 | +1.08 % / +0.0043 | +0.33 % / +0.0055 |
| orion_hct116 | +1.16 % / +0.0090 | +1.76 % / +0.0076 | −0.68 % / +0.0068 |
| orion_hek293t | +0.73 % / +0.0064 | +1.35 % / +0.0101 | −1.36 % / +0.0075 |
| kolf2 | +0.65 % / +0.0031 | +1.75 % / +0.0047 | −3.17 % / −0.0023 |

The per-target sharpening gain holds with either key, at the same size, and the targeting again shows in direction.

## Addendum 492f (00:05Z, EXPLORATORY): the strongest-response third's dose under the strength assignment
`arms492.py <fold> --f`. Direction members vs 446b's bulk, realised:
| arm | q (weak, mid, strong) | rel cos cancer (folds +) | kolf2 | Δpds_bd cancer / kolf2 |
|---|---|---|---|---|
| S_P2 | 2.5, 2.0, **1.75** | **+1.75 % (3/3)** | +2.03 % | +0.0072 / +0.0049 |
| S_P2h15 | 2.5, 2.0, 1.5 | +0.83 % (3/3) | +3.32 % | +0.0045 / +0.0059 |
| S_P2h125 | 2.5, 2.0, 1.25 | −1.65 % (0/3) | +2.98 % | +0.0007 / +0.0062 |
| S_P3h15 | 3.0, 2.25, 1.5 | +0.56 % (2/3) | +4.13 % | +0.0064 / +0.0090 |

Keeping the strong third at 1.75 is right on the cancer folds. kolf2 prefers less sharpening there, as it did for the
global q curve (414: kolf2 transfers sharpening at ~1/3). S_P2 stays as defined.

## Mechanism note (00:25Z): why terciles, and why a continuous rule on these features won't do better
From 420's saved per-target arrays (`stack|K3` at q 1 … 2, realised), each target's own-cosine-maximising q:
- 66–68 % of targets prefer the grid's top (q 2) on the cancer folds, and 15–25 % prefer little or no sharpening
  (q ≤ 1.25). kolf2: 46 % and 38–45 %.
- **The per-target optimum is not predictable from the key or from strength:** Spearman(q_opt, w) +0.00 … +0.17 and
  Spearman(q_opt, strength) −0.09 … +0.01. The rows that sharpening hurts sit evenly across strength terciles.
- What IS predictable is the average gain. From q 1.75 to 2, the weak-strength tercile's own-cosine rises +0.0042 /
  +0.0037 / +0.0027 (k562 / hct116 / hek293t), against +0.0030 / +0.0017 / +0.0007 for the strong tercile.
- So the tercile rules are average-case denoising. A finer continuous q(w) or q(strength) should not beat them. The
  next lever would be a predictor of the ~20 % of rows that sharpening damages, and the key and strength are not it.

## Addendum 492g (04:40Z): re-priced on 446bj2 (the LODO-key stack), board-calibrated, for 479's preregistration

**Answer.**
- **P2g and S_P2 still price above +0.001 on the 446bj2 base: P2g +0.0012 and S_P2 +0.0013 overall (+0.0008 each at
  the pessimistic mse transfer).** The figures are net of the ballast's nmae cost.
- The LODO key's 0.5 floor did not absorb the per-target sharpening. P2g's PLOLO price on the LODO stack (+0.0035 at
  k 1, 3/3 cancer folds) equals its price on the K3 stack (+0.0034).
- The targeting is still real: P2g's q dealt to random targets loses direction (relative cos −0.5 %, kolf2 −2.2 %) and
  nets +0.0000.

**What changed vs 492.** The base bulk is `dump_h1x2_cd4x15_lodo05.npz` (sha `00ccfcba…`) with 446bj2's factors
(A 1.0962, B 1.1461, C 1.1681). The tercile key is the LODO artifact's clipped w (`lodo415_lo0.5_hi3`, `8c82bef3`), with
the floor guard at 0.5. The anchor is 446bj2's board raw mse, 0.6264081539. Scripts: `abc492g.py` (A/B/C energies, ballast
K, level maps), `arms492g.py` (PLOLO, lib420), `ballast492g.py` (the ballast nmae model, below), `analyse492g.py` (pricer),
`run492g.sh` (gated runner). Peak RSS 5.4 GiB per fold, run after 04:30Z at 37 GiB free+inactive.

### Known answers (all pass)
| row | measured |
|---|---|
| LODO dump `carriers` == K3 dump `carriers` (the same unkeyed stack) | identical |
| builder's apply_pert_scale(carriers, LODO w) vs the LODO dump's `bulk` | \|d\| 0 |
| LODO bulk at 446b's factors (= 446bj) vs ntclevel-446bj.json's measured A_own / B_own | −1.4 / −1.5 / −2.4 % (446b's documented route offset; only ratios are used) |
| emulated ballast K at 446bj2's factors vs 446bj2's meta | median 7, max 16 (meta 7 / 16) |
| per-target path at q 1.75 for every row vs comp420 `stack\|V1_5_3+sh1.75` (fold weights = max(418's `V1__PP`, 0.5)) | \|d\| 0, 4/4 folds |
| per-target path at q 2 for every row vs comp420 `stack\|V1_5_3+sh2` | \|d\| 0, 4/4 folds |
| cliff model at 446bj2's solved X reproduces its board raw | < 1e-12 |

### The pricer retrodicts the board's own key step (446b → 446bj2), no fit
| member | instrument | board |
|---|---|---|
| Δmse raw | k_X 0.61: −0.0068; k 1: −0.0111; k 1.407: −0.0156 | **−0.0145 (reproduced at k_X 1.31)** |
| Δpds raw | board-emulated × 0.361: +0.0015; **literal × 0.34: +0.00024** | **+0.00021** |

- **Key changes transfer at the composition class** (k_X ≈ 1.3, like 417e's stack), not at sharpening's 0.61.
- For key and floor changes the β-emulated pds over-predicts by 7×, and the literal pds is right. 492c's floor-stacking
  value was partly that emulation (F5 +0.0006 → +0.0002, F5P2 +0.0021 → +0.0016 with literal pds).
- The per-target arms are a sharpening-class change, so they stay at k_X 0.61 (0.12 pessimistic). Their literal and
  emulated pds agree (P2g +0.0026 either way), so the pds calibration choice does not move them.

### Board-calibrated predictions vs 446bj2 (pds = 0.34 × literal PLOLO Δpds; mse at k_X 0.61 / 0.12, re-levelled; nmae = ballast)
| arm | Δpds raw | Δmse raw (0.61 / 0.12) | Δnmae raw (ballast: D / K model) | **Δoverall net (0.61 / 0.12)** | PLOLO k 1 cancer (folds +) [pooled CI] | kolf2 k 1 | rel cos cancer [pooled CI] | A/B/C energy |
|---|---:|---:|---:|---|---|---:|---|---|
| **P2g** (LODO tercile, 2.5 / 2.0 / 1.75) | +0.0026 | −0.0031 / −0.0006 | +0.0010 / +0.0010 | **+0.0012 / +0.0008** | +0.0035 (3/3) [+0.0023, +0.0049] | +0.0020 | +1.26 % [+0.94, +1.66] | ×0.989 |
| **S_P2** (strength tercile) | +0.0025 | −0.0035 / −0.0006 | +0.0009 / +0.0006 | **+0.0013 / +0.0008** | +0.0038 (3/3) [+0.0026, +0.0051] | +0.0032 | +1.43 % [+0.89, +2.04] | ×0.94 |
| L30g (low tercile only, q 3) | +0.0020 | −0.0017 / −0.0003 | +0.0009 / +0.0009 | +0.0008 / +0.0005 | +0.0021 (2/3) | +0.0027 | +0.70 % | ×0.998 |
| P2g_K3key (the old K3 tercile on the LODO bulk) | +0.0026 | −0.0023 / −0.0005 | +0.0010 / +0.0009 | +0.0011 / +0.0008 | +0.0035 (3/3) | +0.0018 | +0.94 % | ×0.985 |
| G2 (q 2 for every row) | +0.0023 | −0.0019 / −0.0003 | +0.0009 / +0.0006 | +0.0009 / +0.0006 | +0.0029 (3/3) | −0.0003 | +0.84 % | ×0.90 |
| CTRL (P2g's q permuted) | +0.0023 | **+0.0027** / +0.0017 | +0.0012 / +0.0008 | +0.0000 / +0.0002 | +0.0018 | −0.0021 | −0.53 % [−2.1, +1.0] | ×0.875 |

- Per fold, P2g's relative cos is +1.72 / +1.28 / +0.79 % (k562 / hct116 / hek293t) and kolf2 +0.68 %. S_P2's is
  +1.15 / +1.80 / +1.34 % and kolf2 +1.72 %. Every cancer fold's target-bootstrap interval of the k 1 price is above
  zero for both arms.
- Where the ballast cost goes: P2g raises the low tercile's emulated K from 2.5 to 4.9 on this base (2.4 → 3.0 on K3).
  The 0.5 floor gives those rows 1.67× the amplitude, so their sharpened bulk needs more ballast. That doubles P2g's
  ballast cost vs the K3 base (−0.0003 vs −0.00015 overall). It is inside the net figures above.
- The LODO floor holds 29–33 % of each fold's targets and 97/300 (32.3 %) on A/B/C, all under a third. The floor guard
  is therefore a no-op here (P2 == P2g, bit-identical rows). It matters on October-like panels, where the LODO floor
  holds 20–45 % (438's `keys`): 493g checks that next.
- The two arms' A/B/C predictions differ by +0.0001, inside the calibration's own spread, so the same-day pair reads
  the assignment variable (key vs response strength) rather than a pre-decided winner. S_P2's larger energy change
  (×0.94) is refunded by the mandatory re-level (inside the price).

### Suggested bins for 479 (per arm, vs 446bj2; bulk members reproduce exactly, nmae noise ≈ ±0.0002 raw)
| reading | rule |
|---|---|
| **transfers** | Δoverall ≥ +0.0008 and Δpds raw ≥ +0.0015 |
| partial | Δoverall +0.0003 … +0.0008, or Δpds raw +0.0005 … +0.0015 |
| fails | Δoverall < +0.0003 or Δpds raw < +0.0005 |
| ballast check | Δnmae raw predicted +0.0006 … +0.0010. Above +0.0015 means the model under-prices the low tercile's ballast. fid / reach / jac should stay pinned inside ordering noise. |

### The ballast's nmae cost, tested on the board's bulk-only steps (`ballast492g.py`)
k492 had priced the ballast from one board step. Five board steps changed only the bulk, so two models were tested on
them. The K model prices the emulated mean ballast count. The D model prices the per-cell arithmetic mean's
displacement away from the carriers' design lfc (top-100 design genes per block, the scorer's epsilon, the −5
zero-repair floor).

| step | board Δnmae raw | ΔK mean | D (displacement from design) |
|---|---:|---:|---:|
| 410b → 417a (pre-MDE) | +0.000094 | +1.15 | +0.00085 |
| 417a → 417e (pre-MDE) | −0.000147 | −0.10 | +0.00013 |
| 433h → 433k | +0.000814 | +1.52 | +0.00115 |
| 433k → 446b | +0.001346 | +1.37 | +0.00129 |
| 446b → 446bj2 | +0.000060 | +0.18 | +0.00000 |

- Fits on the three steps of the current DE design (MDE era), slope through the origin: D 0.894 (residuals ±0.0002),
  K 0.00073 per unit K (±0.0003). Both put P2g at +0.0010 raw on this base. S_P2 gets +0.0009 (D) and +0.0006 (K).
- Pre-MDE the models do not apply: the D proxy assumes the −5 zero repair, and 410b → 417a cost nothing.
- **Mechanism caveat:** 490 measured, on H1 truth, that removing the ballast's lfc footprint altogether moves nmae only
  −0.0001. So the step costs more likely come through the K-class emitter designing on 400 − K carriers than through
  the footprint. D is then a proxy that happens to track it. For P2g the two models agree, so the price does not depend
  on which is right.
- A carrier pre-distortion that cancels the footprint is therefore not worth building (490's −0.0001).
