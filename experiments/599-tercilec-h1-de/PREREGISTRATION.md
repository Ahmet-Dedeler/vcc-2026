# 599 — Does tercilec's ballast rise cost DE members on a thin-coverage panel? (preregistration)

Written 2026-10-02 ~17:35Z by the session "Explore contribution areas" ([36a81a], formerly [57577c]), at the coordinator's
request (652's F4: 598's DE side is unpriced and it can re-tier half of an October panel). It was written BEFORE any arm
run; only the known-answer run (mode ka) may precede it. The design numbers below come from `prep599.py`. That step
reads no truth: no h1-test perturbed cells, no scorer. Light: each run ~6 min at ~7 GB through gate602.

## Question
Under 598's `tercilec`, the non-carried targets of a thin-coverage panel move tier, and the ballast K moves with them.
On the board, the bulk-only steps that raised mean K cost nmae (492g: 0.00073 raw per unit of mean K, ±0.0003), and S_P2
failed on its low-tier ballast. What does tercilec cost the four DE members on the production emitter, against plain
`tercile`?

## Instrument
The H1 real emitter: run496's 446bj2z construction.
- Production emitter 0f56209d, CallSet(400, 2200, 1.9), CLEAN4 q80 block prior, coarse MDE (c 0.1, cap 0.1, min 0.3);
  the October zs rider (`--zs dn_k250_zs`) is scored as primary and pre-zs as secondary.
- It is driven run612-style (`run599.py` monkeypatches run496 and copies nothing). 611 validated this instrument for
  delivery levers; a tier rule moves K, a delivery lever.
- DE members are scored as in 612: lib521's exact member functions, regimes f 1 / 0.5 / 0.35 / cmp, 192's INV for the DE
  overall-equivalent. **f 0.35 is primary** (611's validated regime).
- pds and mse are the runs' own raw bulk members. Their overall-equivalent uses the BRIEF's mid anchors and is
  approximate.

## Design (fixed now; `common599.py`, `work/prep599.json`, `work/prep599.npz`)
- **Bulk.** 410's whitened shared-energy key on the 85 H1 targets: all 85 carried (5 by 2 donors, 80 by 3).
  - Coverage-thinned: a seeded subset (seed 599) of the carried targets loses its key, so **43 of 85 are non-carried**
    (w = 1), as on efp (~153–160 of ~300).
  - Weights are key410.weights at an October-like clip (lo 0.5, hi 6), re-normalised over the 42 remaining carried
    targets. 14 sit at the floor and 3 at the cap.
  - The scaling is energy-preserving over non-panel genes (c 0.5946), run496.ps_bulk's arithmetic.
- **Render.** The ballast's bulk profile is the builder's own
  `profile_from_effect(eff_bk[i], b_ctrl, min_cpm=10, sharpen_q=q_i, protect=panel)`, production's ballast call. Before
  replacing it, run599 asserts that the profile run496 passed equals the unsharpened one recomputed, bit for bit.
- **Arms** (the q per target comes from 598's patched `sharpen_q_by_tercile`, 5cf0e012's text):

  | arm | rule | LOW / MID / HIGH | non-carried tiers | edges |
  |---|---|---|---|---|
  | s175 | scalar q 1.75 (sharpen_pt off) | — | — | — |
  | **t** | `tercile:2.5,2.0,1.75` (DECIDED) | **21 / 0 / 64** | all 43 HIGH | [0.5946, 0.5946]: both cuts sit on the w = 1 block, MID is empty |
  | **tc** | `tercilec:2.5,2.0,1.75` (598) | **14 / 57 / 14** | all 43 MID | [0.3175, 0.8168] over the 42 carried |

  t → tc changes q for 57 targets:
  - the 43 non-carried targets and 7 carried targets go 1.75 → 2.0 (sharper);
  - 7 carried targets go 2.5 → 2.0;
  - 28 are unchanged.
- **Seeds:** 0, 1, 2 (emission). Same picks across arms (the pick stream is the builder's, seeded).

## The ballast K every run will use (exact: K is fixed before emission)
| | s0 | s1 | s2 |
|---|---:|---:|---:|
| mean K: s175 / t / tc | 6.365 / 7.000 / 7.941 | 6.329 / 6.976 / 7.929 | 6.329 / 6.976 / 7.965 |
| ΔK (tc − t), all 85 | +0.941 | +0.953 | +0.988 |
| ΔK on the 43 non-carried / the 42 carried | +2.02 / −0.17 | +2.02 / −0.14 | +2.09 / −0.14 |
| K at the cap 16 (share): t → tc | 0.28 → 0.34 | 0.28 → 0.35 | 0.28 → 0.35 |

33 targets change K: 28 up, 5 down.

## Prediction (tc − t, post-zs, 3-seed mean)
Two hypotheses are registered, and the read discriminates them as far as the noise allows.
- **H_K (492g's K model):**
  - Δnmae raw = 0.00073 × 0.96 ≈ **+0.0007**. Allowing for the slope's ±0.0003 and 611's board ≈ H1 × 1.2–1.6, the band
    is **+0.0003 … +0.0010** on H1.
  - DE overall-eq at f 0.35 ≈ **−0.0002** (band −0.0005 … 0.0000).
  - fid, jac and reach: small decreases or null (400 − K carriers means slightly less power); |Δ| < 0.002 raw each.
- **H_0 (416 / 490: on H1, K itself was not the effect):** |Δnmae raw| < 0.0002 and DE overall-eq within ±0.0002.
- **Bulk (pds, mse):** no sign is registered. The rule's bulk value is priced elsewhere (493h, the 438 instrument: efp
  +0.84 % rel cos vs +0.48 %). H1's bulk delta is a secondary single-panel read. There is no re-level here, so mse reads
  the un-levelled render.
- **Anchors:** t − s175 has ΔK +0.64 and tc − s175 has ΔK +1.60. Under H_K, Δnmae raw ≈ +0.0005 and +0.0012.

## Reading (fixed now; the decision itself is the coordinator's, at APPLY DAY v2 step 9)
On the primary read, DE overall-eq (tc − t) at f 0.35, post-zs, 3 seeds:
- **R1, clean:** ≥ −0.0003, and Δnmae raw ≤ +0.0010. The DE side does not block adoption; tercilec's net on a thin panel
  is its bulk gain less ≤ 0.0003.
- **R2, bounded cost:** −0.0010 < DE overall < −0.0003. A real cost, to be set against the bulk gain; the net is reported.
- **R3, blocking:** ≤ −0.0010, or Δnmae raw > +0.0015. The K rise costs more than the K model says. Shipping tercilec
  as-is for thin panels is not supported. The variant to price next is carried-only cut points with the non-carried
  block kept at HIGH (its K unchanged); it is named here, not run.

**Noise:** per-seed paired differences, reported as mean ± sd. If the f 0.35 mean is within 2 standard errors of 0, and
gate602 slots allow, seeds 3 and 4 are added for t and tc only and the read is made on 5 seeds. This extension is
registered now.

## Known answers
- **Before this file (prep599):**
  - KA-bulk: mode ka's bulk (cached key, clip 0.3 / 3) == run496.ps_bulk, bit for bit. **PASS**.
  - KA-K: mode ka's emulated K == 435's recorded ballast_K, all 85 targets, seeds 0–2. **PASS**.
- **Before the arms (chain599):** mode ka seed 0 must reproduce 435's base run (stats npz arrays and raw metrics,
  pre- and post-zs) bit for bit, or nothing runs.
- **In score599:** each arm's recorded K == the K table above.

## Files
- `common599.py`, `prep599.py`, `run599.py`, `ka599.py`, `score599.py`, `chain599.sh`.
- `work/prep599.json` and `work/prep599.npz` (the design and the exact K).
- `work/key599.npz` (the cached key).

## Amendment, 2026-10-02 17:35Z (after the KA run, BEFORE any arm run)
The KA run (mode ka, seed 0) reproduced every array of 435's stats npz bit for bit, pre- and post-zs, and the post-zs
raw metrics exactly. The pre-zs raw nmae differed in the last bit: 0.9427040517722616 vs 435's 0.9427040517722618.
- The scorer itself is not bit-deterministic on identical inputs. 435's own saved arrays, scored in three fresh
  processes, give 0.9427040517722617 (`hashseed599.py`, `work/hashseed599.out`). That is three different last bits for
  one input: the order of cell-eval2's DE rows, averaged by `np.nanmean`.
- The KA's purpose is the construction. The criterion is therefore amended:
  - the stats arrays must still match bit for bit;
  - the raw metrics must match to 1e-12.
- The read is unaffected: score599 recomputes every DE member from the stats arrays, and the arms are compared within one
  environment.
