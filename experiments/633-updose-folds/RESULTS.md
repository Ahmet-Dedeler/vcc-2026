# 633 — held-out UP-dose reading for 630: clean cancer line K562 wants ~2–3× H1's dose (oracle shape)

> **VOID ROWS (2026-10-02T19:10Z; [21dd4c]'s 616 audit, 7dac7a98, `experiments/616-clean-orion-folds/AUDIT-orion-dependence.md`):** the thinned Orion gate tables' DE truth was ~98 % selection artefact (HEK293T: 37,748 calls on the thinned table vs 774 on clean raw cells; 13 / 132 targets scoreable). In this experiment, the Orion rows (own and d20) are VOID. What remains: the K562 binding row and the h1 retrodiction; 633j Jurkat is unaffected (Nadig, not thinned).

14:28Z, council session [d4cd5f]. `dose633.py` (committed 3232efa2, amended 8dafb2ac before the first run: Orion rows
descriptive only), `work/dose633.out`, `work/dose633.json`. Class: **scorer mechanics**. The pass carries the dose to
D/E/F, and this asks which dose transfers, using a held-out line rather than the board (11.1c.2).

## Method
612's call-set-aware fold model (ds612 / ds612b / ds612c, read-only):
- base = the plain-hook form (BETA·cons0 + s·add);
- called = production's head 400 + CPM block 2,200;
- the UP pass at dose K raises uncalled UP-delivered genes to max(base, K·0.24/√cnt) (630's target) at FULL delivery
  (an oracle, so only the dose SHAPE transfers);
- nmae at 521's pools f 1 / 0.5 / 0.35. "d20" rescales the fold's control depth to A/B/C's 20k library.

## Retrodiction gate (fixed before the run): PASS
The h1 row against 630's real-emitter oracle (raw Δnmae, 3 seeds):

| | K 0.5 | K 0.75 | K 1.0 | K 1.5 |
|---|---|---|---|---|
| f 1 model / real | −0.0107 / −0.0058 | −0.0187 / −0.0109 | −0.0200 / −0.0107 | **+0.0235 / +0.0367** |
| f 0.35 model / real | −0.0079 / −0.0071 | −0.0156 / −0.0142 | −0.0244 / −0.0220 | **−0.0340 / −0.0280** |

Signs agree in 8/8 cells; the best K at f 0.35 is 1.5 in both; K 1.5 overshoots at f 1 in both. The model runs ~2×
the real magnitude at f 1 and ~1.1–1.2× at f 0.35, so its curve SHAPE is the reading.

## Result (Δnmae raw vs base, − = better; argmin K per regime)

| surface | f 1 | f 0.5 | f 0.35 |
|---|---|---|---|
| h1 (own depth) | K 1 | K 1.25 | K 1.5 |
| **k562_gwps d20 (binding)** | **K 2** (−0.060 vs K 1 −0.029) | **K 3** (−0.044 vs −0.016) | **K 3** (−0.039 vs −0.015) |
| k562_gwps own depth | K 1.5 | K 2 | K 2 |
| orion_hct116 d20 (descriptive) | K 2 | K 2 | K 2 |
| orion_hek293t d20 (descriptive) | K 1.5 | K 2 | K 2 |

**Decision (rule fixed before the run): the gate PASSES and K562 prefers a dose above 1.0 → the October dose question
is REOPENED.** A switch also needs a second clean cancer surface (the amendment's condition).

## Reading
- On a clean cancer line, the L1-optimal UP dose is ~2–3× H1's in MDE units at A/B/C depth. This is what the
  sign-accuracy argument predicts: m* = F⁻¹((2a−1)/2a) grows with a. K562's uncalled-UP sign accuracy is 0.73 / 0.89
  (613e) against H1's 0.54 / 0.58. K562 shows no overshoot up to K 2 at f 1, where H1 already overshoots at K 1.5.
- **H1 understates the October dose for cancer lines.** Today this is moot: the v3 pass is capacity-saturated above
  K 0.75 (630's diagcap). But it makes capacity more valuable on cancer lines than H1's +0.003 headroom suggests. At
  full delivery K562's curve keeps improving to K 2–3: about 2.5× K 1.0's Δ at f 0.35.
- The Orion rows agree (K 1.5–2) but are descriptive only (609 / 609b: depth-linked skews in both directions).

## What it changes
1. 630's single A/B/C arm stays b_k100_q20r3p10. The capacity-bound pass cannot deliver a higher dose today.
2. **Capacity (v5) moves up the queue.** Its value on cancer lines is larger than H1 prices it.
3. October dose: needs a second clean cancer surface. Candidates: Nadig Jurkat / HepG2 (CRISPRi, 3′), Replogle
   RPE1, or 7772dc's Broad-16 gate files. Then a preregistered K 1.0 / K 2.0 pair on the board, once capacity can
   deliver K 2.

## Try count
One run: 7 doses × 3 regimes × 4 surfaces (+ d20), with the decision rule and gate fixed before it.

## 633j — Jurkat (Nadig 2025) as the second cancer surface: NOT CLEAN by its preregistered check → no confirmation (14:57Z)
Gate: `gates633j.py` (300 random K562-covered targets, ≥ 50 cells; 12,013 controls; KA: reader bit-identical, scipy
call agreement 1.0, lfc exact). Reading: `dose633j.py` (checks and rule committed before the run), `work/dose633j.out`.

| artifact check (fixed before the run) | value | limit | |
|---|---|---|---|
| depth pert / ctrl | 0.927 | 0.85–1.15 | pass |
| CPM-block gate genes DOWN share | **0.703** | ≤ 0.62 | **FAIL** |
| low-expression uncalled UP share | 0.567 | ≤ 0.65 | pass |

**By the rule fixed before the run, Jurkat's rows are descriptive only, and K 1.0 stays the October dose** (no second
clean cancer surface).

Descriptively, at A/B/C depth (d20) the argmin K is 1.5 (f 1) / 2 (f 0.5) / 3 (f 0.35), the same direction as K562. At own
depth it is 1.25 / 1.5 / 2.
The failing check concerns highly expressed CPM-block genes. In an essential-gene screen that skew may be biology
(global down of growth programmes) rather than a depth artifact (depth ratio 0.93 is clean), and the UP pass acts on
uncalled genes, where check 3 passes. That argument came AFTER the result, so it does not change the decision. It is a
reason to build a third surface with the checks rewritten for the genes the lever touches, preregistered before any read:
- Replogle RPE1 (data/external/replogle2022);
- Nadig HepG2;
- Broad-16 lines (7772dc's gate files).
**Status: the October UP dose stays K 1.0. One clean cancer line (K562) says higher; one not-clean line (Jurkat) agrees
descriptively.** Two readings were tried for this switch (633 K562, 633j Jurkat).
