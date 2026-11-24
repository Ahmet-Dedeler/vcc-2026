# 624 — A learned per-target bulk MAGNITUDE (Testh-style ridge) adds +1.3 % over 522's strength rider on October-like folds: small, unreadable on the board

2026-10-02, council session [7772dc]; council bulk item (assigned by [4d0207]). Analysis only: no build, no `vcc`, no
`cell-eval2`. `screen624.py` (stage 1, first-order screen from 522's saved per-target arrays, < 1 s) and `exact624.py`
(stage 2, fold522's exact scorer: q 1.75, render, ballast KMAX 16; 522's known answers re-run on load). Peak RSS
2.60 GB under gate602. Outputs `work/screen624.json`, `work/exact624.json`, `work/exact624.log`.

## Construction
- **Truth-free per-target features** (from 522's dumps; all exist on the day):
  - log key/median; log LODO terms (min, max, |diff|); cross-donor cosine; log donor norms (min, max); log Fnorm;
  - pnorm and cosPd (min, max); 518/522's strength z_src; carried.
- **Target:** the ideal relative weight along our direction on the percentile-cap arm, w* = w0 · (x/a)/(X/A) (522's
  u·w), in logs and clipped to [0.1, 20].
- **Fit:** ridge (α 100, chosen inside the training panels) leave-one-PANEL-out over 438's October-like panels: train on
  two panels' six folds, predict the third panel's three folds.
- **Prediction quality:** ρ(prediction, ideal) is 0.40–0.60, above the key alone (0.34–0.58) on 9/9 folds.
- **Arms:** a full replacement of the key is erratic (4/9 folds positive, first order). The **tilt** form is
  w = w_pct · exp(0.25 z_src above the floor) [522's rider] · exp(0.1 · z_res), where z_res is the ridge prediction
  residualised on log key, standardised and clipped at ±2.

## Result: EXACT, vs the percentile cap, relative aggregate cosine

| fold | 522 rider S | S + ridge tilt | increment over S | stage-1 approximation (S / S+tilt) |
|---|---:|---:|---:|---|
| efp k562 | +1.89 % | +4.18 % | +2.25 % | +2.14 / +4.63 |
| efp hct116 | +0.06 % | −0.63 % | **−0.69 %** | −0.07 / −0.80 |
| efp hek293t | +1.40 % | +2.09 % | +0.68 % | +1.58 / +1.99 |
| efr k562 | +3.97 % | +9.37 % | +5.19 % | +4.92 / +11.39 |
| efr hct116 | +1.15 % | +2.57 % | +1.40 % | +1.41 / +3.22 |
| efr hek293t | +0.32 % | +0.96 % | +0.63 % | +0.70 / +1.65 |
| real k562 | +0.95 % | +2.60 % | +1.64 % | +1.13 / +2.81 |
| real hct116 | +0.70 % | +0.94 % | +0.24 % | +0.84 / +0.99 |
| real hek293t | +0.46 % | +0.76 % | +0.30 % | +0.45 / +0.63 |
| **mean (folds > 0)** | **+1.21 % (9/9)** | **+2.54 % (8/9)** | **+1.29 % (8/9)** | |

- pds moves by ≤ 0.0006 everywhere.
- The ridge tilt's gain is concentrated: the top 5 targets carry 55–86 % of each fold's numerator change. The k562
  folds dominate.
- **The stage-1 first-order screen is validated.** It reproduces the exact ranking with ~10–20 % relative error.
  Per-target weight ideas can now be screened in seconds from 522's saved arrays before any exact run.

## Reading
- A Testh-style learned magnitude does carry signal beyond our key. As a tilt on 522's rider it adds +1.3 % aggregate
  cosine on 8/9 October-like folds. By 522's price, that is ≈ +0.0003 … +0.0012 overall on D/E/F (k 1–1.31, retention
  0.3–1.0). It is not readable on the A/B/C board, because the A/B/C panel's strength mix makes the tilt ≈ 0 there (as for
  522's rider).
- **Not adopted alone.** Its natural home is bundled with 522's optional strength rider, which needs a def_day
  pass-through and a rehearsal first. Before any adoption it still needs 522's own bar: ≥ cap on 8/9 Arc-like re-draws
  (`subpanel522.py`) and A/B/C PLOLO folds ≥ 0. That is about 30 min of work if the coordinator wants the rider bundle.
- Class: transfer (truth-free features, fitted on public October-like panels, leave-one-panel-out).

## 624c — 522's adoption bar: the ARC-LIKE re-draws (`arc624.py`, `work/arc624.json`, `work/arc624.log`)
Essential-like targets are removed and each fold is re-drawn to A/B/C's strength mix (20 draws per fold; 522's
`subpanel522` logic and seeds). Exact fold522 scores with the mask. Peak < 3 GB.

| fold | 522 rider S vs pct (frac draws > 0) | S + ridge tilt | increment of the tilt over S |
|---|---|---|---|
| efp k562 | +1.49 % (1.00) | +3.52 % | +2.03 % (1.00) |
| efp hct116 | +0.78 % (0.90) | −0.27 % | **−1.05 % (0.30)** |
| efp hek293t | +1.94 % (1.00) | +1.92 % | −0.01 % (0.60) |
| efr k562 | +2.95 % (1.00) | +6.66 % | +3.71 % (1.00) |
| efr hct116 | +0.14 % (0.65) | −0.33 % | −0.47 % (0.25) |
| efr hek293t | +0.38 % (0.80) | +1.01 % | +0.64 % (0.95) |
| real k562 | +1.74 % (1.00) | +5.65 % | +3.91 % (1.00) |
| real hct116 | +0.26 % (0.90) | −0.97 % | **−1.22 % (0.00)** |
| real hek293t | −0.04 % (0.30) | −1.20 % | **−1.15 % (0.00)** |
| **mean (folds ≥ 0)** | **+1.07 % (8/9)** | +1.78 % (5/9) | +0.71 % (4/9) |

**624 FAILS 522's bar** (≥ the cap on 8/9 Arc-like re-draws). The tilt's gain lives on the K562 folds and reverses on the
Orion folds once the panels look like Arc's. **Closed: do not adopt the ridge tilt, bundled or alone.** 522's strength
rider by itself still passes this bar (8/9), as 522 reported. The learned magnitude's held-out signal (ρ 0.40–0.60) is
real, but it does not survive the Arc-like panel mix as a weight.
