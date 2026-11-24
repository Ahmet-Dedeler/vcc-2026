# 478 — CLEAN4 signs on uncalled genes ARE deliverable through the coarse MDE (426's over-fill closure predates it); analytic H1 prize nmae −0.0015 … −0.0025. Needs an emitter hook and an H1 run

2026-09-27 00:15Z, third session, at the coordinator's request.
- Analytic only, through 490's exact emitter reconstruction (`recon_h1.py`), on 435's three coarse-MDE H1 seeds, with
  calls / p / reach order held and cell-eval2's own DE functions (`score_de.py`).
- Script `price478.py`, output `work/price478.{json,log}`, 31 s, < 2 GB.
- Not preregistered: this is the coordinator's exploratory "any analytic route" request, reported as such.

## Why this was worth re-asking
- 426 closed prior signs through the A-channel ("class over-fill": ~2,200 flips per target over-filled a class by 8,324
  counts against 14,430 sub-gate counts).
- That was measured with the ORIGINAL split. Its single-count polish drains the middle class.
- 426h then shipped the **coarse** channel: non-block moves only between the extreme library classes, in whole counts.
  With it, an A-balanced move set is count-balanced per class, and the MDE's own up/down balance (fu, fd) is the class
  capacity. Prior signs were never re-priced on it.

## Known answers (all pass)
- The reconstruction's per-target diag known answers (85 targets × 3 seeds).
- The delivered tables reproduce the harness metrics (≤ 1e-9).
- The MDE re-implemented with an arbitrary sign vector reproduces recon's A_int exactly at the shipped signs
  (`allclose` 1e-9, every block).

## Result (Δ vs the shipped coarse MDE with design signs; H1, 85 targets; overall-equiv with 475's slopes; seed 0 CI)
| complement sign source | flips / target | fd (shipped 0.79) | nmae (3 seeds) | reach (3 seeds) | overall-equiv (3 seeds) [seed-0 CI] |
|---|---:|---|---|---|---|
| CLEAN4 where covered (τ 0) | 2,080 | 0.76 | −0.0024 / −0.0024 / −0.0026 | +0.0004 / −0.0011 / +0.0002 | +0.0008 / +0.0005 / +0.0008 [−0.0022, +0.0035] |
| CLEAN4 where \|P−0.5\| ≥ 0.05 | 810 | 0.77 | −0.0015 / −0.0014 / −0.0015 | **+0.0035 / +0.0035 / +0.0017** | **+0.0011 / +0.0011 / +0.0007** [−0.0001, +0.0027] |
| \|P−0.5\| ≥ 0.1 | 166 | 0.78 | −0.0001 | 0 | +0.00003 |
| payers (sub-gate up-mass, ρ 0.25 / 0.5) | — | 0.78–0.80 | +0.0001 more at most | — | — |

- **The balance does not collapse.** fd moves from 0.79 to 0.76 with 2,080 flips per target, because the coarse
  channel's balance already absorbs the prior's down-lean. Sub-gate "payer" mass is not needed: ρ 0.5 adds nmae −0.0001.
- **The prize is small but consistent.** All covered genes give nmae −0.0025 raw (the same size as 426's analytic
  −0.003 with the old delivery). The confident subset (τ 0.05) gives −0.0015 plus reach +0.002 … +0.0035 on all seeds.
  fid and jac do not move, because calls are held.
- **Board translation is uncertain.** This is a DE sign-prior change, on uncalled genes, at MDE magnitude. Per
  `docs/H1-BOARD-CALIBRATION.md` rules 3–4, H1 cannot price its nmae or reach reliably: 392's block prior worsened
  nmae on both H1 and the board.
  - Read literally with the MDE's own nmae transfer (1.76×) it would be −0.003 … −0.004 raw, about +0.001 overall.
    The reach part should get no credit before a board read.

## What it needs before any slot
1. An emitter hook: the coarse path's sign source for complement genes, e.g. an `--mde-prior-sign`-style option
   restricted to non-block eligible genes with a confidence threshold τ. It is a production edit, so it belongs to the
   coordinator.
2. A real H1 emitter run (the 435 harness, 3 seeds), to measure leaked calls (unsuppressible, clip) and the |z|
   realisation, which this analytic path holds fixed. 426's c 0.15 leak was a magnitude effect, and the magnitude is
   unchanged here. That is the main unknown.
3. If 458d reads positive on 09-30, re-run this with Arc's programme as the complement sign source (a one-line change
   of `prior_P`). The prize scales with (2a − 1), where a is the source's accuracy on uncalled gate genes.

## Addendum 478b (00:25Z): 491c's recipient-side programme as the complement sign source (`price478.py --recip`)
The sign source is z(CLEAN4) + z(proj_k2 of CLEAN4 on H1's control PCs) + z(−identity vs A/B/C), mapped to a pseudo-P
(0.5 − 0.5·tanh(z/2)) so that the same |P − 0.5| thresholds apply. It covers every kept gene, so it flips more.

| threshold | flips / target | fd | nmae (3 seeds) | reach (3 seeds) | overall-equiv |
|---|---:|---|---|---|---|
| τ 0 | 4,075 | 0.49 | −0.0016 / −0.0016 / −0.0020 | −0.0080 / −0.0087 / −0.0074 | −0.0010 / −0.0012 / −0.0008 |
| τ 0.1 | 3,195 | 0.50 | −0.0033 / −0.0034 / −0.0037 | −0.0019 / −0.0011 / −0.0021 | +0.0006 / +0.0007 / +0.0006 |
| τ 0.2 | 2,325 | 0.53 | −0.0047 / −0.0046 / −0.0050 | −0.0019 / −0.0011 / −0.0021 | +0.0010 / +0.0011 / +0.0010 |

- More nmae than CLEAN4, but a reach cost on every seed. The balance binds harder (fd ≈ 0.5, more down-signs), and wrong
  flips on weak targets' pool genes cost reach. That is the same mechanism as Arc's programme on small-gate targets.
- Net about +0.001, no better than CLEAN4 at |P − 0.5| ≥ 0.05 (+0.001 with a reach GAIN). The source choice stays
  CLEAN4-confident vs Arc's programme after 458d's read. Peak RSS 3.8 GB.
