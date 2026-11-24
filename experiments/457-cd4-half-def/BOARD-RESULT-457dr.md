# 457dr BOARD RESULT — 0.2695096 (+0.0031484 vs 453dr): the CD4 half transfers on the October path; October keeps CD4 ×1

Entry `1FfirZsGsAZvaXYTd39J`, submitted 2026-09-29 00:01Z by `457/arm0929.sh` slot 1 ("Don't compare the outputs—compare
the inputs.", description empty), published 01:50Z (the server took 1 h 49 min), captured 02:41Z, rank 7. Package
`submissions/457dr-october-path-abc-cd4x15-q1.75.vcc` sha256 `c3afbf8e…`. Counterpart 453dr `q3VuTImh4FRGKXH2waCp`
0.2663613. Bins: `PREREGISTRATION.md` (2026-09-26 02:45Z). Reader: `read0929.py cd4` (written 09-28 14:10Z).

457dr = 453dr's October path (`def_day.sh` on 180's relabelled bundle) with the stack dump at CD4 ×1.5 instead of ×1; config
diff vs 453dr exactly the dump and the re-level factors (gate457).

| member | Δ raw vs 453dr | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | +0.001713 | +0.003550 | +0.000592 |
| mse | −0.014313 | +0.015054 | +0.002509 |
| nmae | +0.000200 | −0.000379 | −0.000063 |
| fid | −0.000034 | −0.000118 | −0.000020 |
| reach | +0.000711 | +0.000774 | +0.000129 |
| jac | +0.000004 | +0.000010 | +0.000002 |
| **overall** | | | **+0.0031484** |

Scaled members: pds 0.81556, mse 0.36445, nmae 0.17129, fid 0.00534, reach 0.23205, jac 0.02837.

## Against the bins
- **BIN ≥ +0.0010: the CD4 half transfers as 438 prices it** (predicted +0.0017 at k 1 … +0.0022 at k 1.3; band +0.0005 …
  +0.0035). It landed at the top of the band, 1.83× 438's k 1 price (+0.00172; 1.43× the k 1.3 point; corrected 2026-09-29, 533 D13), all of it bulk (pds +0.0006, mse +0.0025 overall).
- **DE pinned** (|Δ| ≤ 0.0002 overall; reach +0.0007 raw is inside the ordering channel's noise).
- **October decision: `cd4_stack` stays ×1.** The rule fixed on 09-27 re-opens it only at ≥ +0.0040 (438: CD4 ×1.5 loses on
  every October-like panel; an A/B/C gain cannot price D/E/F's CD4 coverage). The CD4 dose remains a final-phase feedback
  knob (docs/FINAL-PHASE-FEEDBACK-PLAN.md row 4).

## What it adds
With 453dr (−0.0036816 vs 446b, bulk −0.00309), the October path at production's CD4 dose sits at 457dr − 446b =
**−0.0005332** overall, bulk ≈ +0.0000 and DE −0.0005 (reach noise). So on A/B/C the October path's whole bulk cost was the
CD4 half; the NTC-proxy re-level, the day's operator and the re-fitted keys and banks cost nothing measurable. The path
reproduces production when its inputs match.
