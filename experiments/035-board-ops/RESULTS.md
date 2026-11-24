# 035 — Board operating points, 2026-08-22 ~12:00 TRT

**Source:** live leaderboard, 113 teams, our v6 at rank 3 / 0.1369. Public model notes treated
as probably-true descriptions of *that team's* file, not as advice.

Each cell is **raw then scaled**. Overall is the mean of the six scaled numbers.

| rank | team | overall | pds | mse | jac | nmae | fid | reach |
|---|---|---|---|---|---|---|---|---|
| 1 | Ibrahim / CellSim | 0.1427 | 0.809 / **0.684** | **0.919 / 0.073** | 0.032 / 0.003 | 0.953 / 0.080 | 0.499 / −0.045 | 0.133 / 0.062 |
| 2 | Vivai m6 | 0.1380 | 0.777 / 0.611 | 3.909 / 0.000 | 0.033 / 0.008 | **0.921 / 0.131** | 0.502 / −0.036 | 0.180 / 0.115 |
| **3** | **us v6** | **0.1369** | **0.816 / 0.701** | 2.102 / 0.000 | 0.031 / 0.002 | 0.961 / 0.065 | 0.502 / −0.037 | 0.158 / 0.090 |
| 4 | IsoCell m4 | 0.1159 | 0.782 / 0.626 | 1.258 / 0.000 | 0.029 / −0.003 | 0.958 / 0.071 | 0.502 / −0.036 | 0.112 / 0.038 |
| 12 | larry ml-1 | 0.0850 | 0.745 / 0.544 | 1.003 / 0.006 | 0.035 / 0.012 | 1.009 / −0.017 | 0.498 / −0.048 | 0.090 / 0.012 |
| 16 | Hongjian effect-transfer | 0.0754 | 0.721 / 0.489 | 16.23 / 0.000 | 0.027 / −0.010 | 1.094 / −0.159 | 0.507 / −0.018 | **0.212 / 0.150** |

Best-on-board vs us (scaled):

| | us | best seen | gap (×1/6 of overall) |
|---|---|---|---|
| pds | **0.701** | 0.701 (us) | — |
| mse | 0.000 | 0.073 (Ibrahim) | +0.012 |
| nmae | 0.065 | 0.131 (Vivai) | +0.011 |
| reach | 0.090 | 0.150 (Hongjian) | +0.010 |
| fid | −0.037 | ~−0.014 (SimpleCellNYC; Claude cited −0.005 off this paste) | +0.004–0.005 |
| jac | 0.002 | 0.015 (Zafrens) | +0.002 |

Sum of best-in-class ≈ **0.18**. We are 0.006 behind #1, 0.001 behind #2.

## What the raw numbers encode

Raw **mse** is a magnitude thermometer (higher = bigger predicted effects, until you shrink
past the baseline ~0.99):

- larry / Vu Mai ~1.00 — almost control-copy on the profile; `reach` collapses
- Ibrahim 0.919 — **shrunk**, and the only top-3 team with positive mse. Cosine good enough
  that shrinking actually helps. Inverting `1−cos² ≈ 0.92` implies cos ~0.28 vs our ~0.13
- us 2.10 — moderate overshoot, mse clamped at 0
- Vivai 3.91 — larger α; pays mse=0, collects best nmae and strong reach
- Hongjian 16.2 — α≈1.6 class; best reach, nmae destroyed

So the field is on the **same α frontier** we measured: bigger effects help reach/nmae-until-
overshoot, hurt mse. Vivai is the one team off the dumb end of that curve (high α *and* best
nmae) — they have a better *direction* on the genes nmae reads, not just more magnitude.

## Public notes (treat as true for *them*)

- **CSBL csbl-v5:** emission-procedure variant on the same transplanted deltas. Same game as
  our v8 — don't copy the note, copy the fact that emission is the live iteration surface.
- **Sidechain SER-1n:** unshrunk log2FCs, min-variance cells. Explicitly probing whether
  shrinkage kills small PDS-useful signals. Our even emitter is already that pole.
- **Hongjian:** Replogle K562-GWPS, 272/300 direct + centroid fill, effect_scale=1.6.
- **Life Code v2:** K562 delta, hybrid log1p/fold reshape, global scale 1.5, STRING k=5 for
  the 26 missing perts. Jurkat harness for scale.
- **Bio_NUDT:** α=2 probe. Confirms people are still sweeping magnitude live.
- **Giovanni Discord control-copy:** overall −0.304, pds 0.500, raw mse **2.02**, nmae 0.999,
  fid 0, reach 0.063, jac 0. mse≠1 because ρ=0 when every prediction is identical.

## Decision implication

Do **not** spend the remaining slot chasing Ibrahim's mse or Hongjian's reach with a bigger α.
v6 already owns pds. The only one-variable change with a construction proof that pds/nmae stay
put is **exact-total spread** (v8). Wait for the somatic confirmation still running
(`out_som_spread`) before submitting it.
