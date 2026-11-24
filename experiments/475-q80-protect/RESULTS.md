# 475 — 484's q80 fix on H1: FAILS its rule; the emitter stays as scored (04:47Z 09-27)

**Result.** Keeping protected (q80) down genes out of the mass-balance shrink makes nmae slightly WORSE on every seed, as
490's exact reconstruction predicted. The small fid gain from the |z| realisation, which 490 held fixed, does not rescue
it. The emitter stays as scored; its comment ("the strongest keep their delivered LFC") is wrong for down genes.
Correcting the comment would change the pinned emitter sha, so it is recorded here and in FINDINGS instead.

## Numbers (fix − 435's production emitter at the same seed, raw; target bootstrap, 3,000 draws)
| seed | subset | fid | reach | nmae | jac | overall-equiv | 95 % CI | P(≤ 0) |
|---|---|---:|---:|---:|---:|---:|---|---:|
| 0 | all 85 | +0.00055 | +0.00199 | +0.00037 | 0.0 | +0.00057 | [−0.00006, +0.00147] | 0.048 |
| 1 | all 85 | +0.00043 | −0.00031 | +0.00039 | +0.00005 | +0.00009 | [−0.00106, +0.00118] | 0.43 |
| 2 | all 85 | +0.00044 | +0.00145 | +0.00038 | −0.00001 | +0.00040 | [−0.00033, +0.00131] | 0.18 |
| 0/1/2 | non-essential 67 | +0.0002 | +0.0026 / −0.0004 / +0.0018 | +0.0001 | ~0 | +0.0006 / −0.0001 / +0.0004 | | 0.10 / 0.53 / 0.25 |

pds and mse are exactly 0 (totals held). Known answer: `orig` == 435 seed 0, stats bit-equal (`work/ka475.txt`).

## Against the preregistration
- Original prediction nmae −0.002 … −0.008: **wrong**. The addendum written before any number (from 490) predicted
  nmae +0.0000 … +0.0010: **right** (+0.0004).
- Rule: P(≤ 0) < 0.10 on every seed (fails: 0.43 on seed 1), and nmae improving on the non-essential subset (fails:
  +0.0001). **Not a board arm.**
- Mechanism (490): kf and the MDE add are substitutes. The MDE was sized on the shrunk magnitudes, so un-shrinking the
  protected down genes overshoots. The fid +0.0005 is the realisation term: a few more correct-sign calls among the
  un-shrunk genes. It is small and seed-dependent in reach.

Note: the first `fix` run failed at import. The copy had a docstring before `from __future__`, fixed as a comment.
`chain475b.sh` then ran the fix seeds.
