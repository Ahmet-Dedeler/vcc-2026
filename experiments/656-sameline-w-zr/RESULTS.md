# 656 — `sameline_w = 0.5` STANDS under production's zero repair (3 seeds; every higher w loses by ≫ 2 sd)

Runs 2026-10-02 21:07–22:24Z (gate602; commit times), council session [7772dc], assigned by the coordinator after 654's
sweep. PREREGISTRATION 28d99181 was committed before any run. Written up 10-03 after the session gap. Class: transfer
(450's same-line rows). This re-prices a DECIDED October value with one construction change: `--zero-repair -5`.

## Result (450's board units: read450's H1→board slopes on the harness's raw members at f 1, vs 435's base; seeds 0–2)
| bank (flip p, gate signs) | w 0.5 | w 1 | w 2 | Δ(w 1 − 0.5) mean (per-seed sd) | Δ(w 2 − 0.5) |
|---|---:|---:|---:|---:|---:|
| p 0.2 (≈ 0.66) | **+0.0323** | +0.0316 | +0.0274 | −0.0007 (0.0002) | −0.0049 (0.0005) |
| p 0.3 (≈ 0.61) | **+0.0168** | +0.0088 | +0.0029 | −0.0080 (0.0003) | −0.0139 (0.0001) |

- **VERDICT (the preregistered rule): 0.5 STANDS.** No higher w beats 0.5 at either bank quality, and every Δ is negative
  on every seed, between 3 and 100 per-seed sd below 0.
- **Zero repair did what it should.** No arm has a scored zero-total cell (non-own, in any pool) on any seed.
  - Seed 0 vs 450 §10's unrepaired numbers: w 2 rises +0.0013 (p 0.2) and +0.0024 (p 0.3); w 0.5 rises +0.0002 at both.
  - So the artifact did penalise the higher doses, as 654's sweep predicted. It was 4–6× smaller than the w 2 gaps it might
    have explained: 0.0013 vs 0.0049 at p 0.2, 0.0024 vs 0.0139 at p 0.3.
- **Prediction check:**
  - "0.5 stands" holds.
  - "w 1 ties 0.5 at p 0.2" is slightly off: w 1 is consistently below 0.5 (−0.0007, about 4 sd), not tied.
- **KA:** all 6 seed-0 arms equal 450's stored unrepaired seed-0 arms (q, p, lfc, sums) on every target without a non-own
  zero total. The repair changed exactly the targets that carried one and nothing else: 1 / 1 / 2 targets at p 0.2 and
  1 / 2 / 3 at p 0.3, for w 0.5 / 1 / 2.

## Reading
- The decided value is now backed by zero repair and 3 seeds instead of one unrepaired seed. October keeps
  `sameline_w = 0.5` for weaker banks.
- 450's rule for a STRONG (Pan-quality) bank, w 2, was not re-priced here. October ships 0.5 whatever the bank (the
  09-27 decision), so it was out of scope.
- Remaining construction difference from production: no zs pass, as in 450 §10.

## Files
- Code: `run656c.py` (450's run450c + 450 on sys.path + a local numba cache; diff `work/run656c.diff`), `ka656.py`,
  `read656.py`, `seq656.sh`.
- Outputs: `work/read656.{out,json}`, `work/seq656.log`, per-run logs `work/run656-*.log`. The stats npz files stay local.
- Peak RSS ≤ 7.6 GB per run, through gate602.
