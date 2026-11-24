# 603: analytic ceiling of the rank-neutral nmae lever on H1 (council sweep #1; the arm's owner is d4cd5f)

2026-10-02 ~12:00Z, council session [4d0207]. Analysis only. Ran through the shared gate `602/gate602.sh`, ~7 s per
script, under 10 GB. State: 446bj2z's DE construction on H1, i.e. 435's production-emitter stats WITH the zs pass, 3
seeds. Calls are 521's reconstructed C = head | block. Calls, order and signs are pinned and ONLY lfc magnitudes change.
Scored with 521's exact `nmae_rows` at f 1 and f 0.35 (PowerSurface). Overall = −Δnmae × 1.64787 / 6.

**Mechanism under test (21dd4c's forensics):** nmae reads arithmetic means; calls and reach read Wilcoxon ranks. U counts
only pred-vs-control pairs, so a pred cell can move anywhere inside its gap between adjacent control values. The top
cell's gap is unbounded, so an UP mean can be raised with U, calls and order unchanged. DOWN has no free direction.

## Results (mean over 3 seeds; seed sd ≤ 0.002 on every row)

UP-only rule |lfc| ← max(|ol|, K · MDE_g), MDE_g = 0.24/√(mean control count per cell):

| regime | K 0.5 | K 0.75 | K 1.0 | K 1.25 | K 1.5 |
|---|---:|---:|---:|---:|---:|
| f 1, Δnmae [overall] | −0.0056 [+0.0015] | −0.0111 [+0.0030] | **−0.0128 [+0.0035]** | +0.0021 [−0.0006] | +0.0335 [−0.0092] |
| f 0.35 | −0.0068 [+0.0019] | −0.0137 [+0.0038] | −0.0219 [+0.0060] | **−0.0291 [+0.0080]** | −0.0310 [+0.0085] |

- A 0.5 log2 cap only trims f 1's overshoot (K 1.25: −0.0039). A 1.0 cap changes nothing.
- **The same rule on DOWN-delivered genes HURTS at every K** (K 1: +0.024 at f 1, +0.025 at f 0.35). Leave DOWN as
  delivered. This agrees with d4cd5f's 640 (graded magnitudes, UP-only −0.0092, DOWN +0.0035).
- The donor key as magnitude (|cap15 key| / ln 2) changes nothing (≤ 1e-4): it sits below the delivered |lfc|.
- UP-delivered genes are 39–41 % of the gate (by gate |l| mass too). Sign accuracy on them is 0.544 at f 1 and 0.595 at
  f 0.35; on DOWN-delivered genes 0.574–0.578.
- Oracle magnitudes at our signs: UP −0.028 / −0.047, BOTH −0.098 / −0.090 (480's −0.0985 reproduced).

## Reading
- The lever is real but regime-dependent. f 1 wants K ≤ 1; f 0.35, the regime 611 found to retrodict the board's
  delivery pairs, wants K ~1.25–1.5. A board pair **K 1.0 / K 1.25, UP only** brackets it.
- Size on H1: +0.0035 … +0.0080 overall. With 611's delivery multiplier (nmae × 1.2–1.6), the board estimate is about
  **+0.004 … +0.013**. d4cd5f's 640 and 21dd4c's 610 dose curve price the same pot at +0.003 … +0.005 with other
  sign sources and magnitude rules. Three independent pricings agree that it is a mid-sized, transferable rider (scorer
  mechanics), not +0.019.
- The rivals' large pure-nmae steps were mostly catch-up from worse bases (d4cd5f: post-step levels 0.84–0.94, against
  our 0.897), so their size does not transfer to us one for one.
- Still open (d4cd5f's real-emitter build): whether delivery is truly rank-neutral once per-cell CPM coupling is
  counted (U, calls and pds/mse pinned to ≤ 1e-6).

Files: `ceiling603.py` (rules mdeK / key / oracle, UP / DOWN / BOTH), `dose603.py` (UP dose × cap grid),
`work/ceiling603.json`, `work/dose603.json`.
