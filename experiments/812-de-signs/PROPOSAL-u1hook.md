# 812 PROPOSAL 2: H-U1PUB, the uncalled-gene hook signed by the per-target mean of the non-H1 public lines (for the judge)

Worker 812, 2026-10-16 ~08:40Z. Trigger: the lab lead's rule (a source that beats CLEAN4 by more than its CI in a paying class
gets a proposal). skill812 (`work/skill/pubmean_raw.json`, ce10e51a…; KA PASS) shows exactly one such class, u1 (uncalled,
truth-significant genes: nmae's bulk and part of reach).

## The reading that triggers it (H1, f 1, overall-eq at H1's rates; pooled = equal-weight mean of the two surfaces)
| source on u1 | vs CLEAN4 (paired cells) | vs the incumbent's delivered hook sign |
|---|---|---|
| pubmean_raw (per-target mean of 12 non-H1 lines' rows, banks averaged within a line) | u1 acc 0.590 vs 0.570 (h1test), 0.613 vs 0.580 (h1tv); **+0.0013 [+0.00007, +0.0025]** | **+0.0022 [+0.0010, +0.0034]** |
| pubmean_rms (each bank RMS-scaled) | +0.0012 [−0.0001, +0.0025] | +0.0022 [+0.0009, +0.0034] |
| pubmean_ctr (per-line centred) | +0.0007 [−0.0008, +0.0019] | +0.0013 [−0.00003, +0.0025] |

- Most of the gain over the incumbent sits at control CPM 5–20 (h1tv +0.0015 of +0.0027): genes where the key is gated to 0
  (< 10 CPM) and CLEAN4 is uncovered or below τ 0.05, so the hook leaves a near-random sign there.
- On a0 the same source LOSES to CLEAN4 (−0.0025 [−0.0043, −0.0007]), so it is a u1-only (hook) candidate and never a block
  prior.
- Tries: 3 public-mean variants plus the bulk sign. Best-of-3 shift e(3) 0.85 sd ≈ −0.0005, so raw's lower bound vs CLEAN4 is
  ≈ −0.0004 after the shift. **Against CLEAN4 it does not survive best-of-N.** Against the incumbent it does
  (+0.0005 lower).

## Why my own prior is NO-GO for October (the judge should weigh this first)
- **516 closed this family.** A Stouffer consensus of six same-target public rows on uncalled genes gained +0.016 [+0.007,
  +0.024] pool-sign accuracy on H1 but ≈ 0 on every held-out cancer fold (HCT116 +0.003, HEK293T −0.000, K562 +0.001). Its H1
  gain was read as the unwhitened-carrier effect that inverted on the board (417c: nmae +0.0089 raw).
- **My reading is the same object.** It gives H1 u1 +0.02–0.03 on paired cells.
- **Class: transfer.** It needs a held-out-line read before any DE-H1 price means anything. H1 is the line where unwhitened
  donor means flattered before.

## What I propose (only if the judge wants the question re-opened)
- **Stage 0, cheap and decisive.**
  - Run skill812's u1 class with the same three variants, leave-one-line-out, on the non-H1 DE truths that 617's audit kept
    valid. The Orion-thinned DE truths are void (616); K562 / GWPS-like truths are not.
  - Prediction: ≤ +0.005 u1 accuracy vs CLEAN4 on every held-out line (516's numbers).
  - GO to stage 1 only if ≥ 2 held-out lines show ≥ +0.01 u1 accuracy, with lower bound > 0, vs the incumbent hook.
  - Cost: < $2 EC2. It needs a DE truth per held line on the 5-CPM axis, which 516 / 613 built; to be located.
- **Stage 1, the DE-H1 price.**
  - The hook swap goes through 803's `de803.py`.
  - The swap needs a per-target hook. run612's `MDE_ADD_SIGN` is one vector for all targets, so it would be set per target `i`
    inside the `carrier_block` wrapper, as run612 already does per target for `c4kg`.
  - Arms, 3 seeds × 2 surfaces:
    - `pubhook`: pubmean_raw's sign on every uncalled gene;
    - `pubfill`: CLEAN4 where |P−0.5| ≥ 0.05, else pubmean_raw. This takes the CPM 5–20 fill-in without displacing CLEAN4.
  - Base: `c4h` (zs + hook). Prediction: Δnmae −0.002 … −0.004 raw at f 1, other members ≤ 1e-5, pds / mse 0.
  - Cost: about $5 (1 h per surface-seed on r7i.2xlarge; 12 runs in parallel on one r7i.8xlarge ≈ 1.3 h).
- **Day-of path:** D/E/F's targets' rows from the same public banks. No D/E/F data is needed, and nothing is fitted to A/B/C.

Nothing has run beyond skill812. No compute is requested until the judge rules.
