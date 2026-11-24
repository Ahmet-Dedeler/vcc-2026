# 364 — Amplitude break-even vs gate sign accuracy (PLAN)

## Hypothesis
Uniform amplitude ×a>1 is net-negative at today's delivered |lfc|-weighted sign
accuracy (p ≈ 0.588) because the wrong-sign ~41% of gate mass pays (1+a) while
the right-sign mass collects at most its overshoot headroom (343/347: e≈0.34–0.38
< 0.5 everywhere). If a sign lever raises p, the wrong-sign tax base shrinks and
a may turn net-positive. This experiment maps overall(p, a) and the break-even
line p*(a) / a*(p) by exact arithmetic on saved arrays — no emission, no scorer.

## Inputs (all read-only, nothing outside `experiments/364-amplitude-breakeven/` written)
- `experiments/323-current-family-validation/work/stats-cap15.npz` (delivered p/q/lfc, 85×10,779)
- `.../reference.npz` (reference q/lfc/keep/ref_lin/genes/targets)
- `.../per-target-cap15.csv`, `.../metrics-cap15.json` (known answers)
- `.../h1-anchors.json` (H1-local baseline/replicate for scaled conversion)
- `experiments/317-operating-point/work/sweep317.json` (recorded amp120/145/200 raw pds/mse deltas)
- Functionals: fid/nmae/reach re-implemented in own code following 361's
  `fid361.py` conventions (read-only reference); cross-checked to 1e-9.

## Estimator
1. **Known answer**: own fid/nmae/reach reproduce 0.49755452 / 0.95083646 /
   0.25434726 (means; per-target max-diff vs `per-target-cap15.csv`).
2. **Sign model**: gate = {ref q<0.05} ∩ finite-nonzero ref lfc ∩ not-own-target;
   eligible targets = gate ≥ 10 (expect 58, as 343). Pooled |l|-weighted accuracy
   p0 measured. For each target p ∈ {0.60, 0.62, 0.65, 0.70, 0.75}: 20 seeded
   draws; per draw shuffle all wrong-sign gate pairs, flip prefix (sign only,
   |lh|/q/p untouched) until pooled accuracy ≥ p; record achieved p.
3. **Amplitude**: on gate genes L = a·s_new·|lh|, a ∈ {1, 1.2, 1.45, 2}; off-gate
   unchanged. Exact recompute of nmae (eligible-target mean), fid (85-mean),
   reach (scored-mean, same NaN set) per (draw, p, a). Verify fid bit-invariant
   across a (sign(a·x)=sign(x)) and reach order-invariance across a.
4. **pds/mse**: 317's measured raw deltas at a (all three grid a's directly
   measured — no interpolation needed): amp120 Δpds≈0.0/Δmse+0.21827,
   amp145 −0.000560/+0.551690, amp200 −0.001261/+1.521028. Assumed independent
   of p (see Assumptions).
5. **Pricing**: H1-local scaled Δ via h1-anchors spans; jac carried at 0 (q
   untouched ⇒ call sets fixed); overall(p,a) = mean of six scaled deltas vs
   current delivered (p0, a=1). Break-even a*(p) and p*(a=1.2, 1.45) by linear
   interpolation; slope = overall pay per +0.01 accuracy at each a.

## Controls / falsifiers
- Known-answer gate (§1): if any of the three misses by >1e-9, STOP and report.
- fid/reach invariance check across a: max |Δ| must be exactly 0 for fid;
  reach reported (expect 0; any nonzero is a finding about tie-breaks).
- Analytic cross-check: per-(p,a) breakeven sign from 188's e=ΔM/ΔU>½ on the
  pooled gate; must agree in sign with the interpolated a*(p).
- Sanity: overall(p0, a=1.2) nmae leg should be ≈ 317's measured −0.0029 raw
  direction (ours is sign-model-free at fixed p0; 317's includes call-set drift
  — same sign, ours cleaner).

## What kills the idea
If overall(p, a) ≤ 0 for all a>1 at every p ≤ 0.75 (i.e. even at 75% weighted
sign accuracy no uniform amplitude pays), the uniform-amplitude route stays
closed and only stratum-targeted (non-uniform) amplitude could matter — which
347 already closed. The deliverable is then the p* number (>0.75 or nonexistent).

## What result reopens it
A p* ≤ ~0.63 for ×1.2 with overall margin meaningfully above the H1-local noise
would say tonight's sign levers (346: +0.033 weighted-sign pilot with own-line
labels; 362/363 pending) are aimed at a live prize. Report which p the running
levers would need, plainly, with the caveat that 346's gain used the
recipient's own training labels.

## Limits (pre-registered)
H1 development panel only; not an A/B/C forecast. overall(p,a) is a LOWER bound
on the true gain (pds improvement from better signs not modeled). Assumes
317's h1-test amplitude deltas transfer to the 323 surface as raw deltas.
