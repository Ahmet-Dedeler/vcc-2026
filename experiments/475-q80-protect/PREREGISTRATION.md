# 475 — 484's q80 finding on H1: protected down genes out of the kf shrink (preregistration, 21:59Z 09-26)
Written before any 475 number.

**Finding (cloud 484; code read here, `src/vcc2026/rc_emit_k.py` ~959–965).**
- With `protect_q` 0.8, block genes with |E0| above the block's q80 get `tgt = E0` restored.
- The mass-balance mask `side = blk & (tgt < 0)` is computed after that, so protected DOWN genes stay in it. Their
  delivered target is `kf · E0`, with kf median 0.60 on 432 b2's audit and below 0.99 on every block. Protected UP genes
  keep E0.
- This contradicts the comment ("the strongest keep their delivered LFC"), and it has been in the scored path since
  394/395.
- The fix (`rc_emit_k_475.py`) keeps protected genes out of `side`. The mass balance then falls on the unprotected down
  genes only: their kf gets smaller.

**Instrument.** 435's H1 harness, production emitter (`orig`, KA = 435 seed 0 bit for bit) vs the patched copy (`fix`),
3 seeds, deltas vs 435 at the same seed. Reported on all 85 targets and the 67 non-essential ones, with a target
bootstrap.

**Prediction (fix − orig, raw):**
- nmae −0.002 … −0.008: the strongest designed down genes get ~1.6× their magnitude. They are the likeliest reference
  DE genes, and our magnitudes run at ~0.17× the truth. The unprotected down genes shrink more, costing a little where
  their signs are right.
- pds and mse exactly 0 (per-gene totals held).
- fid, jac, reach within ±0.003 (class histograms move slightly).

**Rule.**
- A board arm (an opt-in emitter knob, then 446b + the knob with a known-answer rebuild) goes ahead only if:
  - the all-targets overall-equiv improves with target-bootstrap P(≤ 0) < 0.10 on every seed;
  - nmae improves on the non-essential subset too;
  - no ordering member moves by more than 0.005 raw.
- The H1 → board nmae transfer for magnitude-class changes is ~1.76× (433a).
- Otherwise the finding is recorded (the comment corrected) and the emitter stays as scored.

## Addendum 22:51Z 09-26 (before any 475 number): 490 prices the fix analytically, and it loses
The peer session's 490 rebuilt the production emitter's intent exactly (`recon_h1.py`: kf, MDE, flips and block counts,
max |diff| 0.0 on 85 × 3). With calls, p and reach order held at the delivered ones, 484's patch prices at nmae
**+0.0003 raw (worse)** on every seed.
- Overall-equiv −0.00009, [−0.00025, +0.00006], P(≤ 0) 0.86–0.89.
- Mechanism: kf and the MDE add are substitutes. The MDE was sized on top of the shrunk magnitudes, so un-shrinking the
  2,212 protected-down gate genes (sign accuracy 0.656) overshoots. Rebalancing also drops kf 0.53 → 0.32 on the
  unprotected down genes.
- My prediction above (−0.002 … −0.008) is therefore likely wrong. 475's run measures the |z|-realisation term that 490
  holds fixed. Expected now: nmae +0.0000 … +0.0010, ordering members within ±0.003. The rule is unchanged, so it is
  expected to fail.
