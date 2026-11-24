# 535 — the call-propensity DE block on the board: arm 446bj2zcp (preregistration, 2026-09-29 16:45Z)

Adopts 528's `PREREGISTRATION-draft.md` (arm, construction, predictions, bins) as written, before the arm is built and
before any 09-30 entry exists. Slot: fixed later by a day's rule, never before that rule is written. **Amended 2026-09-30
06:38Z, before any 10-01 entry exists:** 10-01 slot 2 (empty after the 09-30 HARM read), via `arm1001cp.sh`, if build535
passes; the bins are unchanged.

## The arm
**446bj2zcp** = 446bj2z (entry `be90smZUDcIV7j5MgQ03`, 0.2728260) + 528's emitter/builder patches +
`VCC_RC_CALLSET_BLOCK=prop:experiments/528-callprop-block/work/callblock528_abc.npz:04fb24ec…:<swap 1>:<steer kept>`: the
CallSet's 2,200-gene CPM block is chosen per context by log10 control CPM + 22 × excess dropout (median-rescaled), from that
context's own controls; only the block's CPM part keeps CLEAN4 steering. DE-only: pds and mse must equal 446bj2z's to ≤ 1e-6
raw. Expected config diff vs 446bj2z: exactly `rc_callset_block` + `de_sign_emitter_sha256`.

## Bins (528's, unchanged; primary Δ(fid + jac) overall vs 446bj2z)
POSITIVE ≥ +0.0012 · PARTIAL +0.0004 … +0.0012 (October rider only if Δjac raw ≥ +0.002 and Δnmae ≤ +0.001 raw) · NULL
−0.0004 … +0.0004 · NEGATIVE ≤ −0.0004. Falsifiers: pds/mse > 1e-6 raw; Δjac ≤ 0 with Δfid at or above prediction;
|Δnmae| > 0.003 raw. Predicted Δ(fid + jac) +0.0020 … +0.0025 overall; reach and nmae reported, not binned.

## Caveats written before the read
- **Forking paths:** 528 chose the `kept` steering form after seeing the `coupled` and `decoupled` forms on H1 (its own
  statement). The H1 CI [+0.0033, +0.0070] is therefore optimistic; the board read is the test.
- **Not additive with 523's K 1,600** (528: the kept swap at K 1,600 prices −0.0012 … −0.0027 at f ≤ 0.5): this arm stays at
  K 2,600. If 523's K pair reads POSITIVE, the call set's size and composition are competing uses of one lever; the combined
  arm is not built without its own H1 price.
- On 446bj2zmz (if the hook transfers first) the same knob prices +0.0085 on H1 (528 §4.5): the variant 446bj2zmzcp is a
  separate arm with its own read.
