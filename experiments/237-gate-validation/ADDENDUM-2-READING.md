# ADDENDUM 2 to the 237 preregistration — how the `nmae` leg is READ. Written 2026-09-04 23:30 UTC, 32 min BEFORE the slot (pid 41699 fires `7a3c2d3e…` at 00:02Z), by the new coordinator session, from an external audit (GPT, 2026-09-05 01:25 Ankara) that was verified line by line.

The armed package, its time, its name and the three pre-committed overall bins at **±0.0015** (§4 of the
preregistration; ADDENDUM-243 §1) are **unchanged**. This addendum fixes what is *inferred* from them.

## 1. The overall bin does not carry the `nmae` sign — read `score_nmae` directly

ADDENDUM-243 §1 says "Δ ≥ +0.0015 means the identity mis-signed `nmae` a third time". **That does not follow.**
The board reports all six members; the other five can outweigh a fully realised `nmae` loss. Concretely, an entry
at **overall +0.0020 with an `nmae` contribution of −0.0048** (the identity's exact prediction, 229 `six229.json`
`overall_nmae_vs_incumbent` = −0.004783, overall units = Δ`score_nmae`/6) would have been declared a mis-sign by the
old `scripts/read_board_entry.py` and used to re-open 228. `scripts/read_board_entry.py --prereg 237` now prints a
**DIRECT read**: realised `nmae` contribution vs −0.004783, decided BEFORE the result:

| realised `nmae` leg (overall units) | reading |
|---|---|
| ≤ −0.0038 (≥ 80 % of the prediction) | the identity's **sign and magnitude hold on this arm**; if overall still cleared +0.0015 the gain came from the other five and the class rule under-priced the direction legs — the identity is **not** retired on this point |
| −0.0038 < leg < 0 | sign holds, magnitude refuted → the identity becomes a **discount**, not a veto (the "between" action of §4) |
| ≥ 0 | mis-signed a **third** time → retire it as a veto; re-derive the leg against the board points (the §4 "≥ +0.0015" action) |

The pre-committed overall bins are still printed and recorded verbatim; the §4 falsifier ("Δ ≤ −0.0015 AND the leg
within 20 % of −0.004783 ⇒ the slot was mis-spent") is evaluated on the direct read.

## 2. What this point can and cannot establish

229's "exact `nmae` identity" (`experiments/229-stack4-build/six229.py`) is exact in its algebra but a **proxy** in
its gene sets and folds: fixed top-89 pools per row (`POOL_N` = 89), not the scorer's Wilcoxon-significant sets,
and its `nmae` pools **all four** held-out folds (KOLF2 included) where the headline `pds` uses the three cancer
folds. So tonight's point tests the instrument **on a six-source union at g = 0.35** — one arm class, 31–44×
outside its calibration range. It does not validate or condemn every direction arm ever priced on it; the §4
actions ("re-open 228", "228 stays dead") apply to arms of this class, and 228 (per-gene gain, raw Δ`nmae` 2.2×
larger) is re-priced on the realised leg, not on a binary verdict. Likewise the "8.2 LOO-σ" quoted for 249 is a
leave-one-out prediction-error unit, not an inferential significance.

## 3. Evidence capture (operational, done before the slot)

The API resolves only the team's current entry (six earlier ids `not_found`, 253 `work/vcc-status.log`), and
`vcc submit --wait` prints the scaled members only. `scripts/arm_submit_after.sh` (exec'd by the running second-slot
waiter, pid 77532, so the edit is safe) now saves this entry's **complete** `vcc --json status` record to
`data/board-entries/<id>.json` (`scripts/capture_board_entry.py`, validated, atomic) **and refuses to submit the
second slot without it**; an independent watcher in the coordinator session does the same capture. h1union's own
record was captured at 23:23Z (`data/board-entries/iUziqE1gewfcPKzYUDaB.json`) before it disappears tonight.
