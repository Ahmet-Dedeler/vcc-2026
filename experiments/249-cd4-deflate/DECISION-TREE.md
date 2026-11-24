# Pre-decided: what happens after `cd4svd5`'s board result (written 2026-09-04 18:20 UTC, ~9 h before the result)

Incumbent `h1union` 0.2118310 (rank 3). `cd4svd5` prereg: overall +0.0026 … +0.0040 (upside +0.0086), board raw
`pds` +0.0086 at the 0.7× transfer, raw `mse` ≈ 0.8191. It fires after the 237 gate test scores; read both with
`scripts/read_board_entry.py <id>` and 249's prereg falsifiers 1–5.

## Branches
**A — raw `pds` ≥ +0.0086 (the 0.7× prediction realised or better) and overall ≥ +0.0026.** The `pds`-first class
transfers at ≥ 0.7×. October: ship the PINNED basis (`VCC_CD4_DEFLATE_BASIS=V5_abc`, 251), update `SHIPPED_ENV` and
re-run 240's row 1 through the driver. Next slot (09-06): the k ladder, **k = 10** (246: `pds` +0.0124, price flat
to +0.0075 at ×1; 250: off A/B/C the plateau sits at k 10–15) — one variable, preregistered on the realised `pds`
transfer, NOT the `+partial` composition (254: its own increment is +0.0009 … +0.0018, under the gate at ≤ 1.0×).
**B — raw `pds` in [+0.0040, +0.0086) and overall ≥ +0.0015.** Transfer 0.33–0.7×. The class is real but weaker than
`h1union`'s 0.70×. October: pin the basis (250: never negative on 4/4 panels). Next slot: hold; re-price k = 10 and
`+partial` on the realised factor; ship only what clears 1.5× the gate on it.
**C — raw `pds` < +0.0040 (falsifier 1 fires) but overall ≥ 0.** The mechanism does not transfer; the arm scored on
its cosine legs by luck. Retire the deflation CLASS on the A/B/C board; October: `VCC_CD4_DEFLATE_BASIS` stays UNSET
(250's pinned rows are PLOLO-only evidence and the board just refuted the instrument's `pds` reading). Do not spend
a slot restoring the incumbent.
**D — overall < 0 or raw `mse` worse by > +0.0030 (falsifier 2).** Same as C, plus: the cosine cost was
under-read; re-derive 243's class factor for a "re-shape occupied rows" class from this point before pricing any
future re-shaping arm.
**E — delivered amplitude moved > 0.5% (falsifier 5, read off the package's own `A_c`).** The arm was not shape-only;
re-classify under the amplitude gate before any reading of A–D.

## What does NOT change on any branch
`VCC_ALPHA` 0.982, `VCC_CD4_TARGET_F` 0.34, `VCC_H1_UNION_G` 0.28, `rc_callset`, the operator. The `cd4svd5+partial`
package (254, `eea5ed18…`) is NOT submitted on any branch without its own preregistration priced on the realised
`pds` factor. Never re-submit a scored package; rank is informational until Oct 22.

## The pairing with the 237 result
If 237 lands Δ ≥ +0.0015 (the `nmae` identity retired as a veto) AND branch A/B fires here, 228's per-gene gain
(+0.0028 … +0.0042 on direction legs, dead only on `nmae`) becomes the second candidate for 09-06 and must be
re-priced against `cd4svd5` as the incumbent (208's rider rule).

## Addendum — precedence and a HOLD branch (written 2026-09-04 23:32 UTC, before either result; packages and times unchanged)

An external audit (GPT, 2026-09-05 01:25 Ankara, verified) found two holes in the branches above, and
`scripts/read_board_entry.py --prereg 249` now applies the tree exactly as written here:

1. **Precedence, hard failures first: E > D > A > B > C.** D ("overall < 0 OR raw `mse` worse by > +0.0030") can
   co-occur with A/B/C's conditions; it wins. E is read off the package (R249 F5: −0.409 %, inside the bound), so it
   does not fire tonight.
2. **HOLD for every uncovered outcome.** A requires overall ≥ +0.0026 and B overall ≥ +0.0015, so "raw `pds` transferred
   but overall under its bin" (e.g. `pds` +0.0090 raw with overall +0.0010) matched nothing. Such outcomes, and any
   outcome where prereg falsifier 3 (overall outside [−0.0010, +0.0100]) fires on an otherwise A/B result, are
   **HOLD**: reassess on all six members, no next-slot escalation, no October decision from this point alone.
3. **k = 10 is NOT automatic on branch A.** 255 measured its increment over `cd4svd5` at `pds` +0.00054 on 1/3 PLOLO
   folds, priced −0.0005 … −0.0002 overall, with a further −0.46 % amplitude (R255 §0, §2.1). Branch A's "next slot:
   k = 10" stands only as a candidate that must clear its own preregistration on the realised `pds` factor; the same
   holds for `+partial` (254, +0.0009 … +0.0018, under the gate).
4. Both entries' complete records are captured to `data/board-entries/<id>.json` before the second slot fires
   (`scripts/arm_submit_after.sh`, `scripts/capture_board_entry.py`); read them from the files, since the API resolves
   only the current entry.
5. **The pairing condition above ("If 237 lands Δ ≥ +0.0015 (the `nmae` identity retired as a veto)") is SUPERSEDED**
   (23:50Z, before either result): the identity is retired only by the 237 entry's DIRECT `nmae` read (realised
   `nmae` contribution ≥ 0 against −0.004783; E237/ADDENDUM-2-READING.md §1), never by its overall total. 228 becomes
   a 09-06 candidate only on that direct verdict (or on a PARTIAL verdict, re-priced with the discount) AND branch A/B here.
