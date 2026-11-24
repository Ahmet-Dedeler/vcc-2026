# 458d / 458g BOARD RESULT (09-30): Arc's H1 programme as the DE block prior HARMS on A/B/C; October runs no programme

Entries: **458d** `b7UI3NRy1rdVNRRyZwmU` **0.2615288** (rank 11; "the way you do one thing is the way you do everything") and
**458g** `JL8llNLHaTyDQzm798TV` **0.2650949** (rank 9; "Don't be afraid, you can do it. (Steve Jobs)"), description empty.
Submitted 2026-09-30 by `459/arm0930.sh`. Slot 1's first upload died at 02:39Z when the Mac slept on battery; its retry
resumed on wake and published 06:01Z; slot 2 published 06:35Z. Reference 446b `UDgjOKL9gI1udAexq6ue` 0.2700428.
Bins: `PREREGISTRATION.md` (458's fid bins + the 19:40Z jac-leak addendum) and `468-prior-gate/PREREGISTRATION.md`.
Reader: `468-prior-gate/read468.py` (+ `517/predict517.py`).

| pair | overall Δ | pds | mse | nmae raw | fid raw | reach raw | jac raw |
|---|---:|---:|---:|---:|---:|---:|---:|
| 458d − 446b (programme on every target) | **−0.0085141** | 0 | 0 | +0.001549 | **−0.005545** | −0.026694 | −0.000042 |
| 458g − 458d (the gate: weak third back to CLEAN4) | +0.0035661 | 0 | 0 | −0.000779 | +0.001432 | +0.013557 | +0.000034 |
| 458g − 446b (programme on the strong two thirds) | −0.0049479 | 0 | 0 | +0.000770 | −0.004113 | −0.013137 | −0.000007 |

## Against the bins
- **458d: HARM** (fid ≤ −0.0010 raw; it is −0.0055). pds/mse bit-identical, jac within ±0.0010 (no leak). Both halves of the
  panel lose: the programme hurts on the design-strong two thirds (−0.0049 overall) and the weak third (−0.0036).
- **458g vs 458d: the gate is POSITIVE** (+0.0036), i.e. removing the programme from the weak third helps, but 458g does not
  carry against 446b. read468's combined October table: **no programme**. Decided: `programme_source off`, `prior_gate off`
  (06:36Z).
- predict517 maps the fid read to an attenuation r ≈ 0.37 of the programme's H1 sign information on A/B/C (the H1 price's
  model: random flips). Its reach loss (−0.027 raw) is larger than a random-flip model implies at the block, which fits 450h's
  lab/line-signature reading: Arc's H1 programme is partly ESC-specific, not only noisier.

## What it closes and what it leaves
- F18's open question on the challenge's own lines: a same-lab, same-chemistry programme from ANOTHER line (H1) does NOT carry
  to these cancer lines as a block prior. With 446bja (another lab's Jurkat programme, null) the block prior stays CLEAN4.
- The Arc-sourced uncalled-gene hook (446bj2zmz, READY) follows 517's HARM branch: it waits for the CLEAN4 hook's board read
  (496a, 10-05) and is read only if that transfers (predict517 at r 0.37: board nmae −0.0011 … −0.0053, still ahead of
  CLEAN4's hook). The 10-01 rule did not add it after a HARM read (464/rule1001b.py amendment).

## Addendum (2026-10-01T08:53Z): the programme was delivered at less than half dose (found by the peer's 554 calibration)
- The K-class emitter bisects the prior's log2 move (λ) down until every library class keeps its fixed 6,000-count margin
  (`rc_emit_k.py` ~1114–1165). Final metas, `output.rc_callset.prior_lambda` pooled over 900 blocks (median): **446b 1.0,
  458 0.243, 458d 0.433, 458g 0.953**, 446bja 1.0 (min 0.243), 446bj2z 1.0 (coordinator re-read, 2026-10-01T08:53Z). On H1 the same programme
  priors were delivered at λ 1.0 on 83/85 targets (MED23 0.85–0.88, PHF12 0.97–0.99; 462, 467, 517 diag jsons; corrected
  2026-10-01T14:11Z after the second peer's 593 audit: the addendum said 85/85).
- So **458d's HARM compares the programme at a pooled ~0.43 dose with CLEAN4 at full dose**: part of the −0.0085 can be lost
  steering, not wrong signs. **Corrected 2026-10-01T13:45Z (second peer's 592 verification):** 458g does NOT show the direction
  at full dose: its programme blocks are byte-identical to 458d's (468) and its pooled median 0.953 is lifted by its 252
  CLEAN4 blocks; ≥ 198 of the 648 programme blocks sat at λ ≤ 0.433 in both entries. The harm is established for the
  programme AS THIS EMITTER DELIVERS IT on A/B/C; full-dose behaviour is open.
- Likely mechanism: the programme's balanced signs (~56 % down vs CLEAN4's ~84 %) need class capacity that A/B/C's thinner
  sub-gate mass cannot supply. October is unaffected (programme off; CLEAN4 is delivered at λ 1.0). 446bja's per-context A
  dose is unknown (pooled medians only): its null is possibly confounded too. A dose-only audit rebuild
  (`VCC_RC_AUDIT_DUMP`, h5ad unchanged) would settle both; not scheduled, since no October decision depends on it.
