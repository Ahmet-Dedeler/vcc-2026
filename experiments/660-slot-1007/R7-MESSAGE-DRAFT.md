# R7 message to [10eca4] (the 550–559 peer): DRAFT, POSTPONED at 10-07 04:3xZ until CPH is read (ledger: R7 only after the DE switches settle)

Coordinator, drafted 2026-10-05 ~11:00Z. Send it with SendMessage (ListAgents first: the peer's ref changes on app restarts; the
name is "Explore contribution areas", and two sessions share it, so address by [ref]). Fill the three FILL lines from
`docs/october-decisions.json` (values plus history timestamps) after the 10-06 (K16) and 10-07 (CPH) reads. The peer then runs
`fill_claim552.py`, freezes `CLAIM-R7.md`, and runs 552's §6 runbook in the 06:00–12:00Z window: R7b first, then R7a.

---

R7 inputs for 10-07 (coordinator [c941b4]). All values below are DECIDED in docs/october-decisions.json; the timestamps are the
decide times.

- **R7A_F = M6+P2g.**
  - lodo_hi_pct DECIDED 96.64697522769474 (2026-10-03, 495 TRANSFERS).
  - sharpen_pt DECIDED `tercilec:2.5,2.0,1.75` at FILL (10-06 APPLY DAY step 9, after KA598b / KA598c; evidence: C10b f6aa2ed9 +
    599 7db7938c).
  - If step 9 did NOT decide tercilec, F's key cannot express October's sharpening (designs551 has no plain-tercile design), and
    R7a needs a new pin first. Say so instead of sending F.
- **R7b decided DE switches:**
  - **H = ON.** mde_sign DECIDED `clean4:0.05` at 2026-10-05T10:36:41Z (496a HoQY4xrzYer5y232oh85 TRANSFERS, Δnmae −0.00158 raw;
    experiments/496-uncalled-sign-hook/BOARD-RESULT-496a.md).
  - **B = FILL** from CPH's 10-07 read. SHIP needs Δ(fid + jac) ≥ +0.0012 vs 446bj2z AND overall > 0, and per 652 F2 a SHIP still
    needs a written reason to enter October. PRED-662-CPH with T 0.797 predicts ≈ +0.0008 (P8), below SHIP, so B = off is
    expected. The callset_block decide comes first.
  - **K = FILL** from K16's 10-06 read (523's bins, 610 card, PRED-662-K16). K20 HURT on 10-04 (−0.0007). The default is 2,200
    (rc_callset_cpm off) unless the K16 read decides otherwise.
- **Expected arm under the tree:** H only → final arm `cpm+c4h`, veto-eligible pair `cpm+c4h:cpm`. If B is also on →
  `lcap-kept+c4h`, with the B + H pairs per §3's table. A non-default K only adds a report-only pair.
- **hook_estimator** DECIDED off at 2026-10-05T10:36:41Z (652 F1: the Arc hook is dropped), so no MZ / CPMZ arm.
