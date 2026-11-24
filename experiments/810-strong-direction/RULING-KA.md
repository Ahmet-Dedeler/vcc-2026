# Lab lead ruling on 810's known-answer tolerance (2026-10-16 08:55Z, before any aggregate of the GO run exists)

**Facts.** On ABC/orion_hek293t the β-0 identity `sgu0 == inc` failed on `a` (⟨P,P⟩ per target and frame, an unnormalised sum)
by 1.364e-12 absolute against the frozen 1e-12 limit. On that fold max|a| is 375.04 (804's inc array), so the relative error is
3.6e-15, i.e. machine precision. It comes from re-adding the four parts (their sum differs from F.fin by 3.6e-14 on this fold,
the largest of 22). All other KA keys on all 22 folds pass (worst 6.3e-13); pds, pdsl and K are exactly 0.

**Ruling.** The β-0 identity KA on the unnormalised sums `x` and `a` is read with a RELATIVE tolerance: |d| / max|ref| < 1e-12.
Normalised and rank quantities (pds, pdsl, K, r) and every 804-reference KA keep the absolute 1e-12. The KA's intent ("at β 0 the
construction IS the incumbent") is unchanged and is met to 3.6e-15. This changes the gate, not the candidate or its price.

**What 810 does.** (1) Patch only the assert in `price810.py` to the rule above (new sha256, recorded in RESULTS.md and
SHA256SUMS; nothing else changes). (2) Re-run only the failed fold ABC/orion_hek293t with the patched file (one small launch). Keep
the 21 completed folds as they are, because their arithmetic did not depend on the assert, and record that all 21 also pass the new
rule. (3) Then aggregate all 22 and apply decide810 exactly as frozen, with G1–G4.
**Disclosure carried into RESULTS:** one fold's arm lines (ABC/orion_hek293t: sgu1 r +2.81 %, dpds −0.0263; sgu1all r +15.4 %,
dpds −0.0533) were seen while diagnosing. No selection depends on them. From now on the KA gate runs before any summary line
is printed.
