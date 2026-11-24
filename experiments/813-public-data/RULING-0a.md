# Lab lead ruling on 813's stage 0a (2026-10-16 10:00Z by the clock)

**Facts (813 RESULTS §7).**
- Stage 0a ran through labcloud.sh: run `20261016T093745Z-s0a813`, $0.094, metadata reads only.
- It wrote `freeze813.json` with sha256 `2f6e5169…ceba`.
- KA-SIM (ii) failed on the injected guide-load design: 0.200 against the ≤ 0.10 bar.
- 813 then reported three defects in its own implementation, before any expression row was read:
  1. `sim813.py` caps the guide-count stratum at 12 or more. PROPOSAL-2 §2 preregisters strata (n_guides_total, pool, rep)
     with no cap.
  2. The single-guide twin keeps one guide per cell, so each target has about 4.7× fewer incidences. PROPOSAL-2 preregisters
     the same incidences per target.
  3. KA-SG's two statistics are dominated by low-count log bias on single-guide rows. On an unbiased simulated design β
     centres at 0.10 and dsh at 8.1, so their bounds cannot transfer to real data.

**Ruling.**
1. **Defects 1 and 2 are departures of the code from the frozen design, not tuning.**
   - C3 forbids re-tuning after a failed gate. It does not forbid running the preregistered estimator that never ran. No real
     row exists, and nothing about the TeloHAEC expression data has been seen.
   - One corrected stage-0a run is allowed (≤ $0.15), as `freeze813-v2.json`. Its code diff from `sim813.py b2c8a20e…` is
     limited to: (a) uncapped strata exactly as PROPOSAL-2 §2; (b) the matched-incidence twin exactly as PROPOSAL-2 specifies.
   - The diff and the new sha256 go in RESULTS.md before the run.
   - The KA-SIM bars stay as frozen: (i) ≥ 0.9, (ii) ≤ 0.10 on both designs, (iii) > 0.9999.
   - If v2 fails any bar, that is the gate's answer: record it and STOP. There is no third 0a.
2. **The failed v1 freeze stays on the record as a negative with its mechanism.** Coarsening guide count (a capped stratum)
   leaks a load-dependent offset into the per-target estimate: 0.200 of a median effect under the injected design. That is a
   finding about the estimator class, and it binds any later estimator here (no coarsened k).
3. **Defect 3: KA-SG takes G5's fallback; no new statistic.**
   - Building KA-SG in count space would be a new statistic, which would need the judge. It is not authorised.
   - KA-SG is recorded as **uninformative at this depth** (G5), and RESULTS states the consequence: on real data only
     G0-BLOB guards against k-dependent offsets. The TeloHAEC folds carry that caveat in their registry entry.
   - The v2 freeze records KA-SG as `"status": "uninformative (G5)"` and drops its bounds.
4. **Then stop.** The lab lead commits freeze813-v2.json's sha before 0b (G1, C3). Every other GO condition is unchanged.
