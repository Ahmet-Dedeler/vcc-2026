# 469 — the October path with Arc's programme AND the strength gate (469dr), end to end on A/B/C

**Result (2026-09-26 19:53Z): PASS.** `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full DEF_DAY_PRIOR_GATE=q33 def_day.sh` (flags `pt`,
q 1.75, w auto) on 180's relabelled bundle ran end to end. Its config differs from 461dr's (the same path without the
gate) only in `de_prior_gate`. This is the end-to-end test of `def_day.sh`'s `DEF_DAY_PRIOR_GATE` switch
(`5f9bfa26`, committed with this result).

## Gates (`gate469.py`, `work/gates469.log`)
- **Re-level:** converged after 2 levelled builds at D=1.0859, E=1.1383, F=1.1662, the same factors as 461dr's. Every
  final ρ is ≤ 1 (0.9793 / 0.9794 / 0.9795). The gate does not touch the bulk.
- **Programme source:** h1_2025_full for D/E/F; flags auto:pt; no LODO.
- **Prior gate:** logged by `def_day.sh` and by every build: threshold 2 in each context.
  - Programme on 214 / 217 / 217 targets, CLEAN4 on 86 / 83 / 83.
  - These are the counts of 458g on the unrelabelled A/B/C.
- **Same-line artifact:** equal to 461dr's.
  - Every array is identical (NaN-aware). Its sha differs, 683f5452 vs 765fa07f, only because the provenance records the
    run's own identity.json path.
- **Stack dump:** arrays identical to 461dr's (the file differs in its env record).
- **Zero-gene audit:** 0 blocks.
- **Config diff vs 461dr:** exactly `['de_prior_gate']`.
- Relabelled to A/B/C, 120,000 cells each. Package `917c217b5d215926661b9f76ff7fca69589a2a1c7510af06f45e3ba63328d249`,
  novelty rc 0.
- The builder was restored to `693fde15` after the run.

## Use
- A 10-01 candidate under `rule1001b.py`, eligible only if 458g reads positive (G = pos), as "the sonicator adopted a
  cat". It is read against 453dr (09-28 s1) and 461dr.
- If the gate is adopted, 469dr is the October configuration's board rehearsal.
