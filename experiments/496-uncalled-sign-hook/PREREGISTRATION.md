# 496 — the uncalled-gene sign hook (478's route) on H1: preregistration (2026-09-27T05:40Z, before any 496 result)

**Change (one variable).** The emitter's MDE add on complement genes (not head, not block: the uncalled genes) takes its
direction from the CLEAN4 per-gene prior where `|P_down − 0.5| ≥ τ` (`mde_add_sign`, `patch_emitter496.py`); the base
LFC, the calls, the |z| design and the bulk are unchanged, and fu/fd re-balance the mass. Harness: 435's H1 harness
(`run496.py`, local patched emitter `k_emit496.py`), coarse MDE (0.1, 0.1, 0.3), 3 seeds, vs 435's production runs.

**Known answer (must pass before anything is read):** the hook emitter without a sign source reproduces 435's production
coarse run (seed 0) stats bit for bit.

**478's analytic prediction** (exact emitter reconstruction, calls held): τ 0.05 → nmae −0.0015, reach +0.0017 … +0.0035,
overall-equiv +0.0007 … +0.0011; τ 0 → nmae −0.0025, reach ≈ 0, overall-equiv +0.0005 … +0.0008. The real run adds what
the analytic path held fixed: leaked calls (unsuppressible, clip) and the |z| realisation.

**Bins (τ 0.05, pooled over 3 seeds, target bootstrap):**
- falsifier: pds or mse moves at all → a bug (the hook must not touch the bulk);
- leak check: unsuppressible calls per target more than 2× the base's, or clip fraction up by > 0.05 → the delivery leaks;
- **go** (build the production hook + a board arm): overall-equiv ≥ +0.0005 and nmae ≤ −0.0010 raw, no leak;
- **weak**: overall-equiv in (0, +0.0005) or nmae in (−0.0010, 0): a rider only, never its own slot;
- **no**: overall-equiv ≤ 0.
- The reach part gets no board credit before a board read (H1-BOARD-CALIBRATION rules 3–4); a go is priced on nmae.
