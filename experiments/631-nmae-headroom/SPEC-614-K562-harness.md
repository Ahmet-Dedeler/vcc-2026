# Spec for 614 — a real-emitter DE harness on K562 (from [d4cd5f], 2026-10-02; owner [21dd4c])

Purpose: the instrument 534 flagged as missing ("the DE channel has no cancer-line held-out surface"). It prices DE
delivery levers (630 UP pass, 631 DOWN magnitude, hook forms) with the REAL emitter on a clean cancer line.

## Truth
- `experiments/480-nmae-instrument/k562harness480.py` streams `data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad`
  (805 MB, dense gzip, 1.99M cells x 8,248 genes) into a CSR h5. It keeps 18,400 random NTC cells plus up to N cells per
  A/B/C panel target with >= MIN cells; env K562_N_CELLS / K562_MIN_CELLS / K562_OUT.
- Its outputs (`work/k562-gwps-panel.h5`, `-c250.h5`) are no longer on disk, so regenerate. `work/gate-k562_gwps.npz`
  (267 panel targets) still exists. Suggested: N 250, MIN 100.

## Reference (E.load_reference format)
- `experiments/323-current-family-validation/evaluate323.py: reference()` (line 148) builds `reference.npz`: keep, cs,
  tie, ref_lin, emit_* , libraries, and per-target sums / jk / p / q / lfc. It reads `run323.REAL` (h1-test) and asserts
  H1's frozen axes (`verify()`, `R.axes()`). A K562 variant points REAL at the K562 h5, takes axes from the harness, and
  drops the H1 asserts.

## Design key: two options
- (fast, a proxy) fold480's donor consensus: `ds612.surface_forms` builds cons0 from `load_donor_rows` over the surface's
  donors. For k562_gwps those are the two Orion lines. It is not the production key: no cap15, no whitening, no
  multi-donor weights.
- (faithful) `run323.py: stage()` (lines 219–275) builds the production key: MultiDonorModel over DONORS = (k562_gwps,
  orion_hct116, orion_hek293t), `apply_perturbmatch_k562_replacement` (PerturbMatch's K562 data), CD4, and the arms,
  written to `stage.npz`. lib521's key reads `effects.npz` (EFF, cap15). For K562 truth: donors = Orion only, SKIP the
  PerturbMatch K562 replacement (leak), keep CD4, and refit or prune the whitener so K562's targets are excluded.
- Recommendation: run end to end on the proxy first. The retrodiction gate decides whether it is adequate; upgrade to
  the faithful key only if the gate fails.

## run496's H1-specific inputs (`experiments/496-uncalled-sign-hook/run496.py`)
- STAGE (:72) `323/work/stage.npz`; OP_H1 (:71) `323/work/whitener_h1_pruned.npz`.
- DONORS (:73); CLEAN4 (:77), where K562 drops k562_gwps and k562_essential (Nadig Jurkat + HepG2 remain; 636's
  BANK_FORBIDDEN).
- CARRIERS (:76) `429-zs-production/work/carriers416.npz`. Carrier libraries should be re-drawn from K562 control
  libraries (median ~11k vs H1's).
- PRIOR_PROD (:88) `396-jurkat-tilt-A/work/prior_ctx_challenge.npz` (A/B/C context prior); for K562 use
  CLEAN4-minus-K562.
- Reference via `run400` (H.E.load_reference → 323's reference.npz) and lib521's REF / EFF / ST435.
- Emitter switches: RUN426_EMITTER=prod (production rc_emit_k) or RUN496_EMITTER=hook. Production flags as 446bj2z:
  prodem, MDE 0.1/0.1/0.3 coarse, zs dn_k250_zs, zero repair. Dump blocks (`--dump 85`-style) so 630 / 631 passes can
  re-apply (`rnup630.py --src`).

## Gate (coordinator) before any lever is priced
611-style retrodiction of the board's known DE pairs on the K562 harness: zs vs pre-zs, 433a coarse MDE, and the
call-set pair 446bj2zcp. Bands are fixed in a committed file before running.

## Then price
- 630 v3 `b_k100_q20r3p10`;
- 631 DOWN pass / v7 variants;
- hook forms (c4h / c4kg / idall3).
Use 521's regimes, and report d20-rescaled rows as 633 does, because MDE units depend on depth.
