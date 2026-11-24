# 616 — the council's Orion folds, rebuilt on CLEAN raw truths

[21dd4c], approved by the coordinator 2026-10-02 ~18:50Z ("start after 615's stage 1 clears").

## Why
615 found that the thinned Orion harnesses' truth-side blob is the old builder's greedy control selection:
- rebuilt from raw cells, it reproduces the thinned HEK truth (plain share 0.437 vs 0.442);
- a pooled raw control set carries h1-test's level (0.060 vs 0.049).

Every Orion fold in the council's instruments reads those thinned truths:
- 480's gate tables (`480/work/gate-orion_*.npz`);
- 636's Inst;
- 639's decomposition;
- 613 / 613c / 613d / 613e;
- 612's ds612b / ds612c;
- 633's dose folds.

A shared DOWN-leaning offset in the truth's gate signs (609: Orion block gate genes 72–74 % DOWN) would favour
DOWN-leaning priors. That is a mechanism for bias, not yet a measured bias. 616 measures it.

## Design: change ONE thing, the truth's cell selection
- **Clean truths** (`truth616.py`):
  - the thinned harness's own 132 targets, 150 cells each;
  - the cells are a seeded random 150 over ALL batches (the thinned build took the first 150 in batch order);
  - 18,400 controls with 615's P′ mixture (the thinned build took the first batches' controls, then thinned them).
- Same gene axis and the same 480 gate rules (`gates480.py`, run unchanged).
- The design keys, call sets and banks depend on the targets only, so they are IDENTICAL between clean and thinned.
  Δ(clean − thinned) on any row is the truth's construction alone.
- **Instruments run UNCHANGED** inside `sandbox616.py`, a process-local sandbox:
  - the two Orion gate-table reads are redirected to the clean tables;
  - every write lands in `work/sandbox/<tag>/`;
  - each run's redirects are logged.
  - No other experiment's file is touched.

## Stages
- **S0, known answers (all must pass before any clean number is read).**
  - (a) `gates480.py` run through the sandbox on the THINNED HEK harness reproduces 480's `gate-orion_hek293t.npz`
    exactly (p, q, lfc, genes, targets).
  - (b) Identity-mode sandbox runs reproduce each instrument's committed output exactly: decomp639 → block613 (its own
    KA vs decomp639 inside), block613c, block613d, ds612b → ds612c, dose633.
- **S1:** clean truths and gate tables, HEK293T now and HCT116 after its raw extraction (615's `extract615.py`, gated).
- **S2:** clean-mode sandbox runs of the same chain.
- **S3, report.** Each row carries the thinned value, the clean value and Δ(clean − thinned):
  - 613d's hook rows (c4h − base, KGID − c4h, ungated v3 − c4h) per fold and regime;
  - 613's block instrument: the 458d retrodiction (harm sign, and size within 0.43 dose) and 392 (fid gain sign);
  - 613c's bank comparison (Nadig Jurkat alone vs CLEAN4);
  - 613e's UP-share premise for 630;
  - 633's Orion d20 rows.
  - The gate tables' DOWN share on block gate genes, clean vs thinned (609's 72–74 %).

## Coordinator's additions (binding)
1. **A CLEAN-FOLD version of 613d's 496a prediction is registered BEFORE 10-04 05:00Z**
   (`PREDICTION-496a-clean.md`). Both predictions are then on record before 10-05's read. The registered
   (thinned) prediction is never touched.
2. 613's block instrument is re-checked on the clean folds: does the 458d retrodiction still pass? It is the only
   validated block-prior instrument. 633's Orion d20 rows are added.
3. **Every row gets Δ(clean − thinned).** A sign flip on an open call goes to its owner AND the coordinator, naming the
   decision it touches. The open calls: idall3's held-out read, KGID's regret table, 613c, and 630's dose.
- Heavy stages finish before 10-04 05:00Z or pause until 10-05 18:00Z.
