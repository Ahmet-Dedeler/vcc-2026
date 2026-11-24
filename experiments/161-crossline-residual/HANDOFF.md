# 161 — cross-line residual: HANDOFF (run transferred 2026-09-02, no result yet)

**Status: scaffolding complete + inventory MEASURED. No model has been fitted. No GO/NO-GO.**
No submission, no `cell-eval2`, no existing file modified. Everything is inside this dir.

## Scripts written (all `uv run python <x>.py`, from this dir)

| file | what it does | run? |
|---|---|---|
| `common.py` | Donor-table loader + line-GROUP map. Key invariant: `sum(expm1(b))==TS` and `sum(expm1(b+delta))==TS` exactly, so a donor table is a legal 156-style delivered-space surface. Groups `k562_gwps`+`k562_essential` as ONE line (K562) and the three CD4 states as one; splitting either leaks (159: same-line pair reads +0.2712 vs +0.03–0.16 cross-line). Flags `rpe1`,`kolf2` as blob-carrying. | n/a |
| `inventory.py` | §1 deliverable: per-table census, pairwise shared targets split panel/non-panel, lines-per-target histogram, candidate training surfaces, leakage audit. | **DONE** → `work/inventory.json` |
| `build.py` | Builds a mass-consistent multi-line surface. Restriction to a shared gene subset is done in PROFILE space and the `log1p(TS*p)` frame rebuilt, so the truth stays exactly on 156's mass manifold (asserted `<1e-9`). CPM gate carried from the FULL axis, not recomputed. | **DONE** → `work/surface-{core5,core4,wideh1}.npz` (13 s, 6.5 GB peak) |
| `residual.py` | The model. NOT YET RUN — plumbing untested, expect small bugs. | **NOT RUN** |

## What the inventory says (real numbers, §1 deliverable)

- **`core5` = K562/RPE1/Jurkat/HCT116/HEK293T: 1,883 shared targets × 4,980 genes, and 0 of the
  1,883 are on the challenge panel.** So the closure in `docs/CURRENT_STATE.md` ("rpe1,
  nadig_jurkat, k562_essential cover zero panel targets") is not merely inapplicable to
  training — **the entire trainable pool is disjoint from the deployment set by construction.**
- `core4` (drop blobby RPE1): 2,006 × 5,390, 0 panel. `wideh1` (adds H1 as a held line):
  226 targets (206 non-panel) × 7,240 genes.
- Pairwise non-panel shared targets: HCT116↔HEK293T 17,409; K562↔Orion ~8,900; K562↔Jurkat
  2,252; K562↔RPE1 2,187; Jurkat↔Orion ~2,050; K562ess↔Orion ~1,800; H1↔Orion ~225.
  **2,129 targets are covered by ≥5 distinct lines, 1,892 of them non-panel.**
- CD4 (3 states) and KOLF2 are **100% panel** — they carry zero non-panel shared targets, so
  they are inference-only donors and cannot train anything.

## Design decisions already made (keep or overturn deliberately)

1. Fit `y = R_L − C`, where `R_L = E_L − mean_t(E_L)` and `C` = equal-weight mean of
   mean-removed donor LINES. Bar = `C` alone.
2. `yhat` is **centred across targets** before use. A target-independent addition is the
   metric's zero *and* a `pds` hub — that is 160's whole result. Uncentred is scored only as
   the demonstration. Consequence: a CONTEXT-ONLY model predicts exactly 0 after centring, so
   159's context claim must be tested through the `cons × context` INTERACTION block, not
   through context main effects. Feature families in `residual.py:FAMILY` are C / T / X.
3. `div` (recipient↔donor control-profile divergence) is one scalar per fold = line identity in
   disguise. It is in the C family and must be ablated separately.
4. Scoring: `deliver()` = gate → production `mass_consistent_effect` → `render` → 156
   `endpoints`. `alpha` is one deployment-available constant per fold, identical across arms.
5. Regimes R1 held-line / R2 held-target-seen-line / R3 both; permutation null + 40-cluster
   target bootstrap wired into R1.

## Next steps

Smoke-test `residual.py` on `wideh1` first (226 targets, seconds), confirm the two known-answer
rows (`zero effect` ≈ baseline-ish, `ORACLE residual` → pds 1.0 / cos 1.0), then run `core5`.
Feature families not yet buildable: per-gene control **variance** and **detection rate** in the
recipient — donor tables store only `b_ctrl`, so these need cell-level h5ads and are the main
untested part of 159's context claim. Laptop is sufficient; no AWS needed to reach a verdict.
