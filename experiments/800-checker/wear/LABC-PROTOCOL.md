# Using LABC from outside the lab: access, registration, wear (lab lead, 2026-10-16 11:17Z by the clock)

LABC = A/B/C's 300 targets built with the fresh-fold (def_day / dump708) recipe on the three cancer held lines, with the
day-recipe operators. It is the lab's best held-out instrument for key and direction changes on an Arc-drawn panel
(`docs/LAB-PANEL-REPORT-2026-10-16.md`). So every configuration scored on it counts, as the coordinator asked (10-16).

## Files (read-only; never write into experiments/801-lab-folds; clone with `cp -c` into your own dir if you need a copy)
| fold | dump708-format npz | sha256 | bytes |
|---|---|---|---|
| LABC/k562_gwps | `experiments/801-lab-folds/work/dumps/dump708-LABC-k562_gwps.npz` | `eba095c118837ae9a881f65d8dca1906d491733ac7d200b0427757963e0d20d3` | 84,136,473 |
| LABC/orion_hct116 | `…/dump708-LABC-orion_hct116.npz` | `e97f1fbd6bd86626533bf733c5d69aef3de2f5bf199943c03c0bc7e90dfb8337` | 101,102,719 |
| LABC/orion_hek293t | `…/dump708-LABC-orion_hek293t.npz` | `34ecd80e26bff4dc147de8d23326181898d3d0bb17381758134ed4638e60fb95` | 102,532,197 |
- Each npz has a `.json` sidecar: spec, operators (ABCp's day fits: `work/ABCp/ops/held_<held>.npz`, sha in the registry),
  rules `def`, sources, n_carried.
- S3 mirror: `s3://<your-s3-bucket>/input/lab/repo/experiments/801-lab-folds/work/dumps/`.
- Registry entry: `experiments/801-lab-folds/work/registry.json` → `LABC/<held>` (status `labc`).
- Standard-format arrays of the lab's own arms on LABC: `experiments/801-lab-folds/work/ka-labc/`.
- ABCp's dumps are in dump518 format under A/B/C rules (`work/dumps/abcp/dump518-<held>.npz`), not dump708 format.

## Registration (before any number is read)
1. Commit, in YOUR experiment dir, a registration json listing:
   - every configuration (name, the exact arm definition / env, the family it belongs to);
   - the scorer code path and sha256 (e.g. score708);
   - the three LABC dump shas above;
   - the decision rule you will apply.
2. Message the lab lead the path and commit id. The lab appends one declared row per configuration to
   `experiments/800-checker/wear/ledger.jsonl` and replies with each family's N and e(N).
3. After reading, send the results path: per configuration, per fold, overall with its target-bootstrap CI. The lab appends
   the result rows. Any configuration you call a winner is quoted at its best-of-N-corrected bound (stats800.best_of_n).

## How tries are counted
- **Surface set:** "Arc-drawn" = {ABC062, ABCp, LABC}. A read of a configuration on LABC is one try for its family on that set.
- **Prior Arc-drawn reads per family (the starting N):**
  - 707 swap (β 0.25): 3 (062 by 804; ABCp by 801/810; LABC by 801/810);
  - 705 θ rider (θ 0.25): 3 (062 by 804; ABCp by 801; LABC by 801);
  - H1/CD4 design (710 = H1 ×2.5 at the A/B/C CD4 dose): 2 (062, ABCp); not yet on LABC;
  - any other family: 0.
- So a dose curve adds its points to N. Example: 707 β 0.125 and 0.5 make the 707 family N = 5 (e(5) = 1.16).

## Caveat on what a LABC dose read can decide
LABC is ONE Arc-drawn panel. A dose chosen on it is chosen on A/B/C's panel (AGENTS 11.1c.3 / 11.1c.5). It can tell the
direction and rough size of a dose effect on an Arc-like panel. The October day doses are set on D/E/F's own panel by 708's
day-of folds.
