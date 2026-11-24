# 376 — PLAN: reusable same-line calibration artifact builder + per-line value table

2026-09-10 UTC. Assignment 376-calibration-playbook. Engineering + numerical on the
source frame; ≤12 GiB RSS; `NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`. Write-only in
`experiments/376-calibration-playbook/`. No git, no builder, no submission, no reserve,
no download (all inputs on disk; nothing over 2 GB fetched).

## Hypothesis (engineering, not a score claim)
372's Jurkat artifact was a hand-adaptation of 346's H1 scripts. The construction —
ridge map from donor-consensus directions to a recipient line's own public CRISPRi
responses, fitted on non-panel targets, applied to panel targets with a full donor
consensus — is line-agnostic. A single parameterized tool can build it for any
recipient bank in under an hour on this Mac, so that on Oct 22 ("context X is line L")
the PI runs one command. The per-line within-bank held-out sign gain (372's two-split
protocol) tells the PI, per candidate line, what the map is worth before October.

## Exact inputs (paths)
- Tool estimator primitives, imported verbatim (never copied/retuned):
  `experiments/341-known-line-calibration/core341.py` (`fit`, `predict`,
  `sign_endpoints`, `interval`), `experiments/334-text-conditioned-transfer/core334.py`
  (`rms_direction`, `source_blocks`, `diagonal`, `derangement`, `endpoints`).
  Loader shape from 346 `fit346.py::bank` (exact-name row/column selection + 3×3
  scalar-index fixture).
- Production donor set: `data/donors/k562_gwps`, `data/donors/orion_hct116`,
  `data/donors/orion_hek293t` (each: `delta.npy`, `perts.txt`, `genes.txt`, `b_ctrl.npy`).
- Panel: `data/challenge/pert_counts.csv` (300 symbols, file order = canonical order).
  Axis: `data/challenge/gene_names.csv` (18,533 symbols).
- Recipient candidates (≥500 perts): `data/donors/k562_gwps` (9,786),
  `data/donors/orion_hct116` (17,766), `data/donors/orion_hek293t` (17,976),
  `data/donors/nadig_jurkat` (2,269), `data/donors/rpe1` (2,203),
  `data/donors/k562_essential` (2,020),
  `experiments/166-donor-acquisition/donors/nadig_hepg2` (2,159).
  Excluded (<500): `experiments/166-donor-acquisition/donors/frangieh_melanoma` (234),
  `data/donors/h1_2025` (254). `data/donors/194-debias-*` are derived copies, not
  independent banks — excluded, recorded.
- 341's frame for the retrieval column (read-only, no refit):
  `experiments/341-known-line-calibration/work/<line>/predictions.npz` (truth,
  crossline_diag, same_line_full) for its six lines
  (k562_gwps, orion_hct116, orion_hek293t, rpe1, nadig_jurkat, nadig_hepg2);
  `k562_essential` → n/a.
- Known answers (reproduce before any new measurement):
  K1. `experiments/375-ctx-a-arm/work/ctx_a_blend_rows_b035.npz` sha256
  `500c24ab1f7c0cc6f81e7a56f1ff361415ac845f76efc17aa2c7e9854627a8bd`
  (keys: targets, genes, matrix, beta, mode, context, norm_rule; matrix =
  372's `prediction_raw`, uncompressed `np.savez`, 10,467,702 B).
  K2. 372's axes counts: 2,006 training / 4,951 covered / 259 query / 41 missing /
  219 panel genes dropped / Jurkat∩panel = 0.
  K3. 341's Jurkat weighted-sign gain vs diagonal +0.03380218 (gate for the
  retrieval-column code path, recomputed from 341's saved predictions as in 374).

## Recipe
1. **Tool `build_calibration_artifact.py`** with CLI
   `--recipient <bank-dir> --donors <bank-dir>... --panel <csv> --axis <csv>
   --beta <float> --context <X> --out <npz>`:
   - Resolve recipient/donor paths; if recipient (resolved) is one of the donors,
     drop it automatically, require ≥2 donors remain, record
     `dropped_recipient_from_donors` in the sidecar.
   - Axes (372 §"Gene axis" generalized): training = sorted((∩ donor perts ∩
     recipient perts) − panel symbols). NOTE deviation from 372's assert: a general
     recipient may carry panel targets (e.g. K562 carries most); they are *excluded*
     from training, never asserted absent; the count is reported. Covered = sorted(
     (∩ donor genes ∩ recipient genes ∩ canonical) − training symbols − panel
     symbols). If `<recipient>/measured_mask.npy` exists with canonical length,
     intersect covered with the measured set (dormant: no in-scope bank has one;
     Pan H1 precedent from 346). Query = panel targets with a full donor consensus,
     in panel-file order; missing listed.
   - Assembly: order = training + query; donor blocks exact-name selected;
     recipient block filled on training rows only, zeros on query rows (372's
     construction; recipient query rows are never read — poison proves it).
     `rms_direction` over the stacked raw. Free the stacked raw before fitting
     (372 peaked 12.22 GiB by holding it; stay under 12).
   - Fit: `source_blocks` + `diagonal` for the cross-line gain; x = stacked
     leave-one-out donor consensuses, y = recipient response tiled over donors;
     `core341.fit` verbatim (λ = n rows); shuffled-label map with
     `derangement(n_train, 3460910)` as null. Query = donor-mean consensus on
     query rows; `matrix` = `C.predict(same_line_full, query)` (= 372's
     `prediction_raw`, the object the 373 builder blends at β itself).
   - Poison replay: recipient query rows → 1e6 noise, refit, require byte identity
     of gains/scales/betas/predictions. (Proves panel exclusion for the shipped
     configuration.)
   - Output `np.savez` (uncompressed — matches K1's 10.47 MB size) with keys in
     order targets, genes, matrix, beta, mode, context, norm_rule; scalars as
     0-d arrays (`np.array(0.35)`, `np.array("blend_rows")`, ...). Sidecar
     `<out-stem>.json`: bank shas, donor-drop note, counts, λ, poison result,
     panel-overlap count, missing targets, dropped panel genes.
   - Loader-validation gate inside the tool: delta shape vs name lists, finite
     block, 3×3 scalar fixture, no NaN/Inf in matrix, covered axis carries zero
     panel symbols (refuse otherwise — the 373 builder refuses panel support).
2. **Known-answer run**: recipient `data/donors/nadig_jurkat`, donors = production
   three in 372's order, `--beta 0.35 --context A`, out into my `work/`. Compare
   sha256 to K1. If mismatch: first compare arrays (values-equal but bytes-differ
   → save-recipe issue, iterate on recipe); if values differ → stop and report
   (do not tune the estimator; it is verbatim).
3. **Value table `value_table376.py`** (imports the tool's axes/assemble/fit as a
   module — no code duplication): for each of the 7 recipients, donors =
   production three minus recipient-if-member (October-faithful: on Oct 22 the
   donors are the production set; this keeps lines comparable and matches 372's
   number for Jurkat). Per recipient: n training, covered genes, n query, 372's
   two-split held-out (even/odd + seeded random half, seed 3760910/3760911…;
   `sign_endpoints` weighted sign, same_line_full − crossline_diag with paired
   10k-target bootstrap CI, frac improved, top-64 gain, shuffled−diag control).
   Retrieval column: blend 341's saved crossline_diag→same_line_full at β=0.35
   (374's verbatim blend), retrieval delta with `core341.interval` target CI;
   n/a for k562_essential. No refits on 341's frame; no emission; no scores.
4. **`PLAYBOOK.md`**: command sequence context-X-is-line-L → validated artifact +
   builder flag (`VCC_CTX_A_PRIOR` pattern per 375's PREREGISTRATION), checks
   (loader validation, panel-gene exclusion, config-key gate), wall time measured
   from the known-answer build.

## Controls / falsifiers
- K1/K2/K3 must pass before the value table is read. K1 is bit-identity; K2 is
  exact-count; K3 tolerance 1e-9 (374's gate).
- Shuffled-label control must be ≤ ~0 gain per recipient (as in 372: −0.023/−0.027);
  a positive shuffled gain invalidates that recipient's row (leakage/bug), not the tool.
- Donor-drop rule: covered/training counts for K562-as-recipient must reflect 2
  donors; sidecar must say so.
- Kill criterion: if K1 fails on *values* (not save bytes) the estimator
  generalization is wrong → stop, report PARTIAL with the diagnosis. If any
  recipient's held-out gain is negative, that is a measurement (reported as-is),
  not a kill — the table's purpose is per-line worth. No board/score claim is made
  from any number here.

## Budget / safety
~4 h total: tool (~30 min) + K1 build (~2 min, ~76 s in 372) + value table
(7 recipients × ~1–4 min, serial, ≤12 GiB) + retrieval (seconds) + PLAYBOOK/RESULTS.
`memory_pressure | tail -1` ≥50% before each fit-class process (observed 77% at
start; a production build ~3 GiB is running — serial fits only). Delete scratch
npz files on completion; keep the K1 reproduction artifact (small, ~10 MB) as
evidence, then delete it too after hashing (record sha in RESULTS).
