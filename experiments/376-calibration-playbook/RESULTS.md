# 376 — Reusable calibration-artifact builder (bit-identical on Jurkat) + per-line value table for October

**Status: COMPLETE — ENGINEERING_ONLY (tool validated bit-identical; 7-line
measurement table delivered; no emission, no score, no production change).**

2026-09-10 UTC. Single worker, ~1 h wall. All writes inside
`experiments/376-calibration-playbook/` (`build_calibration_artifact.py`,
`value_table376.py`, `PLAN.md`, `PLAYBOOK.md`, `work/`). No git, no builder, no
submission, no reserve, no download.

## Known-answer reproduction (own code, before any new measurement)
- **K1 — 375's β=0.35 artifact, sha256
  `500c24ab1f7c0cc6f81e7a56f1ff361415ac845f76efc17aa2c7e9854627a8bd`:
  PASS, byte-identical first attempt.** `work/repro_b035.npz` hashes to the
  recorded value; all 7 keys equal (`targets` 259 `<U9`, `genes` 4,951 `<U10`,
  `matrix` 259×4,951 float64, `beta`/`mode`/`context`/`norm_rule` scalars,
  uncompressed `np.savez`). The tool recomputed the ridge fit from scratch
  (verbatim 341/334 primitives) — values and save recipe both match.
- **K2 — 372's axes counts: PASS exact.** training 2,006 / covered 4,951 /
  query 259 / missing 41 / panel genes dropped 219 / training genes dropped
  1,613 / Jurkat∩panel 0; fit receipts subset-linearity 1.377e-14 and
  normal-equation residual 2.865e-11 reproduce 372's table exactly.
- **K3 — 341's Jurkat weighted-sign gain vs diagonal: PASS** at 0.03380218
  (recomputed 0.0338021810 from 341's saved predictions; gate 5e-9 = half-ulp of
  the 8-decimal record, disclosed — a 1e-9 gate was tried first and correctly
  refused a rounded record).
- Jurkat held-out point estimates reproduce 372 to 5 decimals on the even/odd
  split (gain +0.04204, shuffled −0.02298); the CI endpoint differs in the 4th
  decimal because the bootstrap seed differs (3760913 vs 3720910) — the CI is a
  Monte-Carlo quantity, the gain is deterministic. Random-half gain +0.04340 vs
  372's +0.04198 likewise differs only by split seed. Proceed.

## What was measured (inputs, estimator, controls)
- **Tool** (`build_calibration_artifact.py`, CLI as assigned): axes = (∩ donor
  perts ∩ recipient perts) − panel for training (sorted); covered = (∩ donor
  genes ∩ recipient genes ∩ canonical) − training − panel (sorted; optional
  `measured_mask.npy` intersect, dormant — no in-scope bank has one); query =
  panel targets with full donor consensus in panel-file order. Assembly fills
  recipient rows on training targets only (zeros on query rows); `rms_direction`;
  `source_blocks` + `diagonal` cross-line gain; `core341.fit` verbatim
  (λ = n rows); shuffled map via `derangement(n, 3460910)`; poison replay
  (recipient query rows → 1e6 noise, byte-identity required); output uncompressed
  `np.savez` with exactly targets/genes/matrix/beta/mode/context/norm_rule where
  matrix = unblended ridge prediction; JSON sidecar (bank shas, drop note,
  counts, λ, poison, overlaps). Recipient-in-donor-set → auto-drop (≥2 must
  remain), recorded in sidecar. Panel-gene support on the covered axis refuses
  to write (373-builder rule).
- **Value table** (`value_table376.py`, imports the tool — no duplicated
  estimator): 7 recipients with ≥500 perts; donors = production three
  (k562_gwps, orion_hct116, orion_hek293t) minus recipient-if-member
  (October-faithful; keeps rows comparable; matches 372 for Jurkat). Excluded:
  frangieh_melanoma (234 perts), h1_2025 (254), 194-debias-* (derived copies).
  Per recipient: 372's two-split held-out (even/odd + seeded random half,
  `sign_endpoints` weighted sign, same−diag gain with paired 10k-target
  bootstrap CI, frac improved, top-64 gain, shuffled−diag control) + retrieval
  change of the β=0.35 fixed-norm blend (374's verbatim blend) on 341's saved
  frame (read-only; n/a for k562_essential, not one of 341's six).

## Results tables (units stated)

### Per-line value table (weighted sign accuracy; diagnostic, NOT a score)
| recipient | donors | n train | covered (% canon) | query | recip∩panel | split | same−diag gain [95% paired CI] | frac↑ | top64 gain | shuf−diag |
|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|
| nadig_jurkat | 3 prod | 2,006 | 4,951 (26.7%) | 259 | 0 | even/odd | **+0.04204** [+0.03885,+0.04518] | 0.804 | +0.06418 | −0.02298 |
| | | | | | | random | +0.04340 [+0.04024,+0.04673] | 0.821 | +0.06691 | −0.02646 |
| nadig_hepg2 | 3 prod | 1,906 | 5,135 (27.7%) | 259 | 0 | even/odd | **+0.05769** [+0.05364,+0.06186] | 0.823 | +0.08616 | −0.00502 |
| | | | | | | random | +0.05567 [+0.05177,+0.05950] | 0.826 | +0.09301 | −0.00371 |
| rpe1 | 3 prod | 1,945 | 4,866 (26.3%) | 259 | 0 | even/odd | +0.09965 [+0.09416,+0.10535] | 0.867 | +0.14803 | **+0.05450** |
| | | | | | | random | +0.09813 [+0.09272,+0.10362] | 0.878 | +0.14808 | **+0.05651** |
| k562_essential | 3 prod | 1,754 | 5,840 (31.5%) | 259 | 0 | even/odd | +0.03308 [+0.03044,+0.03575] | 0.804 | +0.02921 | −0.06148 |
| | | | | | | random | +0.03167 [+0.02903,+0.03429] | 0.786 | +0.02641 | −0.06051 |
| k562_gwps | orion pair | 8,814 | 682 (3.7%) | 289 | 269 | even/odd | +0.01499 [+0.01360,+0.01635] | 0.607 | +0.02410 | −0.00051 |
| | | | | | | random | +0.01492 [+0.01357,+0.01628] | 0.602 | +0.02275 | +0.00001 |
| orion_hek293t | k562+hec116 | 8,814 | 682 (3.7%) | 260 | 294 | even/odd | +0.01323 [+0.01160,+0.01489] | 0.598 | +0.02142 | −0.01002 |
| | | | | | | random | +0.01353 [+0.01196,+0.01511] | 0.614 | +0.02096 | −0.00800 |
| orion_hct116 | k562+hek | 8,814 | 682 (3.7%) | 263 | 290 | even/odd | +0.00698 [+0.00557,+0.00844] | 0.524 | +0.01178 | −0.01156 |
| | | | | | | random | +0.00661 [+0.00516,+0.00807] | 0.524 | +0.01056 | −0.00955 |

### Retrieval change at β=0.35 on 341's frame (1,024 targets; diagnostic, NOT pds)
| line | retr diag | retr β=0.35 | Δ [target 95% CI] | reading |
|---|---|---:|---|---|
| nadig_jurkat | 0.8199 | 0.8250 | **+0.00507** [+0.00038,+0.00985] | only resolved-positive row; 0.35 safe here |
| orion_hct116 | 0.7579 | 0.7586 | +0.00064 [−0.00422,+0.00558] | flat |
| k562_gwps | 0.8251 | 0.8238 | −0.00134 [−0.00548,+0.00274] | flat |
| orion_hek293t | 0.7675 | 0.7557 | −0.01181 [−0.01962,−0.00359] | resolved cost |
| nadig_hepg2 | 0.7797 | 0.7697 | −0.00997 [−0.01497,−0.00487] | resolved cost |
| rpe1 | 0.8175 | 0.7867 | −0.03088 [−0.03883,−0.02298] | large cost |
| k562_essential | — | — | n/a (not in 341's six) | — |

## Verdict against the PLAN.md kill criterion
No kill fired: K1 is bit-identical (values AND bytes), K2 exact, K3 within the
record's own precision. The tool is a validated drop-in for 372's hand
adaptation; PLAYBOOK.md gives the Oct-22 command sequence with per-line β
guidance. Two findings constrain how the table may be quoted:
1. **RPE1's +0.100 is ~half shared-component, not per-target information.**
   The shuffled-label control is +0.055 (positive both splits); the defensible
   per-target number is same−shuffled ≈ **+0.045**. Same code reproduces
   Jurkat's negative shuffled exactly, so this is RPE1's response geometry (large
   target-shared direction: even permuted labels teach the mean), not a bug.
   Never quote +0.100 without the floor.
2. **Covered axis collapses to 682 genes (3.7%) whenever training saturates a
   small gene universe.** All three 2-donor rows (K562/HCT116/HEK293T as
   recipient) drop 8,814 training symbols as genes from a K562-bottlenecked
   ~8k universe. The artifact still builds and validates, but touches little of
   the axis — the October value of a K562/HCT116/HEK293T map is narrow by
   construction, independent of its sign gain.

## Limits
- Same-bank, held-target, cross-study logic throughout (as in 372/346): gains
  are within-bank held-out diagnostics, not A/B/C or D/E/F forecasts, and no
  board conversion is claimed. Donor-count confound across rows (3 vs 2 donors)
  is inherent to the drop rule — the orion/K562 rows are not directly comparable
  in magnitude to the 3-donor rows.
- k562_essential calibrates (+0.033, clean null) but is October-irrelevant: the
  real panel contains zero essential genes (§8.1g), so its row validates the
  tool, not a play.
- HepG2's shuffled control (−0.005/−0.004) is near-zero rather than clearly
  negative; its gain stands but with a thinner null margin than Jurkat's.
- Ridge λ/features/architecture are 341's, untuned; β=0.35 is 375's operating
  point, not re-optimized per line (the retrieval column shows per-line β work
  is owed before wiring HepG2/HEK293T/RPE1).
- A production build (~3 GiB) ran throughout; all fits serial, peak 11.2 GiB
  (value table) / 8.4 GiB (K1 build), free ≥77% throughout.

## Exact commands to reproduce
```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # require >=50% free
# K1: tool reproduces 375's artifact bit-identically
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/376-calibration-playbook/build_calibration_artifact.py --recipient data/donors/nadig_jurkat --donors data/donors/k562_gwps data/donors/orion_hct116 data/donors/orion_hek293t --panel data/challenge/pert_counts.csv --axis data/challenge/gene_names.csv --beta 0.35 --context A --out experiments/376-calibration-playbook/work/repro_b035.npz
sha256sum experiments/376-calibration-playbook/work/repro_b035.npz  # expect 500c24ab…8bd
# Value table (44 s)
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/376-calibration-playbook/value_table376.py
uv run python -c "import json;print(json.load(open('experiments/376-calibration-playbook/work/value_table.json'))['recipients'][0]['heldout']['evenodd'])"
```

## Peak memory and wall time
| step | wall | peak RSS | notes |
|---|---|---:|---|
| K1 build (fit + poison, no held-out) | 19.6 s | 8.41 GiB | vs 372's 75.7 s / 12.22 GiB (raw freed before poison) |
| value table (7 recipients, serial) | 43.9 s | 11.16 GiB | slowest bank k562_essential 12.3 s; orion rows ~1.6 s (682-gene axis) |
| retrieval column (341 frame, read-only) | <2 s | <1 GiB | inside the 43.9 s; no refit |
| total worker wall | ~1 h | — | incl. reading, PLAN, tool, gates, PLAYBOOK, RESULTS |
