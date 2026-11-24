# 616 — the council's Orion folds rebuilt on clean raw truths

[21dd4c], 2026-10-02. PLAN.md has the design: the same 132 targets and 150 cells per target. Only the truth's cell
selection changes: the old greedy selection vs random cells over all batches with 615's P′ control pool.

## S0 — known answers: all PASS (`work/ka616.json`)
- `gates480.py` through the sandbox on the THINNED HEK harness reproduces 480's committed gate table bit for bit
  (max |Δ| 0 on p, q, lfc).
- Identity-mode sandbox runs of 14 instruments reproduce all 17 committed outputs byte for byte:
  - decomp639;
  - 613, 613b, 613c, 613d, 613d3–d6, 613e, 613e2;
  - ds612b, ds612c;
  - dose633.

## S1 — the clean gate tables: the thinned Orion DE truths were the selection artefact
| gate table | downstream DE calls | median per target (q75) | targets ≥ 10 | genes called in ≥ 50 % of targets | DOWN share |
|---|---|---|---|---|---|
| HEK293T thinned (480) | 37,748 | 141 (390) | 125 / 132 | 44 | 0.562 |
| **HEK293T clean** | **774** | **1 (3)** | **13 / 132** | **0** | 0.565 |
| HCT116 thinned (480) | 11,309 | 68 (113) | 127 / 132 | 12 | 0.598 |
| **HCT116 clean** | **712** | **2 (4)** | **18 / 132** | **0** | 0.643 |

- **94–98 % of the DE calls every Orion-fold instrument graded against were artefact.**
- Clean, each fold has ~13–18 scoreable targets, below K562's fold. The Orion folds carry no usable DE truth.
- 615's deeper HEK truth (220 cells, 15k) agrees: median 1.

## S2 — HELD (coordinator, 19:06Z)
- The clean-mode instrument re-runs refused to start (`sandbox616.py` returns 3 without `work/RELEASE-S2`), so no clean
  instrument number exists.
- The coordinator re-ruled from the dependence list instead (`AUDIT-orion-dependence.md`; NEXT-SLOT-DECISION 19:10Z).
- With ~13–18 scoreable targets per fold, a clean re-run would mostly measure noise. S2 stays held unless a specific
  question needs it.

## Memory note (AGENTS 9)
- block613 peaked at 5.3 GB in the identity chain, which ran ungated (each instrument takes ~10–20 s).
- The HCT116 extraction went through gate602.
