| section | probe | blind (pre-252) | fire (252) | verdict |
|---|---|---|---|---|
| identity | `Operator(095_union_maxeig_8.036)` | silent | OperatorGuardError | **info**  |
| identity | `Operator(donors_union_maxeig_6.257)` | silent | accepted | **info**  |
| identity | `Operator(donors_union_pruned_maxeig_5.895)` | silent | accepted | **info**  |
| identity | `Operator(TAMPERED_held_orion_hct116)` | silent | OperatorGuardError | **info**  |
| identity | `Operator(union.npz, name='union_pruned')  [SWAP]` | TypeError | OperatorGuardError | **FIRE** blind: TypeError; fire: OperatorGuardError: SWAPPED OPERATOR: data/donors/whitener_union.npz hashes to 5608e807230dea04f |
| identity | `Operator(union_pruned.npz, name='union_pruned')  [PASS]` | TypeError | accepted | **PASS** blind: `name=` not expressible (TypeError); fire: accepted |
| identity | `Operator(TAMPERED, name='held_orion_hct116')  [SWAP]` | TypeError | OperatorGuardError | **FIRE** blind: TypeError; fire: OperatorGuardError: SWAPPED OPERATOR: experiments/252-harness-guards/work/probe/held_orion_hct11 |
| identity | `Operator(FOLDS[k562_gwps].op)  [PASS]` | silent | accepted | **PASS** identical |
| identity | `Operator(FOLDS[orion_hct116].op)  [PASS]` | silent | accepted | **PASS** identical |
| identity | `Operator(FOLDS[orion_hek293t].op)  [PASS]` | silent | accepted | **PASS** identical |
| identity | `Operator(FOLDS[kolf2].op)  [PASS]` | silent | accepted | **PASS** identical |
| e064 | `load_frozen_operator(TAMPERED, expected=sha256_file(TAMPERED))  [TAUTOLOGY]` | silent | ValueError | **FIRE** blind: silent; fire: ValueError: ~/Code/virtual-cell-challenge-2026/experiments/252-harness-guards/work/probe |
| e064 | `load_frozen_operator(genuine, expected=sha256_file(genuine))  [PASS]` | silent | accepted | **PASS** identical |
| e064 | `load_frozen_operator(genuine, expected=<pinned 7cecbf7f...>)  [PASS]` | silent | accepted | **PASS** identical |
| e064 | `load_frozen_operator(TAMPERED, expected=<pinned 7cecbf7f...>)  [FIRE, old check]` | ValueError | ValueError | **FIRE (pre-252 too)** blind: ValueError; fire: ValueError: ~/Code/virtual-cell-challenge-2026/experiments/252-harness-guards/work/p |
| e064 | `load_frozen_operator(genuine, expected='')  [empty expectation]` | silent | ValueError | **FIRE** blind: silent; fire: ValueError: ~/Code/virtual-cell-challenge-2026/experiments/062-cd4-post-whiten/cloud-out |
| axis_trap | `pess axis reconstructed == 233 json axis_sha256 (index digest)` | silent | accepted | **info**  |
| axis_trap | `n_cols` | silent | accepted | **info**  |
| axis_trap | `Operator(233 pess op)` | silent | accepted | **info**  |
| axis_trap | `apply(pess op, ABC wcols)  [TRAP]` | silent | OperatorGuardError | **FIRE** blind: silent; fire: OperatorGuardError: AXIS MISMATCH: self-described:pess__held_orion_hct116 (experiments/233-october- |
| axis_trap | `apply(pess op, PESS wcols)  [PASS]` | silent | accepted | **PASS** identical |
| axis_trap | `apply(pess op, ABC wcols, check_axis=False)  [named bypass, = blind]` | <absent in blind> | accepted | **info**  |
| preserve_endpoints_zero | `operator` | silent | accepted | **info**  |
| preserve_endpoints_zero | `apply(40pct coverage)  [FIRE]` | silent | OperatorGuardError | **FIRE** blind: silent; fire: OperatorGuardError: ON-AXIS PRESERVE UNDEFINED: the median row norm on the operator axis is 0.0 --  |
| preserve_endpoints_zero | `apply(60pct coverage)  [PASS]` | silent | accepted | **PASS** identical |
| preserve_endpoints_zero | `apply(full pre)  [PASS]` | silent | accepted | **info**  |
| preserve_endpoints_zero | `endpoints(real truth)  [PASS]` | silent | accepted | **info**  |
| preserve_endpoints_zero | `endpoints(all-zero truth)  [FIRE]` | ZeroDivisionError | accepted | **FIRE** blind: ZeroDivisionError; fire: sgn/sgn_w = NaN (mean_endpoints skips them) |
| preserve_endpoints_zero | `mean_endpoints([real, zero])` | <absent in blind> | accepted | **info**  |
| preserve_endpoints_zero | `cd4 coverage` | silent | accepted | **info**  |
| preserve_endpoints_zero | `solve_fraction(cd4, 0.34)  [PASS]` | silent | accepted | **PASS** identical |
| preserve_endpoints_zero | `telohaec coverage` | silent | accepted | **info**  |
| preserve_endpoints_zero | `solve_fraction(telohaec, 0.34)  [FIRE]` | silent | SilentZeroError | **FIRE** blind: silent; fire: SilentZeroError: SILENT ZERO REFUSED: solve_fraction(sizing='all_rows') would return fraction 0.0 f |
| preserve_endpoints_zero | `compose218.solve_fraction_carrying(telohaec, 0.34)` | silent | accepted | **info**  |
| preserve_endpoints_zero | `solve_fraction(telohaec, 0.34, sizing='carrying')` | <absent in blind> | accepted | **info**  |
| preserve_endpoints_zero | `solve_fraction(telohaec, 0.34, sizing='legacy_all_rows')` | <absent in blind> | accepted | **info**  |
| preserve_endpoints_zero | `gse295214_h1 coverage` | silent | accepted | **info**  |
| preserve_endpoints_zero | `solve_fraction(gse295214_h1, 0.34)  [FIRE]` | silent | SilentZeroError | **FIRE** blind: silent; fire: SilentZeroError: SILENT ZERO REFUSED: solve_fraction(sizing='all_rows') would return fraction 0.0 f |
| preserve_endpoints_zero | `compose218.solve_fraction_carrying(gse295214_h1, 0.34)` | silent | accepted | **info**  |
| preserve_endpoints_zero | `solve_fraction(gse295214_h1, 0.34, sizing='carrying')` | <absent in blind> | accepted | **info**  |
| preserve_endpoints_zero | `solve_fraction(gse295214_h1, 0.34, sizing='legacy_all_rows')` | <absent in blind> | accepted | **info**  |
| preserve_endpoints_zero | `solve_fraction(cd4, 0.34, sizing='carrying')` | <absent in blind> | accepted | **info**  |
| preserve_endpoints_zero | `solve_fraction(cd4, 0.34, sizing='legacy_all_rows')  [= default where per_unit > 0]` | <absent in blind> | accepted | **info**  |
| cache | `build_fold(unchanged)  [PASS: hit]` | silent | accepted | **PASS** identical |
| cache | `build_fold(sources changed to ('orion_hek293t',))  [FIRE: must miss]` | silent | RebuildAttempted | **FIRE** blind: silent; fire: RebuildAttempted: REBUILD ATTEMPTED: the cache was not used |
| cache | `build_fold(op changed to held_orion_hek293t)  [FIRE: must miss]` | silent | RebuildAttempted | **FIRE** blind: silent; fire: RebuildAttempted: REBUILD ATTEMPTED: the cache was not used |
| cache | `build_fold(PM_BLEND 0.5 -> 0.25)  [FIRE: must miss]` | silent | RebuildAttempted | **FIRE** blind: silent; fire: RebuildAttempted: REBUILD ATTEMPTED: the cache was not used |
| cache | `MIN_CPM alternative that changes the gate` | silent | accepted | **info**  |
| cache | `build_fold(MIN_CPM 10 -> 20)  [FIRE: must miss]` | silent | RebuildAttempted | **FIRE** blind: silent; fire: RebuildAttempted: REBUILD ATTEMPTED: the cache was not used |
| cache | `build_fold(restored)  [PASS: hit]` | silent | accepted | **PASS** identical |
| cache | `fold_cache_status(unchanged)` | <absent in blind> | accepted | **info**  |
| cache | `fold_cache_status(sources changed)` | <absent in blind> | accepted | **info**  |
