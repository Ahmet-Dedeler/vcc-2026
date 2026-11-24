# 481 proposed patches (NOT applied; production files untouched)

Each patch is a unified diff against the file as it was on 2026-09-26 ~21:30Z, plus a `.patched` copy for reading.
Apply only after the running chains have restored the builder / finished with the scripts (AGENTS §10), then commit
with the known-answer lines below.

| file | finding | what it changes | known answer (measured here) | still needed |
|---|---|---|---|---|
| `F2_ntclevel417_cp46.diff` (`experiments/417-bulk-sharpen/ntclevel417.py`) | F2 | each NTC construct's jackknife × n / cells_per_pert (manifest); refuse constructs < 100 cells; print sizes when unequal | `test_cp46_481.py`: A/B/C Cp46 = 14.439751080335085 bit for bit; 800-cell / 200-cell / uneven constructs corrected to 1.0003 / 1.0015 / 0.9998 of the truth (current rule 0.500 / 2.003 / 1.105) | `ntclevel417.py --ka` on an A/B/C artifact (must still PASS: the factor is exactly 1.0 there) |
| `F1_F5_def_day.diff` (`scripts/def_day.sh`) | F1, F5, labels | `DEF_DAY_CD4_THRESHOLDS=relaxed` pass-through + an explicit stop message when 047's CD4 gate fires; contexts read from the manifest (no hard-coded D,E,F); a `def_postbuild.py` step BEFORE the EXIT restore, `DEF_DAY_INCUMBENT_META` (default 446b), allowed diffs + `bulk_load` (+ `sameline`, `de_prior_gate` only when set) | `bash -n` clean; manifest parse gives `D,E,F` / `D=1,E=1,F=1` exactly (so every rehearsal's commands are unchanged); the postbuild call passes 41/41 on 465's artifact with the builder as pinned during that run (40/41 with the driver's current list) | one rehearsal on a panel whose CD4 gate fires (efp: union 186) with `DEF_DAY_CD4_THRESHOLDS=relaxed`, and one on ef-bundle-r that must retrace 452/465 |
| `F7_def_sameline_programme_rows.diff` (`scripts/def_sameline.py`) | F7 | `--programme-source` refuses < 150 programme rows or < 50 % gene coverage (`--min-programme-rows`) | the patched copy rebuilds 465's `sameline.npz` byte for byte (sha256 `b4eb2a6d…`); refuses on a panel inside the 2025 H1 panel (0 rows) | none beyond the chain's usual KA |
| `cp46_fix481.py`, `test_cp46_481.py` | F2 | the fix as a function and its test rows | as above | — |

Not patched (recorded in RESULTS.md): the prior gate (472's `le1` is the fix in flight), the panel-median anchoring of
the carriers' scale (F6), `dayof438.py`'s q grid (1.25 / 1.5 only, shipped q 1.75), the zero-gene audit not gating.
