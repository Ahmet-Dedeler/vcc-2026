# 474 — October-path hardening: 481's patches applied with their known answers (queued 21:57Z 09-26)

Applies `experiments/481-october-audit/patches/` (F2 `ntclevel417.py`, F7 `def_sameline.py`, F1/F5 `def_day.sh`)
after chain471b, with these known answers:

| patch | known answer | pass criterion |
|---|---|---|
| F2 Cp46 × n/cells_per_pert | `test_cp46_481.py`; `ntclevel417 --ka` on 469dr's relabelled artifact | A/B/C Cp46 lines identical to the pre-patch baseline (`work/ntclevel_ka_pre.log`: KA PASS, Cp46 14.4398 / 14.8057 / 16.4397); --ka PASS |
| F7 programme-row minimum | 465's `def_sameline` call | 465's `sameline.npz` rebuilt byte for byte (`b4eb2a6d…`) |
| F1/F5 (a) | efp (CD4 union 186) without the variable | stops at setup with the explicit CD4-gate message |
| F1/F5 (b) | 465 retraced on ef-bundle-r, same switches | final output sha = 465's `6113937f…`; postbuild PASS |
| F1/F5 (c) | efp with `DEF_DAY_CD4_THRESHOLDS=relaxed` | DONE; postbuild PASS; every final ρ ≤ 1 |

Any failure restores the three pre-474 files. The coordinator commits after reading the log.

## Results (02:38Z 09-27): ALL PASS
| check | result |
|---|---|
| F2 KA | `test_cp46_481.py`: A/B/C Cp46 14.4397510803… bit for bit. `ntclevel417 --ka` on 469dr's relabelled artifact: PASS, Cp46 lines identical to the pre-patch baseline |
| F7 KA | 465's `sameline.npz` rebuilt byte for byte (`b4eb2a6d…`) |
| (a) efp strict | stops at setup: "STOP: CD4 extraction gate fired (only 186 panel targets map across any condition; require 290). Check the mapping …, then rerun with DEF_DAY_CD4_THRESHOLDS=relaxed in a fresh work dir" |
| (b) 465 retrace | **byte for byte** (final output `6113937f…`); converged at D=1.2198, E=1.2842, F=1.2496 after 3 levelled builds (= 465); **postbuild PASS** vs 446b |
| (c) efp relaxed | DONE; converged at D=0.9831, E=1.0324, F=1.0663 after 1 levelled build; ρ 0.9808 / 0.9766 / 0.9745; **postbuild PASS** |

Committed: `def_day.sh` `91f350f0`, `def_sameline.py` `19596af1`, `ntclevel417.py` `a4b6b4cc`.
