| id | tree | command | rc | ok |
|---|---|---|---:|---|
| A1 641 --check after 538+554+590 | T0 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python ~/Code/virtual-cell-challenge-2026/experiments/641-depth-row/patch_depth641.` | 1 | PASS |
| A2 642 --check after 641 | T1 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python ~/Code/virtual-cell-challenge-2026/experiments/642-hek293-contingency/patch_` | 1 | PASS |
| A3 642 refuses a second apply | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python ~/Code/virtual-cell-challenge-2026/experiments/642-hek293-contingency/patch_` | 1 | PASS |
| A4 644 skip --check after 642 | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python ~/Code/virtual-cell-challenge-2026/experiments/644-runbook-amendments/patch_` | 1 | PASS |
| A5 644 skip: def_day.sh syntax after the patch | T3 | `bash -n scripts/def_day.sh` | 0 | PASS |
| A6 644 skip: the knob is logged at start | T3 | `grep -c -F log "644: DEF_DAY_SAMELINE_SKIP= scripts/def_day.sh` | 0 | PASS |
| D1 def_day.sh syntax | T2 | `bash -n scripts/def_day.sh` | 0 | PASS |
| D2 def_day.sh logs the WIDE_LIBS override at start | T2 | `grep -c -F log "554: DEF_DAY_ALLOW_WIDE_LIBS scripts/def_day.sh` | 0 | PASS |
| D2 def_day.sh logs the LOW_DOSE override at start | T2 | `grep -c -F log "554: DEF_DAY_ALLOW_LOW_DOSE scripts/def_day.sh` | 0 | PASS |
| D2 def_day.sh logs the DEPTH override at start | T2 | `grep -c -F log "641: the depth row scripts/def_day.sh` | 0 | PASS |
| D2 def_day.sh logs the NONCCLE override at start | T2 | `grep -c -F log "642: the HEK293 contingency scripts/def_day.sh` | 0 | PASS |
| D3 def_day.sh without arguments stops at the bundle argument | T2 | `bash scripts/def_day.sh` | 1 | PASS |
| P1 def_preflight --verify (A/B/C: both rows silent) | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/def_preflight.py --verify` | 0 | PASS |
| P2 def_preflight --help | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/def_preflight.py --help` | 0 | PASS |
| I1 def_identity --verify (prints 590's min top-1) | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/def_identity.py --verify` | 0 | PASS |
| I2 def_identity with DEF_DAY_ALLOW_NONCCLE_REF=hek293 (642 stage 2, A/B/C not called) | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/def_identity.py --bundle data/challenge --contexts A,B,C --out <tmp>/identity_h` | 0 | PASS |
| I3 def_identity refuses an unknown reference | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/def_identity.py --bundle data/challenge --contexts A` | 1 | PASS |
| S1 def_sameline --help shows --min-top1 | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/def_sameline.py --help` | 0 | PASS |
| S2 def_sameline on the stage-2 identity file (no HEK293 call -> A on via Nadig Jurkat, B/C none) | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/def_sameline.py --bundle data/challenge --identity <tmp>/identity_hek.json --ou` | 0 | PASS |
| B1 def_postbuild --help shows --rc-audit | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/def_postbuild.py --help` | 0 | PASS |
| B2 def_postbuild --rc-audit on a healthy dump (433a: no 554 row) | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/def_postbuild.py --meta submissions/433a-cap15-partial50-clean4q80-pertscale3-b` | 1 | PASS |
| K1 october_recipe command emits the decided day (tercilec, LODO pct, zs, clean4:0.05) | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/october_recipe.py command --bundle data/challenge --work <tmp>/recipe-work --na` | 0 | PASS |
| K2 def_day.sh carries the 598 tercilec pass-through | T2 | `grep -c -F tercilec scripts/def_day.sh` | 0 | PASS |
| K2 def_day.sh carries the 612dd KEYGATE pass-through | T2 | `grep -c -F KEYGATE scripts/def_day.sh` | 0 | PASS |
| K2 def_day.sh carries the 630dd DE_UPRN pass-through | T2 | `grep -c -F DEF_DAY_DE_UPRN scripts/def_day.sh` | 0 | PASS |
| K2 def_day.sh carries the 644 skip pass-through | T2 | `grep -c -F ${SLSKIP:+--skip-contexts "$SLSKIP"} scripts/def_day.sh` | 0 | PASS |
| K3 def_sameline --help shows --skip-contexts (644) | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/def_sameline.py --help` | 0 | PASS |
| K4 def_sameline carries 604's identity programme source | T2 | `grep -c -F a.programme_source == "identity" scripts/def_sameline.py` | 0 | PASS |
| K5 the builder carries VCC_DE_MDE_SIGN_KEYGATE | T2 | `grep -c -F VCC_DE_MDE_SIGN_KEYGATE scripts/make_submission_pm_cd4.py` | 0 | PASS |
| K5 the builder carries VCC_DE_UPRN | T2 | `grep -c -F VCC_DE_UPRN scripts/make_submission_pm_cd4.py` | 0 | PASS |
| K5 the builder carries tercilec | T2 | `grep -c -F tercilec scripts/make_submission_pm_cd4.py` | 0 | PASS |
| R1 october_recipe.py show | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/october_recipe.py show` | 0 | PASS |
| R2 october_recipe.py check --allow-pending | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/october_recipe.py check --allow-pending` | 0 | PASS |
| C1 capture_board_entry.py --help (live repo: no bundle patch touches it) | live | `~/Code/virtual-cell-challenge-2026/.venv/bin/python scripts/capture_board_entry.py --help` | 0 | PASS |
| F1 531 readdef.py --help | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python experiments/531-cloud-final-phase-tooling/readdef.py --help` | 0 | PASS |
| F2 637 readdef_factorial637.py (no args prints usage) | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python experiments/637-final-phase-factorial/readdef_factorial637.py` | 1 | PASS |
| F3 637 tests | T2 | `~/Code/virtual-cell-challenge-2026/.venv/bin/python experiments/637-final-phase-factorial/test_factorial637.py` | 0 | PASS |

Not verifiable yet:
- 612 KEYGATE / 630 UPRN / 604 identity prior VALUES: the code paths are verified (K2-K5) but none is in the recipe yet (their reads are 10-08+); a decided value must be re-checked through K1's command
- the D/E/F bundle itself: `vcc datasets download controls` after Oct 22 (and the CLI upgrade)
- a real def_day.sh run on a bundle (~2 h; the coordinator's rehearsals cover it; 644's skip is stub-tested)
- readdef637.diff dispatch into readdef.main (UNAPPLIED until the final-phase tooling freeze)
