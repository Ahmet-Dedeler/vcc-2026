# READ-FIRST — VCC 2026, primary sources only (v2, 2026-09-04 17:00Z)

Key: A=`AGENTS.md` · BL=`docs/BOARD-LEDGER.md` · RB=`docs/DEF-RUNBOOK.md` · DP=`scripts/def_pipeline.py` · MB=`scripts/make_submission_pm_cd4.py` ·
SL=`submissions/SUBMITTED_LEDGER.json` · E<n>=`experiments/<n>-*/` · BR<n>=E<n>`/BOARD-RESULT*.md` · R<n>=E<n>`/RESULTS.md` · P237/AD=E237`/PREREGISTRATION.md`,`/ADDENDUM-243.md` ·
P249=E249`/PREREGISTRATION.md` · LB=E253`/work/ranks.json` (16:53Z) · VS=E253`/work/vcc-status.log`.

## 1. Board position
- `h1union`, id `iUziqE1gewfcPKzYUDaB`: 0.21183096, rank 3/673 (LB:3-7; VS:2). Leader `cqawesome`, model `m0904b`, 0.22454757 at 14:18Z; gap +0.01272 (LB:6,53-57). 2nd `Aginglab.com` 0.21287742 (LB:66).
- Members raw / scaled / world rank (LB:9-49): `pds` 0.83572/+0.74537/8 · `mse` 0.81726/+0.17957/7 · `nmae` 0.93389/+0.11074/58 · `fid` 0.50748/−0.01648/107 · `reach` 0.27655/+0.22294/8 · `jac` 0.04096/+0.02884/1.
- The four 2026-09-03/04 arms, one variable each, from a15 0.1740865 (BL:28; ids in SL):

| arm | overall | Δ | variable |
|---|---:|---:|---|
| `rc_hybrid` | 0.1834190 | +0.0093325 | emitted p-order, pseudobulk bit-identical (BL:92) |
| `rc_callset` | 0.1918638 | +0.0084448 | call set top-400 `|H|` ∪ top-2,200 CPM, `n_pred` 4,480→2,599 (BL:93) |
| `cd4f034` | 0.2049196 | +0.0130557 | `VCC_CD4_TARGET_F` 0.10→0.34 (BR210:5,8) |
| `h1union` | 0.2118310 | +0.0069114 | `VCC_H1_UNION_G`=0.28, all six members up (BR222:4-6) |

- Scorer deterministic: a byte-identical resubmission reproduced every member (BL:8-9); `pds`/`mse` identical to 16 s.f. across a15/`rc_hybrid`/`rc_callset` (BL:100-101).

## 2. The shipped configuration
- Env `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28` (DP:84; RB:118-120); `VCC_UNION_SOURCES` unset = `gse295214_h1, h1_2025` (DP:85; MB:946). Flags `--min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0` (DP:114-115); `--emitter rc_callset` (RB:208); head 400 / cpm 2,200 / z_cap 1.9 (RB:209); γ 0.65, λ 0.0 (RB:205-206). Operator `data/donors/whitener_union_pruned.npz` `e2f65ecf…` (RB:89); D/E/F re-fit `union_pruned_def` (RB:207).
- Arms absent from 191/202: `cd4f034` +0.0130557 and `h1union` +0.0069114, sum +0.0199671 (RB:123-127). Package `submissions/222-h1union-g028.h5ad` `8df4dadd…` 12,648,807,800 B, `.vcc` `f13ce40b…` 2,824,785,920 B (RB:93-95); the driver's defaults rebuild it byte-for-byte (R240:17-21).
- Four inert flags added today; each unset reproduces staging `abeb0749…` and package `8df4dadd…`:
  - `VCC_CD4_SIZING=carrying`+`VCC_CD4_TARGET_G=0.3443655278649641`: 247 `a17b177` (MB:697-702); KA R247:15. October: the driver sets it on every non-A/B/C bundle (DP:98,111; RB:28-29).
  - `VCC_UNION2_SOURCES`+`VCC_UNION2_G`: 248 `2632b92` (MB:1000-1001); KA R248:14-16. October: no runbook row, driver-blind (R251:375-378); its rider NOT BUILDABLE (R248:417-421).
  - `VCC_CD4_DEFLATE_K`: 249 `be029c83` (MB:740); KA R249:14-15; built `cd4svd5` (R249:348). October: never — a day-of re-fit reads ≤ 0 on both low-coverage draws (RB:227-229; R250:5-8).
  - `VCC_CD4_DEFLATE_BASIS=V5_abc`: 251 `c321f76` (MB:793-807), artifact `data/donors/cd4_deflation_basis/V5_abc.npz`; KA R251:16-22 (rehearsal 40/40). October: unset until the `cd4svd5` point; if board raw `pds` rises against 0.83572 ship this, never a re-fit (RB:224-230); driver-blind (R251:375-380).

## 3. Rules that bind every agent
- 2 submissions/day/team, reset 00:00Z, one in flight (A:49). A slot buys information, never rank; never re-submit a scored package (A:639-653,1347-1350). Only the LAST final-round entry counts (RB:8-10).
- §14: `-m` a random joke, `-d` EMPTY, never confirm/deny a guess, every public surface (A:1779-1790).
- One variable per run (A:1343); preregistration committed before the build (BR210:7-8; BR185:7-8).
- Price all six members or state why one cannot move (A:1487-1496). Gate > +0.0019 shape-only, > +0.0089 amplitude-changing; ordering noise sd 0.004502 overall between digest-identical arms (A:1474-1478); no `reach` toll (A:1428).
- Known-answer row in every screen, abort if it does not reproduce (A:1375-1382); diff the built config against the incumbent, exactly the intended differences (A:1360-1362).
- Path + sha, never a description (A:1497-1498). `fid`/`jac`/`reach` normalisation unreadable on harnesses (A:1507-1508).
- Never `git add -A` with agents running (A:1541-1547). Memory gate `memory_pressure` ≥ 40% free AND free+inactive ≥ 20 GiB AND swap free ≥ 400 MiB, heavy jobs serial (A:1552-1557); one `cell-eval2`, ~25 GB (A:601-603).
- `uv run python` (A:1268). Launch via `scripts/claude_subagent.sh <exp-dir> <mandate>` (A:1602).

## 4. Instruments
| instrument | reads | validated on | known failures |
|---|---|---|---|
| board | all six | byte-identical resub (BL:8-9) | only the current id resolves, six earlier ids `not_found` (VS:4-9). Read a package's cosine MEASURED off its own `A_c`/`B_c`, never the floor (R243:289-290) |
| h1-test `data/harness/h1-test.h5ad` (A:599) | DE members, relative | two board points within 0.014 raw (A:552-554) | oracle-leaked with `h1_2025` as donor, 85/85 (A:581-584) = every arm on top of `h1union` (MB:852); a21 2/4 DE signs wrong (A:556-559); real scorer: `mse` sign-inverts on the direction arm, `fid` on the call-set arm (R221:20-23); order-arm `reach` 0.769× (R221:444), not 1.290 (BR185:50) |
| PLOLO `E162/plolo.py` | `pds` of geometry arms | 6/6 signs, LOO RMS 0.001494 (R162:9-10) | PLOLO→board cosine factor is a CLASS constant: re-dose of occupied rows k_cos 1.36 (1.31–1.42), first source into empty rows 0.94 (0.94–1.00), on `cos_agg`, 3 cancer folds, measured cosine (R243:44-47,286-290). Mirror returns 0.0 below 50% coverage (R223:20; A:871-879); 252 owns the fix (E252/MANDATE.md:1) |
| 197 inversion `cos ≥ (1.023611−raw mse)/2.074495` (A:281) | `mse`, cosine floor | `rc_hybrid` 0.8738050 predicted, 0.8738045 board (R197:9) | one-sided; invert only rows with `score_mse` > 0 (A:305) |
| exact `nmae` identity (216/188) | `nmae` | — | mis-signed both board checks: −0.000152 → +0.001762 (R216:555); −0.000108 → +0.00588 scaled (BR222:49). Tested by the 237 slot (§6); 243 does not price it (R243:292) |
| cloud `cell-eval2` (221) | all six on h1-test | 9/10 bit-exact, 10th 3.3e-16 (R221:8) | $1.48 per score (R221:42) |
| Orion thinned / hct116 | relative direction only | — | 26.4%/44.1% of truth energy in a blob vs 3.6% (A:548); hct116 depth-confounded (A:592-598) |
| panel-total digest 900/900 | arm class | `rc_hybrid` (BR185:26-27) | classifies, never prices (A:1448) |

## 5. Closed families
- Global α: 211 `claim` 1.000000/0.999904/1.000000 of cap (R211:19; RB:249).
- Operator γ/λ/exponent/rows: 193 (R193:1); 212 +0.00013/+0.00042 vs +0.0089 (RB:248); 226 +36% rows → +0.57% cosine, 0.49σ (R226:21).
- Ordering: 220 CLOSED (R220:1); 157 (R157:6-7); 164 (R164:6); 232 (R232:1); 234 (R234:1); whole `reach` prize +0.002689 vs +0.0089 floor (R231:19-22).
- Sources: 166 exhausted (R166:1); 205 (R205:1); 218 (R218:1); 225 DLD-1 alone +0.00228 (R225:21); 229 ARM6 −0.00219, ARM4 −0.00589, ARM3 +0.00067 (R229:22-24); 248 DLD-1 rider +0.00216 = 1.14× gate, leave-top-3-out 25.2% (R248:417-421); 204/194/209/203 (R204:1; R194:1; R209:1; R203:9).
- Cross-line/learned: 184 (R184:1); 195 (R195:8); 224 +5.16% cosine (R224:14) killed by `nmae` in 228 (R228:1).
- CD4 dose: 216 `f`=0.34 argmax, turnover 1.018 (R216:1,50); 214 (R214:5-8); 231 Arm A (R231:1). 246 hypotheses 2–5 NOT build-worthy: 3/4 closed, 5 worthless, 2 open-low (R246:23-29,90-92). PerturbMatch −0.003540 but neither channel half clears its bar; ships (R241:13-19,504-508; RB:233-236).
- `nmae` axis: 227 +5.0%/+12.3% cosine vs ±0.2% (R227:13). Transport 213 +1.6% (R213:43). Toll mis-keyed, charge 0; Hamilton +0.002203 inside noise (R230:3-8).

## 6. Open items
- ARMED: pid 26557 (`scripts/arm_submit.sh`) fires `submissions/229-stack6-g035.vcc` `7a3c2d3e…` at 2026-09-05 00:00Z (E237/submit-armed.log:2; P237:82-83): a gate test, 229 measured it −0.00219 (R229:22). Readings, Δ vs 0.2118310: ≤ −0.0015 gate validated; between: veto → discount, re-price 228; ≥ +0.0015 gate wrong, re-open 228 (P237:68-74). 243: +0.0022…+0.0038, differing from 229 by the −0.004783 `nmae` leg (AD:9-16). Raw `mse` bins: ≤0.80768 above all models; –0.81002 dose; –0.81135 coverage; –0.81392 new-source, point 0.81244; ≥0.81392 below all (AD:23-29).
- ARMED second: pid 80453 (`scripts/arm_second_slot.sh`) fires `submissions/249-cd4svd5-f034.vcc` `6e8b41c6…` after the 237 entry is terminal and 00:30Z (E249/submit-armed.log:1). Prereg: `pds` +0.0032/+0.0046/+0.0091 at 0.7/1.0/2.0×, overall +0.0026…+0.0040, board 0.2144–0.2158 (P249:27-40); falsifiers: raw `pds` < +0.0040 retires the class, raw `mse` worse by > +0.0030, Δ outside [−0.0010, +0.0100], amplitude > 0.5% (P249:43-53). Both pids alive 16:58Z (E253/work/arm-procs.log:2-3).
- Running: 252 harness guards #4/5/6/13/14 (E252/MANDATE.md:1); 254 partial rows + composed `cd4svd5+partial` build, no submit (E254/MANDATE.md:1). Pending the 249 point: the October deflation row (RB:224-230; R250:435-448).

## 7. October
- Sequence: 11 stages, one driver, fail-closed (DP:53-54). 244 rehearsal: 29.7 min to a prepped artifact (RB:35-36); build 1302.3 s / 20.00 GiB, prep 22–24 GiB, never concurrent (RB:164,166).
- §0.0 overrides (RB:16-42): CD4 gates fire (union 186 < 290), confirm the mapping with `E244/mapping244.py`, then `--cd4-thresholds relaxed` (RB:19-23); skip-branch hole closed by 247 (RB:24-25); builder `git status` EMPTY before the pin (RB:26-27); six `OK` lines, `def_postbuild` 32/32 (RB:28-31); over-dose 1.170× measured, not delivered under 247's sizing (RB:37-39).
- Guards: 235 `axis_sha256` guard, A/B/C axis `edf6ba88…` (RB:52-63); carrying-row rule 37/37 > 0 where all-rows returns 0.0 (RB:412-414). 240: env are defaults, postbuild diffs against the 222 meta, arms asserted by value (RB:129-135). 247: postbuild operator pins, row count, census gate, novelty ledger (R247:9-16).
- Decision table RB:202-220 + §3.1: PerturbMatch SHIP (RB:233-238); CD4 row 2.07–2.46× understated, `f` 0.34, carrying rows off A/B/C (RB:239-241); deflation row (RB:224-230) and its census trigger (R251:394-395).
- Retention lower bound: random disjoint 300 keeps 93% cosine, 64% `pds`, 60% ordering; 0.2118 → 0.137 / 0.119 (R233:3-11). Essential genes: panel ∩ list 0 on A/B/C, 34/30 on the draws (RB:434); the October-four sources are 70.7–90.4% essential-carried, robust value 16%, `pds` sign flips, 236 retracted (R238:3-10).
- Coverage, pessimistic: `k562_gwps` 269→150, `gse295214_h1` 64→34, `h1_2025` 20→3, `kolf2` 273→0 (RB:326-335); over-dose 1.012–1.017× today, 1.237–1.245× pessimistic (RB:396-397). Day-of known-answer rows RB:268-281.

## 8. Where things are
- Donors `data/donors/<name>/`; union registry MB:840-944 (`gse295214_h1` → `experiments/115-gse295214-h1/work/compact-donor` MB:843; `h1_2025` → `data/donors/h1_2025` MB:852).
- Operators: `data/donors/whitener_union_pruned.npz` `e2f65ecf…` (RB:89); three near-identical operators exist (A:1366-1368).
- CD4 source `experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad` `c355f535…` 16,786,240,107 B (RB:85); banks `data/external/cd4t/cloud-compact/cd4_compact` (RB:279).
- Packages `submissions/`: incumbent `222-h1union-g028.{h5ad,h5ad.meta.json,vcc}` (RB:93-95); armed `229-stack6-g035.vcc` (P237:82-83), `249-cd4svd5-f034.{vcc,h5ad}` (R249:348-349).
- Ledgers: BL; SL 10 entries, five backfilled by 247 (SL:31,46-53; R247:261-275). Driver DP; stages `scripts/def_*.py` (RB:505-507); builder MB; runbook RB. Headline index `experiments/239-knowledge-audit/work/headlines.txt` (R239:13).

Facts checked on 2026-09-04 against the cited files; anything not cited here is not established.
