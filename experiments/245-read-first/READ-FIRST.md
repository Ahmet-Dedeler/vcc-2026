# READ-FIRST — VCC 2026, built from primary sources only (2026-09-04)

Citation key: A=`AGENTS.md` · BL=`docs/BOARD-LEDGER.md` · RB=`docs/DEF-RUNBOOK.md` · DP=`scripts/def_pipeline.py` ·
MB=`scripts/make_submission_pm_cd4.py` · SL=`submissions/SUBMITTED_LEDGER.json` · E<n>=`experiments/<n>-*/` ·
BR<n>=E<n>`/BOARD-RESULT*.md` · R<n>=E<n>`/RESULTS.md` · LB=`experiments/245-read-first/work/ranks.json`
(live board, 2026-09-04 11:49Z) · VS=`experiments/245-read-first/work/vcc-status.log`.

## 1. Board position
- Entry `h1union`, id `iUziqE1gewfcPKzYUDaB`, overall 0.2118310, rank 3/656 at publish (BR222:1-3); live 11:49Z: 0.21183096, rank 3/664 (LB:6-7; VS:2).
- Leader `cqawesome` 0.22431, gap +0.01248 (BR222:21; LB:54,87). 2nd `Aginglab.com` 0.21288 (LB:60-61).
- Members raw / scaled / world rank (LB:8-48; VS:2): `pds` 0.83572/+0.74537/8 · `mse` 0.81726/+0.17957/6 · `nmae` 0.93389/+0.11074/56 · `fid` 0.50748/−0.01648/106 · `reach` 0.27655/+0.22294/7 · `jac` 0.04096/+0.02884/1.
- The four 2026-09-03/04 arms, one variable each (a15 baseline 0.1740865, BL:28):

| arm | entry | overall | Δ | variable |
|---|---|---:|---:|---|
| `rc_hybrid` | `nQUJZhZb3MW2uFJHBX5e` | 0.1834190 | +0.0093325 | emitted p-order, pseudobulk bit-identical (BL:92; BR185:4-5) |
| `rc_callset` | `31oB9yWy8CSNBDEkY4iW` | 0.1918638 | +0.0084448 | call set top-400 `|H|` ∪ top-2,200 CPM, `n_pred` 4,480→2,599 (BL:93; BR202:4-5) |
| `cd4f034` | `QSe3JlNUYZ2BzvlcsL14` | 0.2049196 | +0.0130557 | `VCC_CD4_TARGET_F` 0.10→0.34 (BR210:5) |
| `h1union` | `iUziqE1gewfcPKzYUDaB` | 0.2118310 | +0.0069114 | `VCC_H1_UNION_G`=0.28 (BR222:6) |

- `cd4f034` overall: 0.2049216 (BL:203; RB:200) vs 0.2049195830 (BR210:5).
- Deterministic scorer: a byte-identical resubmission reproduced every member (BL:8-9); `pds`/`mse` identical to 16 s.f. across a15/`rc_hybrid`/`rc_callset` (BL:100-101).
- `vcc status` on 2026-09-04 resolves only `iUziqE1gewfcPKzYUDaB`; six earlier ids return `not_found` (VS:3-13).

## 2. The shipped configuration
- Env: `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28` (DP:84; RB:89-91); `VCC_UNION_SOURCES` unset = `gse295214_h1, h1_2025` (DP:85; RB:183).
- Flags: `--min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0` (DP:92-93; RB:135); `--emitter rc_callset` (RB:135,179); call set head 400 / cpm 2,200 / z_cap 1.9 (RB:180); γ 0.65, λ 0.0 (RB:176-177); donor weights 1:1:1 (RB:187).
- Operator by sha: A/B/C `data/donors/whitener_union_pruned.npz` `e2f65ecf…` (RB:60); D/E/F must be re-fit as `union_pruned_def` (RB:131).
- The two arms absent from the 09-03 recipe (191/202): `cd4f034` +0.0130557 and `h1union` +0.0069114, sum +0.0199671 (RB:93-96). The pre-240 driver dropped both, its post-build check passing 14/14 (R240:22-24).
- Package: `submissions/222-h1union-g028.h5ad` `8df4dadd…` 12,648,807,800 B; `.vcc` `f13ce40b…` 2,824,785,920 B (RB:64-66); the fixed driver's defaults rebuild it byte-for-byte (R240:15-21).

## 3. Rules that bind every agent
- 2 submissions/day/team, reset midnight UTC, one in flight (A:40). A slot buys information, never rank; never re-submit a scored package (A:1319-1322).
- §14: `-m` is a random joke; `-d` EMPTY; never confirm/deny a guess; applies to every public surface (A:1731,1751,1753,1756,1761).
- One variable per run (A:1315). Preregistration committed before the build (BR185:7-8; BR202:5; BR210:7).
- Price all six members or state why one cannot move (A:1459). Gate: > +0.0019 shape-only, > +0.0089 amplitude-changing (A:1446); no `reach` toll (A:1400); ordering noise between digest-identical arms sd 0.004502 overall (A:1406).
- Every screen carries a known-answer row and aborts if it does not reproduce (A:1347); diff the built config against the incumbent and require exactly the intended differences (A:1332).
- Record path + sha, never a description (A:1469). `fid`/`jac`/`reach` normalisation unreadable on harnesses (A:1479).
- Never `git add -A` with agents running (A:1513). Memory gate: `memory_pressure` ≥ 40% free AND free+inactive ≥ 20 GiB AND swap free ≥ 400 MiB, heavy jobs serial (A:1524-1528); one `cell-eval2` at a time, ~25 GB RSS (A:592-593).
- `uv run python` (RB:72-74). Launch via `scripts/claude_subagent.sh <exp-dir> <mandate>` (A:1574).

## 4. Instruments
| instrument | reads | validated on | known failures |
|---|---|---|---|
| board (`vcc status`) | all six | byte-identical resub reproduced every member (BL:8-9) | only the current id resolves via `vcc status` (VS:3-13) |
| h1-test `data/harness/h1-test.h5ad` (A:590) | DE members, relative | a8/a15 `reach` pair, two board points within 0.014 raw (A:543,547-549) | oracle-leaked when `data/donors/h1_2025` is a donor: 85/85 targets (A:573-574); shipped stack carries `h1_2025` (DP:85; MB:699) so every arm on top of `h1union` is leaked. Not validated on emission-law arms: 2/4 DE signs wrong on a21 (A:547-551). `reach` order-arm factor 1.290 (BR185:50) vs 0.769 on the real scorer (R221:23,444); `mse` sign-inverted on the direction arm, `fid` on the call-set arm (R221:20,22) |
| PLOLO `E162/plolo.py` | `pds` of geometry arms | 6/6 signs, LOO RMS 0.001494 scaled (R162:9-10) | blind to the other five, unpriced/priced 7.4× (R162:17). PLOLO→board cosine factor not identified: 1.384 (`cd4f034`) vs 1.000 (`h1union`) (BR222:61-64). Mirror returns `f_per_unit` 0.0 below 50% coverage (R223:20; RB:349-351) |
| 197 inversion `cos ≥ (1.023611−raw mse)/2.074495` (A:272) | `mse`, implied cosine | `rc_hybrid` raw `mse` predicted 0.8738050, board 0.8738045 (R197:9) | one-sided; only invert rows with `score_mse` > 0 (A:295-296) |
| exact `nmae` identity (216/188) | `nmae` | — | mis-signed both board checks: −0.000152 predicted, +0.001762 realised (R216:555); −0.000108 predicted, scaled `nmae` +0.00588 realised (R222:27; BR222:49) |
| cloud `cell-eval2` (221) | all six on h1-test | 9/10 statistics bit-exact, 10th to 3.3e-16 (R221:8) | $1.48 per six-member score (R221:42) |
| Orion thinned / hct116 | relative direction only | — | 26.4%/44.1% of truth energy in a target-independent blob vs h1-test 3.6% (A:539); hct116 depth-confounded: 093's +0.0202 became −0.00291 (A:585-589; BR093:4,16) |
| panel-total digest, 900/900 blocks | arm class | `rc_hybrid` (BR185:27) | classifies, does not price (A:1420) |

## 5. Closed families (closing experiment · headline)
- Global α: 211 · `claim` 1.000000/0.999904/1.000000 of cap (R211:19); α* 0.9808/0.9879 (RB:198). Per-context knee: abc-knee `mse` −0.00226 (BR180:32).
- Operator γ/λ/exponent/rows: 193 (R193:1) · 212 best +0.00013/+0.00042 vs +0.0089 floor (RB:197) · 226 +36% rows → +0.57% macro cosine, 0.49σ (R226:21).
- Ordering keys: 220 CLOSED (R220:1); 157 NO-GO (R157:6-9); 164 call-set arms NO-GO (R164:6-11); 232 sign vote NO-GO (R232:1); 234 flip keys NO-GO (R234:1). Whole `reach` prize +0.002689 vs +0.0089 floor (R231:19-22).
- Donors/sources: 166 universe exhausted (R166:1); 205 artifacts do not exist (R205:1); 218 independence does not predict (R218:1); 225 DLD-1 alone +0.00228 (R225:1,21); 229 ARM6 −0.00219, ARM4 −0.00589, ARM3 +0.00067, NOT RECOMMENDED (R229:17-27); weighting/re-estimation/Orion subspace (R204:1; R194:1; R209:1); 203 estimator −0.001393 vs +0.0089 (R203:9-10).
- Cross-line/learned: 184 (R184:1); 195 (R195:8); 224 per-gene diagonal +5.16% cosine (R224:14) killed by `nmae` in 228 (R228:1).
- CD4 dose: 216 `f`=0.34 argmax, turnover `f*`=1.018 (R216:1,50); 214 PerturbMatch is a complement (R214:5-8); 231 Arm A refuted (R231:1).
- `nmae` as a separate axis: 227 · costs +5.0%/+12.3% cosine vs a ±0.2% gate (R227:13). Transport law 213 · max +1.6% (R213:43). Mass fix a23 −0.00073 (RB:184). Resample emitter a21 `reach` −0.0185 (BL:30). `reach` toll: 230 mis-keyed, charge 0 (R230:3-8). Hamilton rider +0.002203 vs +0.0019, inside noise, not cleared (R230:8; A:1410-1411).

## 6. Open items
- ARMED: `scripts/arm_submit.sh` pid 6813 will `vcc submit submissions/229-stack6-g035.vcc` (`7a3c2d3e…`) at 2026-09-05 00:00Z, name `Reviewer 2 Was Right About Everything` (`E237/submit-armed.log:1`). Purpose: validate the exact `nmae` gate; readings Δ ≤ −0.0015 / between / ≥ +0.0015 (`E237/PREREGISTRATION.md` §4-5). 229 measured that package at −0.00219 and marked it NOT RECOMMENDED (R229:22,27).
- Running (`claude-run.pid` live, 11:50Z): 241 ablation on the shipped renderer (`E241/MANDATE.md:1`); 243 PLOLO→board cosine factor 1.384 vs 1.000 (`E243/MANDATE.md:1`); 246 direction-channel hypotheses (`E246/MANDATE.md:1`).
- Mandated, not started: 244 end-to-end Oct-22 rehearsal on the fixed driver, essential-free pessimistic panel (`E244/MANDATE.md:1`); RB requires the re-run by 2026-10-15 (RB:332-333).

## 7. October
- Sequence: 11 stages, preflight→…→prep, one driver, fail-closed (RB:85; DP:53-54); 1,748 s ≈ 29 min to a prepped artifact (RB:149); build 1302.3 s / 20.00 GiB, `vcc prep` 22–24 GiB, never concurrent (RB:135,137).
- Guards: 235 — operator carries `axis_sha256`, builder recomputes it from the panel and raises `REFUSING TO BUILD`; A/B/C axis `edf6ba88…`; `n_cols` is 17,942 on all four panels so it cannot discriminate (RB:131; R235:19); carrying-row rule sizes 37/37 (source,panel) cells > 0 where the all-rows rule returns 0.0 in all 37 (R235:21-24). 240 — env are defaults, postbuild diffs against the 222 meta, arms asserted by value (RB:136,202; R240:15-24). Trap priced −0.0429 overall on `pds` alone (RB:204). 242 — 21 guards can pass while wrong, 5 on the Oct-22 path: postbuild operator block passes on `data/donors/whitener_union.npz`; `solve_cd4_fraction` over-doses 1.4222× at 58.7% coverage; novelty check passes a scored package (R242:16,43-45).
- Retention lower bound (233): random disjoint 300 keeps 93% cosine, 64% `pds`, 60% ordering; overall 0.2118→0.137 (65%), 0.119 (56%) at pessimistic coverage; lower bound, not forecast (R233:3-11).
- Essential genes: panel ∩ essential list = 0 on A/B/C, 34 realistic draw, 30 pessimistic (RB:383); 238: the October-four sources' gain is 70.7–90.4% essential-carried; robust value 16% on an essential-free panel, `pds` +0.001904→−0.000978; 236's +0.0066…+0.0104 retracted (R238:3-10; R236:3-11).
- Coverage, pessimistic draw: `k562_gwps` 269→150, `gse295214_h1` 64→34, `h1_2025` 20→3, `kolf2` 273→0 (RB:275-284); `panel_support` 248,276→141,167; CD4 union 291→174 (RB:381-382). Over-dose 1.012–1.017× at ~92% coverage, 1.237–1.245× pessimistic (RB:345); 0 of 21 experiments void (R223:1). Day-of known-answer rows: RB:215-231.

## 8. Where things are
- Donors `data/donors/<name>/`; union-source registry MB:687-793 (`gse295214_h1`→`experiments/115-gse295214-h1/work/compact-donor` MB:690; `h1_2025`→`data/donors/h1_2025` MB:699).
- Operators: `data/donors/whitener_union_pruned.npz` `e2f65ecf…` (RB:60); three near-identical files, max eigenvalues 8.036/6.257/5.895 (A:1338-1342); fitter `scripts/build_union_whitener.py` (RB:131).
- CD4: source `experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad` `c355f535…` 16,786,240,107 B (RB:56); banks `data/external/cd4t/cloud-compact/cd4_compact` (RB:228).
- Packages `submissions/`: incumbent `222-h1union-g028.{h5ad,h5ad.meta.json,vcc}` (RB:64-66); previous `210-cd4f034.h5ad` (BL:260).
- Ledgers: BL (every scored arm, six members); SL holds 5 entries, last `IzQVfDpeDWufkPT8Lmd2` 2026-09-02 (SL:12-13) — the five later entries in BL:92-94, BR180:3, BR222:3 are absent (R242:45).
- Driver DP, post-build `scripts/def_postbuild.py`, builder MB, runbook RB (RB:85-137). Harness PLOLO `experiments/162-retrodicting-validator/plolo.py` (RB:351). Headline index `experiments/239-knowledge-audit/work/headlines.txt` (R239:13).

Facts checked on 2026-09-04 against the cited files; anything not cited here is not established.
