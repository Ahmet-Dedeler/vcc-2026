# VERDICT-2: 813 PROPOSAL-2 (TeloHAEC genome-wide CRISPRi as a fourth held-out BULK fold line)

Judge run, 2026-10-16 (clock: `date -u` 2026-10-16T09:27:42Z during the review). Judged the FROZEN copy
`experiments/819-judge/813-public-data-20261016T092345Z/PROPOSAL-frozen.md`, sha256
`a7628506f4dc7b9a6e54aa833dfd75ef1c4bf940a48e30da9366678e65100166` (re-hashed here; it matches the mandate's prefix and the
live `813-public-data/PROPOSAL-2.md` byte for byte). This is the one resubmission after `813-public-data/VERDICT-1.md` (REVISE,
R1–R7, C1–C4). The lab lead's calls are in `800-checker/resume/813-verdict1.md`.

Lookups only. The checks I ran:
- sha256 of the three frozen artifacts: all three match frozen l.9–11;
- `patch --dry-run -p0` of `patches/lib801-telohaec.patch` against `lib801.py` (sha `41035aeb…`, as frozen l.10 says): clean,
  rc 0, nothing written;
- sha256 of `data/donors/whitener_union_pruned.npz`: `e2f65ecf…`, as frozen l.136 says;
- sha256 of the terms text `work/meta/igvf_external_resource_sharing_policy.txt`: `fc17237b…`, and the operative sentence is
  present once.

## Verdict: **GO**, for stage 0 only ($5), with conditions G1–G6 below

Every item of VERDICT-1 is fixed (table below). The idea is not closed. It prices no arm and claims no gain, and its kill
gates are numeric and preregistered. KA-SG is now a noise-unbiased statistic whose bounds come from simulation before any
real row exists. P2 can now fail and is named honestly.

Stage 1 ($6) is not granted by this verdict. It is released by the lab lead after reading the stage-0 RESULTS (C1). The
remaining gaps are small and are made binding conditions here, not a third round:
- an undefined middle zone between P2's two bins;
- a memory estimate that contradicts the fit's own guard flag;
- the unnamed dump invocation for the L folds.

### VERDICT-1 items, one by one
| item | fixed? | where (frozen file) and what settles it |
|---|---|---|
| R1 KA-SG | **yes** | l.70–87. The cell sets are disjoint (n_guides ≥ 2 vs == 1, l.71–72). β = Σ⟨A,S⟩/Σ⟨S^h1,S^h2⟩ and dsh are unattenuated by independent noise (`kasg813.py` l.12–20, sha `b94bf49e…` verified). [lo, hi] is the central 95 % over R = 200 unbiased KA-SIM replicates, set in stage 0a before any real row (l.77–78). A power row turns an underpowered statistic into "uninformative, not a gate" instead of a false STOP (l.79–80). The μ̄ blind spot is stated (l.81–82). The statistic is frozen by sha, committed before 0b, and 0b refuses on a sha mismatch (l.83–87). |
| R2 P2 can fail | **yes** | l.156–177. 702 is dropped with the reason and the citation (l.157–163; `700-pie-transfer/RESULTS.md` l.49 checked: tahoe, jiang, arc_vcc_25, orion; Replogle weight 0). The claim "only 479sp2 is expressible" checks out against `engine801.py` `named()` l.230–235 (cap/pct/p96.6 × P2g/S_P2 only; no 417a, 471c, 458d or DE rows). 479sp2 is reported, not graded, which is correct: every lab surface misses its class (LAB-LOG l.130, l.206–212). The bins are renamed, with no "TRUSTED". SIGN-CONSISTENT is labelled necessary-only, INCONSISTENT rejects the instrument, and the folds stay "uncalibrated for negatives" until a discriminating board row is graded (l.169–177). |
| R3 fold build | **yes** (G3, G4) | l.122–151. Six new full 3-donor day-recipe fits, with the cost basis in `cloud/runs.jsonl` l.92 / l.103 (r7i.8xlarge, $0.697) and `801-lab-folds/RESULTS.md` l.363–365 (322–338 s). The sha recording and the reproducibility rule rest on 801 RESULTS l.30–32. The ABC operator sha `e2f65ecf…` is verified, and the route KA is ABC/kolf2 (801 RESULTS l.326, 8.7e-15). The patch is UNAPPLIED, adds only new keys, and dry-runs clean. Its KA (`ka_patch813.py`, sha `b53d34a3…`) requires every pre-existing registry entry, dump sha and key to be byte-identical, plus a row re-dump of L1/k562_gwps and abc/kolf2. |
| R4 P3 | **yes** | l.178–179. Dropped as a blob sentinel (one of the lead's two allowed choices). |
| R5 role + leakage | **yes** | l.181–199. Reported apart, never pooled into BULK-W, never a veto alone. A cancer-fold gain with a TeloHAEC upper bound < 0 carries a flag with its mechanism (the lead's call, verbatim). Every TeloHAEC row is dropped on a TeloHAEC-held fold. The check800 refusal is fail-closed on a missing `banks` list, and a declared-seen fold is "cost-bound only", which is DESIGN §7 l.93–94 (702 rule 2). Wear goes to `surface=telohaec`. LZ is excluded (C2). |
| R6 novelty / 105 | **yes** | l.29–34, l.202–209. 105 is cited (`105-telohaec-residual/RESULTS.md` l.15–17: −0.141 / −0.127 checked), and P2 is named as the test of whether a fold that far away still orders arms. Use (b) is deleted (l.22–23) and stays 811's pre-counted try (`811-lolo-model/PROPOSAL-2.md` l.24–26, l.156–157 checked). |
| R7 terms | **yes** | l.107. `813/RESULTS.md` §6: the text's sha `fc17237b…` is verified, the operative sentence ("External users may freely download, analyze and publish results … without restrictions") is present, and the citation duty is recorded. It is re-checked as the first step of 0b, with STOP and no byte fetched if the sentence has gone. |
| C1–C4 | **yes** | l.221–229, verbatim from VERDICT-1 l.186–193. |

## The JUDGE.md checks

### 1. Is the mechanism closed already? **No.**
- `grep -ril telohaec experiments/*/RESULTS.md` gives 22 files. Every earlier use is as a source, donor or prediction target
  on the 2,167-target pilot:
  - 105: direct source;
  - 218 l.161 / l.303: VACUOUS / zeroed as a donor;
  - 223 / 225: `telohaec@0.28` −0.00120;
  - 122 / 130 / 135: the atlas, terminal NO-GO;
  - 224 l.66: "105 + 223 NO-GO";
  - 252 l.39: silent zeros as a source;
  - 811 l.294: pilot ablation ≤ −0.0007 as a voting line.
- None uses TeloHAEC as a held-out fold truth. FINDINGS §2's "Direction via a fifth donor line: sources saturate"
  (`docs/FINDINGS.md` l.885) closes adding a DONOR. This lane adds none (frozen l.19–20).
- The proposal is not a re-named closed construction. 105's distance is carried as the risk P2 must survive (l.29–34), not
  waved away.

### 2. Does the prediction have a falsifier and bins written before any number? **Yes, with one gap (G2).**
- **Stage-0 gates are numeric and preregistered:**
  - G0-KD: ≤ 0.60 median, ≥ 60 % < 0.70;
  - G0-POWER: ≥ +0.028, ≤ 34 % at or below zero (271's definition);
  - G0-BLOB: ≤ 0.10 (615's statistic);
  - KA-SIM: three pass rules plus the A/B choice margin;
  - KA-SG: bounds from simulation.
  All are at frozen l.116–120 and l.64–80, with P1's expected values and the likeliest failure at l.154–155.
- **No re-tune after a failure:** frozen l.87 and C3.
- **P2's falsifier is real.** INCONSISTENT means an upper bound < 0 on abc/telohaec_gw, or ≤ 3/6 L folds positive on ≥ 2
  rows (l.172–173). It rejects the line as an instrument. The three graded rows are near-certain positives on CANCER folds
  (LAB-LOG l.127–129), but not on a line 105 puts 0.13–0.14 pds away from every donor. So this is the honest test the
  proposal says it is, and it is only a necessary condition.
- **Gap: the two bins are not exhaustive.** For example, all three rows point-positive with abc/telohaec's lower bound ≤ 0
  and no upper bound < 0, or one row at exactly 3/6 L folds, lands in neither bin. → G2.

### 3. Is it measurable on the checker's surfaces? **Yes, for what it claims.**
- **Bulk members only:** literal pds and mse through r, as on every L fold (DESIGN §4 l.65–72). It is an instrument: it
  prices no arm, so it needs no DE member. DE stays with DE-H1 / DE-CD4.
- **Gene axis:** the 5-CPM gate follows the scorer's `filter_gene_min_cpm_cell` (BRIEF l.83 checked: "5 CPM, control cells
  only") on the line's own background cells. It lives in the bank, so no engine change is needed (l.93–98).
- **L-fold wiring.** It needs the patch's `FOLD_SOURCES["telohaec_gw"] = MAIN` and the `plolo.DONOR_PATHS` entry.
  - The L panels are not in `PANELS`. dump801 binds them from a CLI-built spec (`dump801.py` l.194–198).
  - The L-fold operators therefore enter through `--fold-ops prod --op <operator>`, where `fold_operator` reads
    `ctx.full_op` with MAIN as sources (`lib801.py` l.340–342).
  - The proposal does not name this invocation. → G4.
  - `fold_operator`'s docstring calls 'prod' "leaky: the held library is in its basis" (lib801 l.335–336). That does not
    hold here: TeloHAEC is in no basis.
- **Memory.** The proposal estimates "≤ 18 GiB each" (l.130) but runs `--max-rss-gib 14` (l.125–126), and that guard RAISES
  above the limit (`scripts/build_union_whitener.py` l.105–109).
  - The record says the estimate is pessimistic. Production's own 3-donor fit of the same recipe (29,394 basis rows, as
    `data/donors/whitener_union_pruned.json` records) took 320–437 s with a 12.9–13.7 GB high-water under the same 14 GiB
    guard (`docs/DEF-RUNBOOK.md` l.386).
  - So the fits fit, but on a thin margin. → G3.

### 4. Does it use any seen data as evidence of gain? **No.**
- No gain is claimed. The board reads of 495 / 479p2g / 498p2g enter only as known-answer calibration (DESIGN §7
  l.86–94).
- KA-SIM uses K562-GWPS rows only as a simulated truth for the estimator (l.57–58). That is not evidence about any arm.
- Neither the engine's named arms nor production's stack carries TeloHAEC rows (`813/RESULTS.md` l.61: "Not in any production
  donor mix").
- The leakage rule covers every future candidate (l.187–198).
- No sealed reserve is touched: only L1–L6 and ABC; LZ is excluded by C2 (l.199). R6, R7a and R7b are not named.

### 5. Is the try budget honest? **Yes.**
- 0 arms priced on BULK-W by this lane (l.216).
- The estimator (2 options) is chosen on simulation only, and KA-SG's statistic and bounds are frozen before any real row
  (l.213–214).
- Control, gate and target set: one each, with a seeded, sorted list for stage 0 (l.109, l.215).
- P2's rows are the engine's existing named known answers, so they choose nothing.

### 6. Is the expected gain worth its cost? **Yes.**
- It is an instrument, so there is no board gain to compare with the class noise.
- **What it buys:** a fold line from outside the one-lab Orion pair, the strongest mutual partnership in 811's ablation
  (`811-lolo-model/RESULTS.md` l.291 checked). It covers every lab panel at near-full coverage, and a same-day fold on
  D/E/F's real panel computed from public data only (l.217–219).
- **Its ceiling, stated honestly:** reported apart, never a veto, so it can at most flag (l.183–186). The proposal's own
  P(stage 0 passes) is ≈ 0.6 (l.155).
- **Price:** the kill gate is ≤ $5 and the lane ≤ $11 (l.113), against the lab's $800. The instance prices check out:
  - r7i.xlarge 45 min ≈ $0.20;
  - r7i.2xlarge 120 min ≈ $1.06;
  - r7i.8xlarge 75 min ≈ $2.65.
- **Download:** IGVF serves presigned us-west-2 S3 URLs (`813/RESULTS.md` l.61, l.100). The egress is not ours, and nothing
  lands on the Mac.

## Conditions (binding on the GO; C1–C4 from VERDICT-1 carry unchanged)
- **G1: stage 0 only.** This verdict grants stage 0 (0a, then 0b; ≤ $5, the retry included). Stage 1 (≤ $6) is launched
  only after every stage-0 gate passes and the lab lead has read the stage-0 RESULTS (C1).
  - 0b starts only after the lab lead has committed `freeze813.json`'s sha (C3).
  - Run `labcloud.sh status` after each launch. No `Project=lab` instance may be left running.
- **G2: P2's third bin, fixed now.**
  - Any P2 outcome that is neither SIGN-CONSISTENT nor INCONSISTENT is **UNRESOLVED**. The folds are then reported apart
    under "uncalibrated" and may NOT carry the §6 dossier flag until a later graded board row resolves them.
  - Before stage 1b's P2 rows are read, 813 writes into RESULTS.md its expected bin and its expected sign per row per
    surface (abc/telohaec_gw; the L-fold count).
- **G3: the operator-fit guard.** `--max-rss-gib` only raises (`build_union_whitener.py` l.105–109); it does not enter the
  operator's numbers.
  - A guard trip on a first fit is an infrastructure failure, not a refit. Rerun it on the same instance type and thread
    count with a higher guard, and record it.
  - The "new fold id on refit" rule (l.133–134) applies only once an operator has been stored.
  - Record each fit's wall time and peak RSS.
- **G4: L-fold wiring, explicit.** Stage 1b's job list names, per L panel, the exact dump invocation: dump801 with fold
  `telohaec_gw`, `--fold-ops prod` and `--op` = that panel's stored `full_telohaec_gw.npz`.
  - Each dump asserts the sha of the operator it read against the registry entry.
  - The bank is staged only on S3, at the mirror path the patch reads (`input/lab/repo/experiments/813-public-data/work/
    bank/telohaec_gw/`), with a sha list. It never goes to the Mac.
  - The engine's own fold KAs and `ka_patch813.py` `check` / `rows` must PASS before any P2 row is read. A KA counts as
    passed only on rc 0 AND the matching hash or tolerance (AGENTS 6.2).
- **G5: the KA-SG fallback is stated, not silent.** If KA-SG is recorded as uninformative at this depth, the stage-0
  RESULTS says so, along with the consequence: on real data only G0-BLOB then guards against k-dependent offsets. The
  folds carry that caveat in the registry entry.
- **G6: data and cost hygiene.**
  - Only the resting reps (stimulated conditions dropped, l.114).
  - Nothing > 200 MB to the Mac, and no Mac transfer in a WINDOW (C4). The Mac is under the coordinator's local-job HOLD
    (LAB-LOG l.302–307), which binds local jobs only; EC2 is unaffected.
  - Every launch and its cost go into `docs/LAB-LOG.md`'s cost table through the lab lead. `813/RESULTS.md` cites the IGVF
    accessions and the Engreitz lab, as §6 of that file records the terms require.
