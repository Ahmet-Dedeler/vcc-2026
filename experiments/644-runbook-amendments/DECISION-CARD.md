# Oct 22 day-of decision card (644, draft 2026-10-02; re-verify after the 10-05 / 10-06 patches and on Oct 21)

Read top to bottom on the day. Each row says what you see, what it means, and the exact action.

- **STOP** means: do not ship. Fix the cause, or ask the coordinator.
- **Every rerun** goes in a FRESH work dir. A tripped guard stops `def_day.sh`, and nothing resumes in place.
- **Every override** is exported before `def_day.sh`, and its reason is written into the round's RESULTS.md before the build.

## 0. Before the build

| you see | meaning | action |
|---|---|---|
| `october_recipe.py show` lists a switch as PENDING | its read has not landed | `command` refuses; ship the DECIDED recipe; never hand-type a PENDING value |
| `vcc datasets list` shows no D/E/F id | the release is not up yet | wait; keep A/B/C in `data/challenge-abc` |
| downloads | standing OK from Ahmet (2026-10-02): never ask | download the D/E/F bundle and the `vcc` upgrade; record file, source, size, sha256 |

## 1. Preflight (the driver's first stage)

| you see | meaning | action |
|---|---|---|
| `[X] control library q75/q25 <= 3.0 (554)` FAIL | wide libraries: the emitter's steering dose would collapse silently above ~3.5 | ≤ 3.5 in every context: export `DEF_DAY_ALLOW_WIDE_LIBS=1` and rerun; §6's 554 row must then PASS. Above 3.5 anywhere: **STOP** |
| `[X] control library median within 15% of A/B/C's 20,030 (641)` FAIL | D/E/F released at another depth | all three off by the same factor: price the ratio on h1r7 (553), then export `DEF_DAY_ALLOW_DEPTH=1`. One context alone: **STOP** |
| `obs has 'ntc_id'` FAIL | the release dropped the NTC column | `DEF_DAY_ALLOW_NO_NTC_ID=1` only for a genuine format change |
| `cells_per_pert == 400` FAIL | a declared format change | `DEF_DAY_ALLOW_CELLS_PER_PERT=<n>` only if the manifest really says so |
| gene axis reordered vs A/B/C | the pinned CD4 deflation basis raises | **STOP** (531 F6): do not refit on the day |
| the bundle sha equals A/B/C's | `data/challenge` still holds A/B/C | check `manifest.json` |

## 2. Identity (590, and 642 if needed)

| you see | meaning | action |
|---|---|---|
| any context | — | copy into RESULTS: call, n/8, `min_top1_score`, consensus top 5 with lineages |
| registered, 8/8, min_top1 ≥ 0.70 | a clean same-line match | proceed (same-line rows switch on, w from the bank class) |
| registered, 8/8, min_top1 0.60–0.70 | borderline | read the top 5 by hand. Sister or duplicate lines, or mixed lineages: rerun in a fresh dir with `DEF_DAY_SAMELINE_SKIP=<ctx>` (644; that context takes the builder's unmatched path, the others are untouched). Never `--force` without a written reason |
| `… < 0.65; context unchanged (590)` | a line probably absent from CCLE | expected; read the top 5 |
| …and the top 5 looks kidney / neural / embryonic | possibly HEK293 | run the 642 diagnostic: `DEF_DAY_ALLOW_NONCCLE_REF=hek293 uv run python scripts/def_identity.py --bundle data/challenge --contexts D,E,F --out <work>/identity-hek293.json` |
| `[642] X: … -> CALLED HEK293` (8/8, ≥ 0.75) | the HEK293 reference row beats every CCLE line | coordinator decides in writing. If yes: rerun the day with `DEF_DAY_ALLOW_NONCCLE_REF=hek293 DEF_DAY_NONCCLE_DECISION="<who, when, context, why>"` (orion_hek293t rows, w 2). Orion HEK293T rows remain unpriced on an Arc-Flex HEK293 context |
| `[642] … not called` | stage 2 did not fire | the context stays without same-line rows |
| `642: … not the pinned …` / `… lists Flex …` | the reference file changed | `git checkout -- data/identity_refs/hek293_nonflex_ref642.npz` |

## 3. Setup and build

| you see | meaning | action |
|---|---|---|
| `CD4 extraction gate fired` | Arc-style panels trip 047's gates (481: P 0.29–1.0) | check the mapping (244's mapping244.py), then rerun with `DEF_DAY_CD4_THRESHOLDS=relaxed` |
| the same-line stage dies with `DEF_DAY_PROGRAMME_SOURCE` set | too few rows or low coverage | rerun at once, source unset, flags `noprog` (R0 must not wait) |
| `rc_hybrid_k did not meet the per-cell totals` | same-line w too high for the library capacity | rerun with w halved, then without `:pt` |
| the re-level does not converge, or ρ > 1 | the mse cliff | never ship ρ > 1; `def_day.sh` falls back to the nearest safe artifact |
| `builder on disk … is not HEAD` | a reboot mid-pin | `uv run python scripts/def_pin_panel_artifacts.py --restore --record <previous-work>/pin-artifacts.json` |

## 4. Postbuild and audits (before the builder restore)

| you see | meaning | action |
|---|---|---|
| `postbuild PASS` | config diff and dose guard clean | go on |
| `(554)` row: prior < 0.9, or a dose near 0 | a context lost its DE design | **STOP**: find out why (programme source? library capacity?) |
| `(554)` row: only an MDE median 0.85–0.90 | ~10–15 % MDE loss, ≈ −0.0003 overall | rerun the day in a fresh dir with `DEF_DAY_ALLOW_LOW_DOSE=1` (~2 h); record postbuild.json's `dose_guard` |
| config diff outside the allowed keys | the build is not the recipe | **STOP**: rule 10 |
| zero-gene audit not clean | a gated gene with a zero block total | **STOP**: the repair floor must hold |
| `vcc prep --dry-run` fails | format contract | fix the contract; never upload around it |

## 5. Upload and after

| you see | meaning | action |
|---|---|---|
| the upload byte range has not moved for 30 min | a stall | restart the SAME bytes; a slow but moving range is left alone |
| slot 1 not PUBLISHED by 16:00 UTC | uploads take hours | do not start slot 2 that day |
| the first D/E/F entry: `capture_board_entry.py` writes `members_hidden` | the dark board | the dark plan stands; the last entry = the best-supported scored package |
| …members come back | the final's own truth is readable | follow FINAL-PHASE-FEEDBACK-PLAN. Feedback-use policy (H-631): whatever §6 item 8 says once Ahmet confirms (two readings recorded there). Factorial reads (637: one bulk row 2/3/4 + one prior-only DE row 6/7) are the efficiency default once its A/B/C validation reads PASS |
| Nov 4 16:00 UTC | the final window | upload the chosen scored package under the reserved final name; start by 20:00; nothing new after Nov 5 20:00 |

Overrides are logged at start and recorded in `<work>/def_day.final.env`:
`DEF_DAY_ALLOW_WIDE_LIBS`, `DEF_DAY_ALLOW_LOW_DOSE`, `DEF_DAY_ALLOW_DEPTH`, `DEF_DAY_ALLOW_NONCCLE_REF`,
`DEF_DAY_NONCCLE_DECISION`, `DEF_DAY_SAMELINE_SKIP`.

Older overrides are covered in DEF-RUNBOOK's top block: `DEF_DAY_ALLOW_NO_NTC_ID`, `DEF_DAY_ALLOW_CELLS_PER_PERT`,
`DEF_DAY_CD4_THRESHOLDS`, `DEF_DAY_KEEP_PINS`.
