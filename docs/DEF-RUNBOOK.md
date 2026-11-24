# D/E/F RUNBOOK — Oct 22 (panel lands, board goes dark) → Nov 5 23:59 UTC (last entry counts)

Written by experiment 199 on 2026-09-03 after rehearsing the whole path end to end on two synthetic
D/E/F bundles. Every command below was run; every number is measured on this machine (M-series
Mac, 64 GB, 18 cores). Follow it top to bottom. Do not improvise on the day — if something here
is wrong, fix the runbook first, then act.

**Absolute rules for the window.** Only the LAST final-round entry counts. The board is dark, so
nothing you upload can be read back. `-m` is a random joke, `-d` is EMPTY (AGENTS.md §14). Never
run `cell-eval2` alongside a build. Never re-upload a package you have already uploaded.

---

## THE DAY IN ONE PAGE (2026-09-25; supersedes §2's command block where they differ)

> **2026-10-02 (644: the guards and switches added since 09-30; these override anything older on this page).**
> - **Patch order:** 538 → 554 → 590 on 10-05, 641 at checklist step 6b, then 598 (tercilec) → 612 → 612dd → 630 → 642 → 604
>   → 644 (SAMELINE_SKIP) on 10-06, in 606's tested chain (`experiments/606-applyday-composition/`; NEXT-SLOT-DECISION's 10-06
>   APPLY DAY lists the expected post-apply shas).
> - **Every new guard is silent on PASS.** An override is legitimate only in the case its row names.
> - **To override:** export it in the shell before `def_day.sh`, in a FRESH work dir (a tripped guard stops the day and
>   nothing resumes in place), and write the reason into the round's RESULTS.md before the build starts. `def_day.sh`
>   logs each override at start and records it in `<work>/def_day.final.env`.
>
> | guard (source) | the row that trips | meaning | override, and when it is legitimate | record |
> |---|---|---|---|---|
> | preflight library spread (554) | `[X] control library q75/q25 <= 3.0 (554: …)` | The K-class emitter's library classes are the control quartiles. Above q75/q25 ~3.5 its steering dose (CLEAN4 block prior, MDE, hook) drops to 0 with no error (553); by ~7 the placement crashes. A/B/C sit at 1.57–1.95 | `DEF_DAY_ALLOW_WIDE_LIBS=1`, only if every context reads ≤ 3.5; the 554 dose row below must then PASS. Above 3.5 in any context: STOP (no tested fallback exists for the real day) | the three spreads (preflight.json `library_spread_q75_q25`) |
> | preflight depth (641) | `[X] control library median within 15% of A/B/C's 20,030 (641: …)` | D/E/F were released at another depth than A/B/C's ~20k (Arc's full data is 29,482 UMI/cell). Every per-gene quantity is refitted from the controls, but the DE call-set size was tuned at 20k | All three contexts off by the same factor: price that ratio on H1's regime instrument (553's h1r7 how-to), then `DEF_DAY_ALLOW_DEPTH=1`. One context off alone: STOP (a damaged or mislabelled file) | the three medians and `library_median_ratio`, plus the price |
> | identity floor (590) | `[sameline] X: … min top-1 correlation 0.6xx < 0.65; context unchanged (590)` | A line absent from CCLE is matched, unanimously, to its nearest CCLE line. Correct calls sit at ≥ 0.68: A/B/C 0.772 / 0.719 / 0.733; real-Flex K562 0.693, HAP1 0.747. Absent lines on Flex: HEK293 0.647 (split), H1 0.653–0.665 (split) | none. The abstention is the design; a non-CCLE line switches on only through 642 | every context's call, n/8, `min_top1_score` and consensus top 5 (identity.json) |
> | postbuild dose (554) | `emitter steering dose per context: median prior_lambda >= 0.9, mde_lambda >= 0.9, call_set_missed <= 100 (554)` | A context's DE design was delivered at low or zero dose (553; 458's programme builds ran at λ 0.24–0.44). Healthy: prior 1.0, MDE 0.951–1.0 (A/D lowest), missed 4–21.5 | `DEF_DAY_ALLOW_LOW_DOSE=1`, only when the sole trip is an MDE median of 0.85–0.90 with prior ≥ 0.9 and missed ≤ 100 (≈ −0.0003 overall). Never for a prior below 0.9 or any dose near 0. Postbuild must run on the pinned builder, and the failed day's exit trap has already restored it: rerun the day in a fresh dir with the override exported (~2 h) | postbuild.json `dose_guard` (per-context values and margins) |
> | HEK293 contingency (642, OFF) | `[642] X: with the hek293 reference row: … -> CALLED HEK293` | Stage 2 named the context HEK293: 8/8, min_top1 ≥ 0.75, and the reference row beats every CCLE line | see "The HEK293 path" below; the coordinator's decision only | identity.json `nonccle_ref`, plus the decision text |
>
> **The Oct 22 identity reading rule (590).**
> - Before the build, copy into the round's RESULTS.md, for every context: the call, n/8 methods, `min_top1_score`, and
>   the consensus top 5 with lineages. `def_identity.py` prints all of them.
> - Read any unanimous call with min_top1 in 0.60–0.70 by hand. The top 5 should be one lineage, show no sister-line or
>   duplicate pattern, and the call's tissue should be plausible.
> - If in doubt, do not switch same-line rows on. Abstaining costs at most the same-line gain; a false switch-on puts
>   another line's rows into a context.
> - **How to switch a context off:** rerun in a fresh dir with `DEF_DAY_SAMELINE_SKIP=<ctx[,ctx]>` (644; the last 10-06
>   step).
>   - `def_sameline --skip-contexts` marks those contexts unmatched whatever their call, so each takes exactly an
>     unmatched context's path: no row blend, no `pt` signs.
>   - With a programme source (a bank, or 604's identity prior) a skipped context still gets the source's programme at
>     w 0, like any unmatched context.
>   - Unknown names stop the day.
>   - Do not use sameline w `0` for this: at w 0 the builder skips the blend, but `pt` still applies the bank's per-target
>     signs.
> - Never `--force` an identity on the day without a written reason in RESULTS.md: `--force` bypasses unanimity and the
>   floor.
>
> **The HEK293 path (642; OFF by default, applied 10-06).** It applies only to a context that abstained (not registered, or
> below the 590 floor) and whose top 5 or lineage suggests HEK293 (kidney or neural, embryonic-like).
> 1. **Diagnostic, read-only:**
>    `DEF_DAY_ALLOW_NONCCLE_REF=hek293 uv run python scripts/def_identity.py --bundle data/challenge --contexts D,E,F --out <work>/identity-hek293.json`
> 2. **Read** `nonccle_ref.per_context[X]`: `accepted` (8/8 and min_top1 ≥ 0.75), `row_top1_methods`, and the top 5.
>    Shallow or odd chemistries reached 6/8 at 0.67–0.69 without being HEK293 (642), so the bar is real, not a formality.
> 3. **The coordinator decides, in writing.** The open question stays open: Orion HEK293T's rows have never been priced on
>    an Arc-Flex HEK293 context. 450's same-line value was measured on H1 only.
> 4. **The day:** prefix the usual `def_day.sh` call with
>    `DEF_DAY_ALLOW_NONCCLE_REF=hek293 DEF_DAY_NONCCLE_DECISION="<who, when, which context, why>"`, in a fresh work dir.
>    - `def_sameline` then reads `orion_hek293t` for that context (294/300 on the A/B/C panel, w 2).
>    - A stage-2 HEK293 call vetoes that context's CCLE call, whatever it was.
>    - Without the decision text, the context is unchanged.
> 5. **The reference is pinned:** `data/identity_refs/hek293_nonflex_ref642.npz`, sha256 `1b829de0…`, built from non-Flex
>    controls only. Never rebuild it on the day, and never from Flex data: a Flex-built row pulled Flex HAP1 to 8/8.
>
> **10-06 knobs.** `DEF_DAY_SAMELINE_SKIP=<ctx[,ctx]>` (644, last): a doubtful context takes the unmatched same-line path;
> logged and recorded. 598's `DEF_DAY_SHARPEN_PT=tercilec:2.5,2.0,1.75` replaces the tercile form once its KAs pass and the
> recipe decides it; `october_recipe.py command` then emits it. Re-run 644's `verify644.py` after 10-06, and again on
> Oct 21.
>
> **PENDING switches** (612, 630, 639/604) are not part of the day until their reads decide them. See the PENDING
> section at the end.

> **2026-09-29 (531's tabletop + 526's review; these override anything older on this page):**
> - **Uploads: one new package a day is the norm.** A 2.77 GB upload took ~7.5 h on 09-29. Start a slot's upload in the
>   morning; start slot 2 only if slot 1 has PUBLISHED by 16:00 UTC and slot 2 is already built and gated. Prebuild the next
>   arm while an upload crawls (never with a pin held). **Stall rule:** check the pending upload's byte range (§6 rule 3);
>   a range unchanged for 30 min is a stall: restart the SAME bytes; a slow but moving range is left alone.
> - **Commit the pins before the first D/E/F submit:** run the day with `DEF_DAY_KEEP_PINS=1`, review and commit the builder
>   diff, then `def_pin_panel_artifacts.py --restore --record <work>/pin-artifacts.json`. The default exit trap restores
>   the builder before you can commit (def_day.sh's success line says to re-run the pin stages otherwise).
> - **Context labels come from `manifest.json`.** Every hand-typed `--contexts D,E,F` / `--expect-contexts D,E,F` below
>   means "the manifest's three names"; a label set other than D/E/F is a proceed once the hand commands use it.
> - **A reordered `gene_names.csv` is a STOP, not a proceed** (531 F6): the pinned CD4 deflation basis (`V5_abc`) raises on
>   a gene-axis mismatch and LODO's same-path assert cannot see it. Do not ship, do not refit the basis on the day.
> - **Preflight now fails closed** on a control file without `ntc_id` when the manifest declares NTC guides
>   (`DEF_DAY_ALLOW_NO_NTC_ID=1` overrides a genuine format change) and on `cells_per_pert != 400`
>   (`DEF_DAY_ALLOW_CELLS_PER_PERT=<n>`); `def_day.sh` refuses a builder that is not HEAD (a reboot mid-pin): restore it
>   from the previous work dir's `pin-artifacts.json`, or `git checkout -- scripts/make_submission_pm_cd4.py` if no record.
> - **If R0 includes `DEF_DAY_PROGRAMME_SOURCE` and the same-line stage dies** (fewer than 150 rows, coverage < 0.5, or zero
>   contexts switched on): rerun at once in a fresh dir with the source unset and flags `noprog`. Record it; do not debug
>   the bank on Oct 22 (531 F9). R0 must not be delayed.
> - **Members hidden:** `capture_board_entry.py` now writes a final entry whose scores are all null as `members_hidden`
>   (exit 0). That is the dark branch (§6), not a capture failure to retry. `experiments/531-cloud-final-phase-tooling/readdef.py`
>   reads captured D/E/F pairs for the feedback plan and exits 2 on hidden members.
> - **The last entry** is an already-scored package, uploaded under the reserved name from Nov 4 16:00 UTC (started by
>   20:00). A sealed confirmation asset (HCASMC, §5) may only CHOOSE among scored packages; it never starts a new build.

> **2026-09-26 (cloud 485, F1):** `def_day.sh` treats a pre-existing h5ad as a successful build. **Every rerun goes in a FRESH work
> dir.** Never "lower q and rerun" in place, or the day can ship the previous build. The fix is 476.

> **2026-09-27 (474, 481's patches, all rehearsed, committed 847cced5):**
> - **`DEF_DAY_CD4_THRESHOLDS=relaxed`**: if setup stops with "CD4 extraction gate fired", check the mapping
>   (`experiments/244-october-rehearsal/mapping244.py`), then rerun with it **in a fresh work dir**. 481 puts the gate's
>   firing probability at 0.29–1.0 on an Arc-style panel. The census floors still gate. Rehearsed on efp (union 186):
>   DONE, ρ ≤ 1, postbuild PASS.
> - `def_postbuild.py` now runs inside `def_day.sh` before the builder restore. The pre-upload diff gate is automatic
>   ("postbuild PASS" in def_day.log).
> - Contexts come from the manifest. ntclevel's Cp46 is normalised by construct size. `--programme-source` refuses
>   < 150 rows.
The whole build path is one gated script, `scripts/def_day.sh`. **Rehearsed end to end 2026-09-25 (452): PASS in 1 h 47
min.** Real identity on `ef-bundle-r`; the same-line mode switched on for D alone. It converged after 3 levelled builds
(ρ 0.9795 / 0.9794 / 0.9790 against 0.9793), audits were clean, the prep dry run verified 360,000 cells, and the builder
was restored. 452b (forced heavy blending) is the second record:

```bash
uv tool upgrade vcc-cli && vcc skill install && vcc datasets list   # 2026-09-26: 0.2.2 lists only `controls` (A/B/C);
                                                     # the D/E/F id may differ -- use what `list` shows
cp -R data/challenge data/challenge-abc              # keep A/B/C first (§2)
vcc datasets download controls                        # the D/E/F bundle into data/challenge
uv run python scripts/october_recipe.py show          # every switch: DECIDED value + evidence, or PENDING + its read
uv run python scripts/october_recipe.py command --bundle data/challenge --work experiments/NNN-def-round/work \
    --name <opaque_joke_stem>                         # prints the exact env + def_day.sh call; refuses while PENDING
#   today (2026-09-27, only LODO decided) it resolves to:
#   DEF_DAY_LODO=1 bash scripts/def_day.sh data/challenge experiments/NNN-def-round/work <stem> 1.75 0.5 noprog 4
```

**The switch values live in `docs/october-decisions.json`** (2026-09-27). After each board read that settles a switch,
record it with `uv run python scripts/october_recipe.py decide <switch> <value> --evidence "<read, entry id, file>"`
(it keeps a history). `command` and `check` also verify that the installed `def_day.sh` and builder carry every knob the
recipe needs (e.g. the zs/hook riders need 502's pass-throughs; an empty flags argument needs def_day ≥ e29e793a) and
that dependencies hold (gate only with the programme source; the programme source forces flags `""`). Tests:
`uv run python scripts/test_october_recipe.py`. The prose decisions below remain the reasoning.

**Downloads and upgrades need no go-ahead** (Ahmet, 2026-10-02, standing; AGENTS §11.0). The day needs two:
`vcc datasets download controls` (the D/E/F bundle, ~1.2 GB like A/B/C's) and a `vcc` CLI upgrade (then
`vcc skill install`). Do both the moment the bundle is released, after checking ≥ 60 GiB free.

What it does, in order (each stage gated; the builder is restored on exit):
1. The driver's setup: preflight, identity, CD4 extract, census, the day's whitener, the day's 410 key, and the pins.
2. The same-line artifact (451). A context switches on only if its identity call is unanimous and its line is in the
   registry. Exported as `VCC_SAMELINE=…:<w>:pt:noprog`.
3. The stack dump, H1 ×2 / CD4 ×1 (438).
4. lvl0, then 448's secant re-level loop until CONVERGED (at most 4 levelled builds; never ship ρ > 1).
5. Audits and a `vcc prep --dry-run`.
It stops before the real `vcc prep` and `vcc submit` (by hand, §2's tail).

Open decisions that change the arguments (decide from the board before Oct 22, and write them here):
- q (arg 4): **1.75, DECIDED 2026-09-26 00:37Z.** 446b transferred on the board (+0.0023; pds +0.0032, mse −0.0076 raw;
  `experiments/446-q175-level/BOARD-RESULT-446b.md`). The curve stops at 1.75: 438 prices q 2 at ≤ 0 calibrated on
  October-like panels, and the ballast's nmae cost grows with q.
- **SUPERSEDED 2026-09-27 (487, the peer's robustness review; coordinator's decision): args 5–6 = `0.5` and `noprog`,
  or `0.5` and `""` when `DEF_DAY_PROGRAMME_SOURCE` is on.** Every likely October match (HCT116 via Orion, K562 via
  GWPS) is a bank already in our donor mix. Its direction lift is ×1.16–1.29 against Pan-on-H1's ×4.37, which is below
  the flip-0.3 bank where w 2 already priced ≈ 0 (mse worse) and w 0.5 +0.0166. `pt` on an in-mix bank is ±0.02–0.04
  against CLEAN4 and bypasses the weak-target gate. Plan the mode at ≤ +0.005 overall. Rehearsal: 510. The text below
  is the superseded 09-26 decision.
  - An empty flags argument (`""`) means "no flags" only in `def_day.sh` ≥ e29e793a (497). The older `${6:-pt:noprog}`
    turned `""` into `pt:noprog`, which with the programme source stops the day and without it brings pt back.
- The same-line flags (args 5–6): **`pt:noprog`, DECIDED 2026-09-26 01:04Z.** 446bja (A's block prior from Nadig Jurkat) read
  no transfer on the board: fid +0.0002 raw against +0.0015 … +0.0040 predicted, overall −0.0011, mostly reach
  (`experiments/450-same-line-programme/BOARD-RESULT-446bja.md`). The programme stays off. Only the per-target rows and their
  per-target prior switch on for a matched context.
  - w = 2 if every matched bank is strong (Replogle GWPS/essential, Nadig, Orion, Sunshine: own-target knockdown in
    ≥ 0.9 of rows). This is 450d's peak, +0.093 on H1.
  - w = 0.5 if any matched bank is weaker (the Jiang lines, VIPerturb, Myllia: 0.62–0.83). 450f: the best w falls to
    0.5 as bank quality falls, and w 2 with a weak bank costs mse.
  - **Since builder `4feae27d` (451b), pass `auto`, the default.** `def_sameline.py` writes each matched context's w
    from its bank class (strong 2, weak 0.5) and the builder reads it. A number still forces one global w.
  - Fallback: if a same-line build dies with `rc_hybrid_k did not meet the per-cell totals`, rerun with w halved, then
    without `:pt`. 450's Gaussian-noise simulation hit it at w 1; the real-bank rehearsals (451r, 452, 452b) never
    did.
- Riders that land on the board (zs `VCC_DE_MAG`, LODO key, repair floor −2, DLD-1 union 2) need their own exports.
  Add each one here as it lands.
- **Two switches are prepared, default OFF (460, 2026-09-26).** Set each only when its board read says so, and write
  the decision here:
  - `DEF_DAY_LODO=1 bash scripts/def_day.sh …`: the stack dump uses 415's LODO key, re-fitted on the day
    (`scripts/def_lodo_key.py`; its `--ka` reproduces A/B/C's `8c82bef3`). **Turn on iff 446bj2 reads key_ok** (09-27
    s1). Its October-path known answer is 460dr vs 453dr (`experiments/460-lodo-dayof/PREREGISTRATION.md`).
    - **DECIDED 2026-09-27 00:30Z: ON.** 446bj2 read key_ok: +0.0022338 overall, mse −0.0145 raw, pds +0.0002, DE
      pinned (`446-q175-level/BOARD-RESULT-446bj2.md`). October runs `DEF_DAY_LODO=1`.
    - **PENDING (2026-09-27): the LODO key's cap.** If 495 (cap 6 on 446bj2) reads "transfers", the day-of fit takes the
      cap as the panel's own p96.647 of carried key/median (= cap 6 on A/B/C; 8–16 on heavier-tailed panels; 543/546).
      That needs `def_lodo_key.py --hi-pct` and a rehearsal. If 495 fails, the cap stays 3.
      - 2026-09-27: `def_lodo_key.py --hi-pct` is in (b072b2cc). Use the unrounded **P = 96.64697522769474**: its KA-pct
        reproduces 495's `dc374281` to 3e-14. The rounded 96.647 lands at cap 6.000055 (10 weights off by 5.5e-5), so it
        is not byte-identical to the arm the board reads. The `def_day.sh` pass-through (`DEF_DAY_LODO_HI_PCT`) goes in
        after 476 lands (497).
  - **PENDING (2026-09-27): per-target bulk sharpening (479).** If 479's pair reads positive, the day's lvl builds take
    `VCC_BULK_SHARPEN_Q=stercile:2.5,2.0,1.75` (S_P2, tie-free; the default unless P2g clearly wins on the board). That
    needs a `def_day.sh` pass-through and a rehearsal. 546 prices LODO p96.647 + S_P2 at +14 … +20 % direction on all 9
    October-like folds.
  - `DEF_DAY_PROGRAMME_SOURCE=h1_2025` with sameline_flags **`pt`** (not `pt:noprog`): every context's block prior
    becomes Arc's own H1 programme (`def_sameline.py --programme-source`; KA-2 reproduces 458's artifact). **Turn on iff
    458 reads fid ≥ +0.0015** (09-30 s1). The `def_day.sh` pass-through is in (9c17856f).
  - **Rehearsed (2026-09-26):**
    - LODO on the A/B/C path: 460dr PASS.
    - Programme source on the A/B/C path: 461dr PASS; the re-level retraced 453 exactly.
    - **Both together on a new panel** (efr): 465 PASS. It converged in 3 levels, took 1 h 43 min, and the builder was
      restored.
  - **A third switch, default OFF (468/472, 2026-09-26): `DEF_DAY_PRIOR_GATE=le1`** (NOT `q33`), only together with
    `DEF_DAY_PROGRAMME_SOURCE`. Targets with design strength ≤ 1 keep CLEAN4 instead of the programme (the builder's
    `VCC_DE_PRIOR_GATE`; `le<N>` needs builder ≥ `efb8c9dc`).
    - `q33` is a share cut on an integer strength. On an essential-free Arc-style panel it is predicted to cut at 1 and
      leave the strength-1 class, the source of H1's reach loss, on the programme (481 F3; the review).
    - `le1` selects q33's exact sets on A/B/C: KA-le1 rebuilds 458g byte for byte. So 09-30's read validates le1.
    - Why: 467 on H1 — reach −0.0388 → +0.0015 with most of the programme's fid kept.
    - **Decide by 458g vs 458d on 09-30** (`468-prior-gate/PREREGISTRATION.md`, with the combined table for 458d's read).
      **Amended 2026-09-27, before any 09-30 entry (487 §2):** with any programme source the gate is `le1`, and a NULL read
      keeps it; only a HARM read (≤ −0.0005) closes it. le1 is the one rule positive with its CI above 0 on every H1
      panel (+0.0087 A/B/C, +0.0073 Arc-style, +0.0056 efr), and Arc-style panels carry more weak targets.
    - End-to-end rehearsal on the A/B/C path: **469dr PASS** (2026-09-26 19:53Z; `469-gated-october-path`,
      `gates469.log`).
      - Config diff vs 461dr is exactly `de_prior_gate`, and the re-level retraced 461dr's factors.
      - The gate's counts on D/E/F (86 / 83 / 83 weak) equal 458g's on A/B/C.
      - `def_day.sh` `5f9bfa26` is committed.
    - The source is the denoised bank **`h1_2025_full`**, not `h1_2025`: 458d replaced 458 on 09-30.

---

## 0. Three things 180's runbook did not know — read before anything else

> ### 0.0 — 244's rehearsal corrections (2026-09-04): these OVERRIDE the lines they name below
> 244 ran all eleven stages end to end on 238's disjoint, essential-free bundle at the pessimistic coverage
> draw (`experiments/244-october-rehearsal/RESULTS.md`). Everything printed what §2 says to look for, except:
> 1. **§2 step 3 and §5 "CD4 gates fire":** on such a panel 047's gates fire (CD4 union 186 < 290) and the driver
>    STOPS — the day is two invocations. Confirm the mapping off the source's `obs['target_contrast_gene_name']`
>    with `experiments/244-october-rehearsal/mapping244.py` (A/B/C must read 291 = 283/287/284, every unmapped
>    target with no case-variant hit), then re-launch with `--cd4-thresholds relaxed`. `target_row_map.csv` is
>    written only by the extraction that succeeds; it cannot be the pre-check.
> 2. **§2's ⚠️ box ("one gap the driver still has"):** closed by 247 (`def_pipeline.py` `verify_existing_operator`);
>    244 observed the three outcomes (A/B/C operator refused on D/E/F, no-digest refused, fresh fit skipped).
> 3. **§2 step 7:** `git status --short scripts/make_submission_pm_cd4.py` must be EMPTY before the pin; an
>    uncommitted edit gets pinned and the meta's `hashes.implementation` becomes unreproducible (244 §7.1).
> 4. **§2 step 9:** SIX `OK` lines on D/E/F, not four — `cd4_sizing 'carrying'` and `cd4_target_g
>    0.3443655278649641` are the day's two extra defaults (247), off A/B/C only.
> 5. **§2 step 10:** `def_postbuild.py` **32/32** on a D/E/F build (31/32 on A/B/C), `moved` =
>    `['cd4_sizing','cd4_target_g','contexts','expect_contexts','whitener_id']`.
> 6. **§7:** if the driver refuses at `prep` in a rehearsal the builder is still pinned —
>    `def_pin_panel_artifacts.py --restore --record <work>/pin-artifacts.json`, then `--stages prep`. The gate
>    that refused read 29.5 GiB free+inactive 40 s after a 20 GiB build.
> 7. **Timing (§2):** 1,692.7 s of stages (whitener 360.9, build 1290.1, prep 15.9); 29.7 min launch-to-prepped
>    artifact including the re-launch; prep 16.25 GiB here, 27.11 in 235 — plan for 25+.
> 8. **§5.-1 / §5.0:** the rehearsal is no longer stale on the panel half. On the ef-bundle the H1 union sized on
>    28 carrying rows (all-rows rule 0.0, solved 0.20199); CD4 182/300, over-dose **1.170×** measured (1.196× on
>    199's draw; the 1.24× planning number is an upper figure) — and with 247's sizing it is not delivered.
>    ef-bundle column for §5's table: K562 146, Orion 293/296, carriers 140/155/5/0, CD4 union 186,
>    `panel_support` 133/182, essential 0. §1's builder hash is three commits stale: HEAD `a17b177` is `24a30055…`
>    (and 249 moves it again).


1. **The whitening operator is a function of the PANEL and fails silently on a new one.** Its
   axis is `union(donor coverage) − panel` (`scripts/build_union_whitener.py`). With the A/B/C
   operator on a new 300, every check the builder runs passes (17,942 columns both times, evec
   shape, file sha256) and **98.8% of the coordinates are paired with the wrong gene**. The
   operator MUST be re-fit for the new panel (7.3 min, §2 step 5) and its sha pinned (step 7).
   Known-answer row: re-fitting for the A/B/C bundle must reproduce `e2f65ecf…` (§4).

   > 🟢 **NO LONGER SILENT — experiment 235, 2026-09-04.** `build_union_whitener.py` now writes
   > an **`axis_sha256`** (sha256 over the operator's ordered column gene SYMBOLS) into the
   > operator's metadata json at fit time, and `make_submission_pm_cd4.py` recomputes it from
   > the panel it is building against and **raises `REFUSING TO BUILD: the whitening operator
   > was fitted for a DIFFERENT PANEL`** — in about a second, before the 1.3 GB eigenbasis is
   > read. The A/B/C axis is **`edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f`**.
   > 235 measured why `n_cols` could never do this job: it is **17,942 on all four panels on
   > disk**, and the meta's `columns_sha256` hashes integer indices, so it reads IDENTICAL for
   > `def-bundle` and `def-bundle-badaxis`. The raise is observed, not asserted:
   > `experiments/235-october-hardening/work/fire235.log` fires it on 199's synthetic panel and
   > shows the pre-235 builder accepting the same inputs with no error at all.
   > **The guard does not replace step 5 — it makes skipping step 5 impossible.**
2. **The CD4 banks hold only the panel's rows** (283/287/284 of A/B/C's 300). On a new panel
   `cd4_complete_targets` is 0 and the build aborts loudly at the residual (`VCC_CD4_TARGET_F`
   raises on a zero-norm residual). The banks are re-extracted from the pinned 16.8 GB source in
   2 s (step 3). Known-answer row: re-extracting the A/B/C panel reproduces the shipped cloud
   banks byte-for-byte on all 36 data files (done 2026-09-03).
3. **Both of those artifacts are pinned by sha256 IN THE BUILDER'S SOURCE**, not in
   `scripts/frozen_inputs.json` (180 moved only the bundle half). So after re-making them the
   builder REJECTS them ("frozen model input hash drift") until `scripts/def_pin_panel_artifacts.py`
   re-pins them (step 7). That is a source edit; it is diffed, verified and must be COMMITTED on
   the day so the meta's recorded builder hash is reproducible.

Everything else 180 established still holds: contexts come from the manifest, the bundle half of
the freeze is re-pinned with `--repin-frozen-inputs`, `panel_support` with `--repin-panel-support`,
and **`vcc prep` needs `--contexts D,E,F`** or it rejects the labels before upload.

---

## 1. Before Oct 22 — verify these exist, monthly and on Oct 21

> **Registry defects found 2026-09-26 (third session, 473 notes). Fix before a D/E/F line can match these banks:**
> - the Jiang banks (A549 / HT29 / MCF7 / BXPC3 / HAP1) and `dld1_lowmoi` store the authors' log2FC (Jiang: mean
>   log2FC × ln2), not our log1p(TS·p) differences. `sameline_blend` unit-normalises them, so their direction
>   over-weights low-expression genes by up to ~2.4×. Convert through the bank's b_ctrl;
> - `dld1_lowmoi` and `myllia_a549` rows carry A/B/C's 300 panel coordinates as exact zeros. On a D/E/F panel those are
>   holes: rebuild the rows without the A/B/C zeroing.

| what | path | sha256 | bytes |
|---|---|---|---|
| CD4 atlas source (Marson 2025 `GWCD4i.DE_stats.h5ad`) | `experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad` | `c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62` | 16,786,240,107 |
| CCLE mirror TPM (cBioPortal `ccle_broad_2025`, DepMap 25Q3) | `data/external/ccle_broad_2025/data_mrna_seq_tpm.txt` | `4b57a32973364aaee93226f1895c67c2abe60042857e648cdcb80cb622e43a11` | 312,778,421 |
| CCLE clinical | `data/external/ccle_broad_2025/data_clinical_sample.txt` | `264b1cf2f571cbaba2c5228d7c16404ed465f03184eddd105023e7872270988e` | 322,069 |
| CCLE parsed cache (matcher input) | `data/external/ccle_broad_2025/ccle_cache.npz` | `1dafee293ac40a8aebbf799bf0d17f08fdafcf4b3576808ceb65a1c9fdff4b73` | 126,728,010 |
| A/B/C production operator (rule 11: this sha, not the filename) | `data/donors/whitener_union_pruned.npz` | `e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123` | 1,287,805,498 |
| its metadata | `data/donors/whitener_union_pruned.json` | `838f263c6a4e1bed069d5faeb3bc3bb2b9234c5a5a933219a94241037a801f4d` | 569 |
| A/B/C frozen bundle table | `scripts/frozen_inputs.json` | `a654d66bc6d093140b09cfe368eb6aa2c0ff591403876588e066f260f500b319` | — |
| builder at the time of writing (git HEAD `ab8aa3e`) | `scripts/make_submission_pm_cd4.py` | `6cc9150122d5e0689eee64ea89fc831c20c3bb80726be75cee335bcdb49c51a6` | — |
| **incumbent package** `h1union` (board **0.2118310**, rank 3/656, entry `iUziqE1gewfcPKzYUDaB`, 2026-09-04) | `submissions/222-h1union-g028.h5ad` | `8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` | 12,648,807,800 |
| its meta (**the config every D/E/F build is diffed against** — `def_postbuild.py`'s default since 240) | `submissions/222-h1union-g028.h5ad.meta.json` | — | — |
| its `.vcc`, and 240's known-answer row 1 | `submissions/222-h1union-g028.vcc` | `f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e` | 2,824,785,920 |
| previous incumbent `rc_callset` (board 0.1918638, entry `31oB9yWy8CSNBDEkY4iW`) — **NOT a valid `--incumbent`: it predates both direction arms** | `submissions/202-rc_callset.h5ad` | `2a894aa543db7f574ebc024dbb3ebd5e36f923827c4e9a01f263092fe49692e0` | 12,648,807,800 |
| previous incumbent `rc_hybrid` (board 0.1834190, entry `nQUJZhZb3MW2uFJHBX5e`) | `submissions/191-rc_hybrid.h5ad` | `e0b0677f984a94a1c4a7af2add196908cdab90d7acb463eb1d6e0397c2d870a8` | 13,202,460,120 |
| **460/461 day-of switch input: Arc's H1 from all public cells** (`DEF_DAY_PROGRAMME_SOURCE=h1_2025_full`; rebuild with `experiments/458-arc-lab-programme/denoise458d.py`, 11 s, from `data/external/vcc2025/*/adata_*.h5ad`) | `data/donors/h1_2025_full/delta.npy` | `26662c527164ddd434ae19a13076e915d1f82de9683c75e27e0af1f5b28591c0` | 21,696,128 |
| 460's A/B/C known answer for `scripts/def_lodo_key.py --ka` | `experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz` | `8c82bef37d0a59b3ed1f4917ccad0d8fd24618a89fd9e2700c058b9aeff4339b` | — |

```bash
# 60 seconds. All three known-answer rows must print PASS.
uv run python scripts/def_preflight.py --verify          # AGENTS.md §2's measured A/B/C numbers
uv run python scripts/def_identity.py  --verify          # A=JURKAT B=HELA C=CAL33, 8/8 variants, anchors
uv run python scripts/def_census.py    --verify          # K562 269, pm 248, cd4 276
shasum -a 256 experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad   # c355f535… (31 s)
vcc --version && vcc prep --help | grep -- --contexts    # 0.2.0 (2026-09-04) has it; re-checked, prep dry-run 16 GB; re-check after ANY upgrade
df -h ~/Code | tail -1                        # need ≥ 60 GB free: 19 GB artifact + 4 GB .vcc + 5 GB scratch
```

**Verified 2026-09-26 10:30Z** (builder 4feae27d, driver 62b5110d):
- all seven pinned inputs above match their sha256;
- `def_preflight --verify` 66/66 PASS; `def_identity --verify` PASS; `def_census --verify` GATE PASS;
- the environment is unchanged (numpy 2.5.2, scipy 1.18.0, h5py 3.16.0, anndata 0.13.2);
- `vcc` is **0.2.2**, not the 0.2.0 recorded above, and `prep --contexts` is present. 453's real prep and dry run
  (2026-09-26) ran on 0.2.2.
- Disk has 70 GB free against the ≥ 60 GB the day needs. Sweep scored `.h5ad`s and packages before Oct 21.

Also keep: `uv` env unchanged (numpy 2.5.2, scipy 1.18.0, h5py 3.16.0, anndata 0.13.2 — the meta
records them); do not upgrade `vcc` in the window without re-running the rehearsal (§4).

---

## 2. Oct 22 — the sequence. One driver, eleven stages, fail-closed.

> ### The configuration the day must carry — check the driver against this line, first
>
> **The package is `h1union`'s configuration on the D/E/F bundle: `VCC_ALPHA=0.982`,
> `VCC_CD4_TARGET_F=0.34`, `VCC_H1_UNION_G=0.28`, `VCC_UNION_SOURCES` unset (= the H1 pair
> `gse295214_h1, h1_2025`), `--emitter rc_callset`.** Board **0.2118310**, rank 3/656, entry
> `iUziqE1gewfcPKzYUDaB`. The two arms that are NOT in 191/202 and what the board paid for them:
>
> | arm | the one variable | Δ overall | entry |
> |---|---|---:|---|
> | `cd4f034` | `VCC_CD4_TARGET_F` 0.10 → **0.34** | **+0.0130557** | `QSe3JlNUYZ2BzvlcsL14` |
> | `h1union` | `VCC_H1_UNION_G` = **0.28** (H1 pair) | **+0.0069114** | `iUziqE1gewfcPKzYUDaB` |
> | | | **+0.0199671** | 0.1918638 → **0.2118310** |
>
> `scripts/def_pipeline.py`'s `SHIPPED_ENV` is that line, as **defaults an exported value
> overrides** — and after the build the driver re-reads the arms off the finished artifact and
> refuses the package if they are not there. Until 240 they were **overrides** pinned to
> `build191.sh`, so this path silently rebuilt the two-arms-ago incumbent and `def_postbuild.py`
> passed, because it diffed against `rc_callset` (experiment 239 §2.1, experiment 240).
> **If a future arm lands, both this table and `SHIPPED_ENV` move together, and 240's row 1 —
> "the driver's defaults reproduce the shipped package byte-for-byte" — is re-run.**

Download the bundle into `data/challenge/` (replace the A/B/C files; keep a copy of A/B/C at
`data/challenge-abc/` first — `cp -R data/challenge data/challenge-abc`). Then:

```bash
mkdir -p experiments/2xx-def-round/work
uv run python scripts/def_pipeline.py \
    --bundle data/challenge \
    --work experiments/2xx-def-round/work \
    --name <opaque-name> \
    --whitener-id union_pruned_def \
    --emitter rc_callset           # the incumbent's emitter (§3); the driver's default
```

That runs, in order, stopping at the first failure and writing `<work>/timings.json` plus one log
per stage under `<work>/logs/`. What each stage is, what to look at, and how long it took in the
rehearsal (pessimistic synthetic panel, 2026-09-03):

| # | stage | command (the driver prints the exact one) | look for | rehearsal wall / peak RSS |
|---|---|---|---|---|
| 1 | `preflight` | `scripts/def_preflight.py --bundle data/challenge` | **66/66 checks pass**; contexts `['D','E','F']`; `n_genes 18533`, `cells_per_pert 400`, 3 × 18,400 control cells; gene axis read from the file; both CSVs have a header row (`vcc prep` says "headerless" — it skips it itself) | 4.5 s / 0.3 GiB |
| 2 | `identity` | `scripts/def_identity.py --bundle data/challenge --contexts D,E,F` | the consensus top-1 per context and whether all 8 feature/rank variants agree; calibration anchors K562/HCT116 at rank 1–2, Jurkat 1–3 (041's table). **Write the three names into the round's RESULTS.md.** They are cancer lines (§2.1). What to do with them: §3, row "auxiliary donors" | 14.7 s / 2.7 GiB |
| 3 | `cd4-extract` | `experiments/047-cd4-direct/extract_cd4.py --input <source> --output <work>/cd4_compact --panel-csv data/challenge/pert_counts.csv --genes-csv data/challenge/gene_names.csv` | `panel coverage {Rest, Stim8hr, Stim48hr}` and `union`. Expect ≈ 283–296 per condition, union ≥ 290 (A/B/C: 283/287/284/291). 047's gates (275 / 290) are ON by default; if they fire, that is a finding, not a threshold to weaken — see §5 | 2.2 s / 0.3 GiB (+9–31 s source sha) |
| 4 | `census` | `scripts/def_census.py --bundle data/challenge --cd4-root <work>/cd4_compact` | the carrier histogram over the three production donors (3/2/1/0), `panel_support predicted` = the two numbers step 8 will pin (A/B/C: 248 / 276), and the panel's overlap with the essential list (A/B/C: 0). Responses in §5 | 0.4 s |
| 5 | `whitener` | `scripts/build_union_whitener.py --bundle data/challenge --basis-prune essential --out <work>/whitener_union_pruned_def.npz --scratch <work>/.whitener_scratch --max-rss-gib 14` | `union axis: 17942 genes` (may differ by a few); `basis rows … = 29394`; `dropped_in_panel` (expect ≈ `[]`); a NEW sha256 (the driver aborts if it equals the A/B/C one). **AND THE AXIS HASH, which is the check that actually identifies the panel (235):** the fitter now prints `axis sha256 <64 hex>` and `panel sha256 <64 hex>` right after the union axis line and writes both into `<out>.json` as `axis_sha256` / `panel_sha256`. **It MUST NOT be `edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f`** — that is A/B/C's, and reading it here means the bundle on disk is still A/B/C. `n_cols` cannot substitute: 235 measured it at **17,942 on all four panels**. Step 9's builder re-derives this digest from the panel it is building against and refuses to run unless it matches, so a skipped or mis-pointed step 5 can no longer reach a build. The 14 GiB guard is needed because the 4.2 GB basis memmap is counted while resident; the real anonymous footprint is ~8 GiB at the eigensolve | **437.6 s** / 8.1 GiB at eigh, 12.9 GiB process high-water (235 re-measured 2026-09-04 on A/B/C: **320.0 s and 364.1 s**, 13.7 GB `ru_maxrss`, 5.54 GiB anonymous at the covariance) |
| 6 | `repin-frozen` | `make_submission_pm_cd4.py --repin-frozen-inputs --bundle data/challenge` | the printed old → new diff: 3 contexts removed, 3 added, `manifest.json` and `pert_counts.csv` moved, `gene_names.csv` unchanged (if it moved, STOP: the gene axis changed — the donors' `genes.txt` mapping is still valid, but read §5), `panel_support → None` | 2.1 s |
| 7 | `pin-artifacts` | `scripts/def_pin_panel_artifacts.py --whitener <npz> --whitener-metadata <json> --whitener-id union_pruned_def --cd4-root <work>/cd4_compact --apply` | the unified diff: ONE new `WHITENER_REGISTRY` entry (additive; `union_pruned` untouched) and 15 replaced `cd4/…` hashes, nothing else. Then `git diff scripts/make_submission_pm_cd4.py` and **commit it** | < 1 s |
| 8 | `repin-panel-support` (effects only; 180 called this "4 min", it is 9 s here) | `make_submission_pm_cd4.py --repin-panel-support --bundle data/challenge --expect-contexts D,E,F --whitener … --whitener-id union_pruned_def --cd4-root … --min-cpm 10 --max-rss-gib 24` with `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28` (§2's box; the driver prints the resolved values and marks any the operator supplied) | `[contexts] ['D','E','F'] (manifest, frozen table and --expect-contexts all agree)`; `[panel-support] UNPINNED … Observed …`; then `[repin] panel_support pinned to {…}`. **The two numbers must equal step 4's prediction exactly.** The driver refuses to start this or step 9 if another builder is running or free+inactive memory < 26 GiB | **9.2 s** / 19.8 GiB |
| 9 | `build` | the shipped flags byte-for-byte from `experiments/222-h1-union-build/build222.sh` — which is 210's, which is 202's and 191's, the flags never moved — (`VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28`, `VCC_UNION_SOURCES` unset, `--min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0`) plus `--bundle data/challenge --expect-contexts D,E,F --whitener … --whitener-id union_pruned_def --cd4-root … --emitter rc_callset -o <work>/<name>.h5ad` | **`[operator] axis_sha256 <hex> over N columns MATCHES metadata:<json>`**, printed within the first minute (235) — if instead it raises `REFUSING TO BUILD: the whitening operator was fitted for a DIFFERENT PANEL`, step 5 or step 7 pointed at the wrong file and NOTHING after it is worth reading; then the builder's own audit; `peak_rss_gib` and the new `operator.axis_sha256` / `operator.panel_sha256` in the meta; and last, 240's `[arms]` block — SIX `OK` lines (see §0.0 item 4) reading `cd4_target_f 0.34`, `h1_union_g 0.28`, `h1_union_sources ['gse295214_h1','h1_2025']`, `alpha 0.982` **off the written meta**, or the driver raises `REFUSING THE BUILD` | **1255.8 s** with `rc_callset` (1287.1 s with `rc_hybrid`) / 19.9 GiB — 240 re-measured the SHIPPED configuration on the A/B/C bundle at **1302.3 s / 20.00 GiB**, and `vcc prep` at **26 s / 19.72 GiB** |
| 10 | `postbuild` | `scripts/def_postbuild.py --meta <name>.h5ad.meta.json --bundle data/challenge --cd4-root <work>/cd4_compact --incumbent submissions/222-h1union-g028.h5ad.meta.json --allow-config-diff whitener_id` (that `--incumbent` is now the default; **never point it back at 202/191** — 240 checks the reference itself carries the arms and fails if it does not) | **20/20 PASS**: operator sha ≠ A/B/C's, operator axis recomputed from THIS panel matches, pm/cd4 counts match the set-intersection prediction, config differs from the incumbent in exactly `contexts`, `whitener_id` (and `targets` if 300 changed), **and the five shipped arms are present by value** (`alpha` 0.982, `cd4_target_f` 0.34, `h1_union_g` 0.28, `h1_union_sources` the H1 pair, `emitter` rc_callset) | 1 s |
| 11 | `prep` | `vcc prep <name>.h5ad -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv --contexts D,E,F -o <work>/<name>.vcc --dry-run` | `✓ dry run`, `cells 360000`, `D=120000, E=120000, F=120000`, `targets verified`, `counts-preserved` | **11.8–13.9 s** / **22–24 GiB** — the CLI loads `X`; the driver guards it like a build: never alongside one |

> ✅ **CLOSED (531 F13, 2026-09-29): the driver now refuses a stale operator** (`verify_existing_operator` in
> `scripts/def_pipeline.py` compares the stored axis sha with this bundle's and refuses the skip). The historical note below
> is kept for the record only; delete the two files and refit if it fires.
>
> ⚠️ *(historical)* **One gap the driver still had, and 235 could not close it** (`scripts/def_pipeline.py` was
> outside 235's authorised files). Its `whitener` stage **skips the fit entirely if
> `<work>/whitener_<id>.npz` and its json already exist**, with no check on either — so a stale
> operator left by an earlier run is silently reused, and the driver's `sha == ABC_OPERATOR_SHA`
> abort never sees it. Since 235 the build (step 9) catches this, but 20 minutes later. **The
> one-line fix, for whoever owns that file next:** in the skip branch, read
> `json.loads(wmeta.read_text())["axis_sha256"]` and compare it against the axis recomputed from
> **this** bundle (`experiments/235-october-hardening/axis235.py:axis_for` is 15 lines and needs
> only the donors' `genes.txt`), and refuse the skip if it differs or is absent.

**Total for one panel, bundle download to a prepped artifact: 1,748 s ≈ 29 min of wall clock**
(measured 2026-09-03, `work/run1/timings.json`), of which 28.6 min is the operator fit (7.3), the
effects run (0.2), the build (20.9) and `vcc prep` (0.2). Add ~15 min of upload.

Then, by hand (never by the driver):

```bash
vcc prep <work>/<name>.h5ad -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv \
    --contexts D,E,F -o <work>/<name>.vcc                                   # writes the .vcc (2.8–3.7 GB)
git add scripts/make_submission_pm_cd4.py scripts/frozen_inputs.json && git commit -m "D/E/F pins"   # BEFORE submitting
vcc submit <work>/<name>.vcc -m "<next MODEL-NAMES.md queue name>"   # -d EMPTY. The LAST entry before Nov 5 23:59Z uses the reserved final name (docs/MODEL-NAMES.md). Upload ≈ 13–16 min at ~3.8 MiB/s (a22: 3.5 GiB in 15m31s)
```

Record entry id, package sha256, joke name and the meta path in the round's RESULTS.md and in
`submissions/SUBMITTED_LEDGER.json`.

---

## 3. The blind decision table — every shippable flag, pre-decided from A/B/C evidence

`n_real` cannot be measured on D/E/F and the board is dark, so nothing below is tuned on the day.
Each row: what ships, the instrument that validated it, and the pre-decided rule if the D/E/F
controls or census look different. "Board" means a scored entry in `docs/BOARD-LEDGER.md`.

| flag / choice | SHIPS | validated by | if D/E/F look different |
|---|---|---|---|
| **`VCC_ALPHA`** | **`0.982`, global** (a15 / a23 / 191 / 202) | board series a11 → a14 → a15; 187 solved `mse`+`nmae` jointly: the global optimum is 0.982 to four decimals. **`VCC_ALPHA=knee` (per-context) is board-FALSIFIED**: abc-knee (entry `N3ui4hCkKx9bxExG2pxE`) scaled Δ`mse` **−0.00226** against a local +0.005 | The alternative is RECORDED, not shipped: 187 §6 — under a 1% cliff uncertainty the raw knee's expected gain is ≈ 0 and the expected-value optimum is **knee × (1 − 1.25%)**, worth **+0.0009** overall over global 0.982. Pre-decided trigger: run the knee diagnostic (180's `knee180.py` pattern, or `VCC_ALPHA=knee VCC_ALPHA_REF=0.982` on the effects-only stage) and read the per-context knee alphas. If **every** context's knee alpha is within **[0.982, 1.13]** (A/B/C: 0.986 / 1.013 / 1.046), ship 0.982 — the gain is inside the proxy error. Only if a context's knee alpha exceeds **1.13** (≈ ρ < 0.80 under 0.982, where 180 §8 prices the forfeit at ≥ 0.018 scaled `mse`) ship **knee × 0.9875 for that context** and 0.982 for the others; never the raw knee, never above it (10% long costs 0.090 scaled). The coordinator ratifies before the build |
| **`VCC_WHITEN_GAMMA`** | `0.65` | 170: end-to-end V on the challenge panel with its vertex at 0.65 (raw `mse` 0.8813 / **0.8738** / 0.8795 at 0.60 / 0.65 / 0.70), every gamma in [0.35, 0.80] net negative on all six; a19 board −0.0135 (gamma 0.80 + ridge). Mechanism = the ρ cliff, panel-independent | never move blind |
| **`VCC_WHITEN_LAM`** | `0.0` | 198: production has no ridge; a19's `lam 0.007` was 79% of its loss | never |
| **operator** (rule 11: by sha) | A/B/C: `data/donors/whitener_union_pruned.npz` `e2f65ecf…`. **D/E/F: the re-fit `union_pruned_def`, sha from step 5, pinned in step 7. The A/B/C file must NOT be used** (§0.1) | 091: essential-pruned basis +0.002409 avg6 (t = 19.85) vs a size-matched random prune; a15 shipped it | mechanical, no judgement: refit, check the sha changed, check `dropped_in_panel` ≈ `[]` |
| **emitter** | **`rc_callset`** (202; board **0.1918638**, +0.0084448 over `rc_hybrid`, `jac` +0.0276 scaled, `reach` +0.0226 scaled, `pds`/`mse` reproduced to 16 s.f. for the third time) | board, entry `31oB9yWy8CSNBDEkY4iW` | context-independent (exact-marginal emission at bit-identical pseudobulk; the identity is a repeated board measurement, not a construction argument). Never `resample_mix` (a21 board `reach` −0.0185; 177/179/181 dead). `rc_hybrid` (191, 0.1834190) is the fallback if `rc_callset` ever fails a falsifier |
| **call-set policy** | 192's `top 400 by |H| ∪ top 2,200 by control CPM` (2,600 per (context, perturbation)); `rc_callset_head 400`, `rc_callset_cpm 2200`, `rc_callset_z_cap 1.9` | 192 h1-test +0.009 to +0.011 predicted; **board +0.0084448 realised** (202), `jac` positive for the first time | `n_real` is unreadable on D/E/F (rule 9), but §2.1 says the panel is binned by `n_real` BY CONSTRUCTION, so the A/B/C regime transfers by design. Do not re-derive; ship the board winner |
| **`VCC_CD4_TARGET_F`** | **`0.34`** (A/B/C effective fraction **0.32207201067578145**) | **board `cd4f034` +0.0130557, entry `QSe3JlNUYZ2BzvlcsL14`, 2026-09-04** — the first arm ever to move the DIRECTION channel (`cos` 0.07228 → 0.09018); 207 swept 25 doses × 4 PLOLO folds and 212 measured +13.83% PLOLO cosine at this dose; the pre-decided rule in §3.1 fired and 0.34 superseded 0.1. Before it: a15 board +0.00366 vs a14 at `f = 0.1`; 062 LODO gains on HCT116 / HEK293T / K562 folds — i.e. **not** Jurkat-specific, so it ships whatever D/E/F's lineages are | a board-validated dose before Oct 22 supersedes 0.34, and `SHIPPED_ENV` moves with it. If CD4 complete coverage of the new panel is low, the builder already preserves incomplete targets before the global scale (`incomplete_targets_preserved_before_global_scale`); nothing to tune by hand. If step 3's gates fire (< 275 / 290), see §5 |
| **`VCC_H1_UNION_G`** | **`0.28`** | **board `h1union` +0.0069114, entry `iUziqE1gewfcPKzYUDaB`, 2026-09-04** — a fourth donor estimate as a second un-whitened post-whitening residual, sized on its carrying rows. **ALL SIX members rose**, `nmae` included (+0.00588), which is the signature of a SIGN gain and not a magnitude rebalance (docs/FINDINGS.md §8.1d); 218 measured +7.680% macro cosine and `pds` +0.004887 across three cancer PLOLO folds, positive on 12 of 12 (fold × dose) cells | the residual is sized on the rows the sources actually carry, so a panel the H1 pair covers less well shrinks the arm rather than dropping it; 235 made a source that carries NOTHING raise instead of contributing a silent zero. If §5.-1's D/E/F coverage note fires, read it there — but the flag itself is not tuned on the day |
| **`VCC_UNION_SOURCES`** | **unset** = the H1 pair `gse295214_h1, h1_2025` | 229: the four-source stack `stack4` was **−0.00102** marginal and `stack6` **+0.00268** against a +0.0019 floor — inside the ordering noise, so nothing supersedes the pair. Registry membership is inert (229's reg10 row: the ten-entry registry rebuilds the six-entry package byte-for-byte) | never set it blind. A new source is an arm with its own A/B/C validation, not a day-of choice |
| **`VCC_CD4_DEFLATE_K` (A/B/C) / `VCC_CD4_DEFLATE_BASIS` (D/E/F)** | **A/B/C: `K=5` (the re-fit the incumbent was built with); D/E/F: `BASIS=V5_abc` (the PINNED A/B/C basis, 251), NEVER the re-fit** — driver defaults by panel since 256 | **board `cd4svd5` +0.0042634, entry `Xv0ayvuoHpsKYWv2iu42`, 2026-09-05** — raw `pds` +0.01492 (k_pds 1.22), raw `mse` +0.00454 (re-shape cosine factor 3.76); §3.1's 250 trigger fired; 251 KA 2 (pin == re-fit on A/B/C to 2.4e−15); 256 KA 1 (driver defaults rebuild `aba27758…` byte-for-byte) | never change k or the basis on the day (250: a re-fit reads ≤ 0 off A/B/C); record `captured_frac_of_centred` (0.10–0.14 measured; outside 0.05–0.19 is unmeasured) |
| **`VCC_MASS_FIX`** | unset | a23 (`minnorm`) board −0.00073: `pds` +0.0017 and `mse` +0.0015 as predicted, `reach` −0.0079 unpriced | never blind |
| `--overdispersion 4.0` | 4.0 | a8 board (`reach` peak 0.1890 raw) | never |
| `--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0 --min-cpm 10` | as listed | a12 (tiered spread) and a14 (per-context cap) board points; a15 carried them | never |
| donor weights | `k562_gwps : orion_hct116 : orion_hek293t = 1 : 1 : 1` | v5 / §4.2 (equal weighting is the best configuration under union whitening) | coverage-driven behaviour is the builder's own: a target carried by < 2 donors gets an all-zero agreement key and no shrinkage; 0 carriers → panel-mean fallback. Nothing to set. §5 sizes it |
| PerturbMatch | 50:50 K562-slot replacement, only where both support the target | 061 (+0.00009 avg6, non-negative prerequisite); a15 carried it | mechanical; count pinned in step 8 |
| seeds | `seed 0`, `emission_seed 1000003`, `rc_jitter_seed_base 1000`, `library_pick_seed 0` (builder defaults) | re-seeded twins differ by 0.2–4.8% of span (§4) | never |
| **auxiliary donors chosen by identity** (the CCLE call) | none beyond the CD4 residual | 041 → 062 → a15: the one time identity selected a donor it shipped as a board-validated residual | **Pre-decided:** an identity may add a donor ONLY through a mechanism that already has a board point (the post-whiten residual), never through a new mechanism, because no D/E/F instrument can price a new mechanism. If a D/E/F line has public Perturb-seq of its own (HepG2, RPE1, A549, MCF7, HT29, HAP1, K562, HCT116, HEK293T, Jurkat), that is in-domain TRAINING data (§7 boundary: condition a model, never look up an effect) and is a coordinator decision with its own A/B/C-style validation — which must be finished before Oct 22 or not attempted |
| `n_real` regime | a15's (whatever the emitter) | rule 9 | unmeasurable; do not infer it from the dark board |

### 3.2 Updates from 2026-09-22/23: the board-validated bulk/DE separation (rows to add to the table above)

| flag / choice | SHIPS | validated by | on D/E/F |
|---|---|---|---|
| **`VCC_PERT_SCALE_BALLAST=1`** (ballast cells) | on | **board 410b 0.2365722 (+0.0171 over the welded 410)**. The four DE members sat within 0.0007 raw of 395 while mse kept the key's −0.079 | panel-agnostic. Check `ballast_jk_ratio` (median ≈ 1.02) and `ballast_clip_mass` (≤ a few thousand) in the meta |
| **`VCC_PERT_SCALE`** (per-target key) | `wshare_lo0.3_hi3`, **rebuilt for the new panel** | **board 410: raw mse −0.0797** (the key transfers); 410b as above | the artifact is panel- and axis-specific. `build_weights410.py LO HI` with the day's refit operator (its Oct-22 guard refuses the A/B/C axis). Candidates on the ladder: `lodo415_lo0.5_hi3` (417i/417j) |
| **`VCC_BULK_SHARPEN_Q`** | `1.25`; `1.5` if the levelled dose step is board-positive | **board 417a 0.2435334 (+0.0070)**: pds +0.0064, mse −0.0249 raw, DE pinned. Transfer k_sh 0.42–0.65 of PLOLO cancer folds (kolf2's fold predicted it) | no artifact. Applies at render, with the new panel's columns protected automatically |
| **`VCC_BULK_CTX_SCALE`** (re-level) | the EV-optimal factors from `ntclevel417.py` | pending on the board (417g/417h, 09-24). The tool's KA passed: NTC-proxy claims == 187's to ≤ 0.4 % | run `ntclevel417.py` on the day's un-levelled build with `--sd 0.01`, then rebuild with its string |
| **bulk composition via `VCC_BULK_DUMP`/`VCC_BULK_LOAD`** (stack H1 ×2 + CD4 ×1.5) | pending (417e, 09-24) | KA-LOAD: a same-env LOAD reproduces 417a's X byte for byte | a dump is tied to the bundle (LOAD refuses cross-bundle). On the day: dump with the stack env, then LOAD into the shipped-env build |
| **`VCC_CARRIER_GAMMA`** | **unset** (γ 0.65 carriers) | board 417c: γ 0 carriers −0.0017 (fid/reach up, nmae +0.0089 raw) | never set it blind |

### 3.1 Updates from 2026-09-03/04

> **UPDATE 2026-09-04 (250, `experiments/250-deflation-october/RESULTS.md` §6) — `VCC_CD4_DEFLATE_K` / the pinned basis:**
> unset until the 2026-09-05 `cd4svd5` board point is read. If board raw `pds` rises against 0.83572, ship the
> deflation on D/E/F **as the PINNED A/B/C basis, k = 5** (`U_5` sha `f06e0262…`; the path exists: `VCC_CD4_DEFLATE_BASIS=V5_abc`, artifact `data/donors/cd4_deflation_basis/V5_abc.npz` sha `ba7ade3f…`, 251, rehearsed on the ef-bundle 40/40) —
> **never a day-of re-fit**: on disjoint panels the re-fit costs 4–14× A/B/C's cosine and reads ≤ 0 where CD4
> carries < 250/300 rows, while the pinned basis is never negative (+0.0061 / +0.0009 / +0.0014 / +0.0030 on
> `real` / `pess` / `efp` / `efr`). Never use `VCC_CD4_DEFLATE_K` (the re-fit) on D/E/F. Never change k on the day. Record the fraction of the day's centred
> CD4 energy captured by the pinned `V_5` (0.10–0.14 measured; outside 0.05–0.19 is unmeasured territory).
>
> **UPDATE 2026-09-05 01:20Z (BR249, `experiments/249-cd4-deflate/BOARD-RESULT.md`) — the trigger above FIRED:** board raw `pds`
> 0.83572 → **0.85064** (`cd4svd5` 0.2160944, +0.0042634). **The deflation SHIPS on D/E/F as the PINNED A/B/C basis
> `VCC_CD4_DEFLATE_BASIS=V5_abc`, k = 5**; A/B/C's incumbent carries `VCC_CD4_DEFLATE_K=5` (251 KA 2: the two agree to 2.4e−15 on
> A/B/C). The cosine cost transferred at 3.76× PLOLO's aggregate read (falsifier 2), so k is NEVER raised on this basis (k = 10 dead,
> 255) and every further re-shaping arm is priced on that factor. **Driver wiring landed 01:50Z (256): `def_pipeline.py` defaults carry the re-fit k = 5 on A/B/C and `VCC_CD4_DEFLATE_BASIS=V5_abc` on any other panel, `def_postbuild.py` requires them by panel with `cd4svd5` as the reference, and 240's row 1 re-run reproduces `aba27758…` byte-for-byte (`experiments/256-pinned-basis-wiring/RESULTS.md`).** The table row above still reads "unset until the point": read this update as the row. **Provisional (2026-09-05 03:00Z): BR249's prereg falsifier 2 fired (the cosine cost transferred at 3.76× PLOLO's read) while falsifier 1 did not and both branch-A conditions held; the SHIPS decision is BR249 §3's reasoning, re-decided at the next A/B/C point that re-shapes the CD4 rows.**
>
> **UPDATE 2026-09-04 (241, `experiments/241-ablation-reprice/RESULTS.md`) — the two rows this section said must be re-decided, re-decided:**
> - **PerturbMatch: SHIP IT, and retire the "−0.001615, inside the floor" citation.** On the shipped renderer at the shipped dose the row is
>   **−0.002178** (−0.002579 at `f` 0.10), i.e. outside the ±0.0019 gate — but split by channel (§8.1b rule 3) the DIRECTION half is +0.001403
>   (0.74× its gate) and the ORDERING half +0.000774 (0.17 sd of the 0.004502 ordering noise). Neither half clears its own bar, so dropping it is
>   not a measurable arm; 231's admissibility test passes on 241's number (implied board `reach` 0.28289 ≤ 0.290835), so 231 and 241 agree.
>   Dominant uncertainty: its `pds` credential is `kolf2` (+0.003233 pluripotent vs +0.000342 on the three cancer folds); on cancer folds alone
>   dropping it reads +0.002713. Mechanism: a magnitude arm, `e` = 0.162 (§8.1d).
> - **CD4 residual row: 200's +0.002635 was understated 2.07–2.46×** (its four DE members are a structural zero on h1-test: coverage 4/85);
>   filled from its own PLOLO cosine move it reads **+0.0055…+0.0065** at the shipped dose, direction UP, against the board's +0.0130557.
>   Ship `f` 0.34 (216: argmax). On D/E/F size CD4 on its carrying rows (247) — the all-rows rule over-doses 1.24–1.42× at the pessimistic draw.
> - 200's ten-component ORDERING is unchanged on the shipped renderer (no rank move larger than its error bar); the six censored `mse` rows are
>   each worth **+0.029925** on the shipped stack (1.6× what 200 could see), so the load-bearing components are load-bearing by more, not less.
 — read these ON TOP of the table above

| row | change | evidence |
|---|---|---|
| **`VCC_WHITEN_GAMMA`** | **Closed a THIRD time, now on the right endpoint.** 212 re-swept it with **cosine** primary (170 and 195 both optimised `pds`, the member least sensitive to direction). The two maximisers do differ — cosine wants γ = 0, `pds` wants γ = 0.80 — and the composite barely moves: best **+0.00013** on `rc_callset`, **+0.00042** on `cd4f034`, most generous anywhere **+0.00321** at γ = 0.40 against a +0.0089 floor. **Structural reason it can never pay:** the operator trades cosine for retrieval at **96.2%** of the board's own exchange rate (`d mse*/d cos` +0.3655 vs `d pds/d cos` −0.3517, cancelling to 3.8%), and that ratio does not move with the emission, the dose, or α. **Never move it blind, and now: never move it at all.** | 212 |
| **`VCC_ALPHA`** | **Closed a SECOND and final time, on the current emission.** 211 re-solved on `rc_callset`: **α\* = 0.9808** (worth +0.000076) and on `cd4f034` **α\* = 0.9879** (worth +0.00039) — both far inside noise. Clearing +0.0089 on amplitude alone needs **α 1.147**; the binding context permits **−0.12%**. **`claim` is at 100.0000% of its cap in the one context that reads it** (context A, `q_A` = 1.000000), which the board already proved: `rc_hybrid`'s and `rc_callset`'s Δ`mse` of exactly 0 to sixteen significant figures **is** that proof. **Ship 0.982.** The knee rule in row 1 stands unchanged as the only conditional. | 211 |
| **standing rule (new)** | **No emission arm may be re-priced on amplitude grounds.** 208's "the `mse` cliff is 2.6× further away" measures the **claim-saturation** margin (Ĉ/C_p), not the **ρ = 1** amplitude cliff — different sides of the same `min`. That re-read is correct **only** for the depth knob (196 §2g's `m = 2, φ = 8`). | 211 |
| **`VCC_CD4_TARGET_F`** | 207 swept 25 points × 4 PLOLO folds (744 configs): the cosine maximiser is `f = 0.80` (**+31.5%**), `q_w`'s is 0.80, `pds`'s is 0.20 — **the shipped 0.10 is below all three**, because 062 chose it on LODO `pds`. Ships **`f = 0.34`** by three pre-stated criteria; 212 independently measured **+13.83% PLOLO cosine** at that dose. **`cd4f034` (entry armed for the 2026-09-04 00:00 UTC slot) decided this row: it landed at +0.0130557 and `f = 0.34` SUPERSEDES `0.1` for D/E/F** (board 0.2049196, entry `QSe3JlNUYZ2BzvlcsL14`). §3's row and `scripts/def_pipeline.py`'s `SHIPPED_ENV` both carry 0.34 as of 240. Mechanism that makes it panel-independent: CD4 must enter **un-whitened** (a whitened residual is flat at every dose), the shuffled control is strictly negative at all 12 doses, and the gain is **+7.6% to +15.2% on four unrelated non-T lines — lineage-GENERIC**, so it carries whatever D/E/F's lineages are. | 207, 212, 210 |
| **`VCC_H1_UNION_G`** (new row) | **Ships `0.28`.** 222's arm landed at **+0.0069114** (board 0.2118310, rank 3/656) with **all six members up**. 229 then tested whether other sources should ride with or instead of the H1 pair: the pre-registered four-source list was **−0.00102** and the six-source stack **+0.00268**, both inside the ordering channel's own noise (sd 0.004502, §8.1b), so **the pair ships and `VCC_UNION_SOURCES` stays unset**. | 218, 222, 229 |
| **the driver itself** | 🔴 **239 found `scripts/def_pipeline.py` pinned to `build191.sh` and applying its env as an OVERRIDE**: run as §2 documents it, the Oct-22 path discarded BOTH arms above — **−0.0199671**, more than our whole gap to first — and `def_postbuild.py` PASSED, because it diffed against `rc_callset`. **240 fixed it**: the env is defaults an operator overrides, the defaults are the shipped configuration, the postbuild reference is the current incumbent and the arms are asserted by value, and the driver re-reads them off the finished artifact. | 239, 240 |
| **PerturbMatch** | ⚠️ **200's ablation makes it the only NET-NEGATIVE component of the ten: −0.001615 six-member.** Inside the ±0.0019 floor as a standalone arm, which is why 200 declined to drop it. **Experiment 214 is testing whether it turned negative because 207 tripled the CD4 dose** — 064 measured the two as *synergistic* at the old dose. **If 214 finds them substitutes, this row and the CD4 row must both be re-decided before Oct 22.** Until then: ship it, as the table says. | 200, 214 |
| **operator trap, priced** | §0.1's silent panel-mismatch failure now has a number: **−0.0429 overall on `pds` alone** — roughly a quarter of our entire score, silently, on the only entry that counts. Step 5 is not optional. | 212 |

**Closed for good as of 2026-09-04, do not spend D/E/F time on any of them:** more donor lines
(166, 205), cross-line adaptation (184, on 50 Tahoe lines), learned transforms (195), the operator's
eigenvalue map (193), **the operator's exponent (212)**, donor re-estimation (194), donor weighting
(204), donor–donor correlation (209), the hierarchical estimator (203), **global amplitude (211)**,
ordering keys (189, 190). The live questions are the **transport law** (213) and **component
composition** (214).

---

## 4. Known-answer rows the day must reproduce (rule 10's corollary)

| row | command | must print |
|---|---|---|
| preflight | `def_preflight.py --verify` | AGENTS.md §2: 18533 / 300 / 400 / 18400 / 46 / 360000; A: nnz 109,910,478, library 20,109 [3,275, 52,420] |
| identity | `def_identity.py --verify` | A = JURKAT, B = HELA, C = CAL33; 8/8 variants unanimous; 24 anchor ranks identical to 041's `results.json` |
| census | `def_census.py --verify` | K562 269, pm 248, cd4 276 |
| CD4 extraction | `extract_cd4.py --input <source> --output <tmp> --dry-run` (A/B/C panel) | coverage 283 / 287 / 284, union 291, 10,113 genes |
| operator refit | `build_union_whitener.py --bundle data/challenge-abc --basis-prune essential --out <tmp>.npz --expect-sha256 e2f65ecf…` | `MATCH` — 409.2 s, 10.95 GiB high-water, `e2f65ecf…` reproduced exactly on 2026-09-03; re-run twice under 235's change on 2026-09-04 (320.0 s / 364.1 s), both `MATCH`, both writing `axis sha256 edf6ba88…`, the two jsons byte-identical |
| **operator axis guard** (235) | `uv run python experiments/235-october-hardening/fire235.py` | `ROW 5 PASS` — the A/B/C operator LOADS on the A/B/C panel and RAISES on 199's `def-bundle`, on `def-bundle-badaxis`, and for any operator carrying no axis digest, while the pre-235 builder accepts the same wrong panel silently. ~2 min, 16.5 GiB (the three-donor model) |
| **carrying-row sizing** (235) | `uv run python experiments/235-october-hardening/tests/test_october_hardening.py` | `9/9 passed`, ~2 min, 1.0 GiB. Includes the shipped `h1union` residual frame reproduced to the last bit and every registered source sizing above zero on a panel it covers |
| A/B/C build under the pinned builder | `experiments/202-callset-rider/build202.sh rc_hybrid` | sha256 `e0b0677f…` on 13,202,460,120 bytes (202 §1, row 1) |
| **A/B/C build under the SHIPPED configuration** | `experiments/222-h1-union-build/build222.sh` with `VCC_H1_UNION_G=0.28` (or 229's `build229.sh`, `VCC_UNION_SOURCES` unset) | `.h5ad` sha256 `8df4dadd…` on 12,648,807,800 bytes; its `.vcc` `f13ce40b…` on 2,824,785,920 (222; reproduced by 229 row 2 and 235 row 2) |
| **the DRIVER's defaults reproduce that package** (240 row 1 — the row that proves the Oct-22 path builds what we ship) | `def_pipeline.py --stages build --bundle data/challenge --whitener-out data/donors/whitener_union_pruned.npz --whitener-id union_pruned --cd4-out data/external/cd4t/cloud-compact/cd4_compact --work <tmp> --name <x>`, then `vcc prep <x>.h5ad -g … --perts … -o <x>.vcc` | the same `8df4dadd…` / `f13ce40b…`, from the driver's own `SHIPPED_ENV`, with nothing exported. Re-run this whenever `SHIPPED_ENV` moves |
| **the arms guard fires** (240 row 4) | `uv run python experiments/240-driver-fix/rows240.py` | `REFUSING THE BUILD` on a config missing either arm, no raise on the shipped one, and `def_postbuild.py` **FAIL 16/20** on `rc_callset`'s config whichever `--incumbent` it is given |
| rule 10 | `uv run python experiments/199-def-runbook/rule10_199.py` | `RULE 10 PASS`, meta diff 0 fields of 50 |

---

## 5. Coverage-drop responses — decided now, not on the day

### 5.-3 THE SEALED HCASMC ASSET — keep it sealed, and use it as the FINAL CONFIRMATION set

126 acquired the official IGVF HCASMC source, byte-verified it
(`IGVFFI0206LUDV.tsv.gz`, 1,531,381,676 bytes, sha256 `b87452c2…`) and **sealed it**: no command may
decompress it, inspect its header, or derive any value-dependent statistic until a new hash-bound
authorization exists. 225 respected that boundary and flagged it as *"the one remaining unmeasured
source in the census."*

**Coordinator decision, 2026-09-04: it stays sealed, and this is a decision, not an omission.**

- **As a source it is worth roughly nothing.** 27/300 panel targets. On the measured rate of the
  87-row `igvf_cardiomyocyte` — **+0.00038 = 0.20× the shape-only gate** — 27 rows project to
  **≈ +0.0001**. It cannot clear any gate we use.
- **As a validation set it is the most valuable data we hold**, because it is the *only* dataset in
  this project that has never influenced a single decision. Every other table has been used to
  select, tune, or price something. **Opening it as a source destroys that property permanently and
  buys ≈ +0.0001.** That is a bad trade by three orders of magnitude.

**The intended use, and the one thing worth spending it on:** after the D/E/F stack is frozen and
before the final prize-eligible submission, run the frozen pipeline against HCASMC **once**, as a
clean out-of-sample confirmation that the method transfers to a line that never touched its
construction. **One shot, after freeze, no tuning against the result** — if it disagrees with our
internal estimates, that is information we could get no other way, and it arrives while there is
still time to act.

**Authorising it is a real decision with a cost** (it ends the asset's special status), so record the
freeze hash and the intent in writing before the run, exactly as 126 did on acquisition.

### 🔴 5.-2 SOURCE COVERAGE INVERTS IN OCTOBER (225, 2026-09-04). Read before re-selecting anything.

`experiments/225-source-sweep/` §8.3 projected every source we hold onto **both** §5 synthetic
draws, reusing 199's own bundles rather than a new sampling model, with three known-answer rows
against this section's numbers (`k562_gwps` 269→**150**, Orion **291**/**293**) — **exact on all
three**. Selected rows, pessimistic draw:

| source | today /300 | **PESSIMISTIC /300** | Δ |
|---|---:|---:|---:|
| `orion_hek293t` | 294 | **293** | −1 |
| `orion_hct116` | 290 | **291** | +1 |
| `k562_gwps` | 269 | **150** | **−119** |
| `dld1_lowmoi` | 67 | **73** | **+6** |
| `igvf_cardiomyocyte` | 87 | **36** | −51 |
| `gse295214_h1` *(shipped, in the H1 union)* | 64 | **34** | −30 |
| `h1_2025` *(shipped, in the H1 union)* | 20 | **3** | **−17** |
| `nadig_hepg2` | **0** | **33** | **+33** |
| `nadig_jurkat` | **0** | **33** | **+33** |
| `rpe1` | **0** | **32** | **+32** |
| `k562_essential` | **0** | **30** | **+30** |
| **`kolf2`** | **273** | **0** | **−273** |
| **`feng_ipsc`** | **182** | **0** | **−182** |

**CORRECTION 2026-09-26: the `kolf2` and `feng_ipsc` rows are artefacts of our LOCAL copies, not of the sources.**
- `data/external/kolf2` holds only the 273 A/B/C-panel rows 055 fetched by range request.
- The KOLF2.1J atlas has **10,167 perturbations** (`data/external/kolf2/official/meta.json`). It covers **251 / 181 / 230**
  of 300 on the realistic / pessimistic / efr panels. Feng targeted 7,226 genes.
- Consequence 1 below is wrong as stated: neither source goes to zero in October.
- Nothing shipped changes: both are NO-GO as donors on price (`kolf2` −0.00174).
- A kolf2 PLOLO fold IS buildable on October's panel, by range-fetching that panel's rows (~10 MB, 055's route). It is the
  fold that predicted the board for bulk-shape arms (417a). It over-reads H1-residual arms (iPSC ≈ H1), so it is no
  judge of the H1 dose.
- Rivals using KOLF2 keep it in October.

**Four consequences, all operational.**

1. **`kolf2` and `feng_ipsc` go to EXACTLY ZERO.** They are panel-**fitted** libraries — 273 of 273
   and 182 of 182 of their targets are on today's panel — and October's panel is **disjoint by
   construction**. They rank #3 and #5 by today's overlap and **#28 and #29 in October.** Both are
   already NO-GO on price so nothing shipped is affected, but **never rank a source by today's
   overlap again.**
2. **`k562_gwps` loses 44% of its coverage** (269→150). It is one of our three production donors.
   §5's single-carrier response tables assume this; they are now measured, not assumed.
3. **Our shipped H1 union degrades hard**: `gse295214_h1` 64→34 and `h1_2025` **20→3**. `h1union`
   scored **+0.0069114** today (rank 3) — **do not assume that carries to D/E/F.** Re-price the H1
   dose `g` on the new panel as a first-class step on Oct 22, not as an afterthought.
   **`dld1_lowmoi` by contrast goes 67→73 — it is one of the few sources that IMPROVES.**
4. **Four sources we already hold are worth ZERO today and ~30 targets each in October**
   (`nadig_hepg2`, `nadig_jurkat`, `rpe1`, `k562_essential` — common-essential screens). **PLOLO
   cannot price them at all**, because it can only score today's panel and they carry none of it.
   **Build them into the source registry now and hold them ready**; their first measurable price is
   on the new panel. Not doing this means discovering them on the day.

**And the rule 225 derived: a BROAD stack degrades better than a well-matched narrow one.** Against
a disjoint panel, redundancy that looks wasteful today is October insurance and is nearly free —
225's rank-2 stack (all thirteen sources) is the better October arm even though its price today is
marginally below the four-source rank-1 stack. **Prefer breadth on the final entry.**

### 🔴 5.-1 THE REHEARSAL IS NOW STALE — it must be re-run before Oct 22

> 🟢 **448 (2026-09-25): the re-level is a LOOP on D/E/F, not one pass.** Rehearsed at q 1.75 with the stack
> H1 ×2 / CD4 ×1 on `ef-bundle-r`, builder `b6c3f131` (bound [0.7, 1.4]). See `experiments/448-def-rehearsal-q175/`.
> - The panel needs factors 1.21–1.29 at q 1.75. The old bound would have refused every context.
> - ntclevel's factor assumes ρ ∝ f². The render gives ρ ∝ f^e with e 1.89 / 1.90 / 2.30 (D/E/F), and F's local e near
>   f 1.25 is 2.70. **One pass put F 5.7 % past its cliff (ρ 1.057)**, and D/E came in 2 % short.
> - **The day-of loop:**
>   1. lvl0 = the un-levelled build (unit factors).
>   2. `ntclevel417.py <lvl0.h5ad> --bundle $B --contexts D,E,F --sd 0.01` gives lvl1's factors.
>   3. After each build n ≥ 1, run ntclevel on it, then
>      `uv run python experiments/448-def-rehearsal-q175/relevel_iter.py <lvl0.ntclevel.json>:D=1,E=1,F=1 <lvl1.ntclevel.json>:<f1> …`
>      with every build so far. It takes a per-context secant step (bracketing when it can) and prints the next
>      `VCC_BULK_CTX_SCALE`, or CONVERGED.
>   4. Stop when every context is within 1 % of ρ* = 0.9793 in ρ and none is above 1.
> - Budget 3–4 builds (~25 min each). Never ship a build with any ρ > 1. If time runs out, ship the last build whose
>   contexts are all below their cliffs.
> - In the rehearsal the secant retrodicted lvl2's D and E to ±0.0001 in ρ. On lvl3 one bracketing step put the
>   expansive context F at ρ 0.9784 (target 0.9793), CONVERGED. That is two levelled builds for D/E, three for F.
>
> 🟢 **432 (2026-09-23) rehearsed the September stack end to end on synthetic D/E/F.** The stack is 433h's: stack
> LOAD, sharpening q 1.25, the 410 key re-fitted for the panel, coarse MDE, zero repair, re-level, plus zs. See
> `experiments/432-def-rehearsal-new-knobs/RESULTS.md`: **§8 is the day-of command sequence (~58–68 min, passes
> `bash -n`), §6 the traps, §9 the new `SHIPPED_ENV`.**
> - The builder was restored byte-identical after all 12 calls.
> - The three negative controls refuse within 10 s: an A/B/C dump, the A/B/C key, the A/B/C re-level string.
> - The traps to fix BEFORE Oct 22 (proposals in `432/proposals/`, not applied; the chains in flight pin the current
>   driver sha):
>   - **T6:** 417's CD4 ×1.5 dump recipe (`VCC_CD4_TARGET_F=0.51`) silently does nothing on D/E/F. Use
>     `VCC_CD4_TARGET_G`.
>   - **T2/T3:** a re-fitted key needs a builder source edit, and `build_weights410.py` overwrites the pinned A/B/C key.
>     Use `key432.py`.
>   - **T4/T5:** `def_pipeline.py`'s `SHIPPED_ENV`, arm guard and postbuild reference predate September.
> - **The re-level is mandatory on D/E/F.** On the low-coverage panel the un-levelled stack puts D past its mse cliff
>   (ρ 1.071). Derive the string on the LOADed artifact.
>
> 🔴 **438 (2026-09-23): re-price the stack's WEIGHTS on the day's panel. CD4 ×1.5 turns negative on October-like panels.**
> - On `efp`/`efr` (levelled, 3 cancer folds):
>   - CD4 ×1.5 on top of H1 ×2 costs −0.0011 / −0.0024, on 0/3 folds;
>   - H1 ×2 holds (+0.0026 / +0.0003);
>   - q 1.5 and the LODO key keep their sign at about half their A/B/C value.
> - **Default for October: H1 ×2 / CD4 ×1**, the lowest-regret point. That is `VCC_H1_UNION_G=0.56`, with **no** CD4
>   override in 432 §8 step 6.
> - Then run the day-of re-pricer between 432 §8 steps 1 and 6 (~30 min, ≤ 8.8 GB):
>   `uv run python experiments/438-october-stack-repricing/dayof438.py --bundle $B --op $OP --cd4-root $CD4 --tag def_$N
>   --board-q <q>`. It prints the dump override to use.
> - Full text: `experiments/438-october-stack-repricing/RESULTS.md` §5. Keep `438/work/bankstats/` (committed).
>
> 🟢 **439 (2026-09-24): the installed driver runs end to end on an October-like panel with NOTHING exported.**
> - Panel: `ef-bundle-r`. Stages preflight → build → postbuild.
> - Every September knob lands from the defaults (arms guard 66/66); the day's key is fitted in 3 s and registered by
>   the pin.
> - Build 24 min, 20.2 GiB. MDE dose median 1.0; 0 zeroed genes; builder restored in `finally`.
> - **Pass `--incumbent-meta` = the shipped A/B/C arm's meta on the day.** Without it, postbuild reads 40/41, and the
>   one FAIL is the drift against the September-5 incumbent. See `experiments/439-def-driver-e2e/RESULTS.md`.
>
> 🟢 **437 (2026-09-23): 432's proposals are INSTALLED and proved.**
> - Driver `5463fe28`, pin script `10af47e7`, `scripts/def_pert_scale_key.py` `0e375b77`, builder `63dfd7d0` (T6 guard).
> - **Row 1 PASS:** the driver's defaults plus only the two day-of exports (LOAD, CTX) rebuild 433h byte for byte
>   (`5d9e2e38…`).
> - With the install, 432 §8's steps 2–4 become the driver's own `key` / `pin-artifacts --pert-scale` stages.
> - The pin stages now run with a bare env (T3b): day-of exports no longer break them.
> - `--restore-builder-after` restores in `finally`, and only this run's pin (T10).
> - T6 is refused by the builder: `VCC_CD4_TARGET_F` ≠ 0.34 under carrying sizing.
> - `SHIPPED_ENV` includes the coarse MDE pending 433a's board read. Remove it if the read fails, and re-run row 1
>   (`experiments/437-driver-row1/`).

199 rehearsed the whole path end to end on both synthetic panels on **2026-09-03**
(`experiments/199-def-runbook/work/run1/`). **That rehearsal predates the stack we now ship** and
therefore does not cover:

- **`VCC_H1_UNION_G` = 0.28** — a code path added by 222 on 2026-09-04 and now in our best entry
  (`h1union`, board **0.2118310**, rank 3). The rehearsal never executed it.
- the **coverage-aware sizing rule** 222 added to `scripts/make_submission_pm_cd4.py` (the carrying-
  row median), which is precisely the rule that behaves differently at the low coverage a new panel
  produces — **the case the rehearsal most needs to cover and currently does not.**
- **any source 225 adds.** Branch A fired on 2026-09-04, so the stack is expected to grow.

> 🟢 **The DRIVER half of this staleness is closed (240, 2026-09-04).** `scripts/def_pipeline.py`
> now defaults to the shipped configuration (`SHIPPED_ENV`), takes an operator override instead of
> silently discarding it, asserts the arms off the finished meta, and `def_postbuild.py` diffs
> against `h1union` and checks the five arms by value. Verified on the A/B/C bundle: the driver's
> defaults rebuild the shipped package byte-for-byte. **What is still stale is the PANEL half** —
> the H1-union path and the carrying-row median have still never run against a low-coverage
> synthetic D/E/F panel. That is what the re-run below is for, and 240 does not substitute for it.

**Do not re-run the rehearsal yet — re-run it once the stack is final**, and no later than
**2026-10-15** to leave a week of margin. Re-running it now would only have to be redone. The
rehearsal is cheap relative to the risk: 199's run is fully scripted through
`scripts/def_pipeline.py --dry-run`.

**The specific thing to verify on the re-run**, beyond 199's existing rows: that the H1-union path
and the carrying-row median both behave at the **pessimistic** draw's coverage (§5's table:
`panel_support` 141/167, CD4 union 174), not merely at A/B/C's ~92%. That is where the over-dose
factor reaches **~1.24×** (§5.0) and where the old all-rows median silently returned zero.

### 5.0 The over-dose factor is now MEASURED (222, 2026-09-04)

The residual sizing over-doses in proportion to missing coverage. Measured on four PLOLO folds:
**1.012-1.017x at today's ~92% coverage** (negligible) and **1.237-1.245x at the pessimistic q0.14
draw** -- tightly clustered across folds, so treat ~**1.24x** as the planning number.

Sec 5's pessimistic draw puts CD4 at 174/300 (58%), which is **below** every coverage in that table,
so 1.24x is a floor for that scenario, not a central estimate. **On the day: compute the carrying-row
median, never the all-rows median** -- the harness mirror at
`experiments/162-retrodicting-validator/plolo.py:363-369` silently returns `f_per_unit = 0.0` below
50% coverage, so a real source reads exactly zero and looks dead (223). The builder raises instead of
returning zero, so this is a **measurement** hazard, not a shipping one -- but on Oct 22 we measure
before we ship.

> **Measured for every registered source, 2026-09-04 (235).** `work/carry235.json` and
> `work/essfree235.json` size all ten registry entries on **three** panels, both ways: A/B/C
> (**0 essential genes**), §5's pessimistic draw (**30**), and a third panel 235 drew to be
> **disjoint AND essential-free** (**0**), as §8.1g now requires. **Every source that covers a
> panel at all sizes above zero under the carrying-row rule — 37 of 37 (source, panel) cells —
> and the shipped all-rows rule returns EXACTLY 0.0 in all 37**; none of them clears 50% coverage
> on any panel. The coverage table reproduces 225 §8.3 on 8 of 8 rows *on the same bundle 225
> used*. A selected source that carries **nothing** is now a raise, not a silent zero.
>
> ⚠️ **And the October four look very different once the panel is essential-free** (235, and it
> agrees with the direction of §8.1g): `k562_essential` **30 → 0** (its entire overlap with the
> pessimistic bundle *was* the essential set — the new raise fires on it), `nadig_hepg2` 33 → 12,
> `nadig_jurkat` 33 → 13, `rpe1` 32 → 12; `dld1_lowmoi` goes **73 → 81** and `gse295214_h1`
> 34 → 22. **Re-rank October sources on an essential-free draw only.** 238 is recomputing §5.-2's
> table; 235's numbers here are one independent draw, not a replacement for it.
>
> 🔴 **436 (2026-09-23): expect ~0 essential-screen targets on D/E/F, not 12–13.** A/B/C carries none of the
> essential screens' 2,389 targets. That includes the 368 targets that are **not** on the K562-essential list.
> A random 300 from 238's essential-free pools would carry about 7 of those; P(0) ≈ 8·10⁻⁴. So Arc's filter
> also removes genes that are essential in RPE1, Jurkat or HepG2.
> - Point estimate: 0 covered D/E/F targets, 95 % upper bound 3.
> - 436's October-only DE lever (sign the MDE add by the covering donors' mean, with a capped dose) is therefore
>   worth +0.00000 at the point estimate. It is built only if the real panel has at least 20 covered targets.
> - **Oct 22 check (one command, no build):** `uv run python experiments/436-october-essential-signs/size436.py`,
>   then read its first row's `union` column. The production design is in `experiments/436-.../RESULTS.md` §5.


Measured on the two synthetic panels (both 300 targets disjoint from A/B/C):

| quantity | A/B/C (real) | realistic draw (`def-bundle-r`: K562 ∩ Orion ∩ CD4) | pessimistic draw (`def-bundle`: any donor) |
|---|---|---|---|
| K562-gwps coverage | 269 | 300 | **150** |
| Orion HCT116 / HEK293T | 290–294 | 300 / 300 | 291 / 293 |
| carriers 3 / 2 / 1 / 0 | — | 300 / 0 / 0 / 0 | 141 / 152 / 7 / 0 |
| agreement key defined (≥ 2 donors) | — | 100% | 97.7% |
| CD4 per condition, union | 283 / 287 / 284, 291 | 291 / 296 / 293, 300 | 173 / 173 / 168, 174 |
| `panel_support` (pm, cd4) | 248, 276 | 273, 286 | 141, **167** |
| panel ∩ essential list | 0 | 34 | 30 |

Responses:

- **3 covering donors** (the §2.1 expectation): nothing changes. Ship the table.
- **2 covering donors on a target**: the consensus is the mean of the two, the agreement key is
  defined, shrinkage runs. Nothing to do. On A/B/C this is already the case for ~30 targets.
- **1 covering donor**: the effect is that donor's alone, agreement key all-zero, alpha = 1 for the
  target. Acceptable up to ~10% of the panel (the pessimistic draw has 2.3%). If more than ~30
  targets are single-carrier, record it; do not re-weight donors by hand — there is no instrument
  to price it and equal weighting is the board-validated setting.
- **0 covering donors**: panel-mean direction at the median donor magnitude (`fallback="mean"`).
  Real prediction, not a zero row. If > 5 targets: that contradicts §2.1's selection criterion —
  re-check the gene symbols (case, aliases) before believing it.
- **CD4 gates fire (< 275 per condition or < 290 union)**: run step 3 with `--min-targets 0
  --min-union-targets 0` ONLY after confirming the mapping is right (`target_row_map.csv`); the
  residual then applies to fewer targets, the builder preserves the rest. Do not change the dose.
- **`gene_names.csv` changed** (step 6 shows its hash moved): **STOP (corrected 2026-09-29, 531 F6).** The operator, CD4
  banks and donor mappings key on symbols, but the pinned CD4 deflation basis (`VCC_CD4_DEFLATE_BASIS=V5_abc`) raises on a
  gene-axis mismatch, and LODO's same-path assert cannot see a reorder. Do not ship; do not refit the basis on the day.
  (The earlier "proceed" advice was wrong.)
- **panel ∩ essential list > 0**: the operator drops those basis rows (`dropped_in_panel` lists
  them). Expected ≈ 0 by §2.1. It is informational; no action.

---

## 6. Submission discipline in the dark — the checklist

1. **Submit the best-supported package EARLY**: on Oct 22 or 23, the incumbent recipe (§3) rebuilt
   on D/E/F, after steps 1–11 and the known-answer rows pass. Budget from bundle download to
   "scoring started": **≈ 30 min of compute + 15 min of upload + your own review ≈ 1.5 h**.
   Do not start it after 21:00 UTC.
2. Then **replace it only** for a package that is (a) strictly better on the local instruments
   (PLOLO `pds`, 096's claim `mse`, h1-test for emission-law arms at the 1.29 factor), (b) has a
   preregistered prediction written before the build, (c) passes the composed slot gate (rule 13),
   and (d) passes `def_postbuild.py` with EXACTLY the intended config differences vs the entry it
   replaces. A late unverified change is the one way to lose the round.
3. **Every re-submission costs ~1 h of wall clock before scoring starts** (30 min build + prep,
   15 min upload, review) and the last one must be uploaded with margin: **nothing new after
   Nov 5 20:00 UTC.** Slots: 2/day, one in flight at a time; the dark board still enforces both.
   - **Amended 2026-09-29: budget uploads in hours, not minutes.** On 09-29 the 2.77 GB upload of 446bd2 ran from
     03:08Z to ~10:40Z (148 MB after 3.6 h, then back to ~3.9 MB/s), and the same day's slot 1 took 1 h 49 min from
     upload start to published. One in flight at a time means a crawling upload also blocks the next slot. So: start the
     final entry's upload by **Nov 4 20:00 UTC**, check progress with the resumable session's byte range (a PUT with
     `Content-Range: bytes */<size>` to the pending `upload_url` in `~/.config/vcc/state.json` returns `range: bytes=0-N`),
     and treat Nov 5 as reserve only.
4. `-m` = the next name from `docs/MODEL-NAMES.md` (`scripts/next_model_name.py --claim <arm>`); the final Nov 5 entry uses the
   reserved final name there. `-d` EMPTY, never a number. Record name → config in RESULTS.md and the ledger.
5. Never `cell-eval2` while a build runs (25 GB + 19 GB does not fit with agents running), and
   never `vcc prep` alongside a build either (22–24 GiB measured).
6. Before each upload: `def_postbuild.py` PASS, `vcc prep --dry-run --contexts D,E,F` clean,
   sha256 of the `.h5ad` and `.vcc` written down, git commit of the pins.
7. If `vcc` is upgraded in the window: re-run §4's prep row on the rehearsal artifact first.
8. **Verify the dark-board assumption on the first D/E/F upload (added 2026-09-27).** Nobody has confirmed that
   final-phase scores are hidden from the submitting team: the CLI reference only says `vcc status` shows "scores, once
   published" (checked 2026-09-27; the rules page is a JS app and was not readable). After the first D/E/F entry
   finishes, run `vcc status <id>` and `vcc --json status <id>`, and capture it with `scripts/capture_board_entry.py`.
   - If members come back, the ~27 remaining final-phase slots are informative on the real panel. Stop and write a
     preregistered test plan (dose curves on the switches this runbook marks as uncertain) before the next upload. Rule
     §6.2 still applies: the last entry counts, so the final upload must be the best-supported package.
     **The plan is pre-written (2026-09-28): `docs/FINAL-PHASE-FEEDBACK-PLAN.md`** (order of reads, bins, and the rule
     that the last entry is the best scored package).
   - **Feedback-use policy (H-631): DECIDED. MAXIMUM-INFORMATION use, fully documented.** Ahmet's words, 2026-10-02:
     - ~12:23Z, [537757]'s session, in direct reply to the restriction question (which named the ~0.01 cost and the
       disqualification risk), verbatim: "i dont find that an issue. what we should do is do what helps us learn most. we
       should document and make use of all info. its not gonna disqualify."
     - ~12:40Z, [7772dc]'s session, on 637 §4's decode of D/E/F's own module programme, verbatim: "we can do the risky
       thing".

     His earlier line in [4d0207]'s session ("for 1., yeah we can do that ...") was ambiguous; the two answers above
     supersede it. Nobody re-asks him.
     1. If members come back, design the D/E/F reads for the most information per slot and use all of it. 637's factorial
        reader stays the default for switch reads (after its A/B/C validation PASS). The fused bulk-switch + module-decode
        design (637 PROPOSAL §4) is APPROVED in principle. It is adopted once it meets the acceptance criteria in
        FINAL-PHASE-FEEDBACK-PLAN §6 (cadence, re-roll noise), with 625's A/B/C run as its end-to-end rehearsal.
     2. Document every read: each D/E/F entry goes into `docs/BOARD-LEDGER.md` with what it was designed to read, the read,
        and how it was used. The finalist write-up discloses this.
     3. Unchanged: R0 first, the last entry is the best-supported scored package, `-m` from the queue, `-d` empty
        (`docs/FINAL-PHASE-FEEDBACK-PLAN.md` §1).
     The council's earlier recommendation (637 PROPOSAL §2: preregistered switches only, no truth-bit reads) is superseded
     and kept there for the record.
   - If nothing comes back, the dark-board plan above stands unchanged.
   - **Factorial reads (637, 2026-10-02; reader accepted, validation PENDING).** If members come back, each entry after
     R0 may carry ONE bulk switch AND ONE prior-only DE switch:
     - the bulk switch is a row 2/3/4, read from pds + mse;
     - the DE switch is a row 6/7, read from the DE members after a −0.0003 ± 0.0002 bulk-leak allowance.

     Command: `uv run python experiments/637-final-phase-factorial/readdef_factorial637.py --factorial BULK_ROW DE_ROW ENTRY BASE`.
     - A verdict counts only if it is the same at every leak in the band; otherwise it reads `inconclusive_leak_band`
       and nothing is adopted.
     - Refused: row 5 (moves both channels), rows 8–9 (ballast-moving), rows 10–12 (combinations).
     - Use it only after its A/B/C validation arm reads PASS (`--validate FACT X_ALONE Y_ALONE BASE`, exit 0; the
       all-switches arm, 10-08 … 10-10). Until then, one switch per entry.

---

## 7. If something breaks

| symptom | cause | fix |
|---|---|---|
| `bundle manifest declares contexts [...] but scripts/frozen_inputs.json pins [...]` | step 6 not run | run step 6 |
| `frozen model input hash drift for cd4/...` or `whitener/...` | step 7 not run, or run against different files than the build uses | re-run step 7 with the exact paths the build gets; `git diff` |
| `unknown --whitener-id` | id not registered | step 7 registers it; check the id string |
| `PerturbMatch replacement support changed: N != 248` | `panel_support` still A/B/C's | step 6 sets it to null; step 8 pins it |
| `RSS ... exceeded the 10.0 GiB budget at covariance` (operator fit) | memmap pages counted | `--max-rss-gib 14` as in step 5 |
| builder aborts at the CD4 residual with a zero-norm error | A/B/C banks on a new panel | step 3 + step 7 |
| `vcc prep: Unknown context label(s) in 'context': D, E, F` | missing `--contexts D,E,F` | add it |
| operator sha equals `e2f65ecf…` after step 5, or `axis sha256` equals `edf6ba88…` | the bundle on disk is still A/B/C | check `data/challenge/manifest.json` |
| `REFUSING TO BUILD: the whitening operator was fitted for a DIFFERENT PANEL` | the build is pointed at an operator fitted for another panel — the Oct-22 failure, caught (235) | re-run step 5 against THIS bundle, re-pin with step 7, and check the `--whitener`/`--whitener-id` pair the build actually got. Never widen the guard |
| `REFUSING TO BUILD: this operator carries no axis digest` | an operator predating 235 that is not in `LEGACY_OPERATOR_AXIS_SHA256` (e.g. `whitener_union_full.npz`) | re-fit it with `build_union_whitener.py`, which writes `axis_sha256`. Do not hand-add a digest you have not independently derived |
| `union source <name> carries 0 of 300 panel targets` | a `VCC_UNION_SOURCES` member covers none of this panel and would deliver exactly nothing (225: `kolf2` 273 → 0, `feng_ipsc` 182 → 0 in October) | drop it from the selection, or select a source that covers this panel — the coverage table is `experiments/235-october-hardening/work/carry235.json` |
| the driver refuses a heavy stage | another builder running or < 26 GiB free+inactive | wait; never run two builds |
| source download needed | it is not: the source is pinned locally (§1). If the file is gone: `https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad`, 16.8 GB, then verify `c355f535…` | |
| `FAIL [X] control library q75/q25 <= 3.0 (554: …)` | wide control libraries: the emitter's steering dose would collapse | the top block's 554 row: ≤ 3.5 everywhere → `DEF_DAY_ALLOW_WIDE_LIBS=1` in a fresh dir; above → STOP |
| `FAIL [X] control library median within 15% of A/B/C's 20,030 (641: …)` | D/E/F released at another depth | the top block's 641 row |
| `[sameline] X: … min top-1 correlation … < 0.65; context unchanged (590)` | the line is probably absent from CCLE | expected; read the top 5. Only 642 can switch a non-CCLE line on |
| `postbuild FAILED …` with a `(554)` row in postbuild.json | a context's DE design was delivered at low dose | the top block's 554 row |
| `STOP: 554: the final build wrote no rc_audit dump` | the build ran without `VCC_RC_AUDIT_DUMP` (a pre-554 `def_day.sh`) | check `def_day.sh`'s sha against the post-bundle table (594), then rerun in a fresh dir |
| `642: … is not the pinned …` or `642: reference … lists Flex (or no) sources` | the HEK293 reference file changed | `git checkout -- data/identity_refs/hek293_nonflex_ref642.npz`; never rebuild it on the day |
| `[sameline] X: the non-CCLE reference called HEK293; …; context unchanged (642)` | the knob or the decision text is missing | by design; see the HEK293 path |
| `642: DEF_DAY_ALLOW_NONCCLE_REF='…' is not one of ['hek293']` | a typo in the knob | fix it |
| `STOP: 644: DEF_DAY_SAMELINE_SKIP names X, which is not a manifest context (…)` | a typo, or a label not in the manifest | use the manifest's context names |

---

## 8. Files

- `scripts/def_pipeline.py` — the driver (§2). `--dry-run` prints every command without running.
- `scripts/def_preflight.py`, `def_identity.py`, `def_census.py`, `def_postbuild.py`,
  `def_pin_panel_artifacts.py` — the stages. All new in 199; none is imported by production code.
- `experiments/199-def-runbook/RESULTS.md` — the rehearsal record, timings and the rule-10 proof.
- `experiments/199-def-runbook/work/run1/` — the pessimistic rehearsal (logs, timings, meta, .h5ad).
- `experiments/199-def-runbook/work/run-r/` — the realistic bundle's cheap stages.
- 2026-10-02 additions (644):
  - `experiments/554-dose-guard/` (spread and dose guard);
  - `experiments/590-identity-validation/` (identity floor);
  - `experiments/641-depth-row/`;
  - `experiments/642-hek293-contingency/` and `data/identity_refs/` (the HEK293 contingency);
  - `experiments/637-final-phase-factorial/` (factorial reader);
  - `experiments/644-runbook-amendments/`: `verify644.py` re-runs every command these rows cite against scratch
    post-bundle trees (27/27 on 10-02), and `DECISION-CARD.md` is the one-page day-of card.

## 410 / 410b additions (2026-09-22) — per-target bulk key and ballast cells
- **The per-target weight artifact is panel-specific AND axis-specific.** Rebuild it on Oct 22 AFTER the operator refit:
  `experiments/410-pertscale-revisit/build_weights410.py LO HI` with `OP_PATH` pointed at the refit
  `whitener_union_pruned_def.npz` and `wcols` taken from that operator's own axis (not `fold-kolf2.npz`, which is the A/B/C
  axis). The script now REFUSES if any panel gene sits on the operator axis (the A/B/C axis contains 0 A/B/C panel genes
  and would contain many D/E/F panel genes). Then re-pin the new artifact's sha in `PERT_SCALE_ARTIFACTS` and re-run the
  knob tests (`test_pert_scale.py`, `test_ballast410.py`).
- **Ballast (`VCC_PERT_SCALE_BALLAST=1`) is panel-agnostic** (no artifact), but it assumes the incumbent emitter path
  (`--emitter rc_callset` + `VCC_DE_SIGN_PRIOR`); the builder refuses other paths.

## 417 additions (2026-09-22): bulk sharpening, carrier gamma, per-context bulk scale, bulk dump/load
- **`VCC_BULK_SHARPEN_Q` and `VCC_CARRIER_GAMMA` need no artifact.** Sharpening acts on each context's own 10-CPM render
  gate and restores each row's non-panel norm. It protects the panel columns, so the NEW panel's columns are protected
  on D/E/F automatically, because the mask comes from `spec.perts`. The carrier gamma re-runs whitening with the
  day's refit operator.
- **`VCC_BULK_CTX_SCALE`: the numbers do not carry, but the method does** (corrected 2026-09-23 01:30Z).
  - The cliff for context c sits at bulk scale √(claim_c / B_c).
  - B_c is our artifact's across-perturbation spread; `decomp414` reads it from the built h5ad.
  - claim_c = Σ_p w_p E[min(Ĉ_p, C_p)]. Ĉ_p is our blocks' jackknife. C_p is proxied from the context's own **control
    cells** (the 46-NTC-group jackknife: 096's `mech_a12_audit.cp_proxies()["ntc"]`, `Cp46` in 180/187), times a
    board-fitted factor f_c = 1.00 / 1.00 / 0.9825 on A/B/C.
  - So for D/E/F:
    - compute Cp46 from the D/E/F controls with f = 1;
    - compute A_c/B_c/Ĉ_p from the D/E/F build;
    - set each context's scale at the EV-optimal point below its proxy cliff, with a wider margin. Take
      sd(ln α) ≈ 0.01 for the unfitted f_c, against A/B/C's 0.0004–0.0025.
    - sumD (the truth's effect energy) is unknown, but it only sizes the gain. It does not move the cliff.
  - This needs a small tool, `ctxlevel417.py` driven by NTC-proxy claims instead of 187's fitted model. Build and
    rehearse it on the synthetic bundle before Oct 15.
  - **Verified 2026-09-23 01:38Z:** 187's fitted claims equal the saturated proxy f_c · Cp46_c · Σw_c to 1.0000 / 0.9927
    / 1.0000 (A/B/C). So claim_c ≈ Cp46_c · Σw_c is the D/E/F formula.
    - Cp46 = the mean over the 46 NTC constructs of cell_eval2's `jackknife_correction` on each construct's 400 cells.
      `experiments/180-def-readiness/cpproxy180.py::per_construct_jackknife` computes it from any bundle. A/B/C
      reference values: 14.4398 / 14.8057 / 16.4397.
    - Σw_c = the scorer's row weights (`decomp414`).
    - Ballast's jackknife boost keeps our Ĉ_p ≥ C_p, so the claim stays saturated.
  - **The tool exists and passed its known answer (2026-09-23 02:30Z):**
    `uv run python experiments/417-bulk-sharpen/ntclevel417.py <artifact.h5ad> --bundle <bundle> --contexts D,E,F --sd 0.01`.
    - On 417i's A/B/C artifact with f = (1, 1, 0.9825), Cp46 reproduced 180's values exactly. Its cliffs matched 187's to
      +0.00 % / +0.37 % / +0.00 %.
    - Runtime 20 s, peak 12 GB (it reads the three control files).
    - It prints the `VCC_BULK_CTX_SCALE` string relative to the artifact it reads. For D/E/F, run it on the day's
      un-levelled build, then rebuild with that string (one extra build).
- **The 415 LODO key and 410 floor-0.5 artifacts are panel-specific** (same rule as above). Rebuild them with the day's
  operator if they are in the shipped config.
- **`VCC_BULK_DUMP` / `VCC_BULK_LOAD`**: a dump is tied to the bundle's perts/genes. LOAD refuses on mismatch, so an
  A/B/C dump cannot leak into a D/E/F build.
- **Pre-upload check for any bulk arm:** `experiments/414-bulk-reprice/cliffcheck414.py <artifact.h5ad>`. On A/B/C it
  reads render-space energy against 410b. On D/E/F there is no 410b artifact, so compare against the day's incumbent
  build instead: an arm that changes render-space energy by more than a few percent in any context is an amplitude
  arm. 410d moved context C ×2.16 while its effect-space energy was held fixed.

## 🔴 OPEN (2026-09-24 22:20Z): the re-level factor bound [0.8, 1.2] binds on efr-like panels at q ≥ 1.5
- `VCC_BULK_CTX_SCALE` is refused outside [0.8, 1.2] (builder import-time check).
- 438's own level file (`experiments/438-october-stack-repricing/work/efr/level438.json`) needs these factors on the
  essential-free, realistic-coverage draw:
  - stack q 1.5: A 1.145 / B 1.172 / **C 1.208**;
  - stack q 1.75: **A 1.212 / B 1.213 / C 1.261**;
  - stack q 2: 1.29 / 1.26 / 1.32.
- The A/B/C factors fit: q 1.75 is 1.10 / 1.15 / 1.17 (446b), `real` 1.09 / 1.14 / 1.18, `efp` 1.00 / 1.04 / 1.07.
- 438 §3 says "inside the builder's [0.8, 1.2]" for efr. That holds for q 1.25 only. Its levelled efr prices assumed
  uncapped factors.
- **Decide before Oct 22:**
  - (a) widen the bound (e.g. to [0.7, 1.4]) with a known-answer row: unset and in-range values bit-identical, and an
    A/B/C build at factor 1.3 on one context checked for ballast realisation (KMAX capping) and render integrity; or
  - (b) a day-of fallback rule: if `ntclevel417.py` asks for any factor > 1.2 at the shipped q, rebuild at the lower q
    whose factors fit, or cap at 1.2 and accept the shortfall. Price the shortfall first: 0.31 raw mse per unit scale
    below the cliff.
- Owner: the coordinator after the 09-26 read of 446b (q 1.75), which decides the shipped q.
- **(a) DONE 2026-09-25 01:20Z (447):** the builder bound is now [0.7, 1.4] (builder `b6c3f131`, 56/0/1 suite).
- **Rehearsed 2026-09-25 02:31Z (448):** factors 1.21–1.29 build cleanly on `ef-bundle-r` at q 1.75, and the re-level
  becomes a secant loop (§5.-1, 448 entry). lvl3 converged (F ρ 0.9784). CLOSED.

## RESOLVED (2026-09-26, Ahmet): the same-line rows and the ballast cells are allowed
Both are within the rules; ship them as designed. No help-desk question. Describe them plainly in the finalists'
write-up (same-line rows as a similarity-conditioned donor inside the fitted transfer model; ballast cells as F19).

## 🔴 OPEN (2026-09-25 03:10Z): the SAME-LINE mode (§3 row 220) is being built. It may be the largest October lever
- **Evidence (450, H1, analytic).** Pan 2025's H1 hESC screen is a different lab and chemistry (10x 5′) from the 2025
  VCC H1 Flex truth. Scored against h1-test:
  - block-sign programme 0.625 vs CLEAN4 0.581 (oracle 0.653);
  - per-target gate signs 0.764 vs our 0.551 (26/26 targets);
  - per-target centred cosine 0.260 vs 0.069 (30/30 targets).
  A D/E/F line with a public bank would carry its own biology into our prediction.
- **Likely lines.** The six lines come from different tissues and blood is taken by Jurkat. Banks we hold (451
  registry): HCT116 (Orion, genome-wide), K562 (GWPS), DLD-1, HepG2 / Jurkat (essential: programme only), A549, HT29,
  MCF7, BXPC3, HAP1, Calu-3, THP-1.
- **Pieces (451):**
  - `scripts/def_sameline.py`: identity JSON → `<work>/sameline.npz`. Its known answer passes: A = Jurkat reproduces
    396's `psign_A` bit for bit.
  - Builder knob `VCC_SAMELINE="<npz>:<sha>:<w>[:pt]"`: patch `451/patch_builder451.py`, unit tests
    `451/test_sameline451.py`.
  - Chained known-answer builds (`451/chain451.sh`): unset must rebuild 446b; programme-only must rebuild 446bja.
- **Board test of the programme part:** 446bja (A's block prior = Nadig Jurkat). The per-target part has no A/B/C
  instrument: Nadig covers none of the essential-free panel. Its dose (w, per-target prior) comes from 450b/450c on H1.
- **The emitter-level H1 price (450 §9).** With 30/85 targets covered, rows w 1 + per-target prior + CLEAN4 default
  reads **+0.084 overall (board units)**: mse −0.17 raw, nmae −0.074, reach +0.070, pds +0.011, fid +0.027.
  - Raw rows beat centred rows.
  - The programme default costs reach on uncovered targets.
  - Pending: w 2 / 4 (450d) and the board read of the programme (446bja).
- **Day-of rule (draft, pending 450d and 446bja):**
  1. After the driver's `identity` stage:
     `uv run python scripts/def_sameline.py --bundle $B --identity <work>/identity.json --out <work>/sameline.npz`.
     A context switches on only if its call is unanimous and its line is registered.
  2. If the artifact names any context: export
     `VCC_SAMELINE=<work>/sameline.npz:<sha256>:<w>:pt:noprog` (w = 450d's winner, provisionally 2) for EVERY
     levelled build. The dump does not need it: the knob acts in the writer.
     - Add the programme (drop `:noprog`) only if 446bja's board read is positive.
     - The blend changes the across-perturbation spread, so the re-level loop (§5.-1) runs with the knob set from lvl0.
  3. Audits as usual. `config.sameline.contexts` in the meta says what each context got.
  4. No registered unanimous line: nothing changes.
- **Still to do before Oct 22:**
  - apply the builder patch (`451/chain451.sh`, after its two known-answer builds);
  - the forced-identity mechanics rehearsal (`451/run451r.sh`);
  - optionally a driver stage (the manual export above works without one).

## 426 / 428 additions (2026-09-23): the DE-channel magnitude floor and the zero-gene repair
- **`VCC_DE_MDE="c,cap,min_count"` (426) needs no artifact.**
  - The add is min(c/sqrt(count), cap), with count = the context's own control mean CPM × median control library / 1e6
    (from `rc_control` and `libraries`, already built per context).
  - It runs inside the K-class prior mode, so it needs `VCC_DE_SIGN_PRIOR`.
  - On A/B/C the delivered dose (`mde_lambda`, rc diagnostics in the meta) had median 0.84 because the prior's
    class-capacity check binds. Read that diagnostic on D/E/F: a median far below 0.8 means the arm delivered little.
  - Board validation: 426a on 09-24 (see BOARD-LEDGER); ship it on D/E/F only if its read landed in the "transfers"
    bin.
- **`VCC_DE_ZERO_REPAIR=-5` (428) needs no artifact and cannot hurt.** It repairs gated genes whose block total is 0
  (delivered lfc ≈ −35). Ship it always.
  - On the day, run `experiments/428-zero-gene-audit/audit428.py <artifact.h5ad> --bundle <bundle>` on the built D/E/F
    artifact: it must report 0 zeroed gated genes.
  - The audit also prints the extreme-positive tail.

## PENDING switches: placeholders (644, 2026-10-02). Not part of the day until a read decides them.
If adopted, each one needs:
- `october_recipe.py decide <switch> <value> --evidence …`;
- a `def_day.sh` pass-through;
- an inert KA and a rehearsal.

Only then does its row move into the top block.

- **612 KGID (key-gated identity hook): RETIRED 2026-10-02** (652 F8). The ungated v3 identity hook `idall3` replaces it; its board
  read vs 496a decides the identity hook (612/PREREGISTRATION-draft-IDALL3.md, cf3f8004). The KEYGATE knob and its 612dd
  pass-through stay applied, inert. Text kept for the record:
- **612, the key-gated uncalled-gene hook** (council [21dd4c]).
  - Builder knob `VCC_DE_MDE_SIGN_KEYGATE=<q>`; needs `VCC_DE_MDE_SIGN`.
  - A robustness lever for unseen lines: worst-case regret 0.0020 against 0.0115 for the plain hook, ~0 EV on
    CLEAN4-friendly lines. **VOID:** these regret numbers came from a call-set-blind fold model that failed re-validation
    (612 CORRECTION 2; NEXT-SLOT-DECISION 10-02 12:58Z). On H1's real emitter, KGID ties the plain hook at f 0.35–0.5, and
    KGID now needs a board WIN over 496a (A1's FAIL branch).
  - Read: a 10-08+ arm. The def_day pass-through `DEF_DAY_MDE_SIGN_KEYGATE` is 612dd, applied 10-06.
- **630, rank-neutral UP magnitude** (council [d4cd5f]).
  - Builder knob `VCC_DE_UPRN=<s|b>_k<K×100>`, applied after zs.
  - H1: nmae −0.0038 … −0.0066 raw, with the other five members exact and 0 call flips.
  - Read: 10-08+. No pass-through yet.
- **639 / 604, the identity-axis sign prior** (council; [4d0207] does the plumbing).
  - Delivered through the same-line path (`VCC_SAMELINE=<npz>:<sha>:0`, w 0, no rows) and computed on the day from
    D/E/F's own controls. Ship point λ 0.5.
  - Hard gate before any slot: per-context median prior_lambda ≥ 0.9 (554's dump).
  - Day-of path (604, applied 10-06; no further `def_day.sh` edit): block prior `DEF_DAY_PROGRAMME_SOURCE=identity:0.5`;
    identity hook `DEF_DAY_MDE_SIGN=hookbank:identity:0` + `DEF_DAY_MDE_SIGN_KEYGATE=0.5` (612dd).
