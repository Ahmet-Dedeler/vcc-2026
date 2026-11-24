# 644 — DEF-RUNBOOK amendments for everything added since 09-30, plus the Oct 22 decision card (DRAFT)

Date 2026-10-02, session [537757], at the coordinator's request (due 10-08). Shared docs are not edited; the coordinator
merges. Class: operations.

## Deliverables

| file | what |
|---|---|
| `DEF-RUNBOOK.644.diff` | unified diff against `docs/DEF-RUNBOOK.md` (base sha256 `cfa1873d…`, commit 9ce09a9f); +109 lines; no line-77 or H-631 text (both are the coordinator's, done); regenerate with `make_amendment644.py` if the runbook moves again |
| `DEF-RUNBOOK.644.md` | the amended runbook in full (sha256 `388f2e0f…`) |
| `DECISION-CARD.md` | the one-page day-of card: every stop and flag in day order, what it means, the exact action |
| `make_amendment644.py` | regenerates both from the live runbook (anchored inserts; refuses on an anchor miss) |
| `verify644.py` | runs every command the new rows cite against scratch post-bundle trees; `work/verify644.{json,md}` |
| `patch_skip644.py` + `test_skip644.py` | `DEF_DAY_SAMELINE_SKIP=<ctx[,ctx]>` v2, UNAPPLIED: `def_sameline --skip-contexts` + the `def_day.sh` pass-through; anchored after 642; the LAST 10-06 step (gap 1's fix, approved by the coordinator) |

## What the amendment adds

1. **A top block** above 531's 09-29 notes.
   - The patch order: 538 → 554 → 590 on 10-05, 641 at step 6b, 598 → 612 → 642 on 10-06.
   - A guard table. For each guard: the row that trips, what it means, when an override is legitimate, and what to record.
     - 554's library spread: `DEF_DAY_ALLOW_WIDE_LIBS` only at ≤ 3.5 in every context.
     - 641's depth row: `DEF_DAY_ALLOW_DEPTH` only after pricing a common ratio.
     - 590's identity floor: no override.
     - 554's dose guard: `DEF_DAY_ALLOW_LOW_DOSE` only for an MDE median of 0.85–0.90 with prior ≥ 0.9. Postbuild needs
       the pinned builder, so this is a fresh-dir rerun of about 2 h.
     - 642's HEK293 contingency.
   - The Oct 22 identity reading rule: log the call, n/8, min_top1 and top 5; hand-read any unanimous call at 0.60–0.70;
     never `--force` without a written reason.
   - **How to switch same-line off:** `def_day.sh` arg 5 = `0` AND arg 6 = `noprog`.
     - w 0 alone still applies `pt`'s per-target signs: the builder skips the blend at w 0 but not the signs (~5693).
     - It acts on every matched context.
   - The HEK293 path, step by step: a read-only diagnostic command, the 8/8 ≥ 0.75 bar, the coordinator's written
     decision, the day's env, and the pinned reference.
   - 598's tercilec knob.
2. **§6.8 additions.**
   - 637's factorial reads: the command, the allowed and refused rows, the leak-band rule, and "only after its A/B/C
     validation reads PASS".
   - No H-631 text: §6 item 8 is the coordinator's paragraph and records both readings until Ahmet confirms.
   - The obsolete "get Ahmet's go-ahead for downloads" paragraph is replaced by his standing OK.
3. **§7 rows** for the 8 new failure messages.
4. **§8 files.**
5. **A PENDING section** for 612 (`VCC_DE_MDE_SIGN_KEYGATE`), 630 (`VCC_DE_UPRN`) and 639/604 (`VCC_SAMELINE` w 0
   programme default). None is on the day until its read decides it; each then needs a decide, a pass-through, a KA and
   a rehearsal.

## Verification: 27/27 PASS (`verify644.py`, 37 s, 3.6 GB; trees T0 = 538+554+590, T1 = +641, T2 = +642, T3 = +644 skip v2)

| group | checks |
|---|---|
| patches | 641 `--check` on T0; 642 `--check` on T1; 642 refuses a second apply; 644 skip `--check` on T2, then `bash -n` and its start-log line on T3 |
| `def_day.sh` | `bash -n`; logs all four override knobs at start (WIDE_LIBS, LOW_DOSE, DEPTH, NONCCLE); no-arg run stops at `${1:?bundle}` |
| preflight | `--verify` PASS with the spread and depth rows silent on A/B/C; `--help` |
| identity | `--verify` prints 590's lowest top-1 correlation; with `DEF_DAY_ALLOW_NONCCLE_REF=hek293`, `[642]` lines for A/B/C, none called; an unknown reference is refused |
| same-line | `--help` shows `--min-top1`; run on the stage-2 identity file: A on via Nadig Jurkat (590 gate passed), B/C none |
| postbuild | `--help` shows `--rc-audit`; `--rc-audit` on 433a's healthy dump adds no 554 row |
| recipe / capture / readers | `october_recipe.py show` and `check --allow-pending`; `capture_board_entry.py --help`; 531 `readdef.py --help`; 637 reader usage and its 27-check suite |

**Not verifiable yet** (listed in `work/verify644.json`):
- 598's `DEF_DAY_SHARPEN_PT=tercilec:…` (lands 10-06);
- 612's `VCC_DE_MDE_SIGN_KEYGATE` and 630's `VCC_DE_UPRN` (builder knobs, not applied, no pass-through);
- 639/604's day-of `def_day.sh` patch (not written);
- the D/E/F bundle itself, and a real `def_day.sh` run on a bundle (~2 h; the coordinator's rehearsals cover it).
  Correction: post-538, `def_day.sh` derives REPO from its own location; only the pre-538 file hard-coded the live repo;
- `readdef637.diff`'s dispatch (unapplied).

## Gaps found while writing it (for the coordinator)

1. **No per-context same-line off switch: CLOSED by `patch_skip644.py` v2 (UNAPPLIED; the last 10-06 step, after 606's
   chain 641 → 598 → 612 → 612dd → 630 → 642 → 604).**
   - **v1 was replaced** (11716c46: `--contexts` in `def_day.sh` only). With a programme source (604's identity prior on
     October, or a bank), v1 would also have removed the skipped context's programme.
   - **v2 adds `def_sameline.py --skip-contexts`.** A skipped context is recorded "none (644: …)" before any bank, 642 or
     590 logic, so it is not `matched`, and everything downstream treats it as unmatched:
     - without a source it is absent from the artifact;
     - with a source (a bank, or 604's identity branch, which keys on `matched` too) it gets the source's programme at
       w 0, with no rows and no signs.
   - **The `def_day.sh` side:** it checks the names against the manifest (else STOP), passes `--skip-contexts`, logs at
     start and records the knob in `def_day.final.env`.
   - **v3 (2026-10-02 ~14:00Z, Grok review 607's MINOR + NOTE):**
     - names are whitespace-stripped in both files (`"E, F"` == `"E,F"`; v2 failed closed on `" F"`);
     - the skip is logged as "treated as unmatched (644)", not "context unchanged": a programme source still writes the
       context at w 0.
   - v3 patch sha256 `65fab4e9…`. On 606's current chain (… 630 → 630dd → 642 → 604 with 607's fixes):
     - pre-skip: def_sameline `9549b055…`, def_day `cb2c9e67…`;
     - after v3: def_sameline `5cd3ab07…`, def_day `ca84075d…` (as of 14:00Z). **Superseded for def_day:** 630dd moved
       twice since. The authoritative pin is the last PINNED block of `606-applyday-composition/RESULTS.md` (cfff2927):
       def_day `54fcf040`, def_sameline `5cd3ab07`. 604f anchors on 644 v3's `${SLSKIP:+--skip-contexts "$SLSKIP"}` line,
       and the chain plus 604f composes ALL PASS ([4d0207]);
     - v2 on the same chain gives def_sameline `f147f8bb…` (the coordinator's pin) and def_day `a4372191…`.
   - v2 was 4c0bb12e; v1 11716c46.
   - **`test_skip644.py` v3: 4/4 OK** on 606's chain trees (test598.canon_tree + 641 → … → 630dd → 642 → 604). 74 s, peak RSS **6.2 GB** (v3's run was announced to the coordinator first):
     loading orion_hct116 for the programme-source case put it above the 4 GB announce level, unannounced. No
     coordinator window was open; reported here.
     - **M:** `--check`, apply, second apply refused, a chain without 642 refused, an anchor miss.
     - **(a) unset:** on 538's stub harness, identical calls (argv + VCC_*; 15 and 17 records in two env
       configurations) and an identical `def_day.final.env`; no `--skip-contexts`. The real `def_sameline.py` without the
       flag writes byte-identical artifacts with no source and with `--programme-source orion_hct116`
       (v3 run: `15058635…` / `b8f52926…`).
     - **(b) SKIP:** the day passes `--skip-contexts E`, logs it and records it; `"E, F"` passes `E,F` (v3); `SKIP=Q`
       stops the day; `def_sameline --skip-contexts "B, C"` works (v3). The real
       `def_sameline.py` on an identity matching A/B/C (JURKAT / K562 / HCT116; covered 0 / 299 / 290):
       - no source: A and C byte-identical, B absent;
       - with a source: A and C byte-identical, and B in an unmatched context's form (w 0, no targets, rows or signs,
         `programme_B` == the source programme);
       - a 642 HEK293 call on B, with knob and decision set: the skip wins.
     - **(c) builder 4daae237, no build:**
       - its own `sameline_load` returns {A, C} for the no-source artifact, and `.get('B')` is None;
       - for the source artifact, B's record has no targets and w 0;
       - every same-line effect sits inside `if _slr is not None:` (5685–5732), whose else is the unmatched path, and the
         `pt` loop iterates the record's own (empty) targets.
2. **A tripped postbuild guard costs a full rerun.** Postbuild reads the pinned builder, and the failed day's exit trap
   restores it, so `DEF_DAY_ALLOW_LOW_DOSE` cannot be applied after the fact.
   - **Coordinator's call (10-02): no script.** The procedure is a fresh-dir rerun (~2 h); 14 final-phase days absorb it.
   - A pin-keeping resume path would be a new, unrehearsed path in the dark.
   - The failure row and the decision card both say "rerun in a fresh dir".

## Re-run after the 10-05 / 10-06 patches land (due 10-08)

`verify644.py` follows 606's P1 rule (2026-10-02 ~17:20Z):
- When it builds its trees, a patch whose `--check` says "already" is skipped. A "partially" result still fails.
- Its patch `--check` rows read "already on HEAD" once that patch is applied.
- Proof: on a fully patched tree, all six steps (538, 554, 590, 641, 642, 644) report "already". Today's HEAD is still
  27/27.
- So it runs as-is on 10-08 against the applied live tree. Before this fix it would have stopped at 538's applier.

```bash
uv run python experiments/644-runbook-amendments/make_amendment644.py   # re-anchor on the live runbook
uv run python experiments/644-runbook-amendments/verify644.py          # every row against the applied tree
```

After 10-06 the T-trees should be replaced by the real applied tree, or by 594's composition plus 598 / 612 / 642, and
the tercilec row added.

## The re-verify on the APPLIED tree (2026-10-06, after APPLY DAY v2 76d0b2ed / 5ee58ffc; due 10-08, done early)

`verify644.py`: **37/37 PASS** against the live tree. Its shas equal the coordinator's pins: def_day `54fcf040`,
def_sameline `5cd3ab07`, def_identity `7073406f`, def_preflight `17f0836f`, def_postbuild `557c2b52`.
- A1 / A2 / A4 (patch `--check` on earlier trees) read "already on HEAD": by design, nothing is left to check.
- Every day-of row passes on the applied files (def_day, preflight, identity incl. 642's stage 2, same-line incl. 590's
  floor, postbuild incl. 554's dose guard, recipe, capture, readers).
- **New post-apply rows (K)**, the 10-06 knobs that were unverifiable before:
  - K1: `october_recipe.py command` emits the decided day, `DEF_DAY_LODO=1 DEF_DAY_LODO_HI_PCT=96.64697522769474
    DEF_DAY_SHARPEN_PT=tercilec:2.5,2.0,1.75 DEF_DAY_DE_MAG=dn_k250_zs DEF_DAY_MDE_SIGN=clean4:0.05 bash scripts/def_day.sh
    … 1.75 0.5 noprog 4`;
  - K2: def_day.sh carries 598's tercilec, 612dd's KEYGATE, 630dd's DE_UPRN and 644's skip pass-throughs;
  - K3: `def_sameline --help` shows `--skip-contexts`;
  - K4: 604's identity programme source is in def_sameline;
  - K5: the builder carries `VCC_DE_MDE_SIGN_KEYGATE`, `VCC_DE_UPRN` and tercilec.
- Still not verifiable:
  - the VALUES of 612 KEYGATE / 630 UPRN / 604's identity prior (not in the recipe; reads 10-08+), so re-check K1 once
    one is decided;
  - the D/E/F bundle;
  - a real D/E/F def_day run.
- The pre-apply suites (`test_depth641`, `test_hek293_642`, `test_skip644`) now SKIP with a reason on the applied tree,
  since their pre-patch base cannot be built from live files. To re-run one, use a worktree at the pre-apply commit with
  the gitignored inputs linked (606's P1/P2 proof). `verify644.py` is the post-apply regression check.
