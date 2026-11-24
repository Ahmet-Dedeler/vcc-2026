# Final phase (Oct 22 – Nov 5): the plan IF our own D/E/F entries return their members

Written 2026-09-28 by the coordinator, before any D/E/F data exists. It is the preregistered answer to DEF-RUNBOOK §6.8's
"members come back" branch. If `vcc status <id>` returns nothing for the first D/E/F entry, ignore this file: the dark plan
(DEF-RUNBOOK §6) stands unchanged.

## 0. Why write it now
The rules hide a public final leaderboard ("no leaderboard displayed for final entries"); they do not say our own entry's
members are hidden. If they are returned, the ~27 slots after the first entry are the only measurements we will ever
have on the real test lines and the real final panel. A plan improvised on Oct 22 under time pressure wastes them. Every
switch below already has a PENDING or DECIDED value in `docs/october-decisions.json`; this file only says which of them
to re-read on D/E/F, in what order, and how to decide.

## 1. Ground rules (they override anything below)
1. **R0 first.** The first D/E/F entry is the recipe exactly as `scripts/october_recipe.py command` prints it after the
   last A/B/C read, built by `def_day.sh` (DEF-RUNBOOK). Nothing in this file delays R0.
2. **The last entry counts.** Nov 5's final upload must be the best-supported package: either the best-scoring scored
   entry re-uploaded unchanged, or a combination whose every component was scored as an improvement AND whose
   combination was itself scored at least once. Never end on an unscored combination.
3. **One variable per arm** against the current best (rule 5), with the config diff asserted by `def_postbuild.py`
   (`--allow-config-diff` exactly the intended keys) and a preregistered bin written before the build.
4. **Noise.** pds and mse reproduce to 16 significant figures between re-builds of one configuration, so bulk arms read
   exactly from a single entry. The DE members' noise depends on the arm's class (AGENTS rule 13; 537 §2.5, CONFIRMED
   2026-10-03 on 12 board pairs):
   - same emission law with pass 1 re-drawn (bulk changes that move ballast K): sd ≈ 0.00025 overall;
   - prior-only DE arms that share emission draws with their base (468): ≈ 1e-5;
   - a change of emission LAW: ≈ 0.0045 (230).

   A DE arm needs an effect > ~2 sd of its class before adoption. (Before 10-03 this item said "> ~0.009 overall, or two
   seeds", which was the emission-law figure applied to every arm.)
5. **Every build in a fresh work dir** (485 F1), serial, under the memory gate; build the next day's arms while today's
   upload scores.
6. `-m` from the model-name queue; `-d` empty.

## 2. Order of the reads (highest value of information first)
Each line: the switch, the arms, what decides it. "B" = bulk (exact), "D" = DE (noisy).

| day | arms (vs current best) | reads | adopt if |
|---|---|---|---|
| 1 (Oct 22/23) | R0 | the baseline on D/E/F | — |
| 2 | re-level ρ* ×0.98 and ×1.02 on every context (B) | mse, pds | the better side beats R0 on mse+pds overall; if both lose, R0's ρ* stands. The NTC-proxy target 0.9793 is A/B/C-calibrated; this is the single most D/E/F-specific number in the recipe |
| 3 | q 1.5 and q 2.0 (B), each re-levelled | pds, mse | overall better than the day-2 best; stop the q curve where it turns |
| 4 | LODO cap: the percentile cap vs fixed 6 (B) — whichever R0 did not ship — and CD4 ×1.5 vs ×1 (B) | pds, mse | overall better. 438 says CD4 ×1.5 loses on October-like panels; D/E/F decide it |
| 5 | same-line w: 0 vs 1 (vs R0's 0.5), only if a context matched a registry line (B+D) | all six | overall better; a w that raises the bulk and costs DE beyond 0.009 is NOT adopted |
| 6 | the programme source on vs off (D, prior-only: shares draws with its base) | fid, nmae, reach | fid + nmae + reach overall ≥ +0.0015 (458's primary bin scaled to the board) |
| 7 | the prior gate le1 vs off, if the programme is on (D, prior-only) | reach, fid | read468's rule on D/E/F |
| 8–9 | DE magnitude scale s on called genes 1.25 and 1.6 (D; 513f) | nmae | nmae ≤ −0.002 raw at pinned pds/mse |
| 10–12 | the combined best (every adopted switch at once), then its re-level ρ* re-check (B) | all | the combination must be scored before Nov 5 |
| 13–14 | reserve: re-upload the best scored package as the LAST entry on Nov 4 or 5 (rule 2) | — | — |

Two slots a day means two rows can share a day when their builds are ready. Skip a row whose switch cannot apply (no
same-line match, programme source closed by the 09-30 read). If the D/E/F bundle lands late on Oct 22, day 1 slips and
the reserve days absorb it.

## 3. What NOT to do with feedback
- Do not tune per-context factors by probing one context at a time: the members come back as panel means over three
  contexts, and three-arm designs per factor burn the window.
- Do not chase DE differences under ~0.009 overall from a single entry (ordering noise).
- Do not submit anything whose config diff against its base is not exactly the intended keys.

## 4. Pre-work that can be done before Oct 22
- Each arm above is an env change to `def_day.sh`'s call or a `--allow-config-diff` addition; list the exact commands per
  row once the 10-05 reads settle the recipe (rehearse two of them on efr, the way 497e/510 do).
- A reader script that takes a D/E/F entry pair and prints the bin verdict, like `read0928.py`.

## 5. Amendments 2026-09-29 (531's tabletop; they override §2 where they differ)
- **Cadence.** With uploads of up to ~7.5 h (09-29), plan ONE new arm a day; a second only if the first has published by
  16:00 UTC and the second is built and gated. The full list (~16 arms) does not fit Oct 23 → Nov 3 at one a day. **Drop
  order, from the bottom:** the post-combination ρ* re-check (row 12), then rows 8–9 (called-gene scale), then row 4a. Keep
  R0, the ρ* pair (row 2) and a SCORED combination (rows 10–11) before Nov 4 12:00 UTC.
- **Tooling.** `experiments/531-cloud-final-phase-tooling/readdef.py ROW ENTRY BASE [ALT]` encodes every row's bin (74
  tested verdicts); its RESULTS §2 lists each row's exact command and must-move / may-move config keys. `def_postbuild.py`
  checks "nothing outside the allowed set moved", not "exactly these moved": until it grows a `--require-config-diff`, read
  its `moved=` list by hand and reject an arm whose knob did not fire.
- **Rows that need knobs not yet built:** row 2 and 12 (`DEF_DAY_RHO_SCALE`) and row 4b (`DEF_DAY_CD4_STACK`) are in cloud 530's
  bundle; row 4a (fixed cap 6) needs a `DEF_DAY_LODO_HI` pass-through; rows 8–9 need a builder knob that scales called
  genes' delivered lfc (none exists; 513f's H1 curve puts the likely gain near +0.0006). A row whose knob is missing on its
  day is skipped, never hand-edited.
- **Row 2's ×1.02 arm** aims at ρ 0.99893, inside 1 % of the cliff: if the log says "nearest safe build" instead of
  CONVERGED, that artifact is not the arm; do not submit it.
- **Row 3 only if R0 left per-target sharpening off**: with `DEF_DAY_SHARPEN_PT` set, the positional q is ignored.

## 6. Amendment 2026-10-02: feedback-use policy (H-631) and the read design
- **Policy: maximum-information use, fully documented. DECIDED** (Ahmet, 2026-10-02; his words are in DEF-RUNBOOK §6
  item 8). §3's "do not" list stays as slot-efficiency advice. It does not limit what may be read.
- **Ledger.** Every D/E/F entry goes into `docs/BOARD-LEDGER.md` with what it was designed to read, the read, and how
  the read was used. The finalist write-up discloses it.
- **Read design.**
  - 637 §1's factorial reader (one bulk switch + one DE switch per entry) is the default once its A/B/C validation
    passes (the all-switches arm).
  - 637 §4's fused module decode (bulk switch X_k + Hadamard module pattern h_k per entry) is APPROVED in principle
    (Ahmet ~12:40Z: "we can do the risky thing"). It is adopted once it meets two acceptance criteria. [7772dc] is
    pricing them and will propose the module count: in a Broad-16 sim, 8 CLEAN4 octiles keep most of the value (+0.0165
    vs +0.0172 for 16) at ~10 entries.
    1. **Cadence.** §5 plans one new arm a day (uploads up to ~7.5 h). A full 16-row decode plus R0 and a scored
       combination needs ~18–20 entries; Oct 23 → Nov 4 at one a day gives ~13. Price either a partial decode (least
       squares, rank reported) or fewer, larger modules.
    2. **DE re-roll inside each equation.** A bulk switch that moves any block's ballast count K re-draws pass 1
       (537 §2.6). Each entry's DE delta then carries the re-roll, about σ ≈ 0.00025 overall (547, 637 §4's leak
       −0.0003 ± 0.0002), on top of the module pattern's effect. With the 16-row Sylvester design, a decoded module net
       has noise sd ≈ (2/16)·√15·σ ≈ 0.48σ (≈ 1.1σ for module 0, whose column sums to 16), so ≈ 0.00012 overall.
       Nets near 0.0002 overall are then unreadable. The fix is K-neutral bulk switches, or per-block emission streams
       (an emitter change, after 10-06).
  - 625's A/B/C run measures the real module net sizes and rehearses the decode end to end. It prices both costs.
- **Update 10-02 13:57Z (625 Amendment 3, fe2064b7): the D/E/F module decode is DROPPED.**
  - 625 FINAL (73c0b648): sign improvements reach the score at ≈ +0.002 overall per +0.01 gate sign (the earlier +0.0011 was
    degradation-only). The readout is worth ≈ +0.003 on D/E/F, still too little for an 8–16-entry decode in the final phase.
  - New plan: A/B/C reads 8 single-octile prior-only probes plus 1 combination entry (exact reads; no Hadamard needed).
  - If GO, the mask enters R0 (DEF_DAY_DE_PRIOR_FLIP, 604f) and ONE prior-only D/E/F entry verifies it. If NULL, no D/E/F entry
    is spent.
  - The two acceptance criteria above therefore apply to no current design. The factorial reader stays the default for switch
    reads.

