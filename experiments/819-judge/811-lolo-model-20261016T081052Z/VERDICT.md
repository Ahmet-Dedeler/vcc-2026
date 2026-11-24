# VERDICT: experiments/811-lolo-model/PROPOSAL.md (M811) — **REVISE**

Judge run 819, written 2026-10-16 08:15Z (clock). I judged the frozen copy
`experiments/819-judge/811-lolo-model-20261016T081052Z/PROPOSAL-frozen.md`, sha256
`cc5b8a8fdc2ee07bf790d5e25f4f6b8c3886e0f33155efe1bb303e8e70e2b7f9`. The live `PROPOSAL.md` had the same sha at 08:13Z. This is
the first submission: there is no earlier `VERDICT-*.md` in `experiments/811-lolo-model/`. I computed nothing beyond lookups,
including reading fields from `work/run-base811/agg811-*.json`.

**Summary.**
- **Why not KILL.** No closure covers a per-target, many-line consensus with per-line centring, blended into the composition.
  The baseline falsified its own preregistered NULL for "more lines": all vs src is +0.0018 / +0.0050 / +0.0105 on BULK-W,
  17/18 folds (`work/table-base811.md` rows "all_lX vs src_lX | W"). The run costs about $0.5.
- **Why not GO, in five points:**
  - The line set it would select on (L folds), confirm on (BULK-ABC) and ship on (D/E/F) is not the same set. The KOLF2
    bank that motivates factor C exists only on BULK-ABC.
  - The mechanism paragraph omits the proposal's own per-target diagnostic, which contradicts it on 4 of the 6
    surface × line cells.
  - Factor W is agreement weighting. That has been closed three times, and the proposal misstates the closure.
  - P1 and the first falsifier cannot be computed from the planned run.
  - The try count and the decision rule leave out the family's 6 configurations, already read on every confirmation fold.

---

## The JUDGE.md checks

### 1. Is the mechanism closed already? **Not at this construction. Factor W is close to a closed construction, and the proposal cites that closure wrongly.**
- **The consensus blend itself is not closed.**
  - "Direction via a fifth donor line: sources saturate" (`docs/FINDINGS.md:885`) covers donors inside the whitened stack.
    It does not cover a raw many-line consensus blended at λ.
  - 343–382 "consensus priors … coherence is the price" (`docs/FINDINGS.md:870`) covers SHARED-across-targets priors. U(t)
    is per target.
  - 706/707 (`experiments/706-strong-tier-budget/RESULTS.md:34-71`) are sources only, top tier. The proposal's §5 names
    the difference correctly (`PROPOSAL-frozen.md:62-63`).
- **Factor C (per-line centring) has close precedents, none at this construction:**
  - **706 §2's control: "centring U changes nothing"** (`experiments/706-strong-tier-budget/RESULTS.md:40-41`), on the
    sources-only U, whose shared energy was 0.2–1.8 %. The proposal does not cite it. It is the closest prior result to C.
  - 075/079: reciprocal kNN *local* centering, post-whitening. Closed on KOLF (`experiments/079-kolf-local-centering/RESULTS.md`,
    "Completed successfully and rejected").
  - 393: rescaling the bulk's shared mean-across-targets component, ≤ +0.0015 (`docs/FINDINGS.md:879`).
  - 483: projecting other lines' top-k shared directions out of the bulk lowers pds (`docs/FINDINGS.md:854`).
  - 541–544: a "blob-free" key (`docs/FINDINGS.md:836`).
  - Per-line centring of the *donor rows before a many-line average*, where U_all carries 2–9× the truth's shared energy,
    is a different construction. But the proposal must say why it expects an effect where 706 found none (item B).
- **Factor W (reliability weights) is agreement weighting, which is closed:**
  - MODZ "weights each source by its mean correlation with the *others*". It lost to equal weights on a 4-line LODO
    (`experiments/016-modz/RESULTS.md:8, 26`).
  - 142 nested selection chose η = 0, exact equal weights, "terminal NO-GO"
    (`experiments/142-batch-reliability-weighting/RESULTS.md:3, 15`).
  - 271: "the reliability-weighted mean loses to the equal-weight mean (−0.0037 pds)"
    (`experiments/271-donor-reliability-estimation/RESULTS.md:1`).
- **The proposal's distinction is factually wrong.** `PROPOSAL-frozen.md:58-59` says these closures "re-weighted the three
  production donors inside the whitened stack". 016 was a raw 4-line LODO with H1 as the fourth line; whitening arrived
  the same day, in 017/019. 142 and 271 are not cited at all.
- **What W does differ in:** line-level weights in place of per-signature weights, and 6–15 lines in place of 3–4. That is a
  variant, not a new mechanism, and the proposal's own P2 prices it at ≤ +0.0005 (`:46`).

### 2. Falsifier and bins written before any number? **Partly. P1 and falsifier 1 cannot be computed as planned. No decision rule.**
- Present: P1–P3 and two falsifiers (`PROPOSAL-frozen.md:43-52`).
- **P1 is stated on BULK-ABC at λ 0.25** (`:44`). The protocol prices only the ONE selected configuration on BULK-ABC
  (`:37-38`). If the W-tune winner is λ 0.10 or C off, P1 is never read.
- **Falsifier 1** ("centring removes ≥ 70 % of the r gain", `:49`) names no surface. It does not say whether the r gain is
  measured against `inc` or against the uncentred arm at the same λ.
- **Bins are only referenced** ("as DESIGN §6/§8", `:48`):
  - there is no BULK-W numeric prediction for the winner, member by member;
  - there is no WEAK bin;
  - the proposal does not say what happens if W-conf misses its "±30 % of W-tune" condition;
  - the proposal does not say which result produces a DOSSIER and which closes the family.
- pdsl is not mentioned. PHASEB requires it beside literal pds (`experiments/800-checker/PHASEB.md:9-11, 53`).

### 3. Measurable on the checker's surfaces? **mse yes. The design's key confound cannot be seen on the selection surface. DE-H1 has a leak.**
- **Line-set mismatch (decisive).**
  - The local KOLF2 bank carries 273 targets on A/B/C's panel and 0 on every L panel
    (`experiments/811-lolo-model/RESULTS.md:45`).
  - The KOLF2 ablation (−0.009 / −0.009 / −0.015, `RESULTS.md:159`) and the 4.3–4.6 % shared energy on A/B/C
    (`RESULTS.md:148-150`, against 2.1–2.5 % on L) are both A/B/C-only.
  - So:
    - **C is selected on W-tune (L1–L3), where its motivating culprit is absent;**
    - it is **confirmed on BULK-ABC with a bank the D/E/F day will not have** (full KOLF2 rows are not on disk,
      `PROPOSAL-frozen.md:74-76`);
    - the BULK-ABC read therefore prices a different U from the one the day ships.
- **The DE-H1 route** (`:39-40`) does not exclude H1:
  - The consensus's H1 group includes `h1_vcc2025` = `data/donors/h1_2025_full` (`PREREGISTRATION.md:16`).
  - DE-H1 requires that "the h1_2025 donor never enters" (`experiments/800-checker/DESIGN.md:48`; AGENTS §6.5 oracle
    leakage).
  - The proposal must hold the whole H1 group out on DE-H1, and say which banks then carry h1-test / h1-trainval targets.
- **The unmeasured DE members are pinned acceptably.** Ballast K is reported, and |dK_rel| > 5 % goes to DE-H1
  (`:39-40`; DESIGN §4 `experiments/800-checker/DESIGN.md:70`). The baseline moves K by up to +5.9 % at λ 0.25
  (`RESULTS.md:136`).
- **PARTIAL.** Off-axis change energy reaches 0.30 on some folds (`RESULTS.md:138`), above DESIGN §7's 25 % line
  (`DESIGN.md:91-92`). The proposal's bins do not mention PARTIAL.
- **The BULK-ABC operator.** A/B/C's fold operators are 062's unpruned whiteners, and 801 is refitting them as ABCp
  (`docs/LAB-LOG.md`, 08:05Z entry). The proposal does not say which A/B/C folds it confirms on.

### 4. Does it use seen data as evidence of gain? **No: no A/B/C board reads, no sealed reserve, no held line in U. Not a KILL. The confirmation surfaces are worn by this same family.**
- The held line's whole group is out on every fold (`PREREGISTRATION.md:21`).
- LZ is left to the lab lead (`PROPOSAL-frozen.md:41`).
- Wear on the confirmation surfaces:
  - the baseline's 6 configurations were priced on all 22 folds, W-conf (L4–L6) and BULK-ABC included (`RESULTS.md:166`);
  - C and the KOLF2 culprit were chosen by reading BULK-ABC's fold truths in align811 (`RESULTS.md:141-163`);
  - so neither confirmation set is fresh for this family. Item 5 charges it.

### 5. Is the try budget honest? **No, it is under-counted on the confirmation reads.**
- The proposal counts 8 on W-tune, 1 on W-conf and 1 on BULK-ABC (`PROPOSAL-frozen.md:71`), and says the baseline's 6
  "stay on record".
- DESIGN §8 counts N per family per surface set (`DESIGN.md:99-102`). On L1–L3 that gives 6 + 8 = 14, above 10. Nested
  selection is the right answer to that, and the proposal has it (`:35-38`).
- On W-conf and BULK-ABC the family has already priced 6 configurations, and the ABC truths shaped the grid. A confirmation
  read at N = 1 (no shift) is not honest. Use e(7) ≈ 1.35 sd, or state and defend another count.
- P3's "2.5 % bound > 0" (`:47`) does not say which shift it is after.

### 6. Is the expected gain worth its cost? **The cost, yes. Factor W, no.**
- About $0.5 on EC2, plus about $0.4 for DE-H1 (`PROPOSAL-frozen.md:69-70`). That is far under the $50 cap and the lane's
  $50 (the lane has spent $0.45, `RESULTS.md:196`).
- The P3 bar of +0.0010 on BULK-ABC is above the ~0.0003 readability of bulk direction (JUDGE.md).
- **W's predicted effect is ≤ +0.0005** (`:46`), at the readability floor, on a construction closed three times. It doubles
  the grid, and so the best-of-N shift, for an effect the proposal expects to be unreadable.
- **Calibration.** The class is bulk-external, with one board point (702, negative).
  - ρ is 0.93 on BULK-ABC and 0.19 on BULK-W (`docs/LAB-LOG.md`, Known-answer report part 1).
  - The aggregator's W projections divide by ρ 0.187. For example, it projects all_l25 vs src_l25 to a board +0.023
    (`work/run-base811/agg811-811-consensus-e8e1a491.json`, `projection`). Such numbers are not usable.
  - The proposal does not quote them, which is right, but it must say that no board projection is claimed beyond a range
    that includes 0 (DESIGN §7, `DESIGN.md:89-90`).

### Day-of path
Mostly acceptable. Centring and weights come from D/E/F's own panel, the rows come from public banks, and nothing is
fitted to A/B/C (`PROPOSAL-frozen.md:72-73`). Class: transfer. Two gaps remain:
1. The day's bank list must equal the priced list (item A).
2. "D/E/F's line never in any bank" (`:21`) is asserted, not guaranteed. A D/E/F context could be one of the bank lines
   (K562, HCT116, HEK293T, Jurkat, A549, MCF7, HT29, DLD-1 …). The proposal needs a rule: use the identity call, and hold
   out or keep the matched line's group, with the reason. The oracle_vip probe shows that a shallow same-line bank loses
   badly (−34 … −64 % r, `RESULTS.md:121`), so this case matters.

---

## What is missing (fix all of these, then resubmit once, as a frozen file)

**A. One bank list on every surface.**
- Freeze the line/bank list that the D/E/F day will have, and use exactly that list on W-tune, W-conf, BULK-ABC and the
  day path. Today that list has no local KOLF2 bank.
- Price "with KOLF2" only later, as a separate counted try, once 813 wires the full rows.
- Add the D/E/F identity rule from the Day-of section.

**B. Correct the mechanism paragraph** (`PROPOSAL-frozen.md:7-19`).
1. State that the per-target raw cosine does NOT rise with more lines for held HCT116 / HEK293T. On BULK-ABC it falls
   (0.0381 → 0.0331, 0.0376 → 0.0238). On L it is flat (`RESULTS.md:148-154`). Only the K562 folds rise. The engine gain is
   energy-weighted (`RESULTS.md:161-163`). So "more lines carry more transferable direction" is the hypothesis the falsifier
   tests, not a finding.
2. Correct "2–4×": the ratio is 3.9–8.6× on BULK-ABC and 1.4–4.5× on L (`RESULTS.md:148-154`).
3. Cite 706's "centring U changes nothing" (`706 RESULTS.md:40-41`).
4. Say quantitatively how removing a component that holds 2–5 % of U's energy can cut the BULK-ABC pds cost by ≥ 60 % (P1),
   or revise P1.

**C. Factor W.**
- Recommended: drop it (grid 4, C × λ).
- If it stays, cite 016 (`:8, 26`), 142 (`:3, 15`) and 271 (`:1`) correctly, give the mechanism by which line-level
  agreement weights differ from them, and add a W-only falsifier.

**D. Make every prediction computable from the planned run.**
- State P1 and falsifier 1 on W-tune, where all configurations are priced.
- Alternatively, pre-commit the C-on / C-off pair at the winner's λ on BULK-ABC as mechanism reads, counted as tries.
- Define the "r gain" baseline: against `inc`, or against all_lX at the same λ.
- Give each P with literal pds and pdsl side by side.

**E. Honest tries.**
- State N per surface set: L1–L3 = 14 (6 + 8); W-conf and BULK-ABC ≥ 7.
- Apply e(N) to the confirmation reads, or defend N = 1 explicitly against the align811 read of the ABC truths.

**F. A decision rule written before the run.**
- Suggested DOSSIER condition, all of:
  - BULK-ABC lower bound > 0 after e(N), literal pds, with pdsl reported;
  - paired > all_lX at the same λ;
  - W-conf lower bound > 0 after e(N);
  - the PARTIAL label if off-axis > 25 %;
  - DE-H1 clean if |dK_rel| > 5 %.
- Say which outcome closes the family, and what a W-conf miss of the ±30 % condition means.
- Board projection: a range including 0 (single negative board point for bulk-external).

**G. DE-H1 hygiene.**
- The whole H1 group (`h1_pan2025`, `h1_vcc2025` = `data/donors/h1_2025_full`) goes out on DE-H1.
- List the banks that form U there.

**H. Provenance and instrument.**
- At freeze time (08:10:52Z) the proposal cited "RESULTS.md §1/§2" (`:3, :22`). That file was written at 08:11:53Z
  (now sha256 `f15f0298204a5200…`). Cite it by sha256 in the resubmission.
- Name the M811 arm module and record its sha256 before the run.
- Say whether BULK-ABC means 062's folds or 801's ABCp refit.
- Freeze the file before the judge reads it.
