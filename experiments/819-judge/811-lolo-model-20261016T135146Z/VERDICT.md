# VERDICT: 811 PROPOSAL-3 (sign skill on responding genes) — **KILL** (the D3 diagnostic is declined; the proposal's own "no construction" answer is accepted as the lane's closure)

Judge run 819. Written 2026-10-16 13:54Z by the clock. I judged the FROZEN copy
`experiments/819-judge/811-lolo-model-20261016T135146Z/PROPOSAL-frozen.md`, sha256
`200755d77a9a879a372237e1f5ca3e4911c787dee0af29ff6f1b024204994543`. I re-hashed it, and it matches the mandate. It also equals
the live `experiments/811-lolo-model/PROPOSAL-3.md` and the line in `811/work/SHA256SUMS.txt:8`.

Other files read:
- the request: `experiments/800-checker/resume/811-proposal3.md`;
- 811's `MANDATE.md`, `VERDICT-1.md` and `VERDICT-2.md`;
- `800-checker/PHASEB.md`;
- LAB-LOG's two known-answer reports (`docs/LAB-LOG.md:136`, `:171`);
- FINDINGS §1–§2;
- the RESULTS files the proposal cites: 516, 621 (via FINDINGS), 430, 803 §8, 811, 812, 813 and 702;
- `812/VERDICT-u1hook.md`.

Computed: sha256 lookups (`m811.py` = `e00192e6…`, which matches §5; `wear811.jsonl` = `82663b4b…`, which matches §6) and greps.
Nothing else.

**What is being judged.** PROPOSAL-3 offers no construction (§0). That is a valid, well-supported negative, and it needs no
compute. The one fundable item is **D3** (§4): skill812 reads the sign content of c0_l25's day-list consensus U on H1's
targets. KILL applies to D3, and the objective it would serve is recorded as closed for lane 811. **c0_l25 is untouched.** It
stays the lab's frozen 6a candidate, to be read once on D/E/F-panel folds (LAB-LOG 10:49Z, "Lead plan for 811").

**Why KILL rather than GO at $0.10.**
1. D3 is the third H1 read of a construction the record has already priced and closed:
   - 516 read a same-target public consensus on the essential-free H1 targets;
   - 812 read the same object three ways;
   - the 812 read was KILLed as H-U1PUB.
2. Neither of D3's outcomes can lead to compute:
   - below the bar, it closes what is already closed;
   - above the bar, it returns as H-U1PUB with a different bank list. 812's verdict names that explicitly as not a reopening.
3. Its "no new arithmetic" claim is false. skill812 has no essential-free subset, and the essential list D3 would filter by
   is not pinned.

---

## JUDGE.md checks

### 1. Is the mechanism closed already? **Yes, at its construction → KILL.**
- **D3's U is 812's killed source with small changes.**
  - U = `m811.day_consensus(..., centre=False, hold={"H1"})` (`m811.py:86-90`). That is `base811.consensus` over
    `DAY_GROUPS` = `base811.GROUPS` minus KOLF2 (`m811.py:41`; `base811.py:42-49`). The lines are K562 (gwps, VIPerturb,
    Jiang), HCT116, HEK293T, CD4, TeloHAEC, DLD1, A549, RPE1, Jurkat, BXPC3, HAP1, HT29 and MCF7. Each line gets equal weight.
  - 812's `pubmean_raw` / `pubmean_rms` is the per-target mean over the same 12 non-H1 lines: banks are averaged within a
    line, H1 is never read, and CD4 is absent (`812/RESULTS.md:208-215`).
  - The only differences:
    - CD4 is added. It covers 5 / 11 H1 targets (`803/RESULTS.md:416`).
    - The unit scaling is base811's instead of bank-RMS.
  - On the H1 surfaces, D3's U is therefore 812's `pubmean_rms` plus CD4 on ≤ 16 of 211 targets.
- **That source is closed, with its numbers.**
  - 812 table (`812/RESULTS.md:224-227`): a0 vs CLEAN4 is −0.0023 to −0.0025, with every CI below 0. u1 vs the incumbent is
    +0.0022. a1 is −0.0011 to −0.0014.
  - KILL in `812/VERDICT-u1hook.md:1-30` and `:126-134`. It is 516's closed rider. Its gain is unwhitened signs replacing
    carrier signs at CPM 5–20 (430 / 417c). Its H1 trigger carries the essential-donor banks. About half of it is a line
    programme.
  - FINDINGS `docs/FINDINGS.md:266-269` (621): "public sign sources are closed for this round".
- **D3's one new element, the essential-free (`no_ess`) H1 subset, was read for this object type in 516.**
  - 516: "Panel-like H1 targets (the 64 that no essential screen covers) +0.0118 [−0.0007, +0.0226]"
    (`516-de-sign-sources/RESULTS.md:13-14, 81`).
  - That read did not carry to the held-out K562 fold: +0.0010 [−0.0030, +0.0053] (`:20`). The plain mean is −0.0023 there
    (812 verdict `:24-26`, citing 516 l.204).
- **The reopening path 812's judge left does not run through D3.** `812/VERDICT-u1hook.md:118-124` asks for three things:
  1. a gain that survives on the essential-free subset AND after per-line centring;
  2. a held-out DE truth pinned by path and sha256 before any read;
  3. delivery "as one of lane 811's model outputs".

  D3 meets only (3):
  - it is uncentred by design (`PROPOSAL-frozen.md` §4 "Construction", `centre=False`), and c0 is frozen uncentred;
  - it names no held-out-line DE truth.

  Its above-bar branch is therefore 812's killed hook source with c0's bank list (`PROPOSAL-frozen.md` §4 "At or above the
  bar"). 812's verdict calls that "a new bank list", not "a new mechanism" (`:118`).
- **The proposal's §2 mechanism holds against the record, with one caveat.**
  - 803 §8 shows that c0 as delivered leaves fid / jac ≤ 4e-5 raw and reach inside seed SD (`803/RESULTS.md:358, 442`).
  - But the transform acts on the composition rows only. The DE carriers and priors are unchanged (`803/RESULTS.md:370-375`).
  - So 803 shows that c0 adds no sign change *on the delivered path*. It does not, by itself, show that U has no sign
    content. That second claim rests on 812's direct read of the same source (above), and that read closes it.
- **One citation needs correcting for the log.** §1 quotes 516's HCT116 +0.0033 and HEK293T −0.0002. Those are Orion fold
  rows, VOID per `616-clean-orion-folds/AUDIT-orion-dependence.md:57`, as cited in `812/VERDICT-u1hook.md:37-39`. What
  remains of 516 is K562 plus H1. The conclusion does not change.

### 2. Does the prediction have a falsifier and bins written before any number? **Yes, formally. But the predicted outcome fails the proposal's own bar.**
- Bins and bar are at `PROPOSAL-frozen.md` §4 ("Prediction … before any number"; "Bar (812's own) and decision"). They are
  812's bar (`812/RESULTS.md:244-247`).
- The proposal predicts a0 below CLEAN4 by −0.002 … −0.003. That matches 812's CIs, which exclude 0. The bar requires "a0 not
  below CLEAN4 by more than its CI".
  - So the proposal's central prediction already fails the bar on a0, whatever u1 does.
  - Its expected outcome is the "below" branch (§4: "I expect this outcome"). The run would buy confirmation of a closure
    already on record.

### 3. Is it measurable on the checker's surfaces? **Partly. The H1 read is measurable; its transfer is not, and the subset is not built.**
- skill812 takes H1-free predictions (`812/skill812.py:232-235`: `--name`, `--pred`, `--h1-free`, `--ka`). It has **no
  target-subset or essential-free option**.
  - The `no_ess` report would need new code, so the claim "no new arithmetic" (§4 "Construction (frozen, no new
    arithmetic)") is false.
  - The essential list is also unpinned. Two lists are candidates: 516's three screens (K562-ess, RPE1, Nadig Jurkat;
    `516/RESULTS.md:81`), or the lab's esslib (`814/work/flags/esslib.txt`, 2,394, `24fffbfe…`; DESIGN §3.1 amendment 2,
    `DESIGN.md:73-77`).
- Transfer of sign priors from H1 to the board is 0.07–0.67× on fid (`FINDINGS.md:722`, cited in §3 of the proposal itself).
  - The checker's DE-sign class has one hit (CLEAN4 hook, ρ 0.95) and one unexplained miss (idall3), both on
    `docs/LAB-LOG.md:184-196`.
  - A held-out cancer-line DE truth is what would make a sign read decision-grade. None is on disk:
    - 616 voided the Orion DE truths;
    - 813 left TeloHAEC unvalidated (`813/RESULTS.md:296-298`).

### 4. Does it use seen data as evidence of gain? **No.**
- The H1 group is held out of U (`m811.py:86-88`; `803/RESULTS.md:396-397`).
- No A/B/C board read and no sealed reserve is used.
- DE-H1 is worn, but it is not a training line of U. This check passes.

### 5. Is the try budget honest? **No. The family is under-counted.**
- §4 opens a fresh family, `811-signs`, at N = 1, and assigns 812's four source reads to "812's family".
- The object is the same: a per-target mean of the same 12 non-H1 lines, read by the same instrument on the same two
  surfaces. 516's `no_ess` read is the same object type on the same subset.
- Best-of-N is charged per construction, not per lane (AGENTS 11.1c.4). Honest N for this sign source on DE-H1 is therefore
  ≥ 5:
  - 812's 3 public-mean forms plus the bulk sign (`812/RESULTS.md:206-207`, "tries 4");
  - D3;
  - plus 516's family, if the Stouffer forms count.
- The selection over the class columns a0 / a1 / a1head / u1 also inflates the effective N (`812/VERDICT-u1hook.md`
  check 5).

### 6. Is the expected gain worth its cost? **No, even at $0.10.**
- The proposal's own predicted u1 + a1 is < +0.002 overall-eq at H1 rates (§4). Apply:
  - the sign-prior transfer factor 0.07–0.67× (`FINDINGS.md:722`);
  - 516's fold-to-H1 ratio. K562 +0.0010 vs H1 +0.0157 is about 1/16 (`516/RESULTS.md:11, 20`).
- The expected board-relevant value is then < +0.0003. That is below the readable floor of a prior-only DE arm and inside
  best-of-N's bias at the re-roll sd 0.00025 (AGENTS §11.1 rule 13, 11.1c.4).
- The cash is trivial. The decision value is zero: no outcome changes c0_l25 (§5) or reopens a hook source (check 1). The
  run would add one more read to DE-H1, which is worn, for a family the record has closed four times (516, 621, 525 / 634,
  812).

### Day-of path
Not applicable: D3 is diagnostic only and ships nothing. c0_l25's own day path was judged in VERDICT-2 ("Day-of path:
acceptable") and is unchanged (`m811.py` sha `e00192e6…`, which matches §5).

### The lead's four asks (`800-checker/resume/811-proposal3.md`)
- **Ask 1, mechanism and closures:** answered. §1's table matches the record line for line: 516 `:10-22`; 621
  `FINDINGS.md:266-269`; 430 `FINDINGS.md:889`; 458d `FINDINGS.md:322`; 812 `RESULTS.md:224-240`.
  - The only error is the Orion fold rows quoted from 516 (check 1).
- **Asks 2–3 (instrument, bins, class, tries, cost):** they apply only to a construction. None is offered, and that is
  justified.
- **Ask 4, relation to c0_l25:** answered. There is no swap and no rider.
  - The standing evidence is quoted correctly: ABCp +0.00264 [−0.00003, +0.00561] 3/3 (`811/RESULTS.md:187`); LABC
    +0.00230 3/3 (`:115`); A-only increment +0.00193 14/18 with Δpds +0.00447 (`:52`); DE-H1 no cost (`803/RESULTS.md:356-363`).

---

## Recorded closure (for docs/LAB-LOG.md, in the lab lead's hand)
**811 sign skill on responding genes (PROPOSAL-3 / D3): KILLED (judge, 2026-10-16 13:54Z).**
- **What was proposed:** no construction. Diagnostic D3 was skill812 on c0_l25's day-list consensus U (H1 held out,
  uncentred), with an essential-free report.
- **Why it is closed:**
  - U on H1's targets is 812's `pubmean` (same 12 non-H1 lines, equal line weight) plus CD4 on ≤ 16 targets. That source was
    KILLed as H-U1PUB (`812/VERDICT-u1hook.md`), and its a0 loss (−0.0023 … −0.0025, CI < 0) already fails 812's bar.
  - The essential-free H1 read of a public same-target consensus exists: 516 measured +0.0118 [−0.0007, +0.0226] on 64
    targets, and it did not carry to K562.
  - The above-bar branch would be a new bank list for a killed hook source. 812's judge requires centring and a pinned
    held-out DE truth, and D3 has neither.
- **Mechanism, for the record:**
  - On the Arc-drawn panels, the many-line increment works through target discrimination: pds on A-only, K562-held folds
    (R_tgt +0.0015 … +0.0034).
  - c0 as delivered does not touch the DE carriers (803 §8.1), and the DE members do not move (≤ 4e-5).
  - A per-gene sign source would need new information: a public line near D/E/F's lineages at depth, or a held-out
    non-H1 DE truth. None is on disk.
- **Tries charged to the sign-source family on DE-H1:** 812's 4 reads + 516's forms. D3 was not run.
- **What would reopen it:**
  - a DE truth on a held-out non-H1 line, built and pinned by path + sha256 before any read (for example a K562-GWPS cell
    fold, 805's DE-CD4 with CD4 out of U, or a validated KOLF2.1J cell surface per `813/RESULTS.md:63, 111`), on which a
    centred source's u1 + a1 clears 812's bar;
  - or a public line near D/E/F's lineages identified after the 10-22 release.
- **c0_l25:** unchanged. It is the lab's frozen 6a candidate, read once on D/E/F-panel folds after the window (LAB-LOG 13:46Z
  day plan).
- **Budget:** none granted.
