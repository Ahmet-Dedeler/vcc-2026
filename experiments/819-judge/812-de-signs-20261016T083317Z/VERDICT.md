# VERDICT: 812 PROPOSAL 2 (H-U1PUB, the uncalled-gene hook signed by the public per-target mean): **KILL**

Judge run 2026-10-16 08:36Z. Judged file: the frozen copy
`experiments/819-judge/812-de-signs-20261016T083317Z/PROPOSAL-frozen.md`, sha256
`de2d2295a1c8560e18816e224a37865744a91ccea46c3408d72794959d7f412c`. It is byte-identical to the live
`experiments/812-de-signs/PROPOSAL-u1hook.md` at judging time. No earlier `VERDICT-*.md` exists in `experiments/812-de-signs/`,
so this is a first submission. Nothing was computed; every check below is a lookup.

**The verdict.** The construction is 516's closed one: same-target public rows, averaged, signing uncalled genes through
496's hook. The proposal says so itself (PROPOSAL-frozen l.23–26: "516 closed this family … My reading is the same object").
Three things are new, and none is a new mechanism:
- 12 lines instead of 6 rows, with CD4 dropped;
- a plain mean instead of Stouffer;
- the read is taken against the incumbent's delivered hook rather than against production's pre-hook sign.

The H1 read that triggers the proposal also carries the essential-donor artefact that closed 430 (below). That artefact cannot
exist on D/E/F's panel. The closure is recorded with its mechanism. The question "does a many-line public signal carry DE signs
to a held-out line" stays reachable through lane 811: skill812 scores 811's H1-free per-target predictions, and a held-out
confirmation is already part of that lane's bar (812 RESULTS l.244–247). It does not need a separate 812 compute line.

---

## JUDGE.md checks

### 1. Is the mechanism closed already? **Yes: closed at its construction → KILL.**
- **516 is the same object, already priced as this rider.** `experiments/516-de-sign-sources/RESULTS.md`:
  - l.1 and l.10–21: the best source is a consensus of same-target rows on uncalled genes. It gains +0.016 on H1 but is ≈ 0
    on the cancer folds.
  - l.27–36: the H1 gain is (a) the unwhitened-carrier effect, about +0.017, which inverted on the board as 417c (nmae +0.0089
    raw), plus (b) consensus, which is 0.000 ± 0.005 on the folds.
  - l.49–52: priced as a rider on 496's hook: "Fold-predicted A/B/C value … ≤ +0.001 overall … **Not a lever. It is at most a
    rider for a DE-sign build made for another reason** (496's uncalled-gene hook)". The proposal is exactly that rider,
    with no other build to ride on.
- **The plain mean, which is 812's form, already lost to the production analog on the one valid held-out line.** 516 §6, l.204
  ("donors only" = the plain donor mean): HCT116 −0.0007, HEK293T −0.0004, **K562 −0.0023** against the analog. Stouffer
  standardisation was the only part that ever beat the plain mean on the folds (l.202). 812 drops it.
- **The Orion void does not reopen 516.** `experiments/616-clean-orion-folds/AUDIT-orion-dependence.md`:
  - l.57: 516's Orion fold rows are VOID, and what remains is K562 plus H1.
  - l.93: "516 / 621 / … NEGATIVES | closed | … | H1 / K562 / Broad rows | all still NEGATIVE without Orion".
  - `experiments/617-closed-family-audit/AUDIT.md` l.57 agrees.
- **FINDINGS closes the family three times.** In `docs/FINDINGS.md`:
  - l.381–385 (516): "no target-specific public sign source transfers … 'Sign sources are at the public-data limit' now covers
    target-specific combinations too".
  - l.249–252 (621): "With 516 / 525 / 634, **public sign sources are closed for this round**".
  - §1 l.626–628: "Unwhitened donor means, blends and TF regulons do not beat the whitened composition on cross-study cancer
    folds (430, 108, 190)".
- **The part of the gain that beats the incumbent is the 430 / 417c construction.** The proposal places most of its
  advantage over the incumbent at control CPM 5–20 (PROPOSAL-frozen l.14–15; 812 RESULTS l.237). There the hook has no
  confident CLEAN4 sign, so the swap replaces the emitter's whitened carrier sign with an unwhitened donor-mean sign. That is
  §2's "unwhitened donor-mean or blended carrier signs (430: cross-study +0.001, essential-donor artefact on H1)" (FINDINGS
  l.889). It is also 516's component (a).
- **The H1 trigger includes the essential-donor artefact.**
  - 812's bank includes K562-essential (43 H1 targets), RPE1 (57) and Jurkat (56) (812 RESULTS l.210).
  - `experiments/430-cancer-fold-signs/RESULTS.md` l.26–28: H1's unwhitened-mean sign gain "comes from K562-essential, RPE1
    and Nadig-Jurkat … those targets hold **60 % of H1's real-gate cells**. The same three donors cover **0 of the 300 A/B/C
    targets**". With only the panel-covering donors, the gain reverses to −0.0017.
  - The panel has zero essential genes by construction (FINDINGS l.1622–1625, 8.1g; AGENTS §2.1), so those banks contribute
    nothing on D/E/F.
  - The proposal reports no essential-free (`no_ess`) subset. Its H1 number therefore over-states what its own day-of path
    (PROPOSAL-frozen l.46) can deliver.
- **About half of the gain is a line programme.**
  - Per-line centring removes about half of the u1 gain (812 RESULTS l.239–240; PROPOSAL-frozen l.12: centred +0.0007, CI
    spans 0).
  - The remainder is the shared across-target component of the banks, which is a programme prior.
  - That family is closed: FINDINGS l.671–677 (programmes are lab signatures; a multi-lab consensus adds ~0) and l.305–307
    (458d HARM; "October runs no programme").

### 2. Does the prediction have a falsifier and bins written before any number? **Partly.** This would only matter for a REVISE.
- Stage 0 has a falsifier and a decision rule written before its run (PROPOSAL-frozen l.33–35). The rule is ≥ 2 held-out lines
  at ≥ +0.01 u1 accuracy with a lower bound > 0 against the incumbent hook.
- The trigger itself is post hoc. The u1 class was picked after six class columns × 4 sources had been read (812 RESULTS
  l.222–228), so the trigger is a selected cell, not a preregistered prediction.
- Stage 1's prediction (PROPOSAL-frozen l.44: Δnmae −0.002 … −0.004 raw, other members ≤ 1e-5) has no SHIP / WEAK / NULL /
  HARM bins at the checker's rule and no falsifier.

### 3. Is it measurable on the checker's surfaces? **No, as written.**
- **Stage 0's held-out surface does not exist yet.** PROPOSAL-frozen l.36 says the DE truth "is to be located".
  - The Orion DE truths are void (616 S1: 94–98 % of their calls were artefact; ~13–18 scoreable targets per clean fold).
  - The one valid non-H1 DE truth on record for this question is 516's K562 fold: 3,179 pool cells, with a CI of about
    ±0.004 at 516 l.21.
  - "≥ 2 held-out lines" therefore cannot be met with what is on disk, and no path or sha is named.
  - "vs the incumbent hook" on a non-H1 line has no realized emitter call set or delivered sign. Only DE-H1 has one. 516 had
    to use an analog for this reason.
- **Stage 1's surface (DE-H1) is the one where this family flatters.** The proposal concedes it (PROPOSAL-frozen l.28–29).
  - FINDINGS l.721–722: "DE sign priors 0.07–0.67× on fid; carrier nmae and sign-prior reach not priceable on H1".
  - LAB-LOG l.164: the checker's DE-sign class has one hit, CLEAN4 hook ρ 0.95.
  - LAB-LOG l.167 and l.184: it has one miss, idall3. Identity- and programme-type sign priors are over-read and must carry
    a range spanning zero.

### 4. Does it use seen data as evidence of gain? **No.**
- The sources are H1-free: VCC25 H1 and Pan's H1 are never read (812 RESULTS l.212).
- No A/B/C board read and no sealed reserve is used.
- DE-H1 is worn, but it is not a training line of this source. This check passes. The closure in check 1 is what kills the
  proposal.

### 5. Is the try budget honest? **No, under-counted.**
- RESULTS l.207 records "tries 4": three public-mean variants plus the bulk sign. The proposal uses e(3) = 0.85 (l.18).
  - At N = 4, e(4) ≈ 1.03.
  - The selection was also over classes: u1 was picked among a0 / a1 / a1head / u1, two comparators each. So the effective N
    is larger than 4.
- The proposal's own conclusion still stands, and it is stronger at the honest N: against CLEAN4 the u1 gain (+0.0013
  [+0.00007, +0.0025]) does not survive best-of-N (PROPOSAL-frozen l.19).
- Against the incumbent, +0.0022 is mostly the CPM 5–20 carrier-sign replacement covered in check 1. That is the part 430
  and 417c closed.

### 6. Is the expected gain worth its cost? **No.**
- Stage 1's own prediction, Δnmae −0.002 … −0.004 raw, is +0.0005 … +0.0011 overall-eq at H1's rates.
  - Apply the H1 → board factor for sign priors (0.07–0.67×, FINDINGS l.722). Then apply 516's fold-predicted one third of H1
    (516 l.49–50).
  - The expected board value is ≈ +0.0001 … +0.0004. That is at or below the readable floor of a prior-only DE arm, and inside
    the 0.0003 bias that best-of-4 adds at the re-roll sd (AGENTS 11.1c.4).
  - The checker would have to report it with a range spanning zero (check 3).
- The cash cost is small (< $2 for stage 0, ~$5 for stage 1). The wear cost is not: 12 more emitter runs on DE-H1, which is
  already worn, for a family the record has closed four times.

---

## What would reopen it (for the log; this is not a REVISE path)
A **new mechanism**, not a new bank list. Two examples:
- a source whose H1 u1 gain survives on the essential-free H1 subset and after per-line centring;
- a held-out DE truth (K562-GWPS fold, or 805's DE-CD4 if the bank excludes CD4) that is built and pinned by path + sha256
  **before** any read.

It should reach the lab as one of lane 811's model outputs scored by skill812 under 811's GO, not as a separate 812 compute
line.

## Recorded closure (for docs/LAB-LOG.md, lab lead's hand)
**H-U1PUB (812): KILLED (judge, 2026-10-16 08:36Z).**
- What was proposed: the public per-target mean of 12 non-H1 lines as the uncalled-gene hook sign.
- Why it is closed:
  - it is 516's closed rider, and the plain-mean form lost to the analog on K562 (516 l.204);
  - its gain over the incumbent is unwhitened donor-mean signs replacing whitened carrier signs at CPM 5–20 (430 / 417c);
  - its H1 trigger includes the essential-donor banks that cover 0 / 300 of an Arc panel (430 l.26–28);
  - half the gain is a line programme (812 RESULTS l.239; FINDINGS l.671–677).
- Tries charged to the family: 4 sources × the class table (812 RESULTS l.207–228).
