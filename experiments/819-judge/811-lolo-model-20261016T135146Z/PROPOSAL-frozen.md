# 811 PROPOSAL-3 — sign skill on responding genes: the honest answer is NO construction; one optional diagnostic

FROZEN. Written 2026-10-16 13:48–13:51Z (clock) by the 811 worker for the lab lead's round of 13:47Z. Nothing was computed for
it; every number below is a lookup with its source. Nobody edits this file after it is written; its sha256 is recorded in
`work/SHA256SUMS.txt`.

## 0. Answer
- **I do not propose a new sign-skill construction.** Every per-target public-line sign source my lane could build has
  already been priced and closed:
  - plain mean, Stouffer, gated by |Z|, centred, with CD4 at three votes, conservation-gated, an oracle per-target choice,
    and a supervised stack (516, 621, 525, 634);
  - 812's 12-line mean, judged on 10-16 (KILL, `812/VERDICT-u1hook.md`).
- My lane's one positive finding, the many-line increment, works through **target discrimination** (pds). Under the
  scorer, DE signs are mostly **generic** (versus control). §2 explains why the first does not turn into the second, and
  the measured fact that it does not (803 §8: Δfid ≤ 4e-5).
- **What I offer** is one optional, pre-counted, $0.10 diagnostic (§4), for the judge to accept or decline. It is not a
  candidate.
- **c0_l25 stays the frozen 6a candidate, unchanged (§5).**

## 1. What the record closed (the constructions this lane could have offered)
| construction | where closed | the number |
|---|---|---|
| Stouffer consensus of same-target public rows (K562, HCT116, HEK293T, CD4 × 3) on uncalled genes, \|Z\| ≥ 0.5 / ≥ 1, centred variant, CD4 weighting, conservation gate | 516 (`516/RESULTS.md` l.10–45, 140–160, 198–227) | H1 +0.0157 (essential-donor and unwhitened-carrier artefacts); cancer folds: HCT116 +0.0033, HEK293T −0.0002, K562 +0.0010 [−0.0030, +0.0053] vs the production analog. "Not a lever; at most a rider" |
| a supervised (chill-style) sign model across public lines | 621 (FINDINGS l.266–269) | loses to CLEAN4 held out (h1 −0.013, hct116 −0.099, hek293t −0.205, k562 +0.001) |
| a confident family-balanced consensus (\|Zbal\| ≥ 2) | 621 (same lines) | 0.81–1.00 vs CLEAN4 0.60–0.73 on its few cells, ≈ +0.001 overall: the confidence gate is right but covers too little |
| unwhitened donor-mean / blended carrier signs | 430 (FINDINGS l.889) | cross-study +0.001; the H1 gain is the essential-donor artefact |
| 812's per-target mean of 12 non-H1 lines (raw / RMS / centred) as the u1 hook source | 812 RESULTS §7; KILL in `812/VERDICT-u1hook.md` | u1 +0.0013 vs CLEAN4 (+0.0022 vs the incumbent); a0 −0.0025; centring halves u1 and deepens a0 |
| agreement / reliability weighting of lines | 016, 142, 271 (`VERDICT-1.md` §1) | equal weights win |
| lab programme as every context's prior (Arc's H1) | 458d board read (FINDINGS l.322) | −0.0085 (HARM); the gated form is still −0.0049 vs 446b |
| gene-level CLEAN4 composition (subsets / supersets) | 613f/g (FINDINGS §2) | equal 4-bank CLEAN4 wins |

My own negatives that bear on signs:
- per-line centring keeps 89 % of the bulk increment, but 812 found centring halves the u1 sign gain (above);
- the essential-screen effect: the A-only panels removed the essential targets and the direction over-read stayed
  (LAB-LOG 13:21Z). The essential-donor artefact that inflates H1 sign reads (430) is a separate thing, and it applies to
  any H1 read of a bank set that includes RPE1 / Jurkat / K562-ess.

## 2. Mechanism: why the many-line gain is not a sign-skill gain
- **What the bulk increment is.**
  - On the Arc-drawn panels it holds (LABC +0.0036, 3/3; 062 +0.0031, 3/3). On A-only it is +0.0019, 14/18, carried by pds
    (+0.0045) with r +0.45 % (RESULTS §00000).
  - Per target, it lives on held K562 folds: R_tgt +0.0015 … +0.0034. On HCT116 / HEK293T it is ≈ 0 (RESULTS §00).
  - So the extra lines sharpen what makes target t different from the panel's other targets. That is pds's quantity.
- **What a DE sign on a responding gene is.** It is the sign of the gene's log fold change against the control pool.
  812 §2 shows the movable part:
  - a0 on weak targets is 0.50 for every source;
  - what moves lives on strong targets and near-significant genes, where the GENERIC down-programme (CLEAN4's 84 % down
    lean) already carries most of the information.
  - A target-specific direction (the bulk increment) mostly adds the parts that differ between targets. On a responding
    gene, those parts rarely flip the sign that the generic programme and the composition already give.
- **The measured consequence.** 803 §8 put c0_l25 (the many-line blend) through the real emitter on both H1 surfaces:
  - fid and jac moved ≤ 4e-5 raw and reach stayed inside seed noise, so called-gene signs did not change;
  - nmae improved by −0.0004 … −0.0010 raw (magnitude, not sign).
  So the bulk route delivers no sign skill. A sign route needs a per-gene sign SOURCE, and every source in §1 is closed.
- **What would change this** is new information, not a new combination:
  - a public line close to D/E/F's lineages, at depth;
  - or Arc-lab Flex data on non-H1 lines.
  - TeloHAEC genome-wide failed 813's truth-line gate (lane parked). KOLF2 full rows are pluripotent, like H1, so H1 sign
    reads would flatter them (702 rule 2: a lineage-matched surface cannot supply the gain). No such data is on disk.

## 3. Instrument limits (why I do not ask for a DE-H1 price of anything new)
- skill812 is an H1 instrument. It needs H1-free predictions and reports essential-free subsets only if asked.
- Two of its known inflations apply to my lane's sources. Any H1 read of the extra lines inherits:
  1. the essential-donor artefact: RPE1 / Jurkat cover 21–36 H1 targets and 0 A/B/C targets (430, 812 §7);
  2. the unwhitened-carrier effect: the board inverted it (417c, nmae +0.0089 raw).
- The per-class rates are H1's; board transfer of sign priors on fid is 0.07–0.67×.

## 4. Optional diagnostic (D3), pre-counted; accept or decline
**Question.** Does the DAY-list consensus that c0_l25 actually ships carry any sign skill on responding genes beyond
CLEAN4 and the incumbent? This needs its exact U, H1 group held out, read on the essential-free targets. It is the one
H1 read of c0's own sign content that does not exist yet:
- 812's bank had no CD4;
- 803 priced c0's DE members, not its signs by class.

**Construction (frozen, no new arithmetic).**
- U = `m811.day_consensus(panel = the H1 surface's targets, centre=False, hold={"H1"})` on the H1 targets.
- The banks are 803's `bt/ext803.py` extraction for the H1 targets (KA-passed against 811's manifest). Its CD4 rows cover
  only 5 / 11 targets, because they came from panel compacts.
- If the judge wants CD4 at full coverage, ext811a's CD4 path (047's extractor) runs on the H1 targets first. Its KA:
  re-extracting L1 must be byte-identical to 801's L1 compact.
- pred = U's sign on the 5-CPM axis. Scored with `skill812.py --h1-free yes`.

**Reports.** Pooled and per-surface a0 / a1 / a1head / u1, vs CLEAN4 and vs the incumbent, on all targets and on
`no_ess` (the targets no essential screen covers; 516 used 64 on h1-test).

**Prediction, with bins, before any number** (pooled overall-eq, vs the incumbent, on `no_ess`):
- u1 in [−0.0005, +0.0015];
- a1 in [−0.002, +0.0005];
- a0 below CLEAN4 (as 812's pubmean: −0.002 … −0.003);
- a1 + u1 combined < +0.002.

**Bar (812's own) and decision.**
- Pooled u1 + a1 vs the incumbent ≥ +0.002 with lower bound > 0 on `no_ess`, and a0 not below CLEAN4 by more than its CI.
- **Below the bar:** the sign-skill objective is closed for this lane for October, on c0's own U. I expect this outcome.
- **At or above the bar:** the read goes back to the judge as a hook-source proposal with a held-out cancer-line stage.
  It does NOT change c0_l25.

**Tries and cost.**
- Family `811-signs`, N = 1 on DE-H1 (skill812). This family has not been priced on DE-H1 before; 812's four
  source reads are 812's family.
- Cost: one r7i.xlarge job, ≈ $0.05–0.10 (≈ $0.3 with the CD4 extraction).
- Class: transfer, if it passed.

## 5. Relation to c0_l25
- **c0_l25 stays the lab's frozen 6a candidate exactly as it is** (`m811.py` `e00192e6…`, day list, λ 0.25, C off). No
  swap.
- Its standing evidence:
  - bulk on the Arc-drawn panels: ABCp +0.00264 [−0.00003, +0.00561] 3/3, LABC +0.00230 3/3, both WEAK after e(N);
  - DE on DE-H1 through 803's hook: no cost, nmae −0.0004 … −0.0010 raw;
  - its one read on D/E/F-panel folds (N = 1 on DEF) decides it.
- Nothing in this proposal is a rider on it. D3 is a diagnostic of c0's U, not a change to it.

## 6. Wear (family `811-consensus`, from `work/wear811.jsonl`, sha256 `82663b4b…`)
| surface set | N |
|---|---:|
| W-tune (L1–L3) | 10 |
| W-conf (L4–L6) | 7 |
| BULK-ABC 062 | 9 |
| ABCp | 10 |
| LABC | 11 |
| A-only | 1 |

The lab's ledger lists c0_l25 at N 11 on the Arc-drawn set. D3 opens a separate family on DE-H1 (N = 1).
