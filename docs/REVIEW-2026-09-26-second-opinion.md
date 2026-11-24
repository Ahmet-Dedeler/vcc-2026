# Second-opinion review, 2026-09-26 ~20:10Z

A second Claude Code session reviewed the coordinator session ("Resume previous chat session", Opus 5.5, max effort,
50 % context) at Ahmet's request. The reviewer made no changes except this file: no builds, no git, no edits to live
scripts. Scratch analysis lives in the reviewer's scratchpad; every number below can be re-derived with the commands
given.

## 1. Verdict

The coordinator is not degraded. Its recent work (468's known-answer rebuilds, the gate verifiers, the three-way split,
synthetic dry runs of every slot rule) is as careful as anything earlier in the record. It catches most of its own
slips (waiter patterns matching quoted `END` lines, timestamps written from memory, a few over-claims later corrected).
The blind spots are systematic rather than new:
1. Recent H1 claims (467, C-71) carry seed noise only, not target-sampling noise, though 423 did report both (§3).
2. Constants tuned on A/B/C are carried to October without checking what they mean on an October-like panel (§2).
3. (Compliance of the same-line rows: settled by Ahmet, allowed; see §4.)

## 2. FINDING (October): the prior gate's q33 is a different rule on an October-like panel

`VCC_DE_PRIOR_GATE=q33` (builder 693fde15) marks a target weak when `strength < np.quantile(strength, 0.33)`.
Strength is a coarse integer (count of gated genes with |design row| > 0.2), so the quantile jumps between classes.

- On H1 (467), 28 of 85 targets (32.9 %) have strength ≤ 1: the q33 threshold is 1.72, a hair from 2
  (470's log prints `threshold 1.72`). Any q below 0.33 protects only the 5 strength-0 targets. All four reach-loss
  targets (CNN3, COL6A1, ZNF664, AASS) have strength 1.
- On A/B/C (458g, what the 09-30 pair tests) the threshold is 2 and the weak set is {0, 1}: 86/83/83 targets.
- On the efr synthetic October panel the q33 threshold drops to 1. Recomputed from the stack dump's `carriers`
  (`465/work/run/dump_h1x2_cd4x1.npz`) with each context's control gate (min_cpm 10):

| panel | strength ≤ 1 | strength = 0 | q33 threshold | weak set | strength-1 targets protected |
|---|---|---|---|---|---|
| A/B/C relabelled (453 dump) | 25–28 % | 21–22 | 2 | {0,1}: 76–83 | yes |
| efr synthetic (465 dump) | 54–58 % | 68–73 | 1 | {0}: 68–73 | **no (94–100 per context)** |

  The recomputation reproduces the builder's A/B/C counts to within about 5 % (the builder rescales rows after the
  dump), which is far smaller than the efr shift.
- So `DEF_DAY_PRIOR_GATE=q33` on a panel like efr leaves exactly the strength-1 class, the source of all of H1's reach
  loss, on the programme. If a panel had ≥ 33 % strength-0 targets (efr is at 23–24 %), the threshold would be 0 and
  the gate would protect nobody.
- Suggested fix: an absolute form, e.g. `VCC_DE_PRIOR_GATE=le1` (strength ≤ 1 keeps CLEAN4). On A/B/C it selects the
  same weak set as q33, so it should rebuild 458g byte for byte: a free known answer, and the 09-30 read then
  validates the rule that October actually runs. Whether ~55 % CLEAN4 is right on a weak-design panel is a separate
  question worth one H1/efr look. Not urgent tonight: 09-30 and 10-01 run on A/B/C, where the threshold is 2.
- Do not touch the builder while `chain471.sh` runs (it pins, patches and restores it).

## 3. FINDING (method): the gate's H1 evidence carries target-sampling noise the RESULTS do not report

467 reports the gate's paired gain as +0.0043 with seed sd 0.0002. Seed sd measures emission re-rolls. For
October the question is target sampling. A target bootstrap over h1-test's 85 targets (467's per-target CSVs,
4,000 resamples, 467's spans):

| H1 arm vs CLEAN4 (overall-equiv) | point | 95 % target-bootstrap CI | P(≤ 0) |
|---|---|---|---|
| Arc programme, ungated (462l) | +0.0067 | [−0.0037, +0.0161] | 0.08–0.10 |
| Arc programme, gated q33 | +0.0110 | [+0.0067, +0.0160] | 0.000 |
| gate's own benefit (gated − ungated) | +0.0043 | [−0.0033, +0.0134] | 0.13–0.15 |

- Same picture on seeds 0, 1, 2. Without the four reach-loss targets the order reverses: ungated +0.0153, gated +0.0115.
- Reading: the gated programme is a robust H1 gain; the gate's own benefit is not established on H1. It rests on four
  targets, two with a single DE gene (CNN3, COL6A1: n_real 1–2, so reach is one coin flip each).
- The 09-30 pair measures the gate on ~250 weak (target, context) pairs, far better powered. Scaling the H1 bootstrap,
  its generalisation sd to a new panel is about 0.001 overall, above the rule's ±0.0005 bins. That is acceptable for a
  low-stakes switch, but claim C-71 and 467's RESULTS should carry the target interval, not the seed sd.
- 423 already reported target-bootstrap intervals; 467 did not. Suggest making it standard beside every H1 six-member
  claim, above all for reach.

## 4. Compliance: settled

Ahmet ruled (2026-09-26) that the same-line rows and the ballast cells are allowed. No action.

## 5. Strategy note: nmae is the largest gap to #2–#3

Board 2026-09-26 19:44Z, scaled: us pds 0.815 / mse 0.365 / nmae 0.173 / fid 0.006 / reach 0.227 / jac 0.028;
cqawesome 0.827 / 0.402 / 0.287 / 0.085 / 0.291 / 0.018; Testh 0.801 / 0.408 / 0.249 / 0.053 / 0.282 / 0.016.
- Against cqawesome the member gaps are nmae 0.114, fid 0.079, reach 0.064, mse 0.037, pds 0.012 (jac +0.010 ours).
  nmae alone is ~0.019 overall.
- "chill" posts raw nmae 0.835 with a much worse bulk than ours (mse 0.749 vs 0.641, pds 0.844 vs 0.867). Their DE
  channel beats ours with less bulk information, which fits the record's reading (gate sign accuracy, 434f) and says
  the DE channel, not the bulk, is where the next points are. 458d/458g are the right tests.
- Question for the coordinator: nmae reads arithmetic CPM means, the Wilcoxon members read ranks. The record already
  exploits this: `VCC_DE_MAG` (423/429, exact-U pass with "mass-holder" carriers) moves means with every U unchanged,
  and the 09-28 zs riders put `dn_k250_zs` on the board. (The reviewer's first grep missed it and wrongly said nothing
  was in the record.) The up-gene selections were measured and closed (423: up_k140 costs nmae +0.0005 on H1). So this
  question is answered; nmae's remaining gap is signs on the gate genes, which is what 458d/458g test.

## 6. Automation audit (all PASS at 19:50–20:05Z)

- 09-27 waiters (pids 50989/50990): both packages match their build-log sha256 (446bj2 `c7206c46…`,
  446br2 `9bbc8d2a…`), both novel, `vcc whoami` can_submit true, NOGO files absent, clock caps 12:00Z.
- arm0928/0929/0930 are exact copies of arm0927 except dates and paths. chain0928, chain0929, chain0930c and
  chain1001b wait on lines only their predecessor writes last, plus process exit, with sane time caps. rule1001b.py
  matches its docstring.
- Readiness of every candidate package through 10-01 (package present, PASS, CLIFF OK, novelty rc 0 in its log, sha
  line, novelty check now): 446bz, 446bj2z, 446br2z, 446bd2, 453dr, 457dr, 458d, 458g, 458, 459dr, 461dr, 460dr, 458e,
  469dr all ready.
- 453's relabelled bundle is A→D, B→E, C→F by content (control cells match exactly), so 453dr/457dr/459dr/460dr/
  461dr/469dr map back correctly.
- 468's known answers are in the logs: KA-1 `9559923a` = 446b, KA-2 `1d09cc91` = 458d, builder restored each time.
- Loose ends: `chain471.sh` ends with the builder patched but uncommitted by design (commit after its KAs, or restore
  693fde15). 41 local commits are not yet pushed to the private backup. Caffeinate runs to ~10-01 12:00Z; the Mac must
  stay on power with the lid open through the slots.
