# 318 — Feasibility: is PER-TARGET magnitude redistribution a real lever, separate from the uniform amplitude axis 8.1d closed? Isolate it from the weighting that is confounded with it.

You are a research subagent. Work in `experiments/318-per-target-magnitude/`. **Edit nothing outside it.**
**NEVER run `vcc submit`; never build into `submissions/`; do not edit `scripts/make_submission_pm_cd4.py` or
`scripts/def_pipeline.py`** (they carry two shipped arms and a third session may build from them). `uv run python`.
**One `cell-eval2` instance at a time, ever.** You MAY use `scripts/cloud_score.sh` — **cap total spend at $12**
and report it. **Time box: 5 hours.** A well-supported NO-GO is a valid result and closes the axis.

## Why — three board points now disagree with a closed family, and the disagreement is specific
| entry | what it is | overall |
|---|---|---:|
| `Xv0ayvuoHpsKYWv2iu42` `cd4svd5` | incumbent | 0.2160944 |
| `Q519I2pLWRCW7PAEx5ul` | per-target weighting **+ an accidental per-target norm distortion** | **0.2177189** |
| `2K5ppxkowMSj2ugoZqYb` | the **same weighting**, made exactly shape-only (313/316) | 0.2169733 |

So the weighting alone is worth **+0.00088** and **the distortion adds +0.00075 on top**. The distortion's
delivered per-target row norms run **0.806–2.155** (whole matrix +6.0 %) where the parity arm holds every row at
**1.0000000**. Member-wise the distortion **bought** `mse` +0.00399 and `pds` +0.00225 and **paid** `nmae` −0.00339.

**The reason this is not already closed.** 317 swept **uniform** amplitude on the real six-member scorer and found
×1.20 *pays* `mse` −0.038293 — the **opposite sign**. AGENTS 8.1d ("signs pay, magnitudes do not") and 211's
amplitude closure at 0.982 were both established on the **uniform** axis. **A per-target magnitude rule is a
different operation and has never been tested.** 317 also measured, over 1,961 consecutive same-team board pairs,
that **buying `mse` is what wins here** (mean +0.01538, 83.2 % of teams improve), which is what the distortion did.

## The confound you must break
The shipped arm's distortion is **not** a clean per-target magnitude effect: that arm also applies the
PerturbMatch correction at the **uniform** `k562_share` and outside the norm match (313 §1). So its norm profile
and its effective weighting differ from parity's **together**. **Your first job is to separate them.**

## What to do
1. **Characterise the distortion, from saved artifacts** — no rebuild needed first. Recover the per-target
   delivered norm ratio for the shipped arm from `experiments/316-parity-arm/work/` (its `test316.py` already
   computes both profiles on 313's fixture) and from the two builds' metas. **What predicts a target's norm
   ratio?** Its donor disagreement, its weight entropy, its `pds`, its reliability. State whether the profile is
   structured or effectively arbitrary — **this decides whether there is a rule to learn at all.**
2. **Construct the isolated arm.** A per-target magnitude multiplier `m(t)` applied to the **parity** (exactly
   shape-only) prediction, so the weighting is held fixed and ONLY the magnitude profile varies. Two forms, both
   declared before scoring: (a) **replay** — `m(t)` = the shipped arm's own measured norm ratio, which asks "was
   that particular profile the thing that paid?"; (b) **rule** — `m(t)` monotone in a build-time observable you
   justify from step 1. If step 1 says the profile is arbitrary, run (a) only and say so.
3. **Score on the six-member instrument** (317's path: `scripts/cloud_score.sh`, ≈ $1.56 and ≈ 2.8 h per point,
   workers self-terminate). **Read 317 §1.2–1.4 first and reuse its KA rows** — the recorded `rc_callset`
   six-member row must reproduce before you believe any number. **Read 317 §2.2 before quoting any h1-test
   `mse`**: h1-test is scored with the H1 union OFF (oracle-leaked, 85/85), its raw `mse` is 1.016 against the
   board's 0.81671, and **it cannot locate the board's amplitude optimum.** Do not repeat that error.
4. **Price it** against the +0.0089 gate an amplitude-changing arm must clear (AGENTS §8.4), and say plainly
   whether the +0.00075 the distortion bought is reachable deliberately, or was an accident of one profile.

## Rules
Programme rules 1–9. Name the unit. Controls: the parity arm itself (the baseline), a **shuffled** `m(t)` across
targets at matched entropy (3 seeds — if this matches the arm, the profile carries nothing), and a **uniform**
scale matched to the same whole-matrix norm (this is 317's axis and must lose). **Nothing competition-facing is
selected or tuned; you do not build a submission.**

## Deliverable
`RESULTS.md`: what predicts the norm ratio; the isolated arm's six members against parity, with the shuffled and
uniform controls; the price against +0.0089; and one line — **318 (per-target magnitude redistribution): REAL AND
PRICED / NOT SEPARABLE FROM ITS SHUFFLED CONTROL / NO-GO (…)**. Report spend. Then stop.
