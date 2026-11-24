# 317 — Feasibility: can OUR pipeline trade `mse` for `nmae` + `reach`, and where does the SIX-MEMBER composite peak? The operating point is the variable, not the arm.

You are a research subagent. Work in `experiments/317-operating-point/`. **Edit nothing outside it.**
**NEVER run `vcc submit`.** You may build predictions **only inside your own work directory** — nothing under
`submissions/`, and do not edit `scripts/make_submission_pm_cd4.py` or `scripts/def_pipeline.py` (a submission
built from them is armed for 00:00Z and the files must not move). `uv run python`. **One `cell-eval2` instance at
a time, ever** (25 GB RSS, the machine is shared). Peak RSS under ~24 GB and only one heavy job at a time.
**Time box: 6 hours.** A well-supported NO-GO is a valid result.

## Why — a standing conclusion of ours is contradicted by the live board
`docs/READ-FIRST.md` §4.1c concluded "**there is ONE route, not two: transfer cosine**", and 227 priced the `nmae`
credit at +5.0 % relative aggregate cosine on h1-test / +12.3 % on PLOLO against a ±0.2 % gate — "failed by
25–60×". **The new board leader falsifies the premise that pricing rests on.** `ViRtUaL cElL builders`,
**0.240185**, first of 778, with the inversion **in domain** (`score_mse` +0.07252 > 0):

| raw | us | leader |
|---|---:|---:|
| `expr_mse_unbiased_capped_norm` | 0.81671 | **0.99153** (near the α=1 baseline) |
| `de_wilcoxon_lfc_nmae` | 0.93795 | **0.87166** |
| `de_wilcoxon_direction_reach_raw` | 0.27428 | **0.37115** |

Scaled, leader − us: `mse` **−0.10761**, `nmae` **+0.11257**, `reach` **+0.10990**, `pds` +0.02657 → total
**+0.02247**. **Their implied transfer cosine is 0.01546 against our 0.09973 — six times worse — and they beat us
by +0.0225.** 227's pricing assumed we must KEEP the cosine. This entry gave it away and bought the DE members.

**Treat that as a HYPOTHESIS about what is reachable, never as a target to imitate** — it is one public entry with
a self-reported model name, and §4.1c already carries a *retracted* amendment from exactly this kind of cross-team
inference. **Your job is to measure whether OUR pipeline has such a trade at all, on our own instruments.**

## READ FIRST
`docs/READ-FIRST.md` §1.0, §1.1, §4 (instruments) and §5 (closed families); `AGENTS.md` §4 (the six metrics and how
each is scaled), §4.1/4.1b/**4.1c**, §6.4 (local scoring), §8.1c–8.1f (the ordering channel and why it was closed
"by arithmetic"). Then the two experiments that already touched this axis and must not be re-run blind:
**034-emission-variance** (*"`emit_counts(mode='quasi', spread=…)` changes cell-to-cell variance while preserving
gene totals — a DE lever, not an `mse` lever"*; `reach` +14–18 %) and **094-reach-ordering** (*"the term is already
saturated at v11, so growing the call set further cannot pay"*). **Say explicitly whether 094's saturation is a
property of the metric or of our construction — the leader's `reach` is 0.371 against our 0.274, so it is not
saturated in absolute terms.**

## Gate 1 — the instrument (do this first and report it before any sweep)
**You cannot answer this question with `pds` alone; you need all six members.** Establish what scores six members
locally: `scripts/cloud_score.sh` / `cloud_score_worker.sh` (the coordinator's note says the real six-member
pipeline is ≈ $1.50 per run on AWS and `prep-real-bundle --baseline` wants the prediction h5ad), `cell-eval2`
locally (ONE instance), and the h1-test harness (R273 reproduces the scorer's own `pds` to 5.6e-17 — establish
whether the other five are available there too). **State the instrument, its cost per evaluation, and its
known-answer row** (reproduce a known six-member row — `cd4svd5`'s or `312-normrel15`'s — before trusting it).
**If no six-member instrument is available within budget, say so and stop: that is the answer.**

## Gate 2 — the sweep, with the OPERATING POINT as the variable
Sweep the emission/calibration knobs that 034 and 094 identified, on a FIXED prediction (do not change the
model): `--overdispersion` (shipped 4.0), the quasi `spread`, `--tail-spread` (1.0), `--keep-top` (8000),
`--keep-spread` (0.25), and `VCC_ALPHA` amplitude (0.982 — 211 closed it *for cosine*, which is not the same
question). Report **all six members at every point**, plus the composite. Name the experimental unit for every
interval.

## Gate 3 — the question
Is there an operating point where the composite **beats our current 0.2177** by more than the emitter noise floor
(034 measured sd 0.0004 / 0.0012)? Report the peak, the trade it makes (which members it pays and which it buys),
and how far it is from the leader's shape. **If the composite is flat or falls everywhere, that is a clean NO-GO
and it re-establishes §4.1c on stronger evidence than it had.**

## Rules
Programme rules 1–9 bind. Strong simple baselines: the shipped point, and each knob alone. **No model change, no
new arm, nothing competition-facing selected or tuned** — this is a measurement of a surface we already ship.
Carry a known-answer row for every reproduction. Report cost if you spend cloud money, and **cap total spend at
$25** — stop and report if a sweep would exceed it.

## Deliverable
`RESULTS.md`: Gate 1's instrument with its known-answer row and per-evaluation cost; the six-member sweep table;
the peak and its trade; whether 094's saturation is metric or construction; and one line —
**317 (operating-point trade `mse` → `nmae`+`reach`): EXISTS / DOES NOT EXIST (peak composite … vs 0.2177 shipped;
it pays … and buys …)**. Then stop.
