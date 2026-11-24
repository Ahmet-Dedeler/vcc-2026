# 308 — Competition track, feasibility: weight each donor row by the ON-TARGET KNOCKDOWN ITS OWN CONSTRUCT ACHIEVED. Reliability-based weighting is structurally blind to the dominant error term (F14/293: a FAILED construct still yields a REPRODUCIBLE transcriptome); measured knockdown is the build-time-observable discriminator that separates them.

You are a research subagent. Work in `experiments/308-kd-weighted-donors/`. **Edit nothing outside it.**
**NEVER build or write a submission package; nothing under `submissions/`; no `vcc`; no `cell-eval2`; no download.**
`uv run python`. Heavy work serial; peak RSS under ~10 GB. **Time box: 4 hours.** Siblings (306, 307) run in
parallel; do not read or write their directories. This is a MEASUREMENT. A clearing arm becomes a preregistered
build mandate later — you do not build it.

## The mechanism, stated before the test
The incumbent `cd4svd5` is an **equal-weight** mean of three donor screens (K562-GWPS, Orion HCT116, Orion HEK293T),
whitened, plus CD4 and H1 residuals (`docs/READ-FIRST.md` §2). F14 (285, `docs/FINDINGS-2026-09.md`) measures the
**guide construct as the dominant irreducible term** of a perturbation response — larger than a different study in
the same line. 296 shows that term is largely **guide efficacy** in immortalised CRISPRi lines (K562-GWPS both-≥cut
ρ_guide +0.340 [+0.113, +0.548] vs weaker-<cut +0.001 [−0.082, +0.105]; Nadig Jurkat +0.458 vs +0.096), with library
failure rates of **19–27 %**. 293 measures the consequence directly: where both constructs work ρ_guide **+0.687**,
where the weaker fails **−0.003** — and, decisively for this mandate, **a failing construct still gives a
reproducible transcriptome**.

**Therefore:** every donor-weighting scheme tried so far keys on *reliability*, *norm* or *donor–donor agreement*
(016 MODZ — does not beat fixed weights; 097 route (a) `normrel^−1.5` — +0.0009; 262/263's arbiter scores; 299's
reliability-selected rows — does not clear). All of them are **blind by construction** to a construct that failed to
knock down its target, because that construct's response is reproducible, well-normed and self-consistent — just not
about the target. The discriminator that is NOT blind is the donor's **own measured on-target knockdown**.

**Note the distinction from a closed result.** 097 route (d) ("exploit the on-target knockdown") is dead *as a
prediction target* — the scorer drops each row's own-target coordinate from every member (`_drop_on_target`). Using
|kd| as a **weight on which donor to believe** is a different use of the same quantity and is not closed. State this
distinction in your RESULTS.md and verify the scorer exclusion claim yourself.

## READ FIRST
`docs/READ-FIRST.md` (all; §2 the shipped configuration, §4 instruments, §5 closed families — **do not re-run any
closed family**), `AGENTS.md` §4 (scoring), §8.1–8.4 (the pricing rules; the ≥ +0.0089 pessimistic gate and the class
factors k_pds 1.22 / k_cos 3.76 from BR249), `docs/FINDINGS-2026-09.md` F1, F4, F14, and §3.
`experiments/293-challenge-construct-ceiling/RESULTS.md` (§ on-target knockdown of every target in the three splits —
**it already computed |kd| per target; reuse it and reproduce it as a known-answer row**),
`experiments/296-guide-efficacy/RESULTS.md` (the cuts, the library-relative vs absolute thresholds, the null
correction), `experiments/016-modz/RESULTS.md` and `experiments/097-direction/RESULTS.md` routes (a) and (d) (the
weighting baselines you must beat and the closure you must not re-run), `experiments/262-*`/`263-*` (the PLOLO fold
caches, `plolo.build_fold`, `E062.prepare_main_source`, `build_baseline`; the per-target `pds` arrays and their
known-answer values 0.815384 / 0.801659 / 0.780978 / 0.640294), `experiments/273-*` (the h1-test instrument).

## The feature (build-time observable — prove it)
For every panel target t and every donor screen d, the **on-target knockdown** |kd|(t,d): the target's own gene
coordinate in donor d's own measured response to perturbing t, in d's own units, computed from d's data only.
- **It must not touch any held truth.** Prove this: state the provenance of every array, and carry a leakage row
  showing the feature is identical whether or not the held fold's data is present.
- Report its distribution per donor, its missingness, and how many targets have a **disagreement** between donors
  (one donor's construct working, another's failing) — that disagreeing set is where the arm can act at all. If that
  set is small, say so early; it bounds everything below.
- Use 296's library-relative cut AND 293's absolute cut; report both. Depth-condition per programme rule 9.

## The arms (all on top of the incumbent, one variable each)
1. **Hard gate:** drop donor d for target t when |kd|(t,d) is below the cut; renormalise over the survivors; if all
   three fail, fall back to the equal-weight mean (state the fallback rate).
2. **Soft weight:** w(t,d) ∝ a monotone function of |kd|(t,d) (state it, one parameter, and sweep it coarsely).
3. **The same, restricted to the disagreeing set** (act only where donors disagree; keep the incumbent elsewhere).

## Controls — every arm must beat all of these, and this is the point of the experiment
- **The incumbent** (equal weight), reproduced bit-identically first (KA row).
- **Reliability weights** in 016's and 097(a)'s forms — *the decisive control*: if |kd| adds nothing over
  reliability, the mechanism above is refuted and you say so plainly.
- **Shuffled |kd|** across targets within donor (3 seeds), and **random weights matched to the arm's weight entropy**
  (3 seeds), and **reversed |kd|** (weight the FAILING construct more — must lose if the mechanism is real).
- **Norm-matched output** so no gain is magnitude (AGENTS 8.1d: signs pay, magnitudes do not).

## Instruments and the bar
- **PLOLO** (three cancer folds; `pds`, `cos_agg`, `nmae`) — an arm must clear **3/3 folds** with a paired-target
  bootstrap interval clear of zero, and every control at or below half the arm.
- **h1-test** (the Flex-400 reference, R273 — it resolves every target, so read the transfer cosine there too) with
  an interval.
- **Price it** by AGENTS §8.4's rule (the class factors, the five-member propagation, the coherent-pessimistic
  convention) against the **+0.0089** pessimistic gate. Report the price BEFORE the verdict sentence.
- **`nmae` is the killer** (ARM6 died there, BR237): report the `nmae` leg of every arm explicitly. An arm that buys
  `pds` and loses `nmae` is not a candidate.

## Rules
Programme rules 1–9 bind (name the unit; strong simple baselines; a failed feasibility test closes its construction
at its scope, not the field; multiplicity; depth-conditioning). Known-answer rows before anything else: the fold
caches reproduce 262/263's `pds` arrays at 0.0e+00; 293's |kd| values reproduce; the scorer's own-target exclusion is
what 097 says it is. **A well-supported negative is a valid and expected result** — if reliability-weighting already
captures everything |kd| carries, that is a finding about F14's term and it belongs in FINDINGS.

## Deliverable
`RESULTS.md`: the KA rows; the feature's distribution and the size of the disagreeing set; each arm × each control ×
each fold with intervals; the h1-test reading; the `nmae` leg; the price against the gate; and one line —
**308 (|kd|-weighted donors): CLEARS / DOES NOT CLEAR (best arm … `pds` [lo, hi] on n/3 folds; vs reliability
control …; `nmae` …; priced … vs +0.0089).** Plus `work/kd308.csv` (per target × donor: |kd|, cut calls, weights,
per-target `pds` under each arm). No build, no submission, no production file touched. Then stop.
