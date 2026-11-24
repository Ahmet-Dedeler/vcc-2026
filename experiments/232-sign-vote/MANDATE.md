# 232 — SIGN VOTING across donors. The one construction §8.1d points at that nobody has tried.

You are a research subagent. Work in `experiments/232-sign-vote/`. **No `vcc`, no submission build,
no production change. No `cell-eval2` while another agent runs one — check first.** Peak RSS ≤ 10 GB;
**check `memory_pressure` and `sysctl vm.swapusage` before any job over ~4 GB and hold if swap free
< 2 GB** — 229 is running 20 GiB builds and a job that exceeds memory here is killed **without a
traceback**. `uv run python`. Write only inside your own directory.

## Why this experiment exists

**`AGENTS.md` §8.1d (228, today): a cosine gain earned by reshaping per-gene MAGNITUDES is destroyed
by `nmae`; one earned by improving SIGNS is not.** `nmae`'s matched term is sign-matched
`min(|pred|,|truth|)`, so it **cannot** collect a magnitude rebalance. Four independent arms have now
died on that wall (224's per-gene gain, 158's per-target scale, 169's sharpening, 187's amplitude).

**So the only shippable direction arm is one that moves SIGNS. And there is exactly one construction
in the consensus that has never been examined: how the donors' signs are combined.**

Production predicts `sign(mean of donor effects)`. **That is not the same object as the donors'
majority sign**, and the two differ precisely where it matters — where donors disagree. For a gene
with donor effects `+5, −1, −1`: **the mean says `+`, the majority says `−`.** One large-magnitude
donor can outvote two agreeing ones under a mean, and never under a vote.

**Nothing closed covers this.** 194 swept **eleven estimators** — every one a **per-cell weight
inside a single donor**, a different level entirely. 204 swept **29 donor weightings** and found
equal weighting maximises the transfer cosine — **but every one of the 29 is a weighted MEAN, so all
29 share `sign(weighted mean)`.** 209 closed donor–donor *correlation*. **No experiment has changed
how signs aggregate.**

## READ FIRST
`AGENTS.md` **§8.1, §8.1b, §8.1c, §8.1d** (especially 8.1d — it is the reason for this experiment),
**§4.1c**, §11.1. Then `experiments/228-theta-gain-plolo/RESULTS.md` §0 (the exact `M`/`U`
decomposition and the `e = ΔM/ΔU` efficiency statistic — **reuse that instrument, do not invent
one**), `experiments/204-repression-weighted-mixing/RESULTS.md` (what is closed and why it does not
cover you), and `experiments/162-retrodicting-validator/` (PLOLO, the challenge frame, the only
`pds` instrument validated against a board point).

## Work

### Phase 0 — known-answer rows, mandatory
Reproduce 228's `M`/`U` decomposition and its `e = 0.072` for 224's gain; reproduce our own PLOLO
baseline on the three cancer folds against the shipped `h1union` stack at **Δ = 0.0e+00**; and
reproduce 204's finding that equal weighting maximises the transfer cosine. **If these do not
reproduce, stop.**

### Phase 1 — is the majority sign actually better? (measure before you build)
On PLOLO donor folds, **for each held line, on the significance-selected DE pool**, measure
**|t|-weighted sign accuracy** of:
- `sign(mean)` — the incumbent;
- unweighted **majority sign** across covering donors;
- **confidence-weighted** sign vote (weight by |t|, by per-donor coverage, by 189's `theta` — say
  which you use and why);
- the **oracle** sign (upper bound).

**Report the disagreement rate too** — the fraction of scored (target, gene) pairs where mean and
majority disagree. **If that rate is tiny, or if majority sign is no better than `sign(mean)`, the
experiment ends here and that is a complete result.** Do not proceed to build on a null.

### Phase 2 — the arm, constructed to move signs and NOTHING else
**Keep every magnitude exactly as production computes it; change only the sign** where the vote
disagrees with the mean. State explicitly that `|effect|` is bit-identical to the incumbent's, and
**prove it** — that is what makes this a clean test of §8.1d rather than another magnitude arm.
Sweep the vote's confidence threshold (flip only where the vote is decisive) as the dose.

### Phase 3 — price all six members through PLOLO, challenge frame
Three cancer folds, **`kolf2` reported separately** (pluripotent, gains 3–5×, and it inflated a
coordinator prediction badly enough to force a band reissue on 2026-09-04).
**Lead with `nmae`.** §8.1d predicts a sign arm should move `nmae` the *right* way — that prediction
is falsifiable and it is the point of the experiment. Report `M`, `U` and `e` alongside, so the
result is comparable with 228's line for line.

## Falsifiers, all mandatory
1. **Sign-shuffle null** — flip the same number of signs at random. A gain that survives is not
   measuring the vote.
2. **Wrong-recipient null** — 184 died here (arm +0.0665, control +0.0667, null +0.0663).
3. **Leave-top-3-out retention on cosine and `nmae` separately** (222 §9: 14.6% vs 72.7% on one arm —
   the channels differ).
4. **Magnitude invariance** — prove `|effect|` is unchanged. If it is not, you have built a
   magnitude arm and §8.1d says it is dead; report that.
5. **Coverage honesty** — a gene covered by only one donor has no vote. Report how many pairs have
   2+ and 3+ covering donors, and what the arm does with singletons (**leave them alone**).

## Reporting
**Do not convert PLOLO to board with a scalar** — the factor is **not identified** (1.384 on
`cd4f034`, 1.000 on `h1union`, matched `cos_agg`; and 231 found a 7.84× spread across four pairs).
Quote PLOLO units and name the fold set, functional, gate and tag.

## Deliverable
`RESULTS.md`: the Phase-0 rows, the sign-accuracy table with the disagreement rate, the six-member
price led by `nmae`, the five falsifiers, and a one-line verdict. **"The majority sign is no better
than the sign of the mean" is a complete and valuable result** — it would close the last untested
construction in the consensus and confirm §8.1d's wall is total. **Propose no slot; build nothing.**
