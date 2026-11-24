# 203 — Reopen experiment 125. Its rejection was a pure scale artefact, and its gain was a DIRECTION gain.

You are a research subagent. Work ONLY in `experiments/203-reopen-125/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
**Peak RSS ≤ 10 GB** — five other agents and one submission build are running on this 64 GB box.
`uv run python`.

## READ FIRST
`AGENTS.md` §4.2, **§8.2** (the gap is now entirely direction, and it is sized), §11.1 rules 1–13.
Then, in order:
1. `experiments/198-closure-reaudit/RESULTS.md` — **the diagnosis this mandate rests on. Verify it
   before building anything on top of it.**
2. `experiments/125-hierarchical-panel-posterior/` — `RESULTS.md`, `run.py`, `README.md`,
   `preregistered_protocol.json`. The estimator, its gates, and its frozen artifacts.
3. `experiments/122-panel-response-atlas/RESULTS.md` — the incumbent 125 was scored against
   (`equal_mean_agreement`, macro signed credit 0.029965, top-20% sign precision 0.537409).
4. `experiments/142-batch-reliability-weighting/` — 198 reports it uses the same name
   "signed credit" for the **scale-invariant** version (`develop.py:1379` divides by the
   prediction's L1, not the truth's). That naming collision is the heart of this.
5. `experiments/197-cosine-inversion/RESULTS.md` — the exchange rate you must price against:
   **+0.00365 overall per +0.01 of cosine**.

## The claim under test

125 fitted a strictly nested hierarchical random-effects estimator over the donor panel and was
rejected. Its own result table:

| endpoint | macro delta vs 122 | bootstrap 95% |
|---|---:|---|
| signed credit | **−0.015050** | [−0.020071, −0.010144] |
| top-20% truth-magnitude-weighted sign precision | **+0.012175** | [+0.006403, +0.017913] |

All three held families gained precision; all three lost credit. Its own conclusion: *"The posterior
confidence is useful for ranking signs, but the posterior mean shrinks effect magnitude enough to
lose signed credit."* And, explicitly: *"A future experiment could test a direction-only
confidence/ranking layer while retaining Experiment 122's effect magnitudes, but that is a new
hypothesis and is not authorized by this result."* **You are that experiment.**

**198's diagnosis, which you must verify first.** "Signed credit" is **none of the six scored
members**, and as implemented at `125/run.py:774-776` it is **exactly degree-1 homogeneous in the
prediction** — 198 verified `credit(k·pred)/credit(pred) = k` to **3e-9**. The endpoint that
*improved* is **exactly scale-invariant** (verified: 0.0 difference under rescaling). The estimator
shrinks by **6.78×**. Under a pure-rescale null, the recorded −0.01505 implies
`credit(122) = +0.0177`, which is where a 0.04–0.08 transfer cosine puts it. **So the entire
recorded loss is explained by shrinkage alone**, and shrinkage is precisely what `pds_cosine` cannot
see (it normalises each row) and what `VCC_ALPHA` corrects for free.

**Why it matters now and did not in August.** §8.2: after `rc_hybrid` we are ahead of the leader on
`fid` and `jac`, and `mse` + `pds` + `nmae` account for **110%** of the remaining gap. Everything
left is direction. 184 closed cross-line adaptation on 50 lines, 195 closed learned transforms at the
cell-line boundary, and 197 shows donor count cannot reach the leader. **A per-target hierarchical
shrinkage that borrows strength across targets within the data we already hold is one of the last
two direction levers standing** — and its one measurement improved the direction endpoint by
+0.0122 while being rejected on a scale.

## Work

### Phase 0 — verify the diagnosis, do not inherit it (mandatory, abort if it fails)
1. Reproduce 125's two macro deltas from its saved per-target arrays (its `RESULTS.md` says an
   independent recomputation matched both exactly).
2. Independently confirm the homogeneity claim: `credit(k·pred)/credit(pred) = k` and the precision
   endpoint's invariance, at several `k`. **If "signed credit" is not degree-1 homogeneous, this
   mandate is void and you should say so and stop.**
3. Confirm the 6.78× shrink, and confirm 142 uses the same name for a scale-invariant quantity.

### Phase 1 — re-score 125 on endpoints the competition actually reads
Rescore the frozen 125 candidate against 122 on:
- **`pds` on PLOLO** — 162's four donor folds, challenge frame, pre-committed unweighted mean. This
  is the only `pds` instrument validated against a board point. **`debench.Surface.pds()` must never
  size a `pds` arm** (188 Phase 3, standing rule).
- **the implied transfer cosine**, using 197's route, so the result lands in the units §8.2 is
  written in.
- **|t|-weighted sign accuracy `q_w`** (187 §5: `nmae` moves only through `q_w`).

**Re-solve the amplitude before scoring anything magnitude-sensitive.** 125 shrinks 6.78×; a15's
`VCC_ALPHA` is fitted to the unshrunk scale, and 187 established the joint `mse`+`nmae` optimum is a
sharp ρ cliff. Scoring a 6.78× shrunk prediction at a15's alpha measures the shrink, not the
estimator. Report the arm **at its own re-solved alpha**, and say what the re-solve cost.

### Phase 2 — the direction-only variant 125 itself proposed
Take 122's magnitudes and 125's **posterior ordering/confidence**, i.e. use the hierarchical fit only
where it improved (sign ranking) and not where it lost (scale). Price it the same three ways.

Note what §8.2 implies about where this can and cannot pay: the ordering *key* family is closed
(189, 190 — `reach`'s 0.90 purity bar caps it, and 185/191 already deliver `|H|` at 95% efficiency),
so a confidence layer is **not** worth pursuing for `reach`. It is worth pursuing only if it raises
**cosine** or **`q_w`**. Say which of the two you are measuring and do not conflate them.

### Phase 3 — the six-member price and the gate
Rule 12. `pds` on PLOLO; `mse` through 096's `claim` route (an estimator change moves the effect's
shape, so `mse` is **not** exempt — a19 lost −0.01352 exactly here); `nmae`/`reach`/`fid`/`jac` on
h1-test with a **significance-selected pool** (188: `|lfc|·√cpm` at the reference's BH depth, median
~89 genes/row; a top-N-by-magnitude pool inverts `nmae`'s sign). Donors on h1-test: `k562_gwps` +
`orion_hek293t` only — `data/donors/h1_2025` covers 85/85 of h1-test's targets and is pure leakage.

Rule 13's `reach` toll does **not** apply to a profile-class arm (169 falsified it there; the
`rc_hybrid` board point confirmed the toll belongs to profile-*changing* emission arms). Gate:
+0.0019 shape-only, +0.0089 amplitude-changing — and this arm changes amplitude, so state which.

### Phase 4 — falsifiers
1. **A pure-rescale control**: take 122's prediction, shrink it 6.78×, and score it on every endpoint.
   If that alone reproduces 125's recorded −0.01505, the diagnosis is proven and any residual
   difference is the estimator's real contribution. **This is the decisive experiment; run it first
   after Phase 0.**
2. 125's own frozen nulls (target-label and response-label derangement, 200 draws each) must still
   be beaten on whatever endpoint you promote.
3. 125's precision gain **failed its family-share gate** (Orion supplied 0.5667 of the positive gain
   against a 0.5 maximum) and **did not exceed either complete-procedure null p95** (p = 0.33, 0.13).
   Those two failures are about the *precision* endpoint and they still stand — do not promote on
   precision without clearing them, and say plainly if you cannot.

## Deliverable
`RESULTS.md`: the Phase-0 verification, the pure-rescale control, the re-scored endpoints in `pds` /
cosine / `q_w`, the direction-only variant, the six-member price, and a verdict in overall units.
**"The rescale control explains everything and the estimator adds nothing" is a full result** — it
would close the family properly instead of leaving it closed for the wrong reason.

## Territory
193 operator eigenvalue floor · 194 donor tables · 196 emitted depth · 200 shipped-stack ablation ·
202 the call-set build. **You own the panel-level hierarchical estimator.** Do not tune the operator
or the donor tables; if your result implies a change there, hand it over.
