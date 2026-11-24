# PREREGISTRATION — `h1union`, written 2026-09-04 01:40 UTC, BEFORE any build

Coordinator (lead). Committed before the builder is touched. Not editable after the board result.

## 1. The arm

Add the **union of two H1 sources — Pan-2026 H1 (64/300) and `h1_2025` (20/300), 7 shared, union
77/300** — as an un-whitened post-whitening residual **sized on its carrying rows**, on top of
`cd4f034`. One variable against the incumbent: **a coverage-aware sizing rule.**

Incumbent = `cd4f034`, entry `QSe3JlNUYZ2BzvlcsL14`, board **0.2049195830**, rank 7.

## 2. Why — a silent bug, not a new idea

218's headline is not that it found a new source. It is that **production has been silently
discarding two sources we already hold.**

`solve_cd4_fraction` sizes a residual by the **median row norm over every panel row**. A source
covering ≤ 50% of the panel has a median row norm of **exactly zero**, so its solved fraction is
**exactly zero** and it contributes **nothing** — with no error, no warning, and every downstream
check passing. **Pan-2026 H1 and `h1_2025` are therefore identically absent from every build this
project has ever made.** This is the same failure class 200 and 207 found in CD4's own sizing on
h1-test, where 4/85 coverage silently zeroed the residual.

It also explains a prior closure: **117/118 closed Pan on a `pds`-only price *and* through the very
sizing rule that set its dose to zero** — i.e. it was measured at a dose of nothing.

## 3. What 218 measured, sized on carrying rows instead

| endpoint | value |
|---|---|
| macro cosine, three cancer-line folds | **+7.68%** |
| **`pds`** | **+0.0049 — UP** |
| fold × dose cells positive on cosine *and* `pds` | **12 of 12** |
| six-member price (§8.4 rule, no member carried at zero) | **+0.00919** conservative |
| marginal over the free row-matched-CD4 alternative | **+0.01054** |
| × rule-13 shape-only gate | **5.5×** |
| leave-one-fold-out, dose chosen on the other three | **+0.02830**, 4/4 folds |
| carrying rows vs 159's `1/k` prediction for k 3→4 (+13.7%) | +8.8% to +11.7% |

**`pds` goes UP, which no CD4 dose arm has ever done.** `cd4f034` paid −0.001345 scaled on `pds` to
buy its cosine; this gains both. That is the signature of a genuinely 4th independent estimate *on
the rows it covers*, and it is why this is not a repeat of the dose lever — which 216 has just
proven is exhausted at `f = 0.34` (the argmax over a 31-point grid and all three forms).

**Controls 218 ran, all clean:** shuffled control **strictly negative at every dose**; non-carrying
rows move **|Δ| ≤ 0.018 pp** (the arm does only what it claims); the row-matched CD4 control
reproduces **1/3.9** of the gain *and costs* `pds`, so this is not a row-reweighting artefact.
The `kolf2` fold is 3–5× the others because both sources and that recipient are pluripotent — it is
**excluded from every headline number above**.

## 4. Predictions

**Predicted board Δ vs `cd4f034`: +0.0055 to +0.0105**, i.e. board **0.2104 – 0.2154**.

I am deliberately predicting **below** 218's +0.00919/+0.01054 rather than at it, for two reasons
this project has now measured: 216 showed a proxy-vs-exact `nmae` difference flipped `f = 0.40` from
+0.00025 to −0.00029, and my own `cd4f080` extrapolation failed because I scaled a member past its
turning point. **The conservative end assumes the union's gain is ~60% realised.**

Per-member, using §8.4's corrected propagation (2.51×, range 2.00–2.56) and 219's warning that the
factor is **local**: `mse` up, `nmae` up **or flat** (this arm is not a dose increase, so the `k > 1`
overshoot that killed `f = 0.80` should not apply — **that is a prediction, and falsifier 4 tests
it**), `pds` **up**, `fid`/`reach` up with sign accuracy, `jac` ≈ flat.

## 5. Falsifiers — pre-committed

1. **`pds` does not rise.** 218's central claim, 12/12 cells. If `pds` falls, the arm is behaving
   like a dose increase rather than an independent estimate and the whole 4th-source reading is wrong.
2. **raw `mse` does not improve** on `cd4f034`'s 0.83635.
3. **Overall Δ outside [−0.002, +0.016].**
4. **raw `nmae` rises by more than 0.002.** This is the member that killed `f = 0.80` (−0.004298 at
   dose). If it rises here, the union is acting as an amplitude increase and §3's mechanism is wrong.
5. **Config diff is not exactly the pre-declared sizing-rule change plus its deterministic
   consequences**, with `operator.sha256` unchanged, `alpha` **0.982**, `cd4_target_f` **0.34**, and
   `effects.cd4.*` **untouched** — CD4's own dose and form must not move (216/217's territory).

## 6. What the slot buys even if flat
The first board test of a **coverage-aware sizing rule**, which is a builder-level fix that also
governs D/E/F: `docs/DEF-RUNBOOK.md` §0.2 already records that the CD4 banks are panel-restricted and
that coverage collapses on a new panel. **If the rule is wrong here it is wrong in October**, when
the board is dark and nothing reads back.

## 7. Discipline
Joke model name, **empty description** (§14). This is the second of the two 2026-09-04 UTC slots;
the first was `cd4f034`, and `cd4f080` was cancelled unspent.
