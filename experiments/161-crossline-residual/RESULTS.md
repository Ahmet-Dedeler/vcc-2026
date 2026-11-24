# 161 — the cross-line residual is not predictable from context, and the gain that looked like it was is a re-shaped consensus

**Date:** 2026-09-02 · five cell lines, two surfaces, delivered-space scoring, deterministic
(seed `20260902`) · **no submission was built, no `cell-eval2` was run, no file outside this
directory was modified.**

---

## Verdict

**NO-GO**, and the closure has a mechanism with three independent legs.

1. **Context features carry nothing.** A context-only model is annihilated by the mandatory
   across-target centring — not approximately, measurably: its eval `R²` on the residual is
   **−0.0024**, its cosine to the true residual is **−0.0033**, and its delivered Δpds is
   **statistically identical to its own permutation null** on all five lines.
2. **Context does not beat target-side, it slightly loses to it.** Cosine between the
   predicted and the true line-specific residual: **T (no recipient information at all)
   0.0606**, X (consensus×context) 0.0592, ALL 0.0550. Adding the recipient lowers it. The
   two scalars a donor table cannot express — per-gene control-cell **detection rate** and
   **dispersion**, built here from 419k control cells across five lines — make it **worse**:
   dropping them takes the worst-line Δpds from **−0.00097 to +0.00382**.
3. **The gain that survives is not residual prediction.** A closed-form, zero-learning,
   context-free per-entry re-shaping of the consensus, `C ← sign(C)·|C|²`, is worth
   **+0.0118 worst-line / +0.0206 mean** raw pds — **3.1× the entire 27-feature learned
   model** (+0.00382 / +0.00549). Put that re-shaping in the baseline and the learned model
   adds **+0.0024 worst-line on `core4` but −0.0027 on `core5`**, while costing **−0.007 to
   −0.013 aggregate cosine on every line of both surfaces**.

Inconsistent in sign across two surfaces, and it sells cosine to buy pds. Under the repo's
worst-surface selection rule that is a NO-GO. **Do not reopen cross-line residual prediction
without a mechanism that beats a per-entry magnitude re-shaping of the consensus**, because
that is what a flexible learner finds when pointed at this objective.

This is the fifth route into 159/160's wall and the first that attacked the residual rather
than the consensus, as `docs/CURRENT_STATE.md` asked for.

---

## 0. The bar, and the known-answer rows (AGENTS.md §11.1 rule 10)

**The bar is the equal-weight consensus of the donor LINES, scored in the delivered space**
(gate → production `mass_consistent_effect` → render → `156/screen.py:endpoints`).

`core5` = K562 / RPE1 / Jurkat / HCT116 / HEK293T, **1,883 shared targets × 4,980 genes, and
zero of the 1,883 are on the challenge panel** — the trainable pool is disjoint from the
deployment set by construction (`inventory.py`).

| held line | pds | agg cos | row-row | sse* |
|---|---:|---:|---:|---:|
| HCT116 | 0.83887 | 0.13479 | +0.00890 | 0.963687 |
| HEK293T | 0.83242 | 0.11338 | +0.00945 | 0.977312 |
| JURKAT | 0.88723 | 0.14419 | +0.01022 | 0.960592 |
| K562 | 0.88463 | 0.14040 | +0.00902 | 0.962246 |
| RPE1 | 0.88824 | 0.12901 | +0.00577 | 0.967893 |
| **MEAN** | **0.86628** | **0.13236** | | |

Every fold aborts unless these reproduce. All passed on all folds, both surfaces:

| known-answer row | expected | measured (core5, mean over 5 folds) |
|---|---|---|
| truth is on the mass manifold | `Σexpm1(b+E)=TS` | max\|mass−1\| ≤ **4.2e-9** |
| truth is renderer-invariant | `render(E)=E` | ≤ **1.2e-7** |
| delivered zero effect is zero | ‖·‖ ratio ~0 | **2e-15** |
| zero effect scores at chance | pds 0.5 | **0.50000** |
| ORACLE truth, ungated, α=1 | pds 1, cos 1 | **1.00000 / 1.00000** |
| ORACLE truth through `deliver()` | pds 1 | **1.00000**, cos 0.98314 |
| ORACLE residual `R_L` (this experiment's ceiling) | pds 1 | **0.99996**, cos 0.92118 |
| `par-CONTROL` (component ∥ to C) | Δ = 0 | **+0.00000 / +0.00000** |

Three of these caught real bugs during the run and are recorded in §6.

---

## 1. The architecture, and the two parameterisation traps it hid

`y = R_L − C` where `R_L = E_L − mean_t(E_L)` and `C` = equal-weight mean of mean-removed
donor lines. The model is a `HistGradientBoostingRegressor` over **27 per-(target, gene,
context) scalars with no gene identity, no target identity and no line identity** — the shape
`AgingLab_BI` shipped and lost with (0.1746 → 0.1673). Families are defined by *what they
require*:

- **T (9)** computable **without the recipient**: `cons, acons, dsd, sign, tmag, gresp, bD, bsd, x_cons_gresp`
- **C (10)** recipient context: `bB, bdiff, btgB, btgdiff, div` + cell-level `cdet, cxvar, cfano, cdetdiff, cxvardiff`
- **X (8)** consensus × recipient context: `sat, x_cons_bB, x_cons_bdiff, x_cons_btgB, x_cons_cdet, x_cons_cxvar, x_cons_cdetdiff, x_cons_cxvardiff`

**`yhat` is centred across targets before use** — a target-independent addition is the metric's
zero *and* a `pds` hub, which is 160's whole result. A gene-only feature therefore predicts a
value constant across targets and centring annihilates it **exactly**, so 159's context claim
can only appear through the **interaction**. That is why the X family exists.

**Trap 1 — `sat` was mislabelled.** `sat = C/(1+b_recipient)` sat in the T family in the first
draft. It is definitionally an interaction with the recipient's basal profile, and leaving it
there would have moved credit from context to target — i.e. corrupted this experiment's
headline. Fixed before any core5 run.

**Trap 2 — the arm's parameterisation, which also made the null vacuous.** The fitted residual
comes out at **cos(yhat, C) = −0.97**: the model overwhelmingly learns the *shrinkage*
`y ≈ −C`, because only ~30% of the consensus transfers (159). Every endpoint here
unit-normalises each row, so **the component of `yhat` parallel to `C` is exactly inert**. Left
in:

- the arm `C + w·yhat` becomes `(1−a·w)C + w·u` with `a ≈ 0.94–0.98`, which passes through zero
  just above `w = 1` and then **flips the sign of the consensus** — the −0.4 pds I first
  measured at `w = 2` was that, not a model failure;
- the permutation null then read **−0.05** for a trivial reason (it splices in another
  target's consensus with a negative sign), so it would have looked like a passed control
  **even if the model carried no target-specific information at all**.

So the arm uses only the component **orthogonal to `C`**, at restored row norms, and the
parallel component is scored as the `par-CONTROL` known-answer row above. After the fix the
permutation null reads **−0.0004 to −0.005** against arms of **+0.004 to +0.007** — a null that
can actually fail. This is the third vacuous-control near-miss in this project in three days.

---

## 2. Leave-one-LINE-out, leave-one-TARGET-out, and their interaction

Best worst-line arm in each regime, Δpds vs the consensus. R1 = line held, targets seen on
other lines. R2 = line seen, targets held. **R3 = both held — and R3 is the deployment
regime**, because the challenge's 300 panel targets appear in none of the 1,883 trainable ones.

| regime | core5 worst / mean | core4 worst / mean |
|---|---:|---:|
| **R1** line unseen, targets seen | `ALL-cell w=1` **+0.00382** / +0.00549 | `T w=1` **+0.00368** / +0.00404 |
| **R2** line seen, targets unseen | `T w=1` **+0.00256** / +0.00668 | not run |
| **R3** both unseen (**deployment**) | `ALL w=0.5` **−0.00096** / +0.00258 | `T w=1` **+0.00342** / +0.00406 |

**The interaction is the result.** On `core5` the model needs *either* the line or the targets
to have been seen; hold out both and the worst line goes negative. On `core4` it survives.
Since the architecture has no target identity, the R1→R3 drop on `core5` says the 27-dim
feature vector is acting as a target fingerprint. I did **not** separate that from the smaller
training-target count in R3 (2/3 of targets), so I state the regime gap as measured and the
memorisation reading as the most likely of two explanations, not as established.

Per-line detail for the best R1 arm, with the target-cluster bootstrap (400 draws, 40 clusters)
and the permutation null:

| `ALL-cell w=1.0`, core5 R1 | Δpds | bootstrap 95% | permutation null |
|---|---:|---|---:|
| HCT116 | +0.00647 | [+0.00299, +0.01093] | −0.00373 |
| HEK293T | +0.00705 | [+0.00483, +0.00969] | −0.00512 |
| JURKAT | +0.00458 | [+0.00331, +0.00595] | −0.00035 |
| K562 | +0.00552 | [+0.00374, +0.00700] | −0.00124 |
| RPE1 | +0.00382 | [+0.00191, +0.00604] | −0.00403 |

All five intervals exclude zero and all five nulls are ≤ 0. **The R1 effect is real. It is
just not what it looks like** (§4).

---

## 3. Which feature families carry the signal — CONTEXT LOSES

`cos(orth, TRUE orth)` is the direct measurement: cosine between the model's predicted
line-specific residual and the true one, per target, mean over the five held lines.

| feature set | eval R² on `y` | cos(yhat, C) | \|orth\|/\|C\| | **cos(orth, TRUE orth)** |
|---|---:|---:|---:|---:|
| **T** — no recipient information at all | 0.1782 | −0.976 | 0.137 | **+0.0606** |
| X — consensus × context | 0.1814 | −0.985 | 0.120 | +0.0592 |
| T+X | 0.1809 | −0.976 | 0.147 | +0.0579 |
| T+C | 0.1777 | −0.971 | 0.145 | +0.0554 |
| ALL (27) | 0.1774 | −0.971 | 0.149 | +0.0550 |
| ALL−div | 0.1791 | −0.969 | 0.149 | +0.0521 |
| ALL−cell | 0.1646 | −0.958 | 0.155 | +0.0605 |
| **C — context only** | **−0.0024** | +0.002 | 0.039 | **−0.0033** |

Three readings, and they answer the mandate's question directly:

- **The context-only model predicts nothing and is indistinguishable from its own null.** At
  `w=1` its Δpds is −0.00098 / +0.00002 / −0.00136 / −0.00186 / −0.00091 against nulls of
  −0.00026 / −0.00261 / −0.00174 / −0.00168 / −0.00097. Arm ≡ null, line by line. This is
  exactly what the centring argument predicts and it is the cleanest confirmation in the run.
- **Context features do not beat target-side features.** The best residual cosine belongs to
  the family that cannot see the recipient at all. Adding context moves it **down**
  (0.0606 → 0.0550).
- **The per-gene control-cell detection rate and dispersion — 159's main untested claim, and
  the reason this experiment built `ctrlstats.py` — do not help.** `ALL−cell` beats `ALL` on
  the worst line, **+0.00382 vs −0.00097**. They add noise dimensions, not information.

Cell-level stats were built from the raw sources for all five lines and each was validated by
rebuilding the control pseudobulk and comparing to the shipped `b_ctrl.npy`: RPE1 r=1.000000,
Jurkat 0.999995, HCT116 0.999939, HEK293T 0.999932, K562 0.995280 (K562's controls come from
the *essential* library because gwps is a 1.99M×8,248 dense array chunked (3886,33), so its
scattered control rows cost a 65 GB read; same line, same lab, same chemistry).

---

## 4. The mechanism: it is a re-shaped consensus, not a residual (`recal.py`)

The winning family, `T`, **contains no recipient information**. So whatever it does, it does
identically to every recipient — it cannot be carrying line-specific information. That makes
the honest reading "the equal-weight consensus is the wrong per-entry *shape*", which is a
consensus improvement, not a residual model.

Tested with three closed-form re-calibrations using only `C` and the cross-donor spread, no
learning and no context. Worst-line / mean Δpds over held lines:

| arm | core5 worst | core5 mean | core4 worst | core4 mean |
|---|---:|---:|---:|---:|
| **`power p=2`  `C ← sign(C)·\|C\|²`** | **+0.01178** | **+0.02056** | **+0.01231** | **+0.02174** |
| `power p=1.5` | +0.00809 | +0.01228 | +0.00806 | +0.01242 |
| `soft λ=0.5` (threshold on donor disagreement) | +0.00525 | +0.01320 | +0.00122 | +0.00485 |
| `hard keep-top 25%` | +0.00013 | +0.00506 | +0.00252 | +0.00611 |
| `eb κ=0.25` (per-entry empirical Bayes) | −0.00765 | −0.00497 | −0.00658 | −0.00582 |
| `power p=0.5` | −0.01901 | −0.01441 | −0.01628 | −0.01369 |
| **the entire learned model (`ALL-cell w=1`)** | +0.00382 | +0.00549 | +0.00234 | +0.00374 |

**A one-parameter context-free squaring beats the 27-feature learned residual model 3.1× on
the worst line of `core5` and 3.3× on `core4`**, and replicates across both surfaces
(per-line core5: +0.02404 / +0.03153 / +0.01178 / +0.01428 / +0.02116, every bootstrap
interval excluding zero). EB shrinkage is *negative*, consistent with AGENTS.md §12's closed
list — the GBM did not rediscover EB, it rediscovered magnitude sharpening.

**The decisive test: is the learned model additive on top of it?** Re-run with
`C_base = sign(C)·|C|²` (features still built from the raw `C`, which is what deployment has):

| | baseline mean pds | best learned arm, R1 worst | R3 worst (deployment) | Δ agg cosine |
|---|---:|---:|---:|---:|
| core5 | 0.86628 → **0.88684** | `ALL-cell w=0.5` **−0.00017** | `T w=0.5` **−0.00268** | −0.008 to −0.013 |
| core4 | 0.85655 → **0.87829** | `T w=0.5` **+0.00237** | `T w=0.5` **+0.00237** | −0.007 to −0.010 |

On `core5` the learned model adds nothing on the worst line in either regime. On `core4` it
adds +0.0024. **It disagrees in sign between two surfaces and it lowers aggregate cosine on
every line of both.** That fails the worst-surface rule, and a cosine loss is not free —
`mse` and `nmae` both read the profile directly (rule 12, and a19's −0.01352).

Hubness moves the right way for the sharpening (`row-row` 0.00890 → 0.00784 on HCT116), so
`power p=2` is acting as a crude hubness suppressor. My first reading of that was that it is
therefore **a poor imitation of the whitening production already runs** (AGENTS.md §4.2) —
**that reading is wrong, and prior work already settled it.** Experiment 030 measured the two
as roughly *additive* (whitening alone +0.080, sharpening alone +0.072, both +0.120), and
experiments **168/169**, running concurrently on the production-shaped delivered surface,
reproduce a global sharpening gain of **+0.01145 worst-surface Δpds at q = 1.15 on the already
whitened artifact**. So sharpening is not redundant with whitening; it is a separate lever.

That does not change this experiment's verdict — §4's point stands either way, and it stands
*more* strongly: the thing the learned residual model was partially rediscovering is a real,
independently-confirmed, context-free consensus lever that another experiment already owns.

---

## 5. All six members (AGENTS.md §11.1 rule 12)

Nothing is being shipped, so nothing moves; this prices the arms that were *considered*.
DE-side statistics are per target then averaged, never pooled (rule 8), over each target's top
200 \|truth\| genes, on the delivered effect.

| member | how read here | consensus | `power p=2` | ORACLE |
|---|---|---:|---:|---:|
| `pds` | **measured**, delivered space | 0.86628 | **+0.02056** | 1.00000 |
| `mse` | `sse*` proxy = 1−cos²_agg at best global amplitude (156's proxy, flagged) | 0.966346 | 0.968847 (**worse**) | 0.000 |
| `nmae` | proxy: mean\|pred−truth\| / mean\|truth\| over the pool | 0.9225 | **0.9899 (worse)** | 0.0370 |
| `reach` | proxy: deepest ≥90%-sign-pure prefix of our \|effect\| ordering | 0.1674 | 0.1677 (flat) | 1.0000 |
| `fid` | **NOT COMPUTED** | — | — | — |
| `jac` | **NOT COMPUTED** | — | — | — |

Sign accuracy over the reference pool: 0.5942 → 0.6002; head purity (top decile) 0.7367 →
0.7361. The ORACLE row reads head purity 1.0000, reach 1.0000, nmae 0.0370 — the known-answer
row for this table.

**Why `fid` and `jac` are unmeasured, stated rather than laundered.** This surface is donor-table
pseudobulk: there are no cells, so there is no Wilcoxon and no call set, and `fid`'s cohort *is*
our call set. Independently, both carry `n_real` in a denominator or a union, which rule 9 says
is unreadable on every harness we hold. **I cannot claim they would not move** — an effect-direction
change does move them. This is a genuine gap; it does not affect the verdict only because the
verdict ships nothing. Note `power p=2` visibly costs `nmae`, which is exactly the shape of the
a19 failure and a second reason not to treat §4's sharpening as a free win.

---

## 6. Controls that could have failed, and the three bugs they caught

Two vacuous checks were caught in this project in the two days before this run, so every
control here was designed by asking what it would read **if the thing being tested were absent**.

1. **`par-CONTROL` reads +0.00000/+0.00000** — the component of the prediction parallel to `C`
   is provably inert under row-normalised metrics. It is *not* vacuous: it reads a large
   negative when the parallel coefficient drops below −1 and inverts a row, which is how the
   `w`-sweep pathology in §1 was found.
2. **The permutation null went from −0.05 (trivially bad) to ≈0** once the parallel component
   was removed. Before the fix it could not have distinguished a real target-specific
   prediction from none at all.
3. **`ctrlstats.py`'s pseudobulk round-trip caught a silent gene-axis failure**: Jurkat stores
   symbols as a legacy `__categories` categorical, so 0 of 8,881 genes mapped and every
   statistic would have been NaN-medianed into a constant. The check reads r=nan and aborted.
4. **The `zero effect` row caught my own over-tight tolerance twice** — first on the
   render round-trip (3.7e-9 against a 1e-9 gate), then on the *direction* of a 2e-15 vector,
   which is meaningless. The check now asserts `‖deliver(0)‖/‖truth‖ < 1e-9`, which is the
   quantity that can actually fail, and the pds-at-chance row separately.

**Limits of the controls.** R2/R3 rows carry **no permutation null or bootstrap** in this build
(both are computed only on the ungrouped R1 rows). The R3 verdict therefore rests on the
worst-line point estimate across two surfaces, not on an interval.

---

## 7. What I did NOT measure, and why

- **`fid` and `jac`** — not computed at all (§5). No call set exists on this surface and both
  carry `n_real` (rule 9). They are *not* covered by the "cannot move" exemption.
- **Whitening.** The entire experiment runs on an **unwhitened** equal-weight consensus, while
  production whitens (§4.2). **Every delta here is measured against a weaker baseline than
  production's**, so no absolute here — including §4's +0.0206 — is a board number. What this
  does *not* mean, having checked: sharpening is not merely re-deriving whitening. 030 measured
  them as roughly additive and 168 reproduces the gain on the whitened artifact (§4). The
  residual verdict is unaffected either way, because the residual model loses to sharpening on
  *this* surface, where both are measured against the same baseline.
- **The challenge panel.** `core5`'s 1,883 targets contain **zero** panel targets by
  construction. 159 §5.3 measured panel cosines 3–5× smaller than non-panel ones, so absolute
  levels here are not board levels. Only the ordering of arms is claimed to transfer.
- **Cell depth.** Donor tables carry a median 77–203 cells/target against the challenge's 400.
  159 §1.3 puts the *residual*'s split-half reliability at 0.05–0.12 in that regime, so the
  target being fitted is mostly noise: disattenuating `cos(orth, TRUE orth) = 0.06` by
  `√ρ ≈ 0.32` gives a true residual predictability near **0.19**, which is not zero.
  **This experiment cannot rule out a residual signal that 400-cell references would reveal**,
  and that is the honest limit of the NO-GO.
- **CD4 and KOLF2** as donors — 100% panel, zero non-panel shared targets, so they cannot
  train anything (`inventory.json`). They remain inference-time-only donors.
- **`core6`** (two K562 libraries as one line) and the `wideh1` surface beyond the smoke test.
- **H1 control stats** — `ctrlstats.py H1` was still running at write-up; `core5`/`core4` do
  not use it.

---

## 8. What would change the verdict

A residual model that beats **`power p=2` as its baseline**, on the **worst** line of **both**
surfaces, in **R3**, **without losing aggregate cosine**. Nothing here does. The nearest miss
is `core4 R3 T w=0.5` at +0.0024 worst-line, which loses 0.007–0.010 cosine and reverses sign
on `core5`.

**The sharpening lead is already owned — do not re-run it here.** §4's transform is the same
object as experiments **030 / 168 / 169**, which are pricing it properly on the
production-shaped delivered surface (168 reproduces +0.01145 worst-surface Δpds at q = 1.15
post-whitening; 169 is completing the six-member pricing, and 168's own banner says do not ship
on it until then). This experiment's contribution to that thread is one **independent
corroboration from a different direction**: a flexible 27-feature learner, pointed at the
line-specific residual on five cell lines, converges on magnitude sharpening as the dominant
usable structure and finds nothing else. Sharpening's `nmae` cost, which is what parked
exp 030, reproduces here too (§5: 0.9225 → 0.9899).

**Integration point for this experiment: none.** Nothing here is proposed for the builder.

---

## 9. Files

| file | what it does | ran? |
|---|---|---|
| `common.py` | donor-table loader, line-GROUP map, blob flags | — |
| `inventory.py` | §1 census, panel/non-panel split, leakage audit | ✅ `work/inventory.json` |
| `build.py` | mass-consistent multi-line surfaces, restriction done in profile space | ✅ `work/surface-{core5,core4,wideh1}.npz` |
| `ctrlstats.py` | per-gene control-cell detection rate + dispersion from cell-level h5ad/parquet, each validated against the shipped `b_ctrl` | ✅ `work/ctrl-{K562,RPE1,JURKAT,HCT116,HEK293T}.npz` |
| `residual.py` | the model, families, regimes R1/R2/R3, nulls, bootstrap, `--base-power` | ✅ 6 runs |
| `recal.py` | §4 — closed-form context-free consensus re-calibrations | ✅ `work/recal-core{5,4}.json` |
| `de_price.py` | §5 — DE-side ordering/sign members | ✅ `work/de_price.json` |
| `report.py` | aggregation into the tables above | ✅ |

```bash
cd experiments/161-crossline-residual
uv run python ctrlstats.py                                    # ~6 min, peak RSS ~1 GB
uv run python residual.py --surface wideh1 --regimes R1 --no-ctrl \
    --n-train 300000 --max-iter 80 --out work/smoke-wideh1.json   # smoke first, 30 s
uv run python residual.py --surface core5 --regimes R1 \
    --fsets T,C,X,T+C,T+X,ALL,ALL-div,ALL-cell --out work/residual-core5-R1.json
uv run python residual.py --surface core5 --regimes R3,R2 --fsets T,ALL \
    --out work/residual-core5-R3R2.json
uv run python residual.py --surface core5 --regimes R3,R1 --fsets T,ALL-cell,ALL \
    --base-power 2.0 --out work/residual-core5-pow2.json
uv run python recal.py --surface core5 && uv run python de_price.py --surface core5
uv run python report.py work/residual-core5-*.json
```

Peak RSS 9.7 GB (the R3 folds), ~50 min of compute for everything above.
