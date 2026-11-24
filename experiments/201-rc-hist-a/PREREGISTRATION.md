# PREREGISTRATION — `rc_hist_a`, written 2026-09-03 11:55 UTC, BEFORE any build

Coordinator (lead). Committed before the builder is touched. Not editable after the board result;
corrections go in a separate `BOARD-RESULT.md`.

## 1. The arm

`rc_hist_a` — 185's **full-axis re-key**, against `rc_hybrid`'s **call-set-preserving** re-key.
Identical in every other respect: a15's configuration, a15's per-gene totals, per-cell totals welded
to `T = S // 400`, a15's own per-perturbation `|z|` multiset, `abs` sign policy. **One variable: the
re-key is applied to the whole gene axis instead of inside a15's BH call set.**

The `abs` policy is retained. 185 flagged that it leaves 33% of designed genes with a Wilcoxon `z`
sign opposite their lfc sign, and that a scorer reading the `z` sign would punish it. **The
`rc_hybrid` board point settles that: `cell-eval2` does not read it.**

## 2. Why this slot

`rc_hybrid` proved the delivery converts and refuted rule 13's toll for profile-preserving arms. It
left one thing unmeasured, and it is the thing experiment 192 needs: **what the board charges for
changing the call set.** 185 measured on h1-test that a full re-key buys `reach` and costs `jac`,
because *the reference's DE set is detectability-driven and a15's expression-dominated call set
matches it better than `|H|` does* — but the h1-test → board conversion for `jac` is only known
through 181's 0.44 ratio, which has never been tested on a designed-order arm, and `rc_hybrid`'s own
`jac` move was too small to calibrate it (+0.00013 raw).

**So this slot measures the `jac`↔`reach` exchange rate on the board, for the class 192 is about to
design in.** Every call-set arm after it is priced against this point.

## 3. Predictions

h1-test (185 §3b), `rc_hist_a` against `rc_hybrid`: `reach` 0.23278 vs 0.21681 (**+0.01597 raw**),
`jac` 0.05932 vs 0.06827 (**−0.00895 raw**), `fid` 0.50808 vs 0.50693 (+0.00115), `n_pred` 6040 vs
6043, `nmae` 0.96670 both.

Converted with the factors the `rc_hybrid` point now fixes — `reach` **×1.290** (measured on that
arm, superseding 164's 0.4277), `jac` ×0.44 (181, **untested on this class — that is the point**):

| member | predicted Δ raw board | Δ scaled | Δ overall |
|---|---:|---:|---:|
| `reach` | +0.0206 | +0.0233 | **+0.00388** |
| `jac` | −0.0039 | −0.0107 | **−0.00179** |
| `fid` | ≈ 0 | ≈ 0 | ≈ 0 |
| `pds`, `mse`, `nmae` | 0 by construction | 0 | 0 |
| | | | **+0.0021** |

**Predicted board overall ≈ 0.1855**, against `rc_hybrid`'s 0.1834190. Composed gate (+0.0019
shape-only, and rule 13's toll does **not** apply — profile-preserving, per the amendment): **clears,
narrowly.**

## 4. Falsifiers, pre-committed

1. **`pds` or `mse` differs from `rc_hybrid`'s board value at all.** Both are identities under this
   construction; `rc_hybrid` reproduced a15's to sixteen significant figures. Any movement means the
   weld or the total repair broke. Check the per-gene panel-total digest against `rc_hybrid`'s
   **before** submitting.
2. **raw `nmae` differs from `rc_hybrid`'s 0.94327 by more than 0.0005.**
3. **`jac` falls by more than 0.006 scaled**, i.e. more than 3× the prediction. That would mean the
   call-set channel is far more expensive than h1-test says and 192 must be re-scoped.
4. **`reach` gains less than +0.010 scaled.** That would mean the 1.290 factor does not generalise
   within the class and is specific to `rc_hybrid`, retiring it as a planning constant.

## 5. What this is NOT

Not an attempt at a large gain — the prediction is +0.0021, barely over the gate. It is a designed
measurement of an exchange rate (`docs/SLOT-PROGRAM.md` row 3), spent because experiment 192 is
about to make call-set decisions that this number prices and no local instrument can.

**If experiment 192 lands with a better-supported call-set arm before this is submitted, 192's arm
takes the slot and this one waits.** That is a deliberate option, recorded here so the choice is not
made after seeing a result.

## 6. Discipline
Joke model name, **empty description** (`AGENTS.md` §14).
