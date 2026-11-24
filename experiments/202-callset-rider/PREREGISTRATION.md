# PREREGISTRATION — `rc_callset` (192's rider on `rc_hybrid`), written 2026-09-03 12:10 UTC, BEFORE any build

Coordinator (lead). Committed before the builder is touched. Not editable after the board result.

## 1. The arm — one variable against `rc_hybrid`

```
call set per (context, perturbation) = top 400 by |H|  ∪  top 2,200 by control CPM   (n_pred 2,600)
order                                = 191's rc_hybrid law, unchanged
everything else                      = rc_hybrid, unchanged
```

Against `rc_hybrid`'s challenge `n_pred` of **4,480** and a15's **4,474**.

## 2. Why this and not `rc_hist_a`

`rc_hist_a` was preregistered for this slot at +0.0021 with an explicit option clause: 192's arm
takes the slot if it lands first with better support. It did. See
`experiments/201-rc-hist-a/SUPERSEDED.md`.

## 3. The mechanism, and why the two halves of the union are both necessary

192 measured that **the set maximising `jac` and the order maximising `reach` want different
things, and the conflict is repairable.** The reference's significant set is **detectability**-driven
— control CPM alone enriches for it at fold **2.05**, while our own `|lfc_pred|` reads fold **1.00**,
exactly random, and is *anti*-informative in its head. So a call set keyed to `|H|` matches the
reference's set poorly (`jac` falls), and one keyed to detectability throws away the `|H|` head that
`reach` walks (−0.008 raw `reach` at fixed size).

The union does both: **top 2,200 by control CPM** buys membership overlap, and **protecting the top
400 by `|H|`** inside the set preserves the order `reach` reads and adds **+0.011 raw** on top.

192 also falsified the two things that would have made this arm unnecessary or bigger:
`fid` is **closed** to the call set (over 7 keys × 10 sizes × 4 `n_real` levels, including the
reference's own oracle key, it never exceeds a15's by more than **+0.0041 raw** — it is a
sign-accuracy member), and **`n_real_p` is not predictable** at build time (LOTO Pearson +0.24, 6%
of variance; oracle per-perturbation sizing +0.036, buildable −0.097). So this is a **global**
call-set arm by measurement, not by convenience.

## 4. Predictions

192's h1-test measurement, marginal to `rc_hybrid`, at the retrodicted level: `fid` **+0.0010**,
`jac` **+0.0129**, `reach` **+0.0112** raw.

**Predicted board Δ overall vs `rc_hybrid`: +0.0091 to +0.0111**, i.e. board overall
**0.1925 – 0.1945** against `rc_hybrid`'s **0.1834190095**. 192 reports this is positive at every
`n_real` level it tested and at **both** `reach` transfer factors (0.4277 and the 1.290 the
`rc_hybrid` point now fixes), with and without rule 13's toll — which the `rc_hybrid` board point
refuted for profile-preserving arms anyway.

`pds`, `mse` and `nmae` are **pinned by construction**: the delivered per-gene panel totals are
unchanged, and `rc_hybrid` proved on the board that this construction reproduces both profile
members to **sixteen significant figures**.

**~88% of the predicted gain is `jac`**, the member we rank 169th of 612 on and the one with no
prior board point in its favour. That is the risk and it is stated up front.

## 5. Falsifiers — 192's four, pre-committed, plus two of mine

1. **Per-gene panel-total digests identical to `rc_hybrid`'s on 900/900 blocks.** Else `pds`/`nmae`
   are not pinned and the arm is void. Check **before** submitting.
2. **Challenge `n_pred` within ±3% of 2,600.**
3. **Board `pds` and `mse` reproduce `rc_hybrid`'s to ≥12 significant figures.**
4. **Board `jac` RISES.** It is 88% of the predicted gain. If `jac` falls, the detectability half of
   the union does not transfer and the whole call-set channel is closed on the board.
5. *(mine)* **Board `reach` does not fall by more than 0.005 scaled.** The union changes which genes
   are called; `reach`'s pool is the reference's and 171 proved the call set is a formal no-op for
   it, so a real `reach` loss would mean the order was disturbed by the set change.
6. *(mine)* **Rule 10 config diff**: exactly the pre-declared call-set difference against
   `rc_hybrid`'s written meta; operator sha256 identical; `after_all_median_norm` not drifted.

## 6. What this slot buys even if the arm is flat

The first board measurement of the **call-set channel** with the profile provably held fixed. `fid`,
`jac` and `reach`'s normalisation are unreadable on every local harness (rule 9), and the two
existing board points in this family (093/v14, a21) both confounded the set with the profile or the
emission law. This one does not.

## 7. Discipline
Joke model name, **empty description** (`AGENTS.md` §14).
