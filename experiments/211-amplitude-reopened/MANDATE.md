# 211 — `rc_callset` moved the `mse` cliff 2.6× further away. 187 closed amplitude against the OLD cliff.

You are a research subagent. Work ONLY in `experiments/211-amplitude-reopened/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 12 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §4.1, §4.1b, §8.1, §8.2, §11.1 rules 1–13 with the rule-13 amendment. Then:
1. `experiments/187-joint-amplitude/RESULTS.md` **in full** — the closure this mandate reopens.
   **Its instrument and its known-answer rows are yours.**
2. `experiments/208-free-riders/RESULTS.md` — **its Phase-4 by-product is the reason this exists.**
3. `experiments/202-callset-rider/BOARD-RESULT.md` and
   `experiments/185-rank-controlled-emission/BOARD-RESULT-rc_hybrid.md` — the two board points that
   changed the emission underneath 187's algebra.
4. `experiments/170-gamma-mse/RESULTS.md` (the ρ knee), `096-mse-floor/` (the `claim` route),
   `193-whitener-shrinkage/RESULTS.md` §Phase 3 (`nmae` wants *larger* `k`).

## The argument

187 closed global amplitude with a sharp result: the joint `mse` + `nmae` optimum for a global α is
**0.982 to four decimals**, because *any* increase crosses **context A's ρ cliff at α 0.9823**, where
`mse` loses **30×** what `nmae` gains per unit `ln α`. That was measured on a15's emission.

**Two board submissions have since replaced the emission underneath it.** 208 measured, as a
by-product of an unrelated screen:

> `rc_callset`'s **Ĉ/C_p is 2.940** against a15's **1.125** — the `mse` cliff that vetoed the depth
> knob is now **2.6× further away.**

`mse`'s numerator is `s²A + Kc − min(claim, s²B) − 2sX`. The cliff is where `s²B` overtakes `claim`
and the `min` saturates, after which extra amplitude stops being refunded. **`Ĉ/C_p` rising from
1.125 to 2.940 means our emitted per-cell noise is now far above the reference's, so `claim` sits
hard at its cap with much more room before `s²B` reaches it.** 187's "30× asymmetry at α 0.9823" is
a statement about a15's `B`, not `rc_callset`'s.

**And the other side of the trade now wants more amplitude, from two independent directions.** 193
§Phase 3 derived that with `q = 0.5375` and `k = 0.278`, `nmae = 1 + k(1 − 2q)` **improves with
larger `k`** — and it used that to *forbid* regularisation. 187 §5 measured the |t|-weighted per-gene
`k` at **median 0.12, p99 1.11**, i.e. we are far below the `k = 1` kink for almost every gene, on
the cheap side of an asymmetric loss (slope −(2q−1) below, +1 above).

**So: the member that blocked amplitude has moved its blocker 2.6× away, and the member that wanted
amplitude still wants it. The joint optimum has to be re-solved, and 187's own algebra does it.**

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
Re-run 187's instrument unchanged and reproduce **its** rows first:
- the C_p-free LS solve: `X` **0.019240**, `ΣD` **0.229141**, residuals ≤ 6.1e−5 raw;
- the four pre-registered proxy hypotheses to ≤5e−5, and the per-context factor
  `f = (1, 1, 0.9825)` retrodicting **abc-knee** at +0.00018 scaled;
- `claim(a15 as shipped)` **0.24581**; a15 raw `mse`/`score_mse` **0.88228 / 0.11112**;
- the two out-of-sample predictions: a10 **1.3567** (board 1.3524), a11 **0.9021** (board 0.9100).
- **Then reproduce 208's `Ĉ/C_p` figures independently**: a15 **1.125**, `rc_callset` **2.940**.
  If you cannot, this mandate's premise is unverified and you should say so and stop.

### Phase 1 — re-locate the cliff on the CURRENT emission
Compute, per context, `claim`, `B`, `s²B/claim` and the cliff position `s*` for **each** of a15,
`rc_hybrid` and `rc_callset`. Report where a15 sat (187 measured ρ = 0.992/0.940/0.882, mean 0.938)
and where `rc_callset` sits. **State the new per-context cliff α explicitly** against 187's
0.9823 / +3.1% / +5.8% for A/B/C.

Note 202 measured `rc_callset`'s challenge `n_pred` at 2,599 against a15's 4,474 and 185/191
measured the emitted per-cell law changed while the *pseudobulk* stayed bit-identical — so `A`, `X`
and `Kc` should be **unchanged** and only `B`/`claim` should move. **Verify that decomposition
rather than assuming it**; if `A` or `X` moved, 187's solve does not transfer and you must say so.

### Phase 2 — re-solve the joint optimum
Rebuild 187's `overall(α) = [mse_scaled(α) + nmae_scaled(α)]/6` on `rc_callset`'s `B` and `claim`,
and report the maximiser. Use 187's measured board `nmae` law on the current pipeline —
`raw_nmae(α) = 0.944330 − 0.0416·mean_c ln(α_c/0.982)` — and say where its linear range ends
(187 put the h1-test turnover at `s ≈ 1.29`, and 193 independently says `nmae` improves with larger
`k` throughout our range).

**Then extend the objective beyond two members**, which 187 explicitly did not do: it carried
`fid`/`jac`/`reach` at zero with 180's twin band. We now hold **two** board points whose emission
is `rc_callset`'s, so report what the DE members did across the a14→a15→…→`rc_callset` amplitude
history and whether carrying them at zero is still defensible.

### Phase 3 — retrodict, or the result is worthless
**Every board point on this axis must be retrodicted by the re-solved model**, not just the ones it
was built from: a10, a11, a12, a14, a15, a19 (sign only — operator ridge, out of domain),
**abc-knee**, and now `rc_hybrid` and `rc_callset` (both of which should read Δ`mse` = **exactly
0.000000**, since the board reproduced `mse` to sixteen significant figures — that is a free and very
sharp gate on your `B`/`claim` bookkeeping). If your model cannot reproduce `rc_hybrid`'s and
`rc_callset`'s exact zeros, your `claim` is wrong.

### Phase 4 — the gate, and the interaction you must not ignore
Gate: **+0.0089** for an amplitude-changing arm — and unlike 210's case this arm genuinely *is*
amplitude-changing, so do not argue for the shape-only floor. Price all six (rule 12).

**Critical interaction:** `cd4f034` (`experiments/210-cd4-build/`, built, submitting at the 00:00 UTC
reset) raises the CD4 dose 3.4× and 210 measured its delivered `A` at **0.97970× (−2.03%)** — i.e.
**it already spends some of the amplitude headroom you are measuring, in the favourable direction.**
Read its `PREREGISTRATION.md`, `PREREGISTRATION-ADDENDUM.md` and `RESULTS.md`, and report your
optimum **both** on `rc_callset` and on `cd4f034`'s delivered profile. If the two disagree, the α
that ships depends on whether `cd4f034` lands — say so, and give both.

## Deliverable
`RESULTS.md`: the known-answer rows including 208's Ĉ/C_p, the per-arm cliff table, the re-solved
joint optimum with the DE members no longer assumed dead, the full retrodiction table, and a verdict
in overall units for both base profiles. **"The cliff moved and the optimum did not" is a full
result** — it would close global amplitude for the second and final time, on the current emission
rather than on a15's.

## Territory
210 owns the CD4 dose. 208 owns emission riders (its Hamilton rider is a vetoed but live +0.0029
standalone — do not re-price it). **You own the global amplitude scalar and the ρ cliff.** Do not
touch the operator (193), donor tables (194, 204, 209), or the call set (192, 202).
