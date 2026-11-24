# a12 board result — the cap term delivered 40% of its model, and falsifier 2 fired

Entry `ZRptoejHEjkO3RrAtDw3`, submitted 2026-08-29 15:45 UTC.
**Overall `0.16798741`** against a11's `0.16611503` — **`+0.001872`**, rank 8 of 434.

| member | a11 | **a12** | Δ | predicted Δ |
| --- | ---: | ---: | ---: | ---: |
| `pds` | 0.739889 | 0.737885 | −0.002004 | ~0 |
| **`mse`** | 0.082041 | **0.091163** | **+0.009122** | **+0.0229** |
| `nmae` | 0.096982 | 0.093974 | −0.003008 | ~0 |
| `fid` | −0.031535 | −0.031473 | +0.000062 | −0.0016 |
| `reach` | 0.112220 | 0.116790 | +0.004569 | −0.002 |
| `jac` | −0.002907 | −0.000413 | +0.002494 | ~0, sign only |
| **overall** | **0.166115** | **0.167987** | **+0.001872** | **+0.0038** |

Raw `mse` **0.9099893 → 0.9013492**.

## Falsifier 2 fired

Pre-registered: *"raw `mse` outside `[0.86, 0.90]`"*. Measured **0.90135** — outside, narrowly.
The stated consequence was that the decomposition's point estimates are retired, and they are.

**The cap term delivered −0.0086 raw where −0.0218 was modelled — about 40%.** Falsifier 1
(`score_mse` below `+0.0900`) did *not* fire, at `+0.09116`, so the tiered spread policy is
doing something; it is simply doing far less than `mech_ship_policy.py` measured. The most
likely reading is that `claim` did not reach the `#247` cap on the real panel the way the
re-emission grid said it would — which is checkable on this artifact and should be, before any
further arm is built on a cap-saturation assumption.

## What was right

**Three predictions the decomposition got wrong in our favour**, all in the DE members it does
not model: `reach` gained `+0.0046` against a predicted small loss, `jac` gained `+0.0025`, and
`fid` was flat rather than `−0.0016`. The `n_pred` probe (4,440, mean across contexts) was the
basis for expecting no coverage cliff, and there was none.

**`nmae` cost `−0.0030`**, which the 096 model did not price for this arm at all — the same
member that experiment 097's gate later identified as binding on the shrinkage family. Two
arms in a row where `nmae` has been the unmodelled cost.

## Consequence for a13

a13's predicted raw `mse` of `0.8212` was computed on the same cap-saturation assumption. If the
shortfall carries across identically, a13 lands nearer **`0.834`** than `0.821`. The
harness-measured gate is the more conservative and more trustworthy estimate: **`+0.0215`
(h1-test) / `+0.0040` (orion)** net against a12, and rule 1 takes the smaller. So a13 is
expected between **`0.172` and `0.180`**, not the `0.1866` the board model projected.

It was submitted anyway, and the reasoning is worth stating: it is gated positive on **both**
surfaces on all six members, its downside case is still positive, and the alternative was to
leave a slot unspent on the last full day before the field's own mse work closes the gap — four
teams passed us on this member between yesterday and today.


---

## Correction, 2026-08-29: `claim` DID reach the cap, and the "most likely reading" above is wrong

Measured on the shipped `submissions/a12-tier8000-cap.vcc` itself: `claim` = **0.245917**,
**99.75 % of the `#247` cap**, with the cap binding on **720 of 900 rows**. Substituting the
artifact's own `A`, `B` and `claim` for the grid's assumed ones moves the prediction by
**+0.00037** — *away* from the board. The emission did what it was built to do.

The miss was the model's arithmetic: `min(claim, B)` was applied to the **context means** while
the scorer forms #348's budget per context, and at a12's raised dispersion context C's `ρ` fell
to **0.9418**. Worth **+0.0218 raw** against a total error of `+0.0131`. The pre-registered
consequence stands — those point estimates are retired — but the corrected model is not a patch:
it acquires no parameter, and the three board arms whose `claim` never touches the `C_p` proxy
now fit one `(X, ΣD)` to `≤ 0.0001` raw. Full working: **`A12-ARTIFACT-AUDIT.md`**.

Two further corrections to this file: the `nmae` sentence ("the second consecutive arm where
`nmae` is the unmodelled cost") — `nmae` moved **+0.0026 on a11 and −0.0030 on a12**, opposite
signs, and both are consistent with a single monotone response to `VCC_ALPHA` at
`−0.011 / −0.017` scaled per unit α. It is mis-modelled as invariant, not unmodelled. And the
`reach` gain of `+0.0046` is **five times below** the same-policy re-seeded twin floor the 097
gate measured (`0.023–0.028` scaled), so it should not be read as a favourable surprise at all.

a13's expected score is revised down: **`0.170–0.175` overall**, with a downside below a12, not
the `0.172–0.180` stated above.
