# a13 — pre-registered prediction

Written 2026-08-28 18:10 UTC, while the build was running and before any scoring. Both of
today's slots are spent; this is for the next available one.

**Arm:** a12's emission exactly, plus the two-part donor-agreement shrinkage and its amplitude
compensation. `VCC_ALPHA` at its module default 1.0796 — **not** a12's 1.0108, because the
amplitude is carried post-gate instead.

```
--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0 --overdispersion 4.0
--shrink-key agreement --shrink-toward zero --shrink-gamma 0.25
--shrink-pert-key donor-cos --shrink-pert-gamma 1.0 --shrink-pert-floor 0.3
--shrink-pert-missing damp --effect-scale 2.4004
```

## The intervention is verified, not assumed

Rule 3 says a result whose mechanism was never checked is not a result. The builder printed,
at build time:

| quantity | solved | built |
| --- | ---: | ---: |
| mean `z^2` (per-perturbation) | 0.351 | **0.3513** |
| graded gene fraction, median | 0.106 / 0.109 | **0.107** |
| targets shrunk | 294/300 | **294/300** |
| targets left at `alpha = 1` | 6 | **6** (all-zero key) |
| `z` floor | 0.30 | **0.3000** at min and q25 |

The estimator that was solved is the estimator that was built.

## Predictions, against a11 (`0.16611503`)

| member | predicted Δ | reasoning |
| --- | ---: | --- |
| `mse` raw | 0.9100 → **0.8212** | exact identity, `claim` at the cap, re-solved on the shipped functions |
| `score_mse` | **+0.0934** | |
| `pds` | **+0.0266** scaled | h1-test at three donors — **single surface** |
| `reach` | +0.010 | |
| `nmae`, `fid`, `jac` | −0.007 combined | |
| **overall** | **≈ 0.1866** (`+0.0205`) | |

## Falsifiers

1. **raw `mse` outside `[0.79, 0.86]`.** The decomposition has predicted two board results to
   within 1%; a third miss outside this band retires its point estimates.
2. **`score_pds` falls at all.** The most consequential and the least supported — it rests on
   a single surface. **If it fires, re-run at `--shrink-gamma 0.125` rather than abandoning
   the family**; the `mse` gain survives a smaller gamma and only the `pds` bonus is at risk.
3. **`n_pred` below 4,000.**
4. **`reach` down more than 0.02 scaled.**

## What is NOT established, stated plainly

- **The `pds` gain is one surface.** `orion` cannot test it: `orion_hct116` is itself a
  production donor and must be held out, so that harness is structurally two-donor, and forcing
  h1-test down to two donors reproduces orion *and flips the sign of the pds result*. h1-test
  at three donors is the only surface that matches production's donor set.
- **The `claim`-stays-at-the-cap check is per-row and has not landed.** The panel total is
  uninformative while `claim` is saturated — flat is equally consistent with "nothing moved"
  and "our margin is being eaten and has not run out yet". The diagnostic that decides it is
  the count of perturbations at the cap before and after, plus the margin distribution
  `C_hat_p / (K * C_hat_p^ref)` on the rows most shrunk. **If margin is being consumed, the
  failure mode is over-amplification — which costs `reach` and `nmae`, not `mse`, and would
  therefore not appear anywhere in the mse accounting that produced this prediction.**
- **`effect_scale` 2.4004 is the largest amplitude compensation this pipeline has shipped.**
  Clipping was measured at 0.08-0.25% and `n_pred` rises rather than falls, so both known
  checks pass; that is not the same as it being safe.

## Downside

If the `pds` gain is real, `+0.0205`. If `pds` is flat, `+0.0162`. If `pds` moves by the
two-donor amount instead, `+0.0148`. All three are positive; `pds` would have to fall by 0.12
scaled — about 4.5x the worst figure any surface produced — for the arm to lose. The `mse`
term dominates and is measured on the challenge panel itself rather than on a harness.

---

## The DE gate, run before the slot — and it reversed the agents' own recommendation

`experiments/097-direction/de_gate.py`, a13's exact configuration, both harnesses, production
`emit_counts` through the real rank-sum, 40 perturbations, **every arm with a re-seeded twin**.

| harness | arm | `n_pred` | `reach` | *(twin)* | `fid` | `jac` | `nmae` |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| h1-test | a12 | 3835 | 0.1434 | *0.1227* | 0.5016 | 0.0725 | 0.9630 |
| h1-test | **a13** | 3744 | 0.1442 | *0.1228* | **0.5041** | 0.0745 | **0.9678** |
| orion | a12 | 3500 | 0.0552 | *0.0797* | 0.5706 | 0.0145 | 0.9189 |
| orion | **a13** | 3516 | 0.0743 | *0.0780* | **0.5742** | 0.0143 | **0.9369** |

**Net vs a12: `+0.0095` (h1-test) and `+0.0025` (orion) — positive on both.**

### Four things this changed

**1. `reach` is unreadable in this run and no reach claim should be made from it.** The
same-policy twins differ by `0.0207` / `0.0245` raw — an order of magnitude larger than any
between-arm difference. That retires the *favourable* reach readings as well as the
unfavourable ones.

**2. `nmae` is the binding member, and nobody was pricing it** — not the direction study, not
`MSE-MECHANISM.md`, not the builder. The mechanism is clean: `nmae` is an absolute fold-change
error scored over **the reference's** DE genes, and those are not the genes the donors agree
about, so damping by donor agreement degrades it monotonically (h1 0.9630 → 0.9678 → 0.9768;
orion 0.9189 → 0.9369 → 0.9585). Every arm more aggressive than a13 goes **negative on orion**,
with `nmae` doing it each time.

**3. The more aggressive gammas both agents recommended lose.** `gamma_pert = 2` reads
`−0.0037` on orion; the direction study's own recommendation reads `+0.0174` on h1-test and
`−0.0020` on orion. Rule 1 binds on the worse surface. **a13 as built is the arm.**

**4. The `pds` threshold is not confirmed, but the concern is** — by a stronger route. The
specific "turns negative past `gamma_gene` ~0.7" does not reproduce at three donors. What does
reproduce is that **on orion `pds` falls monotonically from `gamma_gene = 0`**, so no positive
value is free and the conservative setting is the right one.

### Challenge-side `n_pred`, probed

4544 / 4406 / 3972 across contexts, **mean 4307**, against a12's 4440. Falsifier 3 (`n_pred`
below 4,000) does **not** fire on the mean, though context C at 3,972 sits just under it — worth
watching if `fid` disappoints.

### The finding that outlives this arm

Every default in the shrinkage implementation — the shrink target, both gammas, the floor, and
the missing-donor policy — **was chosen against `mse` alone**. `gamma_pert` is simply the one
the gate happened to test, and it was wrong by enough to flip the sign on a harness. The other
three rest on exactly the footing it rested on an hour earlier. They are being left as they are,
because substituting one unvalidated choice for another is not an improvement — but the help
text now says so, per member, rather than reading as though someone had checked.

---

## Re-gated on the estimator that was actually built — and it is better than the one gated first

The first DE gate ran a **power law** `alpha = (rank/G)^0.25`; a13 is the builder's
**exponential** `alpha = exp(-0.25 * rank / 2000)`. Two different functions, one name
`gg = 0.25`, across two codebases. Re-run on `agreement_shrink_vector`'s own form, Δ overall
against a12:

| arm | surface | `mse` | `pds` | `fid` | **`nmae`** | `jac` | **net** | `n_pred` |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **a13 as built** | h1-test | +0.0150 | +0.0044 | −0.0011 | **+0.0026** | +0.0006 | **+0.0215** | −7.7% |
| **a13 as built** | orion | +0.0214 | −0.0086 | +0.0057 | **−0.0145** | +0.0000 | **+0.0040** | −2.8% |
| the old `gamma = 1.0` default | h1-test | +0.0198 | −0.0017 | −0.0031 | −0.0034 | +0.0003 | +0.0120 | −12.5% |
| the old `gamma = 1.0` default | orion | +0.0313 | −0.0133 | +0.0097 | **−0.0615** | +0.0002 | **−0.0336** | −8.3% |

**a13 as built is net positive on both surfaces and better than the arm that was standing in
for it.** The mislabel was conservative — that was luck, not design.

It also retrospectively justifies dropping the old default: at `gamma = 1.0` raw `nmae` reads
**1.0667**, *above* that member's own baseline anchor of ~1.0013 — worse on `nmae` than
predicting nothing at all, which the mse side could not see.

### The caveat that matters if it disappoints

**`nmae` disagrees in sign between the two surfaces** (`+0.0026` on h1-test, `−0.0145` on
orion). Its contribution is not reliably signed, so a13's **worst-member floor is `−0.0092`**,
which is *worse* than the milder arm's `−0.0012` despite the better net. **a13 is the better
bet and not the safer one. If it underperforms on the board, look at `nmae` first.**

### The check that defuses the scale scare

The power law at `0.25` (mean `alpha` 0.834) nets `+0.0095 / +0.0025`; the exponential at
`0.05` (mean `alpha` 0.800) nets `+0.0102 / +0.0025` — within `0.0007` on both surfaces despite
different functional forms. **The amount of shrinkage transfers; the knob that produces it does
not.** Mean `alpha` is the quantity both sides should have been quoting all along.


---

## Correction, then RE-corrected: the gated configuration and the built artifact are identical

It was stated — by me, upstream — that the gate and the shipping configuration are "the same
object". That is right for 294 of 300 targets and wrong for the remaining six.

The gate computes the per-perturbation weight as
`max(where(n_donors >= 2, rank01(donor_cos), 0.0), FLOOR) ** gamma`, which replaces an
uncarried target's rank01 **value** with 0.0 and then floors it, landing at exactly **0.300**.
The builder's `--shrink-pert-missing damp` sets that target's **score** to 0.0 and ranks it with
the panel, so it lands wherever a zero cosine falls among the other 299 — **0.430** on the real
panel, because challenge cosines run `[-0.285, +0.439]` with median `+0.0284` and a zero
therefore sits near the 45th percentile, not at the bottom.

Carried targets are bit-identical (max absolute difference `0.00e+00`). Only the six uncarried
ones move, and they move by 0.13 in `z`.

**Neither harness can measure the difference, by construction.** orion has **zero** uncarried
targets; h1-test's four sit near the bottom of its ranking anyway, because its cosines are
mostly positive. So the gate returned clean and *could not have distinguished the two policies* —
the discrepancy exists only on the panel we actually ship to. It is **unmeasured, not settled.**

**Nothing was changed.** With no surface able to separate the two, swapping one undocumented
choice for another is exactly the move this experiment's own rule forbids. The behaviour is
documented at `--shrink-pert-missing` with the 294/6 split and the 0.13 gap, and the honest
statement of a13's status is: **gated on 294 of 300 rows; six rows unmeasured and unmeasurable
with the surfaces we hold.**

This is a fourth distance in the same family as the other three, and the sharpest yet, because
the name concealed the population rather than the formula: **a rank is meaningful only against
the population it was computed over**, and `damp` named a policy while saying nothing about
which panel it ranked against.


### The six-row divergence does not exist on the real panel

The section above was written from the *code* plus summary statistics, and it is wrong. Run
against the actual object — `work/challenge_donor_cos.npz`, the three production donors against
the real 300 targets, not a simulated distribution — the two branches agree on **all 300 rows**
at a13's settings.

**Why.** Only **17.0%** of the real panel has `donor_cos < 0`, so a zero cosine earns
`rank01 = 0.1800`, not the 0.43 a simulated panel with a fatter left tail produced. And
`0.18 < FLOOR = 0.3`, so `max(0.18, 0.3)` and `max(0.0, 0.3)` are both **0.3**. The floor
dominates both branches and the policies cannot differ.

| `FLOOR` | rows differing | uncarried `z`: gate / builder |
| ---: | ---: | --- |
| 0.00 | 6 | 0.0000 / 0.1800 |
| 0.10 | 6 | 0.1000 / 0.1800 |
| 0.18 | **0** | 0.1800 / 0.1800 |
| **0.30 — where a13 ships** | **0** | **0.3000 / 0.3000** |

**a13's gate and a13's artifact are the same object on all 300 rows.** No builder change is
needed and none was made.

### The lesson, which is not the one the first correction drew

**"These two lines of code differ" and "these two artifacts differ" are separate claims, and
the second needs the data.** Two independent readings — mine and the implementer's — reasoned
correctly from the source and from a true median to a six-row divergence. Running the actual
population through both branches found zero. The code reading was right; the artifact
conclusion was wrong; only the real panel could tell them apart.

**And an actionable bound falls out:** because the floor dominates both branches,
`--shrink-pert-missing` is **unmeasurable at `FLOOR = 0.3`** and only becomes a live question
**below 0.18** — at which point it is unmeasured on every surface we hold. That belongs next to
the flag, and it is worth more than the divergence that turned out not to exist.
