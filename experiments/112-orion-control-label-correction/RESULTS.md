# 112 — Orion control-label correction

**Verdict: NO-GO by exact algebra.**

Including `NEGR1` in the Orion control pool translates every donor effect row
by one common vector. Both `Donor.decompose()` and the union-whitener builder
subtract the donor row mean before the values can affect production:

```text
(delta + c) - mean(delta + c) = delta - mean(delta)
```

Because production uses centered donor effects and `beta=0`, correcting the
control label cannot change the prediction or the fitted whitener except for
floating-point roundoff.

No expression data, held truth, count-level scorer, cloud worker, challenge
build, or submission was used.
