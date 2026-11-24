# 110 — STATE/a14 disagreement arbitration

**Verdict: NO-GO before expensive training.**

The preregistered first gate required row-norm-matched STATE disagreement
routing to have at least `+0.015` oracle directional-credit headroom on both
H1 and HEK293T.

| primary surface | oracle credit delta | gate |
| --- | ---: | --- |
| H1 | `+0.14504` | pass |
| HEK293T | `+0.00686` | **fail** |

The HEK293T value uses the frozen experiment-106 checkpoint and held artifact.
Only about 4.25% of its significant a14/STATE sign disagreements are cases
where STATE has the correct sign and a14 does not. A learned selector cannot
beat this perfect-information upper bound under the frozen replacement
family.

The scaffold also contains a deployment mismatch: it constructs candidate
pairs only on the reference-significant mask, which is unavailable for the
challenge. A future unrelated selector must fit on labeled significant pairs
but apply to every prediction-known mutable disagreement and score all
collateral changes.

No additional checkpoint, count-level score, challenge build, or submission
was run.
