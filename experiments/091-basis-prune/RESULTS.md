# 091 — Pruning common-essential rows from the whitening basis

**Date:** 2026-08-24 · harness `orion-hct116` count-level · donors `k562_gwps orion_hek293t` ·
grid `g=0.65 q=1.0 a=1.08 cpm=5.0 spread=0.25`, seeds 0/1/2 · two workers, one arm each.

## Verdict

**POSITIVE, and it is a composition effect, not a size effect.**

| arm | run | instance | seeds | mean `avg6` |
| --- | --- | --- | --- | ---: |
| **`BASIS_PRUNE=essential`** | `20260824T193642Z` | `i-0bd65df670234f97e` | 0.12237 / 0.12240 / 0.12271 | **0.122493** |
| `BASIS_PRUNE=random` | `20260824T193646Z` | `i-06197bcff9375a2e1` | 0.12019 / 0.12005 / 0.12001 | 0.120084 |

**delta `+0.002409`, se `0.000121`, t `= +19.85`.** Per-seed ranges do not overlap.

Both arms remove the **same number of rows** (3,810 from the 2-donor HCT116 basis). The only
difference is *which* rows. So basis size is held constant and the entire effect is composition.

Read together with experiment 090 — which *added* unrepresentative rows and lost at t = −13.9 —
this is the same conclusion approached from opposite directions: **the whitening basis is polluted
by common-essential perturbations, and cleaning it is a real operator gain.**

The prune set is safe by construction: the common-essential list is 2,021 genes with **zero
overlap** against the 300-gene challenge panel, so it cannot remove signal the submission needs.

## Provenance — read this before trusting the arm labels

`launch-metadata-*.json` records `run_id`, `instance_id`, `harness`, `grid` and `donors`, but
**not `BASIS_PRUNE`**, which is the only thing that differs between the arms. The two S3 source
archives are byte-identical apart from an extra `launch-metadata` file, so the tree cannot
distinguish them either, and the SSM command history ages out.

The mapping above was recovered on 2026-08-25 from the launch command itself:

```
INSTANCE_ID=i-0bd65df670234f97e ... BASIS_PRUNE=essential bash .../launch.sh
INSTANCE_ID=i-06197bcff9375a2e1 ... BASIS_PRUNE=random    bash .../launch.sh
```

joined to `run_id` through `instance_id` in the launch metadata. **Fix before the next arm-split
experiment: write the arm variable into `launch-metadata-*.json`.** A two-arm experiment whose
arm labels live only in a shell command is one cleanup away from being unreadable.

## What it does not yet establish

The comparison is essential-prune against random-prune at equal row count. That isolates
composition cleanly, which is the scientific question. It is **not** a measurement against the
unpruned production basis on this harness — the `0.12162` figure quoted in session notes comes
from a different run and must not be differenced against these numbers. Two runs, two harness
instantiations; that cross-run subtraction is the same class of error that produced a spurious
`+0.046` elsewhere on 2026-08-24.

## To ship

Production consumes a **pre-fitted** whitener (`--whitener whitener_union.npz`, sha
`5608e807…`). Shipping this means refitting that operator with the essential rows removed and
re-freezing it, which is a separate job from any scalar change.
