# 091 board result — the pruned basis is real, and seven times smaller than advertised

Entry `o8rCQvdgjb4I8rSLHKk1`, model name `Sixteen candles`, submitted 2026-08-27 17:15 UTC.
**Overall `0.1503084`** against v11's `0.1501827` — **`+0.0001257`**. Config: v11 exactly,
plus `--whitener-id union_pruned`. Every other parameter identical (α `1.0796`, γ `0.65`,
`min_cpm 10.0`, spread `0.25`, no truncation, coverage 300/300).

## The measurement

The scorer is deterministic and v11 has a byte-identical rebuild, so this is an **exact**
reading of the operator change, not an estimate with error bars.

| member | `a2-prune` | v11 | Δ |
| --- | ---: | ---: | ---: |
| `pds` | 0.736561 | 0.736281 | +0.000280 |
| `mse` | 0.000000 | 0.000000 | 0.000000 |
| `nmae` | 0.094339 | 0.093404 | **+0.000935** |
| `fid` | −0.035087 | −0.035698 | **+0.000611** |
| `reach` | 0.102951 | 0.104056 | **−0.001105** |
| `jac` | 0.003086 | 0.003053 | +0.000033 |
| **overall** | **0.1503084** | 0.1501827 | **+0.0001257** |

Raw: `pds_cosine 0.831773`, `mse 1.5708`, `nmae 0.944139`, `fid 0.502094`,
`reach 0.170083`, `jac 0.031430`.

## Verdict: ship it, and distrust the harness that sized it

**Positive, and it is now our best artefact.** Four of six members improve; `reach` pays.

But experiment 091 measured `+0.000876 avg6` against the unpruned baseline on the Orion
HCT116 harness, and the board delivers `+0.0001257` — **seven times smaller**. It also got
the per-metric attribution wrong: 091 predicted the gain would land in `nmae +0.0032`,
`reach +0.0017`, `fid +0.0012` with `pds` flat. `nmae` and `fid` did improve, `pds` moved
slightly *up* rather than staying flat, and **`reach` moved the wrong way**.

That is exactly the failure mode `docs/HARNESS-DEPTH-CONFOUND.md` predicts: the harness's
control cells are 38% deeper than its perturbation cells, so its reference DE set is 90.8%
down-regulated against ~0.50 on clean data, and any per-metric attribution taken from it —
especially on direction-sensitive members like `reach` — is unreliable. The *sign* survived
the confound, which is what a same-harness A/B with matched row counts should preserve. The
*magnitude* and the *attribution* did not.

## What it cost to find out

The operator itself is not a marginal change — λ₂ falls `1.138 → 0.497` (−56%), trace ratio
`0.9073`, median `|Δλ|/λ` of `0.109` across the spectrum. A 56% change in the second
eigenvalue buying `+0.000126` overall says the whitening operator is close to insensitive to
basis composition at this scale, which is itself worth knowing: further prune-rule tuning is
not where the remaining score is.
