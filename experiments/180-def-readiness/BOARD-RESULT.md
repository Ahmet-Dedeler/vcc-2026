# Board result — `abc-knee` (a23 + per-context knee alpha)

**Entry:** `N3ui4hCkKx9bxExG2pxE`
**Submitted:** 2026-09-02 23:28:07 UTC (Sep-2 slot; description empty)
**Public name:** `We regret to inform the genes` (reused v8 joke; no method in `-m` or `-d`)
**Package:** `submissions/abc-knee.vcc` · h5ad sha256 `f565ce0f5a38ce78b82759f01408aa560c4438d3bd714db0564ed9800647d1d0`
**Preregistration:** `PREREGISTRATION-abc-knee.md` (written 02:25 UTC+3, before submit)

Exactly one mechanism vs a23: `VCC_ALPHA=knee` → A 0.982 / B 1.013 / C 1.048 (a23: global 0.982).

## Scores

| member | a23 raw | knee raw | Δ raw | a23 scaled | knee scaled | Δ scaled | prereg Δ scaled | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| pds | 0.836734 | 0.836786 | +0.000052 | 0.747608 | 0.747726 | **+0.000118** | 0 ± 0.001 | holds |
| mse | 0.872353 | 0.874491 | +0.002137 (worse) | 0.121581 | 0.119322 | **−0.002259** | **+0.0050** | **FALSIFIED** |
| nmae | 0.944330 | 0.943171 | −0.001159 (better) | 0.094133 | 0.096214 | **+0.002081** | ≈ +0.0002 | sign right, 10× large |
| fid | 0.503644 | 0.503819 | +0.000175 | −0.030060 | −0.029477 | +0.000584 | ≈ 0 | noise |
| reach | 0.174219 | 0.177303 | +0.003084 | 0.107435 | 0.110910 | +0.003474 | ≈ 0 | inside seed |
| jac | 0.030105 | 0.030145 | +0.000040 | −0.000538 | −0.000433 | +0.000106 | ≈ 0 | noise |
| **overall** | | | | **0.173360** | **0.174044** | **+0.000684** | +0.0008 to +0.0012 | overall near, **wrong member** |

Displayed rank 24 / 614. a15 remains our best at `0.1740865087` (knee is −0.000043 vs a15).

## What the slot bought

The preregistered falsifier fired:

> mean Δ`mse` < +0.002 or negative ⇒ the per-context knee (or the `C_p` proxy) does not transfer;
> the D/E/F default must revert to a global alpha.

Board Δ`mse` = **−0.00226 scaled**. We raised B and C toward a locally computed knee and the real
scorer made raw `mse` *worse*. `pds` stayed put (scale invariance holds). The +0.00068 overall is
`nmae` + `reach` compensating a `mse` miss — not the mechanism we came to measure.

**D/E/F default: global alpha (a15 / a23's 0.982), not `VCC_ALPHA=knee`.** The 180 local price
(+0.005 scaled `mse`) does not transfer. Do not ship knee blind on Oct 22.

Likely cause, in order of suspicion:

1. The `C_p` NTC-jackknife proxy overestimates remaining claim budget on B/C, so the computed
   knees sit *above* the board's true knees.
2. 096's `claim` route is locally consistent and still not the board functional (rho / gene set /
   jackknife).
3. Mass-fix × per-context scale interaction that the effects-only price missed.

Do not reopen knee as an A/B/C arm. Reopen only with a proxy that retrodicts *this* board point.

## Slot-gate post-mortem

Rule 13 vetoed the arm (`predicted −0.00161` vs `+0.0089` amplitude floor). The overall landed
`+0.00068` — the veto was right about magnitude, wrong about sign, and the priced member (`mse`)
was the one that missed. Another reminder that a small overall can hide a failed mechanism.
