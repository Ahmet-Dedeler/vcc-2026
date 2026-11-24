# Preregistered board prediction — `abc-knee` (a23 + per-context knee alpha), written 2026-09-03 02:25 UTC+3 BEFORE any submit

**Package:** `submissions/abc-knee.h5ad`, sha256 `f565ce0f5a38ce78b82759f01408aa560c4438d3bd714db0564ed9800647d1d0`.
**Exactly one mechanism differs from a23:** `VCC_ALPHA=knee` → per-context alphas **A 0.98234 / B 1.01291 /
C 1.04786** (a23: global 0.982). Emitter, mass fix, gate, whitener, weights all identical (180 §3 meta diff).

**Why this slot buys information:** it is the first direct board test of the ρ-knee mechanism with the
amplitude set PER CONTEXT from that context's own controls (`C_p` NTC-jackknife proxy) and our own
prediction's `B` — the exact computation we intend to run blind on D/E/F. The board scorer is
deterministic, so the per-member deltas below are resolvable even though they are small.

## Prediction (deltas vs the a23 board row, scaled units)

| member | prediction | mechanism | falsifier |
|---|---:|---|---|
| `mse` | **+0.0050 scaled mean** (A +0.0006, B +0.0047, C +0.0096) → **+0.0008 overall** | each context moved onto its knee `s* = √(claim/B)`; B and C were 3% / 6.6% below it | mean Δ`mse` < +0.002 or negative ⇒ the per-context knee (or the `C_p` proxy) does not transfer; the D/E/F default must revert to a global alpha |
| `pds` | **0.0000 ± 0.0010** | within-context scale invariance; PLOLO measured −0.000076 raw; small `expm1`/mass-fix nonlinearity | \|Δ`pds`\| > 0.003 scaled ⇒ the profile is not scale-invariant in the way assumed |
| `nmae` | **≈ +0.0002 overall** (slightly positive) | amplitude +1–5% on B/C moves k toward 1 on the cheap side of the kink | Δ`nmae` < −0.003 scaled |
| `fid`, `jac`, `reach` | **≈ 0** (unresolved locally; amplitude moves them by a few %) | rule 13: no reach toll for profile-class arms (169) | any \|Δ\| > 0.01 scaled would be a surprise worth recording |
| **overall** | **+0.0008 to +0.0012** | | negative overall ⇒ NO-GO for `knee` as the D/E/F default |

**What we learn either way:** whether "compute each context's knee from its controls and land on it"
works on the real scorer. That is the assumption the dark-round build rests on.

**Model name:** a joke, no method. **Description:** empty.
