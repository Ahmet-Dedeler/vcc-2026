# 459 — the H1 ×3 twin of 453 under October rules (build log; board read on 09-30)

- Setup (07:55Z 09-26): operator seeded from 453 (66610be7); same-line D = Jurkat, a no-op.
- **Dump check:** `[h1-union] g=0.84` (×3) against 453's 0.56 (×2); `[cd4] sizing=carrying g=0.344366` in both. One variable,
  as preregistered. Dump `b12ef388…`.

| build | factors | ρ D | ρ E | ρ F |
|---|---|---:|---:|---:|
| lvl0 (08:22Z) | 1 / 1 / 1 | 0.8742 | 0.8039 | 0.7651 |
| lvl1 (08:47Z) | 1.0584 / 1.1037 / 1.1314 | 0.9712 | 0.9655 | 0.9631 |
| **lvl2 (09:12Z)** | **1.0584 / 1.1122 / 1.1416** | **0.9712** | **0.9794** | **0.9794** |

H1 ×3 adds across-perturbation spread (453's lvl0: 0.8407 / 0.7703 / 0.7346), so the first step is smaller: 1.0584 / 1.1037 / 1.1314.

**CONVERGED after 2 levelled builds** (A was already within the loop's 1 % tolerance at lvl1). Package
`submissions/459dr-october-path-abc-h1x3-q1.75.vcc`, sha256 `4ac6a81f…`, novelty rc 0. `gate459.py` passes every gate and
wrote `work/gates459.log`, which 09-30 slot 2 requires.
- Config diff vs 453 is exactly `bulk_ctx_scale`, `bulk_load` and `bulk_load_sha256` (`sameline` content identical).
  Zero-gene audit 0.
- **Execution note for the read:** A sits at ρ 0.9712, 0.8 % below 453's 0.9793 (0.4 % in scale). The cliff model prices
  that at about −0.0001 overall (mse slope 0.31 per unit scale, one context), small against the +0.0004 prediction.
  Subtract it when reading the H1-axis effect.
