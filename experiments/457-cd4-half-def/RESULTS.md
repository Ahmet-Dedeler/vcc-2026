# 457 — the CD4 ×1.5 twin of 453 under October rules (build log; board read on 09-29)

- Setup (06:06Z 09-26): the operator was seeded from 453 (66610be7, arrays = e2f65ecf), so the refit guard was skipped.
  The key refit and pins passed. Same-line: D = Jurkat, a no-op as in 453.
- **Dump check:** `[cd4] sizing=carrying g=0.516548` (×1.5) against 453's `g=0.344366` (×1); `[h1-union] g=0.56` in both.
  One variable, as preregistered.

| build | factors | ρ D | ρ E | ρ F |
|---|---|---:|---:|---:|
| lvl0 (06:32Z) | 1 / 1 / 1 | 0.8305 | 0.7562 | 0.7213 |
| lvl1 (06:57Z) | 1.0859 / 1.1380 / 1.1652 | 0.9675 | 0.9619 | 0.9603 |
| **lvl2 (07:21Z)** | **1.0931 / 1.1490 / 1.1775** | **0.9794** | **0.9793** | **0.9796** |

The ×1.5 CD4 stack has a little less across-perturbation spread at unit scale than 453's ×1 (0.8407 / 0.7703 / 0.7346). Its
first step is therefore larger: 1.0859 / 1.1380 / 1.1652.

**CONVERGED after 2 levelled builds.** Package `submissions/457dr-october-path-abc-cd4x15-q1.75.vcc`, sha256 `c3afbf8e…`,
novelty rc 0. `gate457.py` passes every gate and wrote `work/gates457.log`, which 09-29 slot 1 requires.
- Config diff vs 453 is exactly `bulk_ctx_scale`, `bulk_load` and `bulk_load_sha256`.
- `sameline` differed only in the per-run artifact path and sha (content identical, arrays equal with NaNs). The first
  gate run flagged that; the twin gates now strip path and sha before comparing.
- The ×1.5 stack levels 0.7 / 0.9 / 1.0 % higher in scale than 453 (1.0931 / 1.1490 / 1.1775 vs 1.0859 / 1.1383 / 1.1662).
  Zero-gene audit 0 (17 cells below −5).
