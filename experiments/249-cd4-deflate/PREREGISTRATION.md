# PREREGISTRATION — `cd4svd5`: deflate the shipped CD4 residual's top-5 target-shared directions. Written 2026-09-04 13:05 UTC, BEFORE any build.

Lead. Committed before the builder is touched. Not editable after the board result.

## 1. The arm — one variable
On top of the shipped stack (`h1union`, entry `iUziqE1gewfcPKzYUDaB`, board 0.2118310, rank 3), the CD4
residual is **centred over its carried rows and its top-5 left-singular directions removed** (246
`lib246.deflate(mode="svd", k=5)`: `R_c = R − mean`, `G = R_c R_cᵀ`, remove `U_k U_kᵀ R_c`, support
re-zeroed) BEFORE `solve_cd4_fraction` sizes it — so the existing rule sizes the object actually added
(246's `B2_cd4_svd5_resolved`; PLOLO fraction 0.374 → 0.396). `VCC_CD4_TARGET_F` stays **0.34**,
`VCC_H1_UNION_G` **0.28**, the H1 union untouched, α 0.982, operator `e2f65ecf…`, emitter `rc_callset`.
New flag `VCC_CD4_DEFLATE_K=5`, unset = byte-inert.

## 2. Why
The CD4 residual enters un-whitened at `f` 0.34, so nothing ever removes its target-shared
component — the T-cell-programme hub the whitener strips from the donors (209) and nothing strips
from CD4. Its top-5 directions carry 20.2% of its energy across 251 rows. 246 measured removing them on
162's PLOLO instrument, three cancer folds: **`pds` +0.01226 raw (3/3 folds, 8.2 LOO-σ) for −0.46%
aggregate / −0.81% macro cosine** (resolved variant). Controls: energy-matched random subspaces
(same 21% of energy) buy **no** `pds` (−0.0001 / −0.0050) and cost 3× the cosine; shuffled projectors
≤ +0.0013; a split-half projector fitted on the other half of CD4's rows keeps +0.0104…+0.0118;
non-carrying rows move ≤ 0.026 pp; leave-top-3-out retains 56%; leave-one-fold-out k selection
+0.0017 / +0.0029 / +0.0070. Exchange rate: 29× the dose axis's `pds` per unit cosine. §8.1d: the
removed component is unmatched magnitude (`e` 0.156), exact `nmae` +0.00031 (magnitude bound).
**This is the first arm that buys `pds`, the member where the leader's lead is largest (+0.0139).**

## 3. Predictions (246 §4, resolved row; `pds` transfer 0.7× is the only calibration point, `h1union`)
| member | predicted Δ overall |
|---|---:|
| `pds` | **+0.0032** (0.7×) · +0.0046 (1.0×) · +0.0091 (2.0×) |
| `mse` | −0.0003 |
| `reach` | −0.0002 |
| `nmae` | −0.0001 … +0.0003 |
| `fid` | −0.0001 |
| `jac` | ≈ 0 |
| **overall** | **+0.0026 … +0.0040**, upside +0.0086 · board **0.2144 – 0.2158** |

Board raw `pds` predicted **+0.0086** raw at 0.7× (h1union's realised transfer) — 0.8357 → ≈ 0.844;
board raw `mse` predicted **≈ 0.8183** (cosine −0.46% through 197's identity at the arm's own terms;
243: measure it with the package's `A_c`/`B_c`, never the floor).

## 4. Falsifiers — pre-committed
1. **Board raw `pds` rises by less than +0.0040** (half the 0.7× prediction): the mechanism does not
   transfer; retire deflation as an arm class.
2. **Raw `mse` worsens by more than +0.0030** (three times the predicted +0.0009): the cosine cost is
   not the −0.46% PLOLO read.
3. **Overall Δ outside [−0.0010, +0.0100].**
4. **Config diff is not exactly `cd4_deflate_k` plus its deterministic consequences** — the CD4 frame's
   `residual_median_row_norm_at_unit_fraction`, `f_per_unit_fraction`, `solved_fraction` and the
   merged preserve scalar move; `target_f_requested` 0.34, `h1_union_g` 0.28, `alpha` 0.982,
   `operator.sha256`, every `effects.h1_union.*` field and the union's solved fraction do NOT.
5. **Delivered amplitude moves by more than 0.5%** (246: 1.0002 vs the stack's 1.0033) — the arm is
   claimed shape-only; if amplitude moved it is not.

## 5. What the slot buys even if flat
Whether a `pds`-first arm transfers at all: every prior direction arm paid on `pds`; this is the first
priced to gain it, on the member that separates us from first. A flat result kills the class; a
positive one opens k = 3–10 and the `+partial` composition (246 §3.7, +0.0032 … +0.0100).

## 6. Discipline
Joke model name, **empty description** (§14). Intended for the SECOND 2026-09-05 UTC slot, after the
armed 237 reading; never in flight together with it.
