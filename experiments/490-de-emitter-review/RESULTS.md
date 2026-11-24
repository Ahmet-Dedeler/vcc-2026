# 490 — DE emitter core review: where does the code deviate from its documented intent?

**Answer: no new silent score bug in the DE emitter core. The one known deviation (484's q80 mask) costs nothing: fixing
it is slightly worse on H1.** An exact reconstruction of the production emitter prices 484's patch (475's arm) at
**nmae +0.00031 / +0.00034 / +0.00033 raw on 3 seeds (worse)**, with fid, reach and jac 0.0. Overall-equiv is
−0.00009 (CI [−0.00025, +0.00006], P(≤ 0) 0.86–0.89). It fails 475's preregistered rule, which predicted nmae
−0.002 … −0.008. The kf shrink and the MDE add are substitutes: the MDE (c 0.1, cap 0.1) was sized on top of the shrunk
magnitudes, so un-shrinking the protected down genes overshoots the truth.
- One small real shortfall (F1): a fixed 6,000-count capacity margin cuts the prior's dose below 1 on 105 of 900 A/B/C
  blocks (97 in context A) and leaves 1 block (A/CHKB) fully unsteered. It is worth well under 0.001 overall.
- Everything else checked is exact or immaterial, measured on H1 truth and on real A/B/C blocks: the ballast cells' DE
  footprint, order fidelity, sign sources, the coarse MDE delivery, axis mapping, the target's own gene, and call-set
  composition vs jac.

Status: final, 2026-09-26 23:10Z. Analysis only: no build, no `vcc`, no `cell-eval2` CLI, nothing written outside this
directory. Peak RSS 4.37 GB (`ballast_calls_abc.py`, run under the 33 GiB gate at 36.4 GiB free+inactive); every other
process < 1.5 GB.

## Instruments (each with a known-answer row)

| script | what it measures | known answer |
|---|---|---|
| `recon_h1.py` | Exact per-gene intent of `rc_emit_k.py` (`edb5b8dc`) on 435's H1 runs: pass-1 totals (`emit_counts` quasi is deterministic, `floor(n·L·p)`), head / CPM block, q80 threshold, kf, MDE fu/fd, the coarse split | kf, mde_fu, mde_fd, mde_n_up, mde_n_dn, prior_n_flip_target, prior_n_block equal the run's own diag on 85/85 targets, 3/3 seeds (max \|diff\| 0.0) |
| `score_de.py` | The four DE members through cell-eval2's own `prepare_de` + `dispatch_de_metrics`, as `evaluate323.score` does | Delivered tables reproduce the harness's metrics JSON: nmae \|diff\| 1e-16, fid/reach/jac 0 |
| `price475.py` | 484's patch applied inside the reconstruction; delivered 400-cell mean shifted by (n_car/400)·ΔA_carriers | `orig` = the diag; the delivered tables = the harness |
| `mde_by_class.py` | MDE removed per gene class, first order | Removing all MDE gives nmae +0.0097 / +0.0096 / +0.0093; 435 measured −0.00966 for adding it |
| `order_fidelity.py`, the reach walk | Realised vs intended p-order inside the call set | My reach walk reproduces the harness's 0.2693 |
| `audit_lambda.py` | Per-block prior/MDE dose in every recorded rc audit dump (429a, 433a, 432 b0/b1/b2/ef_b1, 448 efr) | — |
| `ballast_abc.py`, `ballast_calls_abc.py`, `prior_delivery_abc.py` | Real A/B/C blocks of `453dr_final.h5ad` (the October-path rehearsal; its DE config equals 446b's: clean4q80, coarse MDE, zero repair −5, ratio ballast), read one block at a time against `data/challenge/context_*.h5ad` | Class quantiles of the controls equal the build log's `[de_sign]` classes |

## X-475 — the q80 fix prices at nmae +0.0003 on H1 (cross-check of 484's known finding, not a re-derivation)

**Where.** `src/vcc2026/rc_emit_k.py:962` (`side = blk & (tgt < 0)` after the q80 restore); comment at
`scripts/make_submission_pm_cd4.py:627` ("the strongest keep their delivered LFC").

**Reproduction and price.** `price475.py` rebuilds each block's intended carrier means twice: with the production mask,
and with 484's `patches/q80_exclude_protected_from_shrink.patch` (protected genes out of `side`; kf, MDE fu/fd and the
coarse split re-derived). The ballast cells are carried over exactly as delivered. Calls, p and reach order are held at
the delivered ones. That is exact for the lfc channel. It leaves out the second-order change in |z| realisation when the
per-class split moves; 475's own run measures that, and its preregistration expects |Δ| ≤ 0.003 on the ordering members.

| seed | nmae (85 targets) | nmae (67 non-essential) | fid / reach / jac | overall-equiv (475's slopes) | CI 95 % | P(≤ 0) |
|---|---:|---:|---|---:|---|---:|
| 0 | +0.000314 | +0.000081 | 0 / 0 / 0 | −0.000087 | [−0.00025, +0.00006] | 0.86 |
| 1 | +0.000340 | +0.000104 | 0 / 0 / 0 | −0.000094 | [−0.00026, +0.00006] | 0.89 |
| 2 | +0.000326 | +0.000083 | 0 / 0 / 0 | −0.000090 | [−0.00025, +0.00006] | 0.88 |

**Mechanism** (seed 0; 2,212 protected-down gate genes, 26 per target, sign accuracy 0.656):
- Un-shrinking them alone (no rebalance): the right-signed ones gain ΣL1 −36.3 (1,450 genes), the wrong-signed ones lose
  +38.7 (762 genes). Net per-target nmae +0.00023.
- The right-signed gain is only about half the shift (0.025 of 0.052 per gene). The delivered value becomes
  E0 − fd·madd, and the MDE add takes |pred| past |l| for many genes.
- The MDE already overshoots this class slightly as shipped: removing it from protected genes alone moves nmae
  −0.0002 (`mde_by_class.py`). In every other class it pays: head +0.0014, flipped +0.0004, kept +0.0017, complement
  +0.0062 (nmae cost of removing it).
- Rebalancing adds about +0.0001: kf median falls 0.53 → 0.32, so the unprotected down genes (78 gate genes per target,
  accuracy 0.62, |E0| median 0.034 against |l| 0.161) are shrunk further.

**Member and direction.** nmae only, and the fix makes it worse (+0.0003 raw on H1; about +0.0005 on the board at 433a's
1.76× magnitude transfer). pds and mse are pinned by per-gene totals.

**Proposed fix.** Keep the code as scored and correct the comment: `patches/builder_q80_comment_and_dose_counts.patch`,
first hunk. Leave `rc_emit_k.py:959` alone: a comment edit changes the pinned emitter sha `edb5b8dc`. If 475's own run
differs from this by more than its ordering-member noise, the difference is the |z|-realisation term this instrument
holds fixed. A variant that exempts the protected genes from both the shrink and the MDE would be a different arm, and
its first-order ceiling is small (~−0.0002 nmae).

## F1 (small, real) — a fixed capacity margin cuts the prior's dose in context A; 1 block goes unsteered

**Where.** `src/vcc2026/rc_emit_k.py:1051` (`margin_ = 6000`), `:1082` (feasibility `Mk.sum(0) <= nk·Lk − margin_`),
`:1084-1099` (λ bisection; λ = 0 means unsteered).

**Mechanism.** The margin is a fixed count per class, but a class's real slack is its sub-gate mass, which scales with
the class's size. Take context A: 66–70 sub-gate counts per cell, 0.37 % of a 17.7k–18.9k library. The lowest library
class has only ~6,400–6,800 counts of slack, so the margin takes ~90 % of it. Raising a gene's arithmetic mean moves
counts INTO that class, which is the sink for up-moves. On the block with the least slack even λ = 1/1024 fails.
Across A's blocks, λ correlates with sub-gate counts per cell at +0.40 to +0.52.

**Prevalence** (`audit_lambda.py`; 433a, the DE configuration 446b ships):

| context | blocks with λ < 1 | λ q10 / median / max among them | blocks below the down-flip threshold 1/(1+kf) ≈ 0.62 | unsteered (λ = 0) |
|---|---:|---|---:|---:|
| A | 97 / 300 | 0.87 / 0.95 / 0.999 | 2 | 1 (CHKB; also in 429a) |
| B | 5 / 300 | 0.47 / 0.85 / 0.98 | 1 | 0 |
| C | 3 / 300 | 0.86 / 0.98 / 0.99 | 0 | 0 |

- Every October rehearsal repeats this in D, the relabelled A: 94–122 blocks with λ < 1, and 0–1 unsteered (AIRE in
  448 efr, ADAMTS18 in 432 ef_b1).
- MDE λ < 1 on 295 blocks and λ = 0 on 2 (433a). 431/426h already documented A's MDE shortfall (A's median MDE dose
  0.956 after `coarse`).
- 446b's meta cannot show any of this: `rc_diagnostics_summary` prints min/median/max only. The `steer_unsteered` key
  it lists is this one block.

**Member and direction.** fid and nmae on at most ~3 blocks' prior signs, plus a ~5–10 % trim of the flip magnitudes on
A's affected third. That is well under 0.001 overall. H1 cannot see it: λ = 1 on 84 of 85 blocks per seed, because H1's
sub-gate slack is larger. So H1 slightly over-prices A for any lever that uses the A-channel.

**Proposed fix.** No emission change is worth a slot on this evidence. `patches/builder_q80_comment_and_dose_counts.patch`
(second hunk, diagnostics only) adds `n_below_1` / `n_zero` counts for `prior_lambda`, `mde_lambda` and
`dose_delivered`, and the `steer_unsteered` block count, to the meta summary. On 433a's rows it reports prior 105 below
1 and 1 at zero, MDE 295 and 2, and 1 unsteered block. The patch dry-runs cleanly against builder `97d6ba31…`. If A's
dose ever matters, the lever is a margin proportional to each class's sub-gate mass rather than a flat 6,000.
`rebalanced` and `unresolved` are 0 on all 900 of 446b's blocks, so placement has headroom. This is untested here, since
it needs emission.

## Checked and negative (measured)

- **Ballast cells vs the carriers' DE design.** The z design, call set and A-channel targets are computed on the
  n_car = 400 − K carriers, while the scorer reads all 400 cells.
  - H1: using the carriers' intended lfc instead of the delivered one (same calls) moves nmae by −0.0001 and
    fid/reach/jac by 0.0.
  - H1 sign changes: 0 in the head and the CPM block, and 3.7 of ~8,200 complement genes per target. Call count
    changes by a median −2 (range −27 … +86), with median 0 leaks.
  - Real A/B/C (30 blocks per context, K median 9–15): median |Δlfc| is 0.0022–0.0026 and q99 about 0.007. There are 3–5
    sign flips per block of ~10,000 genes, and none among the 2,600 highest-CPM genes.
  - Real A calls: blocks with K ≤ 4 are untouched. Blocks with K = 16 gain 13–37 calls the carriers did not make and
    lose 0–14 (net −4 … +32 on ~2,595). At most ~−0.0003 overall through jac's union.
- **Sign source of uncalled genes.** Delivered sign(E0) and design sign(key) agree on 93.3 % of complement genes. On the
  28 truth-significant genes per target where they disagree, key is right 51.0 % and E0 49.0 %.
- **Reach order.** The realised p-order inside the call set follows the documented |H| order: Spearman median 0.964
  (min 0.79). Forcing the exact order moves reach +0.0013 raw, inside H1's reach noise (~0.015).
- **The prior reaches the scorer as designed.** On H1, non-protected block genes carry the CLEAN4 sign 100 % of the
  time; protected genes agree with it 45.7 % and head genes 52.9 %, which gives 0.878 on the 2,200 highest-CPM genes. On
  real A/B/C it is 0.856 / 0.869 / 0.856 (no block below 0.84).
  - The H1 harness prior equals production's pinned CLEAN4 (`7e532f84`) on all 10,777 shared DE-universe genes.
  - Production priors are strictly ±1 (no zero signs to steer a block gene to lfc 0).
- **Coarse MDE.** Non-block moves deliver 83–87 % of the requested add at every expression level (the documented
  fd ≈ 0.8 balance), and no gene is floored to a zero move.
- **The q80 shrink does not starve the kept-down genes.** They are delivered at a median 1.69× E0, because the MDE add
  dominates weak genes.
- **Call-set composition vs jac.** The CPM block's hit rate on the reference's DE set is 0.116, against the head's 0.101
  and 0.059 for the next 2,200 genes by |key|. Set-level jac is 0.0748 for head + CPM block and 0.0378 with the block
  filled by |key|. 192's design matches the metric.
- **Gates and axes.**
  - The emitter's DE universe equals the evaluator-order gate on A/B/C (0 genes differ, build logs).
  - In the glue (`make_submission_pm_cd4.py:5729-5799`), the prior, MDE add, key, zero-repair gate and target column
    all sit on the right axis.
  - Same-line artifacts are shape- and value-checked on load.
  - The per-block seeds (jitter, jk boost, zero repair) are unique per (context, target). The jitter and the jk boost
    share one seed but use it for unrelated draws.
- **Own target.** It sits in the head for 81 of 85 H1 targets and never in the CPM block. That is 484's item 6: one of
  400 head slots.
- **Zero repair on real A/B/C** (`453dr_final.audit428.json`): 0 zero genes remain; 17 cells sit below −5 across 900
  blocks, 1 below −6 (−7.83). These are 484's cap residual and immaterial.
- **Mass-balance sets.** kf acts on every down-signed block gene. That matches 392's comment ("shrinking down-signed
  magnitudes"); only 394's comment says otherwise (X-475). MDE fu/fd act on the up/down eligible sets as documented.
  kf is interior (0.37–0.91) on every 446b block, so the no-root fallback never fired.

## Not checked (territory of others, or out of reach)

- `src/vcc2026/rank_mag.py` (the zs pass): 489.
- Scorer conformance (488). One pointer: cell-eval2 `de_compute.py:432` says gpudge stages per-cell values in float32.
  On H1 a median 268 of our ~2,570 called genes per target have p < 1.2e-38 (173 below 1.4e-45); all are fine in
  float64. If the board's GPU path forms p in float32, those genes tie at p = 0, and reach orders them by |lfc| instead
  of |H|.
- The October day-of scripts (481/485) and 484's own items 2–5 (prior-gate ties, zero-repair edges, the two CPM gates,
  the ballast remainder crash). I only measured item 3's real-data prevalence and item 4's gate agreement.
- The bulk side of the ballast (pds/mse), e.g. ratio mode copying the carriers' floor-rounded composition. 410 and 428
  priced it.
- The |z|-realisation (second-order) part of 475's arm.

## Files

`recon_h1.py` (with `fix=` for 484's patch), `score_de.py`, `price475.py`, `mde_by_class.py`, `ballast_footprint.py`,
`order_fidelity.py`, `audit_lambda.py`, `ballast_abc.py`, `ballast_calls_abc.py`, `prior_delivery_abc.py`;
`work/*.json` (all numbers above); `patches/builder_q80_comment_and_dose_counts.patch`.
`ballast_footprint.py` and `order_fidelity.py` read `work/recon_<tag>.npz`. That file is 67 MB and was deleted to keep
git light; `uv run python recon_h1.py` regenerates it in ~3 s. `price475.py` and `mde_by_class.py` rebuild it in memory.
