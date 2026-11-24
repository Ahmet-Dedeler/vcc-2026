# 488 — Does the emitter's model of the scorer match cell-eval2?

**No, not everywhere.** Against installed **cell-eval2 0.16.0** (scanpy 1.12.3, scipy 1.18.0; `auto` resolves to **scanpy** because pdex and gpudge are not installed), the rank-controlled DE path matches the scorer. The builder's CPM gate does not, and an exact-U move does not leave Wilcoxon p unchanged.

CONFIRMED, on synthetic matrices run through `cell_eval2` itself (`work/conformance.py`, output `work/report.json`):

- `rc_emit.de_from_block` and `rank_mag.block_lfc` match `compute_de` on p, BH q, and log2 fold change (max |Δp| 1.1e-23, max |Δlfc| 8.3e-16).
- `scorer_control_block` and `control_block_from_counts` use the scorer's control gate: unweighted mean of per-cell CPM, reference cells only, strict `> 5`.
- `emission_rank_axis`, `profile_from_effect`, `control_render_gate`, and `TransferModel.cpm_gate_mask` threshold **pseudobulk** CPM. On a two-cell control they keep the gene the scorer drops and drop the gene the scorer keeps.
- `rank_mag.postpass` keeps the rank sum and moves p through the within-block tie correction. A hand-built tie split moves p by a factor 1.074 at fixed U and flips a BH call at m = 9501. The shipped move on a 4-cell block moves p by relative 3.0e-4.

This is a reproduction, not a new mechanism. Experiment 214 already scored the tilt on the scorer's mean-per-cell universe rather than the builder's pooled gate. Experiment 423 already recorded that the tie term moves p, that H1 had 0 call flips, and that raw mse stayed bit-identical because the prediction sits on the trace cap. gpudge versus scanpy (~1e-5 in `cell_eval2/run.py`) is **UNTESTED** here.

Re-run from the repo root:

```bash
uv run python experiments/488-cloud-scorer-conformance/work/conformance.py
```

No challenge h5ad is required. The suite is Python; it does not depend on macOS bash 3.2 or BSD grep.

## Versions

| Piece | Value |
|---|---|
| cell-eval2 | 0.16.0 at `.venv/lib/python3.12/site-packages/cell_eval2` |
| preset | `configs/vcc2026.yaml`: target_sum 1e6, filter 5.0, wilcoxon, arithmetic mean, epsilon 1e-9, fdr_scope per_pert, bulk_target_sum 5e4, cosine, rank denominator n-1, exclusion panel, tie policy midrank, backend auto |
| DE backend used | scanpy (`de_compute.py` calls `sc.tl.rank_genes_groups(..., method="wilcoxon", tie_correct=True)` at line 872) |
| Not installed | pdex, gpudge. `run.py` lines 1001–1002 say those engines differ by ~1e-5. **UNTESTED** numerically |

## Assumption table

Verdict is CONFORMS / DIFFERS / UNTESTED against cell-eval2. Evidence is CONFIRMED when `work/report.json` holds the number, SUSPECTED when the size is reasoned from a synthetic draw or a prior experiment.

| Assumption | Emitter | Scorer | Verdict | Evidence |
|---|---|---|---|---|
| Control gate is mean per-cell CPM, reference only, strict `>` | `rank_mag.scorer_control_block`, `rc_emit.control_block_from_counts` | `de_compute._ref_cpm_from_cells`, `_apply_cpm_filter` (`ref_cpm > threshold`), then BH on survivors | CONFORMS | CONFIRMED. Two-cell matrix: mean CPM of S is 500, of L is 3. Both helpers and `compute_de` keep S and filler, drop L. Gene at CPM exactly 5 is dropped. A gene at 0 in control and ~5e5 in the perturbation is dropped |
| Same gate, builder render path | `emission_rank_axis`, `profile_from_effect`, `control_render_gate`, `TransferModel.cpm_gate_mask` | same | **DIFFERS** | CONFIRMED. Same matrix: pseudobulk CPM of S is 0.999, of L is 5.994. All four keep L and drop S. One-hot effect through `profile_from_effect`: S's profile does not move (max abs 0); L's does (1.20e-5) |
| How often the two gates disagree on a real panel | — | — | UNTESTED on A/B/C | SUSPECTED size. Synthetic 2000×6000, libraries clipped to the challenge depth range. Planted depth–composition covariance: scorer 2755 genes, pseudobulk 2739, only-scorer 43, only-pseudobulk 27, and all 43 only-scorer genes have mean CPM in (5, 50]. Shared composition: 6 and 2. Context h5ads are not on this VM |
| `count / (lib/1e6)` versus `count * (1e6/lib)` | `rank_mag._normalized` is the first; `rc_emit` is the second | scanpy `normalize_total` after the float64 upcast is `count / (lib/1e6)` | CONFORMS for ranking | CONFIRMED. On a random 40×25 integer matrix, rank_mag is bit-identical to scanpy. rc_emit max abs difference 5.82e-11, 0 equality flips, 0 order flips |
| Wilcoxon: midranks, tie correction, no continuity correction, `p = 2 Φ̄(\|z\|)` | `rc_emit.wilcoxon_z`, `de_from_block` | scanpy wilcoxon, `tie_correct=True` (`_rank_genes_groups.py` 456–500). No continuity term | CONFORMS | CONFIRMED. 8 gated genes: max \|Δp\| 1.10e-23, max \|Δq\| 8.77e-23. scipy `mannwhitneyu` **with** continuity differs by up to 7.31e-3; without it, 1.10e-23 |
| BH is scipy `false_discovery_control(method="bh")` on the post-gate genes, per perturbation | `de_from_block` | `_apply_cpm_filter` drops rows, then `_bh_per_target`. Scanpy's pre-filter `pvals_adj` is discarded | CONFORMS | CONFIRMED. Below-gate g1 has scanpy q 3.73e-18 and is absent from `compute_de`. In-gate g0 moves from scanpy 3.23e-9 to cell_eval2 4.04e-9 because m shrank. Per-pert BH of `[1e-8, 0.04]` and `[0.04, 0.04]` is `[2e-8, 0.04, 0.04, 0.04]`; global BH makes the first `4e-8` |
| log2FC is `log2((mean_pert + 1e-9) / (mean_ctrl + 1e-9))` of arithmetic means of per-cell CPM. A zero pert mean is finite | `de_from_block`, `rank_mag.block_lfc` | `vcc2026.yaml` epsilon 1e-9, mean_calc arithmetic | CONFORMS | CONFIRMED. Zero pert mean gives −43.185065, equal to `log2(1e-9 / (ref + 1e-9))`. rank_mag lfc max \|Δ\| 3.89e-16 |
| Reach orders significant genes first (`p_adj < 0.05`), then p_adj, p, \|lfc\| descending, gene name ascending. Depth is the deepest k with purity ≥ 0.9 | design assumption, not a reimplementation | `metrics/direction.py` `_purity_curve`, `REACH_PURITY_FLOOR = 0.9` | CONFORMS | CONFIRMED. Nine genes, one wrong sign, tied p and p_adj. Wrong gene first by \|lfc\|: reach 0. Wrong gene last: 0.8889. Equal \|lfc\|, name `g0` leads (0) and `g8` trails (0.8889) |
| A call is `p_adj < 0.05`. A significant gene with predicted lfc 0 stays in `n_pred` and is a miss. An adjudicable gene the reference did not call still enters `n_pred` | design assumption | `de_direction_fidelity_yield_raw` = k / max(n_pred, n_real). `_defined` is non-null, non-NaN, non-zero; inf has a sign | CONFORMS | CONFIRMED. p_adj exactly 0.05 on the only reference-significant gene: fid 0. p_adj 0.049: fid 1. n_real = 1, extra call of an adjudicable gene at lfc 0: fid 0.5; without that call: fid 1. The earlier two frames in the suite do not separate the claim (both sides land on 0.5 or on 1/3) |
| nmae fills null, NaN, non-finite, and absent predicted lfc with 0 and keeps them in the gate. Reference non-finite lfc and the target gene leave. Gate under 10 genes is omitted | design assumption | `de_lfc_nmae`, `min_gate_size` 10 | CONFORMS | CONFIRMED. Hand calculation 8/12 = 0.666667, returned 0.666667. A 9-gene gate returns `{}` |
| Jaccard of two empty sets after dropping the target gene is 1 | design assumption | `de_sig_jaccard` | CONFORMS | CONFIRMED. Both sides significant only on the target: 1.0. Prediction calls one extra gene: 0.0 |
| exact-U moves leave the scorer's p unchanged | `rank_mag.postpass` asserts U, not the tie term. The docstring already says p can move | tie correction is inside the Wilcoxon sd | **DIFFERS** | CONFIRMED. See below |
| `block_jackknife` equals `moments.jackknife_correction` on a legal block, and is invariant to how counts are arranged at fixed totals | `scripts/make_submission_pm_cd4.py` `block_jackknife` | `moments._loo_bulk`: `r_i == 0` maps that row's v to 0, then subtracts the full-sample bulk | CONFORMS on legal blocks. **DIFFERS** on `r_i == 0`, and C is not allocation-invariant | CONFIRMED. Legal 25×12: relative difference 1.26e-16. Same column sums, two cells swapped: C 0.4416 → 0.4612. One cell holds every count: cell_eval2 C = 29.267, builder NaN |
| PDS is cosine, scale-invariant, panel exclusion, midrank, denominator n−1. A control paste scores 0.5 | submission geometry, not an emitter formula | `discrimination_score` with the vcc2026 kwargs | CONFORMS | CONFIRMED. Identical composition at 3× counts: PDS 1, 1, 1 on both. Paste: mean 0.5. A prediction whose only effect is on the other panel genes: 0.5 under panel exclusion, 0.0 under row exclusion |
| Official GPU scoring (gpudge) matches these scanpy numbers | — | `run.py` 1001–1002 | UNTESTED | source text only |

## DIFFERS 1 — pseudobulk CPM versus mean per-cell CPM

Reproduction is `cpm-gate-mean-vs-pseudobulk` in `work/report.json`.

Two control cells. Gene S has 1 count in the library-1000 cell and 0 in the library-1e6 cell. Gene L has 6 counts in the large cell only.

| Gene | Mean per-cell CPM | Pseudobulk CPM | Scorer / rank_mag / rc_emit | Builder gates |
|---|---:|---:|---|---|
| S | 500 | 0.999 | keep | drop |
| L | 3 | 5.994 | drop | keep |
| filler | 999497 | 999993 | keep | keep |

`compute_de` features are S and filler. `control_render_gate` was called on a real two-cell h5ad, not just the formula copied out of it. `TransferModel.cpm_gate_mask` at threshold 5 matches the builder. Its docstring says genes below the mask are excluded from the DE metrics because the filter is 5 CPM. That is false for S: mean CPM 500, pseudobulk 0.999. The donor default threshold in `profile()` is 20, which is experiment 006's knob, separate from this mismatch.

`profile_from_effect` returns a renormalized abundance, so applying an effect to every gene moves genes the gate zeroed. A one-hot effect does not: max abs change is 0 for S, 1.20e-5 for L, 4.42e-6 for filler.

**Score effect.** SUSPECTED, not a panel count. Genes the scorer keeps and the pseudobulk gate drops sit in mean CPM (5, 50] in the synthetic draw (the 43 only-scorer genes have mean-CPM quantiles 5.03, 5.50, 5.96, 6.05). That band is below `RC_CALLSET_CPM = 2200`, so a zeroed bulk effect misses the |H| head (top 400) and falls through to the tail law, while the gene remains inside the Wilcoxon universe. The opposite set gets a bulk effect the scorer drops before BH. Direction: lost tested signal on a thin low-CPM band, and untested effect on genes the scorer never scores. Do not flip the default without a score. Experiment 214's 5–10 CPM note is the same family, measured there as a 10-versus-5 threshold rather than as two definitions of "5".

**Fix.** `patches/mean-cpm-gate.patch`. Optional `mean_cpm` on `profile_from_effect`, `emission_rank_axis`, and `cpm_gate_mask`, and `per_cell=False` on `control_render_gate`. Unset arguments keep today's pseudobulk gate, so an unflagged build stays bit-identical. The per-cell formula, checked against `_ref_cpm_from_cells` on this two-cell matrix, matched to 0. `patch --dry-run -p1` succeeds. Not applied.

These functions cannot be fixed by swapping the formula in place. Their input is already a pseudobulk vector (`b_ctrl`, or summed counts). The scorer's vector has to be threaded in. `RcControlBlock.ref_lin` is that vector on the kept genes; the full-axis mean is what the new argument wants.

## DIFFERS 2 — exact-U does not leave p unchanged

`rank_mag.postpass` keeps each carrier inside the same open interval between control values, so the number of controls below it, and whether it ties a control, stay fixed. Cells inside the block can untie from each other. That changes the tie-correction factor in the Wilcoxon standard deviation and therefore z and p. U can stay exactly the same. The module docstring already says this (lines 26–30) and then says the calls are unchanged. The second sentence holds only outside a thin window around q = 0.05.

Reproduction is `exact-u-pvalues`.

Hand-built, six perturbation cells inside one control gap, same total, same mean, same arithmetic-mean lfc (−1.736966):

| | Tied at count 100 | Split 50 / 150 |
|---|---:|---:|
| `compute_de` p | 4.898753e-06 | 5.262684e-06 |
| scipy Mann–Whitney U, no continuity | 168 | 168 |
| BH q at m = 9501 (one real p, the rest 1) | 0.04654 | 0.05000 |

The tied arm is a call. The split arm is not. With one real p and m−1 ones, q = min(p·m, 1), so any pair of p-values that straddle 0.05/m flips. m = 9501 is that stratum, not a claim about our call-set size.

Shipped `postpass` (mode zs, selection dn) on a 4-cell block: moved, `dU_car_abs_sum = 0`, gene 0 counts `[3, 12, 0, 0]` → `[1, 14, 0, 0]`. `compute_de` p 0.639110 → 0.639304 (relative 3.03e-4). lfc −6.05889 → −6.47393, which is the move the pass is for. A call flips only if q sits within about 1.5e-5 of 0.05. That matches 423's 0 call flips.

Jackknife of that same integer move, column sums fixed: at n = 4, C goes from 0.971 to 2.870. At n = 400 with 398 identical ballast cells, C goes from 5.385e-3 to 7.035e-3 (relative +0.306). The relative figure is the toy's: C is small because the other cells match. Experiment 423 measured about +7% median jackknife on the production bulk and **raw mse bit-identical**, because `expr_mse` takes min(prediction trace, k × real trace) and the prediction already sits on the cap. Scored mse was not re-run here. There is no harness on this VM.

**Fix.** `patches/exact-u-tie-term.patch` corrects the docstring. It does not change the pass. A guard that refused a move when the within-block tie multiset changes would also refuse the lfc move. 423's nmae gain is that lfc move (−0.00232 raw on H1). Leave the pass alone.

## DIFFERS 3 — jackknife when one cell holds every count

On any block where every leave-one-out remainder is positive, `block_jackknife` matches `jackknife_correction` to relative 1.26e-16 (`jackknife-matches`). The variance is not a function of the column sums alone: swapping two cells' compositions at fixed margins moved C by 0.0196.

When one cell holds every count, `moments._loo_bulk` sets that row's v to 0 (because the remainder is 0) and subtracts the full-sample bulk. `block_jackknife` divides by the remainder and returns NaN. Measured: cell_eval2 C = 29.267008727749964, builder NaN. Production blocks give every cell a positive library, so this edge is not on the scored path.

**Fix.** `patches/block-jackknife-zero-remainder.patch`. Map `r == 0` to v = 0, subtract the full-sample bulk, then take the same centered sum of squares. Checked in the suite (`jackknife-proposed-formula`): the edge matches cell_eval2 exactly, and the legal 25×12 block matches to relative 1.26e-16, the same residual as today's function. `patch --dry-run -p1` succeeds. Not applied. Score effect on current emission: none.

## What was not claimed

- The synthetic disagreement rate is not a count of genes on contexts A, B, or C.
- The +30% jackknife at n = 400 is not a board mse delta. 423 is the measurement that mse did not move, and it is not re-done here.
- No gpudge or pdex number was produced. Every DE comparison in this directory is the scanpy backend.
- Production files were not edited. The three patches live in `patches/` and were dry-run only.

## Coordinator's disposition (2026-09-27 06:50Z; brought to master from `cloud/488-scorer-conformance`)
Accepted as a conformance record; no action. The emitter's DE statistics match cell-eval2 (p, q, lfc). The builder's
pseudobulk-CPM gates versus the scorer's per-cell-mean gate differ only on near-threshold genes, which is immaterial at
18,400 control cells (214 already scored on the scorer's universe). The exact-U pass's tie-term p movement is 423's known
item (0 call flips on H1).
