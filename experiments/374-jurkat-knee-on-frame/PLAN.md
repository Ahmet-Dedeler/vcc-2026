# 374 — Jurkat knee on the six-line source frame (PLAN)

## Hypothesis
On 341's six-line / 1,024-target / 512-coordinate source frame, blending the
cross-line diagonal baseline `d` with the same-line full-ridge prediction `p`
under 346's fixed-norm rule traces a knee: weighted sign gain saturates early
while cosine-rank retrieval keeps falling. Jurkat's own knee (sign gain and
retrieval loss vs β) sets the context-A arm's pds price via 346's measured
retrieval→pds slope (0.9991 Pearson, 1.009), without any emission or H1 work.

## Exact inputs (read-only, no copies over 2 GB, no downloads)
- `experiments/341-known-line-calibration/selection.json` — frame (1,024
  targets, 512 genes, six lines:
  k562_gwps, orion_hct116, orion_hek293t, rpe1, nadig_jurkat, nadig_hepg2).
- `experiments/341-known-line-calibration/core341.py` — `sign_endpoints`,
  `interval` (SEED 3410909, 10,000 draws). Imported, not copied.
- `experiments/334-text-conditioned-transfer/core334.py` — `endpoints`
  (retrieval with midrank/(n−1) rule, own cosine). Imported, not copied.
- `experiments/341-known-line-calibration/work/<line>/predictions.npz` —
  saved `truth`, `crossline_diag`, `same_line_full` (1024×512 each, ~27 MB
  per line). No refit; the blend acts on these saved rows.
- `experiments/341-known-line-calibration/work/results.json` — known-answer
  record (means table §RESULTS, Jurkat contrasts +.03380218 /
  +.02800431, retrieval −.04829784 / −.01343338).
- `experiments/346-same-line-calibration-h1/RESULTS.md` — slope 1.009,
  fixed-norm blend formula, β=0.25 operating point (context only).

## Estimator
Per line, per target row `i`, with `d = crossline_diag[i]`,
`p = same_line_full[i]`:
```
nd = ||d||2, np_ = ||p||2
m  = (1−β)·d/nd + β·p/np_        (β ∈ {0, 0.1, 0.25, 0.35, 0.5, 0.75, 1})
pred_β[i] = m · (nd / ||m||2)
```
This is exactly 346's rule (unit-blend then rescale to the diagonal norm;
the assignment's "scale ridge to diagonal norm first" wording is algebraically
identical — will assert both forms agree to <1e-14 on one β). β=0 returns `d`
bit-identically by construction (assert max abs diff == 0).
Per β endpoints over all 1,024 targets (both halves, as 341):
- weighted sign + top-64 weighted sign via `core341.sign_endpoints`
- retrieval + own cosine via `core334.endpoints`
Intervals: `core341.interval(Δ)` with Δ = metric_β − metric_0, shape (1,1024)
for Jurkat alone and (5,1024) for the other-five pooled; report mean,
target 95% CI (and crossed where defined). Retrieval Δ × 1.009 → expected
pds raw Δ (346 slope; state as translation, not a score claim).
Efficiency: sign gain per unit retrieval loss = Δwsign / |Δretrieval|.

## Controls / falsifiers
1. **Known-answer gate (must pass before reading any curve):** recompute from
   the saved predictions with the imported functions: Jurkat same_line_full
   weighted-sign mean and retrieval mean; crossline_diag means; the two
   Jurkat Δ means and the pooled retrieval Δ mean. Require absolute
   difference 0.0 (≤1e-12 tolerance stated if any) vs `results.json`.
2. **β=0 identity:** blended matrix at β=0 equals `crossline_diag`
   bit-exactly (max abs diff == 0); β=1 equals the norm-rescaled ridge
   (row norms equal diagonal's to <1e-12).
3. **Shuffled-label sanity (read-only):** report saved `labels_shuffled_full`
   retrieval/sign on Jurkat as a floor reference; not blended.
4. **Kill criterion:** if the known-answer gate fails (>1e-9 on any checked
   mean), stop and report the mismatch instead of any curve.

## What result would kill the idea
- Sign gain on Jurkat at β=1 ≤ 0 vs β=0, or retrieval loss ≥ 0 (no trade-off).
- No knee: efficiency Δwsign/|Δret| monotonic increasing through β=1
  (then the "optimum near small β" premise fails and no β is recommended).
- Otherwise recommend the β with maximal efficiency among β with resolved
  positive sign gain, and give its expected pds raw price; state DE-member
  translation unknown and make no board forecast.

## Budget / safety
≤ 8 GiB RSS, ~6×27 MB + intermediates (<1 GiB); single process,
`NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`; `memory_pressure` checked (76%
free at plan time). Writes only inside `experiments/374-jurkat-knee-on-frame/`
(`run374.py`, `work/curve.json`, `RESULTS.md`). No git, no builders, no
submit, no reserve, no emission/323.
