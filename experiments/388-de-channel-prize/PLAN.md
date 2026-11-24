# 388-de-channel-prize — PLAN

## Hypothesis
The scorer's DE channel (per-cell arithmetic-CPM log2FC + sign → nmae/fid/reach/jac)
is computed from per-cell arithmetic means (`mean_calc: arithmetic`, `target_sum 1e6`),
while pds/mse come from bulk sums (`bulk_lognorm`, `bulk_target_sum 5e4`). Our emitter
welds every cell to one library total, tying the two channels. If a future emitter (389)
decouples them — pinning per-gene totals while steering per-gene arithmetic means via
per-cell library layout — the DE channel could carry better signs/magnitudes at zero
bulk (pds) cost. This experiment PRICES that prize on the H1 panel: for each candidate
signed estimator field, what would nmae / fid-proxy / reach-proxy be if the DE channel
carried it?

## Exact inputs (all read-only, all existing)
- `experiments/323-current-family-validation/work/stats-cap15.npz` — delivered E0
  (`lfc` 85×10779, `q`, `p`, `targets`).
- `experiments/323-current-family-validation/work/reference.npz` — `genes` (18080),
  `keep` (10779-universe mask), `targets`, `q`/`lfc` reference matrices.
- `experiments/323-current-family-validation/work/reference_de.parquet` — long-form
  reference DE (target, feature, p_value, p_adj, log2_fold_change).
- `experiments/323-current-family-validation/work/effects.npz` — `cap15` (85×18533)
  designed effect = E1; `genes`, `targets` axes differ from universe (map by name).
- `experiments/172-sign-confidence-key/work/h1test-feats.npz` — `cons` (85×10779)
  signed un-whitened donor consensus = E2. PROVENANCE (verified in source):
  built by `probe172.py` with `SOURCES = ("k562_gwps", "orion_hek293t")` — TWO donors,
  not the three production donors (`fit172.py: PROD = k562_gwps/orion_hct116/orion_hek293t`).
  No h1_2025 (banned oracle), no k562_essential. Axis: recipient = h1-test 85 targets ×
  10779 DE-axis genes — verify by name/order against stats-cap15 before use; if axes
  mismatch, rebuild mapping by name.
- `experiments/169-sharpening-six-member/work/h1test-effect.npz` — fallback (`H`
  85×18080 is POST-whitening; `sums` pre-whitening candidate) if E2 axis fails.
- `experiments/172-sign-confidence-key/work/donorstats.*.npz` — per-donor stats for
  E4 single-donor effects (check keys at runtime; E4 is best-effort: if per-donor
  signed effects cannot be rebuilt leak-free, report E4 as N/A with reason).
- `experiments/348-reach-ordering/reach348.py::walk_reach` — import VERBATIM (copy
  function source, do not modify); known-answer O1 current-panel mean 0.25434726
  (323 `metrics-cap15.json` raw reach 0.25434726355678233).
- cell_eval2 0.16.0 in repo `.venv` (`de_compute.py`, `norm.py`, `streaming_bulk.py`,
  `prep.pseudobulk`) for the mechanism gate.
- 343 scripts `decompose343.py` / `signs343.py` for KA conventions.

## Estimators (signed 85×G fields on the reference-gate gene axis)
- E0 delivered lfc (`stats-cap15` `lfc`).
- E1 designed cap15 (`effects.npz['cap15']`, name-mapped to universe).
- E2 un-whitened consensus (`h1test-feats.npz` `cons`).
- E3 = E1 magnitudes with signs replaced by E2's where they disagree (173 `flip` analogue).
- E4 single production-donor calibrated effects, if rebuildable from donorstats without
  h1_2025/k562_essential leakage (donors allowed: k562_gwps, orion_hct116, orion_hek293t).
- E5 derangement of E2 (permute target labels; control).
- E6 random signs (control; fixed seed).

## Measurements per estimator
1. |l|-weighted sign accuracy on the reference gate (343 convention: gate = reference
   q<0.05, finite nonzero lfc, own target removed, ≥10 genes/target → 58 eligible;
   weight = |l_ref|), plus on gate ∩ our call set (delivered q<0.05) = fid proxy
   (directional precision on the intersection; yield term saturated).
2. Analytic nmae on the gate: (a) estimator's own magnitudes; (b) estimator signs with
   E0 magnitudes; (c) nmae-optimal maps given estimator signs — global scalar,
   per-|est| decile scalars, per-target scalar (343's grid to 6×). Mean over 58 targets.
3. Reach proxy via verbatim `walk_reach`: (i) delivered p-value order + estimator signs;
   (ii) |estimator| order + estimator signs. Report mean reach + mean k*.
4. Paired-target bootstrap, 10,000 draws, seed 20260913, for every delta vs E0
   (follow `342-hamilton-current-family/report342.py` convention).

## Mechanism gate (do early; small synthetic test, cell_eval2's own functions)
Two 400-cell × 2000-gene integer blocks, IDENTICAL per-gene totals, different library
layouts: (a) welded constant total; (b) three library classes at real-control
25/50/75% quantiles (≈14k/20k/27k) with a probe gene's counts concentrated in the
low class. Show `bulk_lognorm` pseudobulk bit-identical while arithmetic-CPM log2FC
differs; report max |Δlog2FC| vs count-split across classes. If cell_eval2 functions
are awkward to call, reimplement both formulas from source and cite line numbers.

## Controls / falsifiers
- KA reproduction first (rule 3): mean NMAE 0.950836, weighted sign acc delivered 0.588
  / designed 0.600, counterfactuals 0.803 / 0.928 / 0.658, O1 reach 0.25434726 — match
  to 1e-4 or STOP and report.
- E5 derangement must give sign acc ≈ 0.5; E6 likewise. If not, the gate/weight code
  is wrong.
- E3 must lie between E1 and E2 in sign acc (it is a sign mixture of the two); if not,
  the flip code is wrong.

## Kill criterion (from brief)
If no estimator beats E0's weighted sign accuracy by ≥ +0.03 AND the best magnitude
map gains < 0.03 raw nmae → verdict "the decoupled channel's H1 prize is amplitude
only (≈ −0.02 raw nmae)".

## Resources / limits
≤ 8 GiB RSS (expect < 2 GiB; all arrays ≤ 85×18533 float64), threads 4
(`NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`, `NUMBA_CACHE_DIR=experiments/388-de-channel-prize/work/numba-cache`),
no 900-block emission, no cell-eval2 CLI full-panel run, no git/submission. Analytic
numpy/polars only + tiny synthetic cell_eval2 calls.
H1 development evidence, one line, 58 targets — not an A/B/C forecast.
