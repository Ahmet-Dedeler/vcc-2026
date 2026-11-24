# Experiment 140 — batch-matched controls

This experiment tests whether target-specific batch composition contaminates
the incumbent donor effects because every target is currently contrasted
against one globally pooled control profile.

The candidate changes only the three production donors:

- K562-GWPS;
- Orion HCT116;
- Orion HEK293T.

The K562 arm is a **processing-replicate diagnostic**, not an entirely new
production mechanism. The live a15 path already replaces the K562 slot with a
50:50 incumbent/PerturbMatch ensemble on genes supported by PerturbMatch, and
PerturbMatch is a limma-voom reanalysis using matched controls from the same
Replogle assay. Orion batch matching is the genuinely new production lever.
Even a passing development result would not be direct authority to modify a15:
the K562 overlap, PerturbMatch composition, unsupported-row behavior, and
production integration would require a separate frozen gate.

For donor target `t` and batch `b`, let `S_tb` be the retained target-cell
count sum and `U_tb = sum_g S_tbg`. The pooled target composition is

```text
S_t / U_t = sum_b (U_tb / U_t) * (S_tb / U_tb).
```

Therefore the exact counterfactual mixture for a pseudobulk formed by summing
counts and normalizing once is **target total-UMI weighting**, not target-cell
count weighting. Experiment 140 preregisters only that mathematically matched
variant. K562 uses `obs["ncounts"]`, which the schema audit found to equal the
analyzed `.X` row sum exactly; `obs["UMI_count"]` includes additional counts
and is not used. Orion uses `total_counts`, which equals the sparse
`gene_expression` sum.

For each supported target:

```text
q_b                  = control_sum_b / total_control_UMI_b
w_tb                 = target_total_UMI_tb / target_total_UMI_t
matched_control_t    = sum_b w_tb * q_b
delta_matched_t      = b_target_t - log1p(50000 * matched_control_t)
delta_eta_t          = (1 - eta) * delta_global_t + eta * delta_matched_t
eta                  in {0, 0.25, 0.5, 1.0}
```

Orion candidate controls use the frozen exact allowlist
`gene_target == "Non-Targeting"`. In particular, `NEGR1` is a perturbation and
is explicitly excluded. The old broad `^neg` detector incorrectly included
154 eligible HCT116 cells and 233 eligible HEK293T cells targeting NEGR1 in
the incumbent global control. For provenance, preparation scans both the
legacy pool (to reconstruct the frozen `b_ctrl.npy`) and the exact
Non-Targeting pool (to construct every matched `q_b`), and records their
measured discrepancy. NEGR1 remains a targeting row when target batch/UMI
mixtures are accumulated.

`b_target_t` is reconstructed from the frozen incumbent artifact as
`delta_global_t + b_control_global`, which makes `eta=0` exactly the incumbent
row and avoids introducing a second target-pseudobulk implementation.
Unsupported target/donor rows have a zero correction and remain bit-exact
incumbent.

All 613 development targets are present in all three production donor tables
and are expected to have supported control mixtures. Within a donor, define
`r_t = b_control_global - b_control_matched,t` and
`mu_r = mean(r_t over supported source rows)`. The production-shaped centered
candidate uses `s_eta,t = s_global,t + eta * (r_t - mu_r)` for supported rows.
An unsupported row receives zero centered correction. This is a
fallback-compatible zero-mean correction rule, but because the development
set does not exercise unsupported rows, any production use on unsupported
targets requires a separate integration experiment.

## Development contract

- Reuse Experiment 132's frozen 613-target development set and 4,992-gene
  response axis.
- Never load or score the frozen 410-target test effects.
- Score K562-GWPS, K562-essential, RPE1, and Nadig Jurkat.
- Treat both K562 datasets as one biological exclusion group.
- Predictions use only the three production donors and exclude any held
  biological source.
- Truth is raw and uncentered.
- Remove all 300 VCC target coordinates and all 1,023 gallery-target
  coordinates before operator fitting, transformation, and scoring.
- Fit the incumbent `union_pruned`, `gamma=0.65` signal-covariance whitener
  separately inside each outer source exclusion, using the exact
  global-control donor artifacts. Reproduce a15's basis rule: one RNG with
  seed 0 consumed in donor order, cap each donor at 12,000 rows, then remove
  selected rows whose perturbation is in
  `data/donors/k562_essential/perts.txt`. The fixed incumbent operator and
  eta-zero donor renormalizations are reused for every eta and null. This
  isolates the matched-control correction rather than allowing
  candidate-dependent covariance or donor-scale changes. The three-donor
  basis must retain exactly 7,784 K562-GWPS, 10,807 HCT116, and 10,803
  HEK293T rows (29,394 total). The K562-held two-Orion basis must retain
  21,610 rows. Each outer group freezes selected-row index hashes.
- Reproduce production precision as well as row selection:
  `Donor.decompose()` means, centered rows, and donor median-norm scales are
  float32; selected centered basis rows are cast to float64 only when written
  into the covariance basis.
- Use complete-gallery raw PDS and continuous signed credit.
- Five deterministic target folds select eta without reading held-fold
  outcomes. Inner selection uses an inner-only gallery; held queries are
  evaluated against the complete 613-target gallery.
- Eligibility is checked at the three-biological-family level (K562,
  RPE1, Jurkat), matching the frozen gate; the two K562 datasets are not
  treated as independent vetoes.
- After primary eta selection, report frozen all-donor, Orion-only, and
  K562-only component diagnostics at the same selected eta values, without
  retuning. They are attribution evidence only and cannot rescue or reject
  the primary candidate.

## Frozen nulls and uncertainty

- `target_label`: one common derangement of the target-specific correction
  rows across all allowed donors. It preserves cross-donor coherence but
  breaks target identity.
- `batch_mixture`: independent donor-specific derangements of the
  target-to-batch-mixture rows. It preserves each donor's correction
  distribution while breaking both target identity and cross-donor
  alignment.
- 49 draws of each null.
- 20,000 target-block bootstrap draws using the Experiment 132 DepMap
  reciprocal-neighbor block definition.

## Gates

All must pass:

1. eta-zero identity is numerically exact at raw-row, centered-row,
   pre-whitening, and post-whitening stages;
2. exact fallback for every unsupported target/donor row;
3. no scored biological family has negative PDS or signed-credit delta;
4. macro PDS delta at least `+0.003`;
5. macro signed-credit delta at least `+0.003`;
6. target-block bootstrap lower 95% is positive for both endpoints;
7. both observed endpoints exceed both null p95 values;
8. positive-gain concentration is at most `0.5` by biological family and
   target block for both endpoints.

Gates 3–8 are evaluated independently on both the nested out-of-fold result
and the final fixed candidate. Each null draw repeats its own nested and
full-development eta selection before comparison.

Preparation revalidates the frozen raw-input hashes. Identity and development
rehash the exact source subsets they consume. Candidate construction uses
control expression and target batch/depth metadata; K562 also streams all
rows only to verify `ncounts`, and Orion reads NEGR1 expression only to
reconstruct the frozen legacy control. Neither check contributes a targeting
expression vector to the candidate correction.

The 410-target test remains sealed regardless of the development result.
Passing development authorizes only a separate review, not an a15 build or
submission.

## Commands

The pre-outcome stages are:

```bash
uv run python -m unittest \
  experiments/140-batch-matched-controls/test_run.py

uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 18 --min-free-percent 15 -- \
  uv run python experiments/140-batch-matched-controls/run.py plan

uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 18 --min-free-percent 15 -- \
  uv run python experiments/140-batch-matched-controls/run.py freeze

uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 18 --min-free-percent 15 -- \
  uv run python experiments/140-batch-matched-controls/run.py prepare

uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 18 --min-free-percent 15 -- \
  uv run python experiments/140-batch-matched-controls/run.py identity
```

Development outcome access is deliberately blocked until a separate
`DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json` is created after Ahmet's explicit
review. No command in this experiment can access the 410-target test.
