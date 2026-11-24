# 370 — Jurkat block-sign prior (PLAN)

## Hypothesis
Context A is Jurkat. Our 2,200-gene CPM block is a coin flip on signs (361:
delivered block-call accuracy 0.508 on H1). The PI measured that on Jurkat's own
high-expression genes 54.0% of knockdown responses are negative, so an all-down
block would be right ~54% on A. A gene-specific prior P(down|gene) estimated
across Jurkat's 2,270 knockdowns could beat all-down if some highly expressed
genes systematically go up under most knockdowns (compensation) while others go
down. This experiment measures, inside the Jurkat bank by target split, the
held-out accuracy of such a prior on block genes — the number that decides how
much fid (and nmae/reach on gate genes) an A-specific prior can buy.

## Exact inputs (read-only; never written)
- `data/donors/nadig_jurkat/delta.npy` (2270 × 8881, float32, mmap), `genes.txt`,
  `perts.txt`, `b_ctrl.npy`. Sign convention pert-minus-control (own-target
  delta < 0 on 99% of rows per mandate; to be re-verified with own code).
- `data/challenge/context_A.h5ad` (18,400 × 18,533 CSR) + `gene_names.csv`
  (header `gene_name`, 18,533 symbols, order TSPAN6/TNMD/DPM1…).
- `data/donors/h1_2025/delta.npy` (255 × 18080; genes == 323 ref genes, verified
  byte-wise 2026-09-10), `b_ctrl.npy`, `genes.txt`, `perts.txt`;
  `experiments/323-current-family-validation/work/reference.npz` (`ref_lin`
  10779-universe basal CPM, `genes`, `keep`) for the H1 block basal.
- `data/donors/k562_gwps/delta.npy` (? × 8248 — shape asserted at runtime),
  `genes.txt`, `perts.txt`, `b_ctrl.npy` for the cross-line check.
- Known-answer records: `experiments/323-current-family-validation/work/
  metrics-cap15.json` (raw fid 0.49755452, nmae 0.95083646, reach 0.25434726)
  + `per-target-cap15.csv`; AGENTS.md §2 (context_A median library ~20,109).

## Estimator
1. **Block for A.** Mean CPM per gene from context_A controls via h5py on
   `X/data|X/indices|X/indptr`, row chunks of 2000 (≤ ~100 MB live): for each
   cell `c`, `lib_c = row sum`, accumulate `data/lib_c*1e6` into gene sums by
   `indices`; divide by 18,400. Record median library size (known-answer:
   ~20,109). Map to Jurkat genes by name. Block = top-2,200 by A-control CPM
   among mapped genes. Per-target eval set = block minus that target's head
   (top-400 |delta| on the donor axis). Exact zeros of delta are dropped from
   accuracy (count reported; expected ~half a dozen entries bank-wide).
2. **Split-half prior (Jurkat primary).** 5 repeats, `rng =
   np.random.default_rng(20260910)`, sequential `rng.permutation(2270)` draws,
   half1 = first 1135 (fit), half2 = rest (eval). Fit per-gene
   `P_down[g] = mean(delta[half1,g] < 0)` and magnitude-weighted
   `Pw[g] = sum(|d|·1[d<0])/sum|d|`. Evaluate on half2 block entries
   (head-excluded): (i) all-down; (ii) gene-prior `P>0.5 → down else up`
   (no exact ties possible for odd n=1135; weighted ties → up, count reported);
   (iii) thresholded `|P−0.5|>τ → prior else all-down`, τ ∈
   {0.05,0.10,0.15,0.20}, with entry- and gene-coverage; same for Pw;
   (iv) shuffled-gene control (P permuted across genes with a sub-RNG seeded
   9000+repeat, then rule (ii)). Report unweighted accuracy (fid reads this)
   and |delta|-weighted accuracy (nmae-on-gate reads this), mean ± sample-sd
   over the 5 splits.
3. **H1 analogue.** Block = top-2,200 `ref_lin` universe genes; head =
   top-400 |delta| per donor row; same 5-split protocol (127/128 halves, same
   seed scheme, fresh RNG). Compare prior vs all-down; reference 361's
   delivered all-down tilt (~0.516–0.518) as context, not as an identical
   quantity (different sign source: donor delta vs delivered lfc).
4. **Cross-line transfer.** Shared axis = Jurkat ∩ K562 donor genes (7,352).
   Full-data priors both lines. K562 eval entries: block = top-2,200 by K562
   `b_ctrl` among shared genes, head-excluded via K562 delta rows, true signs
   from K562 delta; predict with Jurkat full P_down. Vice versa with Jurkat
   `b_ctrl` block and K562 full P_down. Report cross accuracies (unweighted +
   |delta|-weighted), all-down baselines on the same entries, Pearson/Spearman
   of the two P_down vectors and sign-call agreement on shared genes.

## Controls / falsifiers
- **Known answers first, with own code, before any new number:** (a) context_A
  shape 18400×18533 + median library ≈ 20,109; (b) Jurkat own-target
  negativity ≈ 99% (coverage + fraction); (c) mean of `per-target-cap15.csv`
  nmae rows == `metrics-cap15.json` nmae 0.95083646 (valid-cohort
  convention noted; reach NaN set PCCB/TRAFD1/TSKU per 361).
- Shuffled-gene prior (iv) must sit at ~all-down-or-below; if the real prior
  does not beat it, the gene-specificity claim is dead regardless of the
  absolute number.
- Magnitude-weighted vs unweighted priors are reported separately; if only the
  weighted one wins, the gain is an nmae-shape effect, not a fid effect.
- **Kill criterion (mandate): the prior is worth a board arm on A iff held-out
  block accuracy ≥ 0.56 unweighted vs all-down ~0.54.** Also enforced on the
  thresholded variants: a τ-arm must clear 0.56 on entries it actually flips
  (coverage reported), not just on average.
- H1 reading per mandate: what accuracy does a gene-prior get on H1's block
  vs all-down (~0.516), i.e. does the method reach 361's +0.008
  fid-equivalent there.

## Output
- `RESULTS.md` (mandatory format) with accuracy tables (units: fractions,
  mean ± sd over 5 target-splits).
- `work/jurkat_block_prior.npz` (`genes` U-array, `P_down` float64,
  `Pw_down` float64, `n_targets`, `seed`, `block_A` gene names + `block_rank`
  provenance) built on ALL 2,270 targets — production input iff the PI builds.
- `work/metrics370.json` + `work/per-split370.csv` (machine-readable numbers
  behind the tables).

## Compute / memory plan
`memory_pressure | tail -1` ≥ 50% free before the CPM pass (checked 73% at
plan time). Chunked CSR pass ≤ ~0.5 GiB live; donor tables mmap'd or ≤ 350 MB
resident; cap 12 GiB RSS; `NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`; no
densification of any single-cell matrix; delete chunk scratch. No downloads
(no new data), no reserve paths, no builders/submissions, no git writes.
Wall budget: CPM pass ~2–5 min, splits ~minutes; total well inside 3 h.
