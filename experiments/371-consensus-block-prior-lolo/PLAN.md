# 371 — Consensus block prior, leave-one-bank-out (PLAN)

## Hypothesis
370 showed a same-line gene-level sign prior P(down|gene) reaches 0.577 held-out
on context A's high-expression block (vs 0.515 all-down), while a single K562
prior transfers to Jurkat at only 0.548. Contexts B (HeLa) and C (CAL33) have no
same-line bank. A consensus prior pooled over all donor banks might transfer to
an unseen line better than any single bank (errors average out, shared biology
kept). This experiment measures, by leave-one-bank-out over 7 banks, whether the
consensus beats all-down and best-single-bank on held-out-bank block signs.

## Exact inputs (read-only; never written)
- Banks (all `{delta.npy, genes.txt, perts.txt, b_ctrl.npy}`, pert-minus-control):
  - `data/donors/k562_gwps` (9787×8248)
  - `data/donors/k562_essential` (2021×8563)
  - `data/donors/nadig_jurkat` (2270×8881)
  - `data/donors/h1_2025` (255×18080)
  - `data/donors/194-debias-orion_hct116` (17767×38606, debiased/thinned)
  - `data/donors/194-debias-orion_hek293t` (17977×38606, debiased/thinned)
  - `experiments/166-donor-acquisition/donors/nadig_hepg2` (2160×9623; layout
    verified same 5 files, 2026-09-10)
- Shared axis by name: 7-way intersection = **6346 genes** (sorted, deterministic;
  6-way without hepg2 would be 6683 — hepg2 costs 337 genes / 5.0%, noted).
- Challenge: `data/challenge/context_{A,B,C}.h5ad` + `gene_names.csv` (A-block
  definition reuses 370's method).
- Known answers: `experiments/323-current-family-validation/work/
  metrics-cap15.json` (nmae 0.95083646, fid 0.49755452, reach 0.25434726) +
  `per-target-cap15.csv`; 370's context_A median lib 20109, A-block cutoff
  130.4 CPM, Jurkat own-target negativity 0.99707.

## Estimator
1. **Shared axis + per-bank P(down).** P_b[g] = mean(delta<0) over all targets on
   shared columns (exact zeros count as not-down; zero fractions reported at
   full-axis / shared / block-entry level — the debiased Orion check).
2. **Block_L** = top-2200 genes by L's own `b_ctrl` among shared genes.
   Head_L(t) = top-400 |delta| on the bank's FULL donor axis for row t (370's
   method), intersected with the block. Eval entries = block − head − exact
   zeros. Truth = sign of L's delta.
3. **LOLO priors (fitted WITHOUT L):** (a) all-down; (b) each other bank's own
   P(down) (6 single-bank transfers); (c) consensus = equal-weight mean of the 6
   others' P; (d) target-count-weighted mean; (e) agree-only consensus (genes
   where all 6 others agree on the side of 0.5, else all-down; entry+gene
   coverage reported); (f) shuffled-gene consensus control (5 draws, seeded).
   Ties (P==0.5) → predict up; tie counts reported.
4. **Endpoints:** unweighted accuracy (fid reads this) + |delta|-weighted
   accuracy (nmae-on-gate proxy), per held-out bank and pooled (micro over
   entries + macro over banks).
5. **Uncertainty:** 1,000-draw target bootstrap per bank (resample targets with
   replacement; per-target correct/eval sums precomputed) for the
   consensus-equal-minus-all-down difference, unweighted and weighted: mean,
   2.5/97.5 percentiles, frac>0. Seeds deterministic per bank.
6. **A/B/C + A sanity:** chunked-h5py mean CPM for all three contexts (370 code);
   consensus coverage of the challenge axis (expect 6346/18533 = 34.2%) and of
   each context's challenge-wide top-2200 block. A sanity: full-consensus (all 7,
   in-sample for Jurkat — flagged) and LOLO-consensus (no Jurkat — the clean
   transfer number) scored on (i) Jurkat's LOLO block and (ii) 370's exact A-CPM
   block restricted to covered genes, vs 370's same-line 0.577.

## Controls / falsifiers
- Known answers first with own code (323 means; A shape/median-lib/gene order;
  own-target negativity per bank incl. Orion debiased + hepg2 — verifies the
  sign convention on every bank used).
- Shuffled control (f) must sit at/below all-down; if consensus does not beat
  it, gene-specificity is dead.
- **Kill criterion (mandate): consensus is worth a B/C rider only if its LOLO
  gain over all-down is ≥ +0.02 unweighted on the MEDIAN held-out bank with the
  bootstrap interval above zero; otherwise record the negative.**
- Single-bank matrix (b) is the control for "consensus vs best single": if the
  best single donor ties/beats consensus, pooling adds nothing.
- H1 bank (255 targets) is expected noisy — its CI will show it; no separate
  conclusion drawn from it alone.

## Outputs
- `RESULTS.md` (mandatory format).
- `work/consensus_prior.npz` (`genes`, `P_down` equal-weight all-7 consensus,
  `n_banks=7`, plus `P_down_weighted`, `P_banks`, `bank_names`, `n_targets`,
  provenance) for the PI regardless of verdict.
- `work/metrics371.json` + `work/per-bank371.csv` (machine-readable numbers).

## Compute / memory plan
`memory_pressure` 75% free at plan time. All `delta.npy` mmap'd; big banks
streamed in 1000-row chunks (≤0.5 GiB transient); shared-column copies ≤456 MB
one bank at a time, freed after eval; cap 12 GiB RSS;
`NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`; no densification of any single-cell
matrix (chunked CSR pass as in 370); no downloads, no reserve paths, no
builders/submissions, no git writes. Expected wall: minutes (370 ran in 6.9 s
compute; Orion streaming + bootstraps dominate, budget <1 h).
