# 615 stage 0 — the technical-blob test on Orion RAW HEK293T (fixed BEFORE any raw-data statistic exists)

Written 2026-10-02 ~18:40Z by [21dd4c]. Coordinator's assignment (msg 18:0xZ): build the 614 harness on Orion raw
HEK293T, but test first whether its truth would carry the technical blob that disqualified the thinned Orion
harnesses (FINDINGS H.4). If batch matching does not bring the blob well below the thinned truths' 26–44 %, stop and
report: the gate would be measuring batches.

## The statistic (fixed by known-answer rows, `blobka615.py`, run before this file; no raw data read)
`share = n · ‖mean_t d_t‖² / Σ_t ‖d_t‖²`, with `d_t = b_t − b_ctrl`, `b = log1p(TS · P / ΣP)` of pooled counts (the
scorer's pds space, cell-eval2 bulk TS 5e4), over the gene axis minus the panel genes.

| surface | published (FINDINGS H.4) | this definition (`work/blobka615.json`) |
|---|---|---|
| orion-hek293t-thinned | 44.1 % (mean pair cos +0.423) | **44.2 %** (+0.4233) |
| orion-hct116-thinned | 26.4 % (+0.252) | **26.5 %** (+0.2529) |
| h1-test | 3.6 % (+0.024) | **3.8 %** (+0.0257) |

The a-priori definition (TS 5e4, all genes minus targets) reproduces all three; the variants (TS 1e4, CPM-5 genes) do
not, so they are dropped.

**Decision statistic: the split-half share.** Each target's cells are split at random (seed 615) into halves A and B,
and `share_sh = n · ⟨μ_A, μ_B⟩ / Σ_t ⟨d_tA, d_tB⟩`. Sampling noise inflates the plain share's denominator, so the plain
share depends on cells per target and depth (raw HEK: ≤ 250 cells at 20.7k; thinned: 150 at 14.3k; h1-test: 400 at
20k). The split-half share removes the noise from both terms. The plain share is reported beside it, and the split-half
share is also computed on the three known-answer surfaces for reference.

## Data and target set
- **Source:** all 223 `HEK293T_Batch*.parquet` files, `pass_guide_filter == 1`, controls = `gene_target ==
  "Non-Targeting"` exactly (the builder's `^neg` regex would also catch NEGR1).
- **Gene axis:** Orion tokens mapped to the official 18,533 genes (as `build_orion_harness.py` does).
- **Eligible targets:** A/B/C panel targets with ≥ 100 cells.
- **Harness subset (decides):** a seeded random 100 of the eligible targets (seed 615, drawn from the sorted list),
  each capped at a seeded random 250 cells. This is 614's design; the gate would read this object.
- The all-eligible share is reported beside it.

## Control configurations
| config | control profile for target t | role |
|---|---|---|
| **P** (natural pool) | all controls, all batches | reported |
| **P′** (pool-matched; the buildable harness) | one pool, batch mixture w_b ∝ the harness's perturbed cells per batch, Σ_b w_b · C̄_b | **decides** |
| **BM** (batch-matched) | per target and per half: Σ_b n_tb · C̄_b, the target's own batch mixture of per-batch control means | **decides** |
| G (greedy mechanism) | controls from the first batches in file order until 18,400; the first 150 cells per target in batch order | descriptive: does the old builder's selection reproduce the thinned blob from raw cells? |

C̄_b is batch b's mean control count vector. The scorer takes one control pool, so a buildable harness is P′ (BM is
the per-target ideal).

## Decision (harness subset, split-half share)
- **PASS (build 615):** share_sh(P′) ≤ 0.10 AND share_sh(BM) ≤ 0.10.
  - 0.10 is < 1/4 of the thinned HEK truth's 44 %, < 0.4 of HCT116's 26 %, and ≈ 2.6× h1-test's plain 3.8 %.
- **STOP-and-report, per-target batch composition:** share_sh(BM) ≤ 0.10 but share_sh(P′) > 0.10. A single control
  pool would carry the blob; per-target cell selection would have to be designed first.
- **STOP-and-report, not batches:** share_sh(BM) > 0.10. Batch matching does not remove the blob, so 615 would measure
  a target-independent perturbed-vs-control axis that h1-test does not have.

Reported with no rule attached:
- the plain shares;
- the all-eligible shares;
- mean pairwise cosine;
- G;
- the median-UMI ratio of perturbed to control cells under P′ and BM;
- the top-PC share of the mean-removed energy and that PC's loading on the long-gene axis (the donor tables' dominant
  direction, FINDINGS H.4).

## After a PASS (stage 1, before the gate runs)
The reference is built with `ref614.py` on the P′ truth. Then the 614 low-power diagnostic is reported up front:
- median DE genes per target and the number of nmae-scored targets, against K562's 5 / 34 of 100 and H1's;
- the fraction of the reference's DE calls that are negative (h1-test ~0.50; the old greedy harness 0.908).

If they look K562-like, that is said before the gate's bands are applied. The gate file (TRUSTGATE.md, the 614 bands
copied unchanged) is committed before any gate run.
