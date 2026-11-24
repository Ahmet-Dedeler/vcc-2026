# 345 — Common-program structure: does a target-shared vector carry sign information? (PLAN)

## Hypothesis
The reference DE on a held-out line has a large target-shared component (a line-level common
response: proliferation/ribosome/stress programmes moving the same way under most knockdowns).
Our design is optimised for discrimination (whitening strips shared directions), so it may carry
little of that component — and that component may be most of the |lfc| mass on which
nmae/reach/fid are scored. If a common vector estimated WITHOUT the target's own truth predicts
reference signs materially better than 0.60 weighted, there is a design-stage lever needing no
new biology.

## Inputs (exact paths)
- `experiments/323-current-family-validation/work/reference.npz` — `lfc` 85×10,779 (universe),
  `q`, `genes`+`keep` (18,080 → 10,779 universe), `targets`, `ref_lin` (basal CPM on universe).
- `experiments/323-current-family-validation/work/reference_de.parquet` — long table
  (`target, feature, p_adj, log2_fold_change`); authoritative gate source (same as 343).
- `experiments/323-current-family-validation/work/stats-cap15.npz` — delivered `lfc`, `q`, `p`
  85×10,779.
- `experiments/323-current-family-validation/work/effects.npz` — `cap15` design 85×18,533 on
  `effects.npz['genes']`; map by gene NAME onto the universe; 2 universe genes absent → design 0
  (same as 343 `signs343.py`).
- `experiments/323-current-family-validation/work/metrics-cap15.json` — known-answer NMAE 0.950836.
- Transfer banks (step 4b): `data/donors/k562_gwps/delta.npy` (9,787×8,248, float32, mmap) +
  `genes.txt`/`perts.txt`; Orion `data/donors/orion_hct116/delta.npy` if same layout.
- Prior record: `experiments/343-nmae-decomposition/RESULTS.md` (+`decompose343.py`,
  `signs343.py`, `bins343.py`), `experiments/155-control-covariance/RESULTS.md` §4,
  `experiments/246-direction-hypotheses/RESULTS.md` §1, `experiments/344-null-key-killcheck/RESULTS.md`.

## How this differs from the prior record (mandatory positioning)
- vs 155 §4 (constant all-down 0.5974 unweighted vs incumbent 0.5366; uniform additive down-shift
  closed at 4:1 pds cost): my construction is a GENE-SPECIFIC common vector (per-gene mean/vote
  over other targets), not a uniform shift, and it is measured on the |lfc|-WEIGHTED gate
  (343's object), not the unweighted pool. A uniform shift cannot beat the weighted gate by
  construction if the shared signal is gene-specific; a gene-specific vector can.
- vs 246 §1 / 232/234 (sign-consistent consensus closed; CD4-residual shared subspace deflated
  FOR pds): opposite direction of use. 246 REMOVES the shared component to buy discrimination;
  I ADD a shared component to buy sign accuracy on the reference's own gate. Both can be true:
  shared helps signs, hurts retrieval. The tension is expected and will be stated, not resolved.
- vs 343 (0.588 delivered / 0.600 designed weighted; counterfactual table): 343 is the baseline
  and the known-answer. I reuse its gate construction exactly and ask how much of the missing
  0.40 a target-blind common vector recovers.

## Estimator
- Gate per target (exactly 343): reference `p_adj<.05` & finite nonzero lfc & gene ≠ target;
  keep the 58 targets with ≥10 gate genes. All weighted stats use weights |ref lfc|.
- R: 85×G matrix, R[t,g] = reference lfc if g in gate(t) else 0 (G = 10,779 universe).
- D: cap15 design mapped to universe, then D_sup[t,g] = D[t,g] if g in gate(t) else 0.
- SVD energy fractions ||R_k||_F²/||R||_F² for k ∈ {1,2,3,5,10}, same for D_sup; |cos| between
  first left/right singular vectors of R and D_sup (report signed too).
- LOO common vectors from R: (i) mean over other 84; (ii) sign-vote Σsign(R) over others;
  (iii) per-gene |lfc|-weighted mean c_g = Σ R|R| / Σ|R| over others (0/0 → 0). Weighted +
  unweighted sign accuracy of sign(c_{-t}) on gate(t); constant all-down on the same gate.
- Combination: sign(D[t] + α·s·c_{-t}), s = median|D|/median|c| over pooled gate pairs (so at
  α=1 the medians match), α ∈ {0,.25,.5,1,2,4,∞}. Weighted acc; implied NMAE ≈ 1−m(2p−1) at
  m=0.205 (labelled approximation); EXACT NMAE with delivered magnitudes and combination signs.
- Transfer: (a) design's own LOO mean over other 84 D rows, alone + in combination;
  (b) donor-bank mean response (K562 GWPS; Orion HCT116 if same layout), gene-mapped, alone +
  in combination. Donor deltas are pseudobulk log-scale effects; only their SIGNS are read.
- Best combination: ref-UP vs ref-DOWN split + basal-CPM quintiles (343 bins343.py shape).

## Controls / falsifiers
- Known-answer first: reproduce 343's 0.588 delivered / 0.600 designed weighted to 1e-6; reproduce
  323's mean NMAE 0.950836 from saved stats. If either fails, stop and report.
- Target-shuffle control: permute the 84 donor rows of the common vector (deranged target
  assignment, 200 seeds) — if the shuffled common vector scores the same, the gain is a marginal
  effect, not target-shared structure (155's lesson).
- Sign-ambiguity guard on singular-vector cosines: report |cos|.
- Zero-handling stated: sign(0) counts as wrong (343 convention); report nonzero coverage of
  each common vector on the gate (mass share with c≠0).

## Kill criterion
The idea is DEAD as a design lever if the best TRANSFERABLE combination (step 4, design's own
or donor-bank common vector + D) gains < +0.02 |lfc|-weighted sign accuracy over the design
alone (0.600 baseline). The step-3 oracle ceiling (same-line LOO common vector) is reported
separately as what a same-line prior could give — it does not clear the idea by itself.

## Output
`RESULTS.md` (protocol sections), `work/results.json` (every number), scripts in the exp dir.
No emission, no 323 scoring run, no builds. Closing line: "worth a 323 six-member run: yes/no,
with which construction". No recommendations beyond that.
