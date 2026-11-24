# PLAN.md — 377-l1000-response-subspace-prior

## Hypothesis
A recipient line's own perturbation-response subspace — estimated from LINCS L1000
signatures in that line — is a useful prior on cross-line CRISPRi direction.
Concretely: projecting (or shrinking) a donor-consensus prediction onto the
recipient line's low-rank L1000 response subspace raises |response|-weighted sign
accuracy against the line's held CRISPRi truth, beyond what any same-dimension
subspace achieves. If true, a HeLa subspace artifact could feed a later builder
mode for context B (HeLa). If false, the "response-structure prior" route is dead.

## Bank gap found during planning (deviation from the mandate, recorded before running)
The mandate assumes held level-5 **compound** (`trt_cp`) values ("26,170 HELA
signatures … in `data/external/cmap/`"). What is actually on disk (verified
2026-09-10 with `uv run python`, h5py):
- `level5_beta_trt_sh.gctx`: 238,351 rows × 12,328 genes, COL ids `CGS001_<CELL>_96H:<GENE>:<n>`
  (shRNA consensus signatures). Cell mix: VCAP 32k, PC3 31k, MCF7 29k, A375/A549/HT29/HA1E
  ~21–23k, **HEPG2 18,727**, HCC515 18k, …, **JURKAT 256, HEK293T 97, K562 0, HELA 0,
  HCT116 0**.
- `level5_beta_trt_xpr.gctx`: 142,901 rows, mostly `.311` OE lines; K562 6, HELA 0, HEPG2 0.
- `siginfo_beta.txt` (1.2M rows) and `genetic_sigs.parquet` (177k: trt_xpr + trt_sh.cgs)
  carry metadata only; the `trt_cp` **values** were never downloaded (a cp beta file at
  720k × 12,328 would be ~35 GB).
- Consequence: the as-written K562-compound and HeLa-compound legs are **not executable**
  without a >2 GB download (violates the ≤10 GiB / 4 h budget). The executable analogue
  is the **same-line shRNA-response subspace** — same hypothesis object (a same-line
  perturbation-response subspace as a direction prior), and shRNA knockdown is, if
  anything, *closer* to CRISPRi than compounds are. Inference rule: if the same-line
  shRNA subspace cannot beat wrong-line/shuffled/random subspaces, the farther-modality
  compound variant is not rescued by this experiment (compounds would need their own
  experiment with the downloaded bank, out of scope here).
- Test lines actually runnable: **HEPG2 primary** (18,727 same-line sh rows; truth
  `experiments/166-donor-acquisition/donors/nadig_hepg2/`, 2,160 perts), **JURKAT
  secondary** (256 rows — weak, reported with that caveat; truth
  `data/donors/nadig_jurkat/`, 2,270 perts). K562/HELA/HCT116/HEK293T/RPE1 have no
  usable same-line L1000 response bank on disk → not tested. The mandate's
  "K562 and HepG2 both" gate for producing `work/hela_l1000_subspace.npz` therefore
  **cannot fire** (K562 untestable, HELA bank absent); no artifact will be produced.

## Exact inputs (all read-only)
- L1000: `data/external/cmap/level5_beta_trt_sh.gctx` (row-chunked reads, never full),
  `data/external/cmap/geneinfo_beta.txt` (landmark mask: `feature_space==landmark`, 978).
- CRISPRi donors: `data/donors/{k562_gwps,orion_hct116,orion_hek293t,rpe1,nadig_jurkat}/`
  + `experiments/166-donor-acquisition/donors/nadig_hepg2/` (`delta.npy` mmap,
  `genes.txt`, `perts.txt`, `meta.json`).
- Known-answer: `experiments/323-current-family-validation/work/{metrics-cap15.json,
  per-target-cap15.csv, stats-cap15.npz, reference_de.parquet}`.
- Own code only in `experiments/377-l1000-response-subspace-prior/` (`plan` = this file,
  `run377.py`, `work/`).

## Estimator
1. Gene axis G(L) = landmark symbols ∩ L1000 ROW ∩ L-axis ∩ every consensus-donor axis
   (~600–800 genes expected). All vectors on G(L).
2. Base prediction per target t (donors = all CRISPRi lines except L, carrying t; require
   ≥2 donors): unit-RMS donor directions (`rms_direction` = divide by RMS, core334
   verbatim) → (i) `source_mean` (mean of unit vectors, renormalised to unit RMS);
   (ii) `diagonal` (334-style transfer gain: per-gene gain fit pooled over donor lines
   acting as pseudo-recipients, applied to the source mean). Truth = L's own unit-RMS
   response. Targets = all L-bank targets with ≥2-donor coverage.
3. Subspace: L1000 same-line sh matrix S_L (rows = signatures, cols = G(L)), gene-centred
   → SVD → U_k, k ∈ {10, 20, 50, 100}. Report variance explained at each k.
   Corrections on unit-RMS p with centred projection (no mean add-back, to avoid
   injecting the target-independent drug blob): (a) `p' = U Uᵀ(p−μ)`;
   (b) shrinkage `p' = p + γ(UUᵀ(p−μ) − (p−μ))`, γ ∈ {0.25, 0.5} (keeps the common mean
   out of both terms); (c) sign-only hybrid (projected sign where |projected| above its
   per-target median, else original sign).
4. Controls (identical pipeline): (i) wrong-line subspace (MCF7 sh, largest bank, plus an
   n-matched MCF7 subsample); (ii) same-line gene-shuffled subspace; (iii) random
   orthonormal subspace, same k (3 seeds); (iv, if time) held-target subspace (same-line
   rows whose knocked-down gene ∉ evaluated CRISPRi targets).
5. Endpoints per target: |truth|-weighted sign accuracy (**primary**), unweighted sign
   accuracy, own cosine, cosine-rank retrieval across targets (334 `endpoints` formula,
   the 346 pds proxy). Paired 10,000-target bootstrap of Δ vs the uncorrected base
   (target resamples; report mean + 95% CI), per line; pooled = inverse-variance over
   lines (in practice HEPG2-dominated; JURKAT reported separately).

## Falsifiers / kill criterion (from the mandate)
DEAD if any holds: (1) best configuration's weighted sign gain over uncorrected < +0.01
pooled; (2) wrong-line subspace matches it (within the same bootstrap CI); (3) retrieval
falls by more than the sign gain rises (the 155/362/363 exchange). ALSO DEAD (bank):
K562/HELA legs untestable → no HeLa artifact regardless of HEPG2 outcome.

## Budget & hygiene
≤4 h wall; ≤12 GiB RSS (expect <3 GiB); `NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`;
`.h5ad` untouched (donor `.npy` via mmap); L1000 via row-chunked h5py reads
(≤~100 MB blocks); delete chunk scratch; no writes outside this dir; no commits,
no builders, no submissions.
