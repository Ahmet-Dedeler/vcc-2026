# 379 — PLAN: does recipient-similarity donor weighting beat uniform on the six-line frame?

## Check the record first

Grepped `experiments/*/RESULTS.md` and `docs/*.md` for "similarity-weighted", "per-context donor",
"donor weight", "lineage", "nearest donor", and read 041, 099, 175, 195, 204, 224, 246 (its "H5" is
per-*source* dose for CD4/H1, not donor-similarity weighting — no literal "H5" label exists in the
file; the six named hypotheses are enumerated in `RESULTS.md` §0), 263, 271, 289, 308, 312, 313, 337,
353, 360, plus scanned every experiment 325-378 by directory name for a closer match.

**Nothing has run this exact construction, but three prior results bound expectations hard and must
be cited honestly:**

1. **224 §3, rung R4b "basal-similarity donor weights (010's)"** is the closest analogue: a single
   scalar turning basal-correlation similarity into donor weights, on the S7 seven-donor-bank surface
   (968 shared targets, 5,070 shared genes, macro **cosine** transfer as the endpoint, not sign
   accuracy). Result: **relative cosine +2.81% [+2.28, +3.36], Δpds −0.0051, below the +2.9% gate,
   P(>gate) = 0.376**. Positive but sub-gate, and paid for in `pds`.
2. **224 §8, rung R4a "recipient-basal gain (184's shape)"** — a *learned*, recipient-conditioned
   per-gene gain (12 params) — clears only +0.41% and **the wrong-recipient null beats the honest arm
   outright** (0.17634 vs 0.17535, bit-identical in 30/84 configs). This is 184's Tahoe finding
   (recipient-conditioned parameters read `+0.0665` for the wrong recipient) reproduced on CRISPRi.
3. **353 (2026-09-10, same night, predates the 360 DepMap bank)** is a literature+repo scout, not an
   executed arm. It explicitly flags: "Nobody with ≤6 lines has made a *learned* recipient encoder
   work" (TxPert's own ablation; our 332). It proposes three DIFFERENT constructions (expression-gated
   network propagation, supervised per-gene sign classification, genotype-stratified priors) and does
   not propose or test fixed-form similarity-weighted donor averaging.

**What is actually new here, and why it is not a re-run of 224 R4b:** (a) the six-line/1,024-target/
512-gene frame from `experiments/334-text-conditioned-transfer/selection.json` plus `core341.py`'s
matched estimator, not S7; (b) the primary endpoint is |lfc|-weighted **sign accuracy** and
**retrieval** (343's decisive channel — cosine gain has repeatedly failed to imply a score gain, e.g.
195, 224 itself), not cosine; (c) 360's DepMap/CCLE covariate bank is used as an explicit similarity
source, which did not exist when 224 or 353 ran; (d) a softmax temperature sweep and an explicit
wrong-recipient falsifier on the *weight vector itself* (224's wrong-recipient test was only run on
the recipient-*conditioned* rungs R4a/R5, never on R4b). **Not closed. Proceeding, with 224 R4b's
result as the working prior: expect a small, probably sub-gate effect, and take the wrong-recipient
control especially seriously given R4a's outcome on a mechanistically adjacent construction.**

## Hypothesis

Weighting the five donor lines per held-out recipient by a public-covariate similarity score (basal
expression correlation, dependency correlation, or the banks' own control-pseudobulk correlation)
before forming the consensus, then re-fitting 341's per-gene cross-line diagonal gain on top, improves
pooled |lfc|-weighted sign accuracy over the uniform-weight incumbent without a retrieval collapse.

## Exact inputs

- `experiments/334-text-conditioned-transfer/selection.json`: 6 lines (`k562_gwps`, `orion_hct116`,
  `orion_hek293t`, `rpe1`, `nadig_jurkat`, `nadig_hepg2` = K562, HCT116, HEK293T, RPE1, Jurkat, HepG2),
  1,024 targets, 512 genes.
- `experiments/334-text-conditioned-transfer/run334.py`'s `load()` (bank `delta.npy` + `b_ctrl.npy`,
  indexed to the frozen targets/genes) and `experiments/334-text-conditioned-transfer/core334.py`'s
  `rms_direction`, `order`, `source_blocks`-style LOO structure, `endpoints`.
  `experiments/341-known-line-calibration/core341.py`'s `sign_endpoints`, `interval`.
  All imported read-only via `sys.path` insertion (342's pattern); nothing written there.
- `experiments/360-depmap-three-lines/work/three_lines.npz`: `genes` (18,533 challenge axis),
  `expr_ids`/`expr_logtpm1p` (DepMap 24Q4 log(TPM+1)), `dep_ids`/`dep_chronos` (Chronos gene effect).
  DepMap IDs for our six lines (checked directly against `Model.csv` and the subset CSVs in
  `experiments/360-depmap-three-lines/work/`, this run, read-only): **K562 ACH-000551, HCT116
  ACH-000971, HepG2 ACH-000739, Jurkat ACH-000995** all present in both expression and Chronos
  subsets (510/512 and 502/512 of the 512-gene axis by name, respectively). **HEK293T is not a
  DepMap model at all** (0 hits in `Model.csv`; matches the task brief's warning). **RPE1 is also
  absent** from both subsets — `Model.csv` carries only engineered `RPE1SS*` clones, no parental
  RPE1 row (360 §4b flagged this for dependency; confirmed here it holds for expression too, which
  the task brief did not flag). So CCLE/Chronos similarity has **4/6 lines instrumented
  (K562/HCT116/HepG2/Jurkat)**; HEK293T and RPE1 fall back to the banks' own control-pseudobulk
  correlation for any pair touching them, exactly as the brief prescribes for HEK293T.

## Estimator

1. **Similarity, three sources, one 6x6 matrix each** (`i,j` over the six lines, `i != j`):
   - `S_bank`: Pearson correlation of `basal[i]` vs `basal[j]` (each bank's own `b_ctrl.npy`, already
     restricted to the 512-gene axis by `load()`) — **always defined, all 15 pairs**.
   - `S_ccle`: Pearson correlation of `three_lines.npz` `expr_logtpm1p` rows over the intersection of
     the 512-gene axis with the array's `genes` field, for the 4 DepMap-covered lines; **hybrid**:
     any pair touching HEK293T or RPE1 is filled from `S_bank` (the brief's own prescription,
     generalised to both uninstrumented lines).
   - `S_chronos`: same construction on `dep_chronos`, same hybrid fallback.
2. **Weights**: for recipient line `L` (whichever line is currently held out — either the true outer
   fold or an inner pseudo-recipient during gain-fitting) and its available donor set `D_avail`,
   `w_D ∝ softmax(S[L,D]/tau)` for `tau in {0.05, 0.1, 0.2}`; `tau -> inf` is uniform (**this is
   already the incumbent / the known-answer arm below, not a new run**); `top2` zeroes every donor
   but the two highest-`S[L,:]` and splits weight 50/50.
3. **Consensus and gain, mirroring 341 exactly except for the averaging step**: for the true outer
   line, `wbase = sum_D w_D * direction[D, test_targets]`. For the per-gene diagonal gain, redo
   341/334's inner leave-one-line-out *using weighted, not uniform, consensus*: for each of the outer
   fold's 5 donors acting in turn as an inner pseudo-recipient `R`, weight the remaining 4 by
   `softmax(S[R, donor]/tau)` to build `x_R` (train targets only), pair with `y_R = direction[R,
   train_targets]`, and fit `gain[g] = (sum_R x_R[.,g]*y_R[.,g] + prior) / (sum_R x_R[.,g]^2 + prior)`,
   `prior = median(denominator)` — this is core334.diagonal's exact formula, applied to weighted
   inputs instead of the uniform mean. Report both **consensus-only** (no gain) and
   **consensus+diagonal** (comparable to 341's "source mean" / "cross-line diagonal" pair); the
   diagonal-calibrated version is primary since that is what the known answer is stated on.
4. **Controls**: `uniform` (tau=inf, must equal the known answer); **random** — 20 seeded Dirichlet(1)
   draws of weights over each line's donor set, held fixed across the inner-fit and outer-query steps
   for internal consistency, same pooling; **wrong-recipient** — one fixed derangement `pi` of the six
   lines (`pi(L) != L`), weights for recipient `L` computed from `S[pi(L), :]` instead of `S[L, :]`,
   run for all three similarity sources at all four tau (incl. top2).
5. **Endpoints**: `core341.sign_endpoints` (|lfc|-weighted sign accuracy, i.e. the `weighted_sign` key)
   and `core334.endpoints` (`cosine`, `retrieval`) computed per target per held-out line (2 blocks
   covering all 1,024 targets once each), pooled (macro mean over the 6x1024 rows) and via
   `core341.interval` on the paired per-target-per-line delta array against the uniform incumbent
   (`target_ci95`, `crossed_ci95`).

## Known answer (reproduced before any new number, own driver, fresh data load)

341's pooled table: **source mean weighted-sign .589398 / retrieval .790095; cross-line diagonal
weighted-sign .589398 / retrieval .794614** (identical sign accuracy is expected: `core334.diagonal`
returns a strictly positive per-gene gain here, which cannot flip any sign). My driver reloads the six
banks from scratch via `run334.load()` (not by reading 341's cached `.npz` predictions), rebuilds
`direction`/`rms` via `core334.rms_direction`, replays 341/334's exact two-block LOO with uniform
weights, and must match both numbers before any similarity-weighted arm is interpreted. If it does not
match to float tolerance, STOP and report the mismatch instead of proceeding.

## Controls / falsifiers

- Uniform-weight arm must equal the known answer (structural check, not a result).
- Random weights (20 draws) bound how much of any gain is just "a non-uniform weight vector" rather
  than the specific similarity ordering.
- Wrong-recipient weights must NOT match or beat the honest arm; if they do, the effect is not about
  which line is held out (224 R4a's exact failure mode on the adjacent construction).
- `top2` vs the softmax family separates "pick the single best analog" from "graded blending".

## Kill criterion (from the task brief, verbatim)

**Dead** if pooled weighted-sign gain over uniform < +0.005, OR retrieval falls by more than the sign
gain, OR wrong-recipient weights match (perform comparably to) the honest arm.
**Alive** if >= +0.01 sign gain with retrieval not falling: then report Jurkat's and the epithelial
lines' (HeLa/CAL33 are not in this frame; the epithelial donor lines here are HCT116, HEK293T, HepG2)
best weights, since A/B/C could each get their own donor mix through the existing per-context
`context_effects` path.

## Scope and limits, stated in advance

This is the six-CRISPRi-donor-line sandbox (334/341's frame), not A/B/C. A sign-accuracy or retrieval
gain here is development evidence only, per programme rule 5 — it licenses a *board arm proposal*, not
a score claim. No `vcc`, no build, no submission; ~10 MB of in-memory arrays, nowhere near the 10 GiB
cap for this task.
