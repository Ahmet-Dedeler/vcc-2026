# PLAN.md — 359-expression-gated-propagation

## Hypothesis

A knockdown's effect propagates through regulatory/interaction edges, but an
edge carries response only where its nodes are expressed in the recipient.
Donor lines determine how strongly each edge propagates (prior weights +
donor-fitted per-gene gains); the recipient contributes only an expression
GATE over the same edges — label-free, from its control cells, zero fitted
recipient parameters. Prediction for gene g under target t in line L =
propagation of the donor-consensus response from t through the gated graph.

This is worker 353's proposal 1
(`experiments/353-recipient-conditioning-scout/RESULTS.md` §Exp.1).

## How this differs from the closed lookalikes (must state)

- **108** (signed-regulon transfer, NO-GO): regulon-only prediction /
  sign-preserving reweight with NO recipient gate and no donor-consensus base.
  Ours smooths the donor-consensus response and gates edges by recipient
  expression, with a wrong-recipient control 108 never ran.
- **190** (regulon ordering key, M2 closed): CollecTRI used as a confidence
  key re-ranking `|H|`, no propagation of response values, no gate. Ours moves
  predicted SIGNS via propagation, not ranks via a key; endpoint is weighted
  sign accuracy, not `reach`.
- **340** (neighborhood transfer, DOES_NOT_EARN): fitted response-gene ridge
  couplings (per-gene β over neighbor source values, fitted on donors),
  undirected, ungated — smoothing that bought cosine (+.0057) but lost
  retrieval (−.0052 vs diagonal, all 6 lines). Ours fits NO couplings
  (prior weights only, row-normalized) and adds the recipient gate plus
  wrong-recipient / random-gate / ungated controls 340 never ran.
- **339** (gene-wise marginal maps, +0.0029): no network structure at all.
  Ours is network-structured by construction.

## Frame and inputs (exact paths)

- Selection: `experiments/334-text-conditioned-transfer/selection.json`
  (6 lines, 1,024 targets, 512 coords; targets ∩ genes = ∅, verified).
- Loaders/primitives: `experiments/334-text-conditioned-transfer/run334.py`
  (`BANKS`, `load`, `known`) as `B`;
  `experiments/334-text-conditioned-transfer/core334.py` as `P`
  (`rms_direction`, `order`, `source_blocks`, `diagonal`, `endpoints`);
  `experiments/341-known-line-calibration/core341.py` as `C`
  (`sign_endpoints`, `interval`). Read-only imports; nothing written outside
  `experiments/359-expression-gated-propagation/`.
- Donor banks (read via `B.load`, mmap + copy of the 1024×512 window only):
  `data/donors/{k562_gwps,orion_hct116,orion_hek293t,rpe1,nadig_jurkat}` +
  `experiments/166-donor-acquisition/donors/nadig_hepg2`
  (`delta.npy`, `perts.txt`, `genes.txt`, `b_ctrl.npy`, `meta.json`).
  `b_ctrl` units: log1p(5e4·P/ΣP) (`scripts/build_donor_table.py:150`,
  `scripts/build_orion_donor.py:112`) — same units on all lines; Orion's low
  full-axis median is a zeros artifact of its 38,606-gene axis, so thresholds
  are computed over the 512-gene frame V, per line (medians reported).
- Graph priors:
  - STRING (primary, dense): `experiments/340-neighborhood-transfer/selection.json`
    `graph_neighbors` (512×8, positive-cosine KNN from
    `experiments/100-sign-confidence/work/string-pca.npz` unit vectors —
    same topology 340 used). Binary, zero diagonal, used as-is (directed);
    symmetrized variant as sensitivity.
  - CollecTRI (sparse check): `data/external/regulon/collectri_rescued.csv`
    (43,536 edges; sign = int(weight): +1 ×37,662 / −1 ×5,874; `sign.decision`
    is provenance, ignored). Measured coverage on this frame: **50 G×G edges**
    (13 TFs, 35/512 genes with an in-edge), **56 T×G edges covering 19/1024
    targets**. Consequence pre-registered: the CollecTRI smoothing arm can move
    ≤35/512 genes, so gated≈ungated≈diagonal there by arithmetic; the gate is
    tested on STRING, CollecTRI is a mechanism check + target-seeded
    exploratory (19 covered targets only, descriptive).

## Construction

Node set V = 512 response coords (smoothing arms). Base signal per
(outer line L, test target t): `s0(t,·) = baseline(t,·) × gain`,
`baseline` = mean `direction` over the 5 donors, `gain` = `P.diagonal`
fitted on the 5 donors' TRAIN targets (341's split regime — see below).
Truth: `direction[L]` (unit-RMS). All in 341's complementary target halves
(`P.order(targets,'334-split|')`, blocks [::2]/[1::2]), stitched to full
1024 for endpoints — bit-comparable with 341's diagonal.

Propagation operator: `P_L = D_L^{-1}(A ∘ G_L)`, `A` = prior adjacency
(512×512), `G_L` = recipient gate, `D_L` = row sums. Isolated rows (sum 0)
→ self-loop (P_ii = 1). One step: `s1 = s0 @ P_Lᵀ`; two steps:
`s2 = s1 @ P_Lᵀ`. Mixed: `s(λ) = (1−λ)·s0 + λ·s_k`, λ ∈ {0.25, 0.5, 1}.

- `A_string`: binary KNN (340 topology), unsigned (+1). Smoothing.
- `A_ct`: signed directed, `W[i,j] = +mor(j→i)` for j→i CollecTRI edge
  (co-response sign: an activator's target co-moves with it — +mor, NOT
  190's knockdown-direction −mor, which applies to seeds, not smoothing;
  cf. 190 §1's flipped-convention lesson). 1 step only.
- Exploratory `seeded-ct` (NOT under the kill criterion): seed vector over
  full V′ = 1024 targets ∪ 512 genes with −1 at the knocked-down target,
  one step over TxG+GxG CollecTRI edges with knockdown-direction sign
  (−mor, per 190's N), gene-block readout scaled to s0's per-target RMS,
  mixed at fixed λ = 0.5. Reported on the 19 covered targets only.

Gates (recipient-only, from `b_ctrl` restricted to the 512 frame):

- `hard`: node gate = 1[b ≥ τ_L], τ_L = line L's median over the 512.
  Edge gate = outer product (both endpoints). Density ≈ 25%, matched
  across lines by construction.
- `soft`: node gate = sigmoid((b − τ_L)/s), s = 0.5 (log1p units; per-line
  IQR measured and reported; if IQR differs wildly across lines, noted as
  a limit, not re-tuned).
- `random`: permuted node gates within the line (exactly same density,
  structure destroyed); 3 seeds, mean + range.
- `wrong`: hard gate from each of the other 5 lines' basal (own τ); report
  the BEST of the five (strongest control).
- `ungated`: all-ones (pure smoothing — isolates the gate).

λ choice: inner LOLO over the 5 donors on TRAIN targets only (for each
inner-held d: gain from the other 4, base, propagate, score weighted sign
on d's train truths); λ maximizing the inner mean, per outer line.
Primary arm: hard-gate 1-step at donor-chosen λ (fixed before seeing L).
Full fixed-λ sweep reported for transparency; 2-step and soft are
secondary (donor-chosen the same way, reported alongside, cannot rescue
the primary verdict).

No per-edge/per-pathway ridge is fitted (deviation from 353's sketch,
reason recorded): 5 donor lines cannot identify 4k–50k edge strengths —
that is the TxPert/332 failure mode (≤6-line learned recipient/donor maps
fail); the donor contribution is the per-gene diagonal gain + donor-chosen
λ. Recipient contributes zero fitted parameters.

## Endpoints (as in 358)

- Primary: |response|-weighted sign accuracy (`C.sign_endpoints`), pooled
  over 6×1024 + per-line means.
- Secondary: top-64 weighted sign; retrieval (`P.endpoints`) as GUARD only.
- Uncertainty: paired target bootstrap + crossed line/target intervals via
  `C.interval` (SEED fixed in core341).

## Controls / falsifiers

(a) ungated smoothing — is the gate doing anything;
(b) wrong-recipient gate, best of five — recipient specificity;
(c) random (shuffled-node) gate, same density — structure vs density;
(d) diagonal alone (λ = 0) — any value at all.
Mechanism check: CollecTRI gated vs ungated must be ≈ equal (≤35 genes
movable); if STRING gated fails while ungated smoothing helps, the GATE is
dead but smoothing is alive — reported separately, not smuggled.

## Kill criterion

DEAD if (pooled weighted sign): own-hard-gated − best-wrong-gated ≤ +0.005,
or gated ≤ ungated, or gated ≤ diagonal. ALIVE if ≥ +0.015 over diagonal
with the wrong-recipient control resolved below. (Scale: 341's labeled
+0.033; 353's falsifier: +0.005 is half a meaningful nmae step.)

## Known-answer gate (must pass before any arm is read)

End-to-end replay of 341's cross-line diagonal in the split regime with own
code: pooled weighted sign **.589398** and pooled retrieval **.794614**
(341 RESULTS.md table). Required match: ≤1e-12 (same deterministic code
path). If it fails: stop, report, no arms.

## Resources and limits

No 323 run. No downloads (all inputs on disk). No reserve data. Budgets:
≤12 GiB RSS (expected <2 GiB: direction 6×1024×512 f64 = 201 MB, rest
small; Orion delta.npy 2.7 GB touched via mmap window only),
`NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`, ~5 h. Memory-gate check recorded
before the main run. H1/323-panel numbers are not produced; six-line
source-frame sign/retrieval deltas are development evidence, not A/B/C
forecasts; no retrieval/cosine gain will be called a score gain.
Sibling 358 runs the supervised sign classifier on the same frame — no
shared code, no duplicated arm (this experiment fits nothing per-(t,g)).
