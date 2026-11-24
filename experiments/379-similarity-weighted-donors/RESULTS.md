# 379 — Recipient-similarity donor weighting: real but small signal, an order of magnitude below the alive bar on every source and every temperature

**Status: NEGATIVE.** Best honest arm (`ccle_tau0.1`) reaches pooled weighted-sign gain **+0.00402**
over the uniform incumbent — below the pre-registered **+0.005 dead threshold**, let alone the
**+0.01 alive threshold** — and its own retrieval falls by **more** than that gain (−0.00566 vs
+0.00402). No arm, at any of 3 similarity sources × 4 weight shapes (τ = 0.05/0.1/0.2, top-2), clears
either kill-criterion leg. The construction is not spurious — a fixed wrong-recipient derangement and
20 draws of random reweighting are both substantially *worse* than uniform, so the small honest gain is
real, recipient-specific information — but it is roughly an order of magnitude too small to matter, and
top-2 hard donor selection is actively harmful everywhere it was tried.

## Check the record first (mandatory step, done before any code)

Grepped `experiments/*/RESULTS.md` and `docs/*.md` for "similarity-weighted", "per-context donor",
"donor weight", "lineage", "nearest donor"; read 041, 099, 175, 195, 204, 224, 246, 263, 271, 289, 308,
312, 313, 337, 353, 360 in full or in relevant part, and scanned every experiment directory 325–378 by
name. **246's "H5"** does not exist as a literal label — the file enumerates six unlabelled hypotheses
in its §0; the nearest is "hypothesis 5" (§3.6), which is a per-*source* dose knob for CD4/H1, not
donor-similarity weighting, and is unrelated here.

**Nothing in the record runs this exact construction, but three results bound expectations hard:**

1. **224 §3, rung R4b "basal-similarity donor weights (010's)"** — the closest analogue: a single
   scalar turning basal-correlation similarity into donor weights, measured on the *S7* seven-donor-bank
   surface (968 shared targets, 5,070 shared genes) with macro **cosine** transfer as the endpoint.
   Result: relative cosine **+2.81% [+2.28, +3.36]**, Δpds **−0.0051**, below the **+2.9%** shape-only
   gate, P(>gate) = 0.376. Positive, sub-gate, paid for in `pds`.
2. **224 §8, rung R4a "recipient-basal gain (184's shape)"** — a *learned*, recipient-conditioned
   per-gene gain — clears only +0.41% relative cosine, and its **wrong-recipient null beats the honest
   arm outright** (0.17634 vs 0.17535, bit-identical in 30/84 configs). This reproduces 184's Tahoe
   finding on CRISPRi: recipient-conditioned *learned* parameters carry nothing.
3. **353 (same night, predates the 360 DepMap bank this task uses)** is a literature+repo scout, not an
   executed arm. It states "Nobody with ≤6 lines has made a *learned* recipient encoder work" (TxPert's
   own ablation; our own 332) and proposes three different, unexecuted constructions (expression-gated
   network propagation, supervised per-gene sign classification, genotype-stratified priors) — none of
   which is fixed-form similarity-weighted donor averaging.

**What is new here and why this is not a re-run:** the six-line/1,024-target/512-gene frame
(`experiments/334-text-conditioned-transfer/selection.json`) and `core341.py`'s matched estimator
rather than S7; the primary endpoint is |lfc|-weighted **sign accuracy** and **retrieval** (343's
decisive channel) rather than cosine; 360's DepMap/CCLE bank as an explicit similarity source (it did
not exist when 224 or 353 ran); a softmax-temperature sweep and a wrong-recipient falsifier applied
directly to the *weight vector*, which 224 only ran on the recipient-*conditioned* rungs R4a/R5, never
on the plain-weighting R4b. **Verdict: not closed. Proceeded** — full record check is in `PLAN.md`.

## Known-answer reproduction (rule 3)

Reloaded the six donor banks from scratch via `run334.load()` (fresh data load, not 341's cached
`.npz` predictions), rebuilt `direction`/`rms` via `core334.rms_direction`, and replayed 334/341's
exact two-block leave-one-line-out with **uniform** weights:

| quantity | published (341 RESULTS.md) | reproduced here | \|diff\| |
|---|---:|---:|---:|
| source-mean weighted sign | 0.589398 | 0.5893977594 | 2.4e-07 |
| source-mean retrieval | 0.790095 | 0.7900950661 | 6.6e-08 |
| cross-line-diagonal weighted sign | 0.589398 | 0.5893977594 | 2.4e-07 |
| cross-line-diagonal retrieval | 0.794614 | 0.7946138553 | 1.4e-07 |

All four match to float-precision (the published table is rounded to 6 dp; every reproduced value
rounds to the same digits). **Match confirmed — proceeded to the new measurement.**

Additional internal-consistency check (not in the record, done here as extra rigor before trusting my
own reimplementation of the diagonal-gain formula on weighted inputs): my hand-rolled
`diagonal_gain(xs, ys)` reproduces `core334.diagonal(core334.source_blocks(...))`'s output on the
uniform-weight case **exactly, max\|diff\| = 0.0e+00**, on the first fold. Weighted-consensus runs below
reuse the identical formula.

## What was measured

**Inputs**: `experiments/334-text-conditioned-transfer/selection.json` (6 lines = K562/HCT116/HEK293T/
RPE1/Jurkat/HepG2, 1,024 targets, 512 genes) and its `run334.load()`/`core334.py` (imported read-only,
`sys.path` insertion, nothing written there); `experiments/341-known-line-calibration/core341.py`
(`sign_endpoints`, `interval`, read-only); `experiments/360-depmap-three-lines/work/three_lines.npz`
(DepMap 24Q4 `expr_logtpm1p`, `dep_chronos`, both on the 18,533-gene challenge axis, read-only).

**DepMap ID mapping** (checked directly against `Model.csv` and the CSV subsets in
`experiments/360-depmap-three-lines/work/`, this run): K562 = ACH-000551, HCT116 = ACH-000971, HepG2 =
ACH-000739, Jurkat = ACH-000995 — all four present in both the expression and Chronos subsets, 510/512
and 502/512 of the 512-gene axis by name respectively (498/512 finite after intersecting with the
axis's own missing-gene pattern). **HEK293T is not a DepMap model** (0 hits in `Model.csv`, as the task
brief warned). **RPE1 is also absent** — `Model.csv` carries only engineered `RPE1SS*` clones, no
parental row; the brief did not flag this one, but it is real (confirmed here, independently of 360's
narrower note about RPE1's *dependency* row). So CCLE/Chronos similarity is instrumented for **4 of 6
lines** (K562/HCT116/HepG2/Jurkat); any pair touching HEK293T or RPE1 falls back to the banks' own
`b_ctrl.npy` control-pseudobulk correlation, exactly as the brief prescribes for HEK293T, extended here
to RPE1 by necessity.

**Estimator.** For held-out recipient line *L* (either the true outer fold, or an inner pseudo-recipient
during gain-fitting — both use the identical machinery), and its available donor set: `w_D ∝
softmax(S[L,D]/τ)` for τ ∈ {0.05, 0.1, 0.2} (τ→∞ is uniform, already the known answer above); `top2`
zeroes every donor but the two highest-similarity ones and splits 50/50. The weighted consensus
`wbase = Σ_D w_D · direction[D]` replaces the uniform mean at **both** the outer query and every inner
leave-one-line-out fit used to calibrate 341's per-gene diagonal gain (`gain[g] = (Σ x·y + prior)/(Σ x²
+ prior)`, `prior = median(Σx²)`), so the two arms being compared differ *only* in how the consensus is
formed, not in the downstream calibration. Three similarity sources: `bank` (control-pseudobulk
correlation, always defined, all 15 pairs), `ccle` (DepMap basal log-TPM correlation, hybrid-filled from
`bank` for HEK293T/RPE1), `chronos` (DepMap Chronos dependency correlation, same hybrid fill).

**Controls.** `random`: 20 seeded Dirichlet(1) draws of donor weights, redrawn fresh at every weighting
call (outer query and every inner fit) from one continuous per-draw RNG stream — same pipeline shape as
the honest arms, no similarity information. `wrong-recipient`: one fixed derangement π of the six lines
(`k562_gwps→rpe1, orion_hct116→k562_gwps, orion_hek293t→nadig_jurkat, rpe1→orion_hek293t,
nadig_jurkat→nadig_hepg2, nadig_hepg2→orion_hct116`), weights for recipient *L* computed from `S[π(L),
:]` instead of `S[L, :]`, run at every τ and for `top2`, all three similarity sources.

**Endpoints**: `core341.sign_endpoints` (|lfc|-weighted sign accuracy) and `core334.endpoints`
(retrieval, cosine), per target per held-out line (1,024 × 6 = 6,144 rows, two complementary blocks
covering every target exactly once per line), pooled as the macro mean, and via `core341.interval` on
the paired per-target-per-line delta against the uniform incumbent (10,000-draw target-resample and
crossed line/target-resample 95% intervals).

## Results

**Similarity matrices** (off-diagonal range 0.71–0.86; computed, not assumed: Pearson-of-correlations
between `S_bank` and `S_ccle` over all 15 off-diagonal cells is **0.845**, but that number is inflated by
the 9 cells both share by hybrid construction — restricted to the **6 pairs where `S_ccle` is genuinely
DepMap-derived on both sides** (the K562/HCT116/HepG2/Jurkat submatrix), the agreement drops to
**0.574**. The two sources are only moderately redundant, not near-duplicates of each other):

| | k562 | hct116 | hek293t | rpe1 | jurkat | hepg2 |
|---|---:|---:|---:|---:|---:|---:|
| `S_bank` (own control pseudobulk, always available) | — | .781 | .762 | .794 | .852 | .827 |
| `S_ccle` (DepMap basal log-TPM, hybrid-filled for hek293t/rpe1) | — | .823* | .729* | .794 (=bank) | .797* | .783* |
| `S_chronos` (DepMap Chronos, hybrid-filled) | — | .823* | .729* | .794 (=bank) | .788* | .845* |

(* = genuinely DepMap-derived cell; unmarked cells for hek293t/rpe1 pairs fall back to `S_bank`.)

**Honest arms** (pooled over 6 lines × 1,024 targets; Δ is vs the uniform incumbent's `plus_diagonal`
row, weighted-sign .589398 / retrieval .794614 / cosine .198735):

| arm | weighted sign | Δ sign [target 95% CI] | retrieval | Δ retrieval [target 95% CI] | Δ cosine | sign-positive lines |
|---|---:|---|---:|---|---:|---:|
| `bank` τ=0.05 | .591223 | +.00183 [+.00076,+.00289] | .780846 | **−.01377** [−.01765,−.00989] | +.0047 | 4/6 |
| `bank` τ=0.1 | .593109 | +.00371 [+.00305,+.00437] | .790283 | −.00433 [−.00661,−.00208] | +.0089 | 6/6 |
| `bank` τ=0.2 | .592332 | +.00294 [+.00253,+.00334] | .793356 | −.00126 [−.00246,−5e-5] | +.0063 | 6/6 |
| `bank` top2 | .582052 | **−.00735** [−.00859,−.00610] | .771710 | **−.02290** [−.02775,−.01815] | −.0172 | 0/6 |
| `ccle` τ=0.05 | .591847 | +.00245 [+.00145,+.00345] | .779566 | −.01505 [−.01877,−.01139] | +.0064 | 4/6 |
| **`ccle` τ=0.1 (best)** | **.593418** | **+.00402** [+.00338,+.00467] | .788950 | −.00566 [−.00792,−.00343] | +.0096 | 6/6 |
| `ccle` τ=0.2 | .592527 | +.00313 [+.00274,+.00353] | .792629 | −.00199 [−.00320,−.00077] | +.0067 | 6/6 |
| `ccle` top2 | .581462 | **−.00794** [−.00916,−.00672] | .767678 | **−.02694** [−.03158,−.02218] | −.0172 | 0/6 |
| `chronos` τ=0.05 | .589528 | +.00013 [−.00072,+.00099] | .782806 | −.01181 [−.01490,−.00864] | −.0007 | 2/6 |
| `chronos` τ=0.1 | .591541 | +.00214 [+.00159,+.00270] | .790318 | −.00430 [−.00620,−.00235] | +.0042 | 3/6 |
| `chronos` τ=0.2 | .591195 | +.00180 [+.00145,+.00216] | .793172 | −.00144 [−.00250,−.00035] | +.0037 | 3/6 |
| `chronos` top2 | .575318 | **−.01408** [−.01538,−.01282] | .763494 | **−.03112** [−.03581,−.02634] | −.0352 | 0/6 |

**Controls:**

| control | weighted sign | Δ sign vs uniform | reading |
|---|---:|---:|---|
| `random` × 20 draws | mean .578993, range [.573716, .581959] | **−.01040 mean** | worse than uniform on every draw |
| `WRONG_bank` τ=0.1 | .566007 | −.02339 [target CI −.02467,−.02213] | far worse than uniform |
| `WRONG_ccle` τ=0.1 | .564180 | −.02522 [−.02657,−.02388] | far worse than uniform |
| `WRONG_chronos` τ=0.1 | .566614 | −.02278 [−.02403,−.02154] | far worse than uniform |

(all 12 wrong-recipient cells, all τ and top2, all three sources, are negative vs uniform; full table
in `work/results379.json`.)

**Bridge to 224 R4b, a different surface and a different metric.** The relative *cosine* gain here —
`bank`τ0.1 +4.48%, `ccle`τ0.1 +4.82%, `chronos`τ0.1 +2.13% — lands in the same small-positive range as
224's independently-measured **+2.81%** on the unrelated S7 seven-donor-bank surface. Two different
similarity sources, two different surfaces, two different pipelines, same order of magnitude: this is
mutual corroboration that the effect is real, and that it is small everywhere it has been looked for.

## Verdict against the PLAN.md kill criterion

**Dead** (as pre-registered) fires on the primary leg alone: pooled weighted-sign gain over uniform
**< +0.005** for every one of the 12 honest cells (max **+0.00402**, `ccle` τ=0.1, 80% of the dead
threshold and only 40% of the +0.01 alive bar). It also fires on the second leg for the
higher-sign-gain cells: at `ccle`τ0.1 and `bank`τ0.1 the retrieval fall (−.00566 / −.00433) exceeds the
sign gain (+.00402 / +.00371). Only the τ=0.2 cells keep retrieval fall smaller than the sign gain, and
those have the smallest sign gains of the softmax family. **The alive bar (≥ +0.01 sign, retrieval not
falling) is not approached by any cell.** Wrong-recipient weights do not match or beat the honest arm —
they are dramatically worse than even uniform — so this is a clean, unambiguous **DOES NOT CLEAR**, not
an ambiguous middle case.

**What is true and worth carrying forward, short of clearing the gate:**
1. The signal is real, not noise dressed up as a similarity score: the honest arms (+0.001 to +0.004
   sign) sit clearly above both the random-reweighting control (−0.010 mean) and the wrong-recipient
   control (−0.008 to −0.039 across the 12 wrong-recipient cells) — scrambling either the weight
   *values* or the weight *identity* destroys
   it and goes solidly negative. That rules out "any non-uniform reweighting looks like this."
2. It is nonetheless small: the best cell (+0.00402) sits *inside* the pre-registered dead zone
   (< +0.005), not near it, and would have to grow by **~2.5×** just to reach the +0.01 alive bar —
   there is no cell close to either threshold in a way a bigger sample or a different τ grid point
   would plausibly close.
3. **Hard top-2 selection is a clean NO-GO** everywhere: sign −0.007 to −0.014, retrieval −0.023 to
   −0.031, 0/6 lines positive, on all three similarity sources. Discarding graded information for a
   binary pick costs more than the entire softmax family gains — consistent with 263's finding that
   discrete donor-selection arbiters cannot capture what a continuous oracle can.
4. Per-line heterogeneity exists (`ccle`τ0.1: RPE1 +0.0088, HepG2 +0.0067, but HEK293T only +0.0014) —
   noted for the record since the task brief asks for it under "alive", but the pooled gate is the
   pre-registered decision variable and it did not clear, so this is not read as a per-context
   recommendation.

## Limits

Development evidence on the six-CRISPRi-donor-line sandbox only (334/341's frame), not an A/B/C
forecast (programme rule 5). A cosine or sign-accuracy gain here does not imply a board score gain —
this project's own record (195, 224) shows that repeatedly. RPE1 and HEK293T similarity is
bank-pseudobulk-only, not DepMap-informed, for 2 of 6 lines and 9 of 15 pairs; the `ccle`/`chronos`
arms are therefore a hybrid, not a clean 6-line DepMap test (a pure 4-line DepMap-only submatrix was
not separately re-run given the pooled result was already unambiguous). Chronos-based similarity is
weaker than the other two sources at every τ, consistent with 360 §4b's finding that dependency
correlation among these lines is comparatively flat and less discriminative than expression. The
wrong-recipient control uses one fixed derangement, not an average over all 5 wrong choices per line;
given the effect size (−0.008 to −0.039 across all 12 cells, dwarfing the honest arms' +0.001 to
+0.004) this would not
change the qualitative reading. Same-line training targets are available within the inner LOO fit
(341's regime, not unseen-context validation) — inherited from the frame, not introduced here.

## Exact commands to reproduce

```
memory_pressure | tail -1   # gate: proceed only if free >= 50% (measured 78% before this run)
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python3 experiments/379-similarity-weighted-donors/run379.py
```
Writes only `experiments/379-similarity-weighted-donors/work/known_answer.json` and
`.../work/results379.json`. `run379.py` reads `experiments/334-text-conditioned-transfer/`,
`experiments/341-known-line-calibration/`, and
`experiments/360-depmap-three-lines/work/three_lines.npz` (its hardcoded `LINE_ACH` name-to-ACH-ID map
was verified separately, ad hoc, against `experiments/360-depmap-three-lines/work/Model.csv` and the
CSV subsets before being written into the script — see the DepMap ID mapping note above; the script
itself does not re-parse `Model.csv`). Nothing outside `experiments/379-similarity-weighted-donors/`
was written. No `vcc`, no build, no submission.

## Peak memory and wall time

**18.0 s** wall (`main()` internal clock) / 18.81 s real (`/usr/bin/time -l`). **621 MB** maximum
resident set size — direction/basal/rms arrays are ~25 MB total (6×1,024×512 float64); the entire 26-arm
sweep (12 honest + 12 wrong-recipient + 20-draw random + known-answer + internal check) reruns the full
six-line, two-block LOO each time but stays small throughout. Memory gate: `memory_pressure` 78% free
before the run (gate: proceed if ≥ 50%); well under the 10 GiB task cap and the 12 GiB protocol cap.
