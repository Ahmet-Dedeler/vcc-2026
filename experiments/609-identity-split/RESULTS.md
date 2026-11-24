# 609 — where the identity rider's sign gain sits: the 639d (block) vs 622m (uncalled) conflict, on shared definitions

[4d0207] for [7772dc] (council 600), 2026-10-02 ~13:50Z. Analysis only. Peak RSS 3.1–3.2 GB per run, run directly (< 4 GB).

**Answer: the conflict is REAL in the data, but 639d's block result is an ORION ARTIFACT. On the clean surfaces (H1, and
the Broad-16 under production's sign rule) the gain sits mostly on UNCALLED genes. Centring zc (a definitional
difference) is minor and is not the cause.**

## What was run
- **Part B, `split609.py`.** 639's four folds with 639d's OWN masks (head 400 by the consensus key + the 2,200-gene CPM
  block) and 639d's identity (Flex references). **KA: the `rider` row reproduces `decomp639.json` exactly on every fold
  and region.** Four sources, each against its matching CLEAN4 baseline:
  - `rider`: zc uncentred, production v3's rule;
  - `rider_c`: 622m's mean-centred zc;
  - `rider_expr`: expression E = log1p(own control CPM) in place of identity;
  - `rider_IperpE`: identity with expression partialled out.
  The run also reports each region's gate-gene DOWN fraction.
- **Part A, `block609.py`.** 622m's Broad-16 script verbatim plus one switch. **KA: the `centred` run equals
  `622/work/block622_gated.json` exactly (all lines, regions, sources, CIs).** The `uncentred` run uses production's
  zc rule.

## Part B — 639's folds (λ 0.5; gain = sign-accuracy change on the region's gate genes; sh = share of the total net)

| fold | region | gate pairs | gate DOWN frac | rider | rider_c (centred) | **rider_expr** | **rider_IperpE** |
|---|---|---:|---:|---|---|---|---|
| **h1** (clean Arc Flex) | block | 21,640 | 0.56 | +0.0057 (sh 0.30) | +0.0076 | **−0.0078** | **+0.0080** |
| | uncalled | 39,984 | 0.51 | **+0.0132 (sh 0.69)** | +0.0118 | +0.0006 | +0.0129 |
| **hct116** (Orion) | block | 7,811 | **0.72** | +0.0399 (sh 0.99) | +0.0515 | **+0.0474** | +0.0245 |
| | uncalled | 3,048 | 0.28 | −0.0024 | −0.0088 | +0.0129 | −0.0241 |
| **hek293t** (Orion) | block | 24,077 | **0.74** | +0.0333 (sh 0.92) | +0.0497 | **+0.0842** | **−0.0061** |
| | uncalled | 12,834 | 0.25 | +0.0088 | +0.0040 | −0.0064 | +0.0079 |
| k562 | block | 1,824 | 0.51 | −0.0289 | −0.0233 | −0.0444 | −0.0100 |
| | uncalled | 828 | 0.22 | +0.0540 | +0.0423 | +0.1127 | −0.0023 |

- **On the two Orion folds the block's gate genes are 72–74 % DOWN** (h1 and k562: 0.51–0.56). This is the documented
  depth/technical skew of Orion truth: highly expressed genes are called down (AGENTS 6.5;
  `docs/HARNESS-DEPTH-CONFOUND.md`).
- On those folds, **expression alone** (`rider_expr`) earns as much block "gain" as identity, or far more: HCT116 +0.047
  vs +0.040; HEK293T +0.084 vs +0.033.
- **With expression partialled out, HEK293T's block gain disappears (−0.006).** On HCT116 it halves (+0.025).
- So 639d's "≥ 0.83 of the net on the block" on HCT116 / HEK293T is mostly the Orion truth's expression-linked DOWN
  skew, not line identity.
- **On h1, the clean surface validated against the board, the gain is mostly UNCALLED** (+0.0132, share 0.69). The
  block gain is smaller but real: +0.0057, share 0.30, and +0.0080 with expression partialled out, while expression
  alone loses (−0.0078) there.

## Part A — Broad-16 (12 lines, Cas9 KO, gated consensus head), centred (= 622m) vs uncentred (production)

| zc rule | block mean d (neg / 12) | head | **uncalled** | all | pooled net share block / head / uncalled |
|---|---|---|---|---|---|
| centred (622m; KA = block622) | −0.0022 (7) | +0.0133 | +0.0204 (0 neg) | +0.0081 | −0.11 / 0.28 / 0.83 |
| **uncentred (production)** | **−0.0007 (6)** | +0.0089 | **+0.0231 (0 neg; 8 CIs > 0)** | +0.0090 | **−0.09 / 0.17 / 0.92** |

On the production rule, `rider_expr` is block −0.0168 (10 of 12 lines negative) and `rider_I_perp_E` is block −0.0074.
On the Broad-16, expression in the block is harmful, unlike on Orion. **Centring moves the block number by only
+0.0015 and flips nothing.** It is not the cause.

## What this changes
1. **639's block-form prereg (ID05 / ID10; d4cd5f, d5f7b65b).** Its "held-out reading" (HCT116 +0.040, HEK293T +0.033
   block-gene signs) is artifact-dominated. Honest held-out block readings: **h1 +0.0057** (+0.0080 with expression
   partialled out) and **Broad-16 −0.0007**. Its predicted board size (+0.001 … +0.004) should come down roughly
   5–6× on the block, to about +0.000 … +0.001. The arm is still an exact prior-only read, and cheap, but it tests a
   small effect.
2. **Where the identity gain lives (h1, Broad-16): on UNCALLED genes, which only the uncalled-gene HOOK touches.** H1's
   real emitter gives KGID a TIE with the plain hook (612). So the identity's sign gain does not reach the board through
   the hook either, at today's hook magnitudes. That leaves 630's UP pass. It compounds with better UP signs (d4cd5f:
   −0.0105 vs −0.0080 raw at f 1 on the KGID base), so the route most likely to deliver the identity gain is the hook's
   signs with 630's magnitudes.
3. **Any Orion-fold result about HIGHLY EXPRESSED genes needs an expression-only control** before it is read as
   biology. That applies to the CPM block, to identity, and to the "identity genes go down" story in general.
   SCIENCE §2's sign result rests on all-gate-gene LOLO gains (h1, Broad-16 uncalled, Tahoe, 628), and those survive
   the I⊥E control. The block localisation does not.

Files: `split609.py`, `block609.py`, `work/split609.{log,json}`, `work/block609_gated_{centred,uncentred}.{log,json}`.

## 609b — the identity TAILS by region (`tails609.py`, `work/tails609.{log,json}`; 10.6 s, analysis only)
Question: does 622n's "low-identity tail goes UP" (found on UNCALLED genes) also hold in the CPM block, so that a
"low-I block genes → UP" override would be a cheap prior-only arm?

| surface | region | bottom-10 % I: truth UP / acc CLEAN4 → force UP | bottom-20 % I | top-10 % I: acc CLEAN4 → force DOWN | **expression-only bottom-10 %** (truth UP; acc C4 → force) |
|---|---|---|---|---|---|
| **h1** (clean) | block | 0.603 / 0.606 → 0.603 (−0.003; n 312) | **+0.021** (n 1,079) | 0.601 → 0.633 (+0.032; n 3,727) | |
| | uncalled | **0.639 / 0.557 → 0.639 (+0.082; n 2,562)** | +0.067 | 0.589 → 0.563 (−0.026) | **0.580 / 0.580 → 0.580 (0.000)** |
| hct116 (Orion) | block | 0.143 (n 7) | **−0.220** (n 177) | 0.695 → 0.764 (+0.069) | |
| | uncalled | 0.589 / 0.641 → −0.052 | +0.030 | −0.278 | **0.825**: +0.262 |
| hek293t (Orion) | block | 0.318: **−0.242** (n 223) | −0.450 | 0.640 → 0.805 (+0.165) | |
| | uncalled | 0.587: +0.030 | −0.001 | −0.146 | **0.675**: +0.147 |
| k562 | uncalled | 0.871: +0.306 (n 62) | +0.256 | −0.381 | 0.953: +0.375 (n 64) |

- **No block override.** On clean H1 the low-I tail adds only +0.02 in the block (bottom 20 %). On the Orion folds the
  block's low-I genes are mostly DOWN, the depth skew. Broad-16's block was ~0 (609 Part A).
- **On H1 the uncalled low-I → UP effect is IDENTITY, not expression.** The expression-only bottom decile is flat
  (0.580 = CLEAN4), while I's bottom decile is +0.082. So 622n's 630-route signal is real on the clean surface.
- **New caution, from both of the Orion artifact's tails:** on the Orion folds the LOW-expression uncalled genes skew UP
  (truth 0.83 / 0.68, expression-only "gains" +0.26 / +0.15). That is the mirror of the high-expression DOWN skew.
  CPM-normalising shallower perturbed cells inflates low counts, so low-expression lfc goes up. **Any 630 / UP-pass
  pricing on the Orion folds is therefore inflated by an artifact. Price UP-sign levers on H1, the Broad-16 and K562
  only.**
