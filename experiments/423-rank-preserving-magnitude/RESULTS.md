# 423 — rank-preserving carrier magnitudes: an exact-U pass on down-regulated "mass-holder" genes buys nmae −0.0023 … −0.0026 on H1 with the other five members bit-identical

Status: **COMPLETE** (2026-09-23 02:05–06:00Z; worker, Opus 5.5). H1 six-member harness only. No builder edits, no
`vcc`, no `cell-eval2`, nothing outside `experiments/423-rank-preserving-magnitude/`.

## Answer

The literal rank-free pass from the mandate reaches about 15% of the target LFC change, measured mass-weighted. It
reaches the ~3.5% of DE-universe genes that have any carrier cell above the control max or with interval headroom.
Almost all of the nmae prize sits in one class of those genes.

**The mass-holder genes.** These are strongly down-designed, low-expressed genes.
- The emitter delivers their per-gene bulk total by parking it on 3–9 carrier cells above the control max. Every
  other carrier is a zero tied with the control zeros.
- So the Wilcoxon sees a strong down gene, but the arithmetic-mean LFC reports the bulk total, not the rank design.
- On H1, the moved genes' sign accuracy on the reference DE genes is 0.65–0.74, against 0.56 overall. Their median
  |ℓ|/|r| is 0.33–0.63 (0.46 for the genes zs moves).

Up genes are the opposite case. Scaling them costs nmae: at k 1.4, +0.0005.

**Recommended rule: `dn_k250_zs`.** It is a zero-sum carrier-internal transfer. For every down gene, counts move from
low-library-class carriers to high-library-class carriers. Every cell stays strictly inside its own interval between
adjacent reference-control values; an above-max holder may fall to just above the max. The per-gene block total,
every library, the ballast and every cell's control rank are unchanged, so U, the tie structure, the calls and the
signs are unchanged.

| | Δnmae (3 seeds) | 95% CI over targets | other five members | per-gene U |
|---|---:|---|---|---|
| `dn_k250_zs` | **−0.00258 ± 0.00004** | [−0.0046, −0.0009] | bit-identical | identical for every gene; also under a float32 scorer path |
| `dn_k200` (pc, the literal pass) | −0.00375 ± 0.00014 | [−0.0090, +0.0005] | fid/jac identical, reach −2e-5 | the ballast crosses control values on ~100 genes/target |

- `dn_k250_zs` improves 34 targets and worsens 20. Its top 3 targets carry 55% of the gain.
- `dn_k200` has the larger mean, but 3 targets carry 95% of its gain and 27 targets are worse. It is a lottery on a few
  genes; `dn_k250_zs` is not.
- In overall units, `dn_k250_zs` is **+0.0007** on H1 (nmae slope 1.6532/6).
- **Predicted board nmae: −0.0023 … −0.0053 raw, i.e. +0.0006 … +0.0015 overall**, with every other member
  bit-identical (§6). This uses the production-bulk H1 value, −0.00232 (§7).

The arm is small, but on the board it would be a single-member read of the magnitude channel. No other member moves,
so the nmae change is exactly this pass's H1→board transfer ratio.

## 1. Harness and known answers

`run423.py` is 416's `run416.py`, restricted to the base carriers:
- the shipped cap15 carriers;
- CLEAN4 q80 prior;
- K-class emitter;
- ballast in ratio mode;
- pair-jackknife boost.

It adds one change: every post-pass variant is scored from the SAME emitted base blocks inside one process per seed.
The pass itself is `rp423.py`; its verification uses the scorer's own arithmetic (`rp423.verify`, `rp423.r1_tie`).

| known answer | result |
|---|---|
| `base_s0` vs 416 `metrics-base_s0.json` | max\|d\| 1.1e-16, then 0.0 on later runs (driver aborts above 1e-12) |
| `base_s1` / `base_s2` vs 416 | 2.2e-16 / 0.0 |
| `base_s{0,1,2}_bps` vs 416 `_bps` (production 410b bulk) | 0.0 / 0.0 / 1.1e-16 |
| `test_rp423.py` (real H1 blocks, targets 1 and 6) | 86/86 PASS; see list below |
| `test_rp423_synth.py` (synthetic blocks, K = 1, 3, 7 ballast rows) | 0 failed, 3,420 genes moved |

`test_rp423.py` checks:
- k = 1.0 is a bit-identical no-op in every mode;
- per-gene totals, libraries, non-negativity and no sign flip;
- untouched genes' R1 and tie term are exact;
- zs/ex: R1 is identical for every gene and the calls are identical;
- pc: every carrier keeps its exact control position.

Every run's verify record, per target and per variant (`work/h1/diag-*.json`), shows:
- block sums bit-equal (so pds is bit-identical);
- **untouched genes' rank sums unchanged in 100% of target×variant records**;
- **raw mse bit-identical**, even though the pass raises the block jackknife (the prediction's variance trace):
  - zs by +17% median (max +30%) on the cap15 bulk, and +7% on the production bulk;
  - pc by +2%.
  `expr_mse_unbiased_capped` uses min(pred trace, k × real trace), and on H1 our prediction already sits at or above
  that cap, so a larger trace changes nothing.
- On the board, where the cap may not bind for every target, a larger capped correction can only lower raw mse
  (better) or leave it unchanged; it cannot cost.

## 2. What the rank constraint allows

For each DE-universe gene, `rp423.cell_bounds` computes each carrier's exact count range [x_lo, x_hi]:
- inside that range, log1p(count/(lib/1e6)) stays strictly inside the same open interval between adjacent sorted
  reference-control values;
- the inversion is verified against numpy's log1p;
- a cell tied with a control value is frozen, e.g. a zero against control zeros;
- above the control max, x_hi is open.

**Library compensation.** Every count added to or removed from a tested gene in a carrier cell is taken from or given
to that cell's sub-gate genes: control CPM ≤ 5, not tested, not in nmae. So no other gene's CPM moves.

**Ballast bookkeeping.** When the pass routes counts to the ballast (modes pc/bal/ex):
- the ballast mirrors every gene total exactly (design mode);
- its libraries are restored by moving sub-gate counts between ballast rows.

**Feasible fraction**, i.e. achieved dA / target dA, per touched gene (`feas423.py`; seed 0, all_k140 pc):

| genes (per target) | any headroom | mean fraction | \|dA\|-weighted | full (≥ 0.95) |
|---|---:|---:|---:|---:|
| all 10,779 | 3.5% | 0.031 | 0.147 | 2.8% |
| up / down | 4.1% / 3.0% | 0.037 / 0.027 | 0.158 / 0.136 | |
| \|ℓ\| q0–0.5 | 1.7% | 0.015 | 0.140 | |
| \|ℓ\| q0.5–0.8 | 2.1% | 0.018 | 0.096 | |
| \|ℓ\| q0.8–0.95 | 3.1% | 0.023 | 0.129 | |
| **\|ℓ\| q0.95–1 (539 genes)** | **30.8%** | 0.296 | **0.408** | 28.1% |

Feasibility is all-or-nothing and sits in the top 5% of |ℓ|.

**Mass-holder census** (`holders423.py`, 8 dumped targets):
- about 150 genes per target have ALL their nonzero carriers above the control max (median 3–9 holders); ~40% of them
  are down, and another ~130 genes have some cells above the max;
- on the reference DE sets they carry sign accuracy 0.50–0.92 (5 targets with holder genes in S), against 0.47–0.58
  for the rest;
- median |ℓ| is 0.19–0.29, against ~0.04 for the rest;
- median control CPM is ~24, and ~90% of them sit outside the 2,200-gene CPM block. So the prior's C5b magnitude
  shrink is not what shrinks them. Their arithmetic mean is the design's bulk total, split in proportion to mass,
  while their rank design is a strong down.

Example: MYO6 in ACTN1 has 3 holders at ~2,700 CPM against a control max of 492 CPM. Its ℓ is −0.83; the reference
is −2.74. Pulled to the max, the holders alone would give ℓ −3.2.

**Two structural limits.**
- **The ballast crosses control values.** It is one deep cell (K = 1 on the cap15 bulk) sitting mid-distribution.
  Its gap to the next control value is a median 0.013 CPM, about 0.013 counts. Any count it absorbs crosses control
  values: that is the pc mode's U change.
- **The ballast's sub-gate reservoir.** The ballast holds only ~5,300 sub-gate counts (0.5% of its library), and its
  library headroom to 1e6 is ~5,000. So lowering moves that dump counts into the ballast are capped at ~5k counts per
  target. This binds for down rules: dn_k200 skips ~82 genes per target.
- zs needs neither: it is carrier↔carrier and reservoir-neutral.

## 3. The table (H1, cap15 bulk, 3 emission seeds, Δ vs same-seed base from the same blocks)

| variant | mode | seeds | Δnmae (± sd over seeds) | Δpds | Δmse | Δfid | Δreach | Δjac | total → overall | feasible (\|dA\|-w) | genes moved / target | genes with ΔU / target |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all_k080 | pc | 3 | +0.00222 ± 0.00007 | 0 | 0 | 0 | +0.00001 | 0 | -0.00061 | 0.151 | 385 | 368 |
| all_k120 | pc | 3 | -0.00139 ± 0.00016 | 0 | 0 | 0 | -0.00003 | 0 | +0.00038 | 0.158 | 375 | 363 |
| all_k140 | pc | 3 | -0.00195 ± 0.00015 | 0 | 0 | 0 | +0.00073 | 0 | +0.00068 | 0.147 | 375 | 370 |
| all_k160 | pc | 3 | -0.00196 ± 0.00010 | 0 | 0 | 0 | +0.00074 | 0 | +0.00068 | 0.139 | 375 | 373 |
| all_k200 | pc | 3 | -0.00143 ± 0.00015 | 0 | 0 | 0 | +0.00079 | 0 | +0.00054 | 0.122 | 373 | 373 |
| q50_k140 | pc | 3 | -0.00192 ± 0.00014 | 0 | 0 | 0 | +0.00005 | 0 | +0.00054 | 0.150 | 285 | 279 |
| q80_k140 | pc | 3 | -0.00187 ± 0.00014 | 0 | 0 | 0 | +0.00005 | 0 | +0.00052 | 0.201 | 217 | 211 |
| q80_k200 | pc | 1 | -0.00165 | 0 | 0 | 0 | +0.00010 | 0 | +0.00047 | 0.165 | 216 | 216 |
| up_k140 | pc | 3 | +0.00047 ± 0.00003 | 0 | 0 | 0 | +0.00072 | 0 | +0.00001 | 0.157 | 191 | 188 |
| dn_k140 | pc | 3 | -0.00236 ± 0.00018 | 0 | 0 | 0 | -0.00002 | 0 | +0.00065 | 0.112 | 154 | 151 |
| dn_k160 | pc | 1 | -0.00276 | 0 | 0 | 0 | -0.00002 | 0 | +0.00076 | 0.081 | 126 | 126 |
| dn_k180 | pc | 1 | -0.00325 | 0 | 0 | 0 | -0.00002 | 0 | +0.00089 | 0.062 | 110 | 110 |
| dn_k200 | pc | 3 | -0.00375 ± 0.00014 | 0 | 0 | 0 | -0.00002 | 0 | +0.00103 | 0.051 | 102 | 102 |
| dn_k250 | pc | 1 | -0.00269 | 0 | 0 | 0 | -0.00002 | 0 | +0.00074 | 0.034 | 91 | 91 |
| dn_k300 | pc | 1 | -0.00144 | 0 | 0 | 0 | -0.00002 | 0 | +0.00039 | 0.025 | 87 | 87 |
| dn_k400 | pc | 1 | +0.00123 | 0 | 0 | 0 | -0.00001 | 0 | -0.00034 | 0.017 | 81 | 81 |
| q80dn_k200 | pc | 1 | -0.00346 | 0 | 0 | 0 | -0.00002 | 0 | +0.00095 | 0.088 | 72 | 72 |
| q80dn_k250 | pc | 1 | -0.00272 | 0 | 0 | 0 | -0.00002 | 0 | +0.00075 | 0.072 | 71 | 71 |
| dn_k180_zspc | zspc | 1 | -0.00318 | 0 | 0 | 0 | -0.00002 | 0 | +0.00087 | 0.086 | 147 | 106 |
| dn_k200_zspc | zspc | 1 | -0.00352 | 0 | 0 | 0 | -0.00002 | 0 | +0.00097 | 0.070 | 139 | 104 |
| dn_k250_zspc | zspc | 1 | -0.00279 | 0 | 0 | 0 | -0.00002 | 0 | +0.00077 | 0.048 | 132 | 95 |
| all_k140_zs | zs | 3 | -0.00198 ± 0.00001 | 0 | 0 | 0 | 0 | 0 | +0.00055 | 0.054 | 259 | 0 |
| dn_k200_zs | zs | 1 | -0.00255 | 0 | 0 | 0 | 0 | 0 | +0.00070 | 0.029 | 121 | 0 |
| **dn_k250_zs** | zs | 3 | **-0.00258 ± 0.00004** | 0 | 0 | 0 | 0 | 0 | +0.00071 | 0.020 | 121 | 0 |
| dn_k300_zs | zs | 1 | -0.00255 | 0 | 0 | 0 | 0 | 0 | +0.00070 | 0.016 | 122 | 0 |
| dn_k500_zs | zs | 1 | -0.00255 | 0 | 0 | 0 | 0 | 0 | +0.00070 | 0.008 | 122 | 0 |
| all_k140_ex | ex | 3 | -0.00167 ± 0.00114 | 0 | 0 | 0 | 0 | 0 | +0.00046 | 0.101 | 327 | 0 |
| all_k200u120_ex | ex | 1 | -0.00213 | 0 | 0 | 0 | 0 | 0 | +0.00059 | 0.075 | 296 | 0 |
| dn_k200_ex | ex | 3 | -0.00159 ± 0.00051 | 0 | 0 | 0 | 0 | 0 | +0.00044 | 0.067 | 140 | 0 |
| dn_k250_ex | ex | 1 | -0.00128 | 0 | 0 | 0 | 0 | 0 | +0.00035 | 0.046 | 133 | 0 |

Base (mean of seeds 0–2; = 416 base): pds 0.88950, mse 0.97502, nmae 0.95212, fid 0.50902, reach 0.27495, jac 0.07474. "genes with ΔU" counts genes whose rank sum changed, all through the ballast (pc/zspc); zs and ex have 0 everywhere. "total → overall" sums all six members at board slopes (pds +2.2138, mse −1.051, nmae −1.6532, fid +3.3645, reach +1.1292, jac +2.7263)/6. Feasible = achieved/target dA, |dA|-weighted over touched genes (zs rows are measured against the k 2.5 target, which zs never reaches: it is capacity-bound).

Notes:
- **Global k** peaks at k 1.4–1.6 (−0.0020), well short of the unconstrained analytic ceiling (−0.0043 at k 1.3,
  `ideal423.py`). k 0.8 costs +0.0022.
- **q50 / q80 ≈ all.** The feasible genes are already the top of |ℓ|.
- **Up-only hurts (+0.0005); down-only helps (−0.0024 at k 1.4).** In pc, dn peaks at k 2.0 (−0.0037). dn k 2.5 / 3 /
  4 give −0.0027 / −0.0014 / +0.0012 (seed 0).
- **Reach.** The up-moving pc variants moved reach +0.002 on seed 2. That is the ballast's U change reordering p for
  one or two targets, within reach's 0.006 realisation noise. Down-only pc: −2e-5. zs/ex: exactly 0.
- **zs saturates at k ≥ 2.** dn_k200/250/300/500_zs all give −0.00255 (seed 0): the rule is "all available class-shift
  capacity on down genes", and k does not matter above 2.
- **ex is unstable** across seeds: all_k140_ex gives −0.0025 / −0.0004 / −0.0021 on seeds 0 / 1 / 2. It is zs plus ballast absorption with exact U restored by
  zero-sum crossing pairs or a single holder dropping into the controls' top tail. The exact-U solutions depend on
  integer coincidences in the ballast's crossing function, change with the emission realisation, and often land at A
  values that overshoot.

**Target-sampling uncertainty** (`boot423.py`; paired bootstrap over the 58 scored targets, seed-averaged):

| variant | mean Δnmae | 95% CI | p(≥ 0) | targets better / worse | top-3 share of gain |
|---|---:|---|---:|---|---:|
| **dn_k250_zs** | −0.00258 | [−0.00459, −0.00092] | 0.0001 | 34 / 20 | 0.55 |
| all_k140_zs | −0.00198 | [−0.00366, −0.00061] | 0.0003 | 33 / 23 | 0.59 |
| dn_k140 (pc) | −0.00236 | [−0.00475, −0.00040] | 0.008 | 33 / 22 | 0.71 |
| dn_k200 (pc) | −0.00375 | [−0.00900, +0.00050] | 0.047 | 25 / 27 | 0.95 |
| all_k140 (pc) | −0.00195 | [−0.00436, +0.00010] | 0.032 | 29 / 28 | 0.85 |
| q80_k140 (pc) | −0.00187 | [−0.00432, +0.00020] | 0.041 | 24 / 31 | 0.89 |
| all_k140_ex | −0.00167 | [−0.00364, +0.00009] | 0.032 | 33 / 24 | 0.70 |
| dn_k200_ex | −0.00159 | [−0.00499, +0.00152] | 0.16 | 28 / 27 | 1.14 |

Targets with a loss above 0.005: **0** for dn_k250_zs in each of the three seeds, against 5–7 for dn_k200.

**Seed-0 follow-ups (`queueF`): zs followed by pc on the remainder (`zspc`), and finer down k.**
- zspc does not beat pc: dn_k200_zspc −0.00352 against −0.00358.
- Every ballast-using down rule is lumpy (bootstrap on seed 0): dn_k160 CI [−0.0061, −0.0001], dn_k180
  [−0.0076, +0.0003], dn_k200_zspc [−0.0089, +0.0008], q80dn_k200 [−0.0088, +0.0008]. The top 3 targets carry 84–100%
  of the gain.
- The ballast route buys a larger mean on a few genes with huge reference LFCs and loses a little on many others. zs is
  the only construction whose gain is spread: 34 targets better, 0 worse by more than 0.005.

## 4. Library-size side effects (the mandate's question), measured

- **With the pass as built** (sub-gate compensation): carrier libraries are identical. Untouched genes' R1 and tie
  term are identical in every record.
- **Counterfactual without compensation** (`libfloat423.py`: the same moves, carriers' sub-gate counts restored so
  their libraries float by the kept-gene change, max 0.3–0.5% per cell):
  - 97–98% of UNTOUCHED genes change U;
  - max |ΔU| 1,900–12,000 half-ranks, i.e. Δz up to ~0.06;
  - sum |ΔU| 1–9 × 10^6 per target;
  - no call flips on the two targets tested.
- Compensation is required for a clean pass.

## 5. Robustness of "U identical" to the scorer's float path

H1's evaluator ranks float64 values in scanpy order: count / (lib/1e6), then log1p. cell-eval2's DE backends can
stage float32 values. `f32check423.py` recomputes every gene's rank sum with float32 values:

| variant | genes whose U moves |
|---|---|
| `dn_k250_zs` | **0 on all 8 targets** |
| `all_k140_ex` | 2 per target, by 1 half-rank |
| pc | every touched gene (the ballast), max 12.7k ranks |

## 6. Recommended rule and predicted board range

**Rule: `dn_k250_zs`.** It is k-independent above 2, so in words: "for every down gene of each block, move the
above-floor mass of low-library carriers into higher-library carriers, all inside their control intervals".
- H1 Δnmae −0.00258 (sd 0.00004 over seeds; target CI [−0.0046, −0.0009]).
- pds, mse, fid, jac and reach are bit-identical on H1. On the board mse can only stay or fall (§1).
- It needs no ballast and no reservoir.

**Predicted board nmae: H1 × 1 … 2.3.** On the production bulk (H1 −0.00232, §7) that is −0.0023 … −0.0053 raw → **+0.0006 … +0.0015 overall**. On the cap15 bulk it is −0.0026 … −0.0059. Other members unchanged.

Why the ×1 … ×2.3 window:
- It is a pure magnitude change on genes whose delivered sign the board already scores. Pure magnitude changes are the
  one class where the board has read H1's nmae sign correctly: 410's weld, H1 +0.024 vs board +0.054, about 2.3×.
- The two H1 mis-signings (392, 417c) were carrier-DESIGN changes that moved signs and magnitude structure together.
  This pass moves neither the design nor any sign.
- The floor ×1 assumes the A/B/C holder genes are as sign-reliable as H1's.

Risks:
- The gain sits on ~15 reference-DE genes per target (median), about 5% of the nmae denominator (`share423.py`). A
  context with fewer mass-holder genes, or a less reliable down sign on them, lands below ×1.
- The ×2.3 was measured on a harm, not a gain.

Price it as nmae-only. Built on 417a with the same seeds, the base blocks are identical and the pass moves only
down-gene carrier counts. So pds, fid, jac and reach come out bit-identical to 417a's, and mse identical or lower
(§1). The board's nmae delta is then the pass's effect alone, with no realisation noise: a clean calibration read of
H1→board magnitude transfer for this gene class. The upside variant is `dn_k200` pc, H1
−0.0037, but its target CI crosses 0 and it moves U on ~100 genes per target (see §3); it is not recommended.

## 7. Production bulk (410b's key-scaled bulk, `--bulk ps`, K = 1–16 ballast rows)

These runs use the production bulk: 410b's key-scaled bulk, the construction 417a ships before sharpening.
- Known answers: base_s0/s1/s2_bps vs 416's `_bps` rows: max\|d\| 0.0 / 0.0 / 1.1e-16.
- Base (seed mean): pds 0.89486, mse 0.92529, nmae 0.95232, fid 0.50920, reach 0.27382, jac 0.07480; this equals 416's production-bulk base.

| variant | mode | seeds | Δnmae (± sd over seeds) | Δpds | Δmse | Δfid | Δreach | Δjac | total → overall | feasible (\|dA\|-w) | genes moved / target | genes with ΔU / target |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| dn_k200 | pc | 3 | -0.00369 ± 0.00012 | 0 | 0 | 0 | +0.00001 | 0 | +0.00102 | 0.092 | 156 | 153 |
| dn_k200_zspc | zspc | 3 | -0.00376 ± 0.00011 | 0 | 0 | 0 | +0.00002 | 0 | +0.00104 | 0.101 | 171 | 151 |
| all_k140_zs | zs | 3 | -0.00156 ± 0.00012 | 0 | 0 | 0 | 0 | 0 | +0.00043 | 0.056 | 263 | 0 |
| **dn_k250_zs** | zs | 3 | **-0.00232 ± 0.00005** | 0 | 0 | 0 | 0 | 0 | +0.00064 | 0.021 | 123 | 0 |

| variant (production bulk) | mean Δnmae | 95% CI over targets | p(≥ 0) | targets better / worse | top-3 share |
|---|---:|---|---:|---|---:|
| **dn_k250_zs** | −0.00232 | [−0.00429, −0.00074] | 0.0003 | 32 / 21 | 0.60 |
| all_k140_zs | −0.00156 | [−0.00331, −0.00003] | 0.023 | 30 / 25 | 0.76 |
| dn_k200 (pc) | −0.00369 | [−0.00896, +0.00058] | 0.049 | 26 / 28 | 0.96 |
| dn_k200_zspc | −0.00376 | [−0.00889, +0.00054] | 0.046 | 27 / 28 | 0.94 |

Readings:
- The production bulk reproduces the cap15 picture.
- zs is unchanged in kind, at 90% of its cap15 size. zs never touches the ballast; the carriers differ only because K
  differs.
- With K = 1–16 ballast rows, the pc reservoir is K× larger, so feasibility doubles (0.092 against 0.051). But the pc
  gain does not grow (−0.0037), and it stays lumpy: 3 targets carry 96% of it.
- Every other member is bit-identical for zs. pc moves reach by +1e-5.

## 8. Builder spec (not applied; builder untouched)

**Where.** In `scripts/make_submission_pm_cd4.py`'s block loop (~l. 5494), right after
`block, ballast_diag = ballast_rows(...)` and its shape check, before `rc_emit.to_csr(block)`.
- zs never touches the ballast rows. In principle it works without `VCC_PERT_SCALE_BALLAST` too, but `rp423.postpass`
  as written expects at least one ballast row; a K = 0 path needs a guard around `b_dep`, `ball_sub` and `_subgate`,
  which zs never uses. The arm as specified is 417a + the pass, so K ≥ 1.
- It needs a K-class carrier block: `rc_hybrid_block_k` with ≥ 2 library classes, i.e. `VCC_DE_SIGN_PRIOR` set.
- Call: `block, mag = rank_mag.postpass(block, n_car=cells_per_pert - ballast_k, cc=rc_scorer_control,
  lfc_base=<block lfc>, spec={"sel": "dn", "k": 2.5, "mode": "zs"})`.
  - Pop `_fr_by_gene` / `_touched` from `mag` before merging it into `rc_diag`.
  - `lfc_base` = log2((A + 1e-9)/(ref_lin + 1e-9)), where A is the block's per-gene mean of count/(lib/1e6) over the
    kept genes (the function re-checks it to 1e-9).
- Knob: e.g. `VCC_DE_MAG=dn_k250_zs` (the `rp423.parse_variant` grammar). Unset = no call = bit-identical.

**What it needs.**
- A scorer-order control block per context, `exact_de323.control_block` semantics: keep = control mean CPM > 5, cs =
  sorted log1p(count/(lib/1e6)), ref_lin. Build it at ~l. 5315 before `control_matrix` is released: +1.6 GB per
  context.
- The emitter's `rc_control` uses `count*(1e6/lib)`, which differs by ulps. Reusing it would move ranks only when a
  value lands within an ulp of a control. §5 suggests zs is immune even to float32, but the exact block is the safe
  choice.

**Code.** Vendor `rp423.py`'s `cell_bounds`, `fill`, `_levels`, `_zs_step`, `_take`, `_subgate`, `settle_ballast`,
`postpass` and `select`; zs needs no `_offset` / `_drop_dump`.
- Runtime ~0.7 s per block in numpy, ~10 min per 900-block build.
- Asserts inside: per-gene totals, libraries, non-negativity, U per gene (mode zs), no sign flip.

**Known-answer tests before a build.**
1. `test_rp423.py` and `test_rp423_synth.py` pass against the vendored module: k = 1.0 is a bit-identical no-op, and
   the invariants hold.
2. H1 retrodiction: `run423.py 0 --variants dn_k250_zs` with the vendored function swapped in must reproduce
   `work/h1/stats-dn_k250_zs_s0.npz` bit-for-bit and `metrics-dn_k250_zs_s0.json` to ≤ 1e-12.
3. The builder's scorer-order control block for the H1 controls must equal 323's `reference.npz` `cs`/`keep`
   (sha256).
4. On the built package, per-target block column sums are identical to the pass-off build of the same seed. This is
   exact by construction and asserted per block, so it is cheap to also record in `postbuild.json`. The rc audit
   carries per block: genes moved, `n_genes_U_changed` = 0, and the feasible fraction.
5. The emitted A/B/C blocks' mass-holder census must look like H1's before the arm is priced; the arm is priced on
   holder genes. Use `holders423.py` on a few dumped blocks: ~150 pure-holder genes per target, ~40% down.

## 9. Negatives and closed doors (with mechanism)

- **Up genes: closed.** Scaling them costs nmae on H1 (up_k140 +0.0005 at 3 seeds). The feasible up genes are already
  at |ℓ|/|r| ≈ 0.65 with sign accuracy 0.61, so more magnitude overshoots.
- **Global k is capped by sign accuracy.** Even unconstrained, k 1.3 gives only −0.0043 (412 found the same). Under
  the rank constraint the best global rule is −0.0020.
- **ex (exact U with the ballast absorbing) is closed.** It is unstable across seeds, and its target CI crosses 0. The
  exact offsets exist only through integer coincidences.
- **The literal pc pass is not recommended.**
  - The ballast (one deep mid-distribution cell) crosses control values on every touched gene.
  - Its down rules are also capped by the ballast's ~5k-count sub-gate reservoir.
  - Its gain rests on 3 targets.
- **Library-floating variants** would move 97–98% of untouched genes' ranks (§4).

## 10. Compliance

- Serial H1 arms (`drive423.sh`, chained by `chain423.sh`), each started only at `memory_pressure` ≥ 40% free and
  free+inactive ≥ 20 GiB. The gate never held.
- Peak RSS 5.2–6.3 GB per base-bulk run (`/usr/bin/time -l` in `work/logs/run_*.log`); ps runs in §7.
- Small analysis scripts (`quick423*.py`, `f32check423.py`, `libfloat423.py`, `holders423.py`; RSS not measured
  individually, a few GB each) ran alongside H1 arms. System memory stayed at 70–77% free throughout.
- Once, a heredoc wrote `f32check423.py` into the repo root because the shell's working directory had drifted; it was
  moved into this directory within a minute, and nothing else was written outside it.
- Superseded seed-0 ex results (before the zs-fallback fix) are kept in `work/h1/superseded/`.
- **Two restarts of the production-bulk queue.** Both were rp423 K > 1 bookkeeping bugs, caught by the pass's own
  asserts; neither is a scoring error, and no result was taken from a failed run.
  1. The reservoir row was chosen before the pass and re-chosen after it.
  2. The even-split remainders all landed on ballast row 0, so one row owed more sub-gate counts than it held.
- **The fixes:**
  - one deepest row is fixed for the pass;
  - the reservoir is K × the poorest row's sub-gate holdings, minus 64 per row;
  - rows that gained kept counts refund sub-gate counts first;
  - remainders rotate by gene, as `ballast_rows` does.
- **Checks after the fixes:**
  - `dbgps423c.py`: 60 production-bulk targets × 4 variants, 0 failures;
  - `test_rp423_synth.py` passes;
  - `regress423.py`: the K = 1 (cap15) results reproduce bit-for-bit on sampled targets, so the base-bulk table is
    unaffected.
- The production-bulk runs peak at 7.1–7.3 GB RSS, under the 8 GB cap.

## Files

| file | what |
|---|---|
| `rp423.py` | the pass, modes pc/zs/bal/ex/zspc, `verify`, `r1_tie` |
| `run423.py` | harness |
| `drive423.sh`, `chain423.sh`, `work/queue*.txt` | serial drivers |
| `analyse423.py`, `table423.py` | tables |
| `boot423.py` | target bootstrap |
| `feas423.py` | feasible fraction |
| `contrib423.py`, `cmpfull423.py`, `share423.py` | gene-level decomposition |
| `holders423.py` | mass-holder census |
| `f32check423.py` | float32 scorer path |
| `libfloat423.py` | library-size counterfactual |
| `ideal423*.py` | unconstrained ceilings |
| `test_rp423.py`, `test_rp423_synth.py` | invariant tests |
| `regress423.py` | K = 1 bit-reproduction of the stored runs |
| `dbgps423*.py` | K > 1 checks on production-bulk dumps (regenerate with `run423.py 0 --bulk ps --dump 60`) |
| `explore423.py`, `quick423*.py`, `dbg423*.py`, `prof423.py`, `hang423.py`, `cmp423.py`, `genes423.py` | exploration |

Outputs are `work/h1/{metrics,diag,per-target,stats,genes}-<variant>_s<seed>[_bps].*`, the table in
`work/analysis423*.{md,json}`, and logs in `work/logs/`, `work/drive*.log`.
