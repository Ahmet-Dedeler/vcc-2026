# 345 — A same-line common vector holds +0.056 of sign accuracy; no transferable version carries any of it

**Status: NEGATIVE (as a design lever). COMPLETE as a measurement.** Descriptive analysis on saved
arrays only. No emission, no 323 scoring run, no builds.

## Known-answer reproduction

Gate rebuilt exactly per 343 (`p_adj<.05` & finite nonzero lfc & gene ≠ target from
`reference_de.parquet`; 58 targets with ≥10 gate genes). Design mapped by gene name from
`effects.npz::cap15` with the 2 absent universe genes (`TBCE-1`, `HSPA14-1`) as 0:

| quantity | mine | record | match |
|---|---|---:|---|
| \|lfc\|-weighted sign acc, delivered | 0.58796796 | 0.58796796 (343 `signs-cap15.json` `w_all`) | to 1e-9 |
| \|lfc\|-weighted sign acc, designed | 0.60018510 | 0.60018510 (343 `design_w_all`) | to 1e-9 |
| mean NMAE | 0.95083646 | 0.95083646 (323 `metrics-cap15.json`) | exact |
| magnitude ratio m | 0.20548496 | 0.20548496 (343 `summary-cap15.json`) | exact |

## What was measured

- **R**: 85×10,779 matrix, R[t,g] = reference lfc if g in gate(t) else 0. **D/Dsup**: cap15 design
  on the universe / restricted to the gate support. SVD energy shares; |cos| of first singular
  vectors (sign reported too).
- **LOO common vectors from H1 truth** (mean over the other 84 targets; sign-vote; per-gene
  |lfc|-weighted mean), read by sign on the held-out target's gate; constant all-down on the same
  gate; gene-shuffle null (200 seeds, permute gene labels of the LOO mean).
- **Combination** sign(D[t]+α·s·c_{-t}), s = pooled-median |D| / pooled-median |c|
  (s=2.7549 oracle-mean), α ∈ {0,.25,.5,1,2,4,∞}; weighted acc, implied NMAE ≈ 1−m(2p−1) at
  m=0.20548, and EXACT NMAE with delivered magnitudes and combination signs.
- **Transfer vectors**: design's own LOO mean; donor-bank all-target means (K562 GWPS 9,787 rows,
  Orion HCT116 17,767, HEK293T 17,977; `delta.npy` via mmap, chunked means; sign convention
  verified pert-minus-control in `scripts/build_donor_table.py:157`, so the signs below are not a
  flip bug); target-matched donor means (82–84 rows in the H1 85-panel); common-vector cosines.
- Splits: ref UP vs DOWN, basal-CPM quintiles (343 `bins343.py` shape).

## Results

### 1. The reference signal has a large shared component; the design does not share its direction

| k | R energy ≤k | Dsup energy ≤k |
|---|---|---:|
| 1 | 0.2191 | 0.1401 |
| 2 | 0.3172 | 0.2622 |
| 3 | 0.3829 | 0.3774 |
| 5 | 0.4924 | 0.4870 |
| 10 | 0.6828 | 0.6577 |

Rank-1 carries **21.9%** of reference gate energy. But the shared component is NOT the flat mean:
mean-direction share is 0.0241 (R) / 0.0128 (D). First singular vectors R-vs-D: left cos −0.2171
(|.| 0.2171), **right cos −0.0524 (|.| 0.0524)** — the dominant gene direction of the design is
unrelated to the reference's. Pool down-share on this H1 gate: 0.5248 unweighted / 0.5403 weighted
(less down-biased than 155's 0.597 pool — different surface, as expected).

### 2. The LOO common vector beats the design — gene-specifically, not marginally

| vector (sign read on held-out gate) | weighted | unweighted |
|---|---|---:|
| design alone (baseline) | 0.6002 | 0.5694 |
| LOO mean | **0.6170** | 0.6054 |
| LOO sign-vote | 0.5745 | 0.5666 |
| LOO \|lfc\|-weighted mean | 0.6133 | 0.5982 |
| constant all-down (this gate) | 0.5272 | 0.5047 |

Gene-shuffle null for the LOO mean: 0.4781 ± 0.0101, **z = +13.7** — the gain is gene-specific
structure, not the marginal. (155's 0.5974 was unweighted on a different pool; the comparable
number here is 0.5047. The uniform-shift closure is untouched: this vector is gene-specific.)

### 3. Oracle ceiling: +0.056 weighted, −0.015 exact NMAE (mean variant, α=1–2)

| α | wacc | uacc | NMAE exact | NMAE implied |
|---|---|---:|---:|---:|
| 0 (design) | 0.6002 | 0.5694 | 0.9485 | 0.9588 |
| 0.25 | 0.6418 | 0.6034 | 0.9420 | 0.9417 |
| 0.5 | 0.6527 | 0.6186 | 0.9358 | 0.9373 |
| **1** | **0.6564** | 0.6246 | **0.9338** | 0.9357 |
| 2 | 0.6573 | 0.6306 | 0.9338 | 0.9353 |
| 4 | 0.6529 | 0.6326 | 0.9353 | 0.9372 |
| ∞ (common alone) | 0.6170 | 0.6054 | 0.9639 | 0.9519 |

Vote and wmean variants peak at 0.6589 / 0.6556 (same shape; full grids in `work/results.json`).
Implied identity tracks exact within 0.002. (α=0 exact 0.9485 ≠ delivered 0.9508 because α=0 pairs
DESIGNED signs with delivered magnitudes — stated, not a mismatch.) The oracle gains concentrate
on ref DOWN genes (0.6351→0.7333, +0.098) vs UP (0.5774→0.5885, +0.011); all basal quintiles gain
+0.05–0.08, including the lowest-CPM quintile holding 36% of NMAE mass (+0.051).

### 4. Transfer: every design-time common vector scores ≤ design alone — kill criterion fires

| vector | alone (w) | best combo α (w) | verdict |
|---|---|---:|---|
| design's own LOO mean | 0.4957 | 0.6002 at α=0 (all α>0 worse, monotone) | hurts |
| K562 GWPS all-mean | **0.2753** (anti-predictive) | 0.6002 at α=0 | hurts |
| K562 target-matched mean (82 rows) | 0.2849 | 0.6001 at α=0.25 | hurts |
| Orion HCT116 all-mean | 0.5087 | 0.6002 at α=0 | hurts |
| Orion HCT116 matched (84 rows) | 0.5054 | ≤0.5874 | hurts |
| Orion HEK293T all-mean | 0.5277 | 0.6002 at α=0 | hurts |
| Orion HEK293T matched (84 rows) | 0.5295 | 0.6028 at α=0.5 | **+0.0026, grid-max, 8× below gate** |

Common-vector cosines: donor means correlate with EACH OTHER (K562×HCT116 +0.361,
HCT116×HEK293T +0.406, K562×HEK293T +0.158) but barely with H1's LOO mean (+0.168/+0.104/+0.074);
design's own mean is orthogonal to everything (|cos| ≤ 0.06). Each line has a common program;
none of them is H1's — K562's is sign-inverted on the weighted gate.

## Verdict against the PLAN.md kill criterion

**The idea is DEAD as a design lever.** The best transferable combination gains +0.0026 weighted
sign accuracy over the design alone (HEK293T matched mean, α=0.5 — a grid maximum, uncorrected
for selection), against the +0.02 gate. All other transferable constructions are at or below
α=0. Reported separately as required: the same-line ceiling is +0.056 weighted / −0.015 exact
NMAE — i.e. a same-line prior could buy what the gap needs, but no design-time vector we hold
carries any of it. This is consistent with 246's tension note in the opposite direction:
whitening strips shared directions (good for `pds`), and the shared direction that matters for
signs is line-specific and not transferable from any donor bank, target-matched or not.

**Worth a 323 six-member run: no, with no construction.**

## Limits

H1 development panel, one line, 58 targets; weighted sign accuracy is not the official
`fid`/`reach`/`nmae`. Combination grids are sign-only (magnitudes fixed at delivered). Donor
means average over all bank rows including essentials; matched means over only 82–84 rows.
K562's 0.275 anti-accuracy is descriptive — mechanism (line-specific stress programmes vs
technical axes) not investigated; it is not needed for the kill.

## Exact commands to reproduce

```bash
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/345-common-program-structure/run345a.py  # steps 0-2 + stage345.npz
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/345-common-program-structure/run345b.py  # steps 3, 4a/b, 5-transfer
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/345-common-program-structure/run345c.py  # matched means, cosines, oracle splits, results.json
```

All three deterministic (fixed seeds); reruns reproduce bit-exactly (verified: two consecutive
full runs identical, including NMAE 0.95083646).

## Peak memory and wall time

`/usr/bin/time -l`: run345a RSS 0.29 GB; run345b 3.04 GB; run345c 3.06 GB (mmap'd Orion
`delta.npy` pages counted resident; no array densified; cap was 12 GB). Free ≥60% at every
launch. Wall: low single-digit minutes per script (not clocked precisely). Scratch: `work/`
holds `stage345.npz` (~15 MB), `results*.json`; no other large files created.

## PI re-derivation note (2026-09-10 01:05Z, Claude)
The "K562 all-mean 0.2753 anti-predictive" row does not reproduce. With independent code (sign convention verified:
own-target delta < 0 on 99.3% of rows), the K562 all-row mean scores **0.527 macro / 0.547 pooled** |lfc|-weighted
sign accuracy on the H1 gate; Jurkat common (350) 0.534 / 0.553; HepG2 common 0.529 / 0.530; the design mean
0.518; the H1 LOO truth mean reproduces this report's **0.6170** exactly. The kill verdict is unchanged (no
donor common vector approaches the design's 0.600), but the anti-alignment claim is withdrawn. Likely cause: an
axis/sign mishandling in `run345b.py`'s donor-mean path; not traced further.
