# 421: CD4 dose on the gamma-0 carriers (H1). The shipped dose is the optimum; closed

Date: 2026-09-22 23:25–23:57Z. H1 six-member harness: 416's `run416.py`, copied with a new design set. The bulk is
410b's shipped cap15. Three emission seeds, and each arm is compared with 416's same-seed base. No builder, no vcc, no
cell-eval2.

**Question (416's suggestion):** 416's 2×2 found the CD4 residual worth about +0.006 of the DE channel. Does MORE CD4 in
the carriers pay more? Design: gamma-0 carriers with CD4 target f = 0.34 × {1, 1.5, 2, 2.5}, re-solved against the
unwhitened base. Rows are norm-matched to the shipped carriers, as in 416's g000.

**Known answers:** c034 == 416's g000 (sha 488732a0…). c034 seed 0 reproduces 416's g000_s0 metrics with max|d| = 0.0.

| design | CD4 f | Δnmae | Δfid | Δreach | Δjac | DE-only overall (board slopes) |
|---|---:|---:|---:|---:|---:|---:|
| c034 (= 416's g000, i.e. 417c's carriers) | 0.34 | −0.00683 ± 0.00004 | +0.00253 ± 0.00005 | −0.0021 ± 0.0079 | +0.00025 | **+0.00301** |
| c051 | 0.51 | −0.00658 | +0.00300 | −0.0175 ± 0.0097 | +0.00027 | +0.00032 |
| c068 | 0.68 | +0.00200 | +0.00308 | −0.0259 ± 0.0107 | +0.00032 | −0.00356 |
| c085 | 0.85 | +0.00985 | +0.00277 | −0.0356 ± 0.0150 | +0.00034 | −0.00771 |

**Reading:**
- More CD4 raises fid slightly (+0.0005 at ×1.5). It costs reach at every seed and dose, and past ×1.5 it costs nmae
  as well.
- The DE channel's CD4 optimum is the shipped dose. **417c's carrier design (γ 0, CD4 ×1) stands.**
- The bulk, by contrast, wants CD4 ×1.5 (418). The two channels have different optima, which only ballast lets us ship.
  That is the stack LOAD arm's bulk.
