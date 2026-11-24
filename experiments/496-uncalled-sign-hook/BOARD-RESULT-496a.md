# 496a BOARD RESULT: 0.2728318 (+0.0005552 vs 446bj2). The plain CLEAN4 uncalled-gene hook TRANSFERS (Δnmae −0.00158 raw)

Entry `HoQY4xrzYer5y232oh85`, rank 16. Submitted 2026-10-05 00:01:41Z by `arm1005.sh` slot 1 (`chain1005`, rule1005 `P=P2g,C=TRANSFERS,O=OK`;
-m from the queue, -d empty). Upload done ~00:21Z on Wi-Fi; status "launching" until ~01:03Z; captured 01:03:30Z. Package
`submissions/496a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-as05.vcc`
(sha256 `3f5cc1ec…`, prediction identity `e21ff752…`). Base 446bj2 `YvVru3DAF89Jiyr34olR` 0.2722766. Arm: 446bj2 +
`VCC_DE_MDE_SIGN=clean4:0.05`, config diff exactly `de_mde_sign` + `de_sign_emitter_sha256`. Bins: `PREREGISTRATION-board.md`
(2026-09-27 06:57Z). Reader: `experiments/663-idall3-arm/read1005.py --a496 HoQY4xrzYer5y232oh85` (written 10-04 ~20:15Z, before the
entry existed). Read 2026-10-05 ~01:05Z.

| member | Δ raw vs 446bj2 | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | 0 (bit-identical) | 0 | 0 |
| mse | −4.6e-08 (0.6264081539 → 0.6264081078) | +0.000000 | +0.0000000 |
| nmae | **−0.001577** | +0.002423 | +0.0004038 |
| fid | −0.000024 | −0.000082 | −0.0000136 |
| reach | +0.000901 | +0.000988 | +0.0001647 |
| jac | +0.0000003 | +0.000002 | +0.0000004 |
| **overall** | | | **+0.0005552** |

## Against the bins
- **Falsifier: PASS.** pds is bit-identical. mse moved −4.6e-8 raw, so it is not bit-identical. The prereg's wording is
  "exactly"; it was written on 09-27, before 592's 10-01 correction recorded that DE-only arms move the scorer's jackknife-corrected
  mse by O(1e-8) raw. Five earlier pairs show this: 446bja +2.4e-8, 458d / 458g −5.3e-8, 446bj2zcp −2.3e-8, K20 +4.1e-8. The
  project's DE-arm falsifier is ≤ 1e-6 raw (537 §5.2; 535/BOARD-RESULT-446bj2zcp; the frozen 10-07 rule's R3). −4.6e-8 is
  20× inside it. Mechanism, from the scorer's source (`cell_eval2/metrics/delta.py`, `mse_unbiased`): the mean `mu` is built from
  the group SUM, which a DE-only arm leaves unchanged. The subtracted correction `C_pred` is a delete-1 jackknife over the
  prediction's cells, and the hook redistributes counts across cells. read1005.py first printed FAIL at a 1e-12 tolerance; it was corrected to 1e-6 (still printing bit-identity)
  at 01:10Z, before idall3's read.
- **Leak check: clean.** |Δfid| 2.4e-5, |Δjac| 2.5e-7 (both ≤ 0.001 raw); |Δreach| 0.0009 (≤ 0.01).
- **BIN: TRANSFERS** (Δnmae −0.00158 ≤ −0.0010). Per the prereg, October adds the hook (`VCC_DE_MDE_SIGN=clean4:0.05`, through
  `DEF_DAY_MDE_SIGN`, with its own rehearsal).

## Against the registered predictions
| source | predicted Δnmae raw | read −0.00158 |
|---|---|---|
| 496's prereg (H1 f 1 × 1–2.4) | −0.0015 … −0.0035 | inside, at the small end |
| PRED-662-496a, I-662 on 618 (f\* 0.61 × 1.05–1.30) | −0.0014 (−0.0017 … −0.0011) | **inside** |
| PRED-662-496a, 618 at f 0.35 (619's rule) | −0.0014 … −0.0017 | inside |
| PRED-662-496a, I-662 on h1-test (f\* × 0.92–1.37) | −0.0022 (−0.0027 … −0.0016) | just outside (2e-5 above the band); h1-test over-reads 1.4× |
| 613d (registered 10-02; its discrimination rule WITHDRAWN 10-02 19:10Z) | +0.00017 | missed; recorded only, no conclusion drawn |

- **618 reads the hook; h1-test over-reads it** by about 1.4× (−0.0022 predicted vs −0.00158 read).
- **CPH's T** = Δnmae / −0.001979 = **0.797**, below the I-662 band's mid of 1.12 (0.92 … 1.37). PRED-662-CPH's formula, fixed
  before this read: the I-662 primary is +0.000575 + 0.000273 T. That gives **+0.00079 (P8)**, +0.0017 (P1) and −0.0003 (P9)
  overall. On f 1 (537's denominator −0.001458, T 1.08), P8 gives +0.00067. Every value is below CPH's SHIP bar (+0.0012), so
  the block still does not ship (652 F2).
- reach +0.0009 raw (+0.00016 overall) is unpriced, as the prereg said, and positive.

## Consequences
- **October `mde_sign` → `clean4:0.05`.** The board read is TRANSFERS. The held-out reads are H1 (3 seeds, nmae −0.00146 at
  f 1) and 618 (−0.00123 at f\*). Class: transfer (a cross-line CLEAN4 sign prior for uncalled genes, delivered through the MDE
  channel). Tries: 1 board arm for the plain hook. The `october_recipe.py decide` runs only after recipe_dr524 has resolved its
  recipe: 658's 524dr2 design assumes 524dr ran mde_sign OFF (PREREGISTRATION-524dr2 line 23), and 10-05 CHECKLIST v2 step 10
  decides it after dr524.
- **idall3's anchor:** Δ_496a = Δ(fid + nmae) overall = **+0.000390** (fid −0.0000136, nmae +0.0004038).
- The displayed row is now 496a (0.2728318). The best scored entry is still 498p2g (0.2749190).
