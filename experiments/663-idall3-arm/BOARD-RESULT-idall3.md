# idall3 BOARD RESULT: 0.2746732 (−0.0002458 vs 498p2g). The ungated identity hook is BELOW the plain CLEAN4 hook; October keeps the plain hook

Entry `4R0FMEtcLtvZMmg2X0G5`, rank 16. Submitted 2026-10-05 01:04:17Z by `arm1005.sh` slot 2 (the coordinator's fill through
`fill1005s2.sh`; -m from the queue, -d empty), after 496a's capture. Captured 01:37:03Z. Package
`submissions/idall3-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5-idhook.vcc`
(sha256 `dad0e10e…`). Base 498p2g `JWSghc0o7MN5D6o8Iadc` 0.2749190. Arm: 498p2g + `VCC_DE_MDE_SIGN=file:<604 v3 identity artifact,
λ 0.5, 6cd84dc2…>:0`, no KEYGATE; config diff exactly `de_mde_sign` + `de_sign_emitter_sha256`. Bins:
`experiments/612-key-gated-hook/PREREGISTRATION-draft-IDALL3.md` (FINAL, cf3f8004). Predictions: `experiments/661-618-board-preds/PRED-661-idall3.md`
(618 TIE, plus its I-662 alternative). Reader: `read1005.py --a496 HoQY4xrzYer5y232oh85 --idall3 4R0FMEtcLtvZMmg2X0G5`.
Read 2026-10-05 ~01:38Z.

| member | Δ raw vs 498p2g | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | 0 (bit-identical) | 0 | 0 |
| mse | +9.9e-08 (jackknife term; see 496a's record) | −0.000000 | −0.0000000 |
| nmae | **+0.000274** (worse) | −0.000641 | −0.0001068 |
| fid | −0.000005 | −0.000025 | −0.0000041 |
| reach | −0.000701 | −0.000772 | −0.0001287 |
| jac | −0.000013 | −0.000037 | −0.0000062 |
| **overall** | | | **−0.0002458** |

## Against the bins
- **Falsifier: PASS.** pds is bit-identical; mse moved +9.9e-8 raw (≤ 1e-6, the DE-only jackknife move); |Δjac| 1.3e-5 raw.
- **Primary** Δ(fid + nmae) overall = **−0.000111**. Against Δ_496a = +0.000390 (496a − 446bj2, the plain hook's own effect, read
  the same night): **−0.000501 → BELOW** (< −0.0002). Not HARM: the primary against its own base is −0.000111, above −0.0002.
- Reach −0.0007 raw is inside the ±0.0042 band. That band gates only the WIN bin, so it plays no role here.
- **Bin consequence:** A/B/C and October both use the plain CLEAN4 hook. October's `mde_sign` stays `clean4:0.05` (496a TRANSFERS),
  and the identity hook is not adopted.

## Against the registered predictions (primary − Δ_496a, overall)
| source | predicted | read −0.000501 |
|---|---|---|
| h1-test (the prereg) | WIN +0.0007 … +0.0010 | over-predicted by ≈ 0.0012 |
| 618 at f 0.35 (PRED-661-idall3, registered TIE) | +0.00007 (−0.00002 … +0.00015) | over-predicted by ≈ 0.0006 |
| 618 I-662 alternative (f\* 0.61) | +0.00020 (+0.00006 … +0.00028) | over-predicted by ≈ 0.0007 |

- Every instrument over-predicted the identity hook. 618 was the closest (TIE, against a BELOW read), as it also was for 496a.
- On the board, the identity hook does WORSE on nmae than no hook at all (+0.00027 raw vs 498p2g). The plain hook gains −0.00158
  raw on 446bj2. Both are hooks on uncalled genes through the same MDE channel. Only the sign source differs: CLEAN4 versus
  z(P_C4 − 0.5) + 0.5·z(I), where I is the identity axis.
- 639's identity axis beat CLEAN4 on gate signs leave-one-line-out (H1 +0.016, HCT116 +0.015, HEK293T +0.024). That sign-accuracy
  gain does not reach A/B/C's nmae through the hook. H1's real emitter priced idall3 − c4h nmae at −0.00223 (f 0.35); the board
  went the other way.
- This read cannot separate two explanations: the identity term flips the signs of uncalled genes the board scores by magnitude
  (nmae is MAE on the reference's DE genes), or A/B/C's lines (Jurkat, HeLa, CAL33) sit differently on the identity axis from
  the panels 639 used.
- **Tries (11.1c.4):** 1 board arm for the identity hook family. KGID retired without a board read, and ID05 / ID10 (block form)
  never ran.

## Consequences
- **No recipe change.** mde_sign's decide value is `clean4:0.05`, taken after recipe_dr524 resolves (CHECKLIST step 10).
- The queue drops idall3 (it has been read). ID05 / ID10 (the block form) stay empty-slot only. The identity axis now has no
  positive reading on any delivery route: the hook is BELOW on the board, and 609 expects the block route at about +0.000 … +0.001.
- **Instrument note:** for hook-class levers, the board reads below every H1-family surface: 496a by 1.4× against h1-test,
  idall3 by ≥ 0.0006 overall. 618 is the less biased of the two surfaces, so price future hook levers on 618 with a downward
  margin.
- The displayed row is idall3 (0.2746732). The best scored entry is still 498p2g (0.2749190).
