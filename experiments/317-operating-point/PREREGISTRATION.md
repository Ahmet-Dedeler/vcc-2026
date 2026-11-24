# 317 — PREREGISTRATION (committed before any build or score)

**Question.** Does OUR pipeline have an operating point that trades `mse` for `nmae` + `reach`
and beats the shipped composite? The operating point is the variable; the model is fixed.

**Instrument.** The real `cell-eval2` 0.16.0 (`--preset vcc2026`, DE backend `scanpy`) on
`data/harness/h1-test.h5ad`, run on disposable EC2 workers via `scripts/cloud_score.sh`
(≈ $1.48 and ≈ 2.5 h per six-member score, zero local memory). Deliverable is the RAW
`agg_results.csv` six-member row per operating point; raw deltas are converted to board
overall points through board-fitted per-member spans (094's method, refit here on the
2026-09-07 snapshots) and, where one exists, the arm-class transfer factor.

**Known-answer rows, frozen now.**
- **KA-1 (instrument).** Re-score 221's `rc_callset.h5ad`
  (`s3://<your-s3-bucket>/input/cloud-score/ladder8-20260904/preds/rc_callset.h5ad`)
  and reproduce 221's recorded `agg_results.csv` row. ABORT if any statistic moves.
- **KA-2 (build).** The rebuilt `base` must have sha256
  `048b9286030ff0124514134b55d935fdb79a4d002cc04d590d83676b32553c77` — 221's `rc_callset`,
  byte for byte (gate 5). ABORT otherwise.
- **KA-3 (operator path).** 221's gates 1–3 (a15 + a19 operator anchors to 1e-11 relative;
  170's h1 effect to 1e-12; the CD4 f-linearity to 1e-12) are kept verbatim.
- **KA-4 (board pair, read not run).** The `--overdispersion 4.0` knob has its own board
  point: `a2-prune` → **`a8`** (`docs/BOARD-LEDGER.md`), scaled `reach` +0.021398,
  `fid` +0.006393, `jac` −0.003439, `nmae` −0.000038, `pds` 0.000000, `mse` 0.000000,
  overall **+0.0040523**. The amplitude knob has one too: `a8` → **`a11`**
  (`VCC_ALPHA` 1.0796 → 0.838), `mse` +0.082041, `reach` −0.012129, overall **+0.0117543**.
  Signs measured here must be read against those.
- **Noise floor.** 221's `a15` vs `a15_seed1` re-seeded twin on the same surface.

**Frozen structural prediction (gate 4, decided before the build).** `emit_counts(mode="quasi")`
repairs every gene's panel total exactly, so every EMISSION arm must be panel-total
digest-identical to `base` on 85/85 blocks and every AMPLITUDE arm must differ on 85/85. If
that holds, `pds` and `mse` cannot move on the emission axis and the only knob in this sweep
that can pay the `mse` side of the trade is amplitude.

**Arms — one variable each off the shipped point** (`--tail-spread 1.0 --keep-top 8000
--keep-spread 0.25 --overdispersion 4.0`, `VCC_ALPHA 0.982`):

| arm | knob | value |
|---|---|---|
| `base` | — | shipped (= 221 `rc_callset`) |
| `od1` | `--overdispersion` | 4.0 → 1.0 |
| `ks010` | `--keep-spread` | 0.25 → 0.10 |
| `ts070` | `--tail-spread` | 1.0 → 0.70 |
| `kt16000` | `--keep-top` | 8000 → 16000 |
| `amp070` | delivered amplitude | ×0.70 (α 0.982 → 0.687) |
| `amp085` | delivered amplitude | ×0.85 (α 0.982 → 0.835) |
| `amp120` | delivered amplitude | ×1.20 (α 0.982 → 1.178) |
| `amp145` | delivered amplitude | ×1.45 (α 0.982 → 1.424) |

**Decision rule (frozen).** 317 says a trade EXISTS only if some arm's board-converted composite
beats `base` by more than the emitter noise floor AND the gain is carried by `nmae` + `reach`
against an `mse` loss. Flat or falling everywhere = NO-GO, which re-establishes §4.1c.

**Budget.** ≤ $25 of AWS. 9 on-demand `r7i.2xlarge`, `--max-minutes` capped so no instance can
reach the (unbuilt, ~4 h) `prep-real-bundle` stage; every instance terminated on result.

**Nothing here is competition-facing.** No submission is built, nothing is selected or tuned,
`submissions/` is not touched.
