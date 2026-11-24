# 453 — preregistration (written 2026-09-26 00:25Z, before any 453 build or entry exists)

## The arm
The October path (`scripts/def_day.sh`, HEAD `b517b6d2`, unmodified) on 180's `def-bundle` (A/B/C controls and panel as
D/E/F), relabelled back to A/B/C (`relabel453.py`). q and the same-line flags come from the 09-26 reads (`work/args.env`,
written after them).

**Counterpart** (the tuned A/B/C arm with the same q and the same DE settings):
- 446b if the flags keep `noprog`;
- 446bja if `noprog` is dropped. The same-line stage then gives D (= A, Jurkat) exactly 446bja's programme: 451's KA-2
  rebuilt 446bja byte for byte from the same artifact type.

## What differs from the counterpart, and its price (overall, board units)
| component | counterpart | 453 | price | source |
|---|---|---|---:|---|
| bulk stack | H1 ×2 + CD4 ×1.5 (ABC sizing) | H1 ×2 + CD4 ×1 (DEF carrying sizing) | **−0.0017** (k 1) … **−0.0022** (k 1.3, 417e's realised stack rate) | 438 §3.1 "CD4 ×1.5 on top of H1 ×2" on abc +0.00170, on abcdef +0.00172 |
| DEF vs ABC composition rules | ABC | DEF (carrying g 0.344, V5_abc) | ~0 (≤ 0.00014) | 438 §4.4 |
| re-level | ρ matched to 433k's cliffs (0.9908 / 0.9827 / 0.9873) | NTC-proxy ρ* 0.9793 per context, secant loop | about −0.0002 (cliff slope 0.31 raw mse per unit scale below the cliff) | ntclevel417 |
| operator, key, CD4 banks | pinned A/B/C artifacts | re-fitted on the day | ~0 (452: ρ agrees to 4 dp across refits) | 452 |
| carriers (DE channel) | ABC CD4 sizing f 0.3461 | carrying f 0.3432 (0.8 %), fresh operator | ~0, within emission noise | 438 §4.4 |

**Prediction: overall Δ −0.0024 (band −0.0045 … −0.0005).**
- It sits almost entirely in the bulk: pds Δ ≈ −0.001 … −0.003 raw, mse Δ ≈ +0.003 … +0.010 raw.
- The DE members are pinned up to emission noise.

## Bins (453 − counterpart)
- **Bulk channel** (pds and mse, scaled, summed in overall units; these reproduce exactly on the board):
  - −0.0045 … −0.0005: the October path is validated end to end. The CD4 half is priced right.
  - Below −0.0060: a day-of bulk component is mis-priced. Attribute it by local rebuilds (no slot): swap back the A/B/C
    operator, then the pinned key, then 446b's factors.
  - Above +0.0005: the October stack is not worse on A/B/C. Re-read 438's CD4 half (a finding either way).
- **DE channel** (fid, reach, jac, nmae in overall units; unwelded, so only emission noise, sd ≈ 0.0045):
  - within ±0.009: pinned. The DEF carrier path equals the tuned one.
  - Outside: diagnose the carriers (CD4 sizing, operator) before October.
- **Hard failures** stop the October plan until explained:
  - any context's ρ above 1 in the final build (def_day refuses these);
  - an overall below −0.02.

## Name
A fresh joke, chosen at arming. Description empty.

## 453 verified 02:35Z 09-26 (still before any 09-27 or 09-28 entry)
Package `submissions/453dr-october-path-abc-q1.75.vcc`, sha256 `a76283bd…`, novelty rc 0. `gate453.py` passes every
gate and wrote `work/gates453.log`, which `456/chain0928.sh` requires for slot 1.
- Config diff vs 446b is exactly the preregistered set: `bulk_ctx_scale`, `bulk_load` (+sha), `cd4_deflate_basis` /
  `cd4_deflate_k`, `cd4_sizing`, `cd4_target_g`, `contexts` / `expect_contexts`, `pert_scale` (the id only) and
  `sameline`. `pert_scale_sha256` is identical: the day key is byte-identical to 446b's.
- Converged after 2 levelled builds: ρ 0.9793 / 0.9794 / 0.9795 at A 1.0859 / B 1.1383 / C 1.1662.
- Zero-gene audit 0 (17 cells below −5). Relabel D/E/F → A/B/C at 120,000 cells per context.
- The same-line stage was a no-op as designed: A = Jurkat, 0 rows changed, programme off, w 2.
