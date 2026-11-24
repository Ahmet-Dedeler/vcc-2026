# 402 — Preregistration: context A from same-line (Nadig Jurkat) neighbourhood votes (`VCC_DE_SIGN_PRIOR=knn25o400jaq80`) on 399. Frozen 2026-09-16 ~06:00Z, before the build and before any entry exists.

**The arm.** 399's exact configuration; contexts B and C byte-identical to 399 (same rows, same default). Context A's per-target
sign vector is replaced: the same 25 nearest K562-GWPS perturbations per target, but their votes come from the **Nadig Jurkat
bank alone** (2,270 perturbations, 8,881 genes), shrunk toward the Jurkat marginal with prior weight 50; default = Jurkat
marginal sign (CLEAN4 where no vote). One variable: the vote family for context A. Artifact
`experiments/402-ctxA-jurkat-knn/work/knn25o400_A_jurkat50_prior_challenge.npz` (sha `8a7c592c…`; 269 targets with a row; median
26 % of A's block signs differ from 399's). Expected config diff vs 399's meta: exactly `de_sign_prior`, `de_sign_prior_sha256`,
`de_sign_prior_context_A` (absent → nadig_jurkat_knn25_shrink50). Package `submissions/402-cap15-partial50-knn25o400jaq80.vcc`,
model name "the paperclip demanded a raise". **Reference for Δ: 399's entry, whatever its score.** If 399's falsifiers fail, void.

## Evidence before the entry
| instrument | six-bank kNN25 shrink 400 (= 399) | Jurkat-only kNN25 shrink 50 | Δ |
|---|---|---|---|
| 401, Nadig Jurkat split-half (three seeded splits, sd ≤ 0.001), block accuracy | 0.5806 | **0.6065** | **+0.026** |
| 370 (2026-09-10), same-line marginal vs cross-line marginal on the A-block, Jurkat truth | ≈ 0.55 | 0.577 (marginal) | +0.03 |
Context A is Jurkat (041; CCLE match). 319pA (2026-09-10) blended Jurkat-calibrated ROWS into the bulk and lost pds −0.009; this
arm touches the DE channel only (pds/mse pinned by construction), the failure mode that killed 319pA cannot occur here.

## Predictions (board units; fid slope +0.561)
Block precision on A +0.026 → A's fid +0.022 raw → panel fid +0.0074 raw (+0.004 overall) if challenge A behaves like Nadig
Jurkat; realised transfer so far 64 % (392) → **central Δ ≈ +0.003 overall**, range 0 … +0.006. nmae/reach on A expected to
follow fid (weighted accuracy +0.036 in 401). B/C contribute exactly 0.

## Falsifiers
1. `score_pds`, `score_mse` equal 399's to ≤ 1e-9. 2. jac within ±0.0005 raw of 399's.

## Primary reading: Δ score_fid vs 399
| bin | reading |
|---|---|
| ≥ +0.012 (≥ +0.004 raw, i.e. ≥ half the A-only forecast) | same-line public data is a lever through the DE channel; October: use it wherever D/E/F has a same-line bank |
| +0.004 … +0.012 | partial; incumbent if overall ≥ +0.0005 |
| −0.004 … +0.004 | no transfer: challenge A's block signs are not Nadig-Jurkat-like beyond what the six banks carry |
| ≤ −0.004 | harm; close same-line votes for A |

## Overall bins (Δ score_avg vs 399)
| bin | reading |
|---|---|
| ≥ +0.002 | transfers; incumbent |
| +0.0005 … +0.002 | partial; incumbent if fid in its top two bins |
| −0.0005 … +0.0005 | flat |
| ≤ −0.0005 | harm; 399 stays |

**Stated in advance.** Displayed row not corrected. Upload ~12 min. First slot of 2026-09-17 UTC unless a slot opens today.
