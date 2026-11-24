# 375 — Preregistration: context-A same-line calibration rider on 319p (`VCC_CTX_A_PRIOR`, blend_rows, β=0.35). Frozen 2026-09-10 01:00Z before the build.

**The arm.** 319p's exact configuration (cap −1.5, partial rows 0.34) plus `VCC_CTX_A_PRIOR` pointing at
`experiments/375-ctx-a-arm/work/ctx_a_blend_rows_b035.npz` (SHA-256 `500c24ab1f7c0cc6f81e7a56f1ff361415ac845f76efc17aa2c7e9854627a8bd`):
for context A only, the 259 of 300 panel targets with a three-donor consensus get their design row blended, on the
4,951-gene covered axis, toward the Nadig-Jurkat-calibrated prediction (372: ridge map donor consensus → Jurkat
response fitted on 2,006 non-panel Jurkat targets; 346's estimator verbatim) at β=0.35 with fixed row norm (373's
`blend_rows`, 346's rule). B and C are byte-identical to 319p; genes outside the covered axis and the 41 uncovered
targets are untouched. Builder = HEAD `1318cd3c…` + the 366/373 chain patch, SHA-256 `c6f366cad3d549129df0a7e1789a43264fab2fb7a43b861d677324545cd2fe4a`, unit-tested 100/100 and
byte-inert when the flag is unset; driver `27edc486…`. Config must differ from 319p in exactly
{ctx_a_prior_mode, ctx_a_prior_sha256, ctx_a_prior_context, ctx_a_prior_beta, ctx_a_prior_targets} (the run script stops before prep otherwise).
Package `submissions/319pA-cap15-partial-ctxA035.vcc`. Intended as the second 2026-09-10 slot.

**Reference for Δ:** 319p, `SvzrTOQJAAaNjawzSyrC`, **0.2198189778** (same panel, anchors r4).

**Why β=0.35.** 346 (H1, six members): sign accuracy saturates by β≈0.6 while pds keeps falling; knee at 0.25
(+0.0065 [−0.0002, +0.0124] H1-local, five members beneficial), 0.5 gives +0.0071 [−0.004, +0.016]. 374 (Jurkat's own
frame, retrieval ≡ pds at slope 1.01 per 346): Jurkat pays NO retrieval loss up to β=0.5 (+0.005 at 0.25/0.35,
resolved), keeping 64% of the sign gain at 0.35. 0.35 sits inside H1's positive region and Jurkat's no-loss region.

**Expected size, stated for the record and not as a forecast.** If H1's β=0.25 six-member result transfers to A and
the other two contexts are untouched: ≈ +0.002 overall. Smaller than 319p's +0.0012? No — comparable; it is an
information buy about whether same-line public data transfers to the board (the October playbook), not a jump.

## Bins (Δ overall vs 0.2198190)
| bin | reading |
|---|---|
| **Δ ≥ +0.0005** | same-line calibration transfers to the board. 319pA is the incumbent; next point β=0.5 (the end of Jurkat's no-loss region); the October playbook (per-context calibration from any public same-line bank) is confirmed. |
| **−0.0005 … +0.0005** | no net transfer at 0.35. Read the secondaries: if the DE members rose and pds fell by about as much, the exchange exists on the board at par and a β=0.2 point is earned; if nothing moved, the covered 26.7% of the axis is too small to matter and the family is parked. 319p stays. |
| **Δ ≤ −0.0005** | the board's pds price exceeds the frame's prediction; 319p stays; closed at this construction. |

**Secondaries, no verdict:** pds (frame predicts ≥ 0 on A; a resolved fall is the important negative), mse, nmae,
reach, fid (expected carriers, as on H1), jac (expected flat).

**Stated in advance.** The displayed row becomes this entry's score; not corrected by re-submission. The production
builder file is patched for this build (first production code change of the night); the patch is byte-inert when the
flag is unset (unit-proven) and stays in place afterwards so the arm is reproducible; `def_postbuild`'s config-drift
line will list the five new keys. Upload 13–17 min.
