# 395 — Preregistration: q80-protected CPM-block signs (`VCC_DE_SIGN_PRIOR=clean4q80`) on 392. Frozen 2026-09-16 ~05:40Z, before the build finishes and before any entry exists.

**The arm.** 392's exact configuration (`VCC_NORMREL_EXP=-1.5`, `VCC_NORMREL_MODE=cap`, `VCC_CD4_PARTIAL_G=0.5`,
`VCC_DE_SIGN_PRIOR=clean4`) with the knob value `clean4q80` instead of `clean4`: the strongest 20 % of each target's
called CPM-block genes (by |E0|, the delivered log2FC) keep their delivered sign and magnitude; the other 80 % take the
CLEAN4 prior sign at the delivered magnitude with 390 C5b's sign-preserving mass balance, exactly as in 392. Builder
`05e49916…` (commit 7f1cddd6; knob inert unset), emitter `src/vcc2026/rc_emit_k.py` `6ff0b866…` — its diff against 392's
pinned `a12cf8f9…` is the five `protect_q` lines only, so the `clean4` path is byte-identical and this is ONE variable
(`protect_q` 0.8). Prior file unchanged (`7e532f84…`). Expected config diff vs 392's meta: exactly `de_sign_protect_q`
(absent → 0.8) and `de_sign_emitter_sha256` (a12cf8f9… → 6ff0b866…). Package
`submissions/395-cap15-partial50-designclean4q80.vcc`, model name "the thermostat has a union rep".
**Reference for Δ: 392, `TJTOj9qrbiobJBpqLSa3`, 0.2226221** (raw pds 0.854738 / mse 0.809898 / nmae 0.938837 /
fid 0.516157 / reach 0.273878 / jac 0.040898). First slot of 2026-09-16 UTC.

## Why this arm, and what it should teach
392's fid transferred (+0.0082 raw, top bin) but nmae worsened +0.0041 raw and reach fell −0.0039 raw, where H1 had
predicted +0.0024 (unresolved) and +0.0175 (unresolved). Both members read the strongest block genes: nmae is summed over
the reference's own significant set, reach walks the front of our order. The prior's block precision (0.52) barely beats the
delivered sign's, and a wrong flip on a strong gene costs 2|lfc| on nmae and a purity break on reach. Protecting the strongest
20 % (394) is the priced construction that removes exactly those flips: on H1 vs 392's construction (direct paired contrast)
fid −0.0008 [−0.0028, +0.0012], reach +0.0087 [−0.005, +0.025], nmae −0.0012 [−0.0042, +0.0013] (lower better); vs cap15 it is
the only variant resolved positive on H1-local overall (+0.0071 [+0.0009, +0.0139]).

## Predictions (board units; slopes per raw unit from the 2026-09-13 regression: fid +0.561, reach +0.188, nmae −0.255)
Central: fid −0.0008 raw (−0.0004), reach +0.006 raw (+0.0011), nmae −0.0015 raw (+0.0004), pds/mse/jac 0 → **Δ ≈ +0.0011
overall**; plausible range −0.001 … +0.003.

## Falsifiers (checked first; a failure voids the reading as a one-variable test)
1. `score_pds` and `score_mse` equal 392's to ≤ 1e-9 (per-gene totals identical by construction).
2. `de_wilcoxon_sig_jaccard` within ±0.0005 raw of 392's (call set unchanged by construction).

## Primary reading: the two members 392 lost — Δnmae raw and Δreach raw vs 392
| bin | reading |
|---|---|
| nmae ≤ −0.002 raw AND reach ≥ +0.003 raw | the strong-gene flips were the cost; protection is the right construction; dose 2 = q50 (390 §6: fid +0.0081, reach +0.025, nmae −0.0001 on H1) |
| exactly one of the two recovers by those margins | partial; read which member; q50 only if reach carried it |
| neither (|Δnmae| < 0.002 and Δreach < +0.003) | the nmae/reach shortfall of 392 is not the strong block genes; the protection axis is closed at q80 |
| nmae ≥ +0.002 raw or reach ≤ −0.003 raw (worse) | protection hurts: the delivered signs on the strong block genes are worse than the prior's on A/B/C; 392 stays |

## Overall bins (Δ score_avg vs 0.2226221)
| bin | reading |
|---|---|
| **≥ +0.0015** | protection pays; 395 becomes the incumbent; next one-variable step is the q50 dose |
| +0.0003 … +0.0015 | small gain; incumbent if Δfid raw ≥ −0.002 (the fid transfer is not being given back); dose step only if reach carried the gain |
| −0.0003 … +0.0003 | flat; 392 and 395 indistinguishable; 392 stays the incumbent by its fid; protection axis closed |
| ≤ −0.0003 | harm; 392 stays; do not try q50 |

**Stated in advance.** The displayed row becomes this entry's score and is not corrected. Upload ~12 min. Slot 2 of today is
decided after this reads (candidate: a Jurkat-only per-gene tilt for context A, B/C unchanged — being sized now).

## Addendum 2026-09-16 05:50Z — build note, written after the build and before the entry exists; bins unchanged
The build's config diff vs 392's meta was exactly three keys: `de_sign_prior` (clean4 → clean4q80), `de_sign_protect_q` (absent →
0.8), `de_sign_emitter_sha256` (a12cf8f9 → 6ff0b866). The registration above listed two, forgetting that the knob's own value is
recorded; `run395.sh`'s check therefore stopped before packaging and `finish395.sh` completed prep/novelty/sha with the corrected
three-key expectation. Package sha256 `19435f09…`, uploaded 05:50:21Z as "the thermostat has a union rep". Build diagnostics
(`the-thermostat-has-a-union-rep.h5ad.meta.json`): see BOARD-RESULT.md once read.
