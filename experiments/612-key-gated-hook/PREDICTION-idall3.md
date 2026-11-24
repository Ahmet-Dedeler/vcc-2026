# 612: pre-registered prediction for the UNGATED v3 identity hook (`idall3`) on H1's real emitter

Written 2026-10-02 ~13:47Z by [21dd4c], BEFORE run612's `idall3` mode exists and before any idall3 run.
- Form: VCC_DE_MDE_SIGN = 604 v3's identity sign (score = z(P_C4 − ½) + 0.5·z(I) on CLEAN4-covered genes, τ 0) on every
  uncalled gene. NO key gate.
- Source of the prediction: 613d's add-rule instrument (613/block613d5.py, descriptive run at 13:46Z). Its known-answer
  rows match H1's real emitter within 0.0004 (c4h − base, c4kg − c4h).

## Prediction (raw nmae, idall3 − c4h, mean of seeds 0–2 by score612)
| f 1 | f 0.5 | f 0.35 |
|---|---|---|
| −0.00160 | −0.00232 | −0.00225 |
- PASS iff every regime lies within ±0.0006 of its prediction (A3's tolerance: the largest known error 0.0004 + 2 SE):
  - f 1 ∈ [−0.0022, −0.0010];
  - f 0.5 ∈ [−0.0029, −0.0017];
  - f 0.35 ∈ [−0.0029, −0.0017].
- Reported beside: idall3 − idkg3 (613d predicts +0.0004 / −0.0020 / −0.0023) and idall3 − base.

## What each outcome means
- **PASS:** 613d gains a second out-of-sample retrodiction (a new form, on the real emitter).
  - Its rows then say the ungated v3 hook beats c4h and KGID on H1, HCT116 and HEK293T at f 0.5 / 0.35 by 0.001–0.002
    raw.
  - It loses only on the CLEAN4-weakened K562 fold (+0.0016 / +0.0024 vs c4h). INSTRUMENTS.md reading rule 1 treats K562
    as a stress test.
  - The ungated v3 hook becomes a candidate October hook form, in place of the plain CLEAN4 hook. It needs its own
    board prereg (vs 496a's plain hook) with the K562 tail stated.
- **FAIL:** 613d's rows for new hook forms are not trusted beyond its known-answer contrasts.
- Try count: hook forms priced on 613d today: c4h, c4kg, KGID (v3), ungated v3 (4).
