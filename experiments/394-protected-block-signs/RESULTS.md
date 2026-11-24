# 394 — q80-protected CPM-block signs end to end on H1: resolved positive vs cap15 (+0.0071 H1-local), not distinguishable from 391's full-block flip

Status: **COMPLETE** (PI-run, 2026-09-13 07:35–07:40Z).

## Known answer
`run394.py cap15` reproduced 323's cap15 exactly (six raw members, both digests, mse legs; the driver aborts otherwise).
Comparator `k3c5b` is 391's run, copied (deterministic harness: same picks and emission stream; the report checks the
statistics sha).

## Arm
`k3c5bq80`: 391's k3c5b (CLEAN4 prior signs on the called CPM block at the delivered magnitude, down-signed magnitudes
shrunk for zero block mass change) except that only block genes with |E0| ≤ the block's 80th percentile take the prior
sign; the strongest 20 % keep their delivered LFC. Emitter `k_emit394.py` = 391's + `protect_q`.

## Results (raw, paired 10,000-target bootstrap)
| member | k3c5bq80 − cap15 | k3c5b − cap15 (391) | k3c5bq80 − k3c5b (direct) |
|---|---|---|---|
| pds / mse | 0 / 0 exactly | 0 / 0 exactly | identical |
| fid | **+0.0115 [+0.0070, +0.0162]** | +0.0122 [+0.0067, +0.0182] | −0.0008 [−0.0028, +0.0012] |
| reach | +0.026 [−0.004, +0.060] | +0.017 [−0.019, +0.055] | +0.0087 [−0.005, +0.025] |
| nmae (lower better) | +0.0013 [−0.0009, +0.0035] | +0.0024 [−0.0011, +0.0066] | −0.0012 [−0.0042, +0.0013] |
| jac | −0.00003 [−0.00006, −0.00000] | −0.00001 | −0.00002 |
| **H1-local overall** | **+0.0071 [+0.0009, +0.0139]** | +0.0052 [−0.0029, +0.0132] | — |
Board-slope estimate vs cap15: +0.0110 (391's k3c5b: +0.0095). The analytic prediction (390 §6) was +0.0105.

## Verdict
The q80-protected construction is resolved positive on H1 where the full-block flip was not, with a smaller
intervention (≈ 40 % of block calls flipped vs ≈ 50 %). No member separates it from the full-block flip. It is the
preferred construction for the next board step **if 392's fid transfers**; it needs a builder knob value (e.g.
`VCC_DE_SIGN_PRIOR=clean4q80`) added after 392's build releases the builder (the builder re-hashes the pinned emitter per
context), with the unset path unchanged.

## Limits
H1 development panel; the direct contrast is unresolved on every member.
