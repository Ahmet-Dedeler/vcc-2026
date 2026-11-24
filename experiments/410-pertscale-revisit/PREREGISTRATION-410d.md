# 410d — 228's per-gene reliability gain on the BULK only, on 410c; preregistered 2026-09-22 ~19:20Z (before 410b/410c read)

One variable against 410c: `VCC_BULK_GENE_GAIN=g228reg_t025` — m_g = 1 + 0.25 (g/median − 1) from 228's gain_reg (six fit
lines, 1,758 NON-panel targets, production operator; artifact da937cfa…), multiplied into the BULK matrix after the
per-target key and renormalised to the non-panel Frobenius norm; the carriers (DE channel) are untouched.
Expected config diff vs 410c's meta: exactly `bulk_gene_gain` and `bulk_gene_gain_sha256`.

Evidence: 228 (PLOLO, read-back, production composition) pds +0.00180 and mse +0.00227 overall at t 0.25 (its nmae −0.0104
was the weld, removed by ballast); PLOLO stacking with the per-target key is multiplicative (+0.52 relative aggregate cosine
for key uncapped + gain vs +0.34…+0.39 key alone; stack410.py). Caution recorded by 228: the pds gain concentrates in three
targets (leave-top-3-out retention −0.80) and read-back retention of cosine was 0.5 of the effect level.

Predictions vs 410c: raw mse −0.005…−0.020; pds +0.000…+0.004 raw; DE members within noise. Bins (Δ overall vs 410c):
≥ +0.002 → keep the gain; −0.002…+0.002 → flat, drop it (heavier construction for nothing); ≤ −0.002 → the gain does not
transfer through the whole row (it was validated on the whitened base only); drop it.

## ADDENDUM 2026-09-22 20:10Z — 410d as built is withdrawn (NOGO-410d); 410d2 replaces it (before any 410c/410d number)
Worker 414 (`experiments/414-bulk-reprice/URGENT-410d.md`) showed that 410d's EFFECT-space energy normalisation raises the
energy the scorer reads: the render zeroes genes below each context's 10-CPM control gate, 228's gain favours expressed
genes, and the real build logged gain constant 0.5325 (kolf2 analogue 0.545: gated-out energy share 18 % → 6.3 %) — an
amplitude step of ~+14…20 % past context A's mse cliff (ρ_A 1.043 on 410c), predicted −0.003…−0.021 overall. That is an
amplitude arm, not the gain arm registered above.
**410d2** = the same gain artifact (da937cfa) with knob value `g228reg_t025rg`: energy matched on each context's RENDERED
non-panel genes (control CPM > 10), summed over A/B/C. One variable vs 410c (`bulk_gene_gain` + sha). Pre-upload check
(automatic in the build log? no — run by the PI): `experiments/414-bulk-reprice/cliffcheck414.py` on the artifact must
read A_c(410d2)/A_c(410c) within 0.97…1.03 in every context, else NOGO-410d2.
Predictions vs 410c (414's amplitude-matched price): +0.0025…+0.0046 overall; pds +0.000…+0.004 raw; mse −0.005…−0.015
raw; DE members within noise. Bins unchanged from above.
