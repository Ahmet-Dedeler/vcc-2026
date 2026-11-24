# 443 — 228's per-gene reliability gain re-priced on the CURRENT bulk: it does not stack with sharpening. CLOSED (2026-09-24, coordinator)

## Question
410d/410d2 (228's per-gene cross-line reliability gain in the bulk, via ballast) was shelved on render energy (context C
blew up ×2.16), and 414 priced a bounded, amplitude-matched variant at +0.0025…+0.0046 overall **on 410c's bulk**. Since
then the bulk gained the stack (H1 ×2 + CD4 ×1.5) and sharpening (q 1.25 → 1.5), and re-levelling now fixes per-context
energy. Does the gain still buy direction on top of them?

## Instrument
`comp443.py`: 420's PLOLO harness (fold build, stack composition, 410 key K3, sharpened render, 414's ballast realisation
at KMAX 16/24/32, β-emulated pds at λ 0.5), in the builder's order: key → gain (non-panel Frobenius preserve) →
sharpened render. The gain is each fold's own 228 artifact (`gain-<fold>.npz`, `gain_reg`, fitted without the held
line), m = 1 + t(g/median − 1), optionally clipped. **Known answers: `stack|K3+sh1.25` and `stack|K3+sh1.5` reproduce
comp420's rows with |d| = 0 on all four folds** (design and every KMAX). Runtime 33 s per fold, peak RSS 5.7 GiB.

## Result (vs `stack|K3+sh1.5`, realised KMAX 16; relative aggregate cosine / Δ β-emulated pds)
| arm | k562_gwps | orion_hct116 | orion_hek293t | kolf2 (board predictor for shape arms) |
|---|---|---|---|---|
| gain t 0.25, unbounded (410d) | −10.42 % / −0.0038 | +0.29 % / −0.0030 | −1.55 % / +0.0029 | −3.71 % / −0.0107 |
| t 0.25, clip [0.5, 2] | −1.33 % / −0.0030 | +1.85 % / −0.0012 | +1.35 % / +0.0058 | +0.07 % / −0.0063 |
| t 0.25, clip [0.67, 1.5] | −0.89 % / −0.0025 | +3.11 % / −0.0001 | +2.08 % / +0.0055 | +0.33 % / −0.0037 |
| t 0.125, clip [0.5, 2] | −0.85 % / −0.0030 | +2.89 % / +0.0007 | +2.50 % / +0.0055 | +0.91 % / −0.0049 |
| t 0.5, clip [0.5, 2] | −1.57 % / −0.0041 | +1.32 % / −0.0028 | +0.29 % / +0.0047 | −1.07 % / −0.0063 |
| CONTROL: clip [0.5,2] permuted over genes | −17.5 % / −0.0191 | −6.4 % / −0.0029 | −8.2 % / −0.0066 | −11.5 % / −0.0102 |
The controls lose everywhere, so the instrument sees the gain's alignment. But its value on this bulk is ~+1 % cosine
(cancer mean), not 414's +3.9 % on the unsharpened bulk, with pds mixed on the cancer folds and negative on kolf2.

## Verdict
**Closed on the current bulk: sharpening absorbed the per-gene reliability reweighting.** Both amplify the entries whose
signal transfers; once q ≥ 1.25 is in, the gain adds ≤ 1 % cosine at a kolf2 pds cost of 0.004–0.006. No artifact, no
builder registration, no slot. (Also found: `VCC_BULK_GENE_GAIN` is applied before `VCC_BULK_LOAD`, and a LOAD replaces
the bulk wholesale, so the two together silently drop the gain; any future gain arm must go into the dump.)
