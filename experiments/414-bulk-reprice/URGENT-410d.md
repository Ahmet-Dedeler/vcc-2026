# URGENT (414, 2026-09-22 20:30Z): DO NOT UPLOAD 410d as built. Its artifact is measured, and context C's bulk energy doubled

Artifact: `experiments/410-pertscale-revisit/work/build_d/the-lighthouse-keeper-went-remote.h5ad`, built 23:24 local.
It was streamed read-only (`cliffcheck414.py`, `diag410d.py`).

| context | render-space energy A_c, 410d / 410c | rho_c = B_c/claim_c (410c → 410d) | amplitude-only Δraw mse vs 410b |
|---|---:|---|---:|
| A | **×1.197** | 1.043 → 1.246 | +0.221 |
| B | **×1.182** | 0.971 → 1.147 | +0.155 |
| C | **×2.127** | 0.930 → 1.978 | **+1.124** |

**What happened.**
- In A and B, the shipped gain m moves energy from genes the render drops (ungated) onto gated genes. The rise is
  +18–20 %, exactly what PLOLO's kolf2 fold predicted (+19.1/+19.3/+19.9 %; kolf2's multiplier IS the shipped artifact,
  max |d| = 0).
- In C, one target explodes. **SYVN1's bulk row energy goes from 229 to 4,858 (×21).** m ≈ 4–5 on the unfolded-protein-response
  genes (HSPA5 m 4.14, HSP90B1 4.13, SDF2L1 5.32, MANF 5.08, HERPUD1 4.79, PDIA4 4.72). HSP90B1 sits at 668 CPM in
  context C, and its delivered log effect goes from 3.55 to 6.60 (HSPA5 4.07 → 7.75), so it takes a large share of the
  block's mass. The render's mass balance then pushes every other gene down, and that row's energy multiplies.
- m also ranges from −6.7 to 19.2, and flips the sign of 116 genes.

**Consequence (187's per-context model, exact scorer arithmetic past the cliff).** All three contexts end up over their
cliffs, and context C's raw mse rises by about +1.1 from amplitude alone. SYVN1-C alone adds about 2,800 of squared
error against the context's ΣD ≈ 4,246. score_mse would clip to about 0 in C and fall in A/B: **roughly −0.03 to −0.04
overall** vs 410c, whatever the gain's direction gain is. The preregistration expected raw mse −0.005 to −0.020.

**410d2 (render-gated energy matching, `patch_gain410d2.py`) probably does NOT remove the C blow-up.**
- The blow-up comes from m's relative boost on the UPR genes, not from the global constant. 410d applied 4.13 × 0.5325 =
  ×2.20 to HSP90B1. 410d2's constant emulates to c_rg = 0.4846 (`emul410d2.py`), which still gives ×2.00.
- **The emulation cannot settle context C.** It uses 410c's RENDERED bulk as a stand-in for the effect matrix. Its
  known-answer row reproduces the measured 410d energy in A and B (emulated/measured 0.974, 0.993) but FAILS in C (0.554).
  The render's mass fix compresses SYVN1's true effect, so the stand-in understates it.
- So 410d2 is predicted to pass its own 0.97…1.03 gate in A and B (×0.971/×0.978), with C unknown and at risk.
- A safer variant clips m to [1/2, 2] before the render-gated normalisation. Then c = 0.616, HSP90B1's boost is ×1.23,
  and the emulated energy is ×0.994/×1.002/×0.995 with no row above ×1.1. Direction cost unmeasured.
- Either way, **run cliffcheck414 and `diag410d.py <410c.h5ad> <410d2.h5ad>` on the built artifact.** GO only if every
  context's energy ratio is in 0.97…1.03 AND the max per-target energy ratio is ≲ 1.5.

**The fix is a builder change, not mine to make.** Amplitude-match the gain per context (predicted +0.0025 to +0.0046
overall) AND bound m (or each row's render-space energy). Or drop the gain; its direction gain is +4 % cos.

**Pre-upload checks for ANY bulk arm** (10 s each, < 3 GB RSS, read-only):

    uv run python experiments/414-bulk-reprice/cliffcheck414.py <artifact.h5ad>     # render-space energy and rho_c per context
    uv run python experiments/414-bulk-reprice/diag410d.py <ref.h5ad> <artifact.h5ad>  # per-gene / per-target energy vs a reference

cliffcheck414 is validated on 410c: measured ×1.0099/1.0084/1.0149 vs 410b, against PLOLO's predicted ×1.013/1.011/1.012.
