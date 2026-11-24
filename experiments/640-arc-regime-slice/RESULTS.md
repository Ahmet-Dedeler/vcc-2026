# 640 — Arc's published selection rules vs our October-like panels: 438's CD4 verdict holds

Date 2026-10-02, session [537757] (Ahmet's Arc-blog follow-up thread). Analysis only: no build, no `vcc`, no
`cell-eval2`, no download. It reads 438's saved per-target arrays and 204's residual tables. Peak RSS 135 MB.

**Class:** instrument validity (rule 3, "validate the validator"). No recipe change.

## Why

Arc's post of 2026-10-01 ([Behind the Data 2026](https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge-2026))
gives the selection rules for the D/E/F panel:

- every target has ≥ 80 % median on-target knockdown in every line (Fig. 2: VCC median 89 %);
- 400 cells are subsampled from a median of 500 per perturbation;
- the same sgRNA sequences are used in all six lines. The 2025 library used Replogle's K562 day-8 library protospacers
  (301).

Two more selection rules are in FINDINGS H.1: targets are binned by H1-pilot DE count (negligible / subtle / strong) and
sampled for diversity.

438's October-like draws (`real`, `efp`, `efr`) are uniform random draws from public target pools, scored against public
truths: K562 GWPS, Orion HCT116 and Orion HEK293T, each held out. If the October verdicts depended on targets that Arc's
rules would exclude, the instrument would mis-price D/E/F. The verdict at stake is CD4 ×1 vs ×1.5: ×1.5 won the A/B/C
board on the October path (457dr, +0.0031), but October ships ×1 on 438's draws.

## Results

**1. The draws do carry ~2× more weak-knockdown targets** (`kd_premise640.py`). Residual is the target's own
expression in the held-out truth, perturbed / control (204).

| fold | panel | median residual | residual ≤ 0.2 | residual > 0.5 |
|---|---|---:|---:|---:|
| K562 GWPS | abc / real / efp / efr | 0.135 / 0.137 / 0.140 / 0.167 | 0.73 / 0.68 / 0.69 / 0.59 | 0.04 / 0.07 / 0.08 / 0.14 |
| Orion HCT116 | abc / real / efp / efr | 0.294 / 0.313 / 0.254 / 0.315 | 0.19 / 0.22 / 0.39 / 0.21 | 0.10 / 0.19 / 0.22 / 0.24 |
| Orion HEK293T | abc / real / efp / efr | 0.505 / 0.509 / 0.480 / 0.525 | 0.02 / 0.03 / 0.23 / 0.07 | 0.51 / 0.53 / 0.47 / 0.56 |

The HEK293T truth sits at ~50 % pseudobulk knockdown on every panel. Nothing in that fold is in Arc's regime.

**2. They are only mildly weaker in response energy** (`strength_premise640.py`). In abc's terciles of truth energy d,
the draws put 0.39–0.44 of their targets in the weak tercile on K562 (abc 0.33), 0.30–0.44 on HCT116 and 0.34–0.41 on
HEK293T.

**3. CD4 ×1.5 on top of H1 ×2 loses on the draws under every Arc-matched subset** (`slice640.py`,
`strength_slice640.py`). Levelled, k = 1, cancer-fold mean, n positive folds of 3.

| subset | abc | real | efp | efr |
|---|---|---|---|---|
| all targets (438) | +0.00170 3/3 | −0.00160 1/3 | −0.00111 0/3 | −0.00236 0/3 |
| truth knockdown, stronger half | +0.00139 2/3 | −0.00182 1/3 | −0.00055 1/3 | −0.00218 0/3 |
| truth residual ≤ 0.5 (no weak tail) | +0.00166 2/3 | −0.00147 1/3 | −0.00155 0/3 | −0.00231 0/3 |
| K562-guide strong (K562 residual ≤ 0.2) | +0.00110 2/3 | −0.00033 2/3 | −0.00124 0/3 | −0.00157 0/3 |
| re-weighted to abc's strength mix | +0.00170 3/3 | −0.00157 1/3 | −0.00101 0/3 | −0.00227 0/3 |

- On abc, CD4's gain sits in the WEAK energy tercile (+0.0049). On the draws, that same tercile is the most negative
  (−0.0011 … −0.0032).
- So the panel effect belongs to the specific targets in each panel. Neither knockdown nor strength mix explains it.
- On the most Arc-like cell (the K562 fold, K562-guide-strong targets), the three draws average +0.0002 against abc's
  +0.0034. That is still a panel effect.

**4. Side readings (438's older K3 / sh1.25 designs; not re-priced on the current recipe):**
- H1's value concentrates on strong-knockdown targets. The H1 ×2 half on abc is +0.0033 on the stronger-knockdown half
  and +0.0014 on the weaker half.
- H1 ×2 → ×3 on the stronger-knockdown half: abc +0.0009 3/3, real +0.0044, efp +0.0032, efr +0.0001.
- Every other step keeps its sign on every subset: q 1.25 → 1.5 and the LODO key are positive everywhere.

## Verdict

**438's October CD4 ×1 decision survives Arc's published selection rules.** Neither the knockdown rule nor the strength
stratification turns the draws' CD4 loss into a gain. No October change.

One sentence on what changed: a candidate mechanism for "the October-like panels mis-price D/E/F" (knockdown-regime or
strength-mix mismatch) was tested and is not the reason CD4 ×1.5 loses off-board.

## Reproduce

```bash
uv run python experiments/640-arc-regime-slice/kd_premise640.py        # table 1
uv run python experiments/640-arc-regime-slice/strength_premise640.py  # result 2
uv run python experiments/640-arc-regime-slice/slice640.py             # knockdown subsets, all 438 steps
uv run python experiments/640-arc-regime-slice/strength_slice640.py    # strength terciles + re-weighting
```

- Known answer: `strength_slice640.py` asserts that unit weights reproduce 438's `Panel.price` to 1e-12, on every step
  and panel (passed).
- Inputs:
  - 438 `work/{abc,real,efp,efr}/arrays438-*-held.npz`; abc K562 sha256 `b7e03f15…`, HCT116 `e33d12d6…`;
  - 204 `work/resid-*.npy`.
- Outputs in `work/`.
