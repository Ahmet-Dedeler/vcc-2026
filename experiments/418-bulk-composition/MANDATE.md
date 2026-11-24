# 418: re-tune the BULK composition now that ballast unwelds it (mandate from the PI, 2026-09-22 20:45Z)

**Current state.** Read `docs/FINDINGS.md` §0–§2, `experiments/414-bulk-reprice/RESULTS.md` (all of it) and
`experiments/415-lodo-key/RESULTS.md` §1–§4. Since 410b (board 0.2365722, best), the ballast cells pin each block's
pseudobulk (pds, mse) to a BULK design, while the 400−K carrier cells alone set the DE table (fid, reach, jac, nmae).
The DE members sat within 0.0007 raw of 395 on the board.

Every composition constant in the shipped effect was tuned while bulk and DE were welded, with six-member prices where
the DE members moved with the bulk. Examples:
- the CD4 residual fraction (0.346 effective; 212/216/217: pds turns over near f ≈ 0.20 while mse and the DE signs kept
  improving);
- the H1-union residual fraction (0.195 effective);
- CD4 SVD-5 deflation;
- CD4 partial rows g 0.5;
- normrel cap −1.5.

With ballast, the BULK only has to serve pds + mse, so its optimum may sit elsewhere. The carriers can keep theirs.
414 found sharpening this way (169 had been closed on a weld cost). Your job is the composition.

**Territory.** Others: 416 ranks DE-carrier designs on H1, and 417 builds bulk sharpening arms. Do not do their work and
do not touch the carriers. Price ONE-variable changes of the BULK composition with 414's instruments on 415's
reconstruction of the SHIPPED composition ("PP": `prep415.py`, `build_weights415pp.py`; KA-PP reproduces the 410b meta
constants to 8.9e-6):
- (a) CD4 residual fraction × {0.5, 0.75, 1.25, 1.5, 2.0} of shipped;
- (b) H1-union fraction × {0, 0.5, 1.5, 2.0};
- (c) CD4 deflation k ∈ {0, 3, 8} (shipped 5);
- (d) CD4 partial rows g ∈ {0, 0.25, 0.75, 1.0} (shipped 0.5);
- (e) normrel exponent ∈ {0, −0.75, −3} (shipped −1.5, cap);
- (f) at most two ideas of your own that a unwelded bulk makes newly sensible (justify from FINDINGS §1–§2).

Hold the shipped per-target key fixed: 410's `pert_scale_410_lo0.3_hi3.npz` weights, applied as the builder's
`apply_pert_scale` does, which is what a builder arm would do. Price each change twice: on the 410b bulk (no sharpening),
and with 417a's bulk sharpening (sign|h|^1.25 on each frame's 10-CPM gate, per-row gated non-panel norm restored, panel
columns untouched; `experiments/417-bulk-sharpen/patch_sharpen417.py` has the exact transform). Report a
LODO-recomputed key as a secondary column only if cheap.

**Instruments (414's, unchanged):**
- the a23 read-back in the challenge frames: rel aggregate cosine, pds, render-space A_c/B_c;
- the cliff-aware mse route: 187's model, with A_c/B_c ratios on 410b's measured values and X at 410b's board
  calibration;
- the board-emulated pds (β at λ = 0.5);
- predicted Δoverall at board slopes (pds +2.2138, mse −1.051 per raw unit, /6).

Carry k_key = 1.407 only for key-class arms. For composition arms report k = 1 and k = 1.407 as a range.

**Known answers first** (abort on failure):
- 414's `base`, `Kinf`, `Kinf+sh@1.25` rows on 228's composition, from `arms414-<fold>.json`;
- 415's PP rows, from `lodo415.json` / `sens415.json`;
- the unmodified PP composition must read identically through your code path and 415's.

Cancer folds decide. kolf2 is reported but carries no vote.

**Deliverable:** `experiments/418-bulk-composition/RESULTS.md`, containing:
- the ranked table (per fold, plus 3-fold mean and range, with a sign count);
- the two best one-variable bulk changes, with predicted Δoverall vs 417a's design;
- for each, the builder env knob that controls it today. Note that those knobs act on BOTH carriers and bulk, so a
  bulk-only arm needs a builder change (not yours). State which knobs are pure scale-like and which change direction.

**Constraints.**
- One fold per process, ≤ 8 GB RSS.
- Before each heavy step, check `memory_pressure` (≥ 40 % free) and vm_stat free+inactive ≥ 20 GiB. A ~15 GB
  production build and 416's ≤ 12 GB H1 runs share this Mac tonight.
- No builder edits, no vcc, no cell-eval2. Write only in experiments/418-bulk-composition/.
