# 420: price the bulk STACK for the next arms (mandate from the PI, 2026-09-22 21:30Z)

**Current state.** Read `docs/FINDINGS.md` §0–§1, `experiments/414-bulk-reprice/RESULTS.md`,
`experiments/418-bulk-composition/RESULTS.md`, `experiments/415-lodo-key/RESULTS.md` §3–§5 and
`experiments/417-bulk-sharpen/RESULTS.md` + `PREREGISTRATION.md`.

The shipped bulk (410b = 0.2365722) is the PP composition times 410's key at clip(0.3, 3), delivered through ballast.
417a adds bulk sharpening q 1.25. It is built; its measured render energy is 0.961/0.958/0.961 of 410b's, and it is
uploaded at 00:01Z on 09-23. 417b (q 1.5) and 417c/417d (carriers from an unwhitened pass) are building.

Each lever was priced alone, and every one is small:
- 418: H1-union ×2 (+0.002…+0.003), CD4 fraction ×1.5 (+0.0015…+0.0022);
- 415: the LODO key, +0.0017;
- 414: key floor 0.5, +0.0012; per-context re-levelling, +0.0024 after sharpening.

The next arms will be stacks, so they have to be priced as stacks.

**Territory.** 416/419 do the DE carriers on H1 and 417 does the builds. Do not touch the carriers, the builder or the
board. Using 418's machinery (`lib418.py`, `comp418.py`, `abc418.py`, `analyse418.py`: PP composition, a23 read-back,
cliff-aware mse with the A/B/C real-panel energy route, board-emulated pds at λ 0.5, history-calibrated k), price
these BULK designs vs 417a's bulk (PP + 410 key clip(0.3, 3) + sharpening q 1.25):
- (a) the stack H1 ×2 + CD4 ×1.5, and each alone as known answers against 418's table;
- (b) sharpening dose q ∈ {1.35, 1.5, 1.75, 2.0} on PP|K3. Also run each through 418's `ballast418.py` realisability:
  what fraction of blocks hit K = 16, and the realised-vs-design cos. Add KMAX 24 and 32 as an emulation column: the
  carrier count drops to 376 / 368, so say what that would cost the DE channel if you can bound it, otherwise flag it;
- (c) the key variants on the sharpened bulk: 415's LODO V1_PP at clip(0.3, 3), 410's key at clip(0.5, 3) (floor 0.5),
  and LODO at clip(0.5, 3);
- (d) the best combination of (a)–(c) you find, with per-context re-levelling EV-optimal factors on top
  (`experiments/417-bulk-sharpen/ctxlevel417.py` has the model; here use PLOLO's A/B/C-energy route).

For every row report:
- Δrel cos, Δpds_bd, Δraw mse at k = 1 and 1.407, and Δoverall;
- the per-fold values with a sign count;
- the A/B/C render energy ratio vs 417a (flag > 1.03 or < 0.90 in any context);
- the builder expression. Existing knobs are `VCC_BULK_SHARPEN_Q`, `VCC_PERT_SCALE` artifacts (`wshare_lo0.3_hi3`,
  `wshare_lo0.5_hi3`, `lodo415_lo0.3_hi3`, …) and `VCC_BULK_CTX_SCALE`. Composition changes ride
  `VCC_BULK_DUMP` / `VCC_BULK_LOAD`, a dump made with that composition env.

**Known answers first** (abort on failure): 418's `PP|K3`, `PP|K3+sh`, `b:h1x2|K3+sh`, `a:cd4x1.5|K3+sh` rows
(`work/comp418-<fold>.json`) and 414's `Kinf+sh@1.25`/`Kinf+sh@1.5` rows must reproduce exactly.

Cancer folds decide. kolf2 is reported but has no vote, and it is biased on the H1 arms (418 §3).

**Deliverable:** `experiments/420-bulk-stack/RESULTS.md`:
- the ranked table;
- the single recommended next bulk arm as a one-variable step from 417a;
- the recommended sequence of the next three bulk arms, each one variable from the previous, with predicted Δoverall.

**Constraints.**
- One fold per process, ≤ 8 GB RSS.
- Check `memory_pressure` (≥ 40 % free) and free+inactive ≥ 20 GiB before each heavy step: ~15 GB builds and a 5 GB H1
  harness share this Mac.
- No builder edits, no vcc, no cell-eval2. Write only in experiments/420-bulk-stack/.
