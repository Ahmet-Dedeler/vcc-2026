# v8 — the second submission, and the evidence behind each input

**Config:** v6 plus exact-total quasi emission at `spread 0.25`. Nothing else changed.
**Artefacts:** `submissions/v8.h5ad` (20.3 GB, 2.54B stored entries, 53.5% of cap),
`submissions/v8.vcc` (3.67 GB).

## Why this is an unusually clean submission

The pseudobulk is **bit-identical to v6's** — verified on the 21 GB artefacts, not in
principle: gene totals identical on 18,533/18,533 genes for every perturbation checked, max
|delta log-pseudobulk| exactly 0, identical total counts. So `pds`, `mse` and `nmae` **cannot
move**. Only `fid`, `reach` and `jac` can, and only through the Wilcoxon's view of cell-to-cell
spread.

That makes the board result a clean single-family measurement rather than a mixture.

## Every input, and the harness that governs it

| knob | value | evidence | governing harness |
|---|---|---|---|
| donors | 1:1:1 | equal weighting is best *once whitening is on*; worst without it (019) | LODO |
| whitening | gamma 0.85, full 45,531-perturbation basis | +0.070 macro PDS, replicated on a 4-line LODO (017/019/020/024) | LODO (exact `pds`) |
| renorm | preserve | global rescaling silently down-weights the 23% of tested genes off the operator axis (025) | LODO |
| **emission** | **quasi, spread 0.25** | plateau 0.20–0.30 on both harnesses; `reach` +14%, `fid` +0.004, `jac` −0.004 (034) | h1-fast (400 reference cells) |
| CPM gate | 5 | plateau 2–10; gate 0 costs `pds`, gate 20 costs `nmae` 0.9449 → 0.9613 (034d) | h1-fast |
| alpha | 1.0796 | measured optimum; 1.5 pushes `nmae` past 1.0 (031) | somatic + h1-fast |
| sharpening | **off** | refuted: q=2 takes `pds` 0.763 → 0.832 and `nmae` 0.961 → **1.260** (031) | somatic |
| targeted variance | **off** | confidence-targeted spread lands inside the noise of uniform (034c) | h1-fast |

## Expected effect

`reach` 0.1584 → ~0.180 (+0.031 scaled), `fid` 0.5016 → ~0.506 (+0.012 scaled), `jac` −0.001,
everything else frozen: **≈ +0.007 overall**, from 0.1369 to ~0.144.

The downside is bounded to `jac`, which contributes 0.002 scaled to our total score. There is no
configuration of the emitter that can cost us `pds`.

## Operational note

The `.vcc` is **3.67 GB against v6's 936 MB**. The even emitter's base/base+1 counts are
pathologically compressible and realistic ones are not — which is itself a symptom of how
unnatural minimum-variance emission is, and the same artefact that inflates our DE calls.
Budget ~17 minutes for the upload rather than 4.
