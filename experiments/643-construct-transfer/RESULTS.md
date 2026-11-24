# 643 — Does a construct's seed off-target signature travel with the construct into K562? Weakly yes (fragile)

Date 2026-10-02, session [537757]. Analysis only: no build, no `vcc`, no `cell-eval2`, no download. Peak RSS 5.1 GB, 109 s.
Preregistered before any K562 value was read: `PREREGISTRATION.md`, sha256 `721cc81c…`, commit f8f95ecf.

**Class:** science (F16 track), plus a mechanism check for donor weighting. No recipe change.

## Question

- Arc's 2026-10-01 post: "We used the same sgRNA sequences in every cell line".
- The 2025 library reused Replogle's K562 GWPS dual-guide pairs for 297/317 constructs (301).
- F16 (309) found, in H1, that genes whose promoters match a construct's protospacer seed are repressed in that construct's
  cells but not in its partner construct's cells.
- If that signature belongs to the protospacer and not to H1, it should reappear in K562 cells carrying the same pairs.
  K562 GWPS rows would then carry the construct-specific part of Arc's truth wherever Arc reuses Replogle's pairs.

## Data

- 309's 34 constructs (17 two-construct targets), each with its seed-rule candidate genes. The selection is
  sequence-only.
- K562 GWPS cells (scPerturb h5ad, dense, 2.0 M cells × 8,248 genes):
  - 33 of 34 constructs are found by their exact sgID pair (`301/work/vcc2025_spacers.csv`);
  - controls are 10,000 random `control` cells;
  - Δ = mean log1p(CP50k), construct cells minus controls.
- Only 17 constructs have a candidate expressed in K562 (control CP50k ≥ 1) on that 8,248-gene axis, mostly one gene each.

## Result (`k562_f16_643.py`, `summary643.py`; per construct in `work/k562_f16_643.csv`)

| statistic | value |
|---|---|
| **PRIMARY** median D = own − partner over 17 constructs | **−0.0293**; one-sided Wilcoxon p **0.049**; 11/17 negative |
| bootstrap 95 % CI of the median | [−0.086, +0.038] |
| matched null (expression-decile-matched random sets) | median +0.005, p 0.50: the test is resolvable |
| per-construct correlation, D(K562) vs H1's construct contrast (309 `delta_repaired`) | Pearson **+0.51** (permutation p 0.052), Spearman +0.38 (p 0.14) |
| leave-one-construct-out | Wilcoxon p up to 0.088 (dropping CAST_P2); median D −0.038 … −0.028 |
| without STAG2_P2 (36 K562 cells, own −0.43) | Pearson +0.32, Spearman +0.29, Wilcoxon p 0.088 |

The strongest constructs in both lines:

| construct | candidate | K562 own vs partner | H1 contrast |
|---|---|---|---|
| STAG2_P2 | CMC1 | −0.43 vs +0.10 | −0.31 |
| CAST_P2 | UBL5 | −0.15 vs +0.02 | −0.46 |
| DAXX_P2 | RNF115 | −0.17 vs +0.11 | −0.02 |

## Verdict

**REPLICATES by the preregistered rule, but fragile.** The primary test passes at p 0.049 and does not survive leaving
out one construct. The H1↔K562 agreement across constructs rests mostly on STAG2_P2. The evidence leans toward "a seed
off-target signature travels with its protospacer into another line", at low power. It is not established.

**What it changes for D/E/F: nothing now.**

- **Donor weighting.**
  - If constructs carry their off-target signature across lines, h1-test's K562 preference is partly a property of
    the library: 097 a3 tuned the K562 share to 0.60 on h1-test against 0.35 on Orion, and 175 measured sign accuracy
    K562 0.567 vs Orion 0.53–0.54.
  - It would then reach D/E/F only if Arc's 2026 library is Replogle's pairs again. The 2026 files rename the NTCs
    (`non-targeting-3 … -99`, not Replogle's pair IDs), so that cannot be checked from the data.
  - The off-target part itself is small: ~2.4 genes per construct, Δ ~0.08 log1p on H1.
  - 178 already priced a K562-heavy mix (2,1,1) at −0.0012 after the operator on H1. Equal weights stay.
- **Magnitude.** No lever. Knockdown failure was the large construct term (293), and Arc's ≥ 80 % rule removes it.
- **What a board read would mean.** A donor-weight effect on A/B/C that comes from the guides is a property of the
  library, which all six lines share. It would carry to D/E/F, unlike context-tuned knobs. If anyone ever spends a slot
  on donor weights, preregister that reading.

## Better-powered follow-up (science track, not run)

The power limit is the scPerturb axis: 8,248 genes leave most seed candidates unmeasured. Two cleaner tests:

1. Replogle's full-axis files (figshare, large).
2. K562 GWPS vs Replogle RPE1 essential. Same library, same protospacers, two lines. Own vs partner on their shared
   two-construct essential genes is a direct cross-line test with hundreds of constructs.

Neither changes the October recipe. Both would settle F16's line-transferability.

## Reproduce

```bash
uv run python experiments/643-construct-transfer/k562_f16_643.py     # ~2 min, one sorted pass over the dense file
uv run python experiments/643-construct-transfer/summary643.py
```

The first version read scattered rows per construct and re-decompressed the gzip chunks (3,886 × 33) of the 65 GB dense
matrix many times. It was killed after 13 min. The shipped reader touches each row block once.
