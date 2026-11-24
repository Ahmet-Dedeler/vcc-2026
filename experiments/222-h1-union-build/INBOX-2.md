# INBOX 2 — read before you write the verdict. Two things in `six222.json`, one of them changes the prediction.

Good table. `pds` up, `precision` up, `nmae` up only **+0.00039** against the incumbent (falsifier 4
armed at >0.002 — **passes**, say so explicitly with the number). Two items before BUILDABLE.

## 1. The 3-cancer number is the honest predictor, not the 4-fold number — and you must use it

Your per-fold `pds` deltas vs the incumbent:

| fold | cd4 f=0.34 | h1union g=0.28 | Δ |
|---|---|---|---|
| k562_gwps | 0.79913 | 0.80879 | **+0.00966** |
| orion_hct116 | 0.78958 | 0.79098 | +0.00140 |
| orion_hek293t | 0.75625 | 0.75987 | +0.00362 |
| **kolf2** | 0.62096 | **0.64295** | **+0.02199** |

All four are positive — genuinely good, and the mechanism is coherent: **you added H1, an embryonic
stem line, and the fold that gains most by 2.3× is `kolf2`, the other stem-like fold.**

**That is exactly why the 4-fold mean overstates the board.** `pds_4folds` (+0.009166) is inflated by
`kolf2`. **The challenge contexts are Jurkat (A), HeLa (B) and CAL33 (C) — all cancer lines, none
stem-like** (041). The transferable estimate is `pds_3cancer` = **+0.004891**, and even that leans on
`k562_gwps`, whose +0.00966 is 2.7× the other two cancer folds combined-average.

So:
- **Predict the board from `pds_3cancer`, and state the `kolf2`-driven 4-fold number separately as
  the non-transferable part.** If the preregistration's +0.0055…+0.0105 band was built on the 4-fold
  number, **say so and reissue the band** — a preregistered band derived from an inflated statistic
  should be corrected *before* the result, not after.
- Re-read the three-way pds branch in `DECISION-TREE.md` against the 3-cancer number too, and say
  which branch it selects. If the branch flips, that is a finding, not an inconvenience.

## 2. Your one known-answer FAIL — explain it or it blocks the build

```
FAIL  cd4 f=0.34 <- INCUMBENT  delivered_norm_ratio  got 1.00158932  216 1.00094558  |d| 6.44e-04
```

Every other cell matched **0.00e+00**, so this is not noise. The likely cause is benign and is the
point of your experiment: **216 sized with the all-rows median and you size with the carrying-row
median**, so 216's number is the buggy one and yours is the correction. **If that is it, prove it** —
recompute the incumbent's ratio both ways and show the all-rows path reproducing 216's
`1.00094558` exactly. A known-answer row that is explained is a pass; one that is asserted is not.

If it is **not** that, stop and report. Note this is a **harness** quantity — the shipped package is
already independently verified (config diff exact, the `cd4f034` rebuild hashing byte-identical to
`7036e89f…`, 900/900 blocks moved, carrying:non-carrying 3.3×), so this cannot silently corrupt what
we ship. It can corrupt what we *predict*, which is why it must be closed before the verdict.

## 3. Bank the October number — it is now measured

Your coverage table gives the over-dose factor at the pessimistic draw: **1.237–1.245 across all four
folds**, versus **1.012–1.017** at today's ~92% coverage. That is the D/E/F risk quantified for the
first time. Put it in `RESULTS.md` as a standalone finding with the per-fold table — it belongs in
`docs/DEF-RUNBOOK.md` §5 and outlives this build.

**Still do not submit and do not name the model.** Verdict line, then stop.
