# 354 — Audit of 343/344: all claims MATCH (two conventions pinned down, one rounding note)

Status: **COMPLETE**. Independent code (`audit343.py`, `audit344.py`; 343/344 scripts never
opened, only their RESULTS.md). Full numbers in `work/audit.json`.

## Known-answer reproduction
- 323 recorded raw NMAE `0.9508364630325662` → recomputed mean over eligible targets identical
  **bit-for-bit**; per-target CSV mean identical; 58 eligible rows. Plumbing trusted from here.
- `reference_de.parquet` == `reference.npz` lfc/q to `0.0` max-abs-diff (spot target AASS).
- Design map: exactly **2** universe genes absent from the effects axis (→ design 0), as claimed.

## 343 verdicts (all MATCH)

| # | claim | mine | verdict |
|---|---|---|---|
| 1 | mean NMAE 0.950836, 58 eligible, equals 323 | 0.9508364630325662 exact, 58, exact | MATCH |
| 2 | mag ratio 0.205; wacc 0.588 deliv / 0.600 design | 0.20548 / 0.58797 / 0.60019 (**mean-of-targets**; pooled = 0.20647/0.56192/0.56872) | MATCH |
| 3 | CF all-signs-right 0.803; LS-scalar 1.005 (mean 3.1×); grid 0.928 (mean 1.48×) | 0.80339; **1.0049** (mean a 3.111); 0.92761 (mean best 1.4828, grid 0.25–6 step 0.25) | MATCH |
| 4 | pool 64,931; UP-share 0.475; called 1.000 top / 0.033 bottom basal-Q; mag ratio 0.114 top-\|lfc\|-Q | 64,931; 0.47521 unweighted; 1.000 / 0.03281; 0.11404 | MATCH |
| 5 | UP wacc 0.516 vs 0.558; low-basal(<24)&UP 0.445 vs 0.536 | 0.51604/0.55803; 0.44524/0.53572 | MATCH |
| 6 | bin 5.0–15.7 P=0.405, others ≤0.074; 6.7% calls <10 CPM; 15.1% NMAE mass | edges 5.001–15.733 exact from universe basal quintiles; 0.40494; others 0.0683–**0.0742**; 6.730%; 15.074% denom-mass | MATCH* |

\*Two convention pins: (a) claim-2 figures are **mean-of-per-target**, not pooled — pooled differs
by 0.026/0.031 on wacc, so future citations must state this. (b) "Least-squares scalar" is LS on
**magnitudes** (a=⟨|ĥ|,|l|⟩/⟨ĥ,ĥ⟩); naive signed-LS gives 0.935/mean-a-1.24, not 1.005/3.1×.
Rounding note on claim 6: top bin reads 0.07421, a hair above "≤0.074" (rounds to 0.074; probable
bin-edge tie convention). NMAE-mass claim is the **denominator** share (numerator variant 15.63%).

## 344 verdicts, context A (all MATCH, hand-rolled Wilcoxon/BH/Spearman)

| claim | reported | mine | verdict |
|---|---|---|---|
| 5-CPM universe | 9,929 | 9,929 exact | MATCH |
| constructs w/ any BH hit | 16/46 | 16/46 exact | MATCH |
| median hits | 0 | 0 exact (30/46 constructs zero) | MATCH |
| Spearman(mean\|z\|, ctrl CPM) | 0.097 | 0.09734 | MATCH |

Side flag (no verdict): `non-targeting-11` carries **880** significant genes; next-highest is 14.
Headline unaffected, but a single-construct artefact exists in context A and any construct-pooled
statistic should check robustness to dropping it.

## Limits
- 344: pool convention = construct cells tested against all 18,400 (self included), the literal
  reading; exact headline match means the convention question is moot at reported precision.
- Exact per-gene z/p agreement with the scorer's implementation not checked — only the four
  headline statistics, which is what the mandate asked.
- 343 claim-6 "35% top-decile" and sub-10-CPM 0.508/0.033 addendum figures were out of scope,
  not re-derived.

## Reproduce
```
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/354-audit-343-344/audit343.py
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/354-audit-343-344/audit344.py
```
Peak RSS ≈ 0.5 GiB (343, npz mmap) / ≈ 3 GiB (344, CSC build); wall ≈ 20 s / ≈ 5 min. Two
consecutive 344 runs deterministic (identical outputs).
