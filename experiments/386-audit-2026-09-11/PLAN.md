# 386 — Independent audit of PI numbers (2026-09-11) — PLAN

## Hypothesis
The PI's five stated number-sets (board snapshot, 319p2-vs-319p deltas, mse_alloc
oracle ceiling, call-set arithmetic, dose-closure rule) are exactly reproducible
from the on-disk sources with independently written code. Any deviation is either
a rounding/read-source difference (benign) or a real error (reported as MISMATCH).

## Exact inputs (paths)
- Board: latest-by-mtime `data/leaderboard/*.json` → expect `2026-09-11T1910Z.json`
  (verify by mtime, not by name). Entries keyed by `id`. `score_avg` claimed =
  mean of six `score_*` (verify, don't assume).
- Ours: `data/board-entries/lqWNXX5YaYArdglu6Stn.json` (319p2),
  `data/board-entries/SvzrTOQJAAaNjawzSyrC.json` (319p).
- mse_alloc: `/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/0acea69f-98f2-413e-8a3f-27adb48f049a/scratchpad/mse_alloc.py`
  (copy into `./mse_alloc_pi.py` verbatim, read-only source for comparison) +
  `experiments/323-current-family-validation/work/effects.npz` (`cap15`, genes, targets),
  `metrics-cap15.json` (`mse_legs` numerator/denominator D_p),
  `data/harness/h1-test.h5ad` (truth deltas, TS=5e4, own target excluded, genes by name).
- Call-set arithmetic: pure math from metric defs in
  `docs/reference/vcc2026-metrics-BRIEF.md` §4 (fid yield) and §6 (jaccard);
  board raw numbers n_pred=2,599 (mean? — check provenance; treat as panel-mean
  call-set size), raw jac 0.0409, fid 0.508.
- Dose rule: `experiments/319-cap-variant/PREREGISTRATION-partial-dose2.md` (quote verbatim).

## Estimator / method
- Rows 1–2, 5: fresh `audit_board.py` written by me (no reuse of PI scripts):
  parse leaderboard JSON, find rows by id/team, recompute score_avg as mean of
  six members, rank, count, gap decomposition (Δmember/6), team_id sharing,
  Illumina 09-09→09-11 step from the 09-09 snapshot file.
- Row 3: (a) copy PI script verbatim; (b) write `mse_alloc_mine.py` independently:
  own h5py CSR load, own log1p-pseudobulk profiles, own alignment, own oracle-scale
  algebra S(a)=Σ(a²dd−2a·dt+tt), Δraw=(S(a_or)−S(a_glob))/ΣD. Compare every scalar:
  ΣD·G≈18,077·ΣD conversion, −216 signal term, −0.012 raw, top-10 share 58%,
  Spearman 0.19, D-weighted cosine 0.103. Judge the units question explicitly
  (mse_legs per-gene means vs S as gene-sum).
- Row 4: closed-form algebra, no data load:
  J = I/(n_pred+n_real−I) ⇒ I = J(n_pred+n_real)/(1+J); recall=I/n_real;
  fid=k/max(n_pred,n_real), a=k/I ⇒ a=fid·max(n_pred,n_real)/I per n_real.
  Yield-term: calling exactly 400 with n_real=288 vs 1000.
- Known-answer reproduction FIRST (rule 3): recompute 323 cap15 score_avg from
  `metrics-cap15.json`-independent source — the board-entry files are board-side,
  so instead reproduce: (i) 319p2 score_avg = mean of its six score_* from the
  JSON on disk; (ii) 319p2−319p Δ=+0.000297 from the two JSONs. These use only my
  code + on-disk JSON. The 323 `metrics-cap15.json` members are H1-panel, not
  board, so state scope explicitly and do not conflate.

## Controls / falsifiers
- score_avg≠mean(score_*) for any of the three rows ⇒ PI's "score_avg = mean"
  gloss is wrong; report actual formula.
- Latest-by-mtime file differs from 2026-09-11T1910Z ⇒ audit the actual latest.
- Any board scalar mismatches beyond ±0.00005 (rounding) ⇒ MISMATCH.
- mse_alloc: run PI copy and mine; if mine ≠ PI copy beyond fp noise ⇒ bug-hunt,
  report both; if both agree but ≠ PI-quoted ⇒ MISMATCH with direction.
- Dose rule: full-file read; "flat ⇒ no third dose" must be an actual sentence;
  else MISMATCH (quote what it really conditions the third point on).

## What kills the idea
- A single MISMATCH row does not kill the audit — the audit IS the deliverable.
  What would kill a sub-claim: score_avg not equal to the mean of six (kills row
  premise); h1-test.h5ad unreadable under the memory cap (kills row 3 → PARTIAL
  with the algebra-only verdict on units).

## Resources
- Memory gate already checked: free 86% (2026-09-11 session start). h1-test.h5ad
  is 677 MB on disk; CSR float load ≈ 1–2 GiB peak, one job at a time,
  NUMBA/OMP threads=4. Budget ≤ 8 GiB, no harness/builder/submission binaries.
