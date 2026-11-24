# 593 — Adversarial claims audit of manuscript Draft v16 (cloud agent, Grok 4.7 extra-high)

You are auditing a private scientific manuscript in this repository. Your job is to find errors, not to agree. Work only in
`experiments/593-cloud-v16-audit/` on your branch. Do not edit anything else. Never post, upload or share the manuscript or
any part of it outside this repository: it is private until after Nov 5, 2026. Do not run `vcc`, do not download data. This
VM has no gitignored data, and everything you need is tracked.

## Read first
- `AGENTS.md` §4 (scoring: six members, scaled s = (u − b)/(r − b), overall = mean of six) and §11.1 (method rules).
- `docs/reference/vcc2026-metrics-BRIEF.md` §4–§8 (the DE members' definitions).
- The previous audit you are modelled on: `experiments/503-cloud-claims-audit/RESULTS.md` (its standard: trace every
  number to the primary artifact, recompute the arithmetic, list discrepancies with file:line and the correct value).

## The object
- Draft v15.1 (published): `experiments/294-writeup/PAPER.md`, `CLAIMS.md`, `REFERENCES.md`, `figs/MANIFEST.md`.
- Draft v16 (candidate): `experiments/592-manuscript-v16/paper/` (same files). The edits are in
  `experiments/592-manuscript-v16/apply592.py` (anchored replacements) and the author's own check is `verify592.py`
  (`verify592.out`). The author's account is `experiments/592-manuscript-v16/RESULTS.md`.
- Diff v15.1 → v16 yourself (`diff` / `git diff --no-index`). Audit EVERY new or changed sentence, table row and claim row:
  the Round 15 revision note, §3(h) (the October path; C-82 amended; C-94), §3(j) (the dose amendment; C-89 amended;
  C-97), §3(k) (the DE read-noise correction; C-96), §3(l) (the call-set composition read; C-95), §4.1, §4.5 T12, §4.7,
  §5.9 (line identification; C-98), the abstract and plain-language summary changes, the F-12 caption, the two new
  references.

## Primary artifacts (all tracked)
- Board entries: `data/board-entries/<entry id>.json` (`entry.score_*` scaled, raw members by metric name). Ids used:
  446b `UDgjOKL9gI1udAexq6ue`, 446bj2 `YvVru3DAF89Jiyr34olR`, 446bj2z `be90smZUDcIV7j5MgQ03`, 446bd2
  `R5uE5Rve3R9nWXGbeSLM`, 453dr `q3VuTImh4FRGKXH2waCp`, 457dr `1FfirZsGsAZvaXYTd39J`, 460dr `LrkhsIXWMjEsvRXSGToa`,
  446bj2zcp `K1Q0y2OFYjTNRqh9jKfk`, 458d `b7UI3NRy1rdVNRRyZwmU`, 458g `JL8llNLHaTyDQzm798TV`.
- Board-result records: `experiments/460-lodo-dayof/BOARD-RESULT-460dr.md`, `experiments/535-callprop-arm/BOARD-RESULT-446bj2zcp.md`,
  `experiments/458-arc-lab-programme/BOARD-RESULT-458d-458g.md` (with its 10-01 addendum), their PREREGISTRATION files.
- Build metas: `submissions/{446b,458,458d,458g}-cap15-…h5ad.meta.json`, key `output.rc_callset.prior_lambda`.
- The gate sizes and the exact decomposition: `experiments/468-prior-gate/RESULTS.md`.
- DE noise: `experiments/537-arm-pipeline/RESULTS.md` §2 (and `work/noise537b.json`, `work/redraw537.json`),
  `experiments/547-de-noise-check/RESULTS.md` (read its Correction), `experiments/548-cloud-de-noise-audit/RESULTS.md`,
  `experiments/230-toll-magnitude/work/fit230.json`.
- Dose: `experiments/554-dose-guard/RESULTS.md`.
- Line identification: `experiments/590-identity-validation/RESULTS.md`, `PREREGISTRATION.md`, `work/match590.csv`,
  `work/match590_flexbias_fold{0,1}.csv`, `work/summary590.json`, `analyse590.py`.
- H1 prices: `experiments/528-callprop-block/RESULTS.md`, `experiments/535-callprop-arm/PREREGISTRATION.md`.

## What to do
1. Recompute every number the new text states from the primary artifact (not from another RESULTS file when the primary
   exists). Show the recomputation in your RESULTS (script + output). Write the script as
   `experiments/593-cloud-v16-audit/audit593.py` and run it.
2. Check each inference, not only each number. In particular, attack these three:
   a. §3(j)/C-97's bound "at least 198 of the 648 strong-target programme blocks sat at λ ≤ 0.433 in both entries".
      Is the arithmetic right (900 blocks; a median from min/median/max only; 252 weak-target blocks)? Does it rely on
      anything not shown? Is "458g's programme blocks are byte-identical to 458d's" supported by 468's record?
   b. §3(h)/C-82's re-reading that the October path's DE share is "most likely 446b's own realisation". Recompute the DE
      shares of every pair cited; consider the alternative (a systematic cost of the path) and say which the board
      entries on record support, and how strongly.
   c. §3(k)/C-96's noise figure (≈ 0.00025 overall per re-roll) and the claim that bulk-only pairs are re-rolls. Check it
      against 537 §2.2–§2.6, 547's correction and 548's caveats (e.g. 548's earlier-family rms 0.00321 and its "no seed-only
      pair exists"). Is any statement in v16 stronger than those records allow?
3. Check §5.9/C-98 against 590's files (counts 402 / 390 / 146; accuracy 0.985 / 0.995; the 0.51–0.61, ≥ 0.68 q01,
   0.719–0.772 values; the 2 of 19 cost). Note anything the paper states that 590's CSVs do not support.
4. Check that no new sentence contradicts an unchanged part of the paper (search for stale statements of the same facts,
   e.g. any remaining sentence that still prices DE noise at 0.0045 as a re-roll, or calls 458d's harm "full dose").
5. Style is not your job, except where wording changes the claim (e.g. "identical" vs "within 2e-8").

## Deliverable (on your branch `cloud/593-v16-audit`, only in `experiments/593-cloud-v16-audit/`)
`RESULTS.md`: a verdict line, then a table — location (file:line in the v16 copy), the statement, your recomputed value,
the primary artifact, verdict (CONFIRMED / WRONG / OVERSTATED / UNSUPPORTED), and the exact corrected wording for every
non-CONFIRMED row. Then the three attacks (a–c) with your reasoning. Plus `audit593.py` and its output file. Commit
those files to your branch and stop. Do not open a pull request.
