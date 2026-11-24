# 801 — the lab's fresh held-out folds and the bulk fold engine (worker mandate)

Written 2026-10-15 ~22:46Z by the lab lead ([c61d4c]). You own `experiments/801-lab-folds/` and write only there.
Read first: `docs/LAB-2026-10-16.md`, `experiments/800-checker/DESIGN.md` (your work is its §3 and the bulk half of §4),
then `experiments/522-october-bulk-key/RESULTS.md` (§1, §4b), `experiments/534-generalization-audit/RESULTS-B.md`,
`experiments/551-r6-c10b/RESULTS.md`, `experiments/708-def-panel-repricing/PREREGISTRATION.md`.

## Objective
The lab's checker needs held-out bulk folds nobody has tuned on, and an engine that prices ANY bulk candidate on them
member by member, from per-target arrays, exactly as production's arithmetic would. Deliver:

1. **Panels by the rule "arc800"** (DESIGN §3.1), drawn by ONE script `draw801.py`, in order, seeds 800101 … 800108,
   each draw's targets + pool size + sha256 written to `work/draws.json` at draw time (as `draw534.py` does; reuse its
   exclusion code: every `pert_counts.csv` under experiments/ and data/ at draw time, plus earlier lab draws).
   - Pool: official axis ∩ union of the three cancer donors' perts (k562_gwps, orion_hct116, orion_hek293t) − K562
     essential list − excluded panels.
   - Stratify to A/B/C's own panel shares on two axes, frozen in the script BEFORE the first draw:
     (a) coverage pattern: which of the 3 donors + the CD4 atlas carry the target (A/B/C: K562 269/300);
     (b) strength bin: GWCD4i Rest `n_downstream` (data/external/cd4t/DE_stats.suppl_table.csv) in A/B/C-quintile bins,
         plus one bin for targets the atlas lacks. CD4 is the "pilot line" here (Arc used an H1 pilot); it is never a fold.
   - Basal expression: the target's mean control CPM in A/B/C's context files ≥ the A/B/C panel's 5th percentile.
   - 300 targets each; assert 0 essential, 0 overlap with any excluded panel or earlier lab draw, all on the axis.
   - L1–L6 = WORKING, L7–L8 = SEALED. **For L7–L8 write only the target lists, seeds and sha256: no bundle, no CD4 bank,
     no operator, no dump, no statistic involving them** (they are the sealed reserve; the lab lead opens them once).
   - Diagnostics for L1–L6 only, written before any score: strength-bin shares, coverage, CD4 coverage, basal-expression
     quantiles, vs A/B/C and vs 534's efr/efp rules (`work/diag801.json`).
   - Bundles for L1–L6 the way `draw534.py` builds them (relabelled A/B/C controls as APFS clones, `cp -c`).
2. **A vendored bulk fold engine**, `engine801.py` (+ whatever it needs, copied into 801 with a header naming the source
   file and its git blob sha; arithmetic verbatim, only paths/guards changed, as 534 copied 438 and 551 copied 534):
   - per fold (panel, held-out cancer line), the dump708 arithmetic: the fold's sources = the other two cancer donors,
     the fold's held-out operator, composition h2c1 (H1 ×2, CD4 ×1, DEF rules), carriers h1c1, the fold's LODO-V1 key,
     frames, truth (the held line's rows on the scored genes), betas;
   - an arm = (key weights w, per-target sharpening q vector, optional transform of the composition rows `fin`, optional
     external rows to blend/swap per target). The engine renders (minnorm), realises ballast at KMAX 16 against the
     carriers, and returns PER-TARGET, PER-FRAME arrays: x = ⟨P, T⟩, a = ⟨P, P⟩, d = ⟨T, T⟩ (scored genes), literal pds
     (scorer's rank rule), board-emulated pds λ 0.5, the ballast K, and the arm's change energy split into scored /
     unscored genes and by A/B/C control-CPM bin [0,5), [5,20), [20,100), [100,∞) (DESIGN §7, 702 rule 3);
   - the incumbent arm: percentile cap P 96.64697522769474 on the fold key (522's `pct_hi` / `wclip(0.5, hi)`), tercilec
     sharpening (551's `q_p2g`, carried-only cut points, 2.5 / 2.0 / 1.75);
   - named arms for the known answers: `cap3` (wclip 0.5, 3), `cap6`, `pct`, `pct+P2g` (= incumbent), `pct+S_P2` and
     `cap3+S_P2` (strength terciles, 522/551's S_P2), `cap3+P2g`.
   - **Known answers, asserted on load, abort unless |d| < 1e-12:** 522's dump/arm rows on efp (dump522 + 546o's
     `base q1.75` / `fixed 6` / `p96.6`), 551's efp KA rows (P96.647 + P2g, + S_P2), 708's rehearsal row if present. Run
     them on EC2 (the dumps are ~6 GB) and record them in RESULTS.md.
   - The same engine must also load the **A/B/C panel's own PLOLO folds** (k562, hct116, hek293t, kolf2) through 522's
     `fold522.Abc` / 518's path so the checker has the board's panel as a retrodiction surface; same output arrays.
3. **The heavy steps on EC2** (`experiments/800-checker/cloud/labcloud.sh`; read its header): per WORKING panel, CD4
   extraction (534's `cd4_534.py` route), the three held-out operator fits (`scripts/build_union_whitener.py`, the exact
   commands in `experiments/534-generalization-audit/work/partB/state534.json`, `--max-rss-gib 14`), the fold dumps.
   Find the data closure (a `sys.addaudithook` recording every file the engine opens is the quickest way), stage it
   (`labcloud.sh stage <paths>`), and run everything for L1–L6 in as few launches as possible (an r7i.8xlarge runs ~8
   fits at once; set OPENBLAS/OMP threads so they do not oversubscribe). Keep the operators on S3 (they are the fold's
   identity; record their sha256), download only the dumps (~100 MB each) into `work/dumps/`.
4. **A registry** `work/registry.json`: every fold id (`L1/k562_gwps` …, `ABC/k562_gwps` …), its dump path + sha256, its
   operator sha256, target count, coverage, and status (working / abc / sealed-drawn).
5. **A light pricing smoke** (local is fine if < 2 GB, else EC2): incumbent vs `cap3` and vs `pct+S_P2` on every fold,
   per-target arrays saved as `work/smoke/<fold>-<arm>.npz`, plus the fold-level r (relative aggregate cosine) and Δpds.
   Do NOT aggregate into board units or bins; the lab lead's check800 does that.

## Rules and budget
- Cloud budget for this mandate: **$40**. Prefer one r7i.8xlarge batch for fits + dumps; terminate everything.
- Nothing here reads R6/R7 reserves' scores (S1–S4, T1–T4 are only excluded from the pool).
- If a known answer fails, stop and report; do not price anything on an engine that has not reproduced its source.
- Time: aim to have L1–L6 dumps and the KA table by ~10-16 18:00Z. Write RESULTS.md as you go (state survives restarts).

## Return
RESULTS.md: the draws (sha256, pools, diagnostics vs A/B/C), the KA table, the registry, the cloud runs and cost, the
engine's interface (function signatures + what each array means), and anything that did not work.
