FOLLOW-UP 2 from the lab lead (2026-10-16 08:05Z). Your context is intact. A production session (712, commit 9bf8b6c9) points
out that A/B/C's PLOLO folds (lib438 PANELS abc/abcdef → 438, 522, our BULK-ABC via fold522.Abc / 518's dumps) use 062's cloud
whiteners (`experiments/062-cd4-post-whiten/cloud-output/whiteners/held_<line>.npz`, basis = every donor row, K562 essentials
included, no per-donor cap), while real/efp/efr, our L folds, production e2f65ecf and the 10-22 day path all use
`build_union_whitener.py --basis-prune essential` with the 12k per-donor cap. So "fold cap hi 5.2–7.0 on A/B/C vs 7.8–19.9 on L
folds" mixes operator recipe with panel, and BULK-ABC's key-class ρ (1.46) may be calibrated on an instrument the day path does
not use. 712 also found the raw delta-row-norm tail is the same on A/B/C, H1 2025, our draws and L1–L6 (1.9–2.4; Spearman with
hi +0.08). The decisive test, on EC2 (budget $15):

1. Refit A/B/C's three held-out operators with the DAY recipe (`--basis-prune essential`, the 12k cap, the same flags as your
   L folds / state534's commands, on A/B/C's own bundle `data/challenge`), keep them on S3, record sha256.
2. Rebuild A/B/C's three cancer folds with those operators — call them `ABCp/<held>` — through the path that matches the
   BOARD's A/B/C rules (lib438 `abc` rules, CD4 ×1.5 all-rows sizing, SVD-5 deflation re-fitted, as 518's dumps / fold522.Abc;
   if only dump708's D/E/F-rules path can take new operators, build that too and say which is which). KA: with 062's operators
   your rebuild must reproduce the existing ABC fold arrays (fin, kv, truth) to < 1e-12 before you swap operators.
3. Report per fold: the percentile-cap hi (P 96.64697522769474 of key/median over carried) under 062's vs the day recipe's
   operator, and the same for the L folds' known values. Then run `experiments/804-known-answers/ka804.py`'s arm set on the
   ABCp folds (it loads folds by id through your engine: register `ABCp/<held>` in your registry the way L folds are), and
   leave the standard-format arrays under `work/ka-abcp/` (S3 + small local copies). The lab lead regrades.
4. Rules as before (window800, gate602 / EC2 only, explicit S3 mirror, no data/ writes). Finish with a short message.
