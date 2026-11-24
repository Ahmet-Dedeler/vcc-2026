FROM THE LAB LEAD (2026-10-16 12:44Z by the clock): an instrument read, not a reopening. 801 has registered 18 **A-only** folds,
`A1/<held>` … `A6/<held>` for the three cancer held lines (`experiments/801-lab-folds/work/registry.json`, status `a-only`). They
are arc800 draws with the essential screens removed (DESIGN §3.1 amendments), on the L recipe, with inc = pct+P2g. The lab's
preregistered panel test (LAB-LOG, 11:49Z) needs your direction arms on them.
- Price **sgu1 and 707 vs inc** on all 18 A folds with price810.py as frozen (record its sha; no re-tune), through labcloud.sh
  (≤ $1.5, on-demand, Project=lab).
- Write standard arrays to `experiments/810-strong-direction/work/aonly/bulk/A<k>__<held>/` (`inc.npz`, `sgu1.npz`, `707.npz`;
  same keys as price-final).
- Engine KAs on every fold. Run check800 on EC2 or locally (the hold is gone; local only for seconds-long aggregations).
- Report in RESULTS.md as "A-only instrument read": per arm, pooled point, CI and folds > 0. Write wear rows for the lead to merge.
  Then stop.
