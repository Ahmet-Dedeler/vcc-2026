# 440 — The empty-union Jaccard lever: A/B/C has no null targets. CLOSED (2026-09-24, coordinator, ~20 min, no build)

## Hypothesis
`de_sig_jaccard` returns **1.0 when both sides call nothing** (`cell_eval2/metrics/de.py:de_sig_jaccard`, "an empty union
returns 1.0"; the own-gene anti-join can create such rows). We call 2,600 genes on every target, so any target with an
empty reference set (n_real = 0 after own-gene exclusion) scores J = 0 for us. With ballast carrying the bulk, a target's
carriers could be emptied (z_cap on every gene, order kept) at no pds/mse cost.
- Price per target-context: correct empty +1/300 raw jac = +0.009 scaled = **+0.0005 overall**; wrong empty costs fid_p
  ~0.5 → **−0.0003 overall**. Break-even P(null) ≈ 0.38–0.44.
- H1-test has 3/85 targets with n_real = 0 (PCCB, TSKU, TRAFD1), and 192 measured the board's n_real median at 31–45
  (0.35–0.5× H1), so A/B/C might have had more.

## Measurements
1. **Predictability on H1 (the only local truth):** GWCD4i `n_downstream` (per target, max over Rest/Stim8hr/Stim48hr,
   ~700 cells per target; `199/work/GWCD4i.DE_stats.h5ad` obs) tracks H1 n_real at Spearman +0.47 (78 covered targets).
   All three H1 nulls have CD4 max ≤ 2, but the rule CD4 max ≤ 2 catches 22 targets of which only 3 are exact zeros (14 %;
   median n_real 10.5). Precision is far below the break-even. Basal expression of the target does not predict nulls at
   all (Spearman(n_real, target CPM) = +0.04): H1's nulls are well expressed and knocked down 3–5 log2 units.
2. **The decisive measurement is free, on the board.** About 50 teams submitted plain control copies or resamples
   (pds exactly 0.5, fid 0.0000; e.g. "b1-control-copy", "ntc-copy-v1", "control-paste", "control-resample-v1",
   "OfficialContextMeanNoE…"). They call nothing on any target, so their raw jac is exactly the fraction of target-contexts
   with an empty reference set. **Every one of them reads jac = 0.00000** (the smallest non-zero value would be
   1/900 = 0.0011). Scan: `data/leaderboard/*.json`, 159 snapshots, entries with fid < 0.2.

## Verdict
**A/B/C has zero null targets: every target has ≥ 1 reference-significant non-target gene in every context.** This is
consistent with Arc's stated selection ("all had more than 80 % on-target knockdown and were selected to provide a strong
set of perturbations and responses"). D/E/F uses the same selection, so expect the same there. The lever does not exist;
no slot was spent. By-product: the board's control-copy entries are a free truth instrument for any property of the
reference that a no-call submission reads exactly (here jac; also reach ≈ 0.05–0.10 is the random-order floor).

Also noted (not a lever): on H1, per-perturbation compositional shift (median lfc of null CPM-block genes) is tiny,
median |shift| 0.0028 log2 against per-gene noise sd 0.057, so even an oracle shift sign puts null-gene sign accuracy at
only 0.529. fid's null calls are a coin flip for a structural reason, not a design choice.
