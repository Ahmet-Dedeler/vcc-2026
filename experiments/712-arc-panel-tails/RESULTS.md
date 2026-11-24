# 712 — is A/B/C's light key tail Arc's selection? Result: NO READ (the proxy fails its gate), and the A/B/C operators are a confound

Run 2026-10-16 ~08:00Z (science [d80a10]), `tail712.py`, 5 s, peak RSS 1.3 GB; `work/tail712.json`. Preregistered in
`PREREGISTRATION.md` (committed before the run).

## Result
The donor-strength tail ratio (P 96.647 / median of each covered target's delta-row norm in the line) is about the same
everywhere:
- Arc's A/B/C: 1.88 / 2.04 / 2.05 (K562 / HCT116 / HEK293T).
- Arc's H1 2025: 2.28 / 1.88 / 1.88.
- Our real / efp / efr: 1.93–2.14.
- The lab's L1–L6: 1.89–2.39.

**Validity gate FAILS.** Spearman(proxy, 522's `hi546`) over the 9 real / efp / efr folds is +0.083 (the gate needs ≥ 0.6),
and A/B/C is not below those folds on every line. By the preregistration this is NO READ: raw target strength does not
carry the LODO key's tail. That is a finding in itself: Arc's selection does not produce a distinct strength distribution,
and A/B/C's raw strength tail matches the draws'.

## A confound found while following up (not part of the preregistered reading)
A/B/C's fold operators are not built the way every other panel's are:
- **A/B/C's PLOLO folds** (lib438 PANELS abc / abcdef: 438, 522, 801's BULK-ABC, the lab's ρ_ABC) use 062's cloud whiteners
  (`experiments/062-cd4-post-whiten/cloud-output/whiteners/held_<line>.npz`). Their basis is every donor row: orion_hct116
  17,767 rows and orion_hek293t 17,977 for the held-K562 fold. That includes the 2,021 K562-essential perturbations, with no
  per-donor cap ("production equal-source renorm").
- **real / efp / efr** (233 / 238), **the lab's L folds** (801), **the production operator** e2f65ecf and **the 10-22 day
  path** (`build_union_whitener.py --basis-prune essential`) all drop the essentials and cap each donor at 12,000 rows
  (e.g. 10,807 / 10,803 kept for held-K562).

So "A/B/C's fold cap hi 5.2–7.0 vs 7.8–19.9 everywhere else" compares two operator recipes as well as panels. A/B/C's light
tail may be a property of 062's unpruned basis, not of A/B/C's targets. If it is, two other things follow:
- the lab's same-panel ratios (BULK-ABC), measured on 062-operator folds, are calibrated on a different instrument from
  708's day path;
- D/E/F's day-of fold caps will look heavy for any panel.

**Decisive test (cheap):** refit A/B/C's three held-out fold operators with the day recipe (`--basis-prune essential`, 12k
cap; 3 × ~5 min, 13 GB each) and recompute the A/B/C fold `hi` and the key-class prices. Read it as follows:
- `hi` stays at 5.2–7.0: the light tail is A/B/C's panel.
- `hi` jumps to ≥ 8: the light tail was 062's operator.

rehearse708b runs the same test tonight as a side effect. Its day708 run fits pruned held operators on 180's relabelled
A/B/C bundle, and score708 prints each fold's `hi`.

## Resolution (lab, 2026-10-16T11:16Z entry; docs/LAB-PANEL-REPORT-2026-10-16.md, 569978b9): the PANEL, not the operator recipe
- The lab rebuilt A/B/C's 300 targets with the exact day recipe: day-recipe operators (essential-pruned basis, 12k cap),
  D/E/F rules, the same held truth lines (**LABC**), and also A/B/C rules with day operators (**ABCp**).
- With every fold on that recipe, A/B/C's panel still ranks at the BOTTOM of its held line's seven folds in 14 of 21
  (arm, held line) cells, where 3.0 are expected. Its mean percentile is 0.19, against 0.50 if exchangeable. The key arm
  (cap 6 − cap 3) gains +2.0 … +4.2 % r on LABC, against +5.0 … +14.0 % on L1–L6. The fold rules account for only +0.4–0.5
  points.
- So the confound raised above does not explain the tail: 062's unpruned basis is not why A/B/C reads light. The light key
  tail is a property of A/B/C's panel on the source side. The lab's 814 is looking for the target property.
- **Consequence.** If Arc drew D/E/F's panel the way it drew A/B/C's, 708's day-of folds on D/E/F should behave like LABC,
  without the fresh-panel over-read. Then the same-panel ratios in 708's information column apply. The lab's fresh ÷
  Arc-drawn factors are: key ×3.08 [1.91, 6.76], direction ×≈6.5 (lower bound 2.0), shape ×1.19.
