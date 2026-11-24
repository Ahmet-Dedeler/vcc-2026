# 348-reach-ordering — PLAN

## Hypothesis
`reach` ranks the reference's significant genes by OUR confidence (emitted
significance: significant-first, then p_adj, p, |lfc|). Our ordering is
effectively by designed effect SIZE. If sign RELIABILITY concentrates on a
design-time-identifiable subset (e.g. predicted-down genes — 343 found
reference-up genes are a coin flip at 0.516 vs 0.588 overall — or genes
agreeing with the cross-target common vector), ranking that subset first
should keep the 90%-pure prefix P(k) ≥ 0.9 deeper and raise raw reach.

## Inputs (all read-only, all small — no 323 scoring run)
- `experiments/323-current-family-validation/work/stats-cap15.npz`
  (delivered p/q/lfc, 85×10,779 universe)
- `experiments/323-current-family-validation/work/reference.npz`
  (genes, targets, keep mask, ref_lin)
- `experiments/323-current-family-validation/work/reference_de.parquet`
  (truth DE long table)
- `experiments/323-current-family-validation/work/per-target-cap15.csv`
  (known answer: per-target `de_wilcoxon_direction_reach_raw`)
- `experiments/323-current-family-validation/work/metrics-cap15.json`
  (known answer: aggregate raw reach 0.25434726355678233)
- `experiments/323-current-family-validation/work/effects.npz`
  (`cap15` design effect matrix + `genes` axis, for orderings 4/5/7)
- Scorer source: installed `cell_eval2/metrics/direction.py`
  (`de_direction_reach`, universe='adjudicated', corrected=False, P0=0.9)
  and `evaluate323.py:score()` (frame construction).

## Estimator
1. **Known answer.** Rebuild pred/real frames exactly as `evaluate323.score()`
   does and run the installed scorer (`prepare_de` + `de_direction_reach`
   corrected=False). Then implement an independent polars walk of §5
   (pool = ref-sig genes minus own target gene; rank = sig-first, p_adj↑,
   p↑, |lfc|↓, gene↑; adjudicable-only depth; k* = deepest k with P(k)≥0.9
   else 0; R = k*/N_conf; NaN iff pool empty). Require per-target match to
   the csv to 1e-9 (all 85, incl. NaN positions) before any new measurement.
2. **Alternative orderings of the SAME delivered signs** (pool, signs,
   N_conf fixed; only the rank key changes):
   - O1 current (known answer); O2 random permutation (100 seeds, mean±sd
     of the panel mean); O3 |delivered lfc| desc; O4 |design effect| desc
     (mapped by gene name); O5 agree-with-common-vector first
     (common = mean design effect over 85 targets; agree iff
     design*common > 0), then |design| desc; O6 oracle correct-sign first
     (ceiling); O7 own design-time keys: (a) predicted-down first by |lfc|,
     then predicted-up by |lfc| (343-motivated); (b) disagree-first
     (reverse of O5, diagnostic); (c) significant-first then |lfc| (ablation
     of the p-value refinement inside the current key).
3. **P(k) profile** for O1: mean P_t(k) over targets with adjudicable depth
   ≥ k at k = 5,10,20,50,100,200 (+ n contributing); fraction k*=0, k*<10;
   first-position-miss diagnosis (k*=0 ⟺ pos-1 miss AND no recovery —
   verify empirically; tabulate recovery depths for pos-1 misses).
4. Panel mean throughout = nanmean over scored targets (NaN set is
   reference-determined, identical across orderings — verify).

## Controls / falsifiers
- O2 (random) is the floor; O6 (oracle) is the ceiling. Any design-time
  ordering must sit between them; if O6 ≈ O1, signs (not order) bind reach.
- O5 vs O7b (opposite directions of the same key): if both beat O1, the key
  is noise.
- Gene-name tiebreak determinism: all keys end in gene-asc; random uses
  seeded rng (report seed scheme).
- If the 1e-9 known-answer match fails: stop, report, do not interpret
  ordering deltas.

## Kill criterion (from mandate)
No design-time ordering beats current by ≥ +0.03 raw reach (panel mean) →
ordering lever CLOSED. If one does: state exactly which emitter quantity
(p-values via designed z in `rc_emit.rc_hybrid_block`, key/z_cap) must
change, and estimate fid/jac side effects analytically (pure within-call-set
re-ranking leaves both invariant — they read only set+signs; crossing the
significance boundary changes the set).

## Budget / footprint
< 30 min wall, < 2 GiB RSS (largest object: 85×10,779 float64 ≈ 7 MB × ~6).
No downloads, no cloud, no reserve data. Work in
`experiments/348-reach-ordering/` only: `reach348.py`, `work/*.json/csv`.
