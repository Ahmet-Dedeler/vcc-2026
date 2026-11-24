# 264 — The three-donor MAJORITY rule on the CANCER folds, with `kolf2` as the third donor: does 263's one-fold signal hold on cancer held lines? Measurement only, ≤ 3 hours. DO NOT BUILD, DO NOT SUBMIT.

You are a research subagent. Work in `experiments/264-majority-rule-cancer-folds/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. One PLOLO fold per process, peak RSS
under ~8 GB. **Time box: write up after ~3 hours.**

## Why (R263 §7–§8; R262 §5)
`cd4svd5` loses `pds` where the whitened donor CONSENSUS points away from the held line (R262). 263 measured the
per-target oracle at +0.046 `pds` and found no observable arbiter that captures any of it — except one cell on the one
fold that has three donors: on `kolf2` (donors `k562_gwps_pm50`, `orion_hct116`, `orion_hek293t`, the production
operator) the **majority rule** — for each target carried by all three, keep the two un-whitened rows whose cosine on
the operator columns is highest and drop the odd one out, norm-matched to the three-row mean — reads **+0.0121 `pds`**
(240 targets, +2.4 % `cos_agg`, delivered 0.9986, leave-top-3-out 60 %, reverse −0.0097, random −0.0049, not
reproducible by any global donor weight). It is production-realisable (A/B/C see all three donors) and priced
+0.0051 pessimistic / +0.0076 point, but one fold on an iPSC line cannot clear the 3/3-cancer-fold rule. **Question:
with `kolf2` added as a THIRD donor on each cancer fold, does the majority rule buy `pds` on 3/3 cancer folds, against
the three-donor mean, with its controls below half?** If yes, a slot on it is a preregisterable information purchase;
if no, the `kolf2` reading was line-specific and H2 is closed entirely.

## READ FIRST
`docs/READ-FIRST.md` (v3.2), 263's `RESULTS.md` §1, §7, §8 and `lib263.py` / `sweep263.py` / `kolfctrl263.py` (REUSE the
selection path, the norm matching, the gating, the controls; byte-copy the fold caches), `experiments/162-retrodicting-validator/plolo.py`
(`build_fold`: how donors are prepared per fold via `E062.prepare_main_source`; `kolf2` lives at `data/external/kolf2`
and is prepared for the `kolf2` fold as a TRUTH — prepare it as a DONOR the same way, then confirm its prepared rows
are the object 263's `kolf2` fold used, hashed), `experiments/262-where-pds-is-lost/RESULTS.md` §5, `AGENTS.md` §4.2
(donor weights), §8.1b, §11.1 rule 10.

## Construction (three cancer folds; `kolf2` fold reported apart as the known-answer)
For each cancer fold: donors = the fold's two production donors + `kolf2` as a third, all prepared as `build_fold`
prepares donors, on the fold's own held operator (state this: the operator was fit on the two production donors; do
NOT re-fit it; report as a limitation), then: (a) **3-donor equal mean** — the new reference base; (b) **majority** —
per target carried by all three, keep the two agreeing raw rows, norm-matched to the 3-row mean; (c) 2-donor mean ==
`cd4svd5` (263 R1, the incumbent). Residuals as `cd4svd5` (CD4 deflated k = 5 re-solved at f 0.34; H1 pair at g 0.28).
Read: `pds`, `cos_agg`, macro, per-target `pds`, delivered ratio, `nmae89` by the identity, for (a), (b), and (b) − (a)
and (b) − (c). Gating: `all` and τ = 0 as 263. Controls per fold: reverse (drop the agreeing pair's odd… i.e. keep the
two LEAST agreeing), random pair (3 seeds), always-drop-`kolf2` (= the 2-donor mean = `cd4svd5`), always-drop each
production donor. Leave-top-3-out retention and the top-3 target share for (b) − (a) and (b) − (c).

**Known-answer rows first:** (R1) with the two production donors only, `cd4svd5` reproduces 263 R1 at 0.0e+00;
(R2) on the `kolf2` fold, 263's `majority:raw:all:nm` reproduces +0.0121 / +2.4 % / 0.9986 / retention 60 % through
YOUR code path; (R3) the prepared `kolf2` donor rows hash to what 263's `kolf2` fold consumed as truth-side rows (or
state the exact difference: e.g. panel-column zeroing).

## Price (READ-FIRST §3; three coherent scenarios; class by delivered amplitude). CLEARED = (b) − (c) pessimistic
≥ +0.0019 on the three cancer folds AND `pds` > 0 on 3/3 at ≥ 2 LOO-σ AND retention ≥ 50 % AND every control below
half. ALSO report (b) − (a): whether the rule beats the three-donor mean (the mechanism) separately from whether adding
`kolf2` as a donor helps at all (a confound — say which of the two carries any gain).

## Deliverable
`RESULTS.md`: R1–R3; the per-fold table for (a), (b), (c) and the controls; the two increments with all three
scenarios; retention; the exact production construction (three production donors, majority, norm-matched — no
`kolf2` in production — and the builder change 263 §8 named, `VCC_DONOR_SELECT=majority:all` at
`build_v9_pre_whiten_effects`, not made); one line **CLEARS / DOES NOT CLEAR**; if it clears,
`PREREGISTRATION-DRAFT.md` in 249's form with five falsifiers priced on the three scenarios, and the list of
known-answer rows a build would need (flag-off byte identity with `aba27758…`; the `kolf2` row). No model name, no
description, no submission. Then stop.
