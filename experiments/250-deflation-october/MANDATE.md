# 250 — Does the CD4 deflation TRANSFER to a disjoint panel, and where is the k knee there? The October reading of 246/249, before the board's.

You are a research subagent. Work ONLY in `experiments/250-deflation-october/`. Read-only everywhere else.
**No `vcc`, no build, no `cell-eval2`, no production change, no slot proposal.** `uv run python`. Peak RSS
≤ 10 GB; **MEMORY GATE before every fold load:** `memory_pressure` ≥ 40% free AND `vm_stat` free+inactive
≥ 20 GiB (248 is running builds); strictly serial; `/usr/bin/time -l` on every heavy script.

## Why
`cd4svd5` (249; `docs/LOG.md` 17:45) deflates the shipped CD4 residual's top-5 target-shared directions and
is armed for a 2026-09-05 slot. Its `pds` gain (+0.01226 raw, 3/3 cancer folds, 246) and its k = 5 were
measured on the A/B/C panel's PLOLO folds. **The only entry that counts is on a different 300-target panel
(Oct 22), where the CD4 residual is re-extracted and its shared subspace re-fitted on the day.** 246's
split-half projector test says the subspace is a property of the source, not of the rows — a strong hint
it transfers — but nothing has measured the deflation on a panel disjoint from A/B/C, and the k knee
(3–5 on A/B/C, saturating by 10) may sit elsewhere when the rows change. 233 already holds everything needed:
per-panel operators (`work/ops`), per-panel CD4 re-extractions (`work/cd4`), and (panel × cancer-fold)
caches (`work/folds`) for `abc`, `real` (§5 realistic draw) and `pess` (§5 pessimistic draw), plus its
instrument (`fold233.py`, `score233.py`, `retain233.py`).

## READ FIRST
`docs/READ-FIRST.md`; 246 RESULTS §3.3–§3.7 and `lib246.py` (`deflate`, `merged_many`); 249 RESULTS §1
(the production construction: deflate BEFORE `solve_cd4_dose`, resolved variant) and its `stage249.py`;
233 RESULTS §0–§3 and its `panels233.py`, `fold233.py`, `score233.py`; 238 RESULTS §1 (the essential-free
panel and why `pess`'s 30 essential targets matter); 247 RESULTS §3 (CD4 carrying-row sizing on non-A/B/C
panels: `g_cd4` 0.3443655278649641 — use it on the disjoint panels, since that is what October ships);
`AGENTS.md` §11.1 rule 10, §14.

## Known-answer rows (abort on failure)
1. On 233's `abc` folds, the shipped stack (`cd4@0.34` all-rows rule + H1 union `g` 0.28 carrying rule)
   reproduces 233's own baseline rows for that panel to the digit, and deflation k = 5 (resolved)
   reproduces 246's `B2_cd4_svd5_resolved` three-cancer-fold numbers (macro −0.808675%, aggregate
   −0.461748%, Δ`pds` +0.012263) — 249's `repro` did this through the builder; do it through 233's
   instrument and state any difference between the two instruments.
2. 233's `pess` and `real` baseline retention rows reproduce (its `analyse233.json`).

## Work
### Phase 1 — the k ladder on the disjoint panels
For panel ∈ {`real`, `pess`} (and 238's essential-free `ef` panel if 233's fold builder can take its bundle
within budget — say so either way), on the three cancer folds, on top of that panel's shipped stack (CD4
sized by 247's carrying-row rule at `g_cd4`; H1 union at 0.28 carrying-row; operator = 233's per-panel
fit): deflation k ∈ {0, 2, 3, 5, 7, 10, 15}, resolved sizing. Endpoints: `cos_agg`, macro, `pds`, `q_w`,
energy removed. Report the knee per panel and whether **k = 5 sits inside the plateau** (Δ`pds` within
10% of the k-max) on each. `kolf2` apart, as always.
### Phase 2 — the controls that made 246's claim (each at k = 5 on each disjoint panel)
Energy-matched random subspace (two seeds); shuffled projector; **split-half projector** (basis fitted on
half of CD4's carried rows on that panel, applied to the other half); non-carrying rows unchanged but
for the preserve scalar; leave-top-3-out retention on `pds` and cosine separately.
### Phase 3 — the projector's cross-panel identity
Fit `U_5` on the A/B/C CD4 residual and on each disjoint panel's; report the principal angles between
the two subspaces and the fraction of each panel's top-5 energy captured by the other panel's basis.
If they are nearly the same subspace, October could pin the A/B/C basis (a known object) instead of
re-fitting on the day — say which is safer and why, with the numbers.
### Phase 4 — the `+partial` composition (246 §3.7) on the disjoint panels
CD4's partial-condition rows, carrying-row sized at 0.34, composed with k = 5: Δ`pds`, Δcosine, LT3O.

## Deliverable
`RESULTS.md`: the known-answer rows, the k-ladder table per panel, every control's number, the principal-
angle table, the composition rows, and one line per panel: **TRANSFERS / DOES NOT**, with the k to ship on
Oct 22 and the falsifier the D/E/F census would fire (e.g. a top-5 energy fraction outside the A/B/C
range 0.19–0.21). Propose the `DEF-RUNBOOK` §3 decision-table row text for `VCC_CD4_DEFLATE_K`; do not edit
the runbook. Then stop.
