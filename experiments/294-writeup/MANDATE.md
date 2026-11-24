# 294 — The scientific write-up: a manuscript draft, figures reproducible from the saved tables, and a claims ledger — "Which perturbation responses are conserved across cell contexts, what a single-cell reference can resolve, and what bounds zero-shot prediction of perturbation responses". Sourced ONLY from the findings note (v1.9), the experiments' RESULTS and their saved CSV / JSON. ≤ 6 hours. No new measurement, no build, no submission.

You are a research subagent. Work in `experiments/294-writeup/`. **Edit nothing outside it.** No build, no `vcc`, no
`cell-eval2`, nothing under `data/`, `submissions/` or `docs/` (the coordinator integrates). `uv run python` for the
figure scripts only; read saved CSV / JSON / npz outputs of finished experiments — **compute no new statistic on raw
data** (a figure may re-plot a saved table; it may not re-derive a number). Peak RSS under ~4 GB. **Time box: 6
hours.** Two siblings run in parallel: 293 (the construct ceiling in the challenge's own regime) and 295 (claim C5 on
reserve R3). Do not read their directories until a `RESULTS.md` exists there; poll every 15 minutes; when either
lands, integrate its one-line verdict and headline numbers into the relevant section; if one has not landed when
you finish, leave an explicit `[[PENDING 293]]` / `[[PENDING 295]]` marker where its number belongs. Never touch
reserve data. `docs/RESEARCH-PROGRAMME-2026-09.md` rule 10 binds every sentence: no claim without its provenance
verdict (NEW / EXTENDED / KNOWN / CONTRADICTED, from `experiments/287-novelty-audit/NOVELTY-AUDIT.md`) and its
confirmation status (`docs/CONFIRMATION-LEDGER.md`).

## Sources of truth (in this order; where they disagree, the later amendment wins)
`docs/FINDINGS-2026-09.md` (v1.9: F1–F14 with every 13:40Z–16:20Z amendment — the amended wording, never the original
headlines); `docs/CONFIRMATION-LEDGER.md` (C1 relabelled; C2 partial; C3 / C4 FAIL with what survived; R3 declared;
C5 pending); `experiments/287-novelty-audit/NOVELTY-AUDIT.md` (the provenance table and the "what can be claimed as
new" tiers, plus the coordinator's addendum re-grading C-1 to RECONCILED); `docs/RESEARCH-PROGRAMME-2026-09.md` (the
organising question, the rules, the Round 4 results); `docs/READ-FIRST.md` §6 (the one-paragraph readings of
267–293). The experiments' `RESULTS.md` and saved tables for the numbers and figures: 267 (`work/tables267.md`,
`truth-<fold>.npz` `rel|*`), 273 (`work/tables273.md`, `reliability273.json`), 280 (`work/q1_*.csv`, `q2_curve_*.csv`),
282 (`work/state282.csv`), 285 (`work/bench285.csv`, `analyse285.json`, `equiv285.json`), 276 / 279 (`conservation276.csv`,
`znf279.csv`), 281 / 284 (`pred284.csv`, `fit284.json`, `fam284.json`, `rem284.json`), 290 / 291 (`pred290.csv`,
`fit290.json`, `surfA291.json`, `surfB291.json`), 268 / 283 / 289 (`decomp268_cpm5.json`, `lolo289.csv`, `p1p2_289.json`),
271 / 272 / 286 (`analysis272.json` or its tables, `analysisB286.json`, `tablesB286.md`), 262 / 265 (the knee table
in FINDINGS F4), 288 (`work/check288.json`). Cite every number to its file.

## What to write (`PAPER.md`, ~6,000–9,000 words of main text; British or American spelling, consistently)
1. **Title, abstract (≤ 250 words), one-paragraph summary for a non-specialist.**
2. **Introduction:** the organising question; what "reproducible", "conserved" and "transferable" mean here and
   why the literature conflates three quantities (287 §2's structural point: within-reference split-cell
   repeatability +0.243, across-donor +0.054, across-context ≈ +0.03 on the same cells); the challenge as the
   motivating setting (zero-shot, one construct per target, 400 cells, six lines; the 25/300 library overlap fact).
3. **Results, in this order, each with its provenance tag and confirmation status stated in the text:**
   (a) *What a single-cell perturbation reference resolves* — F1 (KNOWN, re-measured; the Orion NTC-Spearman trap),
   F14 (one removable axis; Spearman–Brown on real cell ladders; the construct as the dominant irreducible term with
   ρ_guide by dataset and modality; batch free; study 0.79; line 0.35–0.45; donor ≈ guide in CD4), the challenge-regime
   ceiling from 293 `[[PENDING 293]]`, F10b (donor marginal; replicates vs cells at matched budget; the iPSC-34
   tension), the retracted depth claim stated as retracted.
   (b) *The cross-line conservation map* — F8 (EXTENDED: the GWCD4i source table already has a score; new: reproducible
   across disjoint line sets, predictive with the line held out, transfers off-map — 291's +0.274 on the authors'
   CD4-vs-K562 statistic and +0.210 on Broad-16 at low power; 288's magnitude decomposition; C1 as split-half
   robustness, not confirmation), C5 on R3 `[[PENDING 295]]`.
   (c) *Family axis and the per-gene remainder* — F9 (preregistered on 17,106 genes; EXTENDED; corroborated on an
   independent estimator with the interval clearing zero in one of three states), F12 as amended by 284 (the leak
   fixed; +0.805 under strict folds; ZNF 57–80 % predicted; the remainder REPRODUCIBLE across disjoint source sets —
   F12's "mostly measurement" withdrawn).
   (d) *Response geometry* — F13 (290: the remainder is carried by geometry at the ceiling; the preregistered
   mechanism inverted: broad, spectrally central responses are MORE conserved; programme identity adds little) and
   its narrowing by 291 (the sign confirmed on two outside surfaces; "beats annotation" specific to the remainder
   target; magnitude dominant off-map; the DNA-binder exception not reproduced) and C5 `[[PENDING 295]]`.
   (e) *State dependence at the granularity where it exists* — F5 as amended (268: a reproducible, high-rank
   line-specific component; C2 partial with the depth explanation the source paper had published; 289: the paper's
   own state features do not predict the residual transcriptome-wide, yet P1 reproduces and the LOLO finds UBE3A /
   HPV blind, p53-axis sign gating beats chance — C-1 RECONCILED), F11c (NEW: conservation does not predict
   within-individual state dependence; the within-individual size not measurable at ~107 cells), F11b (the donor ×
   guide construct effect, corroborated by the source's own table).
   (f) *Why rank retrieval rewards noise directions* — F3a (textbook, cited as such) and F3b (NEW: damping the
   noise-dominated directions loses retrieval on PLOLO and reads the same sign on the Flex reference), F4 (the knee
   as a property of the decoy distribution), and 286 (the cross-product estimator: the mechanism survives an
   unbiased basis; 80 % of the whitened prediction energy where the signal covariance finds no signal).
   (g) *What bounds zero-shot prediction in the challenge* — the construct ceiling binds every entry equally (no
   source carries the 2026 constructs); the marginal band (F4) against the ceilings (F14); what a reference of one
   construct at 400 cells is a measurement of; what would move it (a second construct per target; the same library;
   depth to 76–330 cells per side only).
4. **Discussion:** what is new (287's Tier 1 with caveats in the same sentence), what is re-measurement credited to
   its sources (Tier 2), what is NOT claimed (F7, C1 as confirmation, F11a's ratio, F4a's knee as discovery, the
   spectral-over-pathway sentence, the depth ceiling), the contradictions and tensions (287 §3: T1 shrinkage
   whitening, T2 iPSC-34, T3 magnitude, T4 depth-conditioned rules, T5 Wang 2026's positive training selection),
   limitations (one map, one statistic, shared pipelines; n of lines, donors, constructs; nothing measured on A/B/C).
5. **Methods:** the response object and residual space; the split units; the conservation statistic (DEFINITION-276);
   the scorer statistic (218's `_per_target` = cell-eval2 to 5.6e-17); the reliability ρ construction (285 §6);
   Spearman–Brown; intervals (paired target bootstrap; block / component bootstraps; permutation nulls); the
   preregistration, known-answer-row and confirmation-ledger practice; multiplicity accounting; data and code
   availability (paths; the reserve rule).
6. **Figures (6–8), each generated by `make_figs.py` from saved tables into `figs/` as PNG + SVG, with a caption that
   names the unit, the n and the interval source:** F-1 the cell ladder vs Spearman–Brown and cells-vs-UMIs
   (285 `bench285.csv`); F-2 the unit decomposition ρ by dataset and modality (285 `analyse285.json` `rho`); F-3 the
   map's reproducibility and off-map transfer (276 / 291); F-4 the family axis under strict folds and the remainder's
   cross-half reproducibility (284); F-5 the geometry sign inside and outside the map (290 / 291); F-6 the
   noise-direction energy and the retrieval loss (272 / 286); optionally F-7 the per-target p53 mechanism (289
   `p1_scores289.csv`) and F-8 the knee (F4's table).
7. **`CLAIMS.md`:** one row per claim in the paper — claim · evidence (file) · unit and n · interval · provenance
   verdict · confirmation status (ledger) · what would falsify it. **`REFERENCES.md`:** from 287's `work/refs.md`
   plus anything new, with access dates.

## Rules of the text
Numbers only with their intervals and units; no σ language (READ-FIRST §3); every negative stated at its
construction's scope; every re-measurement credits its source; quotations ≤ 15 words. Do not invent a result; if
a number you need is not in a saved file, write `[[MISSING: what, where it should be]]`. Keep the competition track
to section 3(g) and the Discussion — no board numbers beyond the public leaderboard facts already in READ-FIRST §1.

## Deliverable
`PAPER.md`, `CLAIMS.md`, `REFERENCES.md`, `make_figs.py`, `figs/`, `RESULTS.md` (what was written, what is pending,
the figure list, the word count, every `[[PENDING]]` / `[[MISSING]]` marker listed). Then stop.
