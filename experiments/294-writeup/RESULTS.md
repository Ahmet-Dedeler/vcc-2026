# 294 — RESULTS: the scientific write-up

**Delivered: `PAPER.md` (manuscript draft), `CLAIMS.md` (50 claim rows), `REFERENCES.md`, `make_figs.py`, `figs/` (8 figures, PNG + SVG + `MANIFEST.md`).** No measurement was made: every number is read from a finished experiment's saved CSV / JSON / Markdown table, and the figure script re-plots those tables without deriving a new statistic. No build, no `vcc`, no `cell-eval2`, no download, nothing written outside `experiments/294-writeup/`. Peak RSS of the figure script ≈ 0.3 GiB.

**Both sibling experiments landed inside the time box and both are integrated in full** — there is no `[[PENDING]]` marker left anywhere in the deliverables.

---

## 1. What was written

`PAPER.md`, structured as the mandate specified:

| section | content | words |
|---|---|---:|
| Title, abstract, plain-language summary | abstract **240 words** (cap 250); one-paragraph non-specialist summary, 323 words | 579 |
| **1. Introduction** | the organising question; the three quantities the literature conflates (+0.243 / +0.054 / ≈ +0.03 on the same panel); the challenge as the motivating setting, including the 25/300 library-overlap fact | 653 |
| **2. Data and instruments** | the four measurement surfaces; the three reserves, all now spent | 202 |
| **3. Results** (a)–(g) | as ordered in the mandate; every claim carries its provenance verdict and confirmation status in the text | 7,355 |
| **4. Discussion** | 4.1 what is new (287 Tier 1 + three post-audit items) · 4.2 re-measurement credited to its sources (Tier 2) · **4.3 an amendment this synthesis makes to its own record** (F14 → guide efficacy, forced by 293) · 4.4 what is not claimed · 4.5 tensions T1–T5 · 4.6 limitations | 1,361 |
| **5. Methods** | the response object and scored statistic; split units; DEFINITION-276; the reliability ratio ρ (285 §6); intervals, permutations, multiplicity, and the explicit no-σ rule; preregistration / known-answer / ledger practice; leakage rules; data and code availability with the reserve rule | 1,333 |
| Figures + pointers | the 8-figure table with source files | 32 (plus the table rows) |

**Main text (sections 1–4): 9,571 words** by the convention "strip markdown, drop table rows, count whitespace-separated tokens, count each bracketed interval bound as a token" (`[+0.326,` and `+0.444])` count as two). That is **~6 % above the mandate's `~6,000–9,000`**. Four compression passes took it from 9,686 down; the residual overshoot is the ~800 words that 293 and 295 added when they landed at 16:30–16:36Z and had to be integrated (§3(a), §3(b), §3(d), §3(g), §4.1, §4.3, and two CLAIMS rows). I stopped trimming rather than delete measured content; the sections that could still be cut without losing a number are §4.2 and §4.6.

Spelling: British, consistently. No number is expressed in units of σ (READ-FIRST §3); Methods §5.5 states the rule and records the one "8 σ" cell that was re-measured at z = +2.3. No quotation exceeds 15 words (checked mechanically). Every negative is stated at its construction's scope. Board numbers appear only in §3(g) and are confined to the public-leaderboard facts in READ-FIRST §1.

## 2. What the sibling experiments changed

Both landed while this was being written and both are integrated with their headline numbers and their verdicts:

**295 (claim C5 on reserve R3) — C5a PASS, C5b PASS.** Integrated into §3(b) (a full paragraph), §3(d) (a full paragraph), §4.1, and `CLAIMS.md` C-47. Headline as carried in the paper: on eight untouched scPerturb datasets, the 28-bank conservation map predicts a gene's actual agreement with the map's consensus at **ρ +0.5016** (file-block [+0.1525, +0.5659], 8 files as blocks; 0/200 shuffles, z +6.9); partial given the file's reliability and the map's magnitude **+0.3086 [+0.1389, +0.3964]**; within-file ranks **+0.3508 [+0.0773, +0.4852]**; F13's sign travels (`gen_pc20` **+0.3709 [+0.0635, +0.4148]**). Limits carried with it: 8 blocks, one file holds 78 of 174 genes, and the margin over annotation (+0.0797 [−0.0747, +0.1767]) crosses zero — only the margin over magnitude separates.

**293 (the construct ceiling in the challenge's own regime) — and it forced an amendment.** Integrated into §3(a) (two paragraphs), §3(g) (rewritten), the abstract, §4.3 (a new sub-section), and `CLAIMS.md` C-46, C-49, C-50. Headline: **ρ_guide +0.075 [−0.010, +0.288]** in Arc's own 2025 H1 library at challenge depth; SB to 400 cells gives same-construct **+0.476 [+0.260, +0.717]** against different-construct **+0.046 [−0.006, +0.121]**; the zero-shot cosine ceiling is **+0.213 [+0.000, +0.348]**, knockdown-matched ≈ **+0.572**, against **+0.690 [+0.510, +0.847]** for a construct-aware model.

Two things in the draft changed because of 293, and both are recorded as corrections rather than quietly edited:

1. **"The construct ceiling binds every zero-shot entry equally" was too strong.** It is a real, quantified ceiling (cos 0.21–0.57) with the incumbent's median per-target transfer cosine ≈ 0.04 sitting about five times below even its pessimistic end. §3(g) now says so.
2. **F14's "the guide construct is the dominant irreducible term" needs an amendment.** In the challenge's own library the term is almost entirely *knockdown failure*: ρ_guide −0.003 [−0.085, +0.181] where the weaker construct misses |lfc| 1.0, +0.687 [+0.197, +1.067] where both reach it; and two-construct targets knock down significantly worse than the rest of the library (p = 0.0015), so +0.075 is a lower bound. §4.3 states the safe wording ("**guide efficacy**, of which construct identity is the visible proxy") and names the untested next measurement — 293's stratification applied to 285's five public libraries.

Also carried from 293: the **back-of-envelope in F14** ("ρ_guide × ρ_line is of the order of the incumbent's own cosine") conflates a variance share with a cosine and is superseded by 293's direct measurement. §3(g) records it as an error worth naming.

## 3. Figures

All eight are regenerated by `uv run python make_figs.py` from saved tables into `figs/` as PNG + SVG, with full captions (unit, n, interval source) in `figs/MANIFEST.md`.

| # | what it shows | source tables |
|---|---|---|
| **F-1** | (a) the cell ladder against Spearman–Brown on four real ladders; (b) cells and UMI depth as one axis (thinned points on the cell ladders) | `E285/work/bench285.csv`, `analyse285.json` |
| **F-2** | the unit decomposition: ρ by dataset and modality, grouped by unit, n per row, intervals where stored | `E285/work/analyse285.json` (`rho`), `equiv285.json` (`frac`) |
| **F-3** | (a) the map's reproducibility (independent cells; disjoint line sets; 50 random splits); (b) its prediction of transfer on PLOLO, Flex-400, and two off-map surfaces | `E276/work/{analyse276,rel276,predict276}.json`; `E291/work/{surfA291,surfB291}.json` |
| **F-4** | (a) the family axis under strict folds (78 families, observed vs out-of-fold predicted); (b) the remainder's cross-half reproducibility over 20 partitions against the raw-map ceiling and the permutation null | `E284/work/{fam284,rem284}.json` |
| **F-5** | (a) the geometry features against the other half's remainder, with the two preregistered-negative signs marked; (b) the same features and two ridge models on both off-map surfaces | `E290/work/fit290.json`; `E291/work/{surfA291,surfB291}.json` |
| **F-6** | (a) where the whitened prediction's energy sits on four operators; (b) the retrieval loss from damping those directions, with shuffled and random-direction controls | `E272/work/analysis272.json`; `E286/work/analysisB286.json` |
| **F-7** | the per-target p53 mechanism: the two HPV-E6 lines against the other ten, with the E2F specificity check | `E289/work/p1_scores289.csv` |
| **F-8** | `pds` as a threshold in the transfer cosine at ≈ 0.02, on both references, with the decoy counts that set the two slopes | `E273/work/tables273.md` §C |

Design notes: colour-blind-safe palette; every panel labels its unit and n; no panel plots a quantity the source did not store. Two figures were rebuilt after a first-pass inspection found real defects — F-2's K562 guide row was reading the 8-target cohort (ρ −0.78, uninformative) instead of the 748-perturbation cohort the findings note quotes, and F-8's table parser was matching the section header rather than the table body and produced an empty plot. Both were fixed and the script now asserts the F-8 parse returns ≥ 6 bins.

## 4. Markers left in the deliverables

**`[[PENDING]]`: none.** Both siblings landed at 16:30Z (295) and 16:36Z (293) and were integrated.

**`[[MISSING]]`: three, all the same gap** — the 287 novelty audit covers F1–F12, C1 and C2 and was written at 14:07Z, before 285 (F14, 15:52Z), 290/291 (F13, 15:04/15:30Z) and 284's remainder test landed, so those items have no NEW / EXTENDED / KNOWN / CONTRADICTED verdict to carry under programme rule 10:

| # | location in `PAPER.md` | text |
|---|---|---|
| 1 | §3(a), the F14 heading line | `[[MISSING: a 287-style provenance verdict for F14; the novelty audit covers F1–F12, C1 and C2 and predates 285 by ~1.75 h]]` |
| 2 | §3(d), the F13 heading line | `[[MISSING: a 287-style provenance verdict for F13 (290/291); the audit covers F1–F12, C1 and C2 and predates 290 by ~1 h]]` |
| 3 | §4.1 | `[[MISSING: 287-style provenance verdicts for F14, the 284 remainder test and F13 — the audit was written before these experiments landed]]` |

`CLAIMS.md` carries the same gap as `n/a (post-audit)` in the provenance column of 13 rows (C-06…C-11, C-26…C-31, C-41, C-46, C-47, C-49, C-50). **The fix is a short extension of the 287 audit to F13, F14, 284's remainder test and 293's amendment** — a reading task, not a measurement, and the one thing this write-up needed and did not have.

No `[[MISSING]]` was needed for a number: every figure and every claim resolved to a saved file.

## 5. What is claimed, in one table

`CLAIMS.md` carries 50 rows, each with claim · evidence file · unit and n · interval · provenance verdict · confirmation status · what would falsify it. The confirmation summary:

- **preregistered and confirmed by its own frozen rule**: F9's family axis (C-21); C5a and C5b on untouched reserve data (C-47).
- **partially confirmed**: F5 on the reserved Broad-16 lines — C2a FAIL, C2b PASS, C2c PASS (C-32, C-33, C-35).
- **preregistered and not confirmed, reported as such**: F13's mechanism (sign inverted, C-28); C3 and C4 (C-30); F11c (C-36); the retracted depth ceiling (C-14).
- **relabelled**: C1, from confirmation to split-half robustness (C-20).
- **explicitly not claimed**: the depth ceiling, the DNA-binder exception, F11a's ratio, F4a's knee as a discovery, the board calibration constants, and F14's superseded back-of-envelope.

## 6. Verification done on the write-up itself

- Every file path cited in `PAPER.md`, `CLAIMS.md` and `REFERENCES.md` was checked to exist (one correction: 267's per-target arrays are `pertarget267-<fold>.npz`, not the `truth-<fold>.npz` the mandate named).
- Every number in the figures is read from its source file at plot time, so a figure cannot drift from the table behind it.
- Two numeric slips found and fixed against the saved JSON: the per-fold n's for 288's magnitude decomposition (268 / 289 / **293**, not 294), and the CD4-surface `breadth` upper bound (+0.156, not +0.155).
- Mechanical checks: no "σ" outside the Methods disclaimer; no quotation over 15 words; abstract 240 words; section cross-references normalised to `§3(a)`-style.

## 7. Files

`PAPER.md` (76 KB) · `CLAIMS.md` (34 KB, 50 rows) · `REFERENCES.md` (16 KB; 287's list plus the challenge manifest, the annotation resources and the 293 / 295 artefacts) · `make_figs.py` (`uv run python make_figs.py`, ~40 KB, regenerates everything) · `figs/F-1…F-8.{png,svg}` + `figs/MANIFEST.md`.
