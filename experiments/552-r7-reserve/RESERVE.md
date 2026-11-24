# Reserve R7 — declaration (552, session "Explore contribution areas", 2026-09-30, before any draw or build)

Requested by the coordinator (2026-09-30 ~21:10Z): a fresh sealed reserve for confirming the FINAL October recipe, bulk and
DE switches, once, after R6's last access (C10b). Declared here before any R7a panel is drawn and before any R7b cell is
selected; recorded in `docs/CONFIRMATION-LEDGER.md` under "Reserve R7".

## 1. Why R7

R6 (534) holds one more access (C10b, the last bulk step). It prices the bulk channel only: the DE switches still being
decided (call-propensity block, uncalled-gene hook, call-set size K) have been priced off-board only on H1, which has chosen
many settings and is tuned-on (AGENTS §11.1c.5). R7 gives the final recipe one untouched reading in both channels.

## 2. R7a — four Arc-like bulk panels (the bulk channel; 534 / 551's instrument)

- **Rules:** `draw552.py` = 534's `draw534.py` (bundle238's rules) verbatim, with these seeds and names, drawn in this order:

  | panel | rule (238) | pool | seed |
  |---|---|---|---|
  | T1 | efr | K562 ∩ HCT116 ∩ HEK293T ∩ CD4 atlas | 552101 |
  | T2 | efp | K562 ∪ HCT116 ∪ HEK293T | 552201 |
  | T3 | efr | as T1 | 552102 |
  | T4 | efp | as T2 | 552202 |

- **Exclusion:** every target of every `pert_counts.csv` under `experiments/` and `data/` at draw time (R6's four bundles
  included, asserted), each earlier R7a draw, and the K562 essential list. So R7a is disjoint from every tuned-on panel, from
  R6, and internally.
- **Built now:** the draw and the four bundles (target lists; A/B/C's control cells relabelled D/E/F, stored as APFS clones of
  199's files, as R6). Each draw's list and sha256 go to `work/draws.json` at draw time. No operator is fitted, no fold is
  built and no score is computed until the access.
- **Known answers at access:** on `efp`, as R6 (`551/drive551.py ka`'s rows).
- **Limitation, stated now:** R7a's truths are PLOLO folds of the same three donor lines as R6 (K562, HCT116, HEK293T). It
  is fresh in TARGETS, not in lines.

## 3. R7b — one cell-level DE surface in a cancer line: Sunshine et al. 2023, Calu-3 (lung adenocarcinoma)

- **Source:** `data/external/scperturb/SunshineHein2023.h5ad`, sha256
  `4dd14fb68dc16057e1cf4a05d0bcd2efd26c1a9e4ef8c7d21afccaf2cfc84270` (90,380 cells × 33,538 genes, raw counts, CSC; 8 gem
  groups). CRISPRi (scPerturb labels it `CRISPR-cas9`; 225 and the paper: CRISPRi); 1–2 guides per target label.
- **Caveat, stated now: the cells are SARS-CoV-2-infected** (scPerturb `disease` = "lung adenocarcinoma and SARS-CoV-2";
  the same-line registry classes it as an infection screen whose mean response is the stimulation's, not the line's). Every
  cell, controls included, is in the infected state, so R7b tests the DE recipe on an unseen line in a STIMULATED state. That
  is a harder test than D/E/F (unstimulated), and a null or negative R7b read is weaker evidence against a switch than a
  positive read is for it. The coordinator weighs this before the claim is preregistered.
- **Why this line:** among the cell-level perturbation sets on disk it is the only cancer-line CRISPRi screen with Arc-like
  depth (median 17,626 UMI per cell; Arc 20,000) and enough targets (184 single-target labels; 140 with ≥ 100 cells) that no
  DE design has used. Checked: 528, 430, 516, 517, 521 and 537 never read it; 430/516/517/521 used RPE1, kolf2, Orion, Nadig
  and K562. Rejected: XuCao 2023 HEK293 (median 2,155 UMI; HEK293 ≈ a donor line), Frangieh 2021 melanoma (Cas9 knockout in
  co-culture, immune targets), Myllia / Altair A549 and Jiang's five lines (used as banks, and A549 is a likely October line
  whose bank would become a production input), RPE1 (used by four DE experiments), and the drug screens.
- **Prior touches, disclosed:** 225 (bulk donor source sweep: its 6 A/B/C-panel rows priced as a union member, F1 FAIL,
  closed); 284 / 290 (science-track conservation maps: its pseudobulk responses in an 11-line, 87-pair analysis); 451 (the
  same-line registry lists it as the Calu-3 bank). No DE design, sign prior, fold, call-set rule or emission setting was
  ever fitted, selected or scored on it, and its Wilcoxon DE statistics have never been computed in this repo.
- **Construction (built before the access; recorded, with sha256s, before any statistic):**
  - cells: those with exactly one target label (1–2 guides of the same target) and the 3,536
    non-targeting controls; unassigned cells (nperts 0) and multi-target cells (labels joining two targets) dropped;
  - targets: ≥ 100 such cells; carried by ≥ 2 of `k562_gwps`, `orion_hct116`, `orion_hek293t`; not in
    `data/donors/k562_essential/perts.txt` (Arc's panel is essential-free); on the official gene axis. By metadata:
    ~106–118 targets;
  - genes: mapped to the official 18,533-gene axis by symbol (`data/challenge/gene_names.csv` order); an unmapped official
    gene is all-zero; no cell downsampling; depth as measured;
  - scoring: `cell-eval2` with the `vcc2026` preset and its own baseline and 5-split replicate anchors
    (`prep-real-bundle`), exactly as `h1-test`.
- **Sealed:** until the access, nobody computes a DE statistic, a per-target response or a pseudobulk on these cells, runs
  the emitter or the builder on them, or reads the baseline / replicate scores. Metadata (obs labels, per-cell totals, cell
  counts) only.
- **Known answers at access:** the emitter harness first reproduces a recorded `h1-test` row of the arm family it prices
  (as 537's `ka537.py` does), bit for bit; abort otherwise.
- **Void condition:** if an Oct 22 context matches Calu-3 before R7b has been accessed, Sunshine becomes a production input
  (the same-line bank) and R7b is void. An access made before that stands.

## 4. Access rule

- **One claim, once:** the final October recipe (its bulk configuration on R7a; its DE configuration on R7b) against the
  recipe before the switches decided after C10a, after R6's last access (C10b) and after the recipe's DE switches are decided
  (537's tree: 10-07 … 10-10). The claim — designs, statistics, reading rules and prior predictions — is preregistered and its
  sha256 is recorded in the ledger before the first R7 number exists.
- **No fitting, selection, dose search or debugging may read an R7 score.** R7 confirms; it does not choose.
- Every access is appended to the ledger's R7 section with the claim it served.

## 5. Amendment 1 (2026-09-30 ~21:55Z, coordinator's conditions for keeping R7b; written before the surface is built and before
any perturbation response is looked at)

1. **Uninfected cells, if they can be told apart.** If the Sunshine data carry a mock-infected condition, or cells that a
   viral-read threshold fixed in advance calls uninfected, R7b uses only those. **Checked from metadata only (obs columns and
   gene names; no expression read): neither exists in this file.** `obs` has no condition column and `disease` is the same
   for all 90,380 cells; all 33,538 features are human Ensembl genes (`ENSG…`), so there are no viral reads to threshold.
   **The infected state stands, as declared in §3.** The raw GEO data are not fetched for this: viral load is the screen's
   own readout (host factors change infection), so a per-cell infection filter would select cells differently per target.
2. **Asymmetric reading rule.** R7b can VETO a DE switch only on a preregistered HARM beyond its own noise band, derived
   from the surface's own re-rolls (emission seeds and split-half reference re-draws, fixed in the claim's preregistration
   before any R7b score). A null or small effect never counts against a switch; a positive read is robustness evidence
   only. The claim's preregistration states the HARM bin in overall-equivalent units and the band's construction.
3. **Timing.** The surface is built this week outside the coordinator's windows (not 10-04 12:00Z – 10-05 18:00Z). The access
   happens only after the DE switches settle, on the coordinator's message, not before the 10-07 reads.

## 6. Amendment 2 (2026-09-30 ~22:00Z, before the surface is built): what "sealed" means for the prediction side

The prediction for R7b is built the way the October day builds D/E/F's: from the line's CONTROL cells and public donors only.
So these are allowed before the access, and are not reads of the reserve: any statistic of the 3,536 non-targeting control
cells alone (b_ctrl, control CPM and the CPM gate, library sizes, the emitter's sorted control block, the call-propensity
key, NTC-proxy claims), the target LIST, per-target cell COUNTS (metadata), and emitting predicted counts for the targets.
Still sealed until the access: every statistic that involves a perturbed (target-labelled) cell's expression — DE tests,
fold changes, pseudobulks, per-target responses, the cell-eval2 reference bundle (baseline and replicate anchors) and any
score. The truth file `work/r7b/truth_r7b.h5ad` is written by `build_r7b.py` (construction: counts copied, no statistic) and
is opened again only by the access script.

## 7. Amendment 3 (2026-10-01 ~08:45Z, before any target cell's library size is read): a human-library depth floor

553 found that R7b as built cannot be emitted by the production emitter (`553/RESULTS.md`, "R7b: the production emitter
cannot emit this surface as built"): R7b's control libraries are bimodal (a shallow third at 0.5k–3k human-gene counts,
plausibly heavily infected cells whose reads are viral), q75/q25 7.24 against 1.57–1.95 for A/B/C/H1, so the K-class
placement fails on 107 / 107 targets, and under any class rule the steering dose (CLEAN4 prior, MDE, hook) is 0. A floor
restores the unmodified production emitter; 553's sweep (controls only + an h1-test rehearsal) picks 5,000.
- **The rule:** keep a cell iff its raw counts summed over the 18,533 official-gene columns are ≥ **5,000**, for the
  non-targeting controls AND the target cells alike (the prediction side's control file is the same rule on the same
  cells: 2,378 of 3,536 controls, `553/work/r7bd5000/context_R_5000.h5ad`).
- **BST2 is excluded** (it crashes the emitter's placement at every floor; an emission-side edge case, decided from the
  prediction side). 106 targets.
- **After the floor, a target with fewer than 30 truth cells is dropped** (fixed now, before any count is known).
- **What this reads from perturbed cells, and why it is allowed:** one scalar per target cell (its library size on the
  official axis), used only to apply the fixed threshold above, plus the resulting per-target cell counts (recorded). No
  per-gene statistic, response, DE test or score. The threshold was chosen from the controls and the emitter's behaviour,
  never from target cells. It is a QC step of the kind every real surface goes through.
- **Its cost, stated now:** the floor selects cells by a quantity the infection screen itself moves, so per-target truths
  are those of less-infected cells (closer to D/E/F's uninfected state, but a selection that can differ by target), and
  per-target cell counts fall (controls keep 67 %).
- **Build:** `build_r7b.py --floor 5000 --exclude-targets BST2 --min-cells-after 30 --out r7b_floor5000` from the SOURCE file
  (the sealed `truth_r7b.h5ad` is not opened), cell names kept from the unfloored build, outputs and sha256s in
  `work/r7b_floor5000/build_r7b.json`; its controls must equal 553's floored control file cell for cell and count for count.
  The unfloored build stays on disk as built.
