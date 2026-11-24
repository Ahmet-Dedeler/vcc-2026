# 295 — Claim C5 on reserve R3, the single access: on 174 genes in eight datasets no map, bank or fold ever used, does the 28-bank map's per-gene conservation score predict how well each gene's response transfers from the map's lines to the reserve, and does F13's geometry sign hold? Preregistered in the ledger BEFORE the first expression value is read; one computation; ≤ 4 hours. Analysis only. No build, no submission.

You are a research subagent. Work in `experiments/295-reserve-r3-claim/`. **Edit nothing outside it** except appending
claim C5 and its result to `docs/CONFIRMATION-LEDGER.md` as instructed. No build, no `vcc`, no `cell-eval2`, nothing
under `submissions/` or `data/`; no download. `uv run python`, row-blocked / streaming reads; peak RSS under ~8 GB;
one file at a time. **Time box: 4 hours.** One sibling (293) runs in parallel; do not read or write its directory.
**You are named as the experiment that spends reserve R3** (`docs/CONFIRMATION-LEDGER.md`, "Reserve R3";
`experiments/292-reserve-declaration/RESERVE.md` §3 is the access rule and binds you word for word): metadata reads
are free; `X` of any R3 file may be read ONLY after C5 is written into the ledger with its statistic, pass rule,
baselines, controls and known-answer rows; ONE computation; the result is appended whatever it is; every file you
read is then SPENT and the ledger says so. No pilot, no quick look, no second pass, no parameter chosen on the
reserve. The earlier reserves (HEK293T half-B; the four Broad-16 lines) stay untouched.

## Why (F8, F13 in `docs/FINDINGS-2026-09.md` v1.9; R291; R292 §6)
The write-up's two live claims about cross-line conservation are (F8) that the 28-bank map's per-gene score is a
reproducible property that predicts transfer — confirmed off-map on 291's two surfaces (+0.274 on the GWCD4i authors'
CD4-vs-K562 statistic, +0.210 on Broad-16 at low power) — and (F13) that the SIGN of response geometry is right: broad,
spectrally central responses (`gen_pc20`) are MORE conserved (291: +0.186 and +0.295 on the same surfaces; the fitted
geometry model does NOT beat annotation off-map). Both of 291's surfaces had been inspected before for other claims.
R3 has not: no cell of its eight files is inside any object the programme built. **This is the one confirmation on
untouched data the write-up will carry.** R3's shape (R292 §4, §6) permits ONE pooled claim over ~174 genes with the
file as the bootstrap block, modality (Cas9 KO vs CRISPRi) as a declared stratum, cells per gene as a covariate — not
per-file claims.

## Inputs
R3's manifest `experiments/292-reserve-declaration/work/reserve292.json` (paths, sha256 — verify each digest before
reading), `RESERVE.md` §2 (the files, their `obs` label grammars in `work/parse292.py`, guides per gene, control
labels, cells per gene, the var-axis overlap 2,636–2,739 of the map's 2,739 genes); the map: `experiments/279-znf-state-
dependence/work/{held279.npz, znf279.csv}` (`C_held`, `r_bar`, per-pair `Ccross` and the 28 sources' residual rows via
279's `build279.py` `residual_rows` / `load_unit` on the 2,739-gene axis — the map's CONSENSUS response for a gene is
the mean of its unit residual rows over the banks that carry it, excluding any bank whose line equals the reserve
file's line; state the exclusion per file); 290's feature store `experiments/290-…/work/persrc290.npz` (`gen_pc20`,
`breadth`, `mag` per source × gene — aggregate over the same excluded-line bank subset) and 284's annotation model
predictions `experiments/284-…/work/pred284.csv` (`pred_oof_raw` / the S1 ridge) as the baseline `M_ann`; 291's
`lib291.py` (fit-then-apply, the partial Spearman, the block bootstrap); 285's `lib285.py` (the residual-space
response object from raw cells: `log1p(50,000 · S/ΣS)` minus the side's mean over the file's genes, the file's own
controls; target-gene columns removed).

## The claim (write `PREREGISTRATION.md`, hash it into `work/PREREGISTRATION.sha256`, and append it to the ledger as C5
## — all BEFORE any `X` read; the ledger entry names this experiment, the files it will read, and their digests)
- **Outcome per reserve gene g (computed once):** `A_g` = cosine between the reserve file's residual response of g
  (its perturbed cells vs the file's own controls, on the map's 2,739-gene axis restricted to the file's `var`
  coverage) and the map's consensus residual response of g over banks in OTHER lines. Report also the file's
  split-half reliability of g (random cell halves, 8 draws) so `A_g` can be disattenuated, and cells per gene.
- **C5a (F8 off-map on untouched data):** Spearman(`C_held`, `A_g`) over the pooled genes (≥ 50 cells, in the held
  map; expected ≈ 152 core + 22 Tian), block bootstrap with the FILE as the block (2,000 draws) AND a gene bootstrap;
  **PASS iff ρ ≥ +0.15 and the file-block interval excludes 0.** Baselines: `M_ann` (annotation), magnitude (the map's
  `mag`), the file's own reliability, the gene's DE count in the file (a magnitude proxy computed in the same access);
  a partial Spearman given the file's reliability and magnitude (reported; not in the rule). Null: 200 gene shuffles.
- **C5b (F13's sign on untouched data):** Spearman(`gen_pc20` aggregated over the excluded-line banks, `A_g`) > 0 with
  the file-block interval above 0; `breadth` likewise (reported). **PASS iff `gen_pc20`'s interval is above 0.**
- **Strata, declared:** modality (CRISPRi: Adamson, Dixit, Tian-day7 restricted; Cas9 KO: Datlinger ×2, Papalexi ×2,
  Shifrut), reported separately and pooled; Tian day-7 flagged as same-study-as-a-bank and the pooled result given
  with and without it. **Depth covariate (rule 9):** cells per gene and the file's median UMIs per cell (measured in
  the access), reported next to every cell; the rule is unconditioned, as declared.
- **Known-answer rows (before any reserve number):** every file's sha256 equals the manifest; 279's `c_held` re-read
  (n 17,106, median +0.018605); 291's C4 cell rebuilt from its saved predictions (`work/pred291_gwcd4i.csv`, Stim48hr
  ρ +0.2739 for `C_held` at 0.0e+00); the map consensus for a gene reproduces the mean of `Ccross`'s unit rows where
  the construction overlaps (state the check); the file's controls-vs-controls split reads ≈ 0 (a no-signal row).
- **Multiplicity:** two verdict cells (C5a, C5b). Everything else is secondary and labelled.

## Protocol
1. Metadata + KA rows; write and hash the preregistration; append C5 to the ledger. 2. The single access, all eight
files in one pass each, everything computed in that pass (responses, reliabilities, DE counts, UMIs per cell). 3. The
statistics, the rules, the strata, the baselines, the null. 4. Append the result to the ledger under C5 and mark the
eight files SPENT. 5. Descriptive (post-hoc, labelled): the best- and worst-transferring reserve genes; whether the
Datlinger 2017 → 2021 same-library pair and Papalexi pooled → arrayed pair agree (the construct / assignment readings
R292 §2.3 named), reported as descriptive only.

## Deliverable
`PREREGISTRATION.md` (hashed), the ledger entries (claim before; result after; files marked SPENT), `RESULTS.md`:
KA rows; the per-file table (n genes, cells/gene, UMIs/cell, reliability, median `A_g`); C5a and C5b with both
bootstraps, the strata, the baselines, the null; the descriptive pairs; `work/reserve295.csv` (per gene: file, cells,
reliability, `A_g`, `C_held`, `gen_pc20`, `breadth`, `mag`, `M_ann`); one line: **C5a (F8 on untouched data): PASS /
FAIL (ρ … [file-block lo, hi]; annotation …; magnitude …); C5b (F13's sign): PASS / FAIL (`gen_pc20` ρ … [lo, hi];
`breadth` …); CRISPRi-only … / Cas9-only …**. No model use for the competition is proposed. Then stop.
