# 303 — Claim C7 on reserve R5, the single access: with a SECOND overexpression library measured at bulk depth (Nakatake et al. 2020, 510 inducible TF-overexpression hESC lines), is a transcription factor's overexpression response the reverse of its CRISPRi knockdown response in the same line (hESC, the map's `gse295214_h1` bank) — and does knockdown geometry predict the SIZE of the overlap (302's post-hoc finding, now preregistered on independent data)?

You are a research subagent. Work in `experiments/303-second-overexpression-library/`. **Edit nothing outside it**
except appending claim C7 (before) and its result (after) to `docs/CONFIRMATION-LEDGER.md` under the "Reserve R5"
heading the coordinator wrote at 04:20Z. No build, no `vcc`, no `cell-eval2`, nothing under `submissions/` or
`data/`; no download. `uv run python`; row-blocked h5py reads of
`data/external/gse295214/pseudo_bulk_perturbation_database.h5ad` (read only the rows you need; `X` is dense float32,
20,345 × 18,108); peak RSS under ~6 GB. **Time box: 3 hours.** Two siblings (304, 305) run in parallel; do not read or
write their directories. **Reserve rule (R3's, word for word):** the `Nakatake2020_hESC_TF_screen` rows of `X` (and
their `layers/control` rows) may be read ONLY after C7 is written into the ledger with its statistic, pass rules,
baselines, controls and known-answer rows; ONE computation; the result is appended whatever it is; the rows are then
SPENT. `obs`, `var`, `layers/measure_mask` and every other dataset's rows are free. R1–R4 are spent: the JoungZhang
atlas's `X` is NOT to be opened; 302's saved sums (`experiments/302-orf-vs-crispri-symmetry/work/sums302.npz`, keys
`axis`, `acc_sym`, `acc_code`, `cells_sym`, `lib_sym`, `used_codes` — the pooled per-TF count sums from the spent
access, reconstructed exactly as `lib302.py` does) may be used for DESCRIPTIVE readings only, never a verdict cell.
Nothing competition-facing is selected or tuned.

## Why
302 (`experiments/302-orf-vs-crispri-symmetry/RESULTS.md`, F15 in `docs/FINDINGS-2026-09.md` v2.7) asked whether
overexpression reverses knockdown in hESC and got VOID verdicts (its KA-E clause) and a bounded null: the atlas's
ORF side has split-half reliability +0.025 at a median 110 cells per TF, which caps any cosine at 0.11, so a true
mirror of order −0.05 was invisible. 302 named the data that would decide it: a second overexpression library. It is
on disk and unread: Nakatake et al. 2020 (Cell Reports; 2,135 hESC lines with PiggyBac inducible overexpression of
one TF each, profiled by RNA-seq after induction) sits in the pseudobulk database as 510 one-row-per-TF profiles
(`obs.dataset == 'Nakatake2020_hESC_TF_screen'`, `obs.id` `hESC_TF_screen_<SYMBOL>`, `perturbation_sign` +1,
`measure_mask` coverage 0.872). Bulk depth removes the cell-count ceiling; a different library removes the
isoform/construct confound. Metadata (04:05Z): Nakatake ∩ Pan2025 knockdown = 100 TFs; ∩ 302's cohort = 58 (the
cross-library reading); the database's own 28 Joung rows ∩ Nakatake = 13.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md` (rules 1, 2, 5, 9, 10), `docs/CONFIRMATION-LEDGER.md` (C6's entry and the R5
declaration; C7 goes under R5), R302 (§ the design, the KA rows, "what would decide it"), 302's `lib302.py`
(residual convention, cosine, own-column removal, the TF bootstrap, the shuffled null — REUSE its functions),
115's `acquire.py` (how the `Pan2025_hESC` rows became the map's `gse295214_h1` bank at
`experiments/115-gse295214-h1/work/compact-donor/`: `delta.npy` rows vs `b_ctrl.npy`, `genes.txt`, `perts.txt`,
`measured_mask.npy` — the SAME transform must be applied to the Nakatake rows so both sides sit in one convention;
state it), 279's `build279.py` (`load_unit`), 290's feature store (`persrc290.npz`: `gen_pc20`, `breadth`, `mag` per
source × gene; use the `gse295214_h1` source's own columns), 279's `work/znf279.csv` (`C_held`).

## Metadata pass first (free; no Nakatake `X` row)
The cohort (Nakatake ∩ the bank's `perts.txt`, expected 100), the shared gene axis (database `var.gene_name` ∩ the
bank's `genes.txt` ∩ Nakatake `measure_mask` ∩ the map's held genes; state the size), the 58-TF cross-library
cohort, the database's provenance line (`uns` is empty; use `obs.id` patterns and the file's own field names).

## The claim (write `PREREGISTRATION.md`, hash into `work/PREREGISTRATION.sha256`, append to the ledger as C7 — BEFORE any Nakatake `X` row)
- **Unit:** the 100 TFs (state the exact count after axis filters); one Nakatake row and one `gse295214_h1` residual
  row each, on the shared axis, the TF's own column removed.
- **Objects (one computation):** overexpression response = the Nakatake row in the bank's convention (the same
  transform 115 applied to Pan2025 rows against `layers/control`), residualised by the mean over all 510 Nakatake
  rows on the axis (the map's convention); knockdown response = the bank's residual row; **S = cos(OE residual, KD
  residual)** per TF; magnitudes of both sides; the Nakatake row's `measure_mask` coverage as a covariate.
- **C7a (antisymmetry exists):** median S < 0 with the paired TF-bootstrap 95 % interval (20,000 draws) below 0, AND
  the share of TFs with S < −0.05 exceeds the share with S > +0.05. Null: TF-shuffled pairings (200 draws). **PASS
  iff both hold.**
- **C7b (conservation predicts mirroring):** Spearman(`C_held`, −S) > 0 with the interval above 0. Baselines: KD
  magnitude, OE magnitude, `r_bar`, coverage; partial given both magnitudes (reported, not in the rule).
- **C7c (F13's sign, second modality):** Spearman(`gen_pc20` of the KD row, −S) > 0 with the interval above 0.
- **C7d (302's size finding, preregistered):** Spearman(`gen_pc20` of the KD row, |S|) > 0 with the interval above 0;
  partial given both magnitudes and coverage reported. **PASS iff the interval is above 0.**
- **Declared secondaries (reported whatever they show, no verdict):** (i) the cross-library reading on the 58 TFs:
  cos(Joung ORF residual rebuilt from `sums302.npz`, Nakatake residual) — the overexpression-side study term in
  hESC, F14's structure in a second modality, and the disattenuation it implies for 302's null (state the
  formula); (ii) the same S restricted to TFs whose Nakatake own-transcript is up (the transgene detected); (iii)
  S against 302's `S_g` on the 58 (agreement of two overexpression libraries about the same knockdown row); (iv)
  the top 15 TFs each way, post-hoc, labelled; (v) the tail shares at ±0.10.
- **Known-answer rows.** LOADER rows (abort clauses: a failure VOIDS the verdicts, as 302's KA-E did): the Pan2025 rows
  through YOUR loader reproduce the bank's `delta.npy` rows for 20 random targets at ≤ 1e-6 (or state the transform
  and reproduce to that tolerance); 279's `Ccross` entries involving `gse295214_h1` reproduce at 0.0e+00 through the
  same path; the cohort counts 100 / 58 / 13 reproduce the coordinator's metadata check. DATA rows (reported, and
  they CONDITION THE INTERPRETATION, never the verdict — the lesson of 302): the overexpressed TF's own transcript is
  up in the Nakatake row (share of TFs; the transgene may or may not map to the endogenous gene — say which); the
  Nakatake row's coverage; a no-signal row (a Nakatake row against a random OTHER TF's knockdown row is the shuffled
  null's element — report its distribution next to the real one).
- **Multiplicity:** four verdict cells (C7a–C7d); everything else secondary and labelled. Depth (rule 9): bulk lines
  — cells not applicable; state what the file records (coverage) and that the OE side is not cell-count limited.

## Protocol
1. Metadata pass. 2. Write and hash the preregistration; append C7 under R5 in the ledger. 3. The single access:
one row-blocked read of the 510 Nakatake rows and their `control` rows; everything downstream from those arrays,
saved to `work/nak303.npz`. 4. Statistics, rules, baselines, null, secondaries. 5. Append the result to the ledger;
mark R5 SPENT. 6. Post-hoc tails, labelled.

## Deliverable
`PREREGISTRATION.md` (hashed), the ledger entries, `RESULTS.md` (KA rows; the cohort; the S distribution with the
null; C7a–C7d with intervals and baselines; the secondaries, the cross-library term first), `work/symmetry303.csv`
(per TF: coverage, own-transcript delta, S, |S|, `C_held`, KD mag, OE mag, `gen_pc20`, `breadth`, 302's `S_g`
where present); one line: **C7a (a bulk-depth overexpression library reverses knockdown, hESC): PASS / FAIL
(median S … [lo, hi]; shares …; null …); C7b: PASS / FAIL (ρ …); C7c: PASS / FAIL (ρ …); C7d (geometry predicts
size): PASS / FAIL (ρ …); cross-library OE term on 58 TFs: … [lo, hi]**. No model use for the competition is
proposed. Then stop.
