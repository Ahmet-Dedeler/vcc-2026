# 302 — Claim C6 on reserve R4, the single access: in the SAME cell line (hESC), is a transcription factor's overexpression response the reverse of its knockdown response — and does the map's cross-line conservation score, or the knockdown response's geometry, predict how strongly it mirrors? JoungZhang 2023 TF-ORF atlas (1.15 M cells; `X` never read until the claim is in the ledger) × the map's hESC CRISPRi screen (gse295214_h1); 249 TFs in both. Preregistered; one computation; ≤ 4 hours. Analysis only. No build, no submission.

You are a research subagent. Work in `experiments/302-orf-vs-crispri-symmetry/`. **Edit nothing outside it** except
appending claim C6 (before) and its result (after) to `docs/CONFIRMATION-LEDGER.md`, and converting the "Reserve R4
candidate" note there into a declared reserve spent by this experiment. No build, no `vcc`, no `cell-eval2`, nothing
under `submissions/` or `data/`; no download. `uv run python`, row-blocked / streaming reads of the 5.8 GB atlas
(h5py; read only the rows you need); peak RSS under ~8 GB. **Time box: 4 hours.** **Reserve rule (R3's, word for
word):** `data/external/scperturb/JoungZhang2023_atlas.h5ad`'s `X` may be read ONLY after C6 is written into the
ledger with its statistic, pass rules, baselines, controls and known-answer rows; ONE computation; the result is
appended whatever it is; the file is then SPENT. Metadata (`obs`, `var`, `uns`) is free. R1–R3 are spent and are
not touched. Every panel target's identity is irrelevant here; nothing competition-facing is selected or tuned.

## Why
The programme's organising question is which perturbation responses are conserved and which depend on the recipient's
state. Every measurement so far compares the SAME kind of perturbation (knockdown) across contexts. R4 allows a
question no dataset on disk could ask before: **within one context, does the direction of a gene's response reverse
when the perturbation reverses?** Joung et al. 2023 overexpressed 3,369 TF ORFs in hESCs (SHARE-seq; 1,145,823 cells;
GFP / mCherry controls, 49,539 + 78,274 cells — 301 resolved the ORF ids, `experiments/301-metadata-acquisitions/
work/joung_orf_map.csv`, 3,136 / 3,369 codes, 99.5 % of cells); the map's `gse295214_h1` bank (Pan et al. 2025, CRISPRi
in hESCs, 2,354 perturbations, `experiments/115-gse295214-h1/work/compact-donor/`: `delta.npy` rows vs `b_ctrl.npy`,
`genes.txt`, `perts.txt`, `measured_mask.npy`) knocked down 249 of the same TFs (233 in the 17,106-gene held map;
148 with ≥ 100 ORF cells — measured from metadata only, 2026-09-07 03:25Z). If overexpression mirrors knockdown, the
mirrored TFs are the ones whose response is a dosage-monotone programme; the programme predicts those are the
conserved, core-machinery-like responses (F9, F13) and that the recipient-specific ones are not mirrored.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md` (rules 1, 2, 5, 9, 10), `docs/CONFIRMATION-LEDGER.md` (C5's format; the R4
candidate note), R301 §3 (how the ORF codes were resolved; the controls' codes 3367 = GFP, 3368 = mCherry; isoforms per
TF), 279's `build279.py` (`load_unit`: a bank's residual rows on the map's axis), 276's `DEFINITION-276.md` (residual =
minus the bank's mean over its perturbations; cosine on the common axis; the target's own column removed), 290's
feature store (`experiments/290-remainder-response-features/work/persrc290.npz`: `gen_pc20`, `breadth`, `mag` per
source × gene — use the `gse295214_h1` source's own columns for the knockdown-geometry features, and the map-wide
aggregate as a second reading), 279's `work/znf279.csv` (`C_held`, `r_bar`), 285's `lib285.py` (split-half reliability
from cells).

## The claim (write `PREREGISTRATION.md`, hash into `work/PREREGISTRATION.sha256`, append to the ledger as C6 — BEFORE any `X` read)
- **Unit:** TFs with ≥ 50 assigned ORF cells (pooled across isoforms of the same symbol; state the count) AND a
  `gse295214_h1` row; expected ≈ 200; the primary cohort ≥ 100 cells (≈ 148).
- **Objects (one computation):** ORF response per TF = `log1p(5e4·S/ΣS)` of its pooled cells minus the same statistic
  of the GFP + mCherry control cells, residualised by the mean over all scored TFs' ORF responses (the map's
  convention), on the genes shared by the atlas's `var` and the H1 bank (state the axis size), the TF's own column
  removed; knockdown response = the H1 bank's residual row for the TF on the same axis (279's `load_unit` path);
  **S_g = cos(ORF residual, CRISPRi residual)** per TF; the ORF side's split-half reliability (random cell halves, 8
  draws) for attenuation; cells per TF and the ORF isoform count as covariates.
- **C6a (antisymmetry exists):** median S_g < 0 with the paired TF-bootstrap 95 % interval below 0, AND the share of
  TFs with S_g < −0.05 exceeds the share with S_g > +0.05. Null: TF-shuffled pairings (200 draws; expected median ≈ 0).
  **PASS iff both hold on the ≥ 100-cell cohort.**
- **C6b (conservation predicts mirroring):** Spearman(`C_held`, −S_g) > 0 with the interval above 0 on the same
  cohort (a conserved knockdown response is more strongly reversed by overexpression). Baselines: knockdown magnitude
  (`mag` of the H1 row), ORF magnitude, the ORF side's reliability, `r_bar`; partial Spearman given both magnitudes
  and the ORF reliability (reported, not in the rule). **PASS iff the interval is above 0.**
- **C6c (F13's sign, new modality):** Spearman(`gen_pc20` of the H1 knockdown row, −S_g) > 0 with the interval above 0
  (broad, spectrally central knockdown responses are the mirrored ones). **PASS iff the interval is above 0.**
- **Declared secondaries:** the same with S_g's absolute value (symmetry of either sign); the 50–99-cell cohort;
  per-isoform S_g where a TF has ≥ 2 isoforms with ≥ 50 cells each (isoform agreement as a reliability reading);
  the fraction of TFs whose ORF response is instead POSITIVELY aligned with their knockdown response (a
  dosage-insensitive / stress-like response — describe the top 20 each way, post-hoc, labelled).
- **Known-answer rows before any reserve number:** the atlas's control-code cell counts equal 301's (49,539 /
  78,274); the H1 bank's rows reproduce 279's `Ccross` entries for pairs involving `gse295214_h1` at 0.0e+00 through
  the same loader; the ORF map's TF list equals 301's `joung_orf_map.csv` (sha256 recorded); a controls-vs-controls
  split (GFP vs mCherry cells, same statistic) reads |cos| < 0.05 against every knockdown row (a no-signal row);
  the 249-TF overlap and the 148-TF ≥ 100-cell count reproduce the coordinator's metadata check.
- **Multiplicity:** three verdict cells (C6a, C6b, C6c); everything else secondary and labelled. Depth (rule 9):
  cells per TF and ORF UMIs per cell reported next to every cell; the rules are unconditioned, as declared.

## Protocol
1. Metadata pass (free): cohort counts, control counts, isoform counts, the shared gene axis, KA rows that need no `X`.
2. Write and hash the preregistration; append C6 to the ledger; convert the R4 note to a declared reserve named for
   this experiment. 3. The single access: one streaming pass over the atlas's `X` collecting the pooled sums for every
   cohort TF (and its two random cell halves) and the two control groups; everything downstream from those sums.
   4. The statistics, the rules, the baselines, the null. 5. Append the result to the ledger; mark R4 SPENT.
   6. Descriptive tails (post-hoc, labelled).

## Deliverable
`PREREGISTRATION.md` (hashed), the ledger entries, `RESULTS.md`: KA rows; the cohort; the S_g distribution with the
null; C6a–C6c with intervals and baselines; the secondaries; `work/symmetry302.csv` (per TF: cells, isoforms, ORF
reliability, S_g, `C_held`, `mag`, `gen_pc20`, `breadth`); one line: **C6a (overexpression reverses knockdown, hESC):
PASS / FAIL (median S … [lo, hi]; share < −0.05 … vs > +0.05 …; null …); C6b (conservation predicts mirroring): PASS /
FAIL (ρ … [lo, hi]; partial …); C6c (geometry sign): PASS / FAIL (ρ … [lo, hi])**. No model use for the competition is
proposed. Then stop.
