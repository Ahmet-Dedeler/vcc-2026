# 424: VIPerturb-seq K562 as a donor table (mandate from the PI, 2026-09-23 01:35Z)

**Why.** 413's source scan (`experiments/413-source-scan/RESULTS.md` §2, §5, §6) found one new public item worth building.
VIPerturb-seq (Bradu … Satija 2026, Zenodo 18460279, CC BY 4.0) is a genome-wide K562 CRISPRi screen with a fixed-cell
probe readout (Flex-style, the challenge's chemistry): 18,885 targets, **300/300 of our A/B/C panel in the manifest,
284/300 at ≥ 20 cells** (median 47), and 7,949 `NO-TARGET` control cells.
- For A/B/C it adds a K562 row for the **30 panel targets our K562-GWPS donor lacks**. Five of those (ABCD1, EPHB2,
  KIF21B, NICN1, SPATA20) are among our thinnest-covered.
- For the other ~254 targets it is a second, independent K562 measurement.
- For October's D/E/F panel, it is the only way to get a K562 row for targets outside K562-GWPS's ~9.8k.

Read `docs/FINDINGS.md` §0–§2 and 413's RESULTS first. The data layer is `data/donors/<name>/` (look at
`data/donors/k562_gwps/`: `delta.npy`, `b_ctrl.npy`, `genes.txt`, `perts.txt`, plus whatever manifest the builder reads)
and `src/vcc2026/`. Find how K562-GWPS's donor table was built (grep `experiments/` and `scripts/` for the builder of
`data/donors/k562_gwps`) and **mirror that construction exactly** (pseudobulk, normalisation, gene axis, delta
definition).

**Territory.** Nobody else touches donor data tonight. 417 builds bulk arms (the production builder is pinned most of
the night: do not touch it), and 423 works the DE channel on H1. Do not edit the builder, `data/donors/` of existing
donors, vcc or cell-eval2.
- (a) **Acquire.** Download the genome-wide bins (binA/B/C RDS, ~10.2 GB, plus the manifest), or only what the donor
  table needs.
  - Seurat v5 RDS: R 4.6.1 is installed locally (`Rscript`). Install SeuratObject/Matrix into a user library if
    needed, or try `rds2py`.
  - Check memory before loading a bin (`memory_pressure` ≥ 40 % free and free+inactive ≥ 20 GiB; ~15 GB builds share
    this Mac). If a bin needs more than ~20 GB, do the conversion on AWS instead:
    - the CLI is configured, and `scripts/cloud_score.sh` shows how the repo launches and tears down an EC2 worker;
    - budget ≤ $10;
    - run `scripts/aws_reap.sh` at the end and state in RESULTS.md that no instance is left running.
  - Keep only what the table needs on disk. Delete raw RDS files after conversion unless under 2 GB total.
- (b) **Build** `data/donors/viperturb_k562/` in the exact K562-GWPS format: per-target pseudobulk deltas vs the
  `NO-TARGET` controls, same gene axis rules, ≥ 20 cells per target. Write a sidecar JSON with provenance: URLs,
  file shas, cell counts per target, the guide/knockdown QC you applied and why.
- (c) **Known answers and validation:**
  - effect correlation with K562-GWPS on shared targets (per-target cosine distribution; it should be clearly positive
    and comparable to 413's or older cross-screen numbers if any exist);
  - the knockdown check: each target's own gene goes down;
  - the panel coverage count must reproduce 413's 284.
- (d) **Price it on PLOLO** with 414/418/420's instruments (`experiments/420-bulk-stack/` has the latest machinery).
  Report the relative aggregate cosine and board-emulated pds on the three cancer folds + kolf2, and predicted
  Δoverall at k = 1 and 1.407, for:
  - (i) VIPerturb as an extra donor at equal weight in the base;
  - (ii) VIPerturb averaged with K562-GWPS into one K562 row (replicate denoising);
  - (iii) VIPerturb only where K562-GWPS is missing (the 30-target fill).

  Hold the operator fixed; say what a refit would change. Folds must exclude the held line's own data, as the
  existing folds do.

**Deliverable:** `experiments/424-viperturb-k562/RESULTS.md`:
- acquisition record (with the AWS cost if used);
- the donor table's validation;
- the PLOLO table;
- a recommendation: which of (i)/(ii)/(iii), if any, is worth a builder arm, and what the builder change would be.
  Do not make the builder change.

**Constraints.**
- ≤ 8 GB RSS for local Python; the memory gate before any heavy step.
- Write only in experiments/424-viperturb-k562/ and the new `data/donors/viperturb_k562/`.
- No git.
