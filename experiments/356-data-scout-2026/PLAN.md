# 356-data-scout-2026 — PLAN

## Hypothesis
The bottleneck is cross-line transfer information (343: design-stage |lfc|-weighted sign
accuracy 0.600 on the reference's significant genes; 325–342 all negative). A dataset
released/updated in 2026 — genome-wide CRISPRi/a Perturb-seq in a new line, Perturb-seq in
Jurkat/HeLa/CAL33, a DepMap 2026 release, or an open-weights perturbation foundation model
with zero-shot line transfer — could supply independent source-stage labels the six-line
frame (334/336/339/340/341) lacks. Value is judged only as prospective source information,
not as a score claim.

## Inputs (read-only; no downloads > 2 GB, none without size+license note)
- `data/external/` inventory (35 entries listed 2026-09-10; already-held list in mandate §"We already have").
- `docs/`: grep `external/` for ingestion notes; `docs/READ-FIRST.md`, `docs/RESEARCH-PROGRAMME-2026-09.md`,
  `docs/LOG.md`, `docs/BOARD-LEDGER.md` for context.
- `data/leaderboard/*.json`, `data/board-entries/*.json` — only if a competitor trace names a dataset.
- Web: bioRxiv/medRxiv search, GEO DataSets (2026 series), HuggingFace datasets, Arc Virtual Cell
  Atlas, Tahoe/Vevo, Xaira X-Atlas/Pisces, Parse GIGA atlas, scPerturb, IGVF portal, PerturbBase,
  CZI cellxgene, DepMap portal (26Q1/26Q2), foundation-model repos (scGPT/Geneformer/STATE successors).
- Context identities (041): A = Jurkat, B = HeLa, C = CAL33 — target line-specific queries.

## Estimator / method
No computation on counts. For each candidate record: URL, access date, size, license, format,
n_cells, n_perturbations, cell lines, modality. Verify against primary source (GEO accession page,
paper, HF dataset card, portal). Cross-check against already-held list to avoid re-listing.
Rank ≤ 5 by expected transfer value for A/B/C and October (D/E/F unknown lines), with licensing flags.

## Controls / falsifiers
- Known-answer reproduction: re-read 323's `work/metrics-cap15.json` mean NMAE and 343's headline
  numbers with own code before any new claim (rule 3).
- Negative control: query each "new" candidate against the already-held list + docs ingestion notes;
  anything already held (or a re-release of the same cells) is excluded from the ranked list.
- Kill criterion: if no candidate offers (a) a new cell line with genome-wide CRISPR perturbation +
  transcriptome readout, (b) same-line (Jurkat/HeLa/CAL33) perturbation transcriptomes, or
  (c) open-weights transfer model/embedding with a 2026 release, the verdict is NEGATIVE
  (nothing to acquire) with the search log as evidence.
- What kills the idea: candidates all commercial/restricted license, bulk-only with no usable
  signature count, < 2 GB download of a re-packaged held dataset, or model weights closed.

## What result would change the verdict
A license-clean, downloadable source with ≥ 1 new informative line (or same-line perturbation labels)
at usable scale, with format/size verified on the primary page, earns a ranked acquisition slot.
Otherwise NEGATIVE.

## Budget / limits
3 h, web research only (+ reading local inventory). No downloads > 2 GB; any download < 2 GB only
with size+license in this plan first (none planned a priori). No reserve paths, no paid APIs, no builders.
Peak RSS target < 1 GiB (JSON reads only).
