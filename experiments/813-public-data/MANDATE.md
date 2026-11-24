# 813 — Lane B4: public data the project lacks (idea worker)

Written 2026-10-16 ~03:30Z by the lab lead ([c61d4c]). Read `experiments/800-checker/PHASEB.md` first, then
`experiments/270-data-reconnaissance/RESULTS.md`, `experiments/278-depth-inventory/RESULTS.md`, FINDINGS H.2b (datasets and
links), `scripts/pisces_watch.sh` (the production watcher for X-Atlas/Pisces; do not touch it, read its state only).
Territory: finding, vetting and wiring new public perturbation data. Not yours: models (811), DE ideas (812), strong-tier
direction (810).

## Objective
Every held-out instrument here rests on three cancer screens (K562, HCT116, HEK293T) and one Arc line (H1). More independent
lines make the checker stronger (new folds) and the model richer (new donors). Find what exists in 2026: CRISPRi Perturb-seq
(ideally 10x Flex) in cancer lines, ≥ 100 cells per target, many targets; also any release by Arc, Xaira (Pisces), Replogle,
Weissman, Marson, Satija, the IGVF consortium, Parse, scBaseCount-style aggregates.

## Steps
1. Inventory (no download yet): source URL, license, lines, chemistry, targets, cells per target, depth, size, overlap with
   the A/B/C panel and with typical Arc-like panels, and what it would add (a fold, a donor, both). Write it into RESULTS.md.
2. For the top 2–3 finds, PROPOSAL.md (what wiring them buys the checker or a model, cost). The judge decides.
3. A GO find is fetched on EC2 (S3 under output/lab/, nothing under data/, ≤ the lane budget), checked (knockdown, power,
   615's blob test) and wired as a checker fold with the 801 engine's conventions.

## Budget and return
Cloud $50 before a GO. RESULTS.md (inventory), DOSSIER.md only if a dataset becomes a production donor candidate.

## Start note (2026-10-16 06:20Z, lab lead)
Started before the checker's DE known answers are graded, for STEP 1 ONLY (the inventory: web, catalogues, metadata; no
download of any data file, no cloud spend). Steps 2–3 wait for the lab lead's GO after the judge. Write the inventory into
RESULTS.md and stop with a short message when step 1 is done.
