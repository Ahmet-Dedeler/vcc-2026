You are a research scout for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master). Read AGENTS.md
(§1, §5, §7 rules on what data may be used, §14 naming discipline), then docs/FINDINGS.md §0–§2, and the earlier scans
experiments/306-recon-oct/RESULTS.md, experiments/482c-cloud-lit-data/RESULTS.md and experiments/520-cloud-def-lines/RESULTS.md
so you do not repeat them. Experiment 532: write ONLY under experiments/532-cloud-lit-refresh/; small public files only
(< 200 MB total) into work/; commit to your branch only; never run `vcc`; never contact anyone, post anything, or reveal
anything about this repo's methods anywhere.

## Question
What appeared since 2026-09-15 (preprints on bioRxiv/arXiv/medRxiv, journal papers, datasets on GEO/ArrayExpress/Zenodo/
figshare/Hugging Face, code releases, talks or posts by Arc Institute) that bears on our open questions:
1. predicting the SIGN of weak/downstream expression responses to CRISPRi knockdown in an unseen cancer line (generic or
   "common" perturbation-response programmes, lab/chemistry signatures, 10x Flex probe-based technical effects);
2. predicting per-target response MAGNITUDE across cell lines;
3. which cell lines Arc used for the 2026 challenge's final contexts (D/E/F), or new same-lab (Arc) Flex Perturb-seq data;
4. method write-ups by Virtual Cell Challenge 2026 participants (what they did on the DE metrics: fidelity, reach,
   Jaccard, log-fold-change error).
For each item: title, authors/lab, date, link, what it contains, whether it is downloadable (size, licence), and a concrete
statement of how it could change one of our decisions (or "no action"). Verify each by opening its page; mark unverified.
Deliverable: `RESULTS.md` ranked by expected value, plus a download manifest for anything worth fetching locally.
