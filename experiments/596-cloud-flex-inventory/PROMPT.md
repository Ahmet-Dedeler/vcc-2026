# 596 — Inventory of public 10x Flex (Fixed RNA Profiling) datasets of named human cancer cell lines (cloud agent)

You work only in `experiments/596-cloud-flex-inventory/` on your branch. Do not edit anything else. **Do not download data
files**: read web pages, API listings, metadata, sample tables and file listings only (a file's size from its listing is
fine; its bytes are not). Do not run `vcc`. Nothing about this repository's methods leaves it.

## Why
Experiment 590 (`experiments/590-identity-validation/RESULTS.md`) validated the day-of cell-line identification (a CCLE
matcher on control-cell pseudobulks) on 402 public profiles. Almost all of them are 10x 3′, 5′ or sci-RNA-seq chemistry.
The final test lines will be measured with **10x Flex** (probe-based Fixed RNA Profiling), and the only Flex profiles we
hold are the challenge's own three validation lines and one H1 stem-cell screen. 590's match-quality floor (0.65; correct
Flex calls 0.719–0.772, false calls of lines absent from CCLE 0.51–0.61) was transferred to Flex through a simulated
chemistry bias. Real Flex profiles of other named cancer lines would test it directly.

## What to find
Every public dataset where named human cancer cell lines (or other immortalised human lines, e.g. HEK293, RPE1) were
profiled with **10x Flex / Fixed RNA Profiling / Chromium Flex** (any version: v1 human probe set, Flex v2, Flex Apex, Flex
multiplexed with probe barcodes), or a close probe-based fixed-RNA chemistry (say which). Search at least:
- 10x Genomics' own public datasets (support.10xgenomics.com/datasets and 10xgenomics.com/datasets; filter Flex / Fixed RNA
  Profiling; many demo sets are cell-line mixtures or cell lines);
- GEO (search terms: "Fixed RNA Profiling", "Chromium Flex", "10x Flex", "FRP", "probe-based", with "cell line");
- CELLxGENE Census / Discover, the Human Cell Atlas data portal, Zenodo, figshare, Single Cell Portal, ArrayExpress /
  BioStudies, the Arc Institute Virtual Cell Atlas listings, X-Atlas (Xaira), Parse is NOT Flex (exclude, but you may list
  probe-based fixed alternatives separately);
- Perturb-seq on Flex in cancer lines (any lab), since their non-targeting controls are exactly the needed profile.

## Per dataset, report
URL (landing page and file listing) · accession · named line(s) (as written, and the CCLE/DepMap name if obvious) · chemistry
and version as stated by the source · cells per line (or total and how lines are separated: probe barcodes, genotype
demux, separate libraries) · controls (unperturbed or non-targeting cells; how many) · file formats and sizes · licence /
terms of use · release date · anything that would make it unusable for a pseudobulk identity test (e.g. drug-treated
only, mixed lines without labels, a targeted gene panel instead of the whole-transcriptome probe set).

## Deliverable (branch `cloud/596-flex-inventory`, only in `experiments/596-cloud-flex-inventory/`)
`RESULTS.md`:
1. A verdict line (how many usable datasets, how many distinct lines, any line likely to be among the final test lines);
2. the table above, sorted by usefulness;
3. for every usable dataset with a permissive licence, a proposal (do NOT run it) for re-checking 590's 0.65 floor on
   it: which files, how to build each line's control pseudobulk on the challenge's gene axis
   (`data/challenge/gene_names.csv`), and how to score it with `scripts/def_identity.py`'s matcher the way
   `experiments/590-identity-validation/match590.py` does, including a "line absent from CCLE" case if one exists
   (e.g. HEK293, primary or engineered lines);
4. what you searched and found empty (queries, dates), so the next sweep does not repeat it.
Mark each fact VERIFIED (you opened the page / API record) or UNVERIFIED. Commit those files to your branch and stop. Do
not open a pull request.
