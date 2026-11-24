# 532 sources opened this session

Dates are the dates on the records, not the date of this scout. Scout day: 2026-09-29.

## APIs

- bioRxiv `https://api.biorxiv.org/details/biorxiv/2026-09-15/2026-09-29/{cursor}` — `total` 3443. medRxiv same interval — `total` 1493. Titles and abstracts filtered locally. Full text opened only for the rows in RESULTS.md.
- GEO esearch `gds`, PDAT 2026/09/15–2026/09/30, Perturb-seq OR CROP-seq OR "CRISPR interference" OR CRISPRi. Nine series, esummary titles listed in RESULTS.md.
- Europe PMC search (FIRST_PDATE 2026-09-01 to 2026-09-30, Perturb-seq / CRISPRi / virtual cell) returned 23 hits; the in-window subset is the same papers as the bioRxiv filter.
- Hugging Face dataset API, `Xaira-Therapeutics/X-Atlas-Pisces`: `usedStorage` 0, `lastModified` 2026-03-25T22:58:18.000Z, three siblings. Saved as `hf_pisces.json`.
- Hugging Face `datasets?author=arcinstitute&sort=lastModified`. Newest rows: `opengenome2` 2026-09-02, then July 2026 and older. No new Perturb-seq set.
- Hugging Face `datasets?author=Xaira-Therapeutics`. Pisces 2026-03-25, Orion 2025-09-19.
- Leaderboard `https://virtualcellchallenge.org/api/leaderboard` — 1212 entries. Descriptions quoted in RESULTS.md are from that response.
- GitHub API repository objects: SPECTRA pushed 2026-09-23, MIT, size 102516 KB; SHERLOCK pushed 2026-07-26, size 336683 KB, licence null; IGNITE_public pushed 2026-09-14, MIT, size 10366 KB; cell-eval2 pushed 2026-08-22; state pushed 2026-07-24; X-Cell pushed 2026-03-16. `GabbyKoki/pop` returned 404.
- GitHub search repositories `virtual cell challenge 2026`: 15 hits, newest `MMatinGerami/zero-shot-perturbation` pushed 2026-09-29.
- BioStudies `Perturb-seq AND release_date:[2026-09-01 TO 2026-09-30]`: totalHits 0.
- Figshare articles search `CRISPRi Perturb-seq`, published_date desc: newest 2026-08-26.
- Zenodo records `"Perturb-seq"`, mostrecent: in-window hits are not a human Perturb-seq matrix (see RESULTS).
- arXiv API `id_list=2604.16642,2609.34391,2609.29221`. Directional coherence v4 was updated 2026-08-27, outside this window, so it is not an item. P2P v1 published 2026-09-28. PRAXIS v1 published 2026-09-24.
- bioRxiv details `10.64898/2026.03.18.712807`: X-Cell still version 1, date 2026-03-20.

## Pages and PDFs whose text was read

- https://arcinstitute.org/news — index text. Latest challenge post linked `/news/virtual-cell-challenge-2026`, dated 20 August 2026 in the index.
- https://www.cell.com/cell/fulltext/S0092-8674(26)00931-1 — commentary HTML. No line names in the fetched text.
- SPECTRA HTML full text, doi 10.64898/2026.09.21.752624v1. Code-availability sentence names `www.github.com/MCalabroCode/SPECTRA` and `www.virtualcellchallenge.org`.
- IGNITE HTML full text, doi 10.64898/2026.09.12.748195v1. Data-availability sentence names `https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq`. That CZI URL returned a script shell only; file list unverified.
- SHERLOCK, Pop-Corn, DELPHAI v2, PerturbBridge PDFs from biorxiv.org. P2P and PRAXIS PDFs from arxiv.org. Availability sentences are quoted in RESULTS.md.
- Matin README `https://raw.githubusercontent.com/MMatinGerami/zero-shot-perturbation/main/README.md` — table row “weighted transfer, 6 screens (26 Sep)”, overall 0.079.

## Not opened

Discord, email, and any contact. `vcc` was not run. Full text was not opened for the one-line abstract-only items in RESULTS.md; those abstracts came from the bioRxiv details API.
