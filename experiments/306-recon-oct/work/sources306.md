# 306 — sources

Every source consulted, with the exact sentence relied on. All accessed **2026-09-07, 15:44–17:05 UTC**
from `~/Code/virtual-cell-challenge-2026/experiments/306-recon-oct/`.
Public, unauthenticated access only. Nothing downloaded into `data/`; total bytes fetched into
`work/` ≈ 1.2 MB (JSON/HTML/text only — see `ls -la` at the end of this file's directory).

Convention: **OBSERVED** = I ran the request and read the bytes; **REPORTED** = a page or a
third party asserts it and I did not verify it independently.

---

## A1 — X-Atlas/Pisces and X-Cell

**S1. HuggingFace dataset API, `Xaira-Therapeutics/X-Atlas-Pisces`** — OBSERVED
`https://huggingface.co/api/datasets/Xaira-Therapeutics/X-Atlas-Pisces?full=true&files_metadata=true`
(raw JSON saved as `work/hf_pisces.json`, HTTP 200, 1,158 bytes).
Fields read verbatim:
```
id = Xaira-Therapeutics/X-Atlas-Pisces
lastModified = 2026-03-25T22:58:18.000Z
createdAt    = 2026-03-25T21:53:52.000Z
private = False   gated = False   disabled = False
downloads = 193   likes = 14
usedStorage = 0
siblings (3): .gitattributes | LICENSE.md | README.md
cardData.license = cc-by-nc-sa-4.0
```
There is no data file of any kind in the repo and `usedStorage` is zero.

**S2. Pisces README** — OBSERVED, `work/pisces_README.md` (956 bytes),
`https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces/raw/main/README.md`.
Exact sentences:
> "X-Atlas/Pisces is a Perturb-seq atlas containing 25.6 million perturbed single-cell
> transcriptomes across 16 biologically diverse contexts, including widely used cell lines
> (HCT116, HEK293T, HepG2), induced pluripotent stem cells (iPSCs), resting and CD3/CD28
> activated Jurkat T lymphoma cells, and multi-lineage differentiating iPSCs."

> "**(Coming Soon)** The following data will be uploaded to this dataset:
> 1. All cells from X-Atlas/Orion (cells with a valid dual-guide pair from the **HCT116** and
>    **HEK293T** screens)
> 2. Cells containing a valid dual-guide pair AND perturbations in the 200-gene validation set
>    OR a non-targeting control from the **HepG2** and **iPSC** screens
> 3. All cells with a valid dual-guide pair from the **Jurkat Resting** and **Jurkat Active** screens"

**S3. HuggingFace model API, `Xaira-Therapeutics/X-Cell`** — OBSERVED, `work/hf_xcell.json`.
```
lastModified = 2026-03-16T23:38:37.000Z   createdAt = 2026-03-14T00:06:40.000Z
usedStorage = 272925   downloads = 0   likes = 15
siblings (3): .gitattributes | README.md | x-cell-overview.png
license = cc-by-nc-sa-4.0
```
No weight file. `downloads = 0`.

**S4. All HuggingFace repos owned by `Xaira-Therapeutics`** — OBSERVED
(`/api/datasets?author=…`, `/api/models?author=…`). Exactly three repos exist:
`X-Atlas-Orion` (lastModified 2025-09-19, downloads 36,279), `X-Atlas-Pisces` (2026-03-25,
downloads 193), `X-Cell` (2026-03-16, downloads 0). No new repo, no mirror.

**S5. Xaira GitHub org** — OBSERVED, `https://api.github.com/orgs/xaira-therapeutics/repos`.
One public repo: `Xaira-Therapeutics/X-Cell`, `pushed_at = 2026-03-16T23:57:26Z`, 111 stars.
Nothing pushed since March.

**S6. bioRxiv API, X-Cell preprint** — OBSERVED,
`https://api.biorxiv.org/details/biorxiv/10.64898/2026.03.18.712807` →
`version 1 | date 2026-03-20 | published: NA`. One version only, never revised, never journal-published.

**S7. X-Cell preprint text already on disk** (`data/external/xcell/xcell_biorxiv_full.txt`, read-only) —
OBSERVED. Its Data Availability section reads
> "Data Availability The following data will be uploaded to https://huggingface…"
i.e. the preprint itself promises HuggingFace and nothing else. (Same "will be uploaded" tense as S2.)

**S8. Alternative hosting checked and NOT found** — OBSERVED (negative):
AWS Registry of Open Data search returned no Pisces listing; figshare's `search_for` API returned
newest-item noise rather than matches for "X-Atlas"/"Xaira"/"Pisces" (so figshare is **unchecked**,
not "checked and empty" — recorded honestly). X-Atlas/**Orion**'s processed h5ads are on Figshare+
(`doi.org/10.25452/figshare.plus.29190726`) — REPORTED by the search summary; Orion is already on disk.

**S9. Press coverage of the March launch** — REPORTED (not relied on for state):
Businesswire, "Xaira Therapeutics Launches X-Cell, Its First Virtual Cell Model, Trained on the
Largest-Ever Genome-Wide Perturbation Dataset, X-Atlas/Pisces",
`https://www.businesswire.com/news/home/20260317710096/en/…` (id encodes 2026-03-17).
This is the announcement; it is not evidence of files.

---

## A2 — Arc Institute

**S10. Arc news index** — OBSERVED via WebFetch, `https://arcinstitute.org/news`.
Posts after 2026-08-20: **"September 3, 2026 — The ancient anaerobic pathways still at work inside
our cells"**; **"September 2, 2026 — Why Jingtian Zhou is building high-resolution maps of the 3D
genome"**; **"August 20, 2026 — The 2026 Virtual Cell Challenge…"**. No VCC or data post since 08-20.

**S11. Arc challenge announcement** — OBSERVED via WebFetch,
`https://arcinstitute.org/news/virtual-cell-challenge-2026`, publication date "August 20, 2026".
Exact quotes:
> "Thursday, August 20: Validation data live, leaderboard active, submissions open."
> "Thursday, October 22: Final test set released."
> "Thursday, November 5: Final submissions due (11:59 pm UTC, check your timezone)"
> "Mid-late November: Winners announced"
> "Three of these cell lines will be used for the validation phase and the live leaderboard.
>  The remaining three are held back for final testing."
> "The grand prize is $100,000, with $50,000 and $25,000 for second and third place."
> "final rankings use an aggregate across a broader panel of six metrics"

**S12. The 2026 challenge paper** — OBSERVED via NCBI E-utilities efetch
(`work/pubmed_42648290.txt`): *Virtual Cell Challenge 2026: Benchmarking zero-shot generalization
across cellular contexts.* Arc Virtual Cell Initiative Team; Goodarzi H. **Cell. 2026 Aug 26**:
S0092-8674(26)00931-1. doi 10.1016/j.cell.2026.08.004. Online ahead of print. PMID 42648290.
Full abstract, verbatim:
> "The Virtual Cell Challenge returns in 2026 with a more demanding test of biological
> generalization: zero-shot prediction across multiple independent cellular contexts. Participants
> will build models to predict gene knockdown responses in a new Arc-generated dataset comprising
> unseen cell lines. The goal is to determine whether the best models can meaningfully close the
> gap between preclinical experimental predictions and human biology."
Four sentences; it is a front-matter piece, not a data or methods paper. ScienceDirect returned
**HTTP 403** to WebFetch, so the body text is **unread** — no claim here rests on it.

**S13. Arc GitHub tooling** — OBSERVED, GitHub REST API:
`ArcInstitute/cell-eval` `pushed_at 2026-07-27`, latest release **v0.8.2 (2026-07-27)**;
`ArcInstitute/cell-eval2` `pushed_at 2026-08-22`, 20 stars;
`ArcInstitute/state` `pushed_at 2026-07-24`, latest release **v0.11.1 (2026-06-29)**;
`ArcInstitute/cell-load` `pushed_at 2026-07-18`. Nothing pushed since the 08-30 freeze.

**S14. Arc HuggingFace datasets** — OBSERVED, `/api/datasets?author=arcinstitute`, 15 repos.
Most recently modified: `arcinstitute/opengenome2` 2026-09-02 (a genome corpus, not perturbation).
No perturbation dataset modified after 2026-07-22 (`bolero-data`). No VCC-2026 data repo exists.

**S15. STATE published in Cell** — OBSERVED, Europe PMC core record MED 42673963:
*Predicting cellular responses to perturbation across diverse contexts with State*, **Cell,
first publication date 2026-08-31**, doi 10.1016/j.cell.2026.07.052, `isOpenAccess = N`.
Abstract quote: "We further introduce Cell-Eval, a comprehensive evaluation framework that can be
used to evaluate future models." This is the journal version of the model already on disk
(`data/external/state-pretrained`, `state-vcc`).

---

## A3 — repositories (GEO / ArrayExpress / scPerturb / Zenodo / figshare)

**S16. GEO, Perturb-seq, human, PDAT 2026-08-01 → 2026-12-31** — OBSERVED, E-utilities esearch+esummary
(`work/geo_summ.json`). 19 GDS ids → 3 distinct series, all before the freeze:
`GSE343369` (2026/08/12, 4 samples, "UCSF Study 2- Perturb-seq targeting 225 genes expressed during
embryoid body formation"); `GSE344269` (2026/08/19, 4 samples); `GSE344535` (2026/08/21, 16 samples,
"Single-cell ORF perturbation screening of an enriched viral microprotein library in hTERT RPE-1 cells").

**S17. GEO, CRISPRi/CRISPRa + single cell, human, PDAT ≥ 2026-08-01** — OBSERVED. 2 distinct series:
`GSE343369` (above) and `GSE283877` (2026/08/27, 7 samples, "CRISPR activation reveals SOX5/6/9 as key
transcriptional regulators directing iPSC-derived cells toward a not[ochordal]…").

**S18. GEO, broader sweeps** — OBSERVED: `perturb*[Title] AND gse[Entry Type] AND PDAT ≥ 2026/08/01`
→ **count 11**; `CRISPR[Title] AND gse AND single cell AND PDAT ≥ 2026/08/01` → **count 2**;
`(Jurkat OR HepG2 OR HCT116 OR HEK293T OR K562) AND CRISPR AND gse AND PDAT ≥ 2026/08/01`
→ **count 14** (none is a transcriptomic screen at Perturb-seq scale; sample counts ≤ 16 throughout).
**No GEO series with a public deposit date after 2026-08-30 matched any query.**

**S19. GSE314342 identified** — OBSERVED, esummary: accession `GSE314342`, PDAT **2026/02/06**,
**577 samples**, "Genome-scale perturb-seq in primary human CD4+ T cells reveals genes regulating
T cell programs and human immune traits." Summary quote:
> "…performed genome-scale perturbations of all expressed genes across 22 million primary CD4+ T
> cells from four donors under three stimulation conditions."
This is the accession a rival names on the public leaderboard (S25) and it is the same study as the
repo's own CD4 donor source (GWCD4i).

**S20. ArrayExpress / BioStudies** — OBSERVED,
`https://www.ebi.ac.uk/biostudies/api/v1/search?query=Perturb-seq OR CRISPRi&collection=ArrayExpress`,
22 total hits, newest **E-MTAB-16124 / E-MTAB-16127, release 2026-08-06** (bacterial *whiA*/*whiB2*
CRISPRi RNA-seq). Nothing human, nothing after the freeze.

**S21. Zenodo** — OBSERVED, `https://zenodo.org/api/records?q=(perturb-seq OR CRISPRi) AND
created:[2026-08-01 TO 2026-12-31]`, total 426, newest 20 inspected: no perturbation screen; the
largest 2026-09 dataset is `ENCODE LR-RNA-seq models and expression values` (11.38 GB, unrelated).

**S22. scPerturb** — OBSERVED, Zenodo query `scPerturb`: the data record
`zenodo 13350497 — "scPerturb Single-Cell Perturbation Data: RNA and protein h5ad files", version 1.4,
created 2024-08-26, updated 2024-08-26, 54 files`. **No release since 2024-08.**

**S23. A screen using the 2025 VCC library** — searched and **not found**. There is no public index of
spacer sequences, so this is a *coverage-limited* negative: no new H1-hESC CRISPRi deposit appears in
S16–S22, and no preprint in S26 describes one.

---

## A4 — preprints since 2026-08-01

**S24. Europe PMC, `("virtual cell" OR "perturbation prediction" OR "Perturb-seq" OR "perturbation
response") AND FIRST_PDATE:[2026-08-01 TO 2026-09-08]`** — OBSERVED, hitCount **46**, 40 inspected.
The ones that matter, with abstracts fetched at `resultType=core`:

- **PPR1303976**, 2026-08-24, bioRxiv doi 10.64898/2026.08.21.746243 — *Virtual-cell models compress
  unseen intervention geometry through a target-specific generalization bottleneck.* Quotes:
  > "We show that high state-level similarity can coexist with loss of the relationships that
  > distinguish perturbations, a failure we term Intervention Geometry Compression (IGC)."
  > "Finally, same-target empirical anchoring transfers intervention identity across contexts far
  > more effectively than increasing exposure to other interventions."
  > "…support a design principle: empirically anchor intervention identity, then use models to
  > generalize anchored effects."
  Not present anywhere in this repo (grep over `docs/` and `experiments/*/*.md`: no match).

- **PPR1305187**, 2026-08-25, doi 10.64898/2026.08.24.746802 — *Dissecting context-dependent cancer
  vulnerabilities using Perturb-seq.* Quote: "We produced a proof-of-concept Perturb-seq dataset
  targeting 100 genes in 16 diverse cancer cell lines." **This is the Broad-16 set already on disk**
  (`data/external/broad16line`); the repo already cites it (287, 289).

- **PPR1297419**, 2026-08-12, doi 10.64898/2026.08.11.744177 — *Reliable single-cell perturbations
  explain and improve model performance.* Quote: "Among 7,170 perturbations from 29 datasets, 65% are
  unreliable, 11% shared, and 24% specific." Already cited in the repo (Wang et al. 2026;
  `experiments/287-novelty-audit/NOVELTY-AUDIT.md`, `experiments/294-writeup/REFERENCES.md`).

- **PPR1312644**, first indexed 2026-09-04, posted 2026-09-02, doi 10.64898/2026.09.02.748925 —
  *Systematic CRISPRi perturbation of 1,408 autism risk genes maps multilevel transcriptional
  convergence in human cortical neurons.* Quote: "we perform a CRISPRi Perturb-seq screen targeting
  1,408 ASD risk genes in human embryonic stem cell-derived immature cortical neurons and profile the
  resulting transcriptomic effects by single-cell RNA sequencing." **The only genuinely
  post-freeze item in the sweep.** I fetched `…v1.full-text` (HTTP 200, 98,061 bytes, `work/asd.html`)
  and found **no GEO/ArrayExpress/Zenodo accession and no availability sentence** in the served HTML
  (the page is JS-rendered, so this is "not verified", not "no deposit exists").

- **PPR1299277**, 2026-08-15 — *Cancer Cell Line Heterogeneity Imposes a Primary Bottleneck for Virtual
  Perturbation Screening at Scale.* Quote: "normal and primary cell lines achieve substantially higher
  prediction fidelity (mean per-cell-line PCC = 0.480, median 0.502) than cancer lines (0.296, median
  0.324; Mann-Whitney p = 3.93e-11)". Modelling paper; no new deposit.

- **PPR1305082**, 2026-08-25 — *MultiFlow: coupled flow matching …* — method, RNA+ATAC, no new data.
- Also seen, no new data for this problem: PPR1308038 (HSPC Perturb-seq, 2026-08-29), PPR1301415
  ("Method Choice, Not Biology, Determines In Silico Perturbation Results", 2026-08-19),
  PPR1294995 (leakage-aware benchmark, 2026-08-10), PPR1300526 (Perturb-ME, 2026-08-18),
  MED 42664972 (the CD4 screen formally published in *Cell*, 2026-08-28, doi 10.1016/j.cell.2026.08.002).

**S25. Europe PMC, `"Virtual Cell Challenge"`, 2026** — OBSERVED, **hitCount 3**: MED 42648290 (S12),
MED 42208532 (*Computational blueprints for cell fate programming*, 2026-05-28), MED 41795653
(*AI-enabled multi-scale virtual cell*, 2026-03-01). **No participant write-up is indexed.**

**S26. Europe PMC, `"cell-eval" OR "perturbation discrimination score" OR "pds_cosine"`, 2026** —
OBSERVED, **hitCount 4**: MED 42673963 (S15), MED 42167225, MED 41795653, PPR1255210
(*Elucidating the Design Space of Generative Models for Single-Cell Perturbation Prediction*, 2026-06-18).
**No indexed paper reports a score on the 2026 challenge's metric.**

---

## B — the board's frontier, public sources

**S27. The challenge leaderboard API is public and unauthenticated** — OBSERVED.
`curl https://virtualcellchallenge.org/api/leaderboard` → HTTP 200, `application/json`, 945,747 bytes,
no credentials sent (`work/lb_live_306.json`, fetched 2026-09-07T15:48:55Z). 769 entries. Each entry
carries the six member scores plus free-text `model_name` and `description` written by the team.
Read numbers (`work/board_summary_306.txt`):
```
leader           cqawesome / m0907h  0.2393600   submitted 2026-09-07T15:29:17Z
ours (rank 3)    Ahmet Dedeler / "You can do it too."  0.2160944
gap                                             0.0232656
entries >= 0.220 : 2   (cqawesome 0.23936, Aginglab.com 0.22652)
entries >= 0.200 : 9   of which with ANY description: 1
descriptions non-empty: 292 / 769; >= 120 chars: 198
best-scoring entry with a >=120-char description: rank 42, 0.17543
```
Everything a team writes about itself is self-reported and unverifiable — treated as a **claim**.

**S28. Verbatim rival descriptions** (from S27; team names as the board prints them):
- rank 6, 0.20364, `Slurm Shadi / Cellulation`: "Apparently the simulation has a leaderboard."
  (the only description in the top 30; not a method).
- rank 42, 0.17543, `Sam Duan Insilico HK Biologist / ihkb-zdfl-k20c3-m042-d10s`: "CHAMPION: the
  rank-33 display recipe (kink mag 0.42, conc 1.0, agree 0.6) plus --exact-colsum (largest-remainder
  exact column sums; ~+0.0015 level and emit-seed noise eliminated, 12/12 bit-identical)."
- rank 55, 0.16363, `Evidoa.com / mx21-c7-research`: "Gene-local compositional residual over
  multi-donor baseline. Public-only training. Research transfer to A/B/C unconfirmed."
- rank 78, 0.14954, `Testh / v6loco_plus_H1hESC_nb_a2`: "v6-LOCO transfer (whiten, all-source
  basal-cosine^1.65 mean, alpha 2/4, NB emission, seed 0) with the VCC2025 H1 hESC context added to
  the evidence pool (**21/300 panel targets**). Internal 5-fold selection score 0.1960 vs incumbent
  0.1929."  ← the claim tested in Part C.
- rank 88, 0.14200, `Ideoscopic / cfuse-w050`: "Centred fusion of two genome-wide CRISPRi screens at
  the column level: per-target response columns from Replogle K562 gwps and from **GSE314342** primary
  CD4+ T cells, each centred on its own covered-panel mean, combined as an equal-weight mix of unit
  vectors and restored to the K562 column's norm, then locally decorrelated against the pool of all
  296 real columns."
- rank 92, 0.14008, `Yusif Mukhtarov / Fusion v6 whitened + calibrate`: "Donor-pruned unit-norm fusion
  of public K562 genome-wide CRISPRi, **X-Atlas HCT116**, and **primary CD4 (Rest/Stim8/Stim48)**
  responses; cross-donor family averaging, common-response centering and Gram row-decorrelation;
  library-preserving coherent raw-count transport at footprint-calibrated per-context amplitude."
- rank 109, 0.13085, `liang zuan / AtlasConsensus-2026-09-06-v5`: "Uses public **K562, HCT116,
  HEK293T and primary CD4** responses plus official context controls."
- rank 84, 0.14452, `Life Code / PetaCell-v4.4`: "a same-cell-line screen reproduces itself at r 0.25
  across all genes but 0.85 once restricted to expressed ones, and 44 percent of our delta energy sat
  on genes where even the same line cannot reproduce[…]".
- rank 80, 0.14793, `tydymy / pseudocount-1e5`: "A reference pseudobulk holds only 155 to 204 cells…".

**S29. Public GitHub repos referencing the challenge** — OBSERVED, GitHub search API
(`q=virtual+cell+challenge&sort=updated`), total_count 36, top 25 read. Actively pushed in the last
10 days: `biona001/Virtual-Cell-Challenge-2026` (2026-09-07), `ENNCELADUS/gene-perturbation-prediction`
(2026-09-07), `mmyen/virtual-cell-challenge-2026` (2026-09-06), `saberhq/sidechain` (2026-09-04),
`algaebrown/virtualcellchallenge` (2026-08-29), `forrestsheldon/virtual-cell-challenge-2026`
(2026-08-29), `wuki-creator/magworldlatest` (2026-08-28), `litholight/vcc2026-hpc-platform` (2026-08-15).
READMEs fetched from `raw.githubusercontent.com`. Disclosed scores, verbatim:
- `biona001`: a dated submissions table — "August 30, 2026 | b01A | 213 / 473 | **0.0184**";
  "August 31, 2026 | b02 | 168 / 511 | **0.0766**"; "September 3, 2026 | b03 | 236 | **0.0649**".
  Method line for b02: "calibrated NB + full shrunken K562 fold changes (|z| > 1)".
- `wuki-creator/magworldlatest`: "The published v4 submission `magworld_h1_v4_directed_cv5_ensemble`…
  At the latest recorded status, it ranked 219 with an overall scaled score of **-0.0448**." and
  "None of these switches improved the five-fold aggregate consistently."
- `saberhq/sidechain`: "a solo entry to the Virtual Cell Challenge 2026, built in the open" — a model
  ladder (GLY/ALA/SER/CYS/…), **no score disclosed**.
- `ENNCELADUS`: "Status (2026-09-06): … Historical Exp13 Stage 2 reached test macro per-gene Spearman
  0.0225 versus context-PCA ridge 0.0851 and nearest-line 0.0462, a negative point estimate."
- `mmyen`, `forrestsheldon`, `algaebrown`, `litholight`: scaffolding only, no results.
**No public repository discloses a score above 0.08.**

**S30. Community forum** — REPORTED (search summary): the challenge's participant forum is **Discord**;
"I found no evidence of a Discourse-based forum". Discord requires an account, so it is **out of
scope** under this mandate and was not accessed.

**S31. The one public write-up found** — OBSERVED via WebFetch,
`https://gmdbioinformatics.substack.com/p/arc-virtual-cell-challenge-has-the`,
"Arc Virtual Cell Challenge: Has the Evaluation collapsed?", **published November 6, 2025**, author
`gmdbioinformatics`. It concerns the **2025** challenge ("deadline stated as 'mid November 2025'") and
its claims ("the organizers decided to cap the penalty for a poor MAE score at zero") are about the
2025 metric set. Already covered by `docs/COMPETITIVE-INTEL-2026-08-28.md` §3.3. Not new.

**S32. Organizer-suggested datasets page** — ATTEMPTED, `https://virtualcellchallenge.org/app/datasets`.
WebFetch returned only "Loading…" (client-rendered SPA). **Unread.** A rival's README cites it as
"Data suggested by organizers" (S29, `biona001`), which is REPORTED, not verified.

---

## C — on-disk files used for the falsifiable check

All read-only. Script: `work/partC_306.py`. Runtime 0.03 s, peak RSS 29.9 MB.

| file | sha256[:16] | rows |
|---|---|---|
| `data/challenge/pert_counts.csv` | — (300 unique `target_gene`) | 301 lines |
| `data/external/vcc2025/train/pert_counts_Training.csv` | `633d202be221418b` | 150 targets |
| `data/external/vcc2025/validation/pert_counts_Validation.csv` | `9da5dec8c25aff1b` | 50 targets |
| `data/external/vcc2025/test/pert_counts_Test.csv` | `c4f0ba439d7c3121` | 100 targets |
| `data/donors/h1_2025/perts.txt` | `977b808e723c73ee` | 255 perturbations |
| `experiments/301-metadata-acquisitions/work/vcc2025_spacers.csv` | `d04e2ca81092523e` | 348 guide rows |

Repo statements checked against: `docs/READ-FIRST.md` §6 ("the 2026 panel shares only 25/300 targets
with the 2025 H1 library"), `docs/LOG.md:4085` (same, dated 2026-09-06 19:20 UTC+3).

---

## Sources list (URLs)

- https://huggingface.co/api/datasets/Xaira-Therapeutics/X-Atlas-Pisces
- https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces
- https://huggingface.co/api/models/Xaira-Therapeutics/X-Cell
- https://huggingface.co/Xaira-Therapeutics/X-Cell
- https://api.github.com/orgs/xaira-therapeutics/repos
- https://api.biorxiv.org/details/biorxiv/10.64898/2026.03.18.712807
- https://www.businesswire.com/news/home/20260317710096/en/Xaira-Therapeutics-Launches-X-Cell-Its-First-Virtual-Cell-Model-Trained-on-the-Largest-Ever-Genome-Wide-Perturbation-Dataset-X-AtlasPisces
- https://arcinstitute.org/news
- https://arcinstitute.org/news/virtual-cell-challenge-2026
- https://pubmed.ncbi.nlm.nih.gov/42648290/ · https://www.cell.com/cell/fulltext/S0092-8674(26)00931-1
- https://www.sciencedirect.com/science/article/pii/S0092867426009311 (HTTP 403, unread)
- https://api.github.com/repos/ArcInstitute/{cell-eval,cell-eval2,state,cell-load}
- https://huggingface.co/api/datasets?author=arcinstitute
- https://www.ebi.ac.uk/europepmc/webservices/rest/search (S15, S24–S26)
- https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{esearch,esummary,efetch}.fcgi (S16–S19)
- https://www.ebi.ac.uk/biostudies/api/v1/search (S20)
- https://zenodo.org/api/records (S21, S22)
- https://www.biorxiv.org/content/10.64898/2026.09.02.748925v1.full-text (S24, ASD screen)
- https://virtualcellchallenge.org/api/leaderboard (S27, S28)
- https://virtualcellchallenge.org/app/datasets (S32, unread)
- https://api.github.com/search/repositories?q=virtual+cell+challenge (S29)
- https://raw.githubusercontent.com/{saberhq/sidechain,ENNCELADUS/gene-perturbation-prediction,mmyen/virtual-cell-challenge-2026,wuki-creator/magworldlatest,algaebrown/virtualcellchallenge,forrestsheldon/virtual-cell-challenge-2026,biona001/Virtual-Cell-Challenge-2026,litholight/vcc2026-hpc-platform}/{main,master}/README.md
- https://gmdbioinformatics.substack.com/p/arc-virtual-cell-challenge-has-the (S31)
