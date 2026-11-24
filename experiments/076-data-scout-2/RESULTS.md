# 076 — Data scout 2: new donor and mirror audit

**Checked:** 2026-08-23  
**Scope:** public repository cards, trees, manifests, portal APIs, and small
metadata/DE files only. No large dataset was downloaded, no cloud job was
launched, and no submission was built or sent.

## Bottom line

Under a strict gate requiring all of the following:

1. a source-identifiable human genetic perturbation dataset;
2. a genuinely new somatic donor or cell context;
3. at least 30 of the 300 VCC targets verified from public metadata; and
4. usable public data today;

there is **one immediate new-context candidate but no fully clean new-donor
candidate**.

The immediate candidate is the released **IGVF WTC11-derived cardiomyocyte
CRISPRi screen**. Its compact filtered global-DE file contains **87/300 VCC
targets**, while the shared production guide manifest contains **91/300**.
This is a new differentiated lineage and is cheap to test from a 1.89 MB
compressed DE table. It is not an independent genetic donor: the cells were
differentiated from the WTC11 iPSC line. Treat it as a new biological context,
not a fourth somatic genotype.

The strongest genuinely independent-donor lead is the **Moonen et al. primary
CD4 screen**: three healthy donors and 7,664 CRISPRi promoter perturbations.
However, its Zenodo record is currently `access_right=embargoed`; the H5AD
redirects to login, the record exposes no file bytes, and no compact public
target manifest was found. It is therefore a **high-priority hold**, not a
current acquisition.

The newly reported Hugging Face repositories do not add a qualifying donor:

- `cfy2yue/perturbseq_normalized` is a well-documented normalized mirror, but
  its only genetic sources above the 30-target gate are X-Atlas/Orion
  (**300/300**) and Replogle (**272/300**), both already in the project.
- `altoslabs/primeflow-vcc-datasets` repackages already audited datasets.
- `Xaira-Therapeutics/X-Atlas-Pisces` still contains no data payload despite
  the accompanying X-Cell/preprint claim of 25.6 million cells.

## Ranked recommendation

### 1. GO — IGVF cardiomyocyte compact DE

Official analysis set:
https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/

Verified public artifacts:

| Artifact | Bytes | Contents |
|---|---:|---|
| `IGVFFI0830FXFI.csv.gz` | 1,887,277 | filtered global DE per perturbation |
| `IGVFFI5989UAVX.csv.gz` | 2,096,066 | earlier/unfiltered global DE |
| `IGVFFI7298IERA.csv.gz` | 6,974,330 | local DE |
| `IGVFFI7637STPX.h5ad` | 112,803,666,874 | full cell-by-gene object |
| `IGVFFI8270UPKB.csv.gz` | 348,103 | shared production guide library |

The released guide library has 14,151 guides and 2,062 unique gene symbols.
Its exact VCC intersection is **91 targets**. The filtered global-DE file
retains **1,889 perturbation targets**, of which **87 are VCC targets**.

For those 87 panel targets the compact DE table contains:

- 1,575 target–response-gene rows;
- 1 to 217 retained downstream genes per target, median 6;
- 9 to 1,584 cells per target, median 514;
- only one panel target below 30 cells.

This is sparse differential-expression output, not a full effect matrix. The
correct first experiment is therefore a bounded sparse donor/DE-prior test:
preserve the incumbent for all missing genes and targets, then test only
whether the cardiomyocyte directions improve held-donor retrieval or
cross-donor sign confidence. Do not download the 112.8 GB H5AD first.

The portal marks the analysis set and relevant files as released. The preprint
is CC-BY-NC 4.0; the API response does not attach a separate dataset license,
so citation/reuse terms should be checked before redistributing derived data.

### 2. HOLD — Moonen primary CD4 promoter screen

Official record:
https://zenodo.org/records/18876252

Public metadata describes:

- primary naive CD4+ T cells from three healthy female donors;
- anti-CD3/CD28 stimulation;
- dCas9-ZIM3-KRAB CRISPRi;
- 7,664 promoter targets;
- whole-transcriptome readout restricted to 8,001 variable genes;
- one raw-count AnnData object named
  `promoter_single_guide_hvg_raw_counts_uint16.h5ad`.

This could be an excellent independent-donor complement to the existing
four-donor Zhu/Dann CD4 atlas. It is not usable today: the record is embargoed,
the file list is hidden, the direct file URL redirects to login, and exact
panel overlap cannot be verified from the public metadata. Re-check after the
embargo changes or request access; do not plan around it for the next
submission.

### 3. HOLD — Southard Hs27/RPE1 TF CRISPRa

Official records:

- Hs27: https://zenodo.org/records/15200179
- RPE1: https://zenodo.org/records/15213619

The study activated 1,836 transcription factors in Hs27 fibroblasts and RPE1.
Both records are CC-BY-4.0 and provide processed mean-population H5ADs:

- Hs27 mean population: 1,720,196,368 bytes;
- RPE1 mean population: 1,897,755,232 bytes.

Hs27 is a genuinely new somatic context, but this is CRISPRa rather than
CRISPRi. No compact target manifest was exposed, so exact VCC overlap is not
verified under this metadata-only audit. It is a lower-priority programme or
representation-learning source, not a direct donor until a target list and
direction-transfer protocol are validated.

### 4. HOLD — remaining IGVF production lineages

Repository:
https://github.com/IGVF/tf_perturb_seq

The same 2,062-gene production TF library, with **91/300 panel overlap**, is
used across:

- HUES8 embryonic stem cells — 17.52 GB H5MU;
- HUES8 definitive endoderm — 18.02 GB H5MU;
- WTC11 cardiomyocytes — 16.65 GB H5MU in the processing repository and a
  separate released 112.8 GB portal H5AD;
- WTC11-derived hepatocytes — 30.31 GB H5MU.

These are valuable cross-lineage contexts, but only cardiomyocyte currently
has a compact released DE artifact that was verified here. The anonymous
Google Storage probes for the HUES8 paths returned 404, hepatocyte is linked
through Synapse, and the GitHub repository has no package license. Revisit
after access and reuse terms are confirmed.

## Hugging Face lead audit

### `cfy2yue/perturbseq_normalized`

Live API snapshot:

- repository SHA: `3f96df93440d35d6e2511c1406550e6863b85039`;
- created **2026-07-13**;
- last modified **2026-07-14**;
- 254 H5AD files;
- 12,718,270 cells;
- exact payload: **260,859,318,437 bytes**.

The live metadata therefore does **not** support “70.4 GB” or “updated within
the last day” as of 2026-08-23. Both the Hugging Face API and
`manifest/RELEASE_FILES.tsv` report 260.859 GB, and the repository has not
changed since July 14.

The 7.24 MB `GENE_TARGET_MAPPING.tsv` provides exact source-level overlap:

| Genetic source | Mapped targets | VCC overlap |
|---|---:|---:|
| X-Atlas/Orion HCT116 + HEK293T | 18,319 | **300** |
| Replogle K562 + RPE1 union | 9,590 | **272** |
| Frangieh melanoma | 240 | 5 |
| Norman K562 CRISPRa | 100 | 3 |
| Adamson K562 | 81 | 2 |
| Dixit K562 | 10 | 2 |
| Papalexi THP1 | 24 | 1 |
| Datlinger Jurkat | 32 / 20 by source curation | 0 |
| Nadig HepG2 + Jurkat | 2,317 | 0 |

Thus the repository is useful as a provenance and normalization manifest, but
it contributes no new qualifying somatic donor. The two large Orion H5ADs
(71.37 GB HCT116 and 119.08 GB HEK293T) are derived from the same Orion source
already used by the incumbent. The remaining 23 genetic H5ADs are mostly
scPerturb repackagings already present locally.

License boundaries are explicit and mixed: scPerturb-derived files are
CC-BY-4.0, while Orion-derived files are CC-BY-NC-SA-4.0.

### `altoslabs/primeflow-vcc-datasets`

Live API snapshot:

- repository SHA: `d773a3476d15e5d32d385b13cebf87e62580097e`;
- last modified **2026-08-11**;
- **363,630,070,836 bytes** stored.

The bundle contains standardized copies of Zhu/Dann CD4, VCC 2025 H1,
Replogle K562/RPE1, Nadig HepG2/Jurkat, Jiang, Feng, McFaline-Figueroa, and
Replogle 2020. It adds preprocessing convenience, not a new biological source.
Prior project audits already found:

- Zhu/Dann CD4: 291/300, already under active evaluation;
- Feng iPSC: 182/300, tested and rejected as an incremental donor;
- Jiang six-line pathway screen: 9/300, rejected by coverage gate;
- Nadig Jurkat: 0/300;
- Replogle and VCC 2025: already core development data.

Do not download this 363.6 GB mirror for donor discovery.

### `Xaira-Therapeutics/X-Atlas-Pisces`

Live API snapshot:

- repository SHA: `8ace335dfe4d23d8fb36cf0df55f7a643991475e`;
- created **2026-03-25**, last modified the same day;
- `usedStorage=0`;
- files: `.gitattributes`, `LICENSE.md`, and `README.md` only.

The card says “Coming Soon.” X-Cell materials describe 25.6 million cells,
seven CRISPRi screens, and 16 conditions spanning HCT116, HEK293T, HepG2,
iPSC, resting/activated Jurkat, and differentiating iPSC lineages. None of the
Pisces payload, target manifests, model weights, or inference code is public
in the linked repositories today. The already public X-Atlas/Orion source is
separate and contains only HCT116 and HEK293T.

Result: **watch, but no-go now**.

## Other checked leads

- **Frangieh-Izar melanoma Perturb-CITE-seq:** exact overlap is only
  `E2F1`, `IFNGR2`, `SLC25A13`, `SLIRP`, and `TIMM50` (**5/300**), despite
  240 mapped targets. It is a useful stress/immune-context benchmark, not a
  direct VCC donor.
- **Epidermal organoid Perturb-seq preprint (July 2026):** potentially relevant
  epithelial biology, but no source-cleared expression object or compact
  target manifest was found. Watch only.

## Final recommendation

1. **Run one bounded experiment from IGVF cardiomyocyte
   `IGVFFI0830FXFI.csv.gz`.** It is the only current lead combining a new
   context, CRISPRi, verified panel coverage above 30, and a tiny usable
   artifact.
2. **Do not acquire either large Hugging Face bundle for this purpose.**
   CellClip and PRiMeFlow mostly repackage data already used or already
   rejected.
3. **Monitor/request access to Moonen CD4.** It is the best future genuinely
   independent-donor lead.
4. **Defer Southard CRISPRa and the unreleased IGVF lineages** until exact
   target coverage, access, and direction-transfer value are established.
5. **Do not treat X-Atlas/Pisces as available.** The paper/card claim is real;
   the current public repository remains empty.

## Resource audit

All parsing was streaming or metadata-only. The largest measured process peak
was approximately **74 MB RSS**, safely below the 1 GiB ceiling.
