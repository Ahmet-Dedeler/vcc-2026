# Experiment 131 — HUES8 and hepatocyte metadata feasibility audit

**Audited:** 2026-08-30  
**Scope:** metadata only. No multi-gigabyte object, expression matrix, held VCC
outcome, scorer, model, count generator, submission artifact, commit, or push
was used.

## Verdict

**Scientific feasibility: GO. Immediate acquisition: NO-GO.**

All three lineages retain the same **91/300 VCC targets** after QC, not merely
in the shared guide library. Their response axes are broad enough to be useful:
9,519 genes for the 17.52 GB HUES8 ESC run, 9,584 for HUES8 definitive
endoderm, and 11,138 for WTC11-derived hepatocyte.

The acquisition blocker is access, not coverage:

- the HUES8 GCS object URLs return HTTP 404 to anonymous HEAD requests;
- every relevant Synapse entity returns HTTP 403 for anonymous entity
  metadata, with `canPublicRead=false`, `canView=false`,
  `canDownload=false`, and `isCertificationRequired=true`;
- no public compact signed gene-level effect table was found.

Do **not** download a 17–30 GB H5MU now. The best next move is to request
authenticated, certified Synapse read/download access and specifically ask for
the already-generated `calibrated_all_results.tsv` or
`trans_de_results.tsv.gz` artifacts. If only full objects become available,
prioritize hepatocyte first because it has the broadest response axis and
largest post-QC cell count.

## Exact target coverage

The released shared production library
[`IGVFDS3299AXST`](https://data.igvf.org/construct-library-sets/IGVFDS3299AXST/)
contains file
[`IGVFFI8270UPKB`](https://data.igvf.org/tabular-files/IGVFFI8270UPKB/).
Its published counts are:

- 14,151 guides;
- 2,061 non-NTC intended-target IDs / 2,062 gene names;
- 1,950 biological targeting IDs / 1,951 targeting gene names;
- 100 negative-control target IDs and 11 positive-control target IDs.

The committed post-QC `cells_per_target.tsv` and `guide_var.tsv.gz` summaries
show:

| Dataset/run | Post-QC cells | Biological target IDs | Target symbols | VCC overlap | VCC target cell-count range |
|---|---:|---:|---:|---:|---:|
| HUES8 ESC, current `qc_barcode_filter_none` | 130,789 | 1,950 | 1,951 | **91** | 227–1,228 |
| HUES8 definitive endoderm, `muddy_penguin` | 269,491 | 1,950 | 1,951 | **91** | 299–1,623 |
| WTC11 hepatocyte, `sara_synapse_syn74842722` | 1,053,820 | 1,950 | 1,951 | **91** | 193–5,218 |

Thus every one of the 91 overlapping VCC targets has at least 100 assigned
post-QC cells in all three current summaries.

The exact 91 symbols are in [`vcc_overlap.csv`](vcc_overlap.csv). They are
identical across the three lineages.

### ESC run-version correction

Experiment 076 identified the 17.52 GB ESC `sceptre_v1` object with 189,599
cells. The upstream repository now calls `qc_barcode_filter_none` the newer
canonical ESC run: 19.4 GB and 130,789 cells.

The older `sceptre_v1` energy-distance result independently confirms
**91/300** overlap. It reports 2,163 coordinate-preserving targeting units
mapping to 1,947 target symbols. This is a different target granularity from
the current Ensembl-collapsed QC summary and must not be mixed with its
1,950-ID / 1,951-symbol count.

## Response-axis breadth

A committed author notebook reads `mod/gene/var` directly from the three
H5MU files and records:

| Dataset | Gene response axis |
|---|---:|
| HUES8 ESC `sceptre_v1` | **9,519** |
| HUES8 definitive endoderm `muddy_penguin` | **9,584** |
| WTC11-derived hepatocyte | **11,138** |

The same notebook references broad calibrated gene-level tables and records
their shapes after excluding mitochondrial and ribosomal genes:

| Dataset | Rows × columns |
|---|---:|
| HUES8 ESC | 22,440,592 × 13 |
| HUES8 definitive endoderm | 22,600,320 × 13 |
| WTC11-derived hepatocyte | 26,305,056 × 13 |

This proves that broad author-derived gene-level effect tables exist in the
authors' working environment. It does **not** establish that those tables are
publicly downloadable.

## Access and published large objects

### HUES8 ESC

- Scout-run object:
  `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/2026_04_13/outs/sceptre_v1/inference_mudata.h5mu`
- Published catalog size: **17.52 GB**
- Current canonical object:
  `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/2026_04_13/outs/qc_barcode_filter_none/inference_mudata.h5mu`
- Published catalog size: **19.4 GB**
- Synapse dataset folder:
  [`syn74835010`](https://www.synapse.org/Synapse:syn74835010)
- Synapse cNMF folder:
  [`syn74893846`](https://www.synapse.org/Synapse:syn74893846)
- IGVF measurement sets:
  `IGVFDS0746MYRH`, `IGVFDS0956UUQN`, `IGVFDS2908DIHX`,
  `IGVFDS2940LYGK`, `IGVFDS3348KRMQ`, `IGVFDS5934EZKT`,
  `IGVFDS6247FKDR`, `IGVFDS8623ONYE`

### HUES8 definitive endoderm

- Object:
  `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/muddy_penguin/inference_mudata.h5mu`
- Published catalog size: **18.02 GB**
- Synapse dataset folder:
  [`syn74834952`](https://www.synapse.org/Synapse:syn74834952)
- Synapse cNMF folder:
  [`syn74893844`](https://www.synapse.org/Synapse:syn74893844)
- IGVF measurement sets:
  `IGVFDS0788TUIL`, `IGVFDS1313VGHL`, `IGVFDS1403KZYV`,
  `IGVFDS1437XFJK`, `IGVFDS2129VHBD`, `IGVFDS2520ZEYH`,
  `IGVFDS3434HGEI`, `IGVFDS4079VXPX`

### WTC11-derived hepatocyte

- H5MU:
  [`syn74728027`](https://www.synapse.org/Synapse:syn74728027)
- Processed folder:
  [`syn74842722`](https://www.synapse.org/Synapse:syn74842722)
- Parent project/folder:
  [`syn70518849`](https://www.synapse.org/Synapse:syn70518849)
- Published catalog size: **30.31 GB**

For all three large objects, the public repository publishes rounded GB sizes
but no exact byte size or content checksum. Anonymous access therefore cannot
independently verify bytes or MD5/SHA-256.

The Synapse permissions response means acquisition requires:

1. a signed-in Synapse account;
2. Synapse certification, as marked by `isCertificationRequired=true`; and
3. explicit READ/DOWNLOAD permission from the owning project/team.

No account, credential, ACL, or access request was changed in this audit.

## Compact artifacts

### Public and useful only as target-level metadata

| Dataset | Artifact | Bytes | What it contains |
|---|---|---:|---|
| HUES8 ESC `sceptre_v1` | `pval_edist_full.csv` | 1,200,137 | Per-perturbation energy distance and p-values |
| HUES8 ESC `sceptre_v1` | `target_by_target_matrix.csv` | 621,799 | Pairwise target-level phenotype distances |
| HUES8 endoderm | `pval_edist_full.csv` | 1,146,352 | Per-perturbation energy distance and p-values |
| HUES8 endoderm | `target_by_target_matrix.csv` | 572,887 | Pairwise target-level phenotype distances |

These files establish perturbation presence and phenotype strength but do not
contain signed downstream gene-response vectors, so they cannot replace the
H5MU for donor-effect extraction.

### Not publicly available

The author notebook names broad
`*_calibrated_all_results.tsv` files for ESC, endoderm, and hepatocyte, and a
committed builder exists for compact `trans_de_results.tsv.gz`. Neither class
of gene-level file is present in the public Git tree at pinned commit
`589f5ff606e7f9d09ad55e938b10727d206d5a60`.

The hepatocyte public Git tree contains only an energy-distance config, not
the corresponding target-level results.

## Acquisition recommendation

| Dataset | Recommendation | Reason |
|---|---|---|
| HUES8 ESC | **NO-GO now; GO after compact-table access** | 91 retained VCC targets and 9,519-gene axis, but anonymous GCS/Synapse access fails. Prefer the current run only after its response axis and compact effects are supplied. |
| HUES8 definitive endoderm | **NO-GO now; GO after compact-table access** | 91 retained VCC targets and 9,584-gene axis; public files are target-level only. |
| WTC11-derived hepatocyte | **NO-GO now; highest-priority after access** | 91 retained VCC targets, 11,138-gene axis, and 1.05M post-QC cells; no public effect artifact and Synapse download is denied anonymously. |

This resolves the source-census uncertainty: retained coverage is adequate in
all three, but none is presently a clean public acquisition. Request the
compact calibrated gene-level tables before considering any full H5MU
download.

## Evidence

- [`summary.json`](summary.json) — computed counts, overlap, response breadth,
  and compact-artifact inventory.
- [`vcc_overlap.csv`](vcc_overlap.csv) — exact overlap symbols by dataset/run.
- [`access_probes.json`](access_probes.json) — anonymous Synapse, GCS, and
  IGVF API responses.
- [`source_manifest.json`](source_manifest.json) — URL, byte count, and
  observed SHA-256 for every downloaded metadata snapshot.
- [`evidence/raw/`](evidence/raw/) — bounded source snapshots pinned to the
  upstream commit.
- [`audit.py`](audit.py) — reproducible audit; rejects responses above 5 MB.

The upstream code repository has no declared repository license. Dataset reuse
terms should be confirmed when access is granted.
