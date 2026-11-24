# 058 — PerturbMatch K562-GWPS reanalysis

## Question

Does Genentech's PerturbMatch reanalysis of the Replogle K562 genome-wide
CRISPRi screen improve VCC-panel effect direction, either as a replacement for
the incumbent K562 donor or as an ensemble with it?

This is a **processing replicate, not a biological replicate**. Both K562
tables derive from the same Replogle assay. Therefore the valid held-out lines
are Orion HCT116 and HEK293T only; scoring PerturbMatch against the incumbent
K562 table as held-out truth would be leakage.

## Authoritative source

- Zenodo record: https://zenodo.org/records/21232999
- DOI: https://doi.org/10.5281/zenodo.21232999
- File: `lfc_matrix_ALL.csv.gz`
- Exact size: `748,423,850` bytes
- MD5: `0017076a5cc7b9979febf3bd19bd2157`
- License: **CC BY 4.0**
- Analysis repository:
  https://github.com/Genentech/perturbmatch_analysis
- Pinned analysis commit:
  `41e1040ddefa381133b293911c55bc9218a2e7dc`
- Raw source assay: GEO `GSE146194`
- Paper: Yeung et al. (2026), DOI
  https://doi.org/10.64898/2026.07.10.737863

The source is a gzip CSV with one row per measured gene and one column per
perturbation. The first field is `gene`; all 9,471 perturbation columns use the
form `<SYMBOL>_K5`. The values are limma-voom log-fold changes from
PerturbMatch-selected controls under the combined singlet+multiplet `ALL`
condition. A complete decompression count gives **9,761 gene rows** (the
separate processed SCE is the artifact described as top-3,000 HVGs), and the
matrix covers exactly **251/300** VCC targets.

## Extraction

`extract.py` streams the compressed object directly from Zenodo, verifies its
compressed byte count and MD5 while reading, and writes:

```text
data/external/perturbmatch/
  delta.npy       float32 Fortran-order, perturbation x gene
  genes.txt
  perts.txt
  panel_perts.txt
  meta.json
```

It does not retain the 748 MB gzip or materialize the decompressed CSV.

## LODO protocol

- Full retrieval gallery: every VCC target shared by PerturbMatch, incumbent
  K562, Orion HCT116, and Orion HEK293T.
- Score/operator axis: official VCC genes shared by all four banks, minus
  **all 300 panel target coordinates**.
- Valid outer folds: hold out HCT116, then hold out HEK293T.
- Held truth: raw, uncentered held-line effects.
- Source centering and whitening basis: source perturbations excluding all 300
  panel targets.
- Candidate-independent geometry: gamma-0.65 whitening fitted only on the
  incumbent K562 donor and the non-held Orion donor. PerturbMatch never changes
  the operator in the primary comparison.
- Comparisons: direct K562 vs direct PerturbMatch; K562 replacement in the
  two-donor baseline; and fixed 10/25/50/75/90% K562-slot blends.
- Every query is ranked against the complete gallery; there is no query-only
  or top-k shortcut.

No challenge submission, count emission, or `cell-eval2` run is part of this
experiment.

## Result

Complete. See `RESULTS.md`. The direct PerturbMatch estimate is weaker than the
incumbent K562 estimate, but a 50:50 processing ensemble improves both Orion
holdouts. Under source whitening, the full two-donor model moves from
`0.794525` to `0.799693` macro PDS (`+0.005168`). This is a follow-up candidate,
not authority to change v9.
