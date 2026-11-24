---
license: other
license_name: mixed-source-licences
license_link: https://github.com/T0-RSI/ai4sci-tasks/blob/main/virtual-cell-zeroshot/environment/data/README.md
pretty_name: ai4sci virtual-cell-zeroshot release v1
tags:
  - single-cell
  - perturbation
  - crispri
  - virtual-cell
  - zero-shot
size_categories:
  - 1M<n<10M
---

# ai4sci virtual-cell-zeroshot — release v1

Prepared data for the [`virtual-cell-zeroshot`](https://github.com/T0-RSI/ai4sci-tasks/tree/main/virtual-cell-zeroshot)
task: predict single-cell CRISPRi knockdown responses in cellular contexts a model has never seen
perturbed (the Arc Virtual Cell Challenge 2026 zero-shot setting), rebuilt from public data.

```
dev/        what the agent sees (mount read-only at /workspace/data)
  train/{k562,jurkat,hct116,hek293t}/   cells.h5ad, pseudobulk.h5ad, se_embeddings.npy, dev_split/
  test/{ctx_near,ctx_mid,ctx_far}/public/   basal.h5ad, requests.jsonl, context.json
  priors/     ESM-2 gene embeddings, GO annotations, DepMap gene effects (training lines)
  manifest.json, splits.json, gene_panel.csv, gene_panel_intersection.json
private/    grader only: responses of the test contexts, request subsets/strata, cached no-effect + DE tables
manifest.json
```

`dev/se600m/` is **not** included: fetch SE-600M (`se600m_epoch16.ckpt`, `config.yaml`,
`protein_embeddings.pt`, licences) from `arcinstitute/SE-600M` at revision `5a9a80f` into
`dev/se600m/` (the task's `materialize.sh` / `fetch_release.sh` do this).

Schemas, split design, sampling caps and QC: the task's `environment/data/README.md`.
Every source file, size and SHA-256: `manifest.json`.

## Sources and licences

Derived (harmonised, subsampled, re-split) from: Replogle et al. 2022 (Figshare 20029387, CC BY 4.0);
Nadig et al. 2025 (GEO GSE264667); X-Atlas/Orion, Xaira Therapeutics (HF `Xaira-Therapeutics/X-Atlas-Orion`
@ 53a5bc98, **CC BY-NC-SA 4.0** — the `hct116`/`hek293t` trees inherit that licence: non-commercial,
share-alike, attribution); Arc Institute Virtual Cell Challenge 2025 H1 hESC data (Arc VCC data terms);
DepMap 24Q4 (CC BY 4.0); Gene Ontology annotations (CC BY 4.0). Please cite the original datasets.
