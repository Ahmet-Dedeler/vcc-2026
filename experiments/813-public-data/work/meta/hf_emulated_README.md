---
license: cc-by-4.0
pretty_name: K562 genome-wide CRISPRi Perturb-seq, training knockdowns
---

Per-knockdown summaries of part of the K562 genome-wide CRISPRi Perturb-seq screen of Replogle et al.,
"Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq", Cell 185:2559 (2022),
as redistributed by scPerturb (Peidli et al., Zenodo record 13350497, CC BY 4.0), with the screen's control
cells as raw counts, Reactome pathways (release 97, CC0) and Gene Ontology annotations and ontology (GO release
2026-08-05, CC BY 4.0).

The per-knockdown arrays are means (and variances) over cells of log1p(counts / cell total over the 8,248
measured genes * 10,000), derived from the scPerturb file. Changes: knockdowns were summarised per gene and a
subset of the knockdowns is included; GO annotations with NOT qualifiers and the evidence codes HTP, HDA, HMP,
HGI, HEP, IEP and ND were removed.
