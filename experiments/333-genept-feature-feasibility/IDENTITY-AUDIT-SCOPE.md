# 333 — representation-identity audit before any biological use

September 9, 2026. The canonical schema/coverage conversion passes, but
independent readback found **734 exact-vector groups covering 1,545
distinct canonical genes**, including **34 challenge targets**. Examples:
PRKACA/PRKCA, PDPK1/PDK1, POR/VDAC2, BRD1/BRPF1, EZH2/EZH1.
These are exact float32-vector duplicates, not high cosine similarities.
They do not alone prove wrong identities, but require source verification.

The feature artifact is **quarantined for biological fitting** pending that
check. Its original coverage gate remains a pass under the disclosed
canonical conversion; a data-quality concern does not retroactively change
that number.

Permitted next diagnostic: inspect only the two description JSON members
already present in the verified archive (`NCBI_summary_of_genes.json` and
`NCBI_UniProt_summary_of_genes.json`), their schemas, primary identifiers and
description hashes for duplicate-vector groups. First inspect schema only.
Determine whether the release distinguishes primary symbols/IDs from aliases
and whether identical vectors correspond to identical supplied text.

This extends the earlier embedding-member-only acquisition for **identity
verification**, not for hand-revising a model. No biological response or
accuracy metric, manual gene-direction annotation, external embedding API,
alias-based substitution, model fitting or submission. Retain ambiguous
identities as unresolved. A correction, alternative artifact or text feature
construction needs its own frozen identifier policy before use.
