# 411 — CD4 T cells do NOT sharpen Jurkat's DE signs; another cancer line (HepG2) helps more. NEGATIVE for a CD4-heavy context-A design

Status: **NEGATIVE** (2026-09-22, ~19:10Z). Question: with ballast making each context's DE carriers designable, should
context A (Jurkat, 041) carry a CD4-heavier design for its signs?

Off-panel test (the panel has no essential genes; Nadig 2025 Jurkat is the only same-line truth): 1,174 Jurkat essential
targets carried by ≥ 2 of the three production donors, by CD4 (GWCD4i, `extract_cd4.py` on the Jurkat target list, union
coverage 1,310) and by Nadig HepG2; 5,813 genes common to all lines; every line mean-removed over the evaluated targets
and divided by its median row norm; own target dropped. Metric: |t|-weighted sign accuracy on each target's top-K Jurkat
|delta| genes, mean over targets; paired 1,000× target bootstrap.

| predictor | sign acc @100 | Δ vs cancer consensus [95% CI] | cos |
|---|---:|---|---:|
| cancer consensus C (K562-GWPS, HCT116, HEK293T) | 0.6384 | — | 0.1328 |
| C + 0.25·CD4 | 0.6410 | +0.0026 [+0.0007, +0.0046] | 0.1356 |
| C + 1·CD4 | 0.6311 | −0.0074 | 0.1158 |
| CD4 alone | 0.5625 | −0.0759 | 0.0478 |
| C + 0.25·HepG2 (control line) | 0.6480 | +0.0095 [+0.0075, +0.0117] | 0.1454 |
| HepG2 alone | 0.6019 | −0.0365 | 0.0881 |

**Verdict.** Primary CD4 T-cell responses (DESeq2 LFCs) predict Jurkat's signs WORSE than the cancer lines do (0.56 vs
0.64 alone) and add less than an unrelated hepatoma line when mixed in. "Same lineage" loses to "more cancer lines" — the
1/k donor law again. No CD4-heavy carrier design for A. Limits: essential targets (strong responses) are not the panel's
regime; top-|delta| genes stand in for Jurkat's DE set. Reproduce: `uv run python experiments/411-jurkat-cd4-signs/run411.py`
(25 s, 2.9 GiB) after the extraction command in the commit message.
