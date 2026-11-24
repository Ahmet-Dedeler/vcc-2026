# 176 — LINCS / CMap 2020 as additional donor lines: inventory (coordinator, 2026-09-02 16:35)

The agent launched for this experiment never ran: Claude Code's OAuth session had expired
(`Failed to authenticate: OAuth session expired and could not be refreshed`, seven launches). The
coordinator did the gating step directly. Metadata: `data/external/cmap/{siginfo,geneinfo,cellinfo}_beta.txt`
from `https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/`.

## Gene axis
978 landmark genes (953 on our 18,533-gene axis); 11,350 inferred (10,808 on our axis).
**Only 22 of the 300 panel targets are landmarks; 235/300 are measured somewhere on L1000.**

## Panel coverage per cell line × perturbation type (distinct target genes)

| cell | pert_type | targets | **panel /300** | h1-test /85 | sigs |
|---|---|---:|---:|---:|---:|
| VCAP | trt_sh | 4,576 | **118** | 55 | 25,215 |
| PC3 | trt_sh | 4,044 | 113 | 57 | 23,343 |
| A375 | trt_sh | 3,837 | 111 | 57 | 16,645 |
| HT29 | trt_sh | 3,766 | 110 | 56 | 16,476 |
| HA1E | trt_sh | 3,812 | 109 | 57 | 14,992 |
| A549 | trt_sh | 3,771 | 109 | 55 | 16,061 |
| MCF7 | trt_sh | 3,699 | 105 | 56 | 20,786 |
| HCC515 | trt_sh | 3,532 | 100 | 51 | 13,992 |
| HEPG2 | trt_sh | 3,573 | 98 | 55 | 12,868 |
| (same 9 lines) | trt_sh.cgs (consensus) | 2,400–3,800 | 76–108 | 33–57 | one per gene |
| 9 lines (A375, A549, HT29, PC3, AGS, BICR6, ES2, U251MG, YAPC) | trt_xpr (CRISPR-KO) | ~5,115 each | **82** each (fixed library) | 41 | 10–15k |

**Union panel coverage across all lines with ≥100 targets: 127/300.** No line reaches the 150-target
gate set in the brief. Per covered target, CMap would take the donor count from 3 to ~12 on the
landmark axis (953 genes, ~9% of the DE universe) and — only if the inferred genes carry sign — on
~10,800 genes.

## Level-5 files (found by probing; the bucket listing is blocked)
- `level5/level5_beta_trt_sh_n238351x12328.gctx` — **11.8 GB**, downloading to `data/external/cmap/level5_beta_trt_sh.gctx`
- `level5/level5_beta_trt_xpr_n142901x12328.gctx` — 6.5 GB
- `level5/level5_beta_ctl_n58022x12328.gctx` — 2.9 GB

## Plan for the measurement (once the file lands)
1. Subset to `trt_sh.cgs` (consensus, seed effects cancelled) for the 9 core lines; z-scores → sign/ordering carriers.
2. For h1-test's 85 targets (50–57 covered per line): sign agreement of each CMap line with each production donor on shared (target, landmark gene) pairs — the platform control. If ~0.5, stop.
3. Sign accuracy vs the H1 reference pool on landmark genes (175's instrument), each line alone and added to the trio; then the same on inferred genes.
4. Only if 3 is positive: `pds`/cos in delivered space via 173's pipeline.

Expected ceiling if everything works: ~1/3 of targets × (9% landmark or ~100% inferred) of genes at
+0.05 sign accuracy → **+0.0025 (landmark only) to +0.017 (if inferred genes carry sign)** pooled
sign accuracy. The landmark-only case is negligible; the whole value rests on the inferred genes.
