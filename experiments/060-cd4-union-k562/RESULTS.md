# 060 — CD4 lineage routing fails under source-only union whitening

Date: 2026-08-22

## Verdict

Do not add the primary-CD4 donor to context A or to any v9 successor.

The frozen lineage-matched test was strongly negative on held K562:

| candidate | PDS | delta versus HCT116 + HEK293T |
|---|---:|---:|
| baseline HCT116 + HEK293T | **0.810157** | — |
| CD4 mean3, fraction 0.50, gene-normalized | 0.667072 | **-0.143085** |
| CD4 mean3, fraction 0.50, raw diagnostic | 0.717191 | **-0.092967** |

The primary gene-normalized candidate improved only 63 of 242 targets, hurt
164, and tied 15. Its median per-target PDS change was `-0.109959`.

## Interpretation

Before whitening, the same fixed CD4 arm was directionally useful:

- baseline: `0.746905`;
- gene-normalized CD4: `0.769024` (`+0.022119`);
- raw CD4: `0.763383` (`+0.016477`).

Union whitening reverses that gain catastrophically. The baseline's operator
norm grew `4.836 -> 28.803`, while the CD4 candidates grew to `112.923` and
`121.074` before preserve renormalization. The CD4 response therefore loads
heavily on low-variance directions learned from the two source donors; those
directions do not transfer to held K562.

This closes the remaining identity-routed CD4 hypothesis. The CCLE match
`context A = Jurkat` remains useful biological interpretation, but it is not a
model-selection signal under the production geometry.

## Protocol and leakage controls

- held truth: centered raw K562-GWPS effects;
- fit sources: Orion HCT116 and Orion HEK293T only;
- fixed candidate: three-condition CD4 mean, fraction `0.50`;
- operator: source-only union covariance, gamma `0.65`;
- gallery: 242 VCC targets shared by all required sources;
- retrieval axis: 7,312 genes;
- operator axis: 17,806 genes;
- all 300 VCC target coordinates excluded;
- K562 values loaded only after operator fitting and prediction transformation;
- no K562 selection, count emission, `cell-eval2`, or challenge submission.

The run completed successfully on AWS `i-0da7fd79427a86d6c`. Peak RSS was
`19.64 GiB`; the process guard was `55 GiB`.

