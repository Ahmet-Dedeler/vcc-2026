# 336 — biomedical semantic features do not improve transfer

September 9, 2026. **DOES_NOT_EARN_FURTHER_DEVELOPMENT.**
Freeze **5d1578af**. All 12 source-only replays are bit-exact.
No alternate encoder, parameter search, count test, package or submission.

This changes only the representation in 334's fixed residual model:
pretrained biomedical semantic vectors instead of TF-IDF/SVD. Targets,
coordinates, descriptions, source features, sample rows, labels, donor
counts, splits, tree capacity and seeds remain matched. No new held
context or independent replication is supplied.

## Measured outcomes

Six lines, 1,024 targets, 512 genes; unit-RMS source directions.
These are **not competition MSE/NMAE or A/B/C scores**.

| arm | absolute error | squared error | cosine | retrieval | decoy cosine |
|---|---:|---:|---:|---:|---:|
| source mean (334) | .807217 | 1.063200 | .195532 | .790095 | .046059 |
| diagonal (334) | .770353 | .971669 | .198735 | .794614 | .046296 |
| source-only (334) | .758690 | .944498 | .205205 | .768977 | .056127 |
| TF-IDF text (334) | .759036 | .945086 | .205626 | .755822 | .062608 |
| **semantic** | .758957 | .944691 | .206780 | .746792 | .067085 |
| shuffled semantic | .758990 | .944827 | .206021 | .752725 | .064119 |
| semantic, shuffled labels | .770477 | .976039 | .150784 | .637260 | .083364 |

Semantic retrieval contrasts, 10,000 paired draws:

| comparator | difference | target 95% interval | crossed line/target 95% interval |
|---|---:|---|---|
| source mean | **-.0433032** | [-.0485103,-.0380140] | [-.0613587,-.0256995] |
| source-only | **-.0221851** | [-.0247725,-.0197005] | [-.0307290,-.0143998] |
| TF-IDF text | **-.0090306** | [-.0111519,-.0069218] | [-.0138612,-.0042960] |
| shuffled semantic | **-.0059327** | [-.0081034,-.0038278] | [-.0101832,-.0015746] |

The semantic arm loses to source-only and TF-IDF on all six lines.
Source-only-relative retrieval is -.0273199 on the lower-cell-number
half and -.0170502 on the higher half. This is not a UMI-depth adjustment.

Cosine and squared-error gates pass; retrieval and robustness gates fail.
The correct semantic correspondence is worse than its shuffled control
on this construction. This rejects the specified feature/model pairing,
not the possibility of every future use of biomedical language models.
All six competition members remain unpriced and potentially movable.

## Independent verification

Four pure projection/missingness/model/permutation tests pass.
All 12 source-only replays reproduce the old predictions with **0.0**
maximum difference. The eight old arms reproduce their aggregate values.

`audit336.py` verifies:

- 49,152 prior endpoint rows remain numerically exact;
- 18,432 new rows independently reconstructed, 67,584 total;
- absolute error, cosine and retrieval differences exactly zero;
- squared error <=1.34e-15; decoy cosine <=5.56e-17;
- five retrieval contrast interval pairs agree within 4.17e-17;
- every registered input and model/prediction receipt is hash-verified.

The audit does not independently refit the semantic trees or replay their
PCA/inference. The isolated source-feature/label audit is inherited,
hash-bound, from 334.

Run 15:24:19.618–15:25:06.598 UTC, 46.98 seconds before terminal
reporting. Process maximum RSS 1,196,081,152 bytes; sampled tree peak
1.17 GiB. All sprint-owned jobs from 334–336 have exited.

- Input SHA:
  `af234341c1b36d7c9a5b26efdae343afbec7ce9179419cb8b466c86d1caa39dd`
- Result SHA:
  `c46a7f34ee64c6f9f3489f987608a6789039c1ec1a0dacca3eae592f676e0ca5`
- Independent audit SHA:
  `5a6fadd34d17bff7b9541b040da8ede75e9e4d14d00398a95baadce77fa6ff7f`

No reserve, paid API, cloud, worker, production change or submission.
