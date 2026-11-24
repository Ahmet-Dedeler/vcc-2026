# 334 — primary-text residual model loses target retrieval

September 9, 2026. **DOES_NOT_EARN_FURTHER_DEVELOPMENT.**
Original freeze **ba649342**; disclosed pre-fit text-source recovery
**5d42ddf5**. No model, dose, split, selected gene, sample count or gate
was tuned. No count test, package, production edit or submission.

## Engineering qualification

The original NCBI-only artifact failed its 95% target-text coverage gate
before fitting: 971/1,024 target bodies and 482/512 gene bodies.
The failure is preserved in `COVERAGE-FAILURE.md`. An exact-key,
exact-header fallback to the SAME archive's NCBI/UniProt descriptions
fills empty bodies, giving 1,021 targets and 507 genes. It never uses
333's quarantined vectors. Missing rows remain in evaluation.

Thirteen pure tests pass. All five relevant historical S5 raw-loader
arrays reproduce exactly. Whole-recipient and whole-target-block
isolation tests pass, including poisoning the actual held response.
Training and test source subsets both contain four donors.

## Measured comparison

Six reused CRISPRi lines, 1,024 targets, 512 selected coordinates;
12 outer-line/target-block splits, 32,768 matched sampled pairs per fit.
The per-gene diagonal deliberately gets all inner training pairs.
These are unit-RMS raw response directions, **not VCC MSE/NMAE or an
A/B/C score**. No H1 response labels or reserve values were read.

| arm | absolute error | squared error | own cosine | retrieval | decoy cosine |
|---|---:|---:|---:|---:|---:|
| source mean | .807217 | 1.063200 | .195532 | .790095 | .046059 |
| zero | .777202 | 1.000000 | .000000 | .500000 | .000000 |
| diagonal | .770353 | .971669 | .198735 | **.794614** | .046296 |
| source-only tree | **.758690** | .944498 | .205205 | .768977 | .056127 |
| STRING tree | .758659 | **.943577** | **.210309** | .727945 | .078966 |
| **text tree** | .759036 | .945086 | .205626 | .755822 | .062608 |
| shuffled-text tree | .759095 | .945264 | .204874 | .756397 | .061694 |
| shuffled-label tree | .770839 | .977089 | .146243 | .645381 | .074355 |

Primary retrieval contrasts, 10,000 paired draws:

| comparator | text minus comparator | target 95% interval | crossed line/target 95% interval |
|---|---:|---|---|
| source mean | **-.0342726** | [-.0390642,-.0294317] | [-.0491696,-.0194126] |
| source-only | **-.0131545** | [-.0152596,-.0111025] | [-.0189327,-.0080589] |
| diagonal | **-.0387914** | [-.0435665,-.0339615] | [-.0548501,-.0227640] |
| shuffled text | -.0005742 | [-.0027195,+.0015444] | [-.0044145,+.0034642] |

Retrieval worsens versus source mean and source-only on **all six lines**.
Versus source-only, lower-median-cell lines average -.0163957 and higher
lines -.0099133. This coarse cell-number split does not adjust for UMI
depth or separate study from line.

The text arm beats STRING, but not the useful simpler comparators.
Cosine and squared-error gates pass; target retrieval and robustness do
not. Own cosine improves alongside decoy cosine, illustrating why own
alignment alone cannot select a retrieval model. The shuffled-text
contrast supplies no evidence of useful text correspondence here.

This closes the fixed TF-IDF/SVD pair-feature residual construction,
not all semantic representations. All six competition metrics remain
unpriced and potentially movable; this result does not clear a candidate.

## Independent verification and provenance

`audit334.py` independently reconstructs training features/labels from
the raw banks and 49,152 endpoint rows from saved predictions:

- maximum training feature/label difference 3.55e-15;
- baseline difference 1.78e-15;
- retrieval and absolute-error differences exactly zero;
- other endpoint differences <=1.56e-15;
- 15 comparator/endpoint interval pairs reproduced within 1.12e-16;
- all original/recovery inputs and 12 prediction/model/training receipts
  verified.

The audit does not independently refit trees, rerun the text SVD or
replay every model inference. It is a readback/isolation audit, not an
independent biological replication.

Run 15:00:17.998–15:01:42.518 UTC: 84.52 seconds before terminal reporting.
Driver max RSS 1,579,663,360 bytes; sampled tree peak 1.51 GiB.
No heavy process remains.

- Original input SHA:
  `853661cd4656fbc6f9515fb478ebb4f89f07e48b9f078078212be73cfa99b1e7`
- Recovery input SHA:
  `d4b394c5821c3783339bd68ec3e533a563b2ecd1a6cbc4ecd6139c2bbc06f02f`
- Result SHA:
  `4d4f00115f2553970ed8b6881552c41462c9fa5e60b9aab4ae33968ec10f4896`
- Independent audit SHA:
  `55613d3d7c4d04f1bf87340811d755a5403adc641c9f777cf3ca329c5c4ca7e2`

Binary models, inputs and predictions stay locally under `recovery-work/`;
their exact hashes are in committed fit receipts. No paid API or cloud.
