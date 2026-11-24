# 340 — sparse neighboring responses buy cosine, lose retrieval

September 9, 2026. **DOES_NOT_EARN_FURTHER_DEVELOPMENT.**
Frozen at **96b7ef6a**. No graph/penalty/axis retuning or count-stage
advancement. No production change, package or submission.

The sparse STRING-embedding-neighborhood model improves squared error
and own-target cosine, but loses retrieval to the diagonal and own-only
ridge on all six lines. Its advantage over random neighborhoods does not
resolve under the crossed line/target bootstrap.

## Results

Six reused CRISPRi lines, 1,024 targets, 512 coordinates. All metrics below
describe unit-RMS raw response directions, not VCC MSE/NMAE or A/B/C scores.

| arm | absolute error | squared error | own cosine | retrieval | decoy cosine |
|---|---:|---:|---:|---:|---:|
| source mean | .807217 | 1.063200 | .195532 | .790095 | .046059 |
| diagonal | .770353 | .971669 | .198735 | **.794614** | .046296 |
| own refit | .762490 | .953757 | .202149 | .794309 | .049088 |
| graph | .761807 | .952226 | .204452 | .789407 | .052640 |
| random graph | .761576 | .951203 | .206058 | .787824 | .054447 |
| correlation graph | **.760390** | **.948726** | **.210180** | .785325 | .055610 |
| shuffled labels | .762679 | .956446 | .192503 | .778123 | .054206 |

Graph minus comparator retrieval, 10,000 paired draws:

| comparator | difference | target 95% interval | crossed line/target 95% interval |
|---|---:|---|---|
| source mean | -.00068779 | [-.00242599,+.00103784] | [-.00406751,+.00288302] |
| diagonal | **-.00520658** | [-.00649882,-.00395619] | [-.00748176,-.00282211] |
| own refit | **-.00490206** | [-.00575775,-.00406261] | [-.00641720,-.00323448] |
| random graph | +.00158353 | [+.00070942,+.00246162] | [-.00041013,+.00382958] |

Graph cosine versus diagonal rises +.00571670, crossed interval
[+.00289218,+.00878641], while decoy cosine also rises. Cosine improves
and retrieval worsens on all six lines. The model therefore does not
improve the discrimination endpoint it was required to improve.
Graph beats the correlation-neighborhood and label-shuffled controls,
but neither is a useful substitute for the failed stronger comparators.

The result rejects this fixed sparse ridge construction. It does not
prove that every graph model fails. All six competition metrics remain
unpriced and potentially movable; no candidate is cleared.

## Integrity

- Ten pure tests pass.
- Metadata coverage passes: 510/512 output genes have eight positive-cosine
  STRING-embedding neighbors. These are embedding neighbors, not claimed
  physical interactions.
- Five historical raw-loader arrays and all twelve source-mean/diagonal
  prediction replays are exact.
- Outer-line poisoning leaves fits/predictions unchanged. The linear
  four-source-subset/five-source-mean identity passes.
- Independent audit checks **60 models / 30,720 scalar ridge optima** from
  raw source reconstruction. Maximum normal-equation residual per
  training row is **3.55e-16**.
- Independent dense-operator prediction replay agrees within **3.55e-15**.
- **43,008** endpoint rows reproduce, with exact retrieval and
  absolute/squared errors; cosine discrepancy <=2.89e-15.
- Eighteen endpoint/comparator bootstrap pairs reproduce with independent
  counts-matrix reductions, maximum difference **8.33e-17**.
- Frozen inputs and model/prediction receipts verify.

The audit does not regenerate the upstream STRING embedding or independently
redo its static metadata graph selection. It validates fit optimality,
prediction implementation, endpoint arithmetic and uncertainty, not an
independent biological replication.

Run **17:28:25.154–17:28:41.110 UTC**, 15.96 seconds before reporting.
Driver peak RSS 1,474,150,400 bytes; sampled tree peak **1.32 GiB**.
Audit sampled peak .70 GiB. Both processes exited.

- Inputs SHA:
  `14b62e6668c612a8556638cd8762814e491f83d5f573de41ed9dbcd03b5d1194`
- Result SHA:
  `3bc5083a9df19cb1c52979dfa41dc6cb934eefb5ecfe14c1586109a1c7fe267d`
- Audit SHA:
  `194d6b67a322689a469440011065ef235c53650da855da4573a7046c458534bb`
