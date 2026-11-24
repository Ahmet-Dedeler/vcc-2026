# Experiment 117 results

Status: **NO-GO after one-shot held evaluation**

The repaired revision-2 evaluator was authenticated and independently reviewed
before truth access:

- protocol SHA-256:
  `32a4433cf66eb131a3dad10b2da215b48c209c2af1bda388cbf467f2e74cfd9c`;
- evaluator SHA-256:
  `3cbc673bababe27ae9e38ff737351fb571033e477c67b11b2c73203ca5a6d958`;
- 52 pretruth tests passed;
- authentication reported `truth_paths_touched: false`;
- the access receipt and result directory were absent before execution.

The one authorized run created the permanent access receipt and immutable
result directory:

```text
work/held-evaluation-v1-access-receipt.json
work/held-evaluation-v1/
```

Result binding SHA-256:

```text
65e6d3a5f9d10a8dd54c570d1320a4cd438e1bae828af53b2f6abebafb8074b3
```

## Primary result

The primary comparison was the frozen `pan_novel_f100` candidate versus exact
identity. H1 was diagnostic only and used the separately frozen
`pan_all_f100` arm.

| Surface | Raw-PDS delta | Directional-precision delta |
| --- | ---: | ---: |
| HEK293T | +0.0000578302 | 0.0000000000 |
| Jurkat gallery 1 | +0.0019620959 | +0.0001813212 |
| Jurkat gallery 2 | +0.0015050167 | +0.0002344025 |
| Jurkat gallery 3 | +0.0015719064 | +0.0003359168 |
| Jurkat gallery 4 | +0.0037681159 | +0.0002414494 |
| Jurkat gallery 5 | +0.0019397993 | +0.0003635450 |
| Jurkat gallery 6 | +0.0032441472 | +0.0002220426 |
| Jurkat macro | +0.0023318469 | +0.0002631129 |
| Cross-context macro | **+0.0011948386** | **+0.0001315565** |

The 20,000-draw target-cluster bootstrap interval was
`[+0.0008382064, +0.0015724834]`. The observed macro exceeded all 99 frozen
target-label nulls; their maximum was `+0.0004315281`.

The same-context H1 diagnostic was stronger:

- raw-PDS delta: `+0.0194677871`;
- directional-precision delta: `+0.0042635659`.

## Decision

Promotion failed because:

1. the frozen cross-context minimum was `+0.003`, above the observed
   `+0.0011948386`;
2. HEK293T directional-precision improvement was required to be strictly
   positive and was exactly zero.

The candidate is scientifically real but too weak and too context-specific to
ship. Do not tune its fraction against these held outcomes, generate counts,
or build a challenge submission from this result.
