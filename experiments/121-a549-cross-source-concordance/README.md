# 121 — A549 cross-source concordance

Status: **complete; preregistered source-only NO-GO**

Date: **2026-08-30**

## Question

This experiment tests whether the independently derived Altair A549 and
Myllia A549 same-target perturbation directions reproduce strongly enough to
justify a later held-context donor experiment.

The analysis is intentionally narrow:

- exact-symbol intersection only;
- the complete 20-target shared gallery;
- one fixed source-native representation with no arm grid;
- one global complete-case response mask;
- all 300 VCC panel coordinates and all 20 common target-gene coordinates
  excluded from every metric;
- bidirectional cosine-retrieval PDS as the primary statistic;
- deterministic target-label derangements and a paired whole-target bootstrap
  frozen before either expression array is opened.

The exact metric, zero-row, midrank, sign, null, bootstrap, and GO rules are in
`preregistered_protocol.json`. `AUTHORIZATION_MANIFEST.json` hash-pins both
source artifact manifests, both compact receipts, the two expression arrays,
the two measured masks, and all text axes.

## Result

The frozen global response axis contains 7,685 exact-symbol coordinates after
requiring complete measured support in all 20 rows and excluding the union of
the 300 VCC panel genes and 20 common target genes.

The candidate failed every preregistered GO criterion:

| Statistic | Result |
| --- | ---: |
| Altair → Myllia PDS | 0.455263 |
| Myllia → Altair PDS | 0.450000 |
| Equal macro PDS | **0.452632** |
| Derangement-null p95 | 0.575000 |
| One-sided null p-value | 0.874913 |
| Whole-target bootstrap 95% interval | [0.339474, 0.569737] |
| Altair → Myllia own-target cosine | -0.029933 |
| Myllia → Altair own-target cosine | -0.029933 |
| Altair → Myllia top-20% sign precision | 0.409141 |
| Myllia → Altair top-20% sign precision | 0.485231 |
| Altair → Myllia directional sign credit | -0.143674 |
| Myllia → Altair directional sign credit | -0.031691 |

Conclusion: these independently derived A549 source directions do not justify
a held-context donor experiment under the preregistered protocol. Do not tune
this comparison against the revealed source outcome.

## Safety boundary

This is source-only. It does not access VCC held perturbation truth, existing
outcome directories, `cell-eval2`, count generation, challenge build or
submission paths. It does not commit or push.

## Run

```bash
PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/121-a549-cross-source-concordance/freeze_plans.py build

PYTHONDONTWRITEBYTECODE=1 uv run python -m unittest \
  experiments/121-a549-cross-source-concordance/test_evaluate.py

uv run ruff check \
  experiments/121-a549-cross-source-concordance/freeze_plans.py \
  experiments/121-a549-cross-source-concordance/evaluate.py \
  experiments/121-a549-cross-source-concordance/test_evaluate.py

PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/121-a549-cross-source-concordance/evaluate.py run
```

Both builders refuse to overwrite existing output directories.

Frozen outputs:

```text
work/preregistered-freeze-v1/
work/source-only-evaluation-v1/
```
