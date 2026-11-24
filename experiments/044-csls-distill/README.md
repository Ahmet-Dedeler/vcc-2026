# 044 — CSLS distillation without a hidden gallery

Exact CSLS is useful only as an oracle diagnostic. For cosine similarities
`S = Q T^T`, its candidate-density term is

```text
score_csls(i, j) = S(i, j) - 0.5 * mean_top_k_i S(i, j)
```

up to a positive scale and a query-specific constant. The second term requires
the hidden reference gallery `T`, so it cannot be used directly in a VCC
submission.

This experiment asks whether that candidate penalty has a stable donor-trained
representation that can be baked into submitted query vectors:

```text
affine:   q' = q + eta * b
spectral: q' = q + eta * V diag(a) V^T q
combined: both terms
```

`b`, `V`, and `a` are learned only from public donor effects. The scorer still
uses ordinary cosine on `q'`.

## Leakage controls

- Outer truth is one of K562-gwps, Orion-HCT116, or Orion-HEK293T and is never
  used to fit or tune its transform.
- Transform training uses 600 perturbations shared by K562-gwps, HCT116,
  HEK293T, RPE1, and Nadig-Jurkat.
- Those 600 perturbations have zero overlap with the 2026 VCC panel.
- Hyperparameters are selected by leave-one-donor-out validation inside the
  outer training donors.
- The final outer test uses the 259 VCC targets shared by the three main donor
  lines.
- The low-rank basis uses randomized range finding. No full gene-by-gene
  eigendecomposition, `cell-eval2`, count generation, or submission build is
  involved.

The primary number is the nested outer-test delta over ordinary cosine. The
fixed-grid LODO table is reported only as an exploratory diagnostic because
choosing a configuration from that same table would be optimistic.

```bash
uv run python experiments/044-csls-distill/run.py
```

The script memory-maps donor arrays, builds small experiment-local caches, and
prints its resident set after each cache/model stage.
