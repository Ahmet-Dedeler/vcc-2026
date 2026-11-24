# 084 — Listwise runtime study

Read-only performance analysis of the pinned implementation in
`experiments/080-listwise-metric-adapter/`.

This experiment:

- does not modify experiment 080;
- does not open held panel truth or donor data;
- does not launch cloud work;
- does not generate counts or submit;
- writes only inside this directory;
- fails closed if the audited `run.py` SHA-256 changes;
- runs its synthetic benchmark below 1 GiB RSS.

The study counts the exact model-fit graph and every explicit dense matrix
multiply under the default protocol. It separately reports QR/SVD call counts
because their hardware FLOP count depends on the linked LAPACK
implementation. It then benchmarks a candidate/null-invariant prepared-fit
cache against the unchanged experiment-080 optimizer.

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 \
uv run --no-project --python 3.12 --with numpy python \
  experiments/080-listwise-metric-adapter/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.1 \
  --report experiments/084-listwise-runtime-study/test-memory.json \
  -- \
  uv run --no-project --python 3.12 --with numpy python \
    experiments/084-listwise-runtime-study/test_study.py
```

```bash
PYTHONDONTWRITEBYTECODE=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 \
uv run --no-project --python 3.12 --with numpy python \
  experiments/080-listwise-metric-adapter/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.1 \
  --report experiments/084-listwise-runtime-study/benchmark-memory.json \
  -- \
  uv run --no-project --python 3.12 --with numpy python \
    experiments/084-listwise-runtime-study/study.py \
    --output experiments/084-listwise-runtime-study/benchmark-results.json
```

Measured conclusions and exact evidence are in `RESULTS.md` and
`benchmark-results.json`.
