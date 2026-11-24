# 079 — KOLF confirmation of frozen local centering

## Question

Does experiment 075's already-frozen reciprocal local query-centering rule
improve the exact production-v9 effect query against an independent KOLF2.1J
truth gallery?

This is a one-candidate confirmation, not a sweep:

```text
k = 64
beta = 0.5
minimum reciprocal degree = 2
```

There is no KOLF-facing tuning option in the CLI. The machine-readable
preregistration is `preregistered_protocol.json`.

## Strict holdout order

1. Build all 300 v9 queries from K562-GWPS, Orion HCT116, and Orion HEK293T
   with equal available-donor weights, `alpha=1.0796`, `beta=0`, no gene
   normalization, and no empirical-Bayes shrinkage.
2. Apply the exact cached production union-axis operator
   (`SHA-256 5608e8…`) at `gamma=0.65`, with all 300 VCC target-gene
   coordinates absent from the operator.
3. Remove all 300 panel coordinates from the query view, apply reciprocal
   local centering across the complete 300-query panel, and freeze both the
   identity and candidate matrices.
4. Only then hash and open KOLF. KOLF contributes held truth only: zero source
   rows, zero covariance rows, zero neighbor rows, and zero tuning decisions.
5. Score all 273 VCC targets present in KOLF against the complete 273-target
   gallery on all 17,806 shared non-panel genes.

The KOLF truth representation is its raw decoded quantized z-score, centered
per gene across the 273 compact VCC rows. It is not normalized, whitened,
blended, or used to alter the query.

## Predeclared gate

Promotion requires all three:

1. raw PDS delta is positive;
2. the target-cluster bootstrap 95% lower bound is positive; and
3. raw PDS delta is at least `+0.0015`.

The bootstrap resamples target genes as clusters. It uses 20,000 replicates
and frozen seed `7902026`.

## Local validation

Synthetic tests do not open real donors, KOLF, or the 1.2 GB operator:

```bash
PYTHONDONTWRITEBYTECODE=1 \
python3 experiments/079-kolf-local-centering/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.25 \
  --report experiments/079-kolf-local-centering/local-test-memory.json \
  -- \
  uv run --no-project --python 3.12 --with numpy==2.5.2 \
  python experiments/079-kolf-local-centering/test_run.py
```

The real run is cloud-only. Its guarded launcher can reuse only the named
already-running worker and cannot create or start compute:

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES \
INSTANCE_ID=i-0e05556f09944111f \
bash experiments/079-kolf-local-centering/aws/launch.sh
```

The launch script is intentionally not run as part of implementation.
Neither AWS script contains a VCC submission command.

## Outputs after a future authorized run

- `results.json`
- `per_target_pds.csv`
- `tests.log`
- `test-memory.json`
- `evaluate.log`
- `evaluate-memory.json`
- `code.sha256`
- `infrastructure.json`
- `artifact-sha256.txt`

All output paths remain under this experiment locally and under the dedicated
`input/output/exp079-kolf-local-centering/<run-id>/` S3 prefixes remotely.
