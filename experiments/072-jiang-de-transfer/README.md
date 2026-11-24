# 072 — Jiang pathway-DE transfer audit

## Question

Can the local Jiang 2025 pathway DE archive support a leakage-safe transfer
test or a new donor contribution to the VCC 2026 predictor?

The experiment is intentionally limited to:

- `data/external/jiang2025/DE_results_all_pathway.zip`;
- the existing compact K562-GWPS, Orion-HCT116, and Orion-HEK293T donor banks;
- effect-level PDS;
- one local process below 1 GiB RSS.

It never downloads the 10 GB PRiMeFlow H5AD, starts cloud infrastructure,
builds counts, runs `cell-eval2`, or creates a challenge submission.

## Archive contract

The runner fails closed unless the source archive has:

- MD5 `f077cba680a1affc599f5153d99b0e45`;
- 271 scientific DE tables;
- 218 unique perturbation targets;
- the five pathway batches `IFNB`, `IFNG`, `INS`, `TGFB1`, and `TNFA`;
- the six effect columns for A549, BXPC3, HAP1, HT29, K562, and MCF7;
- finite numeric values or literal `NA` in every `log2FC_*` field.

Every table is read directly through `zipfile.ZipFile.open`; the archive is
never extracted.

## Leakage-safe effect construction

Pathway/stimulation is treated as a batch.

1. For every pathway, cell line, and gene, estimate its mean using only
   non-VCC perturbation targets.
2. Center a VCC-panel target against that non-panel pathway mean.
3. Center a non-panel proxy target against a leave-one-target-out pathway
   mean, preventing its own response from entering its baseline.
4. If a target occurs in multiple pathway tables, average its centered effects
   equally across those pathway observations.

The development gallery consists only of non-panel Jiang targets shared with
all three established donor banks. The feature axis excludes:

- all 300 VCC target-gene coordinates;
- every non-panel proxy-gallery target coordinate;
- genes absent from Jiang or any established bank.

Any representation and blend fraction are selected on that non-panel gallery.
When K562-GWPS is held as truth, Jiang K562 is omitted from the Jiang
consensus. Missing Jiang coordinates preserve the baseline exactly.

## Panel stopping gate

A held-panel result is allowed only when all are true:

- at least 30 unique VCC-panel targets exist in the archive;
- at least 30 are shared with all three established donor banks;
- at least 100 non-panel proxy targets are shared by all banks;
- at least 5,000 leakage-safe VCC genes are mapped.

The verified archive has only nine VCC targets, so the expected scientific
outcome is an auditable negative metadata gate. The runner may complete the
non-panel development diagnostic and freeze a configuration, but it must not
read or score held-donor effects for those nine panel targets.

## Run

```bash
python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 1.0 --min-free-percent 25 --poll-seconds 1 -- \
  uv run python experiments/072-jiang-de-transfer/run.py \
  --max-rss-gib 0.90
```

Outputs are written under `out/`:

- `audit.json`
- `file_inventory.csv`
- `target_pathways.csv`
- `panel_coverage.csv`
- `compact/`
- `development_rows.csv`
- `development_summary.csv`
- `development_per_target.csv`
- `results.json`

## Test and lint

```bash
python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 1.0 --min-free-percent 25 --poll-seconds 0.25 -- \
  uv run python -m unittest \
  experiments/072-jiang-de-transfer/test_run.py -v

python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 1.0 --min-free-percent 25 --poll-seconds 0.25 -- \
  uv run ruff check \
  experiments/072-jiang-de-transfer/run.py \
  experiments/072-jiang-de-transfer/test_run.py
```
