# 050 — CD4 effect breadth as transfer confidence

This experiment uses the official primary-CD4 supplementary tables to ask
which CD4-only metadata predicts cross-line effect transfer to K562.

It does not use challenge outcomes, run `cell-eval2`, emit counts, or alter a
submission. The K562-versus-CD4 correlations are diagnostic truth only; every
candidate confidence feature is computed from CD4.

Run:

```bash
uv run python scripts/run_with_memory_guard.py --max-rss-gib 1 \
  --min-free-percent 15 -- \
  uv run python experiments/050-cd4-transfer-confidence/run.py
```
