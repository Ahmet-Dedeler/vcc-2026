#!/usr/bin/env bash
# 499 reproductions. From the repo root. No challenge data, no vcc.
set -uo pipefail
cd "$(dirname "$0")/../../.."
exec uv run python experiments/499-october-code-review/repro/run_repro.py
