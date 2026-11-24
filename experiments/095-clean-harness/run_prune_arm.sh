#!/usr/bin/env bash
# Replicate experiment 091 on the REPAIRED harness.
#
# exp091 measured the essential-pruned whitening basis at +0.002409 avg6 (t = 19.85)
# against a size-matched random prune -- on data/harness/orion-hct116.h5ad, which we now
# know is depth-confounded (docs/HARNESS-DEPTH-CONFOUND.md). Both arms faced the same
# confound so the ranking survives, but the magnitude and the per-metric attribution do
# not. This re-runs the three arms on the thinned harness, where the reference's DE
# direction balance is 0.630 instead of 0.908.
#
# Waits for any in-flight run to finish first: the whitener fit peaks near 17 GiB and two
# at once on a 64 GB machine is how this repo previously pushed memory pressure past 100 GB.
set -Eeuo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

REAL="${REAL:-data/harness/orion-hct116-thinned.h5ad}"
while pgrep -f "031-somatic/run.py" >/dev/null 2>&1; do sleep 30; done

for arm in essential random; do
  out="experiments/095-clean-harness/out-${arm}"
  [ -s "${out}/results.json" ] && { echo "skip ${arm}, already done"; continue; }
  echo "=== ${arm} $(date -u +%H:%M:%SZ) ==="
  uv run python experiments/031-somatic/run.py \
    --real "$REAL" \
    --outdir "$out" \
    --configs experiments/095-clean-harness/grid.json \
    --whiten-axis union --jobs 1 \
    --basis-prune "$arm" \
    2>&1 | tail -40
done
echo "=== all arms done $(date -u +%H:%M:%SZ) ==="
