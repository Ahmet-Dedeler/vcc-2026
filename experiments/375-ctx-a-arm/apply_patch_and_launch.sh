#!/usr/bin/env bash
# Apply the 366/373 chain patch to the production builder (only when no builder is running and the file is at HEAD bytes), then launch the 319pA build.
set -uo pipefail; cd ~/Code/virtual-cell-challenge-2026
pgrep -f 'scripts/def_pipeline.py' >/dev/null && { echo "a builder is running; refusing"; exit 2; }
B=$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)
[ "$B" = "1318cd3c7dbe8ed129a4f45dbf7e1e4ce7bed75c5f73af7174d3dd3530f118f4" ] || { echo "builder not at HEAD bytes ($B); refusing"; exit 3; }
uv run python experiments/373-blend-rows-patch/patch373.py --check >/dev/null || { echo "patch --check failed"; exit 4; }
cp experiments/373-blend-rows-patch/work/patched/make_submission_pm_cd4.py scripts/make_submission_pm_cd4.py
A=$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1); [ "$A" = "c6f366cad3d549129df0a7e1789a43264fab2fb7a43b861d677324545cd2fe4a" ] || { echo "post-copy hash mismatch $A"; exit 5; }
echo "patched builder in place: $A"; git diff --stat scripts/ | tail -1
nohup experiments/375-ctx-a-arm/run319pA.sh > experiments/375-ctx-a-arm/work/nohup.log 2>&1 < /dev/null &
echo "launched 319pA build pid $!"
