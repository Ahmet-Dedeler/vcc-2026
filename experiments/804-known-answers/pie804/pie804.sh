#!/bin/bash
# 804: PIE (Arc, replogle_xdataset checkpoint: no Replogle rows in training or evidence, so K562 is an unseen context) on the
# K562-covered targets of the lab's WORKING panels L1-L6 (801), for 702's retrodiction on fresh K562 folds. A copy of
# experiments/700-pie-transfer/cloud700_worker.sh's setup + its run_q for one query; run by labworker.sh in /work/repo.
set -euo pipefail
ts(){ date -u +%FT%TZ; }
OUTD=/work/repo/experiments/804-known-answers/pie804/out; mkdir -p "$OUTD"
export HOME=/root PATH=/root/.local/bin:$PATH
uv venv -q --python 3.12 /work/pievenv && . /work/pievenv/bin/activate && uv pip install -q arc-pie
pip_ver=$(python -c "import importlib.metadata as m; print(m.version('arc-pie'))"); echo "$(ts) arc-pie $pip_ver" | tee "$OUTD/versions.txt"
export PIE_DATA_ROOT=/work/pie/data PIE_RUNS_ROOT=/work/pie/runs PIE_CACHE_DIR=/work/pie/cache HF_HUB_DISABLE_PROGRESS_BARS=1
mkdir -p "$PIE_DATA_ROOT" "$PIE_RUNS_ROOT/replogle_xdataset" "$PIE_CACHE_DIR"
hf download arcinstitute/PIE_replogle_xdataset best_auprc.ckpt config.yaml data_stats.json --local-dir "$PIE_RUNS_ROOT/replogle_xdataset" > /dev/null
CK="$PIE_RUNS_ROOT/replogle_xdataset/best_auprc.ckpt"; echo "checkpoint sha256 $(sha256sum "$CK" | cut -d' ' -f1)" | tee -a "$OUTD/versions.txt"
python - <<'PY'
import yaml
from pie.assets import resolve_asset
cfg = yaml.safe_load(open("/work/pie/runs/replogle_xdataset/config.yaml"))["data"]
for p in cfg["preprocessed_dirs"]:
    print("resolved", resolve_asset(p, kind="preprocessed"), flush=True)
for n, p in cfg["source_dirs"].items():
    print("resolved", resolve_asset(p, kind="source"), flush=True)
print("resolved", resolve_asset(cfg["gene_text_dir"], kind="source"), flush=True)
open("/work/replogle_dir.txt", "w").write(str(resolve_asset(cfg["preprocessed_dirs"][0], kind="preprocessed")))
PY
echo "$(ts) assets resolved; disk $(df -h /work | tail -1 | awk '{print $4}') free"
cd /work/repo
/usr/bin/time -v python experiments/700-pie-transfer/infer700.py --ckpt "$CK" --prep "$(cat /work/replogle_dir.txt)" \
  --query experiments/804-known-answers/pie804/query_k562_L.json --out "$OUTD/pred_k562_L.parquet" --device cpu --batch 8 \
  --context-text data/external/pie/local_sources/context_text 2>&1 | tee "$OUTD/infer_k562_L.log" | tail -5
sha256sum "$OUTD/pred_k562_L.parquet" | tee -a "$OUTD/versions.txt"
echo "$(ts) END"
