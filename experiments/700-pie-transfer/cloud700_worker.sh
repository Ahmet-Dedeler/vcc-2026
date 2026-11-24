#!/bin/bash
# 700 cloud worker: build PIE's evidence cache with Arc's unmodified code on a machine with the RAM it needs (~55 GB peak
# by the array sizes: 9 float64/int64 (perts x genes) accumulators + 2 per-dataset int64 cubes), then run the xdataset
# checkpoint on three query sets and upload predictions + the evidence artifact. Runs as root on AL2023 via systemd-run;
# env: AWS_REGION, BUCKET, PREFIX (input/700-pie/<run>), OUT (output/700-pie/<run>).
# Self-termination: user-data arms `shutdown -h +MAX`; the instance terminates on shutdown; this script's EXIT trap powers
# off. Writes s3://$BUCKET/$OUT/STATUS = "DONE 0" | "FAILED <rc> <step>".
set -uo pipefail
WORK=/work; mkdir -p "$WORK/out"; cd "$WORK" || exit 1
LOG="$WORK/worker.log"; exec > >(tee -a "$LOG") 2>&1
STEP=start
finish(){ rc=$?; [ -s "$WORK/STATUS" ] || printf 'FAILED %s %s\n' "$rc" "$STEP" > "$WORK/STATUS"
  aws s3 cp "$LOG" "s3://${BUCKET}/${OUT}/worker.log" --only-show-errors || true
  aws s3 cp "$WORK/STATUS" "s3://${BUCKET}/${OUT}/STATUS" --only-show-errors || true
  shutdown -h now; }
trap finish EXIT
( while true; do sleep 120; aws s3 cp "$LOG" "s3://${BUCKET}/${OUT}/worker.log" --only-show-errors || true; done ) &
ts(){ date -u +%FT%TZ; }
echo "$(ts) worker start on $(nproc) cpus, $(free -g | awk '/Mem:/{print $2}') GiB"
STEP=setup
dnf install -y -q tar gzip time procps-ng > /dev/null 2>&1 || true
export HOME=/root; curl -LsSf https://astral.sh/uv/install.sh | sh > /dev/null 2>&1; export PATH=/root/.local/bin:$PATH
uv venv -q --python 3.12 /work/venv && . /work/venv/bin/activate && uv pip install -q arc-pie || exit 2
aws s3 cp "s3://${BUCKET}/${PREFIX}/inputs.tar.gz" "$WORK/inputs.tar.gz" --only-show-errors && tar xzf inputs.tar.gz || exit 3
export PIE_DATA_ROOT=/work/pie/data PIE_RUNS_ROOT=/work/pie/runs PIE_CACHE_DIR=/work/pie/cache HF_HUB_DISABLE_PROGRESS_BARS=1
mkdir -p "$PIE_DATA_ROOT" "$PIE_RUNS_ROOT/replogle_xdataset" "$PIE_CACHE_DIR"
STEP=checkpoint
hf download arcinstitute/PIE_replogle_xdataset best_auprc.ckpt config.yaml data_stats.json \
  --local-dir "$PIE_RUNS_ROOT/replogle_xdataset" > /dev/null || exit 4
CK="$PIE_RUNS_ROOT/replogle_xdataset/best_auprc.ckpt"
echo "$(ts) checkpoint $(sha256sum "$CK" | cut -c1-16)"
STEP=assets
python - <<'PY' || exit 5
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
CT=/work/inputs/context_text
run_q(){  # $1 tag  $2 prep dir  $3 query
  STEP="infer-$1"; echo "$(ts) START infer $1"
  /usr/bin/time -v python /work/inputs/infer700.py --ckpt "$CK" --prep "$2" --query "$3" --out "/work/out/pred_$1.parquet" \
    --device cpu --batch 8 --context-text "$CT" > "/work/out/infer_$1.log" 2>&1
  rc=$?; grep -E "Maximum resident|Elapsed \(wall" "/work/out/infer_$1.log" | sed "s/^/  $1 /"; tail -2 "/work/out/infer_$1.log"
  aws s3 cp "/work/out/infer_$1.log" "s3://${BUCKET}/${OUT}/out/infer_$1.log" --only-show-errors || true
  [ "$rc" = 0 ] || exit 10
  aws s3 cp "/work/out/pred_$1.parquet" "s3://${BUCKET}/${OUT}/out/pred_$1.parquet" --only-show-errors || exit 11
  echo "$(ts) END infer $1 rc=$rc"; }
run_q k562 "$(cat /work/replogle_dir.txt)" /work/inputs/query_k562.json      # builds the evidence cache first
run_q h1test /work/inputs/prep/h1test /work/inputs/query_h1test.json
run_q abc /work/inputs/prep/vcc26abc /work/inputs/query_abc.json
STEP=evidence-upload
du -sh "$PIE_CACHE_DIR"/evidence/* 2>/dev/null
aws s3 sync "$PIE_CACHE_DIR/evidence" "s3://${BUCKET}/${OUT}/evidence/" --only-show-errors || exit 12
printf 'DONE 0\n' > "$WORK/STATUS"
echo "$(ts) worker done"
