#!/usr/bin/env bash
# 700: work that needs ~16-18 GB and must not run beside a production build, queued through gate602 (FIFO, starts only
# when no builder runs and memory allows): (1) the mock-day hook rehearsal (mockday700.sh's body), (2) the context-specificity
# inference on Replogle (ctxspec700). Serial. Logs under work/.
#   nohup bash experiments/700-pie-transfer/queue_gap700.sh > experiments/700-pie-transfer/work/queue_gap700.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/700-pie-transfer
# coordinator's guard (10-14 23:5xZ): nothing starts while 702's build chain (r1 -> ntclevel -> r2 has builder-free moments)
# or 671's DN chain is alive; then gate602 adds its own builder and memory checks
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
until ! alive experiments/702-pie-bulk-arm/build702.sh && ! alive experiments/671-dnmag-arm/dn671.sh; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610161200 ] && { echo "$(date -u +%FT%TZ) cap reached"; echo "$(date -u +%FT%TZ) END"; exit 9; }
  sleep 60; done
echo "$(date -u +%FT%TZ) 702 and 671 chains done; now waiting for R7a (coordinator: R7a strictly first)"
# R7a's v2 wrapper must have EXITED and its log's last line must be its terminal "END rc=" (AGENTS 6.2: process exit AND line)
until ! alive experiments/552-r7-reserve/work/r7a_access_1015b.sh && tail -1 experiments/551-r6-c10b/work/r7a/access_wrapper_1015.log | /usr/bin/grep -qE "END rc="; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610161200 ] && { echo "$(date -u +%FT%TZ) cap reached waiting for R7a"; echo "$(date -u +%FT%TZ) END"; exit 9; }
  sleep 60; done
echo "$(date -u +%FT%TZ) R7a done: queueing"
bash $X/mockday700.sh > $X/work/mockday/run.log 2>&1
echo "$(date -u +%FT%TZ) mockday done: $(tail -3 $X/work/mockday/run.log | tr '\n' ' ')"
export PIE_DATA_ROOT=$REPO/data/external/pie/data PIE_RUNS_ROOT=$REPO/data/external/pie/runs PIE_CACHE_DIR=$REPO/data/external/pie/cache HF_HUB_OFFLINE=1
RD="$(ls -d $REPO/data/external/pie/data/hf/datasets/arcinstitute/pie_replogle_nadig_essential/*/preprocessed | head -1)"
GATE602_QUIET=1 /usr/bin/time -l bash experiments/602-council-gate/gate602.sh 700-ctxspec 20 -- \
  $X/.venv/bin/python $X/infer700.py --ckpt $PIE_RUNS_ROOT/replogle_xdataset/best_auprc.ckpt --prep "$RD" \
  --query $X/work/query_ctxspec.json --out $X/work/pred_ctxspec.parquet --context-text $REPO/data/external/pie/local_sources/context_text \
  > $X/work/ctxspec_infer.log 2>&1
echo "$(date -u +%FT%TZ) ctxspec infer rc=$?"
$X/.venv/bin/python $X/ctxspec700.py $X/work/pred_ctxspec.parquet > $X/work/ctxspec700.out 2>&1
echo "$(date -u +%FT%TZ) ctxspec: $(cat $X/work/ctxspec700.out | tr '\n' ' ' | cut -c1-400)"
echo "$(date -u +%FT%TZ) END"
