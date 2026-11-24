#!/bin/bash
# 803 batch (the single labworker job): setup, then each phase file's commands N at a time; a phase with any failed command
# stops the batch (known answers first: phase 1 holds the KA runs and the KA comparisons). Logs: work/out/logs/.
#   bash experiments/803-de-engine/cloud/batch803.sh <parallel> <phase file> [<phase file> ...]
set -u
X=experiments/803-de-engine
O=$X/work/out; mkdir -p $O/logs
L=$O/batch.log
log(){ echo "$(date -u +%FT%TZ) $*" | tee -a "$L"; }
P=$1; shift
# progress for the Mac (labworker only pushes logs at the end): batch.log + summary.tsv every 2 min
if [ -n "${BUCKET:-}" ] && [ -n "${OUT:-}" ]; then
  ( while :; do aws s3 cp "$L" "s3://$BUCKET/$OUT/progress/batch.log" --only-show-errors 2>/dev/null
      [ -f $O/summary.tsv ] && aws s3 cp $O/summary.tsv "s3://$BUCKET/$OUT/progress/summary.tsv" --only-show-errors 2>/dev/null
      sleep 120; done ) &
fi
log "batch803 start: parallel $P, phases $*; $(nproc) vCPU, $(awk '/MemTotal/{printf "%.0f GiB", $2/1048576}' /proc/meminfo)"
if [ ! -f $X/work/SETUP-OK ]; then
  bash $X/cloud/setup803.sh > $O/logs/setup.log 2>&1 || { log "setup FAILED"; tail -20 $O/logs/setup.log | tee -a "$L"; exit 1; }
  touch $X/work/SETUP-OK; log "setup OK"
fi
export NUMBA_NUM_THREADS=${NT803:-3} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 POLARS_MAX_THREADS=2
cat > $O/one.sh <<'EOS'
#!/bin/bash
ph=$1; i=$2; shift 2; c="$*"; o=experiments/803-de-engine/work/out
s=$(date +%s); /usr/bin/time -v -o $o/logs/$ph-$i.time bash -c "$c" > $o/logs/$ph-$i.log 2>&1; rc=$?
e=$(date +%s); r=$(awk -F': ' '/Maximum resident/{print $2}' $o/logs/$ph-$i.time); cat $o/logs/$ph-$i.time >> $o/logs/$ph-$i.log; rm -f $o/logs/$ph-$i.time
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$ph" "$i" "$rc" "$((e-s))" "${r:-NA}" "$c" >> $o/summary.tsv
echo "$(date -u +%FT%TZ) [$ph-$i] rc=$rc wall=$((e-s))s rss_kb=${r:-NA} :: $c" >> $o/batch.log
EOS
chmod +x $O/one.sh
for f in "$@"; do
  ph=$(basename "$f" .txt)
  n=$(awk 'NF && $1 !~ /^#/' "$f" | wc -l)
  log "phase $ph: $n commands"
  awk 'NF && $1 !~ /^#/{printf "%s\t%d\t%s\n", PH, ++i, $0}' PH="$ph" "$f" | while IFS=$'\t' read -r p i c; do printf '%s\0%s\0%s\0' "$p" "$i" "$c"; done \
    | xargs -0 -n 3 -P "$P" bash $O/one.sh
  nf=$(awk -F'\t' -v p="$ph" '$1==p && $3!=0' $O/summary.tsv | wc -l)
  log "phase $ph done: $nf failed"
  [ "$nf" -eq 0 ] || { log "STOP: phase $ph had failures"; exit 3; }
done
log "batch803 END OK"
