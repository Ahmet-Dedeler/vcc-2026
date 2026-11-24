#!/usr/bin/env bash
# 450d H1 queue (re-armed 04:45Z after 450c: the per-target prior with the CLEAN4 default is best, +0.0840; raw rows beat
# centred): the full mode at w 2 and w 4 with pt + CLEAN4 default on 3 seeds, the blend dose beyond w 2 (rows only, w 4 / 8 /
# 32), and w 2 pt with the Pan default (seed 0). Serial, after 451r (the last heavy job).
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/450-same-line-programme; LOGD=$X/work/logs; LOG=$X/work/drive450d.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -qE " END$|stopping|failed" experiments/451-sameline-mode/work/run451r.log 2>/dev/null && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /run451r\.sh/' | grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609252200 ] && { log "451r never ended by 22:00Z"; exit 1; }; sleep 30; done
A="--bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
pt(){ local w=$1 p=$2 s=$3; local lf="$LOGD/d_w${w}_${p}_pt_s${s}.log"; log "START full mode w=$w default=$p seed=$s"
      RUN450_CENTRE=0 RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run450c.py "$w" "$p" base "$s" $A > "$lf" 2>&1 || { log "failed"; exit 1; }
      log "END full mode w=$w default=$p seed=$s"; }
for w in 2 4; do for s in 0 1 2; do pt $w clean4 $s; done; done
pt 1 clean4 1; pt 1 clean4 2
for w in 4 8 32; do lf="$LOGD/d_w${w}_clean4_s0.log"; log "START rows only w=$w"
  RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run450b.py "$w" clean4 base 0 $A > "$lf" 2>&1 || { log "failed w=$w"; exit 1; }
  log "END rows only w=$w"; done
pt 2 pan 0
log "QUEUE DONE"
