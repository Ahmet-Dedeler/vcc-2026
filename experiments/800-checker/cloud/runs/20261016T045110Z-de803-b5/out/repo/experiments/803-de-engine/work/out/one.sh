#!/bin/bash
ph=$1; i=$2; shift 2; c="$*"; o=experiments/803-de-engine/work/out
s=$(date +%s); /usr/bin/time -v -o $o/logs/$ph-$i.time bash -c "$c" > $o/logs/$ph-$i.log 2>&1; rc=$?
e=$(date +%s); r=$(awk -F': ' '/Maximum resident/{print $2}' $o/logs/$ph-$i.time); cat $o/logs/$ph-$i.time >> $o/logs/$ph-$i.log; rm -f $o/logs/$ph-$i.time
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$ph" "$i" "$rc" "$((e-s))" "${r:-NA}" "$c" >> $o/summary.tsv
echo "$(date -u +%FT%TZ) [$ph-$i] rc=$rc wall=$((e-s))s rss_kb=${r:-NA} :: $c" >> $o/batch.log
