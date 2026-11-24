#!/bin/bash
# Terminate each idle VCC worker as soon as ITS root-volume snapshot is completed (Ahmet 10-02: stop idle AWS spend).
# Never terminates an instance whose snapshot is not "completed"; stops on any snapshot "error"; cap 8 h.
L=~/Code/virtual-cell-challenge-2026/experiments/600-breakthrough-council/work/aws-cleanup-1002.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
PAIRS="us-east-1:i-08d90ae88c7b210ed:snap-0a8446d4abe1e3dc7 us-east-1:i-0bd65df670234f97e:snap-03c1a4a7e85a027e3 us-east-1:i-06197bcff9375a2e1:snap-01e6d6763a9557cab us-east-1:i-08300b395f88f6f55:snap-0cfaf592c887fd5ed us-west-2:i-0af19b9a537e08d4b:snap-04d7b42149d401573"
DONE=""; END=$(( $(date +%s) + 28800 ))
log "per-instance terminator armed (pid $$)"
while :; do
  left=0
  for p in $PAIRS; do
    r=${p%%:*}; rest=${p#*:}; i=${rest%%:*}; s=${rest#*:}
    case " $DONE " in *" $i "*) continue;; esac
    st=$(aws ec2 describe-snapshots --region "$r" --snapshot-ids "$s" --query 'Snapshots[0].State' --output text 2>/dev/null)
    if [ "$st" = "error" ]; then log "$s is in ERROR -- NOT terminating $i"; DONE="$DONE $i"; continue; fi
    if [ "$st" = "completed" ]; then
      out=$(aws ec2 terminate-instances --region "$r" --instance-ids "$i" --query 'TerminatingInstances[0].CurrentState.Name' --output text 2>&1)
      log "$s completed -> terminate $i ($r): $out"; DONE="$DONE $i"
    else left=$((left + 1)); fi
  done
  [ "$left" -eq 0 ] && break
  [ "$(date +%s)" -ge "$END" ] && { log "cap reached; $left instance(s) not terminated"; break; }
  sleep 60
done
sleep 90
for r in us-east-1 us-west-2; do log "== $r instances: $(aws ec2 describe-instances --region $r --query 'Reservations[].Instances[].[InstanceId,State.Name]' --output text | tr '\n\t' '; ')"; log "== $r volumes: $(aws ec2 describe-volumes --region $r --query 'Volumes[].[VolumeId,Size,State]' --output text | tr '\n\t' '; ')"; done
log "END"
