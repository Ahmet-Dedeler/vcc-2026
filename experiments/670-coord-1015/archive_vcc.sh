#!/usr/bin/env bash
# 670 (coordinator [2d6f46], 2026-10-17): move scored, superseded submission packages to S3, then free the local copy.
# Disk stays near the 35 GiB build floor during def_day builds until 10-21, and AGENTS §10 wants ≥ 60 GiB through the final phase.
# Every package listed here is scored (SUBMITTED_LEDGER.json maps its package path to an entry id) and is not needed by a live
# driver. For each file: sha256 → manifest line (git-tracked) → aws s3 cp (bandwidth-capped like 600's backup_s3.sh, STANDARD_IA,
# sha256 in object metadata) → head-object must show the same ContentLength AND sha256 → only then rm the local copy.
# Any mismatch keeps the local file. Restore: aws s3 cp s3://$BUCKET/$PREFIX/<name> submissions/<name>, then compare the manifest sha.
#   nohup bash scripts/detach.sh bash experiments/670-coord-1015/archive_vcc.sh <listfile> <STOP YYYYmmddHHMM> > …/archive_vcc.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
LIST="${1:?list of submissions/*.vcc}"; STOP="${2:?stop YYYYmmddHHMM}"
W=experiments/670-coord-1015/work; LOG=$W/archive_vcc.log; MAN=experiments/670-coord-1015/ARCHIVED-VCC.tsv
BUCKET=<your-s3-bucket>; PREFIX=archive-20261017/submissions
ts(){ date -u +%FT%TZ; }; log(){ printf '%s %s\n' "$(ts)" "$*" >> "$LOG"; }
CFG="$W/aws-config-capped"; sed 's/^s3 =$/s3 =\
    max_bandwidth = 4MB\/s/' ~/.aws/config > "$CFG"
/usr/bin/grep -q "max_bandwidth = 4MB/s" "$CFG" || { log "FAIL could not write the bandwidth cap"; log "END"; exit 1; }
export AWS_CONFIG_FILE="$CFG"
[ -s "$MAN" ] || printf 'utc\tpackage\tbytes\tsha256\ts3\n' > "$MAN"
echo $$ > "$W/archive_vcc.pid"; log "START (pid $$): $(wc -l < "$LIST" | tr -d ' ') packages, cap 4 MB/s, stop $STOP"
while read -r f; do
  [ -n "$f" ] || continue
  [ "$(date -u +%Y%m%d%H%M)" -ge "$STOP" ] && { log "stop time reached before $f"; break; }
  [ -s "$f" ] || { log "skip $f: missing"; continue; }
  case "$f" in *695rabc2*|*694rabc*|*h1s710*|*ku38*|*sw707*|*lstr705*) log "skip $f: on the keep list"; continue;; esac
  pend="$(python3 -c 'import json,os,sys; s=json.load(open(os.path.expanduser("~/.config/vcc/state.json"))); p=s.get("profiles",{}).get("default",{}); print(len(p.get("pending_uploads", s.get("pending_uploads")) or {}))' 2>/dev/null)"
  [ "${pend:-0}" = 0 ] || { log "a vcc upload is pending: stopping before $f"; break; }
  SZ="$(stat -f %z "$f")"; SHA="$(shasum -a 256 "$f" | cut -d' ' -f1)"; KEY="$PREFIX/$(basename "$f")"
  left=$(( $(date -j -u -f %Y%m%d%H%M "$STOP" +%s) - $(date +%s) )); [ "$left" -gt 60 ] || { log "no time left for $f"; break; }
  log "copy $f ($SZ bytes, sha $SHA)"
  timeout "$left" aws s3 cp "$f" "s3://$BUCKET/$KEY" --storage-class STANDARD_IA --metadata "sha256=$SHA" --only-show-errors >> "$LOG" 2>&1 \
    || { log "copy FAILED for $f: local kept"; continue; }
  H="$(aws s3api head-object --bucket "$BUCKET" --key "$KEY" --query '[ContentLength, Metadata.sha256]' --output text 2>>"$LOG")"
  RL="$(echo "$H" | awk '{print $1}')"; RS="$(echo "$H" | awk '{print $2}')"
  if [ "$RL" = "$SZ" ] && [ "$RS" = "$SHA" ]; then
    printf '%s\t%s\t%s\t%s\ts3://%s/%s\n' "$(ts)" "$f" "$SZ" "$SHA" "$BUCKET" "$KEY" >> "$MAN"
    rm -f "$f" && log "verified (size + sha) and removed local $f"
  else log "VERIFY MISMATCH for $f (remote $RL / $RS): local kept"; fi
done < "$LIST"
log "END (disk $(df -g /System/Volumes/Data | awk 'NR==2{print $4}') GiB free)"
