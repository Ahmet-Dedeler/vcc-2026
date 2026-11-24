#!/usr/bin/env bash
# chain1004.sh:15-19 eid(), against synthetic slot logs. No challenge data.
set -u
ROOT=$(mktemp -d)
mkdir -p "$ROOT/experiments/456-zs-riders/work" \
         "$ROOT/experiments/479-pertarget-sharpen/work" \
         "$ROOT/experiments/495-lodo-cap6/work"
# copy of chain1004.sh:15-19, ROOT-relative
eid(){ local f; for f in experiments/456-zs-riders/work/submit0928s[12].log experiments/479-pertarget-sharpen/work/submit1002s[12].log \
                         experiments/495-lodo-cap6/work/submit1003s[12].log; do
         [ -f "$f" ] || continue
         /usr/bin/grep -qE "ARMED slot [12]: $1 " "$f" || continue
         /usr/bin/grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' "$f" | tail -1 | awk '{print $4}'; return; done; echo none; }
cd "$ROOT" || exit 1
ID495A=aaaaaaaaaaBBBBBBBBBB
ID495B=ccccccccccDDDDDDDDDD
IDP2G=eeeeeeeeeeFFFFFFFFFF
# 10-02 slot 1 armed 495 as a fallback, submit failed: ARMED, no CAPTURED
cat > experiments/479-pertarget-sharpen/work/submit1002s1.log <<EOF
2026-10-01T20:00:05Z ARMED slot 1: 495 pkg=submissions/495.vcc name='taking action'
2026-10-02T00:01:10Z NOT SUBMITTING: sha mismatch
2026-10-02T00:01:10Z END
EOF
# 10-03 slot 1 actually scored 495
cat > experiments/495-lodo-cap6/work/submit1003s1.log <<EOF
2026-10-02T20:00:05Z ARMED slot 1: 495 pkg=submissions/495.vcc name='taking action'
2026-10-03T00:20:00Z CAPTURED own entry $ID495B -> data/board-entries/$ID495B.json
2026-10-03T00:20:00Z END
EOF
echo "=== case A: earlier ARMED without CAPTURED hides the later success ==="
printf 'eid 495 -> [%s]\n' "$(eid 495)"
echo "    (empty means the 10-03 id was dropped; rule1004 then sees none)"

# reboot re-arm appends a second submission into the same slot log
cat > experiments/479-pertarget-sharpen/work/submit1002s1.log <<EOF
2026-10-01T20:00:05Z ARMED slot 1: 479p2g pkg=submissions/479p2g.vcc
2026-10-02T00:27:00Z CAPTURED own entry $IDP2G -> data/board-entries/$IDP2G.json
2026-10-02T00:27:00Z END
2026-10-02T08:00:00Z ARMED slot 1: 495 pkg=submissions/495.vcc
2026-10-02T08:12:00Z CAPTURED own entry $ID495A -> data/board-entries/$ID495A.json
EOF
echo "=== case B: appended log, last CAPTURED is a different tag ==="
printf 'eid 479p2g -> [%s]  (true p2g id is %s)\n' "$(eid 479p2g)" "$IDP2G"
printf 'eid 495    -> [%s]  (true 495 id in THIS file is %s; 10-03 id %s is never read)\n' "$(eid 495)" "$ID495A" "$ID495B"

echo "=== case C: clean single-tag log still works ==="
cat > experiments/479-pertarget-sharpen/work/submit1002s1.log <<EOF
2026-10-01T20:00:05Z ARMED slot 1: 479p2g pkg=submissions/479p2g.vcc
2026-10-02T00:27:00Z CAPTURED own entry $IDP2G -> data/board-entries/$IDP2G.json
EOF
rm -f experiments/495-lodo-cap6/work/submit1003s1.log
printf 'eid 479p2g -> [%s]\n' "$(eid 479p2g)"
printf 'eid 495    -> [%s]\n' "$(eid 495)"

echo "=== fixed eid (patches/f2-eid.patch) on the same fixtures ==="
# rebuild case A, then B, then B2. The function is the patch.
eid_fixed(){ local f id
  for f in experiments/456-zs-riders/work/submit0928s[12].log experiments/479-pertarget-sharpen/work/submit1002s[12].log \
           experiments/495-lodo-cap6/work/submit1003s[12].log; do
    [ -f "$f" ] || continue
    id="$(awk -v tag="$1" '
      $0 ~ /ARMED slot [12]: / { hit = ($0 ~ ("ARMED slot [12]: " tag " ")) }
      hit && match($0, /CAPTURED own entry [A-Za-z0-9]{20}/) {
        s=substr($0, RSTART, RLENGTH); split(s, a, " "); saved=a[4]
      }
      END { if (saved != "") print saved }
    ' "$f")"
    [ -n "$id" ] && { printf '%s\n' "$id"; return; }
  done
  echo none; }
cat > experiments/479-pertarget-sharpen/work/submit1002s1.log <<EOF
2026-10-01T20:00:05Z ARMED slot 1: 495 pkg=submissions/495.vcc name='taking action'
2026-10-02T00:01:10Z NOT SUBMITTING: sha mismatch
EOF
cat > experiments/495-lodo-cap6/work/submit1003s1.log <<EOF
2026-10-02T20:00:05Z ARMED slot 1: 495 pkg=submissions/495.vcc name='taking action'
2026-10-03T00:20:00Z CAPTURED own entry $ID495B -> data/board-entries/$ID495B.json
EOF
printf 'A fixed eid 495 -> [%s]\n' "$(eid_fixed 495)"
cat > experiments/479-pertarget-sharpen/work/submit1002s1.log <<EOF
2026-10-01T20:00:05Z ARMED slot 1: 479p2g pkg=submissions/479p2g.vcc
2026-10-02T00:27:00Z CAPTURED own entry $IDP2G -> data/board-entries/$IDP2G.json
2026-10-02T08:00:00Z ARMED slot 1: 495 pkg=submissions/495.vcc
2026-10-02T08:12:00Z CAPTURED own entry $ID495A -> data/board-entries/$ID495A.json
EOF
printf 'B fixed eid 479p2g -> [%s]\n' "$(eid_fixed 479p2g)"
printf 'B fixed eid 495    -> [%s]\n' "$(eid_fixed 495)"

echo "=== vcc submit text: primary regex on a real captured log (417/446 shape) ==="
# the line the CLI prints is "submitted — entry <20>", not "entry: <20>"
SAMPLE=$(mktemp)
# excerpt of experiments/446-q175-level/work/submit0927s1.log lines 11 and 30, ids only
sed -n '11p;30p' /workspace/experiments/446-q175-level/work/submit0927s1.log > "$SAMPLE"
echo "--- sample lines ---"; cat "$SAMPLE"
echo "--- primary 'entry [A-Za-z0-9]{20}' ---"
grep -o 'entry [A-Za-z0-9]\{20\}' "$SAMPLE" | head -1 | awk '{print $2}'
echo "--- 'entry:' colon form does not match the primary pattern ---"
printf '  entry: YvVru3DAF89Jiyr34olR\n' | grep -o 'entry [A-Za-z0-9]\{20\}' || echo '(no match)'
rm -rf "$ROOT" "$SAMPLE"
