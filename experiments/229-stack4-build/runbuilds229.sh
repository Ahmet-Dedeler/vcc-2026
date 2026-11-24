#!/bin/bash
# 229 — the three builds, STRICTLY SERIAL, in known-answer-first order.
#
# `vcc prep` peaks at 22.33 GiB and must never run alongside a build (222's measurement), so
# this script does builds ONLY.  The preps run afterwards, one at a time, from the shell.
#
# Order is deliberate: if known-answer row 1 fails, rows 2 and 3 are not worth the disk.
set -uo pipefail
cd "$(dirname "$0")/../.."
HERE=experiments/229-stack4-build
LOG=$HERE/work

run () {  # run <tag> <out> [env assignments...]
  local tag="$1"; shift
  local out="$1"; shift
  echo "=== $tag  $(date -u +%FT%TZ) ==="
  env "$@" "$HERE/build229.sh" "$out" > "$LOG/build-$tag.log" 2>&1
  local rc=$?
  echo "$tag exit=$rc"
  if [ -f "$out" ]; then
    echo -n "$tag sha256 "; shasum -a 256 "$out" | awk '{print $1}'
    echo -n "$tag bytes  "; stat -f %z "$out"
  fi
  grep -E "maximum resident set size" "$LOG/build-$tag.log" | tail -1
  return $rc
}

run ka-off      submissions/229-ka-cd4f034-f034.h5ad  VCC_DUMMY=0        || exit 1
run ka-h1pair   submissions/229-ka-h1union-g028.h5ad  VCC_H1_UNION_G=0.28 || exit 1
run stack4g035  submissions/229-stack4-g035.h5ad \
    VCC_H1_UNION_G=0.35 \
    VCC_UNION_SOURCES=dld1_lowmoi,jiang_5line_noK562,igvf_cardiomyocyte,myllia_a549 || exit 1
echo "=== ALL THREE BUILDS DONE $(date -u +%FT%TZ) ==="
