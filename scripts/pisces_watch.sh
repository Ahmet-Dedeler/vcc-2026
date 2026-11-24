#!/usr/bin/env bash
# Poll HuggingFace for the X-Atlas/Pisces data release (Xaira: genome-scale CRISPRi in Jurkat rest/active, HepG2,
# HCT116, HEK293T, iPSC -- same-line data for context A and likely October lines). As of 2026-09-13 the repo holds
# only README/LICENSE. Writes data/pisces-watch.log; drops data/PISCES-RELEASED.flag the moment data files appear.
#   nohup bash scripts/pisces_watch.sh > /dev/null 2>&1 &
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
INTERVAL="${PISCES_WATCH_INTERVAL:-21600}"   # 6 h
while true; do
  N="$(curl -sf --max-time 60 'https://huggingface.co/api/datasets/Xaira-Therapeutics/X-Atlas-Pisces?full=true' \
       | python3 -c 'import json,sys; d=json.load(sys.stdin); s=[x["rfilename"] for x in d.get("siblings",[])]; print(len([f for f in s if not f.endswith((".md",".gitattributes"))]))' 2>/dev/null || echo NA)"
  printf '%s data_files=%s\n' "$(date -u +%FT%TZ)" "$N" >> data/pisces-watch.log
  if [ "$N" != "NA" ] && [ "$N" -gt 0 ]; then printf '%s X-Atlas/Pisces has %s data files\n' "$(date -u +%FT%TZ)" "$N" > data/PISCES-RELEASED.flag; fi
  sleep "$INTERVAL"
done
