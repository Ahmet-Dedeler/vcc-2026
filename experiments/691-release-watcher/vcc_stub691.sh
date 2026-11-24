#!/usr/bin/env bash
# 691 test stub for `vcc`: the first STUB_BASE list calls return the baseline (A/B/C controls only), later calls add a
# final dataset; `datasets download` copies $STUB_ZIP into -d's dir. Never touches the network.
ST="${STUB_STATE:?}"; n=$(cat "$ST" 2>/dev/null || echo 0)
case "$*" in
  "--json datasets list")
    echo $((n + 1)) > "$ST"
    if [ "$n" -lt "${STUB_BASE:-2}" ]; then
      echo '{"datasets": [{"id": "controls", "filename": "vcc_2026_controls.zip", "description": "A/B/C", "available": true, "unavailable_reason": null}]}'
    else
      echo '{"datasets": [{"id": "controls", "filename": "vcc_2026_controls.zip", "description": "A/B/C", "available": true, "unavailable_reason": null}, {"id": "final_controls", "filename": "vcc_2026_final_controls.zip", "description": "D/E/F", "available": true, "unavailable_reason": null}'"${STUB_EXTRA:-}"']}'
    fi ;;
  "--version") echo "vcc 0.2.2 (691 stub)" ;;
  "skill install") echo "stub skill install" ;;
  datasets\ download*) k=$(cat "$ST.dl" 2>/dev/null || echo 0); echo $((k + 1)) > "$ST.dl"
    [ "$k" -lt "${STUB_DL_FAILS:-0}" ] && { echo "stub: connection refused" >&2; exit 1; }
    d=""; prev=""; for a in "$@"; do [ "$prev" = "-d" ] && d="$a"; prev="$a"; done
    mkdir -p "$d"; cp "${STUB_ZIP:?}" "$d/vcc_2026_final_controls.zip" ;;
  *) echo "stub: unhandled $*" >&2; exit 9 ;;
esac
