#!/bin/zsh
cd "$(dirname "$0")"
score () {
  panel=$1; held=$2; op="work/ops/${panel}__held_${held}.npz"
  for i in {1..300}; do [[ -f "$op" ]] && sleep 25 && break; sleep 20; done
  [[ -f "$op" ]] || { echo "TIMEOUT $op"; return; }
  echo "=== scoring $panel/$held ==="
  uv run python score233.py --panels "$panel" --folds "$held" 2>&1 | tail -11
}
score real orion_hct116
score real orion_hek293t
score real k562_gwps
score pess orion_hct116
score pess orion_hek293t
score pess k562_gwps
echo "=== pipeline done ==="
