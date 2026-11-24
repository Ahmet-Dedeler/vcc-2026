#!/usr/bin/env bash
# 811 C1 board-read dump hook (690g's DEF_DAY_DUMP_HOOK interface; UNAPPLIED: def_day uses it only if DEF_DAY_DUMP_HOOK points here).
#   bash experiments/811-lolo-model/dumphook811.sh <bundle> <work> <dump.npz> <identity.json> <hookdir>
# Writes <hookdir>/<basename dump.npz>: the stack dump with
#   (1) 707's strong-tier swap FIRST: experiments/707-strong-tier-direction/swap707.py, unchanged (beta 0.25), exactly as
#       dumphook707.sh calls it, into <hookdir>/swap707/<basename>;
#   (2) then 811's c0_l25 blend on the swapped bulk rows: hook811.py (lambda 0.25; HOOK811_LAMBDA=0 returns the swapped rows,
#       which is the KA-off path: the output then equals dumphook707.sh's byte for byte).
# Carriers, perts, genes and env are copied unchanged. Needs the donor banks, 811's compact banks (work/banks, sha-pinned) and
# the bundle's pert_counts.csv; no network, no identity (A/B/C board read: no held line).
set -euo pipefail
B="$1"; W="$2"; DUMP="$3"; ID="$4"; HD="$5"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
sha16() { if command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -c1-16; else sha256sum "$1" | cut -c1-16; fi; }
mkdir -p "$HD/swap707"
OUT="$HD/$(basename "$DUMP")"; SW="$HD/swap707/$(basename "$DUMP")"
echo "dumphook811: bundle $B; work $W; identity $ID (unused); in $DUMP ($(sha16 "$DUMP")); lambda ${HOOK811_LAMBDA:-0.25}"
uv run python experiments/707-strong-tier-direction/swap707.py --dump-in "$DUMP" --dump-out "$SW" --bundle "$B" --beta 0.25
echo "dumphook811: 707 swap wrote $SW ($(sha16 "$SW"))"
uv run python experiments/811-lolo-model/hook811.py --dump-in "$SW" --dump-out "$OUT" --bundle "$B" --lam 0.25
echo "dumphook811: wrote $OUT ($(sha16 "$OUT"))"
