#!/usr/bin/env bash
# 811 C0 board-read dump hook (690g's DEF_DAY_DUMP_HOOK interface; UNAPPLIED: def_day uses it only if DEF_DAY_DUMP_HOOK points here).
#   bash experiments/811-lolo-model/dumphook811c0.sh <bundle> <work> <dump.npz> <identity.json> <hookdir>
# Writes <hookdir>/<basename dump.npz>: the UNSWAPPED stack dump with 811's c0_l25 blend on its bulk rows, i.e. C0 = R-ABC2's recipe
# + 811 (the swap-OFF branch of SIXA rule 3). It is dumphook811.sh WITHOUT the swap707 step: hook811.py (unchanged, lambda 0.25)
# applied straight to <dump.npz>. HOOK811_LAMBDA=0 returns the input's rows unchanged (the KA-off path: the output's content
# then equals the input dump's). Carriers, perts, genes and env are copied unchanged. Needs 811's compact banks
# (work/banks, sha-pinned) and the bundle's pert_counts.csv; no network, no identity (A/B/C board read: no held line).
set -euo pipefail
B="$1"; W="$2"; DUMP="$3"; ID="$4"; HD="$5"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
sha16() { if command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -c1-16; else sha256sum "$1" | cut -c1-16; fi; }
mkdir -p "$HD"
OUT="$HD/$(basename "$DUMP")"
echo "dumphook811c0: bundle $B; work $W; identity $ID (unused); in $DUMP ($(sha16 "$DUMP")); lambda ${HOOK811_LAMBDA:-0.25}; no swap"
uv run python experiments/811-lolo-model/hook811.py --dump-in "$DUMP" --dump-out "$OUT" --bundle "$B" --lam 0.25
echo "dumphook811c0: wrote $OUT ($(sha16 "$OUT"))"
