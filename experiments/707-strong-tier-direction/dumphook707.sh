#!/usr/bin/env bash
# 707 day-of dump hook (690g's DEF_DAY_DUMP_HOOK interface; UNAPPLIED: def_day uses it only if DEF_DAY_DUMP_HOOK points here).
#   bash experiments/707-strong-tier-direction/dumphook707.sh <bundle> <work> <dump.npz> <identity.json> <hookdir>
# Writes <hookdir>/<basename dump.npz>: the stack dump with swap707.py's strong-tier direction swap (top z_src tertile of the
# bundle's panel, beta 0.25, norm-preserving toward the production donors' unit mean; carriers, perts, genes, env unchanged).
# This is exactly the construction 707's A/B/C board read tested (KA: on A/B/C's dump it reproduces 77e9de2d byte for byte,
# work/ka_hook707.log). It needs only the donor banks on disk and the bundle's pert_counts.csv (no network, no identity).
# Ships for October only if 707 reads TRANSFERS (>= +0.0010) AND readiness rehearses it before the freeze.
set -euo pipefail
B="$1"; W="$2"; DUMP="$3"; ID="$4"; HD="$5"
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
mkdir -p "$HD"
OUT="$HD/$(basename "$DUMP")"
echo "dumphook707: bundle $B; work $W; identity $ID (unused); in $DUMP ($(shasum -a 256 "$DUMP" | cut -c1-16))"
uv run python experiments/707-strong-tier-direction/swap707.py --dump-in "$DUMP" --dump-out "$OUT" --bundle "$B" --beta 0.25
echo "dumphook707: wrote $OUT ($(shasum -a 256 "$OUT" | cut -c1-16))"
