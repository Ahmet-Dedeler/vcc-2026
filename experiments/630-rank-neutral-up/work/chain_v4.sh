#!/bin/bash
# 630 v4 capacity screen: b_k100_q20r3p10 with each capacity flag alone (a, z, r) and all three, on H1 seeds 0-2.
# KA: the flags-off variant (_f) must reproduce b_k100_q20r3p10's stats bit for bit on every seed.
set -u
R=~/Code/virtual-cell-challenge-2026
W=$R/experiments/630-rank-neutral-up/work
cd $R
V=b_k100_q20r3p10
for s in 0 1 2; do
  /usr/bin/time -l uv run python experiments/630-rank-neutral-up/rnup630.py $s ${V}_f ${V}_fa ${V}_fz ${V}_fr ${V}_fazr > $W/rnup-v4-s$s.log 2>&1
  rc=$?
  ka=$(uv run python -c "
import numpy as np
a=np.load('$W/h1/stats-rnup_s${s}_${V}.npz'); b=np.load('$W/h1/stats-rnup_s${s}_${V}_f.npz')
print('KA PASS' if all(np.array_equal(a[k],b[k],equal_nan=True) for k in ('sums','jk','p','q','lfc')) else 'KA FAIL')" 2>&1 | tail -1)
  echo "$(date -u +%FT%TZ) v4 s$s rc $rc $ka"
done
uv run python experiments/630-rank-neutral-up/price630.py ${V} ${V}_fa ${V}_fz ${V}_fr ${V}_fazr > $W/price630-v4.log 2>&1
echo "$(date -u +%FT%TZ) price rc $?"
echo "$(date -u +%FT%TZ) END"
