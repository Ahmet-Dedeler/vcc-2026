#!/bin/zsh
# Parallel range download of GWCD4i.pseudobulk_merged.h5ad (44,566,657,140 B) into work/.
set -e
URL="https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.pseudobulk_merged.h5ad"
OUT="pseudobulk_merged.h5ad"
SIZE=44566657140
N=8
CHUNK=$(( (SIZE + N - 1) / N ))
mkdir -p parts
for i in $(seq 0 $((N-1))); do
  LO=$(( i * CHUNK ))
  HI=$(( LO + CHUNK - 1 )); [[ $HI -ge $SIZE ]] && HI=$(( SIZE - 1 ))
  curl -sS --fail --retry 5 --retry-delay 2 -r ${LO}-${HI} "$URL" -o parts/p$i &
done
wait
cat parts/p0 parts/p1 parts/p2 parts/p3 parts/p4 parts/p5 parts/p6 parts/p7 > "$OUT"
rm -rf parts
ls -l "$OUT"
