#!/usr/bin/env bash
# Public cross-line CRISPR Perturb-seq datasets. These are the only material from
# which cross-context transfer can be learned. Ordered smallest-first.
# Written for bash 3.2 (macOS default) — no associative arrays.
#
#   ./scripts/download_external.sh replogle_rpe1 replogle_k562
#   ./scripts/download_external.sh all
set -uo pipefail
DEST="data/external"
ORDER="replogle_rpe1 replogle_k562 nadig_hepg2 nadig_jurkat replogle_k562_gwps"

meta() {  # meta <key> -> "<url>\t<relpath>\t<size>"
  case "$1" in
    replogle_rpe1)      printf '%s\t%s\t%s' "https://plus.figshare.com/ndownloader/files/35775606" "replogle2022/RPE1_essential.h5ad" "8.1 GB" ;;
    replogle_k562)      printf '%s\t%s\t%s' "https://plus.figshare.com/ndownloader/files/35773219" "replogle2022/K562_essential.h5ad" "9.9 GB" ;;
    replogle_k562_gwps) printf '%s\t%s\t%s' "https://plus.figshare.com/ndownloader/files/35775507" "replogle2022/K562_gwps.h5ad" "61.3 GB" ;;
    nadig_hepg2)        printf '%s\t%s\t%s' "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE264667&format=file&file=GSE264667%5Fhepg2%5Fraw%5Fsinglecell%5F01%2Eh5ad" "nadig2025/hepg2_raw_singlecell.h5ad" "5.2 GB" ;;
    nadig_jurkat)       printf '%s\t%s\t%s' "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE264667&format=file&file=GSE264667%5Fjurkat%5Fraw%5Fsinglecell%5F01%2Eh5ad" "nadig2025/jurkat_raw_singlecell.h5ad" "8.7 GB" ;;
    *) return 1 ;;
  esac
}

[ "${1:-}" = "all" ] && set -- $ORDER
[ $# -eq 0 ] && { echo "usage: $0 [all|$ORDER]"; exit 1; }

for key in "$@"; do
  info=$(meta "$key") || { echo "SKIP  unknown key: $key"; continue; }
  url=$(printf '%s' "$info" | cut -f1)
  rel=$(printf '%s' "$info" | cut -f2)
  size=$(printf '%s' "$info" | cut -f3)
  out="$DEST/$rel"
  if [ -f "$out" ]; then echo "have  $key"; continue; fi
  mkdir -p "$(dirname "$out")"
  echo "get   $key ($size) -> $out"
  if curl -fL --retry 5 --retry-delay 5 -C - -s -S -o "$out.part" "$url"; then
    mv "$out.part" "$out"; echo "done  $key ($(du -h "$out" | cut -f1))"
  else
    echo "FAIL  $key — leaving .part for resume"
  fi
done
du -sh "$DEST"/*/ 2>/dev/null
