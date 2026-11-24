#!/usr/bin/env bash
# Fetch the 2025 Virtual Cell Challenge H1 hESC dataset (~34 GB).
# Same 10x Flex chemistry and knockdown criteria as 2026, and it includes the
# previously held-out test split — this is the local ground-truth harness.
# Public bucket, no auth required.
set -euo pipefail

BASE="https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025"
DEST="${1:-data/external/vcc2025}"
mkdir -p "$DEST"/{train,validation,test}

fetch() {  # fetch <remote-path> <local-path>
  if [[ -f "$DEST/$2" ]]; then echo "have  $2"; return; fi
  echo "get   $2"
  curl -fL --retry 5 --retry-delay 5 -C - -o "$DEST/$2.part" "$BASE/$1"
  mv "$DEST/$2.part" "$DEST/$2"
}

fetch "gene_names.csv"                        "gene_names.csv"
fetch "train/pert_counts_Training.csv"        "train/pert_counts_Training.csv"
fetch "validation/pert_counts_Validation.csv" "validation/pert_counts_Validation.csv"
fetch "test/pert_counts_Test.csv"             "test/pert_counts_Test.csv"
fetch "validation/adata_Validation.h5ad"      "validation/adata_Validation.h5ad"   #  6.9 GB
fetch "test/adata_Test.h5ad"                  "test/adata_Test.h5ad"               # 12.0 GB
fetch "train/adata_Training.h5ad"             "train/adata_Training.h5ad"          # 15.5 GB

echo
du -sh "$DEST"
