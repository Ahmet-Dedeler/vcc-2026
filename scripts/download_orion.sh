#!/usr/bin/env bash
# X-Atlas/Orion (Xaira): genome-wide FiCS Perturb-seq in HCT116 (colorectal) and
# HEK293T (kidney epithelial). 8M cells, all 18,903 protein-coding genes targeted.
#
# Why this matters: Replogle K562-gwps is our only donor covering the VCC panel
# (272/300), so every prediction leans on one cell line. Orion adds two more lines
# AND covers 25 of the 28 targets gwps misses. It takes us from 1 donor to 3.
#
# https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion  (CC-BY-NC-SA-4.0)
set -uo pipefail
LINE="${1:-HCT116}"
DEST="data/external/orion"
mkdir -p "$DEST"
hf download Xaira-Therapeutics/X-Atlas-Orion \
  --repo-type dataset --local-dir "$DEST" \
  --include "data/${LINE}_*.parquet" --include "metadata/*" \
  2>&1 | tail -3
du -sh "$DEST"
