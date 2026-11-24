#!/usr/bin/env bash
# lzaccess811.sh — the ONE override for the sealed reserve LZ, keyed to LZ-CLAIM-811.md's sha256 (claim §"What is priced against
# what": "The kit refuses L7/L8 by design. The access adds ONE override keyed to this file's sha256 ... and nothing else changes").
#
# No guard, engine or kit file is edited. This wrapper checks, in order, and stops on the first failure:
#   1. LZ-CLAIM-811.md's sha256 == 3a8f7252…c258c2 and docs/CONFIRMATION-LEDGER.md records that sha (the coordinator's line);
#   2. the sealed target files' sha256 == SEALED.md's (L7 33bb3256…, L8 c8d9f1a3…);
#   3. no earlier access: ACCESS-LZ-811.md absent, and no work/LZ7, work/LZ8 kit dirs or ka-lz7/ka-lz8 in 801's work.
# Then it writes the access record's first line and runs 801's kit, UNCHANGED, on the sealed target lists under the tags LZ7 (= L7)
# and LZ8 (= L8), with A/B/C frames (the L-fold recipe; no --frames). The fresh tags are the override: the kit and engine guards
# refuse the literal names L7/L8, and this script is the only place the sealed lists are handed to the kit.
#   bash experiments/800-checker/lzaccess811.sh
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$ROOT"
H=experiments/800-checker; W=experiments/801-lab-folds/work
CLAIM=$H/LZ-CLAIM-811.md; WANT=3a8f725216075e793a9968b53f2ff0f28f845593ba49e48d3c659f009bd258c2
REC=$H/ACCESS-LZ-811.md
die() { echo "lzaccess811: STOP: $*" >&2; exit 2; }
got=$(shasum -a 256 "$CLAIM" | cut -d' ' -f1); [[ "$got" == "$WANT" ]] || die "claim sha $got != $WANT"
/usr/bin/grep -qF "$WANT" docs/CONFIRMATION-LEDGER.md || die "the ledger does not record the claim sha"
s7=$(shasum -a 256 $W/sealed/L7/pert_counts.csv | cut -d' ' -f1); s8=$(shasum -a 256 $W/sealed/L8/pert_counts.csv | cut -d' ' -f1)
[[ "$s7" == 33bb3256ce92f058d7d71c3b7b5e1a3bd15799d6b1e12d4c2af5e0341d8723dc ]] || die "L7 file sha $s7"
[[ "$s8" == c8d9f1a335a66e9ec9252bd0d5983bc104525b88b128be758d3c0bef77e3ab65 ]] || die "L8 file sha $s8"
[[ ! -e $REC ]] || die "$REC exists: LZ was already accessed (one access in total)"
for p in $W/LZ7 $W/LZ8 $W/ka-lz7 $W/ka-lz8; do [[ ! -e $p ]] || die "$p exists"; done
{ echo "# LZ access record (one access; claim $CLAIM sha256 $WANT)"; echo
  echo "- $(date -u +%FT%TZ) access opened by the lab lead. Override: this wrapper (\`lzaccess811.sh\` sha256 $(shasum -a 256 "$0" | cut -d' ' -f1));"
  echo "  tags LZ7 = L7 (file $s7), LZ8 = L8 (file $s8); A/B/C frames; kit801 $(shasum -a 256 experiments/801-lab-folds/kit/kit801.sh | cut -c1-16) unchanged."
} > $REC
for t in 7 8; do
  bash experiments/801-lab-folds/kit/kit801.sh LZ$t $W/sealed/L$t/pert_counts.csv data/challenge/gene_names.csv --no-wait 2>&1 | tail -3
done
echo "- $(date -u +%FT%TZ) kits launched: $(grep -E '^LZ[78] ' experiments/801-lab-folds/kit/kit-runs.txt | awk '{print $1"="$2}' | tr '\n' ' ')" >> $REC
cat $REC
