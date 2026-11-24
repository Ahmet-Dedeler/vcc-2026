#!/usr/bin/env bash
# make-jobs.sh <TAG> -- after kit801.sh (+ kit801-collect.sh) built <TAG>'s folds: write the 810 / 811 job and sync files that
# read them, from 810's / 811's own A-only lists with the A panels swapped for <TAG> (bundle, ops mirror, registry fragment).
# Launch (the lead, through 810 / 811):
#   bash experiments/800-checker/cloud/labcloud.sh launch --name lab810-kit-<tag> --jobs experiments/801-lab-folds/kit/jobs810-<TAG>.txt \
#     --type r7i.4xlarge --minutes 60 --parallel 4 $(sed 's/^/--sync /' experiments/801-lab-folds/kit/sync810-<TAG>.txt) \
#     --collect experiments/810-strong-direction/work/kit-<TAG> --code-extra experiments/801-lab-folds/work/manifest
#   bash experiments/800-checker/cloud/labcloud.sh launch --name lab811-kit-<tag> --jobs experiments/801-lab-folds/kit/jobs811-<TAG>.txt \
#     --type r7i.4xlarge --minutes 60 --parallel 1 $(sed 's/^/--sync /' experiments/801-lab-folds/kit/sync811-<TAG>.txt) \
#     --collect experiments/811-lolo-model/work/run-m811kit-<TAG>
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"; cd "$ROOT"
TAG="${1:?TAG}"; K=experiments/801-lab-folds/kit; W=experiments/801-lab-folds/work
[[ -f $W/registry-$TAG.json ]] || { echo "no $W/registry-$TAG.json: run kit801 first"; exit 2; }
{ echo "bash $K/kit810.sh prep $TAG"; for h in k562_gwps orion_hct116 orion_hek293t; do echo "bash $K/kit810.sh fold $TAG $h"; done; } > $K/jobs810-$TAG.txt
echo "bash $K/run-kit811.sh $TAG" > $K/jobs811-$TAG.txt
swap() {  # their sync files are one repo path per line: drop the A panels' bundles / 810's ops copies, add <TAG>'s inputs
  grep . "$1" | grep -v -E "work/A[1-6]/|810-strong-direction/work/aops"
  printf '%s\n' $W/$TAG/bundle $W/$TAG/ops $W/registry-$TAG.json $W/dumps
}
swap experiments/810-strong-direction/cloud/sync810aonly.txt > $K/sync810-$TAG.txt
swap experiments/811-lolo-model/cloud/sync-aonly.txt > $K/sync811-$TAG.txt
echo "wrote $K/jobs810-$TAG.txt $K/sync810-$TAG.txt $K/jobs811-$TAG.txt $K/sync811-$TAG.txt"
