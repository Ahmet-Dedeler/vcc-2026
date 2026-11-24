#!/usr/bin/env bash
# make-jobs811c.sh <TAG> -- write the job and sync files for run-kit811c.sh on a kit panel <TAG> (after kit801.sh + kit801-collect.sh
# built <TAG>'s folds). A copy of the 811 half of 801's kit/make-jobs.sh, with 811's composite inputs added. Launch (the lead):
#   bash experiments/800-checker/cloud/labcloud.sh launch --name lab811c-kit-<tag> --jobs experiments/811-lolo-model/kit/jobs811c-<TAG>.txt \
#     --type r7i.4xlarge --minutes 45 --parallel 1 $(sed 's/^/--sync /' experiments/811-lolo-model/kit/sync811c-<TAG>.txt) \
#     --collect experiments/811-lolo-model/work/run-comp811kit-<TAG>
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"; cd "$ROOT"
TAG="${1:?TAG}"; K=experiments/811-lolo-model/kit; W=experiments/801-lab-folds/work
[[ -f $W/registry-$TAG.json ]] || { echo "no $W/registry-$TAG.json: run kit801 first"; exit 2; }
echo "bash $K/run-kit811c.sh $TAG" > $K/jobs811c-$TAG.txt
{ grep . experiments/811-lolo-model/cloud/sync-aonly.txt | grep -v -E "work/A[1-6]/"
  printf '%s\n' $W/$TAG/bundle $W/$TAG/ops $W/registry-$TAG.json $W/dumps \
    experiments/811-lolo-model/work/run-m811s1 experiments/811-lolo-model/work/refcomp
  TG=$(echo "$TAG" | tr "A-Z" "a-z"); [[ -d $W/ka-$TG ]] && echo $W/ka-$TG          # kit801 names ka dirs by the lower-case tag
} | sort -u > $K/sync811c-$TAG.txt
echo "wrote $K/jobs811c-$TAG.txt $K/sync811c-$TAG.txt ($(wc -l < $K/sync811c-$TAG.txt) sync paths)"
