#!/usr/bin/env bash
# kit801-collect.sh <TAG> <RUN_ID> -- the kit's second half (Mac, shell only): labcloud wait (collect, terminate, cost), copy the
# outputs into work/ (APFS clones), mirror them server-side to the S3 input mirror, and merge the registry fragment into
# work/registry.json (a sub-second json merge; skipped and left as work/registry-<TAG>.json while window800 prints WINDOW).
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$ROOT"
E=experiments/801-lab-folds; W=$E/work
LC=experiments/800-checker/cloud/labcloud.sh
BUCKET="${BUCKET:-<your-s3-bucket>}"
TAG="${1:?TAG}"; RID="${2:?run id}"; tg="$(echo "$TAG" | tr 'A-Z' 'a-z')"
bash $LC wait "$RID" | tail -3
D=experiments/800-checker/cloud/runs/$RID/out/repo/$W
[[ -d $D ]] || { echo "kit801: nothing collected for $RID"; exit 3; }
mkdir -p $W/dumps $W/ka-$tg $W/cloud/kit-$TAG $W/$TAG/ops
for f in $D/dumps/dump708-$TAG-*; do cp -c "$f" $W/dumps/; done
cp -R $D/ka-$tg/. $W/ka-$tg/
cp -R $D/cloud/. $W/cloud/kit-$TAG/
cp $D/$TAG/ops/held_*.json $W/$TAG/ops/
[[ -f $D/$TAG/cd4_801.json ]] && cp $D/$TAG/cd4_801.json $W/$TAG/
cp $D/registry-$TAG.json $W/
S="s3://${BUCKET}/output/lab/${RID}/repo/$W"; M="s3://${BUCKET}/input/lab/repo/$W"
aws s3 sync "$S/ka-$tg/" "$M/ka-$tg/" --only-show-errors
for f in $D/dumps/dump708-$TAG-*; do aws s3 cp "$S/dumps/$(basename "$f")" "$M/dumps/$(basename "$f")" --only-show-errors; done
aws s3 cp "$S/registry-$TAG.json" "$M/registry-$TAG.json" --only-show-errors
# the operators, server-side, into the mirror (810's / 811's kit jobs sync them from there; 3 x ~1.3 GB, S3-internal)
aws s3 sync "s3://${BUCKET}/output/lab/ops801/$TAG/" "$M/$TAG/ops/" --exclude probe.txt --only-show-errors
echo "kit801: $TAG outputs in $W (dumps, ka-$tg, registry-$TAG.json) and mirrored to $M"
grep -h "cloud801" experiments/800-checker/cloud/runs/$RID/out/jobs/*.log | sort | cut -c1-200 | tail -20
if bash experiments/800-checker/window800.sh >/dev/null; then
  python3 - "$W/registry.json" "$W/registry-$TAG.json" <<'PY'
import json, sys
reg, frag = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
clash = [k for k in frag if k in reg["folds"]]
assert not clash, f"registry already has {clash}"
reg["folds"].update(frag)
open(sys.argv[1], "w").write(json.dumps(reg, indent=1))
print(f"registry: +{len(frag)} folds -> {len(reg['folds'])}")
PY
  bash $LC stage $W/registry.json >/dev/null
else
  echo "kit801: window open; registry fragment left at $W/registry-$TAG.json (merge later)"
fi
