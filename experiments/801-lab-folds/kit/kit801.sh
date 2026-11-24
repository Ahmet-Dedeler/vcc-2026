#!/usr/bin/env bash
# kit801.sh -- the lab's one-command PANEL KIT (follow-up 7, for 10-22). Builds, on ONE EC2 instance, the three cancer held-out
# folds of any 300-target panel by the L-fold recipe, unchanged:
#   CD4 banks (cd4_801.py), a FRESH day-recipe held-out operator per held line (build_union_whitener.py --basis-prune essential,
#   12k cap; fitted at 8 threads, kept at s3://<bucket>/output/lab/ops801/<TAG>/), dump801 folds <TAG>/<held> (D/E/F rules;
#   frames = A/B/C's relabelled controls, as every lab fold), ka804's L arm set (unmodified) -> work/ka-<tag>/{bulk,meta},
#   a registry fragment work/registry-<TAG>.json; with --ref <REF>, a per-array KA of <TAG> against <REF>'s dumps and ka804 arrays.
# The Mac side is shell only (cp, cmp, aws, the launcher): it runs under the lab's local-job hold. It READS the panel files and
# never writes beside them (on 10-22 they are data/challenge-def/*, which production owns).
#
#   bash experiments/801-lab-folds/kit/kit801.sh <TAG> <pert_counts.csv> <gene_names.csv> [--ref REF] [--frames DIR] [--no-wait]
#   --frames DIR (repo-relative bundle dir, e.g. data/challenge-def): the folds' control frames are DIR's own context files
#   (manifest.json's contexts; sha256 recorded on the Mac, re-checked on EC2); default: A/B/C's controls relabelled D/E/F.
#   e.g. 10-22:  bash experiments/801-lab-folds/kit/kit801.sh DEF data/challenge-def/pert_counts.csv data/challenge-def/gene_names.csv
#   dry run:     bash experiments/801-lab-folds/kit/kit801.sh LABCk data/challenge/pert_counts.csv data/challenge/gene_names.csv --ref LABC
# Cost and time per panel (dry run, 10-16): see RESULTS.md section 17. After it returns (or `labcloud.sh wait <run>` with
# --no-wait), the outputs are under work/ and mirrored to s3://<bucket>/input/lab/repo/experiments/801-lab-folds/work/.
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$ROOT"
E=experiments/801-lab-folds; W=$E/work; K=$E/kit
LC=experiments/800-checker/cloud/labcloud.sh
BUCKET="${BUCKET:-<your-s3-bucket>}"
TAG="${1:?TAG}"; PC="${2:?pert_counts.csv}"; GN="${3:?gene_names.csv}"; shift 3
REF=""; WAIT=1; FR=""
while [[ $# -gt 0 ]]; do case "$1" in --ref) REF="$2"; shift 2;; --no-wait) WAIT=0; shift;; --frames) FR="${2%/}"; shift 2;;
  *) echo "unknown $1"; exit 2;; esac; done
die() { echo "kit801: $*" >&2; exit 2; }
[[ "$TAG" =~ ^[A-Za-z][A-Za-z0-9]{0,15}$ ]] || die "TAG must be [A-Za-z][A-Za-z0-9]{0,15}"
case "$TAG" in L[1-8]|A[1-6]|LABC|LABCd|ABC|ABCp|ABC062|efp|efr|real|sealed) die "TAG $TAG is taken";; esac
[[ ! -e $W/$TAG ]] || die "$W/$TAG exists: one kit run per tag"
tg="$(echo "$TAG" | tr 'A-Z' 'a-z')"
# ---- the panel files: read, checked, copied (never modified in place)
head -1 "$PC" | tr -d '\r' | tr ',' '\n' | grep -qx target_gene || die "$PC has no target_gene column"
[[ "$(tail -n +2 "$PC" | grep -c .)" -eq 300 ]] || die "$PC does not list 300 targets"
cmp -s "$GN" data/challenge/gene_names.csv || die "$GN differs from data/challenge/gene_names.csv (the engine's official axis); stop and ask"
mkdir -p $W/$TAG/bundle
cp "$PC" $W/$TAG/bundle/pert_counts.csv
cp "$GN" $W/$TAG/bundle/gene_names.csv
echo "kit801: $TAG <- $PC ($(shasum -a 256 "$PC" | cut -c1-16)), $GN ($(shasum -a 256 "$GN" | cut -c1-16))"
if [[ -n "$FR" ]]; then                 # --frames <bundle dir>: its manifest's contexts are the control frames (repo-relative dir)
  [[ "$FR" != /* && -f "$FR/manifest.json" ]] || die "--frames needs a repo-relative bundle dir with manifest.json"
  CTX="$(jq -r '.contexts[]' "$FR/manifest.json" | tr '\n' ' ')"
  [[ "$(echo $CTX | wc -w)" -eq 3 ]] || die "$FR/manifest.json does not list 3 contexts"
  echo "$FR" > $W/$TAG/bundle/frames_src.txt
  ( cd "$FR" && shasum -a 256 manifest.json $(for c in $CTX; do printf 'context_%s.h5ad ' $c; done) ) | sed 's/ \*/  /' > $W/$TAG/bundle/frames.sha256
  for f in manifest.json $(for c in $CTX; do printf 'context_%s.h5ad ' $c; done); do   # the frames' files in the mirror
    loc=$(stat -f %z "$FR/$f"); rem=$(aws s3 ls "s3://${BUCKET}/input/lab/repo/$FR/$f" 2>/dev/null | awk '{print $3}')
    [[ "$loc" == "$rem" ]] || bash $LC stage "$FR/$f" >/dev/null       # (EC2's prep re-checks every sha256)
  done
  echo "kit801: frames = $FR's own controls ($CTX)"
fi
bash $LC stage $W/$TAG/bundle >/dev/null
# ---- jobs (cloud801.sh subcommands; every one already validated on L1-L6 / A1-A6 / LABC)
J=$W/cloud/jobs-kit-$TAG.txt
{ echo "bash $E/cloud801.sh verify"
  echo "bash $E/cloud801.sh prep $TAG"
  for h in k562_gwps orion_hct116 orion_hek293t; do echo "FIT_THREADS=8 T801_THREADS=3 bash $E/cloud801.sh fold $TAG $h"; done
  for h in k562_gwps orion_hct116 orion_hek293t; do echo "KA_DIR=$tg T801_THREADS=2 bash $E/cloud801.sh ka804 $TAG $h"; done
  echo "bash $E/cloud801.sh regfrag $TAG"
  if [[ -n "$REF" ]]; then echo "bash $E/cloud801.sh kacmp $TAG $REF"; fi
} > $J
SY="$(cat $K/sync-base.txt) --sync $W/$TAG/bundle"
if [[ -n "$FR" ]]; then SY="$SY --sync $FR/manifest.json"; for c in $CTX; do SY="$SY --sync $FR/context_$c.h5ad"; done; fi
if [[ -n "$REF" ]]; then SY="$SY --sync $W/ka-$(echo "$REF" | tr 'A-Z' 'a-z')"; fi
OUTL="$(bash $LC launch --name "lab801-kit-$tg" --jobs $J --type r7i.8xlarge --minutes 60 --parallel 9 --root-gb 300 $SY \
  --code-extra $W/manifest --collect "$W/dumps/dump708-$TAG-*" --collect "$W/ka-$tg" --collect $W/cloud \
  --collect "$W/$TAG/ops/held_*.json" --collect "$W/$TAG/cd4_801.json" --collect "$W/registry-$TAG.json" 2>&1)" || { echo "$OUTL"; exit 3; }
echo "$OUTL" | tail -4
RID="$(echo "$OUTL" | sed -n 's/^run_id=//p')"
[[ -n "$RID" ]] || die "no run id"
echo "$TAG $RID $(date -u +%FT%TZ)" >> $K/kit-runs.txt
[[ $WAIT -eq 1 ]] || { echo "kit801: launched $RID; later: bash $K/kit801-collect.sh $TAG $RID"; exit 0; }
bash $K/kit801-collect.sh "$TAG" "$RID"
