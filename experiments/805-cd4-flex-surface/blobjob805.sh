#!/usr/bin/env bash
# 805 blob-direction diagnostic on EC2 (labworker runs it from /work/repo). Args: the s0 run id (record only; the truths come from
# experiments/805-cd4-flex-surface/work/truths/<cond>/ in the input mirror),
# then the conditions. KA first: blobdir805 on h1-test must reproduce 615's split-half 0.04869 (|d| <= 1e-4), and on the
# thinned HEK truth cos(mu, hekblob) must be 1. Writes experiments/805-cd4-flex-surface/work/blobdir/*.json.
set -Eeuo pipefail
RID="$1"; shift
X=experiments/805-cd4-flex-surface; O="$X/work/blobdir"; mkdir -p "$O" /work/truths
echo "e37f39bb6aae9ce6334f133416b4780ac4226b5a407c4f1a511a5a09ef75d6be  data/harness/orion-hek293t-thinned.h5ad" | sha256sum --check -
echo "06c6baea8c62ef3d38cacdbf58e104908eb93f0ef3b25aca56b86d28240f9692  data/harness/h1-test.h5ad" | sha256sum --check -
uv run python $X/blobdir805.py data/harness/h1-test.h5ad "$O/blobdir_h1test.json"
uv run python -c "import json,sys; d=json.load(open('$O/blobdir_h1test.json')); sys.exit(0 if abs(d['share_sh']-0.04869282269403374)<=1e-4 else 1)" || { echo "KA FAILED (h1-test split-half)"; exit 9; }
uv run python $X/blobdir805.py data/harness/orion-hek293t-thinned.h5ad "$O/blobdir_hekthinned.json"
uv run python -c "import json,sys; d=json.load(open('$O/blobdir_hekthinned.json')); sys.exit(0 if d['cos_mu_hekblob']>0.999 else 1)" || { echo "KA FAILED (hek self-cosine)"; exit 9; }
for c in "$@"; do
  # the instance role cannot read output/: the truths were copied (S3->S3) to the input mirror and synced by labworker
  mv "$X/work/truths/${c}/truth_P.h5ad" "/work/truths/${c}_truth_P.h5ad"
  cp "$X/work/truths/${c}/select.json" "$O/select_${c}.json"
  want="$(jq -r '.truth_P.sha256' "$O/select_${c}.json")"
  echo "${want}  /work/truths/${c}_truth_P.h5ad" | sha256sum --check -
  uv run python $X/blobdir805.py "/work/truths/${c}_truth_P.h5ad" "$O/blobdir_cd4_${c}.json"
done
