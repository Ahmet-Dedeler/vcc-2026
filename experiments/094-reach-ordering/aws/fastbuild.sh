#!/usr/bin/env bash
# Experiment 094 fast build path.
#
# The exp089 ceremony runner re-downloads the donor tree, the controls bundle and the
# whole Python environment on every launch: ~3h20m wall for a build whose own
# `/usr/bin/time -v` reads 16:55. Everything it fetches is already on this worker and
# hash-pinned, and the arms below differ ONLY in emission flags, so this path reuses the
# extracted tree and re-runs the two steps that actually depend on the change.
#
# Provenance is preserved where it matters: the same pinned whitener sha256 is asserted
# by make_submission_pm_cd4.py itself, the test suite runs before every build, and each
# arm writes its own build-report.json carrying the full flag set.
set -Eeuo pipefail

ARM="${ARM:?ARM is required}"
# Per-arm repo. Arms run concurrently and each refreshes its own sources from S3, so a
# shared tree would let one arm swap the script another is mid-import on. The base tree is
# copied once per arm; it is ~200 MB of source next to a 20 GB h5ad, so the copy is free.
BASE_REPO="${BASE_REPO:-/work/exp089-20260826T002404Z/repo}"
REPO="${REPO:-/work/exp094/repo-${ARM}}"
INPUT="${INPUT:-/work/exp089-20260826T002404Z/input}"
EXP061_REPO="${EXP061_REPO:-/work/exp061/repo}"
OUT="/work/exp094/${ARM}"
MIN_CPM="${MIN_CPM:-10.0}"
# The operator. Default is v11's shipped union whitener. WHITENER=pruned selects the
# experiment-091 basis, which drops the 4,393 common-essential rows -- exp091 measured
# +0.002409 avg6 (t=19.85) against a SIZE-MATCHED random prune, so the effect is
# composition and not basis size, and its metadata records dropped_in_panel: [].
WHITENER="${WHITENER:-shipped}"
if [ "$WHITENER" = "pruned" ]; then
  WHITENER_ID="union_pruned"
  WHITENER_SHA256="e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
  WHITENER_NPZ="${INPUT}/whitener_union_pruned.npz"
  WHITENER_JSON="${INPUT}/whitener_union_pruned.json"
else
  WHITENER_ID="union"
  WHITENER_SHA256="5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3"
  WHITENER_NPZ="${INPUT}/whitener_union.npz"
  WHITENER_JSON="${INPUT}/whitener_union.json"
fi

mkdir -p "$OUT"
if [ ! -d "$REPO" ]; then
  mkdir -p "$(dirname "$REPO")"
  cp -al "$BASE_REPO" "$REPO" 2>/dev/null || cp -a "$BASE_REPO" "$REPO"
  # break the hard links for the two files each arm rewrites
  for f in scripts/make_submission_pm_cd4.py scripts/test_make_submission_pm_cd4.py; do
    cp --remove-destination "$(readlink -f "${REPO}/${f}")" "${REPO}/${f}" 2>/dev/null || true
  done
fi
exec > >(tee -a "${OUT}/cloud.log") 2>&1
echo "=== arm ${ARM} starting $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "flags: VCC_ALPHA=${VCC_ALPHA:-<module default>} MIN_CPM=${MIN_CPM} EXTRA=${EXTRA_ARGS:-<none>}"

# Per-arm lock. The exp089 runner serialised every heavy VCC job behind one global lock,
# which was right when a build re-downloaded 100+ GB and peaked at 110 GiB in `vcc prep`.
# This worker has 32 vCPU and 247 GiB, the build itself is single-threaded and peaks at
# ~15 GiB, and `vcc prep` peaks at ~74 GiB for well under a minute. Two arms fit; the lock
# is kept per-arm only to stop the same arm being launched twice.
exec 9>"/work/vcc-${ARM}.lock"
flock -n 9 || { echo "arm ${ARM} is already running" >&2; exit 70; }

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR=/work/uv-cache
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd "$REPO"

if [ "$WHITENER" = "pruned" ] && [ ! -s "$WHITENER_NPZ" ]; then
  echo "--- fetching the pruned operator ---"
  aws s3 cp "s3://${BUCKET}/input/exp094/whitener_union_pruned.npz" "$WHITENER_NPZ" --only-show-errors
  aws s3 cp "s3://${BUCKET}/input/exp094/whitener_union_pruned.json" "$WHITENER_JSON" --only-show-errors
fi
printf '%s  %s\n' "$WHITENER_SHA256" "$WHITENER_NPZ" | sha256sum --check -

echo "--- refreshing the changed sources from S3 ---"
# make_submission_pm_cd4.py verifies FROZEN_IMPLEMENTATION_SHA256 over its helper modules,
# so a build fails closed unless every changed module travels with it. multi_donor.py is
# refreshed here because donor_agreement() lives in it; it is optional so that an older
# SRC_ID without that file still builds.
aws s3 cp "s3://${BUCKET}/input/exp094/${SRC_ID}/make_submission_pm_cd4.py" \
  scripts/make_submission_pm_cd4.py --only-show-errors
aws s3 cp "s3://${BUCKET}/input/exp094/${SRC_ID}/test_make_submission_pm_cd4.py" \
  scripts/test_make_submission_pm_cd4.py --only-show-errors
for m in multi_donor donor_transfer; do
  aws s3 cp "s3://${BUCKET}/input/exp094/${SRC_ID}/${m}.py" \
    "src/vcc2026/${m}.py" --only-show-errors 2>/dev/null \
    || echo "no ${m}.py at this SRC_ID; keeping the tree's copy"
done
sha256sum scripts/make_submission_pm_cd4.py scripts/test_make_submission_pm_cd4.py \
  src/vcc2026/multi_donor.py src/vcc2026/donor_transfer.py

echo "--- test suite (fail-closed before any build) ---"
uv run python scripts/test_make_submission_pm_cd4.py 2>&1 | tee "${OUT}/tests.log"

H5AD="${OUT}/${ARM}.h5ad"
VCC="${OUT}/${ARM}.vcc"
rm -f "$H5AD" "$VCC"

echo "--- build ---"
# shellcheck disable=SC2086
/usr/bin/time -v uv run python scripts/make_submission_pm_cd4.py \
  --bundle data/challenge \
  --k562 "${EXP061_REPO}/data/donors/k562_gwps" \
  --hct116 "${INPUT}/orion_hct116" \
  --hek293t "${EXP061_REPO}/data/donors/orion_hek293t" \
  --perturbmatch "${EXP061_REPO}/data/external/perturbmatch" \
  --cd4-root "${INPUT}/cd4/cd4_compact" \
  --whitener "$WHITENER_NPZ" \
  --whitener-metadata "$WHITENER_JSON" \
  --expected-whitener-sha256 "$WHITENER_SHA256" \
  --whitener-id "$WHITENER_ID" \
  --output "$H5AD" \
  --min-cpm "$MIN_CPM" \
  --max-rss-gib 120 \
  ${EXTRA_ARGS:-} \
  2>&1 | tee "${OUT}/build.log"

sha256sum "$H5AD" > "${H5AD}.sha256"

echo "--- vcc prep ---"
/usr/bin/time -v vcc prep "$H5AD" \
  -g data/challenge/gene_names.csv \
  --perts data/challenge/pert_counts.csv \
  -o "$VCC" 2>&1 | tee "${OUT}/vcc-prep.log"

sha256sum "$VCC" > "${VCC}.sha256"
ls -la "$OUT"

echo "--- uploading package ---"
aws s3 cp "$VCC" "s3://${BUCKET}/output/exp094/${ARM}/$(basename "$VCC")" --only-show-errors
aws s3 cp "${VCC}.sha256" "s3://${BUCKET}/output/exp094/${ARM}/$(basename "$VCC").sha256" --only-show-errors
# The meta.json carries the full audited config -- including experiment 099's measured `f`
# under effects.cd4.residual_frame -- and was previously left on the worker, where it died
# with the instance. It is 17 KB next to a 3.7 GB package; upload it.
for f in cloud.log build.log tests.log vcc-prep.log "${ARM}.h5ad.meta.json"; do
  aws s3 cp "${OUT}/${f}" "s3://${BUCKET}/output/exp094/${ARM}/${f}" --only-show-errors || true
done
echo "=== arm ${ARM} DONE $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
