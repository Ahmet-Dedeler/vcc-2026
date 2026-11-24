#!/usr/bin/env bash
set -Eeuo pipefail

# Resume only the packaging stage of the preserved experiment-071 H5AD.
# This script never rebuilds predictions and never invokes `vcc submit`.

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing resume: set I_HAVE_EXPLICIT_CLEARANCE=YES" >&2
  exit 64
fi

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:-20260823T083438Z}"
EXP065_RUN_ID="${EXP065_RUN_ID:-20260823T040540Z}"
EXPECTED_H5AD_SHA256="${EXPECTED_H5AD_SHA256:-e6a539faa2fe0cc1169c7a6083813c4aaa299d2c609f302749bf8809b071f809}"
EXPECTED_H5AD_BYTES="${EXPECTED_H5AD_BYTES:-20358005624}"
PREP_MAX_RSS_GIB="${PREP_MAX_RSS_GIB:-110}"
WORK="/work/exp071-${RUN_ID}"
REPO="${WORK}/repo"
OUT="${WORK}/output"
H5AD="${OUT}/v10-pm-cd4.h5ad"
VCC="${OUT}/v10-pm-cd4.vcc"
PARTIAL="${OUT}/v10-pm-cd4.partial-${RUN_ID}.vcc"
DRY_RUN_VCC="${OUT}/.v10-pm-cd4.dry-run-${RUN_ID}.vcc"
LOG="${OUT}/resume-prep.log"
EXIT_FILE="${OUT}/resume-prep-exit-code.txt"
OUTPUT_PREFIX="output/exp071-production-pm-cd4/${RUN_ID}"

if [[ "$RUN_ID" != "20260823T083438Z" ]]; then
  echo "resume script is pinned to run 20260823T083438Z" >&2
  exit 65
fi
if [[ "$EXP065_RUN_ID" != "20260823T040540Z" ]]; then
  echo "resume script is pinned to experiment-065 run 20260823T040540Z" >&2
  exit 66
fi
if [[ ! -d "$WORK" || ! -d "$REPO" || ! -d "$OUT" ]]; then
  echo "preserved experiment-071 work tree is missing" >&2
  exit 67
fi

exec 8>/work/exp061.lock
if ! flock -n 8; then
  echo "experiment 061 lock is busy" >&2
  exit 70
fi
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "global VCC heavy lock is busy" >&2
  exit 71
fi
if pgrep -af \
  '[c]ell-eval2|[b]uild_candidate.py|[b]uild_joint_candidate.py|[m]ake_submission_pm_cd4.py|[v]cc prep' \
  >"${WORK}/resume-preexisting-heavy-processes.txt"
then
  echo "refusing to overlap an existing VCC build, scorer, or prep" >&2
  cat "${WORK}/resume-preexisting-heavy-processes.txt" >&2
  exit 72
fi
if [[ ! -f "$H5AD" || -L "$H5AD" ]]; then
  echo "preserved H5AD is missing or is a symlink" >&2
  exit 68
fi
if [[ "$(stat -c %s "$H5AD")" != "$EXPECTED_H5AD_BYTES" ]]; then
  echo "preserved H5AD size mismatch" >&2
  exit 69
fi
printf '%s  %s\n' "$EXPECTED_H5AD_SHA256" "$H5AD" | sha256sum --check -

mem_kib="$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
if (( mem_kib < 120000000 )); then
  echo "worker has less than 120000000 KiB RAM required for vcc prep" >&2
  exit 73
fi

exec > >(tee -a "$LOG") 2>&1

finish() {
  local status=$?
  trap - EXIT
  printf '%s\n' "$status" >"$EXIT_FILE"
  aws s3 cp "$LOG" "s3://${BUCKET}/${OUTPUT_PREFIX}/resume-prep.log" \
    --only-show-errors || true
  aws s3 cp "$EXIT_FILE" \
    "s3://${BUCKET}/${OUTPUT_PREFIX}/resume-prep-exit-code.txt" \
    --only-show-errors || true
  exit "$status"
}
trap finish EXIT

cd "$REPO"
export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"

rm -f "$DRY_RUN_VCC" "$PARTIAL"
uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib "$PREP_MAX_RSS_GIB" \
  --min-free-percent 0 \
  --poll-seconds 2 \
  -- \
  /usr/bin/time -v vcc prep "$H5AD" \
    -g data/challenge/gene_names.csv \
    --perts data/challenge/pert_counts.csv \
    -o "$DRY_RUN_VCC" \
    --dry-run \
  2>&1 | tee "$OUT/vcc-prep-dry-run-resume.log"
rm -f "$DRY_RUN_VCC"

if [[ -e "$VCC" ]]; then
  echo "refusing to overwrite existing final package" >&2
  exit 74
fi
uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib "$PREP_MAX_RSS_GIB" \
  --min-free-percent 0 \
  --poll-seconds 2 \
  -- \
  /usr/bin/time -v vcc prep "$H5AD" \
    -g data/challenge/gene_names.csv \
    --perts data/challenge/pert_counts.csv \
    -o "$PARTIAL" \
  2>&1 | tee "$OUT/vcc-prep-resume.log"
if [[ ! -s "$PARTIAL" ]]; then
  echo "vcc prep returned without a nonempty package" >&2
  exit 75
fi
mv "$PARTIAL" "$VCC"
sha256sum "$VCC" >"$OUT/v10-pm-cd4.vcc.sha256"

(
  cd "$OUT"
  find . -type f \
    ! -name '*.h5ad' \
    ! -name 'artifact-sha256.txt' \
    ! -name 'PACKAGE-SUCCESS.json' \
    ! -name 'resume-prep.log' \
    ! -name 'resume-prep-exit-code.txt' \
    ! -name '*.partial-*' \
    -print0 |
    sort -z |
    xargs -0 sha256sum >artifact-sha256.txt
)

python3 - "$OUT" "$RUN_ID" "$EXP065_RUN_ID" <<'PY'
import hashlib
import json
import os
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
run_id, exp065_run_id = sys.argv[2:4]
package = root / "v10-pm-cd4.vcc"
manifest = root / "artifact-sha256.txt"

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()

marker = {
    "status": "package_complete",
    "authorization_scope": "exploratory_packaging_only",
    "challenge_submission_authorized": False,
    "resumed_from_preserved_h5ad": True,
    "run_id": run_id,
    "exp065_run_id": exp065_run_id,
    "source_h5ad": {
        "bytes": (root / "v10-pm-cd4.h5ad").stat().st_size,
        "sha256": sha256(root / "v10-pm-cd4.h5ad"),
    },
    "package": {
        "path": package.name,
        "bytes": package.stat().st_size,
        "sha256": sha256(package),
    },
    "artifact_manifest_sha256": sha256(manifest),
}
temporary = root / f".PACKAGE-SUCCESS.json.partial-{os.getpid()}"
temporary.write_text(json.dumps(marker, indent=2, sort_keys=True) + "\n")
os.replace(temporary, root / "PACKAGE-SUCCESS.json")
PY

(
  cd "$OUT"
  sha256sum --check v10-pm-cd4.vcc.sha256
)
aws s3 sync "$OUT/" "s3://${BUCKET}/${OUTPUT_PREFIX}/" \
  --exclude '*.h5ad' \
  --exclude '*.vcc' \
  --exclude '*.vcc.sha256' \
  --exclude '*.partial-*' \
  --only-show-errors
aws s3 cp "$VCC" \
  "s3://${BUCKET}/${OUTPUT_PREFIX}/v10-pm-cd4.vcc" \
  --only-show-errors
aws s3 cp "$OUT/v10-pm-cd4.vcc.sha256" \
  "s3://${BUCKET}/${OUTPUT_PREFIX}/v10-pm-cd4.vcc.sha256" \
  --only-show-errors
