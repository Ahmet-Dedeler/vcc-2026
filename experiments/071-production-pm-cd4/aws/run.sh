#!/usr/bin/env bash
set -Eeuo pipefail

atomic_vcc_prep() {
  local final_output="$1"
  shift
  local run_token="${RUN_ID:-test}"
  local partial_output="${final_output%.vcc}.partial-${run_token}.vcc"
  if [[ -e "$final_output" ]]; then
    echo "refusing to overwrite existing package ${final_output}" >&2
    return 90
  fi
  rm -f "$partial_output"
  if "$@" -o "$partial_output"; then
    :
  else
    local status=$?
    rm -f "$partial_output"
    return "$status"
  fi
  if [[ ! -s "$partial_output" ]]; then
    rm -f "$partial_output"
    echo "vcc prep returned success without a nonempty package" >&2
    return 91
  fi
  mv "$partial_output" "$final_output"
}

if [[ "${EXP071_LIBRARY_ONLY:-}" == "1" ]]; then
  return 0 2>/dev/null || exit 0
fi

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:?SOURCE_PREFIX is required}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:?OUTPUT_PREFIX is required}"
EXP065_RUN_ID="${EXP065_RUN_ID:?EXP065_RUN_ID is required}"
EXP065_PREFIX="${EXP065_PREFIX:?EXP065_PREFIX is required}"
SOURCE_ARCHIVE_SHA256="${SOURCE_ARCHIVE_SHA256:?SOURCE_ARCHIVE_SHA256 is required}"
RUNNER_SHA256="${RUNNER_SHA256:?RUNNER_SHA256 is required}"
SOURCE_MANIFEST_SHA256="${SOURCE_MANIFEST_SHA256:?SOURCE_MANIFEST_SHA256 is required}"
AUTHORIZATION_SHA256="${AUTHORIZATION_SHA256:?AUTHORIZATION_SHA256 is required}"
EXP061_REPO="${EXP061_REPO:-/work/exp061/repo}"
WORK="/work/exp071-${RUN_ID}"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
AUTH="${WORK}/authorization"
OUT="${WORK}/output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"
SUCCESS_MARKER="${OUT}/PACKAGE-SUCCESS.json"
MAX_RSS_GIB="${MAX_RSS_GIB:-55}"
PREP_MAX_RSS_GIB="${PREP_MAX_RSS_GIB:-110}"
MIN_FREE_GIB="${MIN_FREE_GIB:-44}"
CONTROLS_SHA256="329a22bc29cac4f43ad0b846b9bb1a383f7a945c83a8cb0450c465f6dc5cc4b5"
WHITENER_SHA256="5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3"
WHITENER_JSON_SHA256="a77f9425ac524cdcf3a53a1744520e1e63cb8b925ef77fae2864c4d8b72337b6"

if [[ ! "$RUN_ID" =~ ^[0-9]{8}T[0-9]{6}Z$ ]]; then
  echo "invalid RUN_ID: ${RUN_ID}" >&2
  exit 64
fi
if [[ ! "$EXP065_RUN_ID" =~ ^[0-9]{8}T[0-9]{6}Z$ ]]; then
  echo "invalid EXP065_RUN_ID: ${EXP065_RUN_ID}" >&2
  exit 65
fi
if [[ "$SOURCE_PREFIX" != "input/exp071-production-pm-cd4/${RUN_ID}" ]]; then
  echo "SOURCE_PREFIX is not bound to RUN_ID" >&2
  exit 66
fi
if [[ "$OUTPUT_PREFIX" != "output/exp071-production-pm-cd4/${RUN_ID}" ]]; then
  echo "OUTPUT_PREFIX is not bound to RUN_ID" >&2
  exit 67
fi
if [[ "$EXP065_PREFIX" != "output/exp065-pm-cd4-hct-count/${EXP065_RUN_ID}" ]]; then
  echo "EXP065_PREFIX is not bound to EXP065_RUN_ID" >&2
  exit 68
fi

rm -rf "$WORK"
mkdir -p "$REPO" "$INPUT" "$AUTH/exp065" "$OUT"
rm -f "$EXIT_FILE"

# Hold both locks for the complete build and packaging pass. The exp061 lock
# protects the reused donor tree; the global lock prevents any other scorer or
# memory-heavy VCC task from sharing this worker.
exec 8>/work/exp061.lock
if ! flock -n 8; then
  echo "experiment 061 still holds /work/exp061.lock" >&2
  exit 69
fi
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 70
fi

if pgrep -af \
  '[c]ell-eval2|[b]uild_candidate.py|[b]uild_joint_candidate.py|[m]ake_submission_pm_cd4.py|[v]cc prep' \
  >"${WORK}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap an existing VCC build or scorer" >&2
  cat "${WORK}/preexisting-heavy-processes.txt" >&2
  exit 71
fi

exec > >(tee -a "$LOG") 2>&1

upload_outputs() {
  local status="$1"
  aws s3 cp "$LOG" "s3://${BUCKET}/${OUTPUT_PREFIX}/cloud.log" \
    --only-show-errors || true
  if [[ -s "$EXIT_FILE" ]]; then
    aws s3 cp "$EXIT_FILE" "s3://${BUCKET}/${OUTPUT_PREFIX}/exit-code.txt" \
      --only-show-errors || true
  fi
  if [[ -d "$OUT" ]]; then
    aws s3 sync "$OUT/" "s3://${BUCKET}/${OUTPUT_PREFIX}/" \
      --exclude '*.h5ad' \
      --exclude '*.vcc' \
      --exclude '*.vcc.sha256' \
      --exclude '*.partial-*' \
      --only-show-errors || true
  fi
  if [[ "$status" == "0" && -s "$SUCCESS_MARKER" ]]; then
    local package="${OUT}/v10-pm-cd4.vcc"
    local checksum="${OUT}/v10-pm-cd4.vcc.sha256"
    if [[ -s "$package" && -s "$checksum" ]]; then
      (
        cd "$OUT"
        sha256sum --check "$(basename "$checksum")"
      ) || return 0
      aws s3 cp "$package" \
        "s3://${BUCKET}/${OUTPUT_PREFIX}/$(basename "$package")" \
        --only-show-errors || true
      aws s3 cp "$checksum" \
        "s3://${BUCKET}/${OUTPUT_PREFIX}/$(basename "$checksum")" \
        --only-show-errors || true
    fi
  fi
}

finish() {
  local status=$?
  trap - EXIT
  printf '%s\n' "$status" > "$EXIT_FILE"
  upload_outputs "$status"
  exit "$status"
}
trap finish EXIT

available_kib="$(df -Pk /work | awk 'NR==2 {print $4}')"
required_kib="$((MIN_FREE_GIB * 1024 * 1024))"
if (( available_kib < required_kib )); then
  echo "only $((available_kib / 1024 / 1024)) GiB free; need ${MIN_FREE_GIB} GiB" >&2
  exit 72
fi
mem_kib="$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
if (( mem_kib < 120000000 )); then
  echo "worker has less than 120000000 KiB RAM required for vcc prep" >&2
  exit 76
fi

required_artifacts=(
  "exit-code.txt"
  "artifact-sha256.txt"
  "results.json"
  "build-report.json"
  "baseline-agg_results.csv"
  "candidate-run/agg_results.csv"
  "candidate.h5ad.sha256"
  "score-memory-guard.json"
  "runtime.json"
)

aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source-manifest.json" \
  "${WORK}/source-manifest.json" --only-show-errors
printf '%s  %s\n' "$SOURCE_MANIFEST_SHA256" "${WORK}/source-manifest.json" |
  sha256sum --check -
aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "${WORK}/source.tar.gz" --only-show-errors
aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/authorization/authorization.json" \
  "${AUTH}/authorization.json" --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/${SOURCE_PREFIX}/authorization/terminal-s3-manifest.json" \
  "${AUTH}/terminal-s3-manifest.json" --only-show-errors
for relative in "${required_artifacts[@]}"; do
  mkdir -p "$(dirname "${AUTH}/exp065/${relative}")"
  aws s3 cp \
    "s3://${BUCKET}/${SOURCE_PREFIX}/authorization/exp065/${relative}" \
    "${AUTH}/exp065/${relative}" \
    --only-show-errors
done

printf '%s  %s\n' "$RUNNER_SHA256" "${BASH_SOURCE[0]}" | sha256sum --check -
printf '%s  %s\n' "$SOURCE_ARCHIVE_SHA256" "${WORK}/source.tar.gz" |
  sha256sum --check -
printf '%s  %s\n' "$AUTHORIZATION_SHA256" "${AUTH}/authorization.json" |
  sha256sum --check -

python3 - \
  "${WORK}/source-manifest.json" \
  "${WORK}/source.tar.gz" \
  "${BASH_SOURCE[0]}" \
  "$AUTH" \
  "$RUN_ID" \
  "$SOURCE_PREFIX" \
  "$EXP065_RUN_ID" \
  "$EXP065_PREFIX" <<'PY'
import hashlib
import json
import pathlib
import sys

manifest_path = pathlib.Path(sys.argv[1])
archive = pathlib.Path(sys.argv[2])
runner = pathlib.Path(sys.argv[3])
authorization = pathlib.Path(sys.argv[4])
run_id, source_prefix, exp065_run_id, exp065_prefix = sys.argv[5:9]
manifest = json.loads(manifest_path.read_text())

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()

expected_header = {
    "schema": 1,
    "experiment": "071-production-pm-cd4",
    "run_id": run_id,
    "source_prefix": source_prefix,
    "exp065_run_id": exp065_run_id,
    "exp065_prefix": exp065_prefix,
    "challenge_submission_authorized": False,
}
for key, expected in expected_header.items():
    if manifest.get(key) != expected:
        raise SystemExit(f"source manifest {key} mismatch")

required = {
    "source.tar.gz",
    "run.sh",
    "authorization/authorization.json",
    "authorization/terminal-s3-manifest.json",
    "authorization/exp065/exit-code.txt",
    "authorization/exp065/artifact-sha256.txt",
    "authorization/exp065/results.json",
    "authorization/exp065/build-report.json",
    "authorization/exp065/baseline-agg_results.csv",
    "authorization/exp065/candidate-run/agg_results.csv",
    "authorization/exp065/candidate.h5ad.sha256",
    "authorization/exp065/score-memory-guard.json",
    "authorization/exp065/runtime.json",
}
entries = {}
for entry in manifest.get("files", []):
    relative = pathlib.PurePosixPath(entry.get("path", ""))
    if (
        relative.is_absolute()
        or not relative.parts
        or any(part in ("", ".", "..") for part in relative.parts)
    ):
        raise SystemExit(f"unsafe source manifest path: {relative}")
    name = relative.as_posix()
    if name in entries:
        raise SystemExit(f"duplicate source manifest path: {name}")
    entries[name] = entry
if set(entries) != required:
    raise SystemExit(
        f"source manifest file set mismatch: "
        f"missing={sorted(required - set(entries))}, "
        f"unexpected={sorted(set(entries) - required)}"
    )

local_paths = {
    "source.tar.gz": archive,
    "run.sh": runner,
}
for name in required - set(local_paths):
    local_paths[name] = authorization / pathlib.PurePosixPath(name).relative_to(
        "authorization"
    )
for name, path in local_paths.items():
    entry = entries[name]
    if not path.is_file() or path.is_symlink():
        raise SystemExit(f"missing regular staged source file: {name}")
    if entry.get("size") != path.stat().st_size:
        raise SystemExit(f"source size mismatch: {name}")
    if entry.get("sha256") != sha256(path):
        raise SystemExit(f"source hash mismatch: {name}")
PY

if tar -tzf "${WORK}/source.tar.gz" |
  awk 'BEGIN {bad=0} /^\// || /(^|\/)\.\.(\/|$)/ {bad=1} END {exit bad}'
then
  :
else
  echo "source archive contains an unsafe path" >&2
  exit 73
fi
tar -xzf "${WORK}/source.tar.gz" -C "$REPO"

python3 - "${WORK}/source-manifest.json" "$REPO" <<'PY'
import hashlib
import json
import pathlib
import sys

manifest = json.loads(pathlib.Path(sys.argv[1]).read_text())
root = pathlib.Path(sys.argv[2]).resolve()
entries = manifest.get("archive_members")
if not isinstance(entries, list) or not entries:
    raise SystemExit("source manifest archive_members is missing")
expected = {
    "experiments/071-production-pm-cd4/validate_authorization.py",
    "scripts/make_submission_pm_cd4.py",
    "scripts/test_make_submission_pm_cd4.py",
    "scripts/run_with_memory_guard.py",
}
seen = set()
for entry in entries:
    relative = pathlib.PurePosixPath(entry.get("path", ""))
    if (
        relative.is_absolute()
        or not relative.parts
        or any(part in ("", ".", "..") for part in relative.parts)
    ):
        raise SystemExit(f"unsafe archive member path: {relative}")
    name = relative.as_posix()
    if name in seen:
        raise SystemExit(f"duplicate archive member: {name}")
    seen.add(name)
    path = root.joinpath(*relative.parts)
    if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(root):
        raise SystemExit(f"missing regular archive member: {name}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != entry.get("sha256") or path.stat().st_size != entry.get("size"):
        raise SystemExit(f"archive member provenance mismatch: {name}")
if seen != expected:
    raise SystemExit(
        f"archive member set mismatch: "
        f"missing={sorted(expected - seen)}, unexpected={sorted(seen - expected)}"
    )
PY

python3 "$REPO/experiments/071-production-pm-cd4/validate_authorization.py" \
  --artifact-root "$AUTH/exp065" \
  --terminal-manifest "$AUTH/terminal-s3-manifest.json" \
  --exp065-run-id "$EXP065_RUN_ID" \
  --expected-prefix "$EXP065_PREFIX" \
  --output "$OUT/authorization-revalidated.json"
cmp "$AUTH/authorization.json" "$OUT/authorization-revalidated.json"
cp "$AUTH/authorization.json" "$OUT/authorization-staged.json"
cp "$AUTH/terminal-s3-manifest.json" "$OUT/exp065-terminal-s3-manifest.json"
cp "${WORK}/source-manifest.json" "$OUT/source-manifest.json"

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

dnf install -y gzip tar time unzip
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

aws s3 cp \
  "s3://${BUCKET}/input/production-challenge/vcc_2026_controls.zip" \
  "${INPUT}/vcc_2026_controls.zip" --only-show-errors
printf '%s  %s\n' "$CONTROLS_SHA256" "${INPUT}/vcc_2026_controls.zip" |
  sha256sum --check -
mkdir -p "$REPO/data/challenge"
unzip -q "${INPUT}/vcc_2026_controls.zip" -d "$REPO/data/challenge"

aws s3 cp \
  "s3://${BUCKET}/input/production-challenge/whitener_union.npz" \
  "${INPUT}/whitener_union.npz" --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/input/production-challenge/whitener_union.json" \
  "${INPUT}/whitener_union.json" --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/input/clean-hct/orion_hct116/" \
  "${INPUT}/orion_hct116/" --recursive --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/input/exp060-cd4-union-k562/cd4_compact.tar.gz" \
  "${INPUT}/cd4_compact.tar.gz" --only-show-errors
mkdir -p "${INPUT}/cd4"
tar -xzf "${INPUT}/cd4_compact.tar.gz" -C "${INPUT}/cd4"
CD4_ROOT="${INPUT}/cd4/cd4_compact"

printf '%s  %s\n' "$WHITENER_SHA256" "${INPUT}/whitener_union.npz" |
  sha256sum --check -
printf '%s  %s\n' "$WHITENER_JSON_SHA256" "${INPUT}/whitener_union.json" |
  sha256sum --check -

cd "$REPO"
uv sync --frozen
uv tool install --force "vcc-cli==0.1.0"
vcc version --json | tee "$OUT/vcc-version.json"

uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib "$PREP_MAX_RSS_GIB" \
  --min-free-percent 0 \
  --poll-seconds 2 \
  -- \
  uv run python scripts/test_make_submission_pm_cd4.py \
  2>&1 | tee "$OUT/tests.log"

H5AD="$OUT/v10-pm-cd4.h5ad"
VCC="$OUT/v10-pm-cd4.vcc"
DRY_RUN_VCC="$OUT/.v10-pm-cd4.dry-run-${RUN_ID}.vcc"

uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib "$MAX_RSS_GIB" \
  --min-free-percent 0 \
  --poll-seconds 3 \
  -- \
  /usr/bin/time -v uv run python scripts/make_submission_pm_cd4.py \
    --bundle data/challenge \
    --k562 "${EXP061_REPO}/data/donors/k562_gwps" \
    --hct116 "$INPUT/orion_hct116" \
    --hek293t "${EXP061_REPO}/data/donors/orion_hek293t" \
    --perturbmatch "${EXP061_REPO}/data/external/perturbmatch" \
    --cd4-root "$CD4_ROOT" \
    --whitener "${INPUT}/whitener_union.npz" \
    --whitener-metadata "${INPUT}/whitener_union.json" \
    --expected-whitener-sha256 "$WHITENER_SHA256" \
    --output "$H5AD" \
    --max-rss-gib "$MAX_RSS_GIB" \
  2>&1 | tee "$OUT/build.log"

sha256sum "$H5AD" > "$OUT/v10-pm-cd4.h5ad.sha256"

rm -f "$DRY_RUN_VCC"
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
  2>&1 | tee "$OUT/vcc-prep-dry-run.log"
rm -f "$DRY_RUN_VCC"

atomic_vcc_prep "$VCC" \
  uv run python scripts/run_with_memory_guard.py \
    --max-rss-gib "$PREP_MAX_RSS_GIB" \
    --min-free-percent 0 \
    --poll-seconds 2 \
    -- \
    /usr/bin/time -v vcc prep "$H5AD" \
      -g data/challenge/gene_names.csv \
      --perts data/challenge/pert_counts.csv \
  2>&1 | tee "$OUT/vcc-prep.log"

sha256sum "$VCC" > "$OUT/v10-pm-cd4.vcc.sha256"
(
  cd "$OUT"
  find . -type f \
    ! -name '*.h5ad' \
    ! -name 'artifact-sha256.txt' \
    ! -name 'PACKAGE-SUCCESS.json' \
    ! -name '*.partial-*' \
    -print0 |
    sort -z |
    xargs -0 sha256sum > artifact-sha256.txt
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
    "run_id": run_id,
    "exp065_run_id": exp065_run_id,
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
