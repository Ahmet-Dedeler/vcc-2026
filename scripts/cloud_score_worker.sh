#!/usr/bin/env bash
# cloud_score_worker.sh — the remote half of `scripts/cloud_score.sh`. Runs ON the EC2 worker.
#
# Never invoke this by hand. `cloud_score.sh` uploads it to S3, and the instance's user-data
# fetches it, checks its sha256 against the value the launcher recorded, and runs it.
#
# CONTRACT
#   in   s3://$BUCKET/$PREFIX/job.json      the job manifest (real, preds, flags)
#        s3://$BUCKET/$PREFIX/payload/...   optional files, restored at the same relative path
#        s3://$BUCKET/$PREFIX/source.tar.gz optional repo source, only when a prologue is set
#        s3://$BUCKET/$PREFIX/prologue.sh   optional build step, run in $REPO before scoring
#   out  s3://$BUCKET/$OUT/<name>/run/      one `cell-eval2 run` output dir per prediction
#        s3://$BUCKET/$OUT/<name>/score.csv the scaled six from `cell-eval2 score`
#        s3://$BUCKET/$OUT/bundle/          the real bundle (baseline + replicate anchors)
#        s3://$BUCKET/$OUT/sha256.txt       sha256 of every artifact returned
#        s3://$BUCKET/$OUT/worker.log       the full transcript
#        s3://$BUCKET/$OUT/STATUS           `DONE <rc>` or `FAILED <rc>` — the launcher's gate
#
# SELF-TERMINATION. Three independent guards, because the 13-day orphan of experiment 061 had
# none of them:
#   (a) the instance is launched with `instance-initiated-shutdown-behavior=terminate`;
#   (b) user-data arms `shutdown -h +$MAX_MINUTES` BEFORE this script starts, so a hang, an
#       OOM kill or a crashed coordinator still ends in a terminated instance;
#   (c) this script's EXIT trap calls `shutdown -h now` on every path, success or failure.
set -Eeuo pipefail

: "${BUCKET:?}" "${PREFIX:?}" "${OUT:?}" "${AWS_REGION:?}"
WORK="${WORK:-/work/cloud-score}"
REPO="${WORK}/repo"
LOG="${WORK}/worker.log"
STATUS="${WORK}/STATUS"

mkdir -p "$WORK"
exec 9>"${WORK}.lock"
flock -n 9 || { echo "another cloud_score worker holds ${WORK}.lock" >&2; exit 70; }
exec > >(tee -a "$LOG") 2>&1

push() {
  aws s3 cp "$LOG" "s3://${BUCKET}/${OUT}/worker.log" --only-show-errors || true
  [[ -s "$STATUS" ]] && aws s3 cp "$STATUS" "s3://${BUCKET}/${OUT}/STATUS" --only-show-errors || true
}

# `COMPLETED` is set to 1 by the LAST line of this script and by nothing else.
#
# Do not infer completion from `$?` in the EXIT trap. A spot interruption sends SIGTERM;
# the trap then runs with `$? == 0` because the last command it interrupted had succeeded,
# and the worker writes `DONE 0` for a run that got two thirds of the way through a
# `cell-eval2 baseline`. That is exactly what happened on run `ka-20260904` (AWS reclaimed
# the instance 42 min in, "Service initiated (2026-09-04 01:49:50 GMT)") — the marker said
# DONE and there was not a single output file. A launcher that trusts the marker would have
# reported a passing known-answer row from a run that never scored anything.
COMPLETED=0

finish() {
  rc=$?
  trap - EXIT INT TERM
  if [[ "$COMPLETED" -eq 1 && "$rc" -eq 0 ]]; then
    printf 'DONE %s\n' "$rc" >"$STATUS"
  else
    printf 'FAILED %s (completed=%s)\n' "$rc" "$COMPLETED" >"$STATUS"
  fi
  push
  # Guard (c). `poweroff` is the fallback for images where `shutdown` is a symlink that
  # refuses without a controlling terminal.
  shutdown -h now || poweroff || true
  exit "$rc"
}
trap finish EXIT
trap 'echo "[worker] caught a signal — this run did NOT complete" >&2; exit 143' INT TERM

echo "[worker] $(date -u +%FT%TZ) starting on $(hostname); $(nproc) vCPU, $(awk '/MemTotal/{printf "%.0f GiB", $2/1048576}' /proc/meminfo)"

dnf install -y -q gzip tar time jq >/dev/null 2>&1 || dnf install -y gzip tar time jq
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"

aws s3 cp "s3://${BUCKET}/${PREFIX}/job.json" "${WORK}/job.json" --only-show-errors
cat "${WORK}/job.json"

CE_VERSION="$(jq -r '.cell_eval2_version' "${WORK}/job.json")"
PRESET="$(jq -r '.preset' "${WORK}/job.json")"
PERT_COL="$(jq -r '.pert_col' "${WORK}/job.json")"
CONTROL="$(jq -r '.control' "${WORK}/job.json")"
REAL_URI="$(jq -r '.real' "${WORK}/job.json")"
HAS_PROLOGUE="$(jq -r '.prologue // "" | length > 0' "${WORK}/job.json")"

mkdir -p "${WORK}/in" "${WORK}/out" "${WORK}/preds"
aws s3 cp "$REAL_URI" "${WORK}/in/real.h5ad" --only-show-errors
REAL="${WORK}/in/real.h5ad"
echo "[worker] real: $(sha256sum "$REAL")"

# ---------------------------------------------------------------- optional prologue (a build)
if [[ "$HAS_PROLOGUE" == "true" ]]; then
  mkdir -p "$REPO"
  aws s3 cp "s3://${BUCKET}/${PREFIX}/source.tar.gz" "${WORK}/source.tar.gz" --only-show-errors
  tar -xzf "${WORK}/source.tar.gz" -C "$REPO"
  if aws s3 ls "s3://${BUCKET}/${PREFIX}/payload/" >/dev/null 2>&1; then
    aws s3 sync "s3://${BUCKET}/${PREFIX}/payload/" "$REPO/" --only-show-errors
  fi
  aws s3 cp "s3://${BUCKET}/${PREFIX}/prologue.sh" "${WORK}/prologue.sh" --only-show-errors
  chmod +x "${WORK}/prologue.sh"
  ( cd "$REPO" && WORK="$WORK" REAL="$REAL" PREDS_DIR="${WORK}/preds" \
      /usr/bin/time -v bash "${WORK}/prologue.sh" ) 2>&1 | tee "${WORK}/prologue.log"
  aws s3 cp "${WORK}/prologue.log" "s3://${BUCKET}/${OUT}/prologue.log" --only-show-errors
  # A prologue may also want its own artifacts back (audits, JSON, logs).
  if [[ -d "${WORK}/prologue-out" ]]; then
    aws s3 sync "${WORK}/prologue-out/" "s3://${BUCKET}/${OUT}/prologue-out/" --only-show-errors
  fi
  # The build is the expensive, non-idempotent half. Park it in S3 immediately so a spot
  # reclaim costs only the scoring, which `--bundle` + `--pred s3://…` can then resume.
  if compgen -G "${WORK}/preds/*.h5ad" >/dev/null; then
    # NOT under `input/`: the instance role is least-privilege — GetObject on `input/*`,
    # PutObject on `output/*` only — so a worker cannot write to its own input prefix. A
    # first version parked here and lost an 18-minute build to AccessDenied. `cloud_score.sh`
    # re-stages an `output/` URI into the next run's input prefix on resume.
    aws s3 sync "${WORK}/preds/" "s3://${BUCKET}/${OUT}/built-preds/" \
      --exclude '*' --include '*.h5ad' --only-show-errors
    echo "[worker] built candidates parked at s3://${BUCKET}/${OUT}/built-preds/"
  fi
  push
fi

# ---------------------------------------------------------------- fetch the predictions
NAMES=()
while IFS=$'\t' read -r name uri; do
  [[ -n "$name" ]] || continue
  NAMES+=("$name")
  if [[ "$uri" == "PROLOGUE" ]]; then
    test -s "${WORK}/preds/${name}.h5ad" || { echo "prologue did not write ${name}.h5ad" >&2; exit 74; }
  else
    aws s3 cp "$uri" "${WORK}/preds/${name}.h5ad" --only-show-errors
  fi
  echo "[worker] pred ${name}: $(sha256sum "${WORK}/preds/${name}.h5ad")"
done < <(jq -r '.preds[] | [.name, .uri] | @tsv' "${WORK}/job.json")

CE=(uvx --python 3.12 --from "cell-eval2==${CE_VERSION}" cell-eval2)

# ---------------------------------------------------------------- 1. run, per prediction
#
# RUNS COME FIRST, DELIBERATELY. `run` is what produces `agg_results.csv` — the RAW six
# members — and the repo's existing local->board factors (164's `reach` 0.4277, 185's 1.290,
# `docs/BOARD-LEDGER.md` §4) are all RAW-delta ratios, so the raw table is the deliverable on
# the short path. `prep-real-bundle` is ~10x one `run` (it averages 5 disjoint half-splits of
# the reference, each needing DE on both halves) and only the SCALED table needs it. Putting
# the bundle first meant a 6-hour wait before the first readable number, and an interruption
# in that window returned nothing at all.
for name in "${NAMES[@]}"; do
  d="${WORK}/out/${name}"
  mkdir -p "$d"
  echo "[worker] === cell-eval2 run: ${name} ==="
  /usr/bin/time -v "${CE[@]}" run \
    -ap "${WORK}/preds/${name}.h5ad" -ar "$REAL" -o "${d}/run" \
    --preset "$PRESET" --pert-col "$PERT_COL" --control "$CONTROL" \
    2>&1 | tee "${d}/run.log"

  # The DE engine is NOT interchangeable: cell-eval2 warns "DE numbers differ between
  # engines", and the 095 h1-test artifact this worker's known-answer row reproduces was
  # scored with `resolved_de_backend: scanpy`. Fail loudly rather than silently compare
  # numbers from a different engine.
  got_backend="$(jq -r '.resolved_de_backend' "${d}/run/run_meta.json")"
  got_version="$(jq -r '.cell_eval2_version' "${d}/run/run_meta.json")"
  [[ "$got_backend" == "scanpy" ]] || { echo "DE backend is ${got_backend}, not scanpy" >&2; exit 77; }
  [[ "$got_version" == "$CE_VERSION" ]] || { echo "cell-eval2 is ${got_version}" >&2; exit 77; }

  aws s3 sync "${d}/" "s3://${BUCKET}/${OUT}/${name}/" --only-show-errors
  push
done

# --------------------------------------------- 2. the real bundle: reuse, wait for, or build
# It depends ONLY on the reference, so one bundle per surface serves every candidate ever
# scored on it. `--bundle` may point at a bundle another worker is still building: this polls
# for its manifest rather than duplicating six hours of work, and falls back to building its
# own if the deadline passes. `score --real-bundle` REFUSES a submission that does not match
# the bundle's identity, so reuse cannot silently mis-pair.
BUNDLE_URI="$(jq -r '.bundle // ""' "${WORK}/job.json")"
BUNDLE_WAIT_MIN="$(jq -r '.bundle_wait_minutes // 0' "${WORK}/job.json")"
mkdir -p "${WORK}/out/bundle"
if [[ -n "$BUNDLE_URI" ]]; then
  deadline=$(( $(date -u +%s) + BUNDLE_WAIT_MIN * 60 ))
  while :; do
    if aws s3 ls "${BUNDLE_URI%/}/manifest.json" >/dev/null 2>&1; then
      echo "[worker] reusing bundle ${BUNDLE_URI}"
      aws s3 sync "${BUNDLE_URI%/}/" "${WORK}/out/bundle/" --only-show-errors
      break
    fi
    [[ "$(date -u +%s)" -lt "$deadline" ]] || { echo "[worker] bundle wait expired"; break; }
    echo "[worker] $(date -u +%FT%TZ) waiting for ${BUNDLE_URI%/}/manifest.json"
    sleep 120
  done
fi

if [[ ! -f "${WORK}/out/bundle/manifest.json" ]]; then
  # `prep-real-bundle --baseline` wants the baseline PREDICTION h5ad, not the `baseline -o`
  # output directory. The pre-2026-09-22 AGENTS.md §6.4 snippet passed the directory; that fails with
  # `IsADirectoryError` after `cell-eval2 baseline` has already burned 3.5 h — which is a fair
  # part of why the four-step pipeline was never actually run end to end in this project.
  # `baseline --save-pred` writes it, and a previously saved one can be supplied instead.
  BASELINE_PRED_URI="$(jq -r '.baseline_pred // ""' "${WORK}/job.json")"
  BASELINE_PRED="${WORK}/out/baseline-pred.h5ad"
  if [[ -n "$BASELINE_PRED_URI" ]]; then
    echo "[worker] supplied baseline prediction ${BASELINE_PRED_URI}"
    aws s3 cp "$BASELINE_PRED_URI" "$BASELINE_PRED" --only-show-errors
  else
    echo "[worker] === cell-eval2 baseline ==="
    /usr/bin/time -v "${CE[@]}" baseline \
      -ar "$REAL" -o "${WORK}/out/baseline" --save-pred "$BASELINE_PRED" \
      --preset "$PRESET" --pert-col "$PERT_COL" --control "$CONTROL" \
      2>&1 | tee "${WORK}/out/baseline.log"
  fi
  echo "[worker] baseline prediction: $(sha256sum "$BASELINE_PRED")"
  echo "[worker] === cell-eval2 prep-real-bundle ==="
  /usr/bin/time -v "${CE[@]}" prep-real-bundle \
    --real "$REAL" --baseline "$BASELINE_PRED" -o "${WORK}/out/bundle" \
    --preset "$PRESET" --pert-col "$PERT_COL" --control "$CONTROL" --force \
    2>&1 | tee "${WORK}/out/bundle.log"
  [[ -d "${WORK}/out/baseline" ]] && aws s3 sync "${WORK}/out/baseline/" \
    "s3://${BUCKET}/${OUT}/baseline/" --only-show-errors || true
fi
aws s3 sync "${WORK}/out/bundle/" "s3://${BUCKET}/${OUT}/bundle/" --only-show-errors
push

# ---------------------------------------------------------------- 3. score, per prediction
for name in "${NAMES[@]}"; do
  d="${WORK}/out/${name}"
  agg="${d}/run/agg_results.csv"
  [[ -f "$agg" ]] || agg="${d}/run/agg.csv"
  echo "[worker] === cell-eval2 score: ${name} ==="
  /usr/bin/time -v "${CE[@]}" score \
    --user-agg "$agg" --real-bundle "${WORK}/out/bundle" -o "${d}/score.csv" \
    2>&1 | tee "${d}/score.log"
  cat "${d}/score.csv" || true
  aws s3 sync "${d}/" "s3://${BUCKET}/${OUT}/${name}/" --only-show-errors
  push
done

# ---------------------------------------------------------------- 4. hashes, runtime, done
( cd "${WORK}/out" && find . -type f \( -name '*.csv' -o -name '*.json' -o -name '*.parquet' \) \
    -print0 | sort -z | xargs -0 sha256sum ) >"${WORK}/sha256.txt"
cat "${WORK}/sha256.txt"
aws s3 cp "${WORK}/sha256.txt" "s3://${BUCKET}/${OUT}/sha256.txt" --only-show-errors

"${CE[@]}" --help >/dev/null 2>&1 || true
uvx --python 3.12 --from "cell-eval2==${CE_VERSION}" python -c "
import importlib.metadata as m, json, platform
print(json.dumps({'python': platform.python_version(),
                  'cell_eval2': m.version('cell-eval2'),
                  'anndata': m.version('anndata'),
                  'numpy': m.version('numpy'),
                  'scipy': m.version('scipy'),
                  'polars': m.version('polars')}, indent=1))" >"${WORK}/runtime.json" || true
aws s3 cp "${WORK}/runtime.json" "s3://${BUCKET}/${OUT}/runtime.json" --only-show-errors || true
COMPLETED=1
echo "[worker] $(date -u +%FT%TZ) complete"
