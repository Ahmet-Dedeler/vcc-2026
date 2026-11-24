#!/usr/bin/env bash
# labworker.sh — the remote half of labcloud.sh. Runs ON a lab EC2 instance (never by hand).
#
# CONTRACT (all under s3://$BUCKET/$PREFIX/ unless noted)
#   in   job.json        {name, sync:[relpath...], collect:[relpath...], parallel:N, venv:bool}
#        code.tar.gz     repo code (labcloud.sh's tarball), unpacked at /work/repo
#        jobs.txt        one shell command per line, run in /work/repo (bash -c), N at a time
#        s3://$BUCKET/input/lab/repo/<relpath>   the data mirror; each job.json `sync` entry is pulled to /work/repo/<relpath>
#   out  s3://$BUCKET/$OUT/STATUS              DONE 0 | FAILED <rc> (completed=0|1)
#        s3://$BUCKET/$OUT/worker.log          this transcript
#        s3://$BUCKET/$OUT/jobs/<i>.log        each job's output + /usr/bin/time -v (peak RSS, wall)
#        s3://$BUCKET/$OUT/jobs/summary.tsv    i, rc, wall_s, maxrss_kb, command
#        s3://$BUCKET/$OUT/repo/<relpath>      every `collect` path after the jobs ran (dirs synced, globs expanded)
#
# SELF-TERMINATION (three guards, as scripts/cloud_score_worker.sh): launched with shutdown behaviour = terminate; user-data
# armed `shutdown -h +MAX` before this script started; this script's EXIT trap powers off on every path.
set -Eeuo pipefail
: "${BUCKET:?}" "${PREFIX:?}" "${OUT:?}" "${AWS_REGION:?}"
W=/work; REPO=$W/repo; LOG=$W/worker.log; STATUS=$W/STATUS
mkdir -p "$W" "$REPO" "$W/jobs"
exec > >(tee -a "$LOG") 2>&1
COMPLETED=0
push() {
  aws s3 cp "$LOG" "s3://${BUCKET}/${OUT}/worker.log" --only-show-errors || true
  [[ -s "$STATUS" ]] && aws s3 cp "$STATUS" "s3://${BUCKET}/${OUT}/STATUS" --only-show-errors || true
}
finish() {
  rc=$?; trap - EXIT INT TERM
  if [[ "$COMPLETED" -eq 1 && "$rc" -eq 0 ]]; then printf 'DONE 0\n' >"$STATUS"; else printf 'FAILED %s (completed=%s)\n' "$rc" "$COMPLETED" >"$STATUS"; fi
  [[ -d "$W/jobs" ]] && aws s3 sync "$W/jobs/" "s3://${BUCKET}/${OUT}/jobs/" --only-show-errors || true
  push
  shutdown -h now || poweroff || true
  exit "$rc"
}
trap finish EXIT
trap 'echo "[labworker] caught a signal: this run did NOT complete" >&2; exit 143' INT TERM
say() { echo "[labworker $(date -u +%FT%TZ)] $*"; }
say "start on $(hostname): $(nproc) vCPU, $(awk '/MemTotal/{printf "%.0f GiB", $2/1048576}' /proc/meminfo), disk $(df -h / | awk 'NR==2{print $4}') free"

dnf install -y -q gzip tar time jq git >/dev/null 2>&1 || dnf install -y gzip tar time jq git
export HOME=/root
curl -LsSf https://astral.sh/uv/install.sh | sh >/dev/null 2>&1
export PATH="/root/.local/bin:$PATH" UV_CACHE_DIR="$W/uv-cache" UV_LINK_MODE=copy
aws s3 cp "s3://${BUCKET}/${PREFIX}/job.json" "$W/job.json" --only-show-errors
aws s3 cp "s3://${BUCKET}/${PREFIX}/jobs.txt" "$W/jobs.txt" --only-show-errors
aws s3 cp "s3://${BUCKET}/${PREFIX}/code.tar.gz" "$W/code.tar.gz" --only-show-errors
cat "$W/job.json"
tar -xzf "$W/code.tar.gz" -C "$REPO"
say "code unpacked: $(find "$REPO" -type f | wc -l) files"

# ---- data mirror
while IFS= read -r rel; do
  [[ -n "$rel" ]] || continue
  src="s3://${BUCKET}/input/lab/repo/${rel}"
  if aws s3api head-object --bucket "$BUCKET" --key "input/lab/repo/${rel}" >/dev/null 2>&1; then
    mkdir -p "$(dirname "$REPO/$rel")"; aws s3 cp "$src" "$REPO/$rel" --only-show-errors          # a file
  else
    aws s3 sync "${src}/" "$REPO/${rel}/" --only-show-errors                                     # a directory prefix
  fi
  say "synced ${rel} ($(du -sh "$REPO/$rel" 2>/dev/null | cut -f1))"
done < <(jq -r '.sync[]?' "$W/job.json")

# ---- python env (the repo's lock, python 3.12)
if [[ "$(jq -r '.venv // true' "$W/job.json")" == "true" ]]; then
  ( cd "$REPO" && uv python install 3.12 >/dev/null 2>&1; uv sync --frozen --python 3.12 2>&1 | tail -3 )
  say "venv ready: $(cd "$REPO" && uv run python -c 'import numpy, scipy, anndata; print(numpy.__version__, scipy.__version__, anndata.__version__)')"
fi

# ---- jobs, N at a time; each job's rc, wall and peak RSS recorded
P="$(jq -r '.parallel // 1' "$W/job.json")"
cat > "$W/runone.sh" <<'EOS'
#!/usr/bin/env bash
i="$1"; shift; cmd="$*"
cd /work/repo
s=$(date +%s)
/usr/bin/time -v -o "/work/jobs/$i.time" bash -c "$cmd" > "/work/jobs/$i.log" 2>&1; rc=$?
e=$(date +%s); rss=$(awk -F': ' '/Maximum resident set size/{print $2}' "/work/jobs/$i.time" 2>/dev/null)
cat "/work/jobs/$i.time" >> "/work/jobs/$i.log" 2>/dev/null; rm -f "/work/jobs/$i.time"
printf '%s\t%s\t%s\t%s\t%s\n' "$i" "$rc" "$((e-s))" "${rss:-NA}" "$cmd" >> /work/jobs/summary.tsv
echo "[job $i] rc=$rc wall=$((e-s))s maxrss_kb=${rss:-NA}"
exit 0
EOS
chmod +x "$W/runone.sh"
: > "$W/jobs/summary.tsv"
awk 'NF{printf "%d\t%s\n", NR, $0}' "$W/jobs.txt" > "$W/jobs.idx"
say "running $(wc -l < "$W/jobs.idx") jobs, ${P} at a time"
while IFS=$'\t' read -r i cmd; do printf '%s\0%s\0' "$i" "$cmd"; done < "$W/jobs.idx" | xargs -0 -n 2 -P "$P" "$W/runone.sh"
nfail=$(awk -F'\t' '$2!=0' "$W/jobs/summary.tsv" | wc -l)
say "jobs done: $(wc -l < "$W/jobs/summary.tsv") ran, ${nfail} failed"

# ---- collect
while IFS= read -r rel; do
  [[ -n "$rel" ]] || continue
  for p in $REPO/$rel; do
    [[ -e "$p" ]] || { say "collect: nothing at ${rel}"; continue; }
    r="${p#$REPO/}"
    if [[ -d "$p" ]]; then aws s3 sync "$p/" "s3://${BUCKET}/${OUT}/repo/${r}/" --only-show-errors
    else aws s3 cp "$p" "s3://${BUCKET}/${OUT}/repo/${r}" --only-show-errors; fi
    say "collected ${r}"
  done
done < <(jq -r '.collect[]?' "$W/job.json")
[[ "$nfail" -eq 0 ]] || { say "some jobs failed"; exit 3; }
COMPLETED=1
