#!/usr/bin/env bash
# 805 stage 0, one condition, on a lab EC2 instance (labworker runs it from /work/repo). PREREGISTRATION-stage0.md.
#   bash experiments/805-cd4-flex-surface/s0job805.sh Rest
#   bash experiments/805-cd4-flex-surface/s0job805.sh Stim48hr Rest    # runs only if Rest's gate verdict is FAIL
# Steps: NVMe scratch; KA (h1-test sha, ref614 DE sha, blob statistic) - abort on a fail; download the four donor files
# (public bucket, anonymous, same region); select; extract x4 in parallel; merge (thinning, truths); ref614 on P and S;
# the gate. Small results -> experiments/805-cd4-flex-surface/work/s0_<cond>/ (collected); big files -> s3 $OUT/big/<cond>/.
set -Eeuo pipefail
COND="$1"; IFFAIL="${2:-}"
X=experiments/805-cd4-flex-surface
RES="$X/work/s0_${COND}"
mkdir -p "$RES"
say() { echo "[s0job805 $(date -u +%FT%TZ)] $*"; }
H1_SHA=06c6baea8c62ef3d38cacdbf58e104908eb93f0ef3b25aca56b86d28240f9692
push() { aws s3 sync "$RES/" "s3://${BUCKET}/${OUT}/early/s0_${COND}/" --only-show-errors || true; }

if [[ -n "$IFFAIL" ]]; then
  v="$(jq -r '.gate.verdict' "$X/work/s0_${IFFAIL}/gate805.json")"
  say "previous condition ${IFFAIL}: verdict ${v}"
  [[ "$v" == "FAIL" ]] || { say "not running ${COND} (the prereg runs it only after a FAIL)"; echo skipped > "$RES/SKIPPED"; exit 0; }
fi

# ---- scratch: the instance-store NVMe if there is one
DEV="$(lsblk -dpno NAME,MODEL | awk '/Instance Storage/{print $1; exit}')"
if [[ -n "$DEV" ]] && ! mountpoint -q /nv; then
  mkfs.xfs -f -q "$DEV" && mkdir -p /nv && mount "$DEV" /nv
fi
mountpoint -q /nv && NV=/nv || NV=/work/nv
mkdir -p "$NV/src_${COND}" "$NV/out_${COND}" "$NV/ka"
say "scratch $NV ($(df -h "$NV" | awk 'NR==2{print $4}') free), $(nproc) vCPU, $(free -g | awk '/Mem/{print $2}') GiB"
aws configure set default.s3.max_concurrent_requests 32
aws configure set default.s3.multipart_chunksize 64MB

# ---- downloads start now (in the background), the KA runs meanwhile
for d in D1 D2 D3 D4; do
  ( aws s3 cp --no-sign-request --region us-east-1 --only-show-errors \
      "s3://genome-scale-tcell-perturb-seq/marson2025_data/${d}_${COND}.assigned_guide.h5ad" "$NV/src_${COND}/" \
    && echo "$d ok" > "$NV/src_${COND}/${d}.done" ) &
done

if [[ ! -f "$X/work/s0_ka/ka805.json" ]]; then
  mkdir -p "$X/work/s0_ka"
  T805_MOCK=1 uv run python "$X/test805_mock.py" "$NV/mock" > "$X/work/s0_ka/mock.log" 2>&1 \
    || { say "MOCK TEST FAILED"; tail -30 "$X/work/s0_ka/mock.log"; aws s3 sync "$X/work/s0_ka/" "s3://${BUCKET}/${OUT}/early/s0_ka/" --only-show-errors; kill $(jobs -p) 2>/dev/null; exit 8; }
  say "mock: $(tail -1 "$X/work/s0_ka/mock.log")"
  echo "${H1_SHA}  data/harness/h1-test.h5ad" | sha256sum --check - > "$X/work/s0_ka/h1_sha.txt"
  uv run python experiments/614-k562-real-emitter/ref614.py data/harness/h1-test.h5ad \
      experiments/614-k562-real-emitter/work/targets_h1.txt "$NV/ka/ref_h1" > "$X/work/s0_ka/ref_h1.log" 2>&1
  cp "$NV/ka/ref_h1/reference.json" "$X/work/s0_ka/reference_h1.json"
  mkdir -p "$X/work/s0_ka/ref_h1" && cp "$NV/ka/ref_h1/reference.json" "$NV/ka/ref_h1/reference_de.parquet" "$X/work/s0_ka/ref_h1/"
  uv run python "$X/s0stats805.py" ka "$X/work/s0_ka/ref_h1" > "$X/work/s0_ka/ka.log" 2>&1 \
    || { say "KNOWN ANSWER FAILED"; cat "$X/work/s0_ka/ka.log"; aws s3 sync "$X/work/s0_ka/" "s3://${BUCKET}/${OUT}/early/s0_ka/" --only-show-errors; kill $(jobs -p) 2>/dev/null; exit 9; }
  aws s3 sync "$X/work/s0_ka/" "s3://${BUCKET}/${OUT}/early/s0_ka/" --only-show-errors || true
  say "KA PASS: $(cat "$X/work/s0_ka/ka.log" | tail -1 | cut -c1-300)"
fi

wait
for d in D1 D2 D3 D4; do [[ -f "$NV/src_${COND}/${d}.done" ]] || { say "download of $d failed"; exit 4; }; done
say "downloads done: $(du -sh "$NV/src_${COND}" | cut -f1)"
ls -l "$NV/src_${COND}" > "$RES/src_listing.txt"

S0="uv run python $X/s0build805.py"
/usr/bin/time -v $S0 select "$COND" "$NV/src_${COND}" "$NV/out_${COND}" > "$RES/select.log" 2>&1
cp "$NV/out_${COND}/select.json" "$RES/"
push
for d in D1 D2 D3 D4; do
  /usr/bin/time -v $S0 extract "$COND" "$NV/src_${COND}" "$NV/out_${COND}" "$d" > "$RES/extract_${d}.log" 2>&1 &
done
wait
for d in D1 D2 D3 D4; do [[ -f "$NV/out_${COND}/X_${d}.npz" ]] || { say "extract $d failed"; tail -20 "$RES/extract_${d}.log"; exit 5; }; done
rm -f "$NV/src_${COND}"/*.h5ad
/usr/bin/time -v $S0 merge "$COND" "$NV/src_${COND}" "$NV/out_${COND}" > "$RES/merge.log" 2>&1
cp "$NV/out_${COND}/select.json" "$NV/out_${COND}/targets_P.txt" "$NV/out_${COND}/targets_S.txt" "$RES/"
push

# ---- the truths at a repo path (s0stats reads them relative to the repo), refs on scratch
mkdir -p "$X/work/truth_${COND}"
for s in P S; do ln -sf "$NV/out_${COND}/truth_${s}.h5ad" "$X/work/truth_${COND}/truth_${s}.h5ad"; done
for s in P S; do
  /usr/bin/time -v uv run python experiments/614-k562-real-emitter/ref614.py "$NV/out_${COND}/truth_${s}.h5ad" \
     "$NV/out_${COND}/targets_${s}.txt" "$NV/out_${COND}/ref_${s}" > "$RES/ref_${s}.log" 2>&1 &
done
wait
for s in P S; do
  [[ -f "$NV/out_${COND}/ref_${s}/reference.json" ]] || { say "ref $s failed"; tail -20 "$RES/ref_${s}.log"; exit 6; }
  mkdir -p "$X/work/truth_${COND}/ref_${s}"
  for f in reference.json reference.npz; do ln -sf "$NV/out_${COND}/ref_${s}/$f" "$X/work/truth_${COND}/ref_${s}/$f"; done
  cp "$NV/out_${COND}/ref_${s}/reference.json" "$RES/reference_${s}.json"
done
uv run python "$X/s0stats805.py" gate "$COND" "$X/work/truth_${COND}" > "$RES/gate.log" 2>&1 || { tail -30 "$RES/gate.log"; exit 7; }
cp "$X/work/truth_${COND}/gate805.json" "$RES/"
push
say "GATE ${COND}: $(jq -c '.gate' "$RES/gate805.json")"

# ---- big files kept on S3 (sha256 recorded)
( cd "$NV/out_${COND}" && sha256sum truth_P.h5ad truth_S.h5ad ref_P/* ref_S/* > "$OLDPWD/$RES/SHA256SUMS.big.txt" )
aws s3 cp "$NV/out_${COND}/" "s3://${BUCKET}/${OUT}/big/${COND}/" --recursive --exclude 'X_*.npz' --exclude 'rows_*' --only-show-errors
push
say "done ${COND}"
