#!/usr/bin/env bash
# 247 item C, the byte-identity row (240 row 1): the driver's defaults, nothing exported, through
# the modified builder on the A/B/C bundle must rebuild the shipped `h1union` byte-for-byte:
#   .h5ad  8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6  12,648,807,800 B
#   .vcc   f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e   2,824,785,920 B
# Build (22 min, 20 GiB) then `vcc prep` (26 s, 19.7 GiB), strictly serial, each gated; then hash,
# verify, and DELETE the two artifacts (they are byte duplicates of submissions/222-h1union-g028.*).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/ka1"; mkdir -p "$W"; LOG="$W/verify247.log"
NAME="247-ka-h1union"
WANT_H5AD="8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6"; WANT_H5AD_BYTES=12648807800
WANT_VCC="f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e";  WANT_VCC_BYTES=2824785920
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
gate_ok(){
  local pct fi others
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  log "gate: memory_pressure free ${pct}%  free+inactive ${fi} GiB  other builders ${others}"
  [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=30)}' && [ "$others" = "0" ]
}
wait_gate(){ local n=0; until gate_ok; do n=$((n+1)); [ $n -ge 240 ] && { log "GATE never opened; giving up"; exit 3; }; sleep 30; done; }

log "ARMED pid=$$"
wait_gate
log "BUILD start (driver defaults, nothing exported; stages build,postbuild)"
env -u VCC_ALPHA -u VCC_CD4_TARGET_F -u VCC_H1_UNION_G -u VCC_UNION_SOURCES -u VCC_CD4_SIZING -u VCC_CD4_TARGET_G \
  uv run python scripts/def_pipeline.py --bundle data/challenge --work "$W" --name "$NAME" \
  --stages build,postbuild --whitener-out data/donors/whitener_union_pruned.npz --whitener-id union_pruned \
  --cd4-out data/external/cd4t/cloud-compact/cd4_compact --frozen-inputs scripts/frozen_inputs.json \
  > "$W/driver.log" 2>&1
RC=$?; log "BUILD rc=$RC"
[ $RC -eq 0 ] || { log "BUILD FAILED; see $W/driver.log"; exit 1; }
H5="$W/$NAME.h5ad"
GOT_H5="$(shasum -a 256 "$H5" | cut -d' ' -f1)"; H5B="$(stat -f%z "$H5")"
log "H5AD sha256 $GOT_H5 bytes $H5B  (want $WANT_H5AD $WANT_H5AD_BYTES)  $([ "$GOT_H5" = "$WANT_H5AD" ] && [ "$H5B" = "$WANT_H5AD_BYTES" ] && echo MATCH || echo MISMATCH)"
wait_gate
log "PREP start"
/usr/bin/time -l vcc prep "$H5" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$W/$NAME.vcc" > "$W/prep.log" 2>&1
RC=$?; log "PREP rc=$RC peak RSS $(awk '/maximum resident set size/{printf "%.2f", $1/2^30}' "$W/prep.log") GiB"
GOT_VCC="$(shasum -a 256 "$W/$NAME.vcc" | cut -d' ' -f1)"; VB="$(stat -f%z "$W/$NAME.vcc")"
log "VCC  sha256 $GOT_VCC bytes $VB  (want $WANT_VCC $WANT_VCC_BYTES)  $([ "$GOT_VCC" = "$WANT_VCC" ] && [ "$VB" = "$WANT_VCC_BYTES" ] && echo MATCH || echo MISMATCH)"
uv run python - "$W" "$NAME" "$GOT_H5" "$H5B" "$GOT_VCC" "$VB" <<'EOF'
import json, sys
w, name, h5, h5b, vcc, vb = sys.argv[1:]
timings = json.load(open(f"{w}/timings.json"))
meta = json.load(open(f"{w}/{name}.h5ad.meta.json"))
out = {"h5ad_sha256": h5, "h5ad_bytes": int(h5b), "vcc_sha256": vcc, "vcc_bytes": int(vb),
       "h5ad_match": h5 == "8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6" and int(h5b) == 12648807800,
       "vcc_match": vcc == "f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e" and int(vb) == 2824785920,
       "build_seconds": timings["build"]["seconds"], "build_peak_rss_gib": timings["build"]["peak_rss_gib"],
       "postbuild": timings.get("postbuild"),
       "config_cd4": {k: meta["config"].get(k) for k in ("cd4_target_f", "cd4_sizing", "cd4_target_g", "h1_union_g", "alpha", "emitter")},
       "effects_cd4_fraction": meta["effects"]["cd4"]["fraction"],
       "output_sha256": meta["output"]["sha256"],
       "builder_sha256": meta["hashes"]["implementation"]["scripts/make_submission_pm_cd4.py"]}
json.dump(out, open(f"{w}/verify247.json", "w"), indent=1)
print(json.dumps(out, indent=1))
EOF
# The two artifacts are byte duplicates of submissions/222-h1union-g028.* (if MATCH); reclaim 15 GB.
if [ "$GOT_H5" = "$WANT_H5AD" ] && [ "$GOT_VCC" = "$WANT_VCC" ]; then
  rm -f "$H5" "$W/$NAME.vcc"; log "DELETED $H5 and $W/$NAME.vcc (byte duplicates of the shipped package; meta and logs kept)"
else
  log "KEPT the artifacts for inspection (a MISMATCH is a finding)"
fi
log "DONE"
