#!/usr/bin/env bash
# 395 -- finish the build after run395.sh stopped on an over-strict config-diff check (the diff was exactly the THREE keys a
# one-knob change produces: de_sign_prior clean4->clean4q80, de_sign_protect_q absent->0.8, de_sign_emitter_sha256 re-pin).
# Runs the remaining steps verbatim (prep, meta copy, novelty, sha) and appends the same log lines the waiter looks for.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; cd "$REPO"
W="$HERE/work/build"; LOG="$W/run395.log"; BUNDLE="$REPO/data/challenge"
NAME="the-thermostat-has-a-union-rep"; PKG="$REPO/submissions/395-cap15-partial50-designclean4q80.vcc"
META="$W/$NAME.h5ad.meta.json"; REF_META="$REPO/submissions/392-cap15-partial50-designclean4.h5ad.meta.json"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
free_gib(){ df -k "$REPO" | awk 'NR==2{printf "%.1f", $4/1048576}'; }
log "FINISH395: re-checking the config diff with the correct expectation (three keys)"
uv run python - "$META" "$REF_META" <<'PY' | tee -a "$LOG"
import json,sys
built=json.load(open(sys.argv[1]))["config"]; ref=json.load(open(sys.argv[2]))["config"]
d={k:(ref.get(k,"<absent>"),built.get(k,"<absent>")) for k in sorted(set(ref)|set(built)) if ref.get(k,"<absent>")!=built.get(k,"<absent>")}
want={"de_sign_prior","de_sign_protect_q","de_sign_emitter_sha256"}
ok = (set(d)==want and d["de_sign_prior"]==("clean4","clean4q80") and d["de_sign_protect_q"]==("<absent>",0.8)
      and d["de_sign_emitter_sha256"]==("a12cf8f9f7231f4b23317fa1288d3f522c4c429eb634fda2bc84b6ae4fb575de","6ff0b866b9b0b88dd4ce0bcc60a49f38b4451422ac6a9f8226c615f0215637ba"))
print("[395] PASS: exactly the knob value, protect_q and the emitter re-pin" if ok else f"[395] FAIL: got {d}")
sys.exit(0 if ok else 1)
PY
[ ${PIPESTATUS[0]} -eq 0 ] || { log "config diff still wrong -- stopping"; exit 1; }
log "builder is restored: $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
log "INVOCATION 2: vcc prep -> $PKG"
/usr/bin/time -l vcc prep "$W/$NAME.h5ad" -g "$BUNDLE/gene_names.csv" --perts "$BUNDLE/pert_counts.csv" --contexts A,B,C -o "$PKG" > "$W/prep.log" 2>&1
RC=$?; log "prep rc=$RC  disk $(free_gib) GiB"
tail -20 "$W/prep.log" | tee -a "$LOG"
[ $RC -eq 0 ] && [ -s "$PKG" ] || { log "PREP FAILED -- stopping"; exit 1; }
log "vcc: $(ls -la "$PKG" | awk '{print $5, $9}')"
cp "$META" "$REPO/submissions/395-cap15-partial50-designclean4q80.h5ad.meta.json"
log "meta copied beside the package: $(shasum -a 256 "$REPO/submissions/395-cap15-partial50-designclean4q80.h5ad.meta.json" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$PKG" 2>&1 | tee -a "$LOG"
log "novelty rc=${PIPESTATUS[0]}"
log "package sha256: $(shasum -a 256 "$PKG" | cut -d' ' -f1)"
log "git status scripts/ src/: $(git status --short scripts/ src/ | tr '\n' ' ')"
log "DONE"
