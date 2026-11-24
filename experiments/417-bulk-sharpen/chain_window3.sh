#!/usr/bin/env bash
# 417 window 3 (PI, 2026-09-23 01:25Z): after 417h's build restores the builder,
#  1. apply patch_lodo05_417.py (VCC_PERT_SCALE=lodo415_lo0.5_hi3), run the unit tests + builder suite (56/0/1), commit;
#  2. build 417i = 417a + key lodo415_lo0.5_hi3                 (one variable vs 417a: pert_scale + sha);
#  3. dump the stack bulk with that key (H1 x2 + CD4 x1.5 + lodo0.5);
#  4. build 417j = LOAD(stack+lodo0.5) + 417h's re-level        (one variable vs 417h: bulk_load + sha).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
E="$REPO/experiments/417-bulk-sharpen"; W="$E/work"; LOG="$W/chain_window3.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
fail(){ log "STOP: $*"; git checkout -- scripts/make_submission_pm_cd4.py 2>/dev/null; exit 1; }
log "WAITING for 417h's builder restore (pid $$)"
n=0; until grep -qE "builder after restore|stopping|GATE never" "$W/build_417h/run417h.log" 2>/dev/null || grep -qE "NOT launched" "$W/chain417h.log" 2>/dev/null; do n=$((n+1)); [ $n -ge 1440 ] && { log "417h never restored"; exit 1; }; sleep 30; done
sleep 60
B0="fc30e8b0fd884867aa65299ff63b236ceb2f251ea85990a1eeef6a509b6e0eb4"
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$B0" ] || { log "builder is not fc30e8b0 -- not patching"; exit 1; }
[ -z "$(git status --short scripts/make_submission_pm_cd4.py)" ] || { log "builder dirty -- not patching"; exit 1; }
uv run python "$E/patch_lodo05_417.py" --write > "$W/patch_lodo05.apply.log" 2>&1 || fail "patch did not apply"
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; log "patched builder $B"
for t in "$E/test_sharpen417.py" "$E/test_next417.py" "$E/test_bulkio417.py" "$E/test_carrier417.py" experiments/410-pertscale-revisit/test_ballast410.py experiments/410-pertscale-revisit/test_pert_scale.py; do
  uv run python "$t" > "$W/w3_$(basename $t .py).log" 2>&1 || fail "$t FAILED"; log "$(tail -1 "$W/w3_$(basename $t .py).log" | cut -c1-140)"; done
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > "$W/builder-tests-417w3.log" 2>&1
grep -q "56 passed, 0 failed, 1 skipped" "$W/builder-tests-417w3.log" || fail "builder suite: $(grep -E 'passed.*failed' "$W/builder-tests-417w3.log" | tail -1)"
git add scripts/make_submission_pm_cd4.py "$W/builder-tests-417w3.log" && git commit -q -m "Builder window 3 (417): VCC_PERT_SCALE=lodo415_lo0.5_hi3 (420's LODO key at clip(0.5,3), 8c82bef3)

Builder fc30e8b0 -> ${B:0:8}; unset knobs bit-identical (56 passed / 0 failed / 1 skipped); 417 unit tests PASS.
Applied by experiments/417-bulk-sharpen/chain_window3.sh after 417h's build restored the builder.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" || fail "commit failed"
log "committed builder $(git log -1 --format=%h) ($B)"
# 2. 417i = 417a + lodo0.5 key
python3 - "$E" "$B" <<'PY'
import sys
from pathlib import Path
E, B = Path(sys.argv[1]), sys.argv[2]
t = (E / "run417a_build.sh").read_text()
t = t.replace("# 417a -- BUILD 410b + 169 magnitude sharpening (q 1.25) on the BALLAST BULK render only (VCC_BULK_SHARPEN_Q=1.25; 414 section 7).",
              "# 417i -- BUILD 417a + the key lodo415_lo0.5_hi3 (420's best key variant; one variable vs 417a).")
t = t.replace('W="$HERE/work/build_a"; mkdir -p "$W"; LOG="$W/run417a.log"', 'W="$HERE/work/build_i"; mkdir -p "$W"; LOG="$W/run417i.log"')
t = t.replace('NAME="the-pigeon-passed-peer-review"', 'NAME="the-knockdown-knocked-back"')
t = t.replace("417a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125", "417i-cap15-partial50-clean4q80-lodo05-ballast-sharpen125")
t = t.replace('REF_META="$REPO/submissions/410b-cap15-partial50-clean4q80-pertscale-ballast.h5ad.meta.json"', 'REF_META="$REPO/submissions/417a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125.h5ad.meta.json"')
t = t.replace("export VCC_PERT_SCALE=wshare_lo0.3_hi3  # 410b key (unchanged)", "export VCC_PERT_SCALE=lodo415_lo0.5_hi3  # THE ARM (417i): 420's LODO key at clip(0.5, 3)")
t = t.replace("export VCC_BULK_SHARPEN_Q=1.25           # THE ARM (417a): sharpening on the ballast bulk render only", "export VCC_BULK_SHARPEN_Q=1.25           # 417a (unchanged)")
t = t.replace("VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN\n", "VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN VCC_CARRIER_GAMMA VCC_BULK_LOAD VCC_BULK_CTX_SCALE VCC_H1_UNION_G VCC_CD4_TARGET_F\n", 1)
t = t.replace('[ "$B" = "400ef6e307896deb8c0adf694547c756000651f77ffb923ab450ef2f01389e31" ] || { log "BUILDER hash $B != 400ef6e3... (417a) -- stopping"; exit 4; }',
              f'[ "$B" = "{B}" ] || {{ log "BUILDER hash $B != {B[:8]}... (417i) -- stopping"; exit 4; }}')
t = t.replace("[417a]", "[417i]")
t = t.replace('ok = set(d)=={"bulk_sharpen_q"} and d["bulk_sharpen_q"]==("<absent>",1.25)',
              'ok = set(d)=={"pert_scale","pert_scale_sha256"} and d["pert_scale"]==("wshare_lo0.3_hi3","lodo415_lo0.5_hi3")')
t = t.replace('print("[417i] PASS: exactly bulk_sharpen_q" if ok else f"[417i] FAIL: got {d}")', 'print("[417i] PASS: exactly pert_scale + sha" if ok else f"[417i] FAIL: got {d}")')
t = t.replace("CONFIG DIFF is not exactly bulk_sharpen_q", "CONFIG DIFF is not exactly pert_scale + sha")
t = t.replace('log "cliffcheck (414): render-space energy vs 410b must be DOWN, about 0.96 per context"', 'log "cliffcheck (414): 417a with the lodo0.5 key -- energy vs 410b near 417a (0.96)"')
t = t.replace('a=r["A_vs_410b"]; ok=all(0.93<=a[c]<=0.995 for c in "ABC")', 'a=r["A_vs_410b"]; ok=all(0.92<=a[c]<=1.00 for c in "ABC")')
(E / "run417i_build.sh").write_text(t)
PY
bash -n "$E/run417i_build.sh" || fail "run417i_build.sh syntax"
git add "$E/run417i_build.sh" && git commit -q -m "417i build script (417a + lodo415_lo0.5_hi3 key), pinned to builder ${B:0:8}

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
log "2. building 417i"
bash "$E/run417i_build.sh" >> "$LOG" 2>&1
log "417i: $(grep -E '^\[417i\] (PASS|FAIL|CLIFF)' "$W/build_i/run417i.log" | tr '\n' ' ')"
# 3. dump the stack with the lodo0.5 key
export VCC_NORMREL_EXP=-1.5 VCC_CD4_PARTIAL_G=0.5 VCC_NORMREL_MODE=cap VCC_DE_SIGN_PRIOR=clean4q80 VCC_PERT_SCALE=lodo415_lo0.5_hi3 VCC_PERT_SCALE_BALLAST=1 VCC_BULK_SHARPEN_Q=1.25 VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.51
unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN VCC_BULK_CTX_SCALE VCC_BULK_LOAD VCC_CARRIER_GAMMA
log "3. DUMP stack + lodo0.5 key"
bash "$E/run_dump417.sh" h1x2_cd4x15_lodo05 "$W/dump_h1x2_cd4x15_lodo05.npz" "$B" >> "$LOG" 2>&1
[ -s "$W/dump_h1x2_cd4x15_lodo05.npz" ] || fail "stack+lodo0.5 dump failed"
unset VCC_H1_UNION_G VCC_CD4_TARGET_F VCC_BULK_DUMP
export VCC_PERT_SCALE=wshare_lo0.3_hi3
log "   dump sha $(shasum -a 256 "$W/dump_h1x2_cd4x15_lodo05.npz" | cut -d' ' -f1)"
# 4. 417j = 417h with the stack+lodo0.5 dump (one variable vs 417h: bulk_load + sha)
log "4. building 417j"
bash "$E/run_load_arm417.sh" 417j 1.25 "$W/dump_h1x2_cd4x15_lodo05.npz" "$REPO/submissions/417h-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-ctxlevel.h5ad.meta.json" \
  the-reference-cells-filed-an-appeal 417j-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stacklodo05load-ctxlevel "$B" "A=1.0088,B=1.0500,C=1.0736" >> "$LOG" 2>&1 || true
log "417j: $(grep -E '^\[417j\] (PASS|FAIL|CLIFF)' "$W/build_417j/run417j.log" | tr '\n' ' ')"
log "END"
