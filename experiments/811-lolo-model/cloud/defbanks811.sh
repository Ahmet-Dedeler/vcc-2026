#!/bin/bash
# 811 X2 banks fix (EC2 only; coordinator's ruling SIXA 344c4c31, items 1-2). One job, in order, stopping on the first failure:
#   1. 801's kit811ext.py UNCHANGED with KIT811_PANELS=DEF (D/E/F's real panel; DEFo has the same targets): the pinned banks + D/E/F's
#      rows appended -> work/banks_ext (the code that built the rows 811 priced on DEFo / DEF).
#   2. defbanks811.py: banks_ext -> work/banks_def in hook811's manifest schema (byte-identical npz); swap it in as work/banks; re-pin
#      hook811.BANKS_MANIFEST_SHA (the ONLY code change; the diff must be exactly one line).
#   3. KA-same-construction: kasame811c0.py (unchanged) on DEFo x3, DEF x3, ABCp x3, L1/k562: the re-pinned hook's transform of
#      composite811's unswapped rows == composite811's C0 rows (|d| = 0), with the new banks.
#   4. KA-priced: composite811.py --mode day (unchanged) with the new banks reproduces the PRICED DEFo / DEF arrays (811's
#      run-comp811kit-DEFo / -DEF, every arm, every array) exactly.
#   5. A/B/C hook answers (content): dumphook811c0.sh at lambda 0.25 on the 10-20 dump 2dfe327d == the Mac's output (96204332...)
#      array by array; lambda 0 == the input dump's content. (The file-sha KA runs on the Mac, where the npz was written.)
# Outputs: work/defbanks/ (logs, diff, summary) and the new work/banks (collected).
set -uo pipefail
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
M=experiments/811-lolo-model; O=$M/work/defbanks; mkdir -p $O
log() { echo "[defbanks811 $(date -u +%FT%TZ)] $*" | tee -a $O/summary.txt; }
# 1
KIT811_PANELS=DEF uv run python experiments/801-lab-folds/kit/kit811ext.py > $O/kit811ext.log 2>&1 || { log "1 kit811ext FAILED"; tail -20 $O/kit811ext.log; exit 3; }
log "1 kit811ext rc 0 ($(grep -c 'old ' $O/kit811ext.log) banks; $(grep 'KA L1' $O/kit811ext.log))"
cp $M/work/banks_ext/manifest.json $O/banks_ext.manifest.json
# 2
KIT811_PANELS=DEF uv run python $M/defbanks811.py > $O/defbanks.log 2>&1 || { log "2 defbanks811 FAILED"; cat $O/defbanks.log; exit 4; }
mv $M/work/banks $M/work/banks.pre-def && mv $M/work/banks_def $M/work/banks
NEW=$(sha256sum $M/work/banks/manifest.json | cut -d' ' -f1)
cp $M/hook811.py $O/hook811.py.orig
sed -i "s/^BANKS_MANIFEST_SHA = \"[0-9a-f]*\"/BANKS_MANIFEST_SHA = \"$NEW\"/" $M/hook811.py
diff -u $O/hook811.py.orig $M/hook811.py > $O/hook811.diff
nd=$(grep -cE '^[-+][^-+]' $O/hook811.diff)
[ "$nd" = 2 ] || { log "2 hook811 diff is not exactly one line ($nd changed lines)"; exit 5; }
cp $M/hook811.py $O/hook811.py.pinned
log "2 new banks manifest $NEW; hook811 re-pinned, sha256 $(sha256sum $M/hook811.py | cut -d' ' -f1)"
# 3
for f in DEFo/k562_gwps DEFo/orion_hct116 DEFo/orion_hek293t DEF/k562_gwps DEF/orion_hct116 DEF/orion_hek293t \
         ABCp/k562_gwps ABCp/orion_hct116 ABCp/orion_hek293t L1/k562_gwps; do echo $f; done \
  | xargs -P 5 -I{} sh -c 'f={}; KASAME811_OUT=/work/kasame811 uv run python experiments/811-lolo-model/kasame811c0.py $f > '"$O"'/kasame-$(echo $f | tr / _).log 2>&1; echo "$f rc=$?"' \
  | sort > $O/kasame-rc.txt
cp /work/kasame811/c0-*.json $O/ 2>/dev/null
if grep -qv "rc=0" $O/kasame-rc.txt; then log "3 KA-same-C0 FAILED: $(grep -v 'rc=0' $O/kasame-rc.txt | tr '\n' ' ')"; exit 6; fi
log "3 KA-same-C0 PASS on $(wc -l < $O/kasame-rc.txt) folds (DEFo x3, DEF x3, ABCp x3, L1)"
# 4 (as run-kit811c.sh ran them: the kit's ka804 references in refcomp, --ka on the tag folds)
for T in DEFo DEF; do for h in k562_gwps orion_hct116 orion_hek293t; do
  mkdir -p $M/work/refcomp/${T}__$h
  cp experiments/801-lab-folds/work/ka-$(echo $T | tr A-Z a-z)/bulk/${T}__$h/705.npz experiments/801-lab-folds/work/ka-$(echo $T | tr A-Z a-z)/bulk/${T}__$h/707.npz $M/work/refcomp/${T}__$h/ || exit 12
done; done
for T in DEFo DEF; do for h in k562_gwps orion_hct116 orion_hek293t; do echo "$T/$h"; done; done \
  | xargs -P 6 -I{} sh -c 'f={}; COMPOSITE811_OUT=/work/compnew uv run python experiments/811-lolo-model/composite811.py $f --ka --mode day > '"$O"'/comp-$(echo $f | tr / _).log 2>&1; echo "$f rc=$?"' \
  | sort > $O/comp-rc.txt
if grep -qv "rc=0" $O/comp-rc.txt; then log "4 composite811 re-run FAILED: $(grep -v 'rc=0' $O/comp-rc.txt | tr '\n' ' ')"; exit 7; fi
uv run python - "$O" <<'PY' || { log "4 KA-priced FAILED (see $O/ka-priced.txt)"; exit 8; }
import sys, numpy as np
from pathlib import Path
O = Path(sys.argv[1]); out, bad, n = [], 0, 0
for T in ("DEFo", "DEF"):
    for h in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        new = Path(f"/work/compnew/bulk/{T}__{h}"); old = Path(f"experiments/811-lolo-model/work/run-comp811kit-{T}/bulk/{T}__{h}")
        arms = sorted(p.name for p in old.glob("*.npz"))
        assert arms and arms == sorted(p.name for p in new.glob("*.npz")), (T, h, arms)
        for a in arms:
            za, zb = np.load(old / a, allow_pickle=True), np.load(new / a, allow_pickle=True)
            assert sorted(za.files) == sorted(zb.files), (T, h, a)
            for k in za.files:
                x, y = za[k], zb[k]; n += 1
                same = x.shape == y.shape and (np.array_equal(x, y) if x.dtype.kind not in "fc" else np.array_equal(x, y, equal_nan=True))
                if not same:
                    bad += 1; out.append(f"DIFF {T}/{h} {a} {k}")
out.append(f"KA-priced: {n} arrays compared on DEFo x3 + DEF x3, {bad} differ")
(O / "ka-priced.txt").write_text("\n".join(out) + "\n"); print(out[-1])
sys.exit(1 if bad else 0)
PY
log "4 KA-priced PASS: $(tail -1 $O/ka-priced.txt)"
# 5
B=experiments/180-def-readiness/work/def-bundle; IN=$M/work/kamac/in/dump_h1x2_cd4x1_s0.25.npz; REF=$M/work/kamac/c0on/dump_h1x2_cd4x1_s0.25.npz
bash $M/dumphook811c0.sh $B /work/hk $IN /dev/null /work/hk25 > $O/hook25.log 2>&1 || { log "5 hook lambda 0.25 FAILED"; cat $O/hook25.log; exit 9; }
HOOK811_LAMBDA=0 bash $M/dumphook811c0.sh $B /work/hk $IN /dev/null /work/hk0 > $O/hook0.log 2>&1 || { log "5 hook lambda 0 FAILED"; cat $O/hook0.log; exit 10; }
uv run python - "$O" "$IN" "$REF" /work/hk25/dump_h1x2_cd4x1_s0.25.npz /work/hk0/dump_h1x2_cd4x1_s0.25.npz <<'PY' || { log "5 A/B/C hook content KA FAILED"; exit 11; }
import sys, hashlib, numpy as np
from pathlib import Path
O, IN, REF, H25, H0 = sys.argv[1:]
def eq(p, q):
    a, b = np.load(p, allow_pickle=True), np.load(q, allow_pickle=True)
    if sorted(a.files) != sorted(b.files): return False, "keys differ"
    d = [k for k in a.files if not (a[k].shape == b[k].shape and np.array_equal(a[k], b[k]))]
    return not d, f"differing keys {d}" if d else "all arrays equal"
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
r25, r0 = eq(H25, REF), eq(H0, IN)
lines = [f"lambda 0.25 vs the Mac's c0on output ({sha(REF)[:16]}): {r25[1]}; this file's sha {sha(H25)[:16]}",
         f"lambda 0 vs the input dump ({sha(IN)[:16]}): {r0[1]}"]
(Path(O) / "ka-abc-hook.txt").write_text("\n".join(lines) + "\n"); print("\n".join(lines))
sys.exit(0 if (r25[0] and r0[0]) else 1)
PY
log "5 A/B/C hook content PASS: $(tr '\n' ' ' < $O/ka-abc-hook.txt)"
log "ALL PASS"
