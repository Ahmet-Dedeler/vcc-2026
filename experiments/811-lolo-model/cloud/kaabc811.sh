#!/bin/bash
# 811 X2 banks fix, follow-up KA (EC2; the 07:42Z defbanks811 job's step 5 compared EC2's hook output with the MAC's file and differed in
# 'bulk'). Same platform, same dump (2dfe327d), lambda 0.25: (old) the pinned banks + hook811 as committed vs (new) defbanks811's banks
# (work/banks_def_ec2, manifest 4cc6d82a) + the re-pinned hook (work/defbanks/hook811.py.pinned). PASS iff old == new array for array.
# Also reported: max |d| of each against the Mac's 10-20 output (platform arithmetic), and lambda 0 == input for the new path.
set -uo pipefail
for d in experiments/*/; do mkdir -p "${d}work"; done
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
M=experiments/811-lolo-model; O=$M/work/kaabc; mkdir -p $O
B=experiments/180-def-readiness/work/def-bundle; IN=$M/work/kamac/in/dump_h1x2_cd4x1_s0.25.npz; REF=$M/work/kamac/c0on/dump_h1x2_cd4x1_s0.25.npz
bash $M/dumphook811c0.sh $B /work/hk $IN /dev/null /work/hkold > $O/hookold.log 2>&1 || { echo "old hook FAILED"; cat $O/hookold.log; exit 3; }
mv $M/work/banks $M/work/banks.pre-def && cp -R $M/work/banks_def_ec2 $M/work/banks && cp $M/work/defbanks/hook811.py.pinned $M/hook811.py
bash $M/dumphook811c0.sh $B /work/hk $IN /dev/null /work/hknew > $O/hooknew.log 2>&1 || { echo "new hook FAILED"; cat $O/hooknew.log; exit 4; }
HOOK811_LAMBDA=0 bash $M/dumphook811c0.sh $B /work/hk $IN /dev/null /work/hknew0 > $O/hooknew0.log 2>&1 || exit 5
uv run python - $O $IN $REF /work/hkold/dump_h1x2_cd4x1_s0.25.npz /work/hknew/dump_h1x2_cd4x1_s0.25.npz /work/hknew0/dump_h1x2_cd4x1_s0.25.npz <<'PY'
import sys, hashlib, numpy as np
from pathlib import Path
O, IN, REF, OLD, NEW, NEW0 = sys.argv[1:]
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
def cmp(p, q):
    a, b = np.load(p, allow_pickle=True), np.load(q, allow_pickle=True)
    assert sorted(a.files) == sorted(b.files)
    out = {}
    for k in a.files:
        x, y = a[k], b[k]
        if x.shape != y.shape: out[k] = "shape"; continue
        if np.array_equal(x, y): continue
        out[k] = float(np.max(np.abs(x.astype(float) - y.astype(float)))) if x.dtype.kind in "fiu" else "differs"
    return out
L = [f"files: old {sha(OLD)} new {sha(NEW)} mac-ref {sha(REF)} input {sha(IN)}",
     f"old vs new (same EC2 platform): {cmp(OLD, NEW) or 'all arrays equal'}",
     f"old vs mac-ref: {cmp(OLD, REF) or 'all arrays equal'}",
     f"new vs mac-ref: {cmp(NEW, REF) or 'all arrays equal'}",
     f"new lambda 0 vs input: {cmp(NEW0, IN) or 'all arrays equal'}"]
(Path(O) / "kaabc.txt").write_text("\n".join(L) + "\n"); print("\n".join(L))
ok = not cmp(OLD, NEW) and not cmp(NEW0, IN)
print("KA-ABC same-platform:", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
PY
