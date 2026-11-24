#!/bin/bash
# 811 X2 banks fix: the MAC file-sha KA (coordinator item 2, "lambda 0.25 on 2dfe327d -> 96204332; lambda 0 -> input content") with the
# re-pinned hook811 (c7cff26d) and the new banks (manifest 4cc6d82a). Local, ~5 s; run through gate602.
set -euo pipefail
cd "$(dirname "$0")/../.."
M=experiments/811-lolo-model; B=experiments/180-def-readiness/work/def-bundle; IN=$M/work/kamac/in/dump_h1x2_cd4x1_s0.25.npz
O=$M/work/kamac/def1022; mkdir -p $O
bash $M/dumphook811c0.sh $B $O/w $IN /dev/null $O/c0on > $O/c0on.log 2>&1
HOOK811_LAMBDA=0 bash $M/dumphook811c0.sh $B $O/w $IN /dev/null $O/c0off > $O/c0off.log 2>&1
on=$(shasum -a 256 $O/c0on/dump_h1x2_cd4x1_s0.25.npz | cut -d' ' -f1)
[ "$on" = 9620433220459420b02bec4eba4d7eba6cdb97b7f3a21ede122945461025df7a ] && r1=PASS || r1="FAIL ($on)"
uv run python -c "
import numpy as np, sys
a, b = np.load(sys.argv[1], allow_pickle=True), np.load(sys.argv[2], allow_pickle=True)
ok = sorted(a.files) == sorted(b.files) and all(a[k].shape == b[k].shape and np.array_equal(a[k], b[k]) for k in a.files)
print('PASS' if ok else 'FAIL')" $O/c0off/dump_h1x2_cd4x1_s0.25.npz $IN > $O/c0off.cmp
echo "$(date -u +%FT%TZ) Mac KA, hook811 $(shasum -a 256 $M/hook811.py | cut -c1-16), banks $(shasum -a 256 $M/work/banks/manifest.json | cut -c1-16): lambda 0.25 file sha == 96204332: $r1; lambda 0 content == input: $(cat $O/c0off.cmp)" | tee $O/summary.txt
