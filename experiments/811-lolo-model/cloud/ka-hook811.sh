#!/bin/bash
# 811 C1 hook known answers (EC2 only), v2 (v1's shasum calls were empty on Amazon Linux: its KA-off line was VOID).
#   KA-off: dumphook707.sh ITSELF (unchanged; its hard-coded repo path is provided by a symlink) and dumphook811.sh with
#           HOOK811_LAMBDA=0 on A/B/C's dump, same machine: outputs byte-identical (sha256 via python, both files must exist);
#   the lambda 0.25 run + 690g's validator; KA-same-construction / KA-scale on engine folds (kasame811.py).
set -uo pipefail
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 KASAME811_OUT=/work/kasame811
K=/work/hk; mkdir -p $K; O=experiments/811-lolo-model/work/hook811/ka2; mkdir -p $O
DUMP=experiments/811-lolo-model/work/hook811/dump_h1x2_cd4x15_lodo05hi6.npz
NAME=dump_h1x2_cd4x15_lodo05hi6.npz
echo '{}' > $K/identity.json
H() { uv run python -c 'import hashlib,sys; print(hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest())' "$1" 2>/dev/null || echo MISSING; }
mkdir -p ~/Code && ln -sfn /work/repo ~/Code/virtual-cell-challenge-2026
rc=0
bash experiments/707-strong-tier-direction/dumphook707.sh data/challenge $K $DUMP $K/identity.json $K/h707 > $O/h707.log 2>&1 || rc=1
HOOK811_LAMBDA=0 bash experiments/811-lolo-model/dumphook811.sh data/challenge $K $DUMP $K/identity.json $K/off > $O/off.log 2>&1 || rc=1
A=$(H $K/h707/$NAME); Bh=$(H $K/off/$NAME); Cs=$(H $K/off/swap707/$NAME)
echo "dumphook707 output      $A" | tee -a $O/summary.txt
echo "dumphook811 lambda 0     $Bh" | tee -a $O/summary.txt
echo "dumphook811's 707 step   $Cs" | tee -a $O/summary.txt
if [ "$A" != MISSING ] && [ "$A" = "$Bh" ] && [ "$A" = "$Cs" ]; then echo "KA-off (lambda 0 == dumphook707, byte for byte, same machine): PASS" | tee -a $O/summary.txt
else echo "KA-off: FAIL" | tee -a $O/summary.txt; rc=1; fi
s=$(date +%s); bash experiments/811-lolo-model/dumphook811.sh data/challenge $K $DUMP $K/identity.json $K/on > $O/on.log 2>&1 || rc=1; e=$(date +%s)
echo "lambda 0.25 run wall $((e-s)) s; output $(H $K/on/$NAME)" | tee -a $O/summary.txt
uv run python scripts/def_dumphook_check.py $DUMP $K/on/$NAME data/challenge > $O/check_on.txt 2>&1 && echo "690g validator on the lambda 0.25 output: PASS" | tee -a $O/summary.txt || { echo "690g validator: FAIL" | tee -a $O/summary.txt; rc=1; }
cp $K/on/*.json $K/on/swap707/*.json $O/ 2>/dev/null
for f in ABCp/k562_gwps ABCp/orion_hct116 ABCp/orion_hek293t L1/k562_gwps; do
  uv run python experiments/811-lolo-model/kasame811.py $f > $O/kasame-$(echo $f | tr / _).log 2>&1 && echo "KA-same $f: PASS" | tee -a $O/summary.txt || { echo "KA-same $f: FAIL" | tee -a $O/summary.txt; rc=1; }
done
cp /work/kasame811/*.json $O/ 2>/dev/null
cat $O/summary.txt
exit $rc
