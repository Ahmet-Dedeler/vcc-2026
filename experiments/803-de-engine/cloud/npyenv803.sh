#!/bin/bash
# 803: run a command with some of NumPy's runtime-dispatched SIMD paths switched off (platform-difference diagnosis).
#   bash npyenv803.sh <none|avx512|avx2> <cmd ...>
#   avx512: every dispatched AVX512* feature off (NumPy's SVML float64 log1p / expm1 / exp / log live there on x86-64 Linux)
#   avx2:   avx512 + AVX2 / FMA3 / AVX / F16C off (SSE baseline only)
set -e
m=$1; shift
if [ "$m" != none ]; then
  feats=$(uv run python -c "
import numpy as np, sys
from numpy._core._multiarray_umath import __cpu_dispatch__ as D, __cpu_features__ as F
m = sys.argv[1]
off = [f for f in D if f.startswith('AVX512') or (m == 'avx2' and f in ('AVX2', 'FMA3', 'AVX', 'F16C'))]
print(' '.join(off))" "$m")
  export NPY_DISABLE_CPU_FEATURES="$feats"
  echo "[npyenv803] $m: NPY_DISABLE_CPU_FEATURES=$feats"
fi
uv run python -c "import numpy as np; from numpy._core._multiarray_umath import __cpu_features__ as F; print('[npyenv803] enabled:', ' '.join(k for k,v in F.items() if v))"
exec "$@"
