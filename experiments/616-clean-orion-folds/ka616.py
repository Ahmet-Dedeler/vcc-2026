#!/usr/bin/env python3
"""616 S0 known answers. (a) gates480 through the sandbox on the THINNED HEK harness == 480's committed gate table
(every array bit-equal). (b) every identity-mode sandbox output == its committed file, byte for byte.
A file that differs is reported with the first differing key or line: drift between a committed output and the current
code/inputs is a finding about that instrument, and its 'thinned' baseline is then the identity run, not the commit.
    uv run python ka616.py <gates-tag> <chain-tag>   -> work/ka616.json; exit 0 iff (a) passes and (b) has no FAIL"""
import json
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
gtag, ctag = sys.argv[1], sys.argv[2]
res = {}
a = np.load(ROOT / "experiments/480-nmae-instrument/work/gate-orion_hek293t.npz", allow_pickle=True)
b = np.load(HERE / f"work/sandbox/{gtag}/experiments/480-nmae-instrument/work/gate-orion_hek293t.npz", allow_pickle=True)
# labels / counts exact; p, q, lfc to 1e-9 with identical NaN patterns (480 §1: the block size alone moves lfc ~7e-12)
bad, dmax = [], {}
for k in a.files:
    if k not in b.files:
        bad.append(k)
        continue
    x, y = np.asarray(a[k]), np.asarray(b[k])
    if x.dtype.kind == "f" and x.shape == y.shape:
        if not np.array_equal(np.isnan(x), np.isnan(y)):
            bad.append(k + ":nan")
            continue
        d = float(np.nanmax(np.abs(x.astype(float) - y.astype(float)))) if x.size else 0.0
        dmax[k] = d
        if d > 1e-9:
            bad.append(f"{k}:{d:.1e}")
    elif not np.array_equal(x, y):
        bad.append(k)
res["gates480_thinned_hek"] = "PASS" if not bad and set(a.files) == set(b.files) else f"FAIL {bad}"
res["gates480_max_abs_diff"] = dmax
print("gates480 (thinned HEK through the sandbox) == committed:", res["gates480_thinned_hek"])
OUTS = ["experiments/639-identity-axis-lolo/work/decomp639.json",
        "experiments/613-block-prior-validation/work/block613.json", "experiments/613-block-prior-validation/work/block613.out",
        "experiments/613-block-prior-validation/work/block613c.out",
        "experiments/613-block-prior-validation/work/block613d.json", "experiments/613-block-prior-validation/work/block613d.out",
        "experiments/613-block-prior-validation/work/block613b.out",
        "experiments/613-block-prior-validation/work/block613d3.out", "experiments/613-block-prior-validation/work/block613d4.out",
        "experiments/613-block-prior-validation/work/block613d5.out", "experiments/613-block-prior-validation/work/block613d6.out",
        "experiments/613-block-prior-validation/work/block613e.out", "experiments/613-block-prior-validation/work/block613e2.out",
        "experiments/612-key-gated-hook/work/ds612b.json", "experiments/612-key-gated-hook/work/ds612c.json",
        "experiments/633-updose-folds/work/dose633.json"]
for rel in OUTS:
    c, s = ROOT / rel, HERE / f"work/sandbox/{ctag}/{rel}"
    if not s.exists():
        res[rel] = "MISSING (instrument not run or failed)"
    elif c.read_bytes() == s.read_bytes():
        res[rel] = "PASS (byte-identical)"
    else:
        cl, sl = c.read_text().splitlines(), s.read_text().splitlines()
        i = next((k for k in range(min(len(cl), len(sl))) if cl[k] != sl[k]), min(len(cl), len(sl)))
        res[rel] = f"DIFFERS at line {i + 1}: committed {cl[i][:90] if i < len(cl) else '<eof>'!r} vs sandbox {sl[i][:90] if i < len(sl) else '<eof>'!r}"
    print(f"{rel}: {res[rel]}")
(HERE / "work/ka616.json").write_text(json.dumps(res, indent=1))
ok = res["gates480_thinned_hek"] == "PASS" and not any(isinstance(v, str) and v.startswith("MISSING") for v in res.values())
sys.exit(0 if ok else 1)
