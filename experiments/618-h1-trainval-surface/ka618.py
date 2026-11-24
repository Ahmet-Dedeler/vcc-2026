#!/usr/bin/env python3
"""618 known answers for the faithful-key pipeline (coordinator condition 1), on the `ka` panel (h1-test's 85 targets):
  (1) work/ka/effects.npz == 323/work/effects.npz: every arm incl. cap15, plus genes / targets, bit for bit;
  (2) work/ka/whitener_h1_pruned.npz == 323/work/whitener_h1_pruned.npz, every array;
  (3) 047's extractor re-run on h1-test's panel (work/ka/cd4_ka_check) == 273's cd4_compact_h1, every .npy / .txt.
Exit 0 iff all three pass. -> work/ka618.json"""
import json
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W323 = ROOT / "experiments/323-current-family-validation/work"
out, ok = {}, True


def cmp_npz(a, b):
    za, zb = np.load(a, allow_pickle=False), np.load(b, allow_pickle=False)
    bad = sorted(set(za.files) ^ set(zb.files))
    for k in sorted(set(za.files) & set(zb.files)):
        if not np.array_equal(za[k], zb[k], equal_nan=za[k].dtype.kind == "f"):
            bad.append(k)
    return bad, sorted(za.files)


for name, a, b in (("effects", W323 / "effects.npz", HERE / "work/ka/effects.npz"),
                   ("operator", W323 / "whitener_h1_pruned.npz", HERE / "work/ka/whitener_h1_pruned.npz")):
    if not b.exists():
        out[name] = "MISSING"
        ok = False
        continue
    bad, keys = cmp_npz(a, b)
    out[name] = "PASS (bit-equal: " + ", ".join(keys) + ")" if not bad else f"FAIL {bad}"
    ok &= not bad
src = ROOT / "experiments/273-flex-reference-instrument/work/cd4/cd4_compact_h1"
chk = HERE / "work/ka/cd4_ka_check"
bad = []
for c in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr"):
    for f in sorted((src / c).iterdir()):
        g = chk / c / f.name
        if not g.exists():
            bad.append(f"{c}/{f.name}:missing")
        elif f.suffix == ".npy":
            if not np.array_equal(np.load(f), np.load(g), equal_nan=True):
                bad.append(f"{c}/{f.name}")
        elif f.suffix == ".txt" and f.read_text() != g.read_text():
            bad.append(f"{c}/{f.name}")
out["cd4_extraction"] = "PASS (every .npy / .txt equal)" if not bad else f"FAIL {bad}"
ok &= not bad
for k, v in out.items():
    print(f"{k}: {v}")
print(f"KA618: {'PASS' if ok else 'FAIL'}")
(HERE / "work/ka618.json").write_text(json.dumps({"pass": ok, **out}, indent=1))
sys.exit(0 if ok else 1)
