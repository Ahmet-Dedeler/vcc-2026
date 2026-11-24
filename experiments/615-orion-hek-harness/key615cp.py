#!/usr/bin/env python3
"""615 G3: the call-propensity block key (528's `lcap` = lc + 22·x) for the Orion raw HEK293T harness (copy of key614cp.py, paths only), from its OWN controls only.
  lc = log10(max(emit_ref_lin, 1e-3)) on the scorer's 5-CPM universe; x = excess dropout (feat528.features, 517's rule);
  beta 22 = 528's value, selected on h1-trainval capture (keys528.json), not refit here.
KA first: the same function on H1 (323's reference + feat528_H1.npz's x) must reproduce keys528.npz's lcap EXACTLY.
Output: work/keys615cp.npz (targets, lcap tiled per target on the gated axis, in the HEK reference's target order)."""
import json
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
sys.dont_write_bytecode = True
import feat528  # noqa: E402  (read-only)
BETA = json.loads((ROOT / "experiments/528-callprop-block/work/keys528.json").read_text())["info"]["lcap"]["beta"]


def lcap(emit_ref_lin_keep, x_keep):
    return np.log10(np.maximum(np.asarray(emit_ref_lin_keep, float), 1e-3)) + BETA * np.asarray(x_keep, float)


z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
keep = np.asarray(z["keep"], bool)
erl = np.asarray(z["emit_ref_lin"], float)
erl_k = erl[keep] if erl.shape[0] == keep.shape[0] else erl
xh = np.asarray(np.load(ROOT / "experiments/528-callprop-block/work/feat528_H1.npz", allow_pickle=True)["x"], float)[keep]
mine = lcap(erl_k, xh)
ref = np.asarray(np.load(ROOT / "experiments/528-callprop-block/work/keys528.npz", allow_pickle=False)["lcap"], float)
ka = all(np.array_equal(mine, ref[i]) for i in range(ref.shape[0]))
print(f"KA lcap on H1 == keys528 (all {ref.shape[0]} rows, bit-equal): {'PASS' if ka else 'FAIL'}; beta {BETA}")
if not ka:
    raise SystemExit("lcap KA failed")
rk = np.load(HERE / "work/ref_hek/reference.npz", allow_pickle=False)
keepk = np.asarray(rk["keep"], bool)
erlk = np.asarray(rk["emit_ref_lin"], float)
erlk = erlk[keepk] if erlk.shape[0] == keepk.shape[0] else erlk
F = feat528.features(HERE / "work/hek-truth.h5")
assert list(F["genes"]) == list(rk["genes"].astype(str)), "feature / reference gene axes differ"
assert np.array_equal(F["keep"], keepk), "feature keep != reference keep"
dm = float(np.max(np.abs(F["mean_cpm"][keepk] - erlk) / np.maximum(erlk, 1e-12)))
print(f"HEK features: n_ctrl {int(F['n_ctrl'])}, lib median {float(F['lib_median']):.1f}; mean_cpm vs emit_ref_lin max rel err {dm:.1e}")
key = lcap(erlk, F["x"][keepk])
T = rk["targets"].astype(str)
np.savez(HERE / "work/keys615cp.npz", targets=T, lcap=np.tile(key, (len(T), 1)))
print(f"wrote work/keys615cp.npz: {len(T)} targets x {key.size} gated genes")
