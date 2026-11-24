"""172 — assemble the key file experiment 171 (or a successor) drops in as `rank_score`.

Everything here is already computed; this only packages it with its axes and provenance so
that nobody has to re-derive which array is which.  `reach` values are recomputed from the
saved arrays as a final round-trip check, so the file cannot ship a key whose reported number
does not belong to it.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from harness172 import H1Test, check  # noqa: E402

OUT = os.path.join(HERE, "key172.npz")
WANT = {
    "key_M1_lodo_selected": "M1 logistic @ C=0.1 w=top25",
    "key_M4_gbt": "M4 gbt @ it=100 w=invrank",
    "key_M2_step": "M2 step-additive @ C=0.01 w=invrank",
    "key_blend_M1_absH": "50/50 blend with |H|",
    "key_absH": "|H| (whitened)",
    "key_cons": "|cons| (un-whitened)",
}


def main() -> None:
    T = H1Test()
    check(T, verbose=False)
    print("known-answer rows OK (harness172.check)")
    Z = np.load(os.path.join(HERE, "work", "allkeys.npz"), allow_pickle=True)
    arrays = {"genes": T.genes, "targets": T.perts}
    reach = {}
    for name, src in WANT.items():
        a = np.asarray(Z[src], dtype=np.float32)
        assert a.shape == T.H.shape, (name, a.shape)
        arrays[name] = a
        reach[name] = T.reach(a.astype(np.float64))
    meta = dict(
        surface="h1-test, a15 emission",
        de_table="experiments/164-callset-artefact/de164.a15.npz",
        effect="experiments/169-sharpening-six-member/work/h1test-effect.npz",
        donors_used_to_build_the_key=["k562_gwps", "orion_hek293t"],
        donors_used_to_train=["k562_gwps", "orion_hct116", "orion_hek293t"],
        donors_banned=["h1_2025 (covers 85/85 h1-test targets: oracle leakage)"],
        higher_is_more_confident=True,
        reach=reach,
        note=("Counterfactual ordering keys only.  The scorer reads our Wilcoxon p-values, "
              "not this key; 171's expression-compensating spread schedule is the delivery "
              "vehicle.  RECOMMENDED: key_absH -- no learned key beat it (172 RESULTS.md)."))
    np.savez_compressed(OUT, meta=np.asarray(json.dumps(meta)), **arrays)
    sha = hashlib.sha256(open(OUT, "rb").read()).hexdigest()
    print(f"\nwrote {OUT}")
    print(f"sha256 {sha}\n")
    for k, v in reach.items():
        print(f"  {k:<26} reach {v:.5f}")
    print(f"\narrays: {sorted(arrays)}")
    print("meta:", json.dumps(meta["note"]))
    with open(os.path.join(HERE, "work", "key172.sha256"), "w") as f:
        f.write(f"{sha}  key172.npz\n")


if __name__ == "__main__":
    main()
