"""Known answers for stats800 (light, local): its per-fold arithmetic must reproduce 522's exact arms and its bootstrap must
land on 522's intervals (same statistic, different draw order, so close, not equal).
    uv run python experiments/800-checker/tests/test_stats800.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
import stats800 as S  # noqa: E402

W522 = ROOT / "experiments/522-october-bulk-key/work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def arm(z, name):
    return dict(x=z[f"{name}__x"], a=z[f"{name}__a"], d=z[f"{name}__d"], pds=z[f"{name}__pds"])


def test_exact_arms():
    worst = 0.0
    n = 0
    for npz in sorted(W522.glob("arms522-*-*.npz")):
        z = np.load(npz)
        J = json.loads(npz.with_suffix(".json").read_text())["arms"]
        base = arm(z, "cap3")
        for nm in ("cap6", "pctPROD", "S1", "cap1e+09"):
            if f"{nm}__x" not in z.files:
                continue
            raw = S.bulk_fold_raw(base, arm(z, nm))
            worst = max(worst, abs(raw["r"] - J[nm]["r"]), abs(raw["pds"] - J[nm]["dpds"]))
            n += 1
    assert n >= 30 and worst < 1e-12, (n, worst)
    return dict(rows=n, max_abs_diff=worst)


def test_bootstrap_vs_522():
    ref = json.loads((W522 / "boot522.json").read_text())
    out = {}
    for panel in ("abc", "efr", "real", "efp"):
        Z = {f: np.load(W522 / f"arms522-{panel}-{f}.npz") for f in CANCER}
        folds = [S.Fold(f"{panel}/{f}", "bulk", arm(Z[f], "cap3"), arm(Z[f], "pctPROD")) for f in CANCER]
        # r only (522's statistic): price with members=("mse",) makes overall proportional to r
        res = S.price(folds, members=("mse",), n_boot=2000, seed=522, top_out=3)
        r_pt = res["members"]["r"]["raw"]
        r_ci = res["members"]["r"]["ci95"]
        want = ref[panel]["pctPROD vs cap3"]
        assert abs(r_pt - want["r"]) < 1e-12, (panel, r_pt, want["r"])
        w = want["ci95"][1] - want["ci95"][0]
        assert abs(r_ci[0] - want["ci95"][0]) < 0.15 * w and abs(r_ci[1] - want["ci95"][1]) < 0.15 * w, (panel, r_ci, want)
        out[panel] = dict(r=r_pt, ci_800=r_ci, ci_522=want["ci95"])
    return out


def test_small():
    assert abs(S.e_max_normal(5) - 1.1629644736) < 1e-9
    assert 1.5 < S.e_max_normal(10) < 1.6
    cal = S.calibrate([dict(**{"class": "k"}, checker=0.002, board=0.001), dict(**{"class": "k"}, checker=0.003, board=0.0015),
                       dict(**{"class": "k"}, checker=-0.001, board=0.001)])
    assert cal["k"]["rho"] == 2.0 and cal["k"]["sign_fail"] == 1
    assert S.decide(0.0001, 0.001, 0.002) == "SHIP" and S.decide(-0.003, -0.002, -0.001) == "HARM"
    # 702 vs ALLSW reproduces the ledger's primary from raw deltas
    assert abs(S.overall_of(dict(pds=-0.0041918, mse=0.0048823)) - (-0.0023631)) < 2e-6
    return "ok"


if __name__ == "__main__":
    print("exact arms:", test_exact_arms())
    print("bootstrap vs 522:", json.dumps(test_bootstrap_vs_522(), indent=1))
    print("small:", test_small())
