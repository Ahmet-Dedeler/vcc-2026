"""318 — read the six raw members off every arm's `agg_results.csv`, price them in board
overall points through 317 §1.4's fitted map, and report the paired per-target intervals.

`base` is 221's `rc_callset`, gated to sha256 `048b9286…` in every arm's build (gate 5). It is
NOT re-scored here: `cell-eval2 run` reads only (prediction, reference, version, preset), and
317 §1.2 proved the recorded row reproduces bit for bit from the same object on a fresh worker.
317's `ka1-317/base/run/` is therefore this experiment's baseline AND its known-answer row.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = HERE / "cloud-output"
BASE_AGG = REPO / "experiments/317-operating-point/cloud-output/ka1-317/base/run/agg_results.csv"
BASE_RES = REPO / "experiments/317-operating-point/cloud-output/ka1-317/base/run/results.csv"

ARMS = ["ptrule", "ptshuf1", "ptshuf2", "ptshuf3", "ptamp6", "ptunif"]
MEMBERS = ["pds", "mse", "nmae", "fid", "reach", "jac"]
RAWKEY = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm",
          "nmae": "de_wilcoxon_lfc_nmae",
          "fid": "de_wilcoxon_direction_fidelity_yield_raw",
          "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}

# 317 §1.2's recorded known-answer row for `rc_callset` (221's, reproduced on a fresh worker).
KA_221 = {"pds": 0.873809523810, "mse": 1.015977751768, "nmae": 0.960046422899,
          "fid": 0.495058277089, "reach": 0.226704287761, "jac": 0.074203737608}


def load(mod_path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, mod_path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


P317 = load(REPO / "experiments/317-operating-point/price317.py", "price317_318")
PAIR = load(REPO / "experiments/317-operating-point/paired317.py", "paired317_318")


def agg_row(path: Path) -> dict[str, float]:
    rows = list(csv.DictReader(path.open()))
    mean = next(r for r in rows if r["statistic"] == "mean")
    return {m: float(mean[RAWKEY[m]]) for m in MEMBERS}


def main() -> None:
    slopes = json.loads((REPO / "experiments/317-operating-point/work/board_slopes.json"
                         ).read_text())["slopes"]

    base = agg_row(BASE_AGG)
    ka = {m: abs(base[m] - KA_221[m]) for m in MEMBERS}
    print("KA-1 (reused from 317 §1.2) — `base` is 221's recorded rc_callset row")
    for m in MEMBERS:
        print(f"   {m:6s} {base[m]:.12f}  vs recorded {KA_221[m]:.12f}   |diff| {ka[m]:.2e}")
    ka_ok = max(ka.values()) <= 1e-9
    print(f"   KA-1 {'PASS' if ka_ok else 'FAIL'}, max |diff| {max(ka.values()):.2e}\n")
    if not ka_ok:
        raise SystemExit("KA-1 FAILED — the baseline row is not 221's. STOP.")

    res = {"ka1": {"base": base, "recorded": KA_221, "max_abs_diff": max(ka.values()),
                   "pass": ka_ok},
           "arms": {}}
    print(f"{'arm':9s} " + " ".join(f"{m:>12s}" for m in MEMBERS) + "   d overall")
    print(f"{'base':9s} " + " ".join(f"{base[m]:12.6f}" for m in MEMBERS))
    for arm in ARMS:
        p = OUT / f"s318-{arm}" / arm / "run" / "agg_results.csv"
        if not p.exists():
            print(f"{arm:9s} (not landed yet)")
            continue
        row = agg_row(p)
        d = {m: row[m] - base[m] for m in MEMBERS}
        _per, _ov = P317.price(d, slopes, factor=1.0)
        pr = {"per_member": _per, "overall": _ov}
        rp = (OUT / f"s318-{arm}" / arm / "run" / "results.csv")
        pair = PAIR.paired(BASE_RES, rp) if rp.exists() else {}
        bj = list((OUT / f"s318-{arm}" / "prologue-out").glob("build318-*.json"))
        build = json.loads(bj[0].read_text()) if bj else {}
        res["arms"][arm] = {"raw": row, "delta_raw": d, "priced": pr, "paired": pair,
                            "multiplier": (build.get("multipliers", {}) or {}).get(arm, {}),
                            "gates": {k: v for k, v in (build.get("gates", {}) or {}).items()
                                      if k in ("gate4_panel_total_digest",
                                               "gate5_base_is_221_rc_callset",
                                               "gate6_shuffle_multiset",
                                               "gate7_delivered_magnitude")}}
        print(f"{arm:9s} " + " ".join(f"{row[m]:12.6f}" for m in MEMBERS)
              + f"   {pr['overall']:+.6f}")

    print(f"\n{'arm':9s} " + " ".join(f"{'d ' + m:>12s}" for m in MEMBERS) + "   d overall")
    for arm, v in res["arms"].items():
        print(f"{arm:9s} " + " ".join(f"{v['delta_raw'][m]:+12.6f}" for m in MEMBERS)
              + f"   {v['priced']['overall']:+.6f}")

    print(f"\n{'arm':9s} " + " ".join(f"{'$ ' + m:>12s}" for m in MEMBERS) + "   d overall")
    for arm, v in res["arms"].items():
        print(f"{arm:9s} " + " ".join(f"{v['priced']['per_member'][m]:+12.6f}"
                                      for m in MEMBERS)
              + f"   {v['priced']['overall']:+.6f}")

    # ------------------------------------------- the uniform-corrected (pure-redistribution)
    # `effect_scale` acts before `profile_from_effect` renormalises, so the DELIVERED
    # whole-matrix norm of the arms is not perfectly equal even though the effect-level norm
    # is (gate 7). `ptunif` measures the price of a pure uniform excess, so the residual is
    # subtracted rather than assumed away: for arm a,
    #     corrected(a) = delta(a) - (excess(a)/excess(ptunif)) * delta(ptunif)
    # with excess = delivered whole-matrix norm - 1. corrected(ptunif) is exactly 0 and
    # corrected(a) is what redistribution alone bought.
    if "ptunif" in res["arms"]:
        gw = {a: res["arms"][a]["gates"]["gate7_delivered_magnitude"][a]["delivered_whole"]
              for a in res["arms"]}
        eu = gw["ptunif"] - 1.0
        print(f"\nUNIFORM-CORRECTED (delivered whole-matrix excess removed; ptunif "
              f"excess {eu:+.6f} = the unit)")
        print(f"{'arm':9s} {'excess':>10s} {'k':>7s} " +
              " ".join(f"{'d ' + m:>12s}" for m in MEMBERS) + "   d overall")
        slopes_l = slopes
        for a, v in res["arms"].items():
            k = (gw[a] - 1.0) / eu if eu else 0.0
            d = {m: v["delta_raw"][m] - k * res["arms"]["ptunif"]["delta_raw"][m]
                 for m in MEMBERS}
            _per, _ov = P317.price(d, slopes_l, factor=1.0)
            pr = {"per_member": _per, "overall": _ov}
            v["delivered_whole"] = gw[a]
            v["uniform_corrected"] = {"k": k, "delta_raw": d, "priced": pr}
            print(f"{a:9s} {gw[a] - 1.0:+10.6f} {k:7.3f} "
                  + " ".join(f"{d[m]:+12.6f}" for m in MEMBERS)
                  + f"   {pr['overall']:+.6f}")

    # ---------------------------------------------------------------- the separability test
    shufs = [a for a in ("ptshuf1", "ptshuf2", "ptshuf3") if a in res["arms"]]
    if "ptrule" in res["arms"] and len(shufs) >= 2:
        print("\nSEPARABILITY — is the replayed profile distinguishable from its own shuffles?")
        print(f"{'member':8s} {'ptrule':>12s} {'shuf mean':>12s} {'shuf sd':>12s} "
              f"{'z':>8s}")
        sep = {}
        for m in MEMBERS + ["overall"]:
            def val(a):
                return (res["arms"][a]["priced"]["overall"] if m == "overall"
                        else res["arms"][a]["delta_raw"][m])
            s = np.array([val(a) for a in shufs])
            r = val("ptrule")
            z = (r - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else float("nan")
            sep[m] = {"rule": r, "shuf_mean": float(s.mean()),
                      "shuf_sd": float(s.std(ddof=1)), "z": float(z)}
            print(f"{m:8s} {r:+12.6f} {s.mean():+12.6f} {s.std(ddof=1):12.6f} {z:+8.2f}")
        res["separability"] = sep

    (HERE / "work/read318.json").write_text(json.dumps(res, indent=1))
    print("\n-> work/read318.json")


if __name__ == "__main__":
    main()
