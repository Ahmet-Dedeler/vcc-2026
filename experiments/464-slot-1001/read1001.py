"""Read the 10-01 entries against their preregistered bins (reader written 2026-09-30 11:11Z, before any 10-01 entry exists;
the bins themselves are older and unchanged).

- lodo: 460dr (the October path + the LODO stack key, day-of fit) vs 453dr (q3VuTImh4FRGKXH2waCp).
        Bins: experiments/460-lodo-dayof/PREREGISTRATION.md (2026-09-26 10:22Z): pds and mse both at 0.7-1.3x of
        446bj2 - 446b (YvVru3DAF89Jiyr34olR - UDgjOKL9gI1udAexq6ue), same sign.
- cp  : 446bj2zcp (446bj2z + 528's call-propensity block) vs 446bj2z (be90smZUDcIV7j5MgQ03).
        Bins: experiments/535-callprop-arm/PREREGISTRATION.md (528's, 2026-09-29 16:45Z): primary d(fid + jac) overall.

usage: uv run python experiments/464-slot-1001/read1001.py lodo <460dr entry> [cp <446bj2zcp entry>]
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
BASE = {"lodo": "q3VuTImh4FRGKXH2waCp", "cp": "be90smZUDcIV7j5MgQ03"}
KEY_NEW, KEY_OLD = "YvVru3DAF89Jiyr34olR", "UDgjOKL9gI1udAexq6ue"   # 446bj2, 446b: the A/B/C key's own board delta


def _find(x):
    if isinstance(x, dict):
        if "score_avg" in x:
            return x
        for v in x.values():
            y = _find(v)
            if y:
                return y
    if isinstance(x, list):
        for v in x:
            y = _find(v)
            if y:
                return y


def load(i):
    return _find(json.loads((ROOT / f"data/board-entries/{i}.json").read_text()))


args = sys.argv[1:]
while args:
    kind, e, args = args[0], args[1], args[2:]
    a, b = load(e), load(BASE[kind])
    d = {m: a[RAW[m]] - b[RAW[m]] for m in RAW}
    ds = {m: a[f"score_{m}"] - b[f"score_{m}"] for m in RAW}
    dov = a["score_avg"] - b["score_avg"]
    print(f"\n== {kind} {e} vs {BASE[kind]}: overall {a['score_avg']:.7f}  d {dov:+.7f}")
    for m in RAW:
        print(f"  {m:5s} raw d {d[m]:+.6f}  scaled d {ds[m]:+.6f}  -> overall {ds[m] / 6:+.6f}")
    if kind == "lodo":
        k1, k0 = load(KEY_NEW), load(KEY_OLD)
        kd = {m: k1[RAW[m]] - k0[RAW[m]] for m in ("pds", "mse")}
        r = {m: (d[m] / kd[m] if kd[m] else float("nan")) for m in ("pds", "mse")}
        print(f"  446bj2 - 446b: pds {kd['pds']:+.6f}  mse {kd['mse']:+.6f}  ->  ratio pds {r['pds']:.3f}  mse {r['mse']:.3f}")
        leak = [m for m, lim in (("nmae", 0.0015), ("reach", 0.010)) if abs(d[m]) > lim]
        print("  DE within ballast noise:", "yes" if not leak else f"NO {leak} -- read by hand")
        if any(v < 0 for v in r.values()):
            print("  BIN: opposite sign on pds or mse: a PATH DEFECT. Keep DEF_DAY_LODO off until diagnosed")
        elif all(0.7 <= v <= 1.3 for v in r.values()):
            print("  BIN: both in [0.7, 1.3] x, same sign: the October LODO path REPRODUCES the board key. DEF_DAY_LODO=1 "
                  "is the day's default (DEF-RUNBOOK)")
        else:
            keep = "keep it ON (overall >= +0.001)" if dov >= 0.001 else "overall < +0.001: read by hand before keeping it"
            print(f"  BIN: same sign, outside [0.7, 1.3]: it transfers at a different size; {keep}; record the ratios")
        print("  predicted: pds and mse at 0.7-1.3x of 446bj2 - 446b, same sign; DE |nmae| <= 0.0015, |reach| <= 0.010")
    elif kind == "cp":
        prim = (ds["fid"] + ds["jac"]) / 6
        bad = []
        if abs(d["pds"]) > 1e-6 or abs(d["mse"]) > 1e-6:
            bad.append("pds/mse moved > 1e-6 raw (a DE-only knob leaked into the bulk)")
        if abs(d["nmae"]) > 0.003:
            bad.append("|d nmae| > 0.003 raw")
        if d["jac"] <= 0 and prim >= 0.0004:
            bad.append("jac did not rise while fid+jac did: not the predicted (jac-led) mechanism")
        print(f"  primary d(fid + jac) overall {prim:+.7f}   (fid {ds['fid'] / 6:+.7f}, jac {ds['jac'] / 6:+.7f})")
        print("  FALSIFIERS:", "none" if not bad else "; ".join(bad))
        if prim >= 0.0012:
            print("  BIN: POSITIVE (>= +0.0012): the propensity block transfers; DEF_DAY_CALLSET_BLOCK=prop joins the October "
                  "bundle and the 10-07 rule considers 446bj2zmzcp")
        elif prim >= 0.0004:
            rider = d["jac"] >= 0.002 and d["nmae"] <= 0.001
            print(f"  BIN: PARTIAL (+0.0004 .. +0.0012): October rider {'YES' if rider else 'NO'} "
                  f"(needs d jac raw >= +0.002 and d nmae <= +0.001 raw: jac {d['jac']:+.4f}, nmae {d['nmae']:+.4f})")
        elif prim > -0.0004:
            print("  BIN: NULL (|d| < 0.0004): no October rider")
        else:
            print("  BIN: NEGATIVE (<= -0.0004): close the propensity block on this construction")
        print("  predicted d(fid + jac) +0.0020 .. +0.0025 overall; reach and nmae reported, not binned")
    else:
        raise SystemExit(f"unknown kind {kind}")
