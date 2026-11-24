"""Read the Arc-programme family against its preregistered bins (written 2026-09-26 10:50Z, before any 09-27…09-30 entry
exists; bins in experiments/458-arc-lab-programme/PREREGISTRATION.md and 461-arcprog-october-path/PREREGISTRATION.md).

  prog  <458d or 458 entry> <446b entry>      458's bins: falsifiers pds/mse <= 1e-6 raw, jac within +-0.0010; fid primary
  dose  <458e entry> <446b entry> <458d entry>   458e's bins: fid ratio to 458d's
  path  <461dr entry> <453dr entry> <458d entry> <446b entry>   461's bins: fid ratio to 458d's, same sign; pds/mse pinned

usage: uv run python experiments/458-arc-lab-programme/read458.py <kind> <entries...> [<kind> <entries...>]
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
NARGS = {"prog": 2, "dose": 3, "path": 4}


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


def delta(a, b, tag):
    d = {m: a[RAW[m]] - b[RAW[m]] for m in RAW}
    ds = {m: a[f"score_{m}"] - b[f"score_{m}"] for m in RAW}
    dov = a["score_avg"] - b["score_avg"]
    print(f"\n== {tag}: overall {a['score_avg']:.7f}  d {dov:+.7f}")
    for m in RAW:
        print(f"  {m:5s} raw d {d[m]:+.6f}  scaled d {ds[m]:+.6f}  -> overall {ds[m] / 6:+.6f}")
    return d, dov


def falsifiers(d):
    bad = [m for m in ("pds", "mse") if abs(d[m]) > 1e-6]
    j = abs(d["jac"]) > 0.0010
    print("  FALSIFIERS:", "OK (bulk pinned, jac inside +-0.0010)" if not bad and not j
          else f"FIRED: {bad or ''} {'jac outside +-0.0010' if j else ''} -- a bug or a leak; stop and diagnose")


def fid_bin(f):
    return ("POSITIVE: the same-lab programme carries across lines (October: DEF_DAY_PROGRAMME_SOURCE=h1_2025_full, flags pt)"
            if f >= 0.0015 else "PARTIAL: build/test the CLEAN4 + Arc blend (458e) before October" if f >= 0.0005
            else "NULL: no (the lab component is too small for this construction)" if f > -0.0010
            else "HARM: close the Arc programme")


args = sys.argv[1:]
while args:
    kind = args[0]
    n = NARGS[kind]
    ids, args = args[1:1 + n], args[1 + n:]
    if kind == "prog":
        d, dov = delta(load(ids[0]), load(ids[1]), f"programme arm {ids[0]} vs 446b {ids[1]}")
        falsifiers(d)
        print("  BIN fid (primary):", fid_bin(d["fid"]))
        print("  BIN nmae (secondary):", "confirms" if d["nmae"] <= -0.0015 else "contradicts fid" if d["nmae"] >= 0.0015
              else "inside +-0.0015")
        print("  BIN overall:", "adopt for A/B/C riders" if dov >= 0.0010 else "close" if dov <= -0.0010 else "flat")
        if d["fid"] >= 0.0015 and d["reach"] <= -0.005 and dov <= 0:
            print("  CONTINGENCY (addendum 11:31Z): fid carries but reach pays -> do NOT set the October switch as is;"
                  " board-test a reach-protecting variant (larger q protection) first")
    elif kind == "dose":
        d, dov = delta(load(ids[0]), load(ids[1]), f"blend 458e {ids[0]} vs 446b {ids[1]}")
        d8, _ = delta(load(ids[2]), load(ids[1]), f"(reference) 458d {ids[2]} vs 446b {ids[1]}")
        falsifiers(d)
        r = d["fid"] / d8["fid"] if d8["fid"] else float("nan")
        print(f"  fid ratio 458e/458d = {r:.2f}")
        print("  BIN:", "the consensus ADDS: October uses the 0.5 blend" if d["fid"] > d8["fid"]
              else "monotone: October uses w 1 (Arc alone)" if r >= 0.4 else "non-linear: keep w 1 and record it")
    elif kind == "path":
        d, dov = delta(load(ids[0]), load(ids[1]), f"October path 461dr {ids[0]} vs 453dr {ids[1]}")
        d8, _ = delta(load(ids[2]), load(ids[3]), f"(reference) 458d {ids[2]} vs 446b {ids[3]}")
        bad = [m for m in ("pds", "mse") if abs(d[m]) > 1e-6]
        print("  FALSIFIER (bulk):", "OK" if not bad else f"FIRED {bad}")
        r = d["fid"] / d8["fid"] if d8["fid"] else float("nan")
        print(f"  fid ratio 461dr/458d = {r:.2f}")
        print("  BIN:", "the October programme path REPRODUCES the board: set the switch for Oct 22" if 0.6 <= r <= 1.4
              else "same sign, different size: adopt if d fid >= +0.0010" if r > 0
              else "opposite sign: a path defect; keep the switch off")
