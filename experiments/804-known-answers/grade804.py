"""804 — grade the checker on the known answers, exactly as PROTOCOL.md (incl. amendments 1-2) fixes it. Light (MB arrays).

Reads work/run/bulk/<fold>/<arm>.npz (ka804.py on EC2) and, when present, work/run/de/<surface>/<arm>/s<k>.npz (803's
engine), prices every row with 800's stats800 at the frozen settings (paired target bootstrap 4,000 draws, seed 800, N = 1),
and grades it against the board (work/board_pairs.json). Writes work/grade804.json and prints the table.
    uv run python experiments/804-known-answers/grade804.py [--surfaces W,ABC]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/800-checker"))
import stats800 as S  # noqa: E402

RUN = HERE / "work" / "run"
BOARD = json.loads((HERE / "work" / "board_pairs.json").read_text())["pairs"]
# row -> (board pair key, class, {surface: (cand, base)}, folds filter)
ROWS = {
    "495":    ("495-vs-446bj2", "bulk-key", {"ABC": ("cap6", "cap3"), "W": ("cap6", "cap3")}, None),
    "495pct": ("495-vs-446bj2", "bulk-key", {"W": ("pct", "cap3")}, None),
    "479p2g": ("479p2g-vs-446bj2", "bulk-shape", {"ABC": ("cap3+P2g", "cap3"), "W": ("cap3+P2g", "cap3")}, None),
    "479sp2": ("479sp2-vs-446bj2", "bulk-shape", {"ABC": ("cap3+S_P2", "cap3"), "W": ("cap3+S_P2", "cap3")}, None),
    "498p2g": ("498p2g-vs-495", "bulk-shape", {"ABC": ("cap6+P2g", "cap6"), "W": ("pct+P2g", "pct")}, None),
    "702":    ("702-vs-ALLSW", "bulk-external", {"ABC": ("702", "inc"), "W": ("702", "inc")}, "k562_gwps"),
    "705":    (None, "bulk-key", {"ABC": ("705", "inc"), "W": ("705", "inc")}, None),
    "707":    ("707-vs-ALLSW", "bulk-direction", {"ABC": ("707", "inc"), "W": ("707", "inc")}, None),   # board read 10-17
    "710":    (None, "bulk-composition", {"ABC": ("710", "inc710"), "W": ("710", "inc710")}, None),
}
PROSPECTIVE = {"705", "710"}   # 707 graded 10-17 (board pair 707-vs-ALLSW)


def surface_of(fold_dir: str) -> str:
    return "ABC" if fold_dir.startswith("ABC__") else "W"


def load_bulk(surface: str, cand: str, base: str, held_filter=None, skip_kolf=True):
    folds = []
    for d in sorted((RUN / "bulk").glob("*")):
        if surface_of(d.name) != surface:
            continue
        held = d.name.split("__")[1]
        if held_filter and held != held_filter:
            continue
        if skip_kolf and held == "kolf2":           # kolf2 is reported separately (stem-cell line, not a cancer fold)
            continue
        b, c = d / f"{base}.npz", d / f"{cand}.npz"
        if b.exists() and c.exists():
            zb, zc = dict(np.load(b, allow_pickle=True)), dict(np.load(c, allow_pickle=True))
            folds.append(S.Fold(d.name.replace("__", "/"), "bulk", zb, zc, meta=dict(
                share_unscored=float((zc["dE_unscored"].sum() + zc["dE_panel"].sum()) / max(zc["dE_total"].sum(), 1e-300))
                if "dE_total" in zc else None)))
    return folds


def board_raw(pair_key):
    """board raw member deltas from the scaled ones (scaled = raw x INV, sign-oriented)."""
    if not pair_key or pair_key not in BOARD:
        return None
    p = BOARD[pair_key]
    d = p["d_scaled"]
    return dict(overall=p["d_overall"], pds=d["pds"] / S.INV["pds"], mse=-d["mse"] / S.INV["mse"],
                nmae=-d["nmae"] / S.INV["nmae"], fid=d["fid"] / S.INV["fid"], reach=d["reach"] / S.INV["reach"],
                jac=d["jac"] / S.INV["jac"])


def grade(board_overall, decision):
    if board_overall is None:
        return "prospective"
    if abs(board_overall) < 3e-4:
        return "not graded"
    if board_overall > 0:
        return "HIT" if decision == "SHIP" else "MISS"
    return ("HIT (HARM)" if decision == "HARM" else "HIT") if decision != "SHIP" else "MISS"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--surfaces", default="W,ABC")
    ap.add_argument("--boot", type=int, default=4000)
    a = ap.parse_args()
    out = {}
    for row, (pk, cls, sides, hf) in ROWS.items():
        brd = board_raw(pk)
        for surf in a.surfaces.split(","):
            if surf not in sides:
                continue
            cand, base = sides[surf]
            F = load_bulk(surf, cand, base, hf)
            if not F:
                continue
            res = S.price(F, members=("pds", "mse"), n_boot=a.boot, seed=800)
            dec = S.decide(res["ci95"][0], res["point"], res["ci95"][1])
            m = res["members"]
            rec = dict(row=row, surface=surf, cls=cls, cand=cand, base=base, n_folds=len(F), point=res["point"],
                       ci95=res["ci95"], p_le0=res["p_le0"], folds_pos=res["folds_pos"], decision=dec,
                       pds_raw=m["pds"]["raw"], mse_raw=m["mse"]["raw"], r=m["r"]["raw"],
                       pds_ci=m["pds"]["ci95"], mse_ci=m["mse"]["ci95"],
                       share_unscored=float(np.nanmean([f.meta.get("share_unscored") or np.nan for f in F])),
                       board=brd, grade=grade(brd["overall"] if brd else None, dec),
                       per_fold=[dict(fold=p["fold"], overall=p["overall"], r=p["raw"]["r"], pds=p["raw"]["pds"]) for p in res["per_fold"]])
            if brd:
                rec["rho_overall"] = res["point"] / brd["overall"] if brd["overall"] else None
                rec["rho_pds"] = m["pds"]["raw"] / brd["pds"] if brd["pds"] else None
                rec["rho_mse"] = m["mse"]["raw"] / brd["mse"] if brd["mse"] else None
            out[f"{row}|{surf}"] = rec
            b = f"{brd['overall']:+.5f}" if brd else "   (prosp.)"
            print(f"{row:7s} {surf:3s} n={len(F):2d} checker {res['point']:+.5f} [{res['ci95'][0]:+.5f}, {res['ci95'][1]:+.5f}] "
                  f"{dec:5s} | board {b} -> {rec['grade']:12s} | pds {m['pds']['raw']:+.5f} mse {m['mse']['raw']:+.5f} "
                  f"r {100 * m['r']['raw']:+.2f}% folds+ {res['folds_pos']}/{len(F)}", flush=True)
    (HERE / "work" / "grade804.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
