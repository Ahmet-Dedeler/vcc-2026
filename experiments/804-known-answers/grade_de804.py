"""804 — grade the checker's DE half on the known answers (PROTOCOL.md: DE members fid + jac + nmae at regime f 0.35, reach
reported unpriced; paired target bootstrap 4,000, seed 800, N = 1; both DE-H1 surfaces pooled with equal weight; seeds averaged
per target first). Inputs: 803's per-target members (work/out/<surface>/<arm>/s<k>/members.npz, S3 mirror `outmerged/`).
Also reports f 1 and f 0.61 (619's amendment for call-set size) beside the frozen f 0.35. Light (MB arrays).
    uv run python experiments/804-known-answers/grade_de804.py --root experiments/803-de-engine/outmerged
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

BOARD = json.loads((HERE / "work/board_pairs.json").read_text())["pairs"]
# row -> (cand arm, base arm, board overall, class, surfaces)
ROWS = {
    "446bj2z zs":   ("zs", "g1pre", BOARD["446bj2z-vs-446bj2"]["d_overall"], "de-mag", ("h1test", "h1tv")),
    "496a hook":    ("c4hpre", "g1pre", BOARD["496a-vs-446bj2"]["d_overall"], "de-sign", ("h1tv",)),
    "496a hook (618 form)": ("c4h", "zs", BOARD["496a-vs-446bj2"]["d_overall"], "de-sign", ("h1test", "h1tv")),
    "DN":           ("updn", "upr", BOARD["DN-vs-630arm"]["d_overall"], "de-mag", ("h1test", "h1tv")),
    "630arm UPRN":  ("upr", "zs", BOARD["630arm-vs-446bj2z"]["d_overall"], "de-mag", ("h1test", "h1tv")),
    "idall3 vs hook": ("idall3", "c4h", -0.000501, "de-sign", ("h1test", "h1tv")),
    "K20":          ("k20", "zs", BOARD["446bj2zk20-vs-446bj2z"]["d_overall"], "de-callset", ("h1test", "h1tv")),
    "K16":          ("k16", "zs", BOARD["446bj2zk16-vs-446bj2z"]["d_overall"], "de-callset", ("h1test", "h1tv")),
}


def load(root: Path, surface: str, arm: str):
    d = root / surface / arm
    seeds = sorted(p.name for p in d.glob("s*") if (p / "members.npz").exists())
    return seeds, {s: dict(np.load(d / s / "members.npz", allow_pickle=True)) for s in seeds}


def grade(board, dec):
    if abs(board) < 3e-4:
        return "not graded"
    if board > 0:
        return "HIT" if dec == "SHIP" else "MISS"
    return ("HIT (HARM)" if dec == "HARM" else "HIT") if dec != "SHIP" else "MISS"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(ROOT / "experiments/803-de-engine/outmerged"))
    ap.add_argument("--boot", type=int, default=4000)
    a = ap.parse_args()
    root = Path(a.root)
    out = {}
    for row, (cand, base, board, cls, surfaces) in ROWS.items():
        for regime in ("f0.35", "f1", "f0.61"):
            folds = []
            for sf in surfaces:
                sb, B = load(root, sf, base)
                sc, C = load(root, sf, cand)
                seeds = sorted(set(sb) & set(sc))
                if not seeds:
                    continue
                folds.append(S.Fold(sf, "de", [B[s] for s in seeds], [C[s] for s in seeds], regime=regime,
                                    meta=dict(seeds=seeds)))
            if not folds:
                continue
            res = S.price(folds, members=("fid", "jac", "nmae"), n_boot=a.boot, seed=800, top_out=3)
            dec = S.decide(res["ci95"][0], res["point"], res["ci95"][1])
            m = res["members"]
            rec = dict(row=row, regime=regime, cls=cls, cand=cand, base=base, surfaces=[f.fid for f in folds],
                       point=res["point"], ci95=res["ci95"], p_le0=res["p_le0"], decision=dec, board=board,
                       grade=grade(board, dec) if regime == "f0.35" else None,
                       rho=res["point"] / board if board and np.sign(res["point"]) == np.sign(board) else None,
                       members={k: v["raw"] for k, v in m.items()},
                       per_fold=[dict(fold=p["fold"], overall=p["overall"], raw=p["raw"]) for p in res["per_fold"]])
            out[f"{row}|{regime}"] = rec
            mm = " ".join(f"{k} {v['raw']:+.5f}" for k, v in m.items())
            print(f"{row:22s} {regime:5s} checker {res['point']:+.5f} [{res['ci95'][0]:+.5f}, {res['ci95'][1]:+.5f}] {dec:5s} | "
                  f"board {board:+.5f} -> {rec['grade'] or '':10s} | {mm}", flush=True)
    (HERE / "work/grade_de804.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
