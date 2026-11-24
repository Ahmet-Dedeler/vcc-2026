#!/usr/bin/env python3
"""708 floor sensitivity (science [d80a10], 10-16): score708's arms priced on all fold targets vs on the targets that clear
Arc's basal-expression floor IN THE HELD-OUT LINE (floor708.py's F_abs: held CPM >= 13.4122, unmeasured counts as below).

Inputs: 522's fold dumps (real / efp / efr; the same fields dump708 writes, so score708.Fold reads them) and floor708.json.
Arms and price exactly as score708.py (r, s, rs, s|r; 522 board units), each fold scored with mask = None and with
mask = floor-kept. No leave-top-3-out (a sensitivity read, not a decision). Folds run one at a time (< 1 GB peak each).
usage: uv run python experiments/708-def-panel-repricing/floorprice708.py [panel ...]     # default: real efp efr
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import score708 as S7  # noqa: E402
import fold522 as FD  # noqa: E402  (on sys.path through score708)
import lib438 as L  # noqa: E402

W522 = ROOT / "experiments/522-october-bulk-key/work"


def arms_for(S, panel):
    """score708.main's arm construction for one fold (identical lines)."""
    wp = FD.wclip(S.kv, S.carried, 0.5, FD.pct_hi(S.kv, S.carried, FD.P_PROD))
    z = S7.zsrc_fold(S.targets, panel, S.sources)
    wr = np.where(wp <= 0.5 + 1e-12, wp, np.maximum(0.5, wp * np.exp(0.25 * z)))
    U = S7.unit_mean(S.targets, S.sources, S.nonpanel)
    have = np.linalg.norm(U, axis=1) > 0
    tier = np.full(S.n, -1)
    r_ = np.argsort(np.argsort(z[have], kind="stable"), kind="stable")
    tier[np.flatnonzero(have)] = np.minimum((r_ * 3) // max(int(have.sum()), 1), 2)
    Fs = S7.swap_fin(S.fin, U, tier == 2, S.nonpanel)
    return wp, {"r": (wr, None), "s": (wp, Fs), "rs": (wr, Fs)}


def main(panels):
    fl = json.loads((HERE / "work" / "floor708.json").read_text())
    out = {}
    for p in panels:
        bundle = Path(L.PANELS[p]["bundle"])
        panel = [r["target_gene"] for r in csv.DictReader(open(bundle / "pert_counts.csv")) if r["target_gene"] != "non-targeting"]
        res = {"all": {}, "floor": {}}
        for f in S7.FOLDS:
            t0 = time.time()
            S = S7.Fold(W522 / f"dump522-{p}-{f}.npz")
            below = set(fl["folds"][f"{p}/{f}"]["below_abs_targets"])
            keep = np.array([t not in below for t in S.targets])
            wp, arms = arms_for(S, panel)
            for ver, mask in (("all", None), ("floor", keep)):
                si = S.score(wp, mask=mask)[0]
                sc = {name: S.score(w, fin=fin, mask=mask)[0] for name, (w, fin) in arms.items()}
                pr = {name: S7.price(s, si) for name, s in sc.items()}
                pr["s|r"] = S7.price(sc["rs"], sc["r"])
                res[ver][f] = pr
            print(f"{p}/{f}: n {S.n}, floor-kept {int(keep.sum())} | " + " | ".join(
                f"{a} {res['all'][f][a]:+.5f} -> {res['floor'][f][a]:+.5f}" for a in ("r", "s", "rs", "s|r"))
                + f" ({time.time() - t0:.0f}s)", flush=True)
            del S
        summ = {}
        for a in ("r", "s", "rs", "s|r"):
            row = {}
            for ver in ("all", "floor"):
                v = np.array([res[ver][f][a] for f in S7.FOLDS])
                row[ver] = dict(mean=float(v.mean()), n_pos=int((v > 0).sum()), n_neg=int((v < 0).sum()),
                                move=bool(v.mean() >= 0.0010 and (v > 0).all()),
                                veto=bool(v.mean() <= -0.0005 and (v < 0).sum() >= 2))
            summ[a] = row
            print(f"{p} {a:4s}: all {row['all']['mean']:+.5f} ({row['all']['n_pos']}/3 pos) move {row['all']['move']} "
                  f"veto {row['all']['veto']} | floor {row['floor']['mean']:+.5f} ({row['floor']['n_pos']}/3 pos) "
                  f"move {row['floor']['move']} veto {row['floor']['veto']}", flush=True)
        out[p] = dict(per_fold=res, summary=summ)
    (HERE / "work" / "floorprice708.json").write_text(json.dumps(out, indent=1))
    print("wrote work/floorprice708.json")


if __name__ == "__main__":
    main(sys.argv[1:] or ["real", "efp", "efr"])
