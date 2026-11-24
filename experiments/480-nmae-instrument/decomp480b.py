"""480 step 2b: the same oracle decomposition on every surface, for the one construction that exists everywhere.

Construction: a_analog = sign(cons) x shipped MDE law + 0.575 * cons (fold480's H1-calibrated analog; cons = the
available panel donors' consensus, never the held-out line's own table). Its "base" for the class-constant add scan is
0.575 * cons. On h1 the real (a) decomposition is decomp480.py; this script shows whether the class structure of the
gap (where the magnitude and sign pots sit) is the same on the folds and on H1 at the folds' power.

    uv run python decomp480b.py     # -> work/decomp480b.json
"""
from __future__ import annotations

import json

import numpy as np

import fold480 as F
from common480 import HERE, load_donor_rows
from decomp480 import CIS_NAMES, decompose, tss_classes

SURF = ["h1", "h1_thin", "orion_hct116", "orion_hek293t", "k562_gwps"]


def main():
    out = {}
    for name in SURF:
        S = F.load_surface(name)
        T, cnt, genes, targets = S["T"], S["cnt"], S["genes"], S["targets"]
        G = np.isfinite(T)
        cons = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, genes, pc_cpm=1.0)[0]
                                                  for d in S["donors"]]), 0))
        s = np.sign(cons)
        s = np.where(s == 0, 1.0, s)
        add = F.mde_add(cnt)[None, :]
        base = F.K_H1_MDECONS * cons
        pred = np.sign(cons) * add + base
        updown = np.where(np.nan_to_num(T) > 0, 0, 1).astype(np.int8)
        q1, q2 = np.quantile(np.broadcast_to(cnt, T.shape)[G], [1 / 3, 2 / 3])
        expr = np.broadcast_to(np.where(cnt < q1, 0, np.where(cnt < q2, 1, 2)).astype(np.int8), T.shape)
        nreal = G.sum(1)
        n_cls = np.broadcast_to(np.digitize(nreal, [30, 100, 500])[:, None].astype(np.int8), T.shape)
        mag = np.digitize(np.abs(cons), np.quantile(np.abs(cons)[G], [1 / 3, 2 / 3])).astype(np.int8)
        classes = {
            "truth direction": (updown, ["up", "down"]),
            "predicted direction (cons)": (np.where(s > 0, 0, 1).astype(np.int8), ["say up", "say down"]),
            "control expression tertile": (expr, ["low", "mid", "high"]),
            "cis position (TSS distance)": (tss_classes(targets, genes), CIS_NAMES),
            "gate size n_real (oracle)": (n_cls, ["10-29", "30-99", "100-499", ">=500"]),
            "donor |cons| tertile": (mag, ["small", "mid", "large"]),
        }
        out[name] = decompose(T, pred, s, base, cnt, classes, f"{name}: a_analog")
    json.dump(out, open(HERE / "work/decomp480b.json", "w"), indent=1)


if __name__ == "__main__":
    main()
