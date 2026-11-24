#!/usr/bin/env python3
"""599 extension prep: the exact ballast K of seeds 3 and 4 for t and tc (the preregistered extension: "if the f 0.35 mean
is within 2 standard errors of 0 ... seeds 3 and 4 are added for t and tc only"). Same arithmetic as prep599.py (its
functions, imported); writes work/prepx599.npz (K_<mode>_s<seed>, and q/carried copies checked equal to prep599's).
No truth is read. ~6 GB (the reference): run through gate602.
    uv run python experiments/599-tercilec-h1-de/prepx599.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prep599 as P  # noqa: E402

C, H, R, E = P.C, P.H, P.R, P.E
SEEDS_X = (3, 4)


def main() -> int:
    b = H.builder389()
    eff, genes, targets = H.native_cap15(b)
    kz = C.fit_key(targets)
    r, cc, emit_c = E.load_reference()
    eff_de = R.map_axis(eff, genes, r["genes"])
    panel = set(map(str, targets))
    panel_mask = np.array([str(g) in panel for g in r["genes"]])
    prof_de = [b.profile_from_effect(eff_de[i], r["b_ctrl"], min_cpm=10.) for i in range(len(targets))]
    reg = np.load(C.WORK / "prep599.npz", allow_pickle=False)
    arrays = {}
    for mode in ("t", "tc"):
        eb, w, carried, c = C.mode_bulk(mode, eff, genes, targets, kz)
        qt, _ = C.mode_q(mode, eb, eff, genes, targets, carried)
        assert np.array_equal(qt, reg[f"q_{mode}"]) and np.array_equal(carried, reg[f"carried_{mode}"])
        eff_bk = R.map_axis(eb, genes, r["genes"])
        pb = [b.profile_from_effect(eff_bk[i], r["b_ctrl"], min_cpm=10., sharpen_q=float(qt[i]), protect=panel_mask)
              for i in range(len(targets))]
        for s in SEEDS_X:
            pr = b.make_rngs(s)[0]
            picks = [pr.choice(len(r["libraries"]), P.NCELL, replace=False) for _ in range(len(targets))]
            arrays[f"K_{mode}_s{s}"] = np.asarray([P.k_of(prof_de[i], pb[i], r["libraries"][picks[i]])[0] for i in range(len(targets))])
            print(f"{mode} s{s}: K mean {arrays[f'K_{mode}_s{s}'].mean():.3f}", flush=True)
    for s in SEEDS_X:
        d = arrays[f"K_tc_s{s}"] - arrays[f"K_t_s{s}"]
        print(f"s{s}: dK (tc - t) mean {d.mean():+.3f}, changed {int((d != 0).sum())}", flush=True)
    np.savez(C.WORK / "prepx599.npz", **arrays)
    return 0


if __name__ == "__main__":
    sys.exit(main())
