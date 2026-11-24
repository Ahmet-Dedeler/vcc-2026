#!/usr/bin/env python3
"""633j — the UP-dose reading on a second clean cancer surface: Nadig 2025 Jurkat (gates633j.py's gate file).

Same model, doses and regimes as dose633.py (its h1 retrodiction PASSED, so the model's SHAPE is readable). The surface
is plugged into 612's fold machinery: the gate file is read from this experiment's work dir, CLEAN4 drops the Jurkat bank
(LOLO: the surface's own data never signs it), and the donor consensus is K562 GWPS + the two Orion lines.

ARTIFACT CHECKS FIRST (609 / 609b; fixed before the run). The surface counts as CLEAN only if all three hold:
  (1) depth: median perturbed library / median control library in [0.85, 1.15];
  (2) CPM-block gate genes' DOWN share <= 0.62 (H1 0.56, K562 0.51; Orion 0.72-0.74);
  (3) uncalled gate genes in the bottom control-expression decile: UP share <= 0.65 (H1 ~0.58; Orion 0.68-0.83).
If any fails, the Jurkat rows are descriptive only.
DECISION (fixed before the run): if CLEAN and the d20 argmin K at both f 0.5 and f 0.35 is > 1.0 AND beats K 1.0 by
> 0.002 raw, Jurkat confirms K562: the October dose question goes to a preregistered K 1.0 / K 2.0 board pair once
capacity can deliver K 2. Otherwise K 1.0 stays the October dose.
  uv run python experiments/633-updose-folds/dose633j.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/612-key-gated-hook"))
import ds612 as D  # noqa: E402
import ds612b as B  # noqa: E402
import ds612c as C  # noqa: E402
from common480 import boot_delta  # noqa: E402

NAME = "nadig_jurkat"
KS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
REG = (1.0, 0.5, 0.35)


def main() -> int:
    gate = np.load(HERE / "work/gate-nadig_jurkat.npz", allow_pickle=True)
    D.F.HERE = HERE                                  # fold480.load_surface reads HERE / work/gate-<name>.npz
    _pv = B.pvals
    B.pvals = lambda name, T: (np.asarray(gate["p"], float) if name == NAME else _pv(name, T))
    D.BANKS[NAME] = ["k562_gwps", "k562_essential", "nadig_hepg2"]
    D.PROFKEY[NAME] = NAME
    BETA, prof = D.setup()
    prof[NAME] = dict(zip(gate["genes"].astype(str), gate["ctrl_cpm"].astype(float)))
    X = D.surface_forms(NAME, BETA, prof)
    T, S = X["T"], X["S"]
    base = X["P"]["plain"]
    called = C.callset(X["cons0"], S["cnt"], 400, 2200)
    head = C.callset(X["cons0"], S["cnt"], 400, 0)
    block = called & ~head
    unc = ~called
    p = B.pvals(NAME, T)
    lines = []

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    # artifact checks
    depth = float(gate["lib_pert_med"]) / float(gate["lib_med"])
    G = np.isfinite(T)
    dn_block = float((T[G & block] < 0).mean())
    cnt = np.asarray(S["cnt"], float)
    low = cnt <= np.quantile(cnt, 0.1)
    up_low = float((T[G & unc & low[None, :]] > 0).mean())
    clean = (0.85 <= depth <= 1.15) and (dn_block <= 0.62) and (up_low <= 0.65)
    say(f"ARTIFACT CHECKS: depth pert/ctrl {depth:.3f} (0.85-1.15), block gate DOWN share {dn_block:.3f} (<= 0.62), "
        f"low-expression uncalled UP share {up_low:.3f} (<= 0.65) -> {'CLEAN' if clean else 'NOT CLEAN (descriptive only)'}")
    say(f"targets {T.shape[0]}, gate genes/target median {np.median(G.sum(1)):.0f}, ctrl lib {float(gate['lib_med']):.0f}")
    sel = unc & (base > 0)
    res = {}
    for dname, cK in (("own", cnt), ("d20", cnt * 20000.0 / float(gate["lib_med"]))):
        mde = 0.24 / np.sqrt(np.maximum(cK, 1e-12))
        rows0 = {f: D.nmae_cells(base, B.pool_T(T, p, f))[0] for f in REG}
        rec = {}
        for K in KS:
            P = np.where(sel, np.maximum(base, K * mde[None, :]), base)
            rec[K] = {f: boot_delta(rows0[f], D.nmae_cells(P, B.pool_T(T, p, f))[0])[:3] for f in REG}
        res[dname] = rec
        say(f"\n== {NAME} [{dname}] uncalled UP-base cells/target {sel.sum(1).mean():.0f}; Δnmae raw vs base (− = better)")
        for f in REG:
            say(f"  f {f:<4g} " + "  ".join(f"K{K:g} {rec[K][f][0]:+.5f}" for K in KS)
                + f"  | argmin K {min(KS, key=lambda k: rec[k][f][0]):g}")
    r = res["d20"]
    pref = True
    for f in (0.5, 0.35):
        kb = min(KS, key=lambda k: r[k][f][0])
        pref &= (kb > 1.0) and (r[kb][f][0] < r[1.0][f][0] - 0.002)
    say(f"\nDECISION (fixed before the run): clean {clean}, prefers a dose above 1.0 (d20) {pref} -> "
        f"{'Jurkat CONFIRMS K562: October dose -> preregistered K 1.0 / K 2.0 pair once capacity delivers K 2' if (clean and pref) else 'no confirmation: K 1.0 stays the October dose'}")
    (HERE / "work/dose633j.out").write_text("\n".join(lines) + "\n")
    (HERE / "work/dose633j.json").write_text(json.dumps({"clean": clean, "depth": depth, "dn_block": dn_block,
                                                          "up_low": up_low, "prefers_higher": bool(pref), "rows": {
        d: {str(K): {str(f): v for f, v in fr.items()} for K, fr in rr.items()} for d, rr in res.items()}},
        indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
