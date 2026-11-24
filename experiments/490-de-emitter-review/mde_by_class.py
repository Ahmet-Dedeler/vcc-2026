#!/usr/bin/env python3
"""490: where does the MDE add pay, and where does it overshoot? First-order, exact on the lfc channel.

The MDE add (c 0.1, cap 0.1, coarse) is applied to every eligible gated gene in the sign of its base target. For each
gene class, remove that class's MDE move from the delivered mean (in linear space over all 400 cells, via the exact
reconstruction: delivered + (n_car/400)(A_noMDE - A_int) on that class only) and re-score nmae with calls and order held.
The other classes keep their MDE and their fu/fd (first-order: the mass balance is not re-solved). A positive delta means
the MDE HELPS that class (removing it costs nmae). Classes: head (top 400 |key|), block protected, block flipped, block
kept (prior agrees with E0), complement.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import recon_h1 as RC  # noqa: E402
import score_de  # noqa: E402

ROOT = HERE.parents[1]


def a_nomde(o, Lk_unused=None):
    """The carriers' mean without the MDE move: block genes at A0*2^(lam*move); non-block genes at the massprop A."""
    return o["A_nomde"]


def main():
    out = {}
    for seed in (0, 1, 2):
        tag = f"base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
        diag = json.loads((RC.H1 / f"diag-{tag}.json").read_text())
        st = np.load(RC.H1 / f"stats-{tag}.npz", allow_pickle=False)
        zr = np.load(RC.REF, allow_pickle=False)
        r = {k: zr[k] for k in ("genes", "targets", "keep", "ref_lin", "emit_ref_lin", "libraries", "b_ctrl")}
        zr.close()
        keep = np.asarray(r["keep"], dtype=bool)
        ref_emit = np.asarray(r["emit_ref_lin"], dtype=np.float64)
        ref_score = np.asarray(r["ref_lin"], dtype=np.float64)
        with RC.R.model_env("uniform"):
            b = RC.R.load_module("builder490m", ROOT / "scripts/make_submission_pm_cd4.py")
        with np.load(RC.EFF, allow_pickle=False) as z:
            eff = RC.R.map_axis(np.asarray(z["cap15"], dtype=np.float64), z["genes"], r["genes"])
        prior_g = RC.prior_signs_full(r["genes"])[keep]
        cnt = ref_score * float(np.median(r["libraries"])) / 1e6
        madd = np.where(cnt >= RC.MDE_MIN, np.minimum(RC.MDE_C / np.sqrt(np.maximum(cnt, 1e-3)), RC.MDE_CAP), 0.0)
        madd0 = np.zeros_like(madd)
        pick_rng, _ = b.make_rngs(seed)
        T = r["targets"].astype(str)
        G = np.asarray(r["genes"])[keep].astype(str)
        lfc_del = st["lfc"]
        m_del = (ref_score + 1e-9) * np.power(2.0, lfc_del) - 1e-9
        dA = np.zeros_like(lfc_del)          # (n_car/400) * (A_noMDE - A_int), per gene
        cls = {k: np.zeros(lfc_del.shape, bool) for k in ("head", "protected", "flipped", "kept", "complement")}
        for i, t in enumerate(T):
            picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
            o = RC.reconstruct(i, diag[i], eff[i], r, keep, ref_emit, ref_score, prior_g, madd, b, picks)
            dg0 = dict(diag[i]); dg0["mde_lambda"] = 0.0           # same prior dose, no MDE move at all
            z0 = RC.reconstruct(i, dg0, eff[i], r, keep, ref_emit, ref_score, prior_g, madd, b, picks)
            dA[i] = (o["n_car"] / 400.0) * (z0["A_int"] - o["A_int"])
            blk, prot, flip = o["blk"], o["protected"], o["flipped"]
            cls["head"][i] = o["head"]
            cls["protected"][i] = blk & prot
            cls["flipped"][i] = blk & flip & ~prot
            cls["kept"][i] = blk & ~flip & ~prot
            cls["complement"][i] = ~o["head"] & ~blk
        base, _ = score_de.score(T, G, st["p"], st["q"], lfc_del)
        res = {}
        for name, m in list(cls.items()) + [("all", np.ones(lfc_del.shape, bool))]:
            mm = np.maximum(m_del + np.where(m, dA, 0.0), 0.0)
            l2 = np.log2((mm + 1e-9) / (ref_score + 1e-9))
            raw, _ = score_de.score(T, G, st["p"], st["q"], l2)
            res[name] = round(raw["nmae"] - base["nmae"], 6)
        out[f"seed{seed}"] = res
        print(seed, res, flush=True)
    (HERE / "work" / "mde_by_class.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
