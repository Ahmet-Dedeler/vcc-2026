"""414 extra arms: controls for sharpening (q < 1 flattening; within-row permuted multipliers), the key floor (lo 0.4 /
0.5 / 0.7), sharpening dose continuation (q 1.35 / 1.5), and the final-row gain again. Same instrument as arms414.score (a23 read-back, challenge frames) PLUS the board-emulating pds of board414
(artifact beta added at lam 0.5 and 1.0). Known answer: `base`, `Kinf`, `Kinf+sh@1.15` must reproduce arms414 exactly.

    uv run python extra414.py
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arms414 import (ALPHA, CD4_F, E062, FOLDS, GAMMA, LAM, PROD_OPS, UNION, UNION_G, Operator, apply_gain,  # noqa: E402
                     build_fold, cd4_for, deliver, endpoints, frames, gene_gain_final, key_scale, merged, multiplier,
                     score, sharpen, shared_energy_key, solve_fraction, solve_fraction_carrying, source_residual,
                     union_rows, weights)

ROOT = HERE.parents[1]
WORK = HERE / "work"
LAMS = (0.5, 1.0)


def sharpen_perm(D, q, gate, seed):
    """CONTROL: the same per-entry multipliers |h|^(q-1) as sharpening, permuted within each row over its gated
    coordinates (keeps the multiplier distribution, destroys its alignment with |h|), row norm restored."""
    rng = np.random.default_rng(seed)
    G = np.where(gate[None, :], D, 0.0)
    out = G.copy()
    on = np.flatnonzero(gate)
    for i in range(G.shape[0]):
        h = G[i, on]
        nz = np.flatnonzero(h != 0)
        f = np.abs(h[nz]) ** (q - 1.0)
        v = h.copy()
        v[nz] = h[nz] * f[rng.permutation(nz.size)]
        n0, n1 = np.linalg.norm(h), np.linalg.norm(v)
        out[i, on] = v * (n0 / n1 if n1 > 0 else 1.0)
    return out


def run_fold(held):
    t0 = time.time()
    fold = build_fold(held)
    fold["scored"] = np.asarray(fold["scored"])
    wcols = np.asarray(fold["wcols"])
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    pre = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production"
                  else Path(FOLDS[held]["op"]))
    cd4 = cd4_for(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    z = np.load(ROOT / f"experiments/228-theta-gain-plolo/work/gain-{held}.npz", allow_pickle=True)
    live, greg = np.asarray(z["live"], bool), np.asarray(z["gain_reg"], float)
    whitened = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    f_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
    f_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
    base, _ = E062.add_post_whiten_residual(whitened, merged(cd4, f_cd4, h1, f_u), 1.0, wcols)
    base[:, fold["panel_positions"]] = 0.0
    key, carried, _ = shared_energy_key(targets, fold["sources"], op, wcols, cols=None)
    K = {f"K{lo:g}_{hi:g}": key_scale(base, weights(key, carried, lo, hi)[0])[0]
         for lo, hi in ((0.3, 1e9), (0.4, 1e9), (0.5, 1e9), (0.5, 3.0), (0.7, 1e9))}
    Ki = K["K0.3_1e+09"]
    ref = json.loads((WORK / f"arms414-{held}.json").read_text())["arms"]
    bz = np.load(WORK / "beta414.npz")
    dec = json.loads((WORK / "decomp414.json").read_text())
    nonpanel = np.ones(base.shape[1], bool)
    nonpanel[fold["panel_positions"]] = False
    betas = {}
    for ctx in "ABC":
        bb = np.asarray(bz[f"beta_{ctx}"], float).copy()
        bb[~nonpanel] = 0.0
        betas[ctx] = bb * np.sqrt(ref["base"]["AB"][ctx]["A"] / (dec["395"][ctx]["A"] * 18532 / 300))

    arms = {"base": (base, None), "Kinf": (Ki, None), "Kinf+sh@1.15": (Ki, ("h", 1.15)),
            "K0.4_inf": (K["K0.4_1e+09"], None), "K0.5_inf": (K["K0.5_1e+09"], None),
            "K0.5_3": (K["K0.5_3"], None), "K0.7_inf": (K["K0.7_1e+09"], None),
            "K0.5_inf+sh@1.15": (K["K0.5_1e+09"], ("h", 1.15)), "K0.5_inf+sh@1.25": (K["K0.5_1e+09"], ("h", 1.25)),
            "K0.5_3+sh@1.15": (K["K0.5_3"], ("h", 1.15)),
            "Kinf+sh@0.87": (Ki, ("h", 1 / 1.15)), "Kinf+shperm@1.15": (Ki, ("P", 1.15)),
            "Kinf+sh@1.35": (Ki, ("h", 1.35)), "Kinf+sh@1.5": (Ki, ("h", 1.5)),
            "Kinf+gf@0.25": (gene_gain_final(Ki, multiplier(greg, live, 0.25)), None)}
    out = dict(fold=held, arms={})
    for name, (E, sharp) in arms.items():
        if sharp is not None and sharp[0] == "P":
            # the permuted control is frame-specific (gate), so score it frame by frame here
            per, pl = [], {lam: [] for lam in LAMS}
            AB = {}
            for i, (ctx, (b, gate)) in enumerate(frames().items()):
                dl = deliver(sharpen_perm(E, sharp[1], gate, 4140 + i), b, gate, "minnorm")
                e = endpoints(dl[:, fold["scored"]], truth)
                e["cos_agg"] = float(np.sqrt(max(0.0, 1 - e["sse"])))
                per.append(e)
                H = dl[:, nonpanel]
                AB[ctx] = dict(A=float((H * H).sum(1).mean()))
                for lam in LAMS:
                    pl[lam].append(endpoints((dl + lam * betas[ctx][None, :])[:, fold["scored"]], truth)["pds"])
            s = dict(pds=float(np.mean([e["pds"] for e in per])), cos_agg=float(np.mean([e["cos_agg"] for e in per])),
                     AB=AB, **{f"pds_lam{lam:g}": float(np.mean(pl[lam])) for lam in LAMS})
        else:
            s, _ = score(E, fold, truth, sharp if sharp is None else (sharp[0], sharp[1]))
            pl = {lam: [] for lam in LAMS}
            for ctx, (b, gate) in frames().items():
                F = E if sharp is None else sharpen(E, sharp[1], gate, mode="row")
                dl = deliver(F, b, gate, "minnorm")
                for lam in LAMS:
                    pl[lam].append(endpoints((dl + lam * betas[ctx][None, :])[:, fold["scored"]], truth)["pds"])
            s.update({f"pds_lam{lam:g}": float(np.mean(pl[lam])) for lam in LAMS})
            s.pop("per_frame", None)
        if name in ref:
            d = max(abs(s["pds"] - ref[name]["pds"]), abs(s["cos_agg"] - ref[name]["cos_agg"]))
            assert d < 1e-12, f"KA FAILED {held} {name} {d}"
        out["arms"][name] = s
        print(f"[{held}] {name:18s} cos_agg {s['cos_agg']:.6f} pds {s['pds']:.5f} lam0.5 {s['pds_lam0.5']:.5f} "
              f"lam1 {s['pds_lam1']:.5f} A {np.mean([s['AB'][c]['A'] for c in 'ABC']):.3f} ({time.time() - t0:.0f}s)",
              flush=True)
    (WORK / f"extra414-{held}.json").write_text(json.dumps(out, indent=1))
    del op


def main():
    for held in sys.argv[1:] or ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"):
        if not (WORK / f"extra414-{held}.json").exists():
            run_fold(held)


if __name__ == "__main__":
    main()
