"""414: a board-emulating pds read on PLOLO — add the shipped artifact's NON-SCALING shared pattern beta (shared414.py)
to every delivered row, and ask whether PLOLO then retrodicts the board's pds cost of the key (410b - 395 = -0.003582 raw).

Why: PLOLO's read-back and H1 both missed 410's board pds (PLOLO -0.0006, H1 +0.0046, board -0.0036/-0.0038). The shipped
395 -> 410b rows show a pattern common to all rows that the key did NOT scale (|beta|^2 ~ 0.05 = 10% of a floor row's
energy; floor rows' cosine with beta 0.30 vs 0.06 in 395). A per-target scale shrinks the signal of 42% of targets to 0.19x
while that pattern stays, so those rows point at a common direction and lose retrieval rank. PLOLO's deliver() has the
render but not the emitter, so this script adds beta (rescaled to PLOLO's row energy per frame) to every row of every arm.

KNOWN ANSWERS: lam = 0 must reproduce arms414's pds exactly (1e-12); PLOLO's own render-only beta_r is reported so
the emulation is not double counting it.

    uv run python board414.py
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
import arms414 as M  # noqa: E402  (imports the whole 228/410/169 stack, sets plolo.WORK to ours)
from arms414 import (ALPHA, CD4_F, E062, FOLDS, GAMMA, LAM, PROD_OPS, UNION, UNION_G, Operator, apply_gain,  # noqa: E402
                     build_fold, cd4_for, deliver, endpoints, frames, gene_gain_final, key_scale, merged, multiplier,
                     sharpen, sharpen_lfc, shared_energy_key, solve_fraction, solve_fraction_carrying, source_residual,
                     union_rows, weights)

ROOT = HERE.parents[1]
WORK = HERE / "work"
LAMS = (0.0, 0.5, 1.0, 2.0)
ARMS = ("base", "K3", "Kinf", "Kinf+sh@1.05", "Kinf+sh@1.1", "Kinf+sh@1.15", "Kinf+sh@1.25", "Kinf+shL@1.15",
        "K3+sh@1.15", "Kinf+gf@0.1", "Kinf+gf@0.25", "gi@0.25+Kinf", "sh@1.15", "Kinf@x0.97")


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

    def compose(m_pre=None):
        whitened = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
        if m_pre is not None:
            whitened, _ = apply_gain(whitened, m_pre, wcols)
        f_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
        f_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
        final, _ = E062.add_post_whiten_residual(whitened, merged(cd4, f_cd4, h1, f_u), 1.0, wcols)
        final[:, fold["panel_positions"]] = 0.0
        return final

    base = compose()
    key, carried, _ = shared_energy_key(targets, fold["sources"], op, wcols, cols=None)
    W3, Wi = weights(key, carried, 0.3, 3.0)[0], weights(key, carried, 0.3, 1e9)[0]
    K3, w3 = key_scale(base, W3)
    Ki, _ = key_scale(base, Wi)
    bulks = {"base": (base, None), "K3": (K3, None), "Kinf": (Ki, None),
             "Kinf+sh@1.05": (Ki, ("h", 1.05)), "Kinf+sh@1.1": (Ki, ("h", 1.10)), "Kinf+sh@1.15": (Ki, ("h", 1.15)),
             "Kinf+sh@1.25": (Ki, ("h", 1.25)), "Kinf+shL@1.15": (Ki, ("L", 1.15)), "K3+sh@1.15": (K3, ("h", 1.15)),
             "Kinf+gf@0.1": (gene_gain_final(Ki, multiplier(greg, live, 0.1)), None),
             "Kinf+gf@0.25": (gene_gain_final(Ki, multiplier(greg, live, 0.25)), None),
             "gi@0.25+Kinf": (key_scale(compose(multiplier(greg, live, 0.25)), Wi)[0], None),
             "sh@1.15": (base, ("h", 1.15)), "Kinf@x0.97": (Ki * 0.97, None)}
    ref = json.loads((WORK / f"arms414-{held}.json").read_text())["arms"]
    bz = np.load(WORK / "beta414.npz")
    dec = json.loads((WORK / "decomp414.json").read_text())
    nonpanel = np.ones(base.shape[1], bool)
    nonpanel[fold["panel_positions"]] = False
    out = dict(fold=held, arms={}, beta={})
    delivered = {}
    for name in ARMS:
        E, sharp = bulks[name]
        per = {lam: [] for lam in LAMS}
        for ctx, (b, gate) in frames().items():
            F = E
            if sharp is not None:
                F = sharpen(E, sharp[1], gate, mode="row") if sharp[0] == "h" else sharpen_lfc(E, sharp[1], b, gate, mode="row")
            dl = deliver(F, b, gate, "minnorm")
            if name in ("base", "K3"):
                delivered[(name, ctx)] = dl
            beta = np.asarray(bz[f"beta_{ctx}"], float).copy()
            beta[~nonpanel] = 0.0
            # rescale beta to PLOLO's row energy: artifact mean row energy = A_c * 18532 / 300
            e_art = dec["395"][ctx]["A"] * 18532 / 300
            e_plo = ref["base"]["AB"][ctx]["A"]
            beta *= np.sqrt(e_plo / e_art)
            for lam in LAMS:
                e = endpoints((dl + lam * beta[None, :])[:, fold["scored"]], truth)
                per[lam].append(e)
        s = {f"lam{lam:g}": dict(pds=float(np.mean([e["pds"] for e in per[lam]])),
                                 cos_agg=float(np.mean([np.sqrt(max(0.0, 1 - e["sse"])) for e in per[lam]])))
             for lam in LAMS}
        d0 = abs(s["lam0"]["pds"] - ref[name]["pds"])
        assert d0 < 1e-12, f"KA lam=0 FAILED {held} {name}: {d0}"
        out["arms"][name] = s
        print(f"[{held}] {name:16s} " + " ".join(f"lam{lam:g} pds {s[f'lam{lam:g}']['pds']:.5f}" for lam in LAMS)
              + f"  ({time.time() - t0:.0f}s)", flush=True)
    # PLOLO's own render-only non-scaling common part (K3 vs base, same regression as shared414)
    for ctx in "ABC":
        D0, Db = delivered[("base", ctx)][:, nonpanel], delivered[("K3", ctx)][:, nonpanel]
        br = ((1 - w3)[:, None] * (Db - w3[:, None] * D0)).sum(0) / float(((1 - w3) ** 2).sum())
        rowE = np.einsum("ij,ij->i", Db, Db)
        bart = np.asarray(bz[f"beta_{ctx}"], float)[nonpanel]
        out["beta"][ctx] = dict(render_beta_energy=float(br @ br), render_beta_over_floor=float(br @ br) / float(
            np.median(rowE[w3 < w3.min() + 1e-9])), cos_render_vs_artifact=float(br @ bart / np.sqrt((br @ br) * (bart @ bart))))
        print(f"[{held}] {ctx}: PLOLO render-only beta |b|^2 {br @ br:.4f} ({out['beta'][ctx]['render_beta_over_floor']:.3f} x "
              f"floor row), cos with artifact beta {out['beta'][ctx]['cos_render_vs_artifact']:+.3f}", flush=True)
    (WORK / f"board414-{held}.json").write_text(json.dumps(out, indent=1))
    del op


def main():
    for held in sys.argv[1:] or ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"):
        if not (WORK / f"board414-{held}.json").exists():
            run_fold(held)


if __name__ == "__main__":
    main()
