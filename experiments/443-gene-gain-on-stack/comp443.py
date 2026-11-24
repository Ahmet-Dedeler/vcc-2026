"""443 — 228's per-gene cross-line reliability gain (410d's bulk arm) re-priced ON THE CURRENT BULK: the stack composition
(H1 x2 + CD4 x1.5), the 410 key K3 (clip .3,3), sharpening q 1.25 / 1.5, on 420's PLOLO instrument (design read-back and
414's ballast realisation at KMAX 16/24/32). Builder order is reproduced: key -> gene gain (non-panel Frobenius preserve,
= `apply_bulk_gene_gain` without render gating) -> sharpened render.

The gain is the fold's own 228 artifact (`228/work/gain-<fold>.npz`, gain_reg; fitted WITHOUT the held line), transformed
m = 1 + t (g/median - 1) (228's `multiplier`), optionally clipped to [lo, hi] (414: bounded = the safe variant; the
unbounded m reaches 19 and blew up 410d's context C).

KNOWN ANSWERS FIRST (abort on failure): `stack|K3+sh1.25` and `stack|K3+sh1.5` must equal 420's comp420-<fold>.json rows
(design and every realised KMAX, every field cmp() reads) to 1e-12.

    uv run python comp443.py <fold>     # k562_gwps | orion_hct116 | orion_hek293t | kolf2
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))
sys.path.insert(0, str(ROOT / "experiments/228-theta-gain-plolo"))

from lib420 import (ALPHA, COMPS, E062, FOLDS, GAMMA, KMAXES, LAM, PROD_OPS, UNION, WORK as W420, Donors, Operator,  # noqa: E402
                    betas_for, build_fold, cd4_frozen_and_partial, compose_arm, de_profiles, deflate, design_and_realised,
                    key_scale, pre_normrel, shared_energy_key, source_residual, union_rows, weights)
from comp420 import cmp, score_rows  # noqa: E402
from arms228 import multiplier  # noqa: E402

TOL = 1e-12
OUT = HERE / "work"


def frob_preserve(E, M):
    return M * (np.linalg.norm(E) / max(np.linalg.norm(M), 1e-300))


def run(held: str) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    fold = build_fold(held)
    fold["scored"] = np.asarray(fold["scored"])
    scored = fold["scored"]
    pp_pos = fold["panel_positions"]
    wcols = np.asarray(fold["wcols"])
    truth_full = np.asarray(fold["truth"], dtype=np.float64)
    truth = truth_full[:, scored]
    pre_cached = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    sources = list(FOLDS[held]["sources"])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production" else FOLDS[held]["op"])
    nonpanel = np.ones(truth_full.shape[1], bool)
    nonpanel[pp_pos] = False
    betas = betas_for(held, nonpanel)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, _ = deflate(frozen, 5)
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    pre = pre_normrel(dn, sources, "cap", -1.5, fallback_rows=pre228.copy())
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    finals = {}
    finals["PP"], _ = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
    m_cd4, g_h1 = COMPS["stack"]
    finals["stack"], _ = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=m_cd4, g_h1=g_h1)
    key410, carried410, _ = shared_energy_key(targets, sources, op, wcols, cols=None)
    wK3 = weights(key410, carried410, 0.3, 3.0)[0]
    P_DE = de_profiles(finals["PP"])
    z = np.load(ROOT / f"experiments/228-theta-gain-plolo/work/gain-{held}.npz", allow_pickle=True)
    live, greg = np.asarray(z["live"], bool), np.asarray(z["gain_reg"], float)
    assert greg.shape[0] == truth_full.shape[1], (greg.shape, truth_full.shape)
    ref = json.loads((W420 / f"comp420-{held}.json").read_text())["arms"]
    res = dict(fold=held, n_targets=len(targets), ka={}, arms={})

    def rec(name, bulk, q, check_ref=None):
        dr = design_and_realised(bulk, P_DE, q=q, kmaxes=KMAXES)
        sd, _ = score_rows({c: dr[c]["design"] for c in "ABC"}, scored, truth, betas, nonpanel)
        entry = dict(q=q, design=sd, realised={})
        for k in KMAXES:
            sr, _ = score_rows({c: dr[c]["realised"][k] for c in "ABC"}, scored, truth, betas, nonpanel)
            K = np.concatenate([dr[c]["stats"][k]["K"] for c in "ABC"])
            sr["frac_capped"] = float(np.mean(K >= k))
            sr["K_mean"] = float(np.mean(K))
            entry["realised"][str(k)] = sr
        if check_ref is not None:
            r = ref[check_ref]
            d = cmp(sd, r["design"])
            for k in entry["realised"]:
                d = max(d, cmp(entry["realised"][k], r["realised"][k]))
            res["ka"][name] = d
            log(f"KA {name} vs comp420 {check_ref}: worst |d| {d:.1e}")
            assert d < TOL, f"KNOWN ANSWER FAILED {held} {name}: {d}"
        res["arms"][name] = entry
        rr = entry["realised"]["16"]
        log(f"{name:40s} design cos {sd['cos_agg']:.6f} pds {sd['pds']:.6f} pds@.5 {sd['pds_lam0.5']:.6f} "
            f"| R16 cos {rr['cos_agg']:.6f} pds@.5 {rr['pds_lam0.5']:.6f} capped {100 * rr['frac_capped']:.0f}% "
            f"| A {np.mean([rr['AB'][c]['A'] for c in 'ABC']):.4f}")
        (OUT / f"comp443-{held}.json").write_text(json.dumps(res, indent=1, default=float))
        return entry

    base = key_scale(finals["stack"], wK3)[0]
    rec("stack|K3+sh1.25", base, 1.25, check_ref="stack|K3+sh1.25")
    rec("stack|K3+sh1.5", base, 1.5, check_ref="stack|K3+sh1.5")
    log("ALL KNOWN ANSWERS PASS")
    variants = {
        "gg0.25": (0.25, None, None),
        "gg0.25c0.5-2": (0.25, 0.5, 2.0),
        "gg0.25c0.67-1.5": (0.25, 2 / 3, 1.5),
        "gg0.125c0.5-2": (0.125, 0.5, 2.0),
        "gg0.5c0.5-2": (0.5, 0.5, 2.0),
    }
    ms = {}
    for tag, (t, lo, hi) in variants.items():
        m = multiplier(greg, live, t)
        if lo is not None:
            m = np.clip(m, lo, hi)
        ms[tag] = m
        res.setdefault("m_stats", {})[tag] = dict(q=np.quantile(m[live], [0, .01, .05, .5, .95, .99, 1]).tolist(),
                                                   frac_clip_lo=float(np.mean(m[live] <= (lo or -1e9) + 1e-12)),
                                                   frac_clip_hi=float(np.mean(m[live] >= (hi or 1e9) - 1e-12)))
        bulk = frob_preserve(base, base * m[None, :])
        rec(f"stack|K3+{tag}+sh1.5", bulk, 1.5)
        if tag == "gg0.25c0.5-2":
            rec(f"stack|K3+{tag}+sh1.25", bulk, 1.25)
    # CONTROL: the bounded multiplier permuted over the live genes (same distribution, alignment destroyed) -- must lose
    rng = np.random.default_rng(443)
    mp = ms["gg0.25c0.5-2"].copy()
    idx = np.flatnonzero(live)
    mp[idx] = mp[idx][rng.permutation(idx.size)]
    rec("CTRL:stack|K3+gg0.25c0.5-2perm+sh1.5", frob_preserve(base, base * mp[None, :]), 1.5)
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (OUT / f"comp443-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    run(sys.argv[1])
