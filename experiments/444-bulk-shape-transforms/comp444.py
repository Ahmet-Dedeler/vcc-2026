"""444 — bulk SHAPE transforms other than 169's power law, on the current bulk (stack|K3), 420's PLOLO harness.

Sharpening sign(h)|h|^q (row norm restored on each context's render gate) is a crude posterior-mean shrinker: it
downweights small entries relative to large ones. If the true per-target effect is sparse (a few large genes, many ~0),
the posterior mean under a sparse prior is closer to SOFT THRESHOLDING, sign(h) max(|h| - tau, 0), than to a power law.
This screen prices soft thresholds (per-row tau at a quantile of the gated |h|), and power+threshold hybrids, against the
shipped q 1.5 and q 1.75 -- same rows, same realisation (414's ballast at KMAX 16/24/32), same beta-emulated pds.

KNOWN ANSWERS FIRST: power q 1.5 and q 1.75 through THIS script's transform path must equal comp420's
`stack|K3+sh1.5` / `stack|K3+sh1.75` rows to 1e-12 (design and every KMAX).

    uv run python comp444.py <fold>
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

import lib420  # noqa: E402
from lib420 import (ALPHA, COMPS, E062, FOLDS, GAMMA, KMAXES, LAM, PROD_OPS, UNION, WORK as W420, Donors, Operator,  # noqa: E402
                    betas_for, build_fold, cd4_frozen_and_partial, compose_arm, de_profiles, deflate, key_scale,
                    pre_normrel, shared_energy_key, source_residual, union_rows, weights)
from comp420 import cmp, score_rows  # noqa: E402

TOL = 1e-12
OUT = HERE / "work"
_power = lib420.sharpen


def transform(spec):
    """spec = ('pow', q) | ('soft', p) | ('powsoft', q, p) | ('hard', p); every variant restores the gated row norm."""
    def f(D, q_unused, gate, mode="row"):
        kind = spec[0]
        if kind == "pow":
            return _power(D, spec[1], gate, mode="row")
        A = np.atleast_2d(np.asarray(D, dtype=np.float64))
        G = np.where(np.atleast_1d(gate)[None, :], A, 0.0)
        if kind == "powsoft":
            G2 = _power(G, spec[1], gate, mode="row")
            p = spec[2]
        else:
            G2 = G
            p = spec[1]
        absG = np.abs(G2)
        tau = np.array([np.quantile(r[gate & (r != 0)], p) if np.any(gate & (r != 0)) else 0.0 for r in absG])
        if kind == "hard":
            out = np.where(absG >= tau[:, None], G2, 0.0)
        else:
            out = np.sign(G2) * np.maximum(absG - tau[:, None], 0.0)
        n0 = np.linalg.norm(G, axis=1)
        n1 = np.linalg.norm(out, axis=1)
        return out * (n0 / np.maximum(n1, 1e-300))[:, None]
    return f


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
    finalPP, _ = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
    m_cd4, g_h1 = COMPS["stack"]
    final, _ = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=m_cd4, g_h1=g_h1)
    key410, carried410, _ = shared_energy_key(targets, sources, op, wcols, cols=None)
    base = key_scale(final, weights(key410, carried410, 0.3, 3.0)[0])[0]
    P_DE = de_profiles(finalPP)
    ref = json.loads((W420 / f"comp420-{held}.json").read_text())["arms"]
    res = dict(fold=held, ka={}, arms={})

    def rec(name, spec, check_ref=None):
        lib420.sharpen = transform(spec)
        try:
            dr = lib420.design_and_realised(base, P_DE, q=1.0 if spec[0] != "pow" else spec[1], kmaxes=KMAXES)
        finally:
            lib420.sharpen = _power
        sd, _ = score_rows({c: dr[c]["design"] for c in "ABC"}, scored, truth, betas, nonpanel)
        entry = dict(spec=list(spec), design=sd, realised={})
        for k in KMAXES:
            sr, _ = score_rows({c: dr[c]["realised"][k] for c in "ABC"}, scored, truth, betas, nonpanel)
            K = np.concatenate([dr[c]["stats"][k]["K"] for c in "ABC"])
            sr["frac_capped"] = float(np.mean(K >= k))
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
        log(f"{name:22s} design cos {sd['cos_agg']:.6f} pds@.5 {sd['pds_lam0.5']:.6f} | R16 cos {rr['cos_agg']:.6f} "
            f"pds@.5 {rr['pds_lam0.5']:.6f} capped {100 * rr['frac_capped']:.0f}% A {np.mean([rr['AB'][c]['A'] for c in 'ABC']):.3f}")
        (OUT / f"comp444-{held}.json").write_text(json.dumps(res, indent=1, default=float))

    rec("pow1.5", ("pow", 1.5), check_ref="stack|K3+sh1.5")
    rec("pow1.75", ("pow", 1.75), check_ref="stack|K3+sh1.75")
    log("ALL KNOWN ANSWERS PASS")
    for p in (0.25, 0.5, 0.7, 0.85):
        rec(f"soft{p:g}", ("soft", p))
    for p in (0.5, 0.85):
        rec(f"hard{p:g}", ("hard", p))
    for q, p in ((1.25, 0.25), (1.25, 0.5), (1.5, 0.25), (1.5, 0.5)):
        rec(f"pow{q:g}+soft{p:g}", ("powsoft", q, p))
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (OUT / f"comp444-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    run(sys.argv[1])
