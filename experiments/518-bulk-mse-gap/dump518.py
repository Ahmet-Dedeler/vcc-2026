"""518 (PLOLO) — one heavy pass per fold: dump everything a per-target KEY arm needs, so every later arm runs without the
whitening operator (<~2 GB instead of ~5.5 GB).

Exactly arms540.py's setup (lib420, 420's path): the fold's unkeyed PP and stack compositions (fPP, fST; panel zeroed), the
LODO-V1 key vector kv1 and its carried mask (KA: weights(kv1, 0.3, 3) == 418's saved `V1__PP`), the scored columns, the
truth, panel positions. A key arm is then  key_scale(fST, w) -> sharpen q 1.75 per frame (row mode on the frame's gate)
-> render (minnorm) -> realise at KMAX 16 against PP's unkeyed carriers (P_DE from fPP).
KNOWN ANSWER (abort unless |d| < 1e-12, run here before saving): recomputing V1_0.5_3 and V1_0.5_6 from the dumped arrays
through lib420's own functions reproduces arms540-<fold>.json's realised cos_agg and pds.

    uv run python experiments/518-bulk-mse-gap/dump518.py <fold>        (~5.5 GB peak, like arms540; one fold at a time)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import lib420  # noqa: E402
from lib420 import (COMPS, E062, FOLDS, PROD_OPS, Donors, Operator, betas_for, build_fold,  # noqa: E402
                    cd4_frozen_and_partial, compose_arm, compose_pp, deflate, key_scale, pre_normrel, source_residual,
                    union_rows, weights)
import lodo415  # noqa: E402
from comp420 import score_rows  # noqa: E402

W418 = ROOT / "experiments/418-bulk-composition/work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
WORK = HERE / "work"
TOL = 1e-12


def run(held: str) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    fold = build_fold(held)
    scored = np.asarray(fold["scored"])
    pp_pos = np.asarray(fold["panel_positions"])
    wcols = np.asarray(fold["wcols"])
    truth_full = np.asarray(fold["truth"], dtype=np.float64)
    pre_cached = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    sources = list(FOLDS[held]["sources"])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production" else FOLDS[held]["op"])
    nonpanel = np.ones(truth_full.shape[1], bool)
    nonpanel[pp_pos] = False
    betas = betas_for(held, nonpanel)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in lodo415.UNION])
    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, _ = deflate(frozen, 5)
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    fb = pre228.copy()
    pre = pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fb)
    f415, _, _, radd415, _ = compose_pp(op, pre, cd4, partial, h1, wcols, pp_pos)
    whit = op.apply(pre * lib420.ALPHA, wcols, lib420.GAMMA, lib420.LAM)
    fPP, pPP = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
    assert np.array_equal(fPP, f415) and np.array_equal(pPP["radd"], radd415), "KA-PP compose"
    mc, gh = COMPS["stack"]
    fST, pST = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=mc, g_h1=gh)
    import lib418 as _l418
    d_rows = {s_: _l418.donor_rows(s_, targets)[0] for s_ in sources}
    A_ = {s_: dn.A[s_] for s_ in sources}
    carried_all = sum(A_[s_].astype(int) for s_ in sources) >= 2
    bases = _l418.lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
    kv1 = np.nan_to_num(_l418.lodo_key(fPP, pPP, bases, d_rows, A_, carried_all, sources, pp_pos), nan=0.0)
    assert np.array_equal(weights(kv1, carried_all, 0.3, 3.0)[0], np.load(W418 / f"arrays418-{held}.npz")["V1__PP"]), "KA V1"
    # per-source per-target features the key does not use directly (row norms on the scored axis, whitened-base agreement)
    src_norm = np.stack([np.linalg.norm(d_rows[s_][:, scored], axis=1) for s_ in sources])
    src_avail = np.stack([np.asarray(A_[s_], bool) for s_ in sources])
    whit_norm = np.linalg.norm(whit[:, scored], axis=1)
    radd_norm = np.linalg.norm(pST["radd"][:, scored], axis=1)
    del bases, whit, f415, pre228, fb
    log("KA-PP and KA V1 PASS; composing the known-answer rows from the dump")

    # ---- KA: V1_0.5_3 and V1_0.5_6 from the dumped arrays == arms540's realised rows
    import rows518
    frames = rows518.challenge_frames()
    P_DE = rows518.de_profiles(fPP, frames)
    ref = json.loads((W540 / f"arms540-{held}.json").read_text())["arms"]
    ka = {}
    for lo, hi in ((0.5, 3.0), (0.5, 6.0)):
        w = weights(kv1, carried_all, lo, hi)[0]
        rea = rows518.realised_rows(key_scale(fST, w)[0], frames, P_DE, q=1.75)
        sr, _ = score_rows(rea, scored, truth_full[:, scored], betas, nonpanel)
        tag = f"V1_{lo:g}_{hi:g}"
        want = ref[tag]["realised"]["16"]
        d = max(abs(sr[k] - want[k]) for k in ("cos_agg", "pds", "pds_lam0.5"))
        ka[tag] = d
        log(f"KA {tag}: cos_agg {sr['cos_agg']:.12f} vs arms540 {want['cos_agg']:.12f}; max |d| {d:.1e}")
        assert d < TOL, f"KA FAILED {held} {tag}: {d}"
    np.savez_compressed(WORK / f"dump518-{held}.npz", targets=np.asarray(targets), sources=np.asarray(sources),
                        fPP=fPP, fST=fST, kv1=kv1, carried=carried_all, scored=scored, panel_positions=pp_pos,
                        truth=truth_full.astype(np.float32), src_norm=src_norm, src_avail=src_avail, whit_norm=whit_norm,
                        radd_norm=radd_norm, **{f"beta_{c}": betas[c] for c in "ABC"})
    (WORK / f"dump518-{held}.json").write_text(json.dumps(dict(fold=held, n=len(targets), sources=sources, ka=ka,
                                                               peak_rss_gib=E062.peak_rss_gib()), indent=1, default=float))
    log(f"dumped; peak RSS {E062.peak_rss_gib():.2f} GiB")


if __name__ == "__main__":
    run(sys.argv[1])
