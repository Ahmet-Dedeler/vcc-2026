"""156 — is the geometry we design the geometry the scorer reads?

THE QUESTION
------------
Production builds a log-space effect `D`, then renders a profile
    p = expm1(b_ctrl + D) / sum(expm1(b_ctrl + D))
and emits counts whose per-gene totals match `p` exactly (the quasi emitter's
exact-total repair, verified bit-identical on the board three times).

The scorer then reads
    delta_hat = log1p(TS * p) - b_ctrl.

That renormalisation is NOT the identity on `D`. Measured on the frozen HCT116
artifact: ||delta_hat - D|| is 22% of ||D||, and the residual's row-row cosine is
0.303 while the intended effects' is 0.0040. Whitening exists to make the rows
mutually dissimilar (hubness is what `pds` punishes); the renderer puts a shared
blob back.

Every effect-level screen in this project -- 017 through 155, every gamma sweep,
every donor comparison -- scored `D`. The board scores `delta_hat`. This script
scores both, on two independent held cell lines, in the scorer's own space.

ARMS
----
`intended`  : score D                      (what prior screens measured)
`delivered` : score render(D)              (what the board measures)
then a family of panel-space geometry transforms, each applied EITHER to D
(pre-render, where prior work put it) OR to the delivered effect with the
renderer inverted by fixed point (`*-post`), which is where it belongs if the
mechanism is real.

Reported endpoints, all on the same gene mask (all 18,533 genes minus the
surface's own target coordinates, matching `pds_cosine`'s panel exclusion):
  pds     exact rank statistic, ties at half
  cos     macro mean cos(pred_p, truth_p) -- the direction quality
  rr      mean off-diagonal row-row cosine of the PREDICTIONS -- the hubness load
  sse*    min over a global scalar of sum||g*pred-truth||^2 / sum||truth||^2,
          i.e. the uncorrected `mse` ratio at its own best amplitude. `mse`'s
          jackknife and rho terms are emitter properties held fixed across arms,
          so this is the part an arm can move. It is a PROXY, flagged as such.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

TS = 5.0e4
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
FROZEN = os.path.join(REPO, "experiments/104-context-heldout-credit-adapter/work/%s/frozen.npz")


# ---------------------------------------------------------------- the renderer
def render(D: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Exactly what the builder + scorer do, end to end, to a log-space effect."""
    lin = np.clip(np.expm1(b[None, :] + D), 0.0, None)
    p = lin / lin.sum(1, keepdims=True)
    return np.log1p(TS * p) - b[None, :]


def invert_render(H: np.ndarray, b: np.ndarray, iters: int = 40) -> np.ndarray:
    """Find D with render(D) ~= H, by fixed point. The map is near-identity, so
    D <- D + (H - render(D)) converges geometrically. Returns D."""
    D = H.copy()
    for _ in range(iters):
        r = H - render(D, b)
        if np.abs(r).max() < 1e-9:
            break
        D = D + r
    return D


# ------------------------------------------------------------------- endpoints
def _unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def endpoints(P: np.ndarray, T: np.ndarray) -> dict:
    A, B = _unit(P), _unit(T)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    pds = float(np.mean((beaten + 0.5 * ties) / (n - 1)))

    GG = A @ A.T
    rr = float(GG[~np.eye(n, dtype=bool)].mean())

    # min over a global scalar g of sum||g P - T||^2 / sum||T||^2
    a = float((P * P).sum())
    x = float((P * T).sum())
    d = float((T * T).sum())
    sse = (1.0 - x * x / (a * d)) if a > 0 else 1.0
    return {"pds": pds, "cos": float(own.mean()), "rr": rr, "sse": sse,
            "gstar": (x / a) if a > 0 else 0.0}


# ------------------------------------------------------------------ transforms
def restore_norms(Q: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Put each row of Q back on the row norm it had in P (every competitor that
    reports this family does it; it keeps alpha's meaning and leaves nmae's units
    alone)."""
    nq = np.maximum(np.linalg.norm(Q, axis=1, keepdims=True), 1e-300)
    npp = np.linalg.norm(P, axis=1, keepdims=True)
    return Q * (npp / nq)


def t_identity(P, **kw):
    return P


def t_centroid(P, dose=1.0, **kw):
    """AgingLab_BI: L <- L - dose * mean_t(L), row norms restored.
    Centering the queries makes every candidate's MEAN similarity to the query
    set exactly zero, which is first-order hubness removal and needs no truth."""
    return restore_norms(P - dose * P.mean(0, keepdims=True), P)


def t_rank(P, k=1, dose=1.0, **kw):
    """Heisenberg: strip the top-k right singular directions of the panel's own
    predictions, per context, preserving row norms."""
    Pc = P
    U, s, Vt = np.linalg.svd(Pc, full_matrices=False)
    V = Vt[:k]
    return restore_norms(Pc - dose * (Pc @ V.T) @ V, P)


def t_rowgram(P, lam=0.01, center=True, **kw):
    """Zhou Li 'Bevel': unit-normalise the rows, remove the panel centroid, apply
    a Tikhonov inverse square root of the ROW Gram (n x n, so this is cheap),
    then restore original row norms. lam -> 0 orthonormalises the rows."""
    R = _unit(P)
    if center:
        R = R - R.mean(0, keepdims=True)
    G = R @ R.T
    ridge = lam * float(np.mean(np.diag(G)))
    w, Q = np.linalg.eigh(G + ridge * np.eye(len(G)))
    w = np.maximum(w, 1e-12)
    W = (Q * w ** -0.5) @ Q.T
    return restore_norms(W @ R, P)


ARMS = [("identity", t_identity, {})]
for _d in (0.25, 0.5, 0.75, 1.0):
    ARMS.append((f"centroid d{_d}", t_centroid, {"dose": _d}))
for _k in (1, 2, 4, 8):
    ARMS.append((f"rank{_k} d1.0", t_rank, {"k": _k, "dose": 1.0}))
ARMS.append(("rank1 d0.5", t_rank, {"k": 1, "dose": 0.5}))
for _l in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0):
    ARMS.append((f"rowgram lam{_l}", t_rowgram, {"lam": _l}))


# ------------------------------------------------------------------------ main
def load(surface: str):
    z = np.load(FROZEN % surface, allow_pickle=True)
    t = np.load(os.path.join(WORK, f"truth-{surface}.npz"), allow_pickle=True)
    pg = np.asarray([str(x) for x in z["genes"]])
    tg = np.asarray([str(x) for x in t["genes"]])
    assert np.array_equal(pg, tg), "gene axes differ"
    ptar = np.asarray([str(x) for x in z["targets"]])
    ttar = np.asarray([str(x) for x in t["targets"]])
    order = {s: i for i, s in enumerate(ttar)}
    assert set(ptar) <= set(ttar), "prediction has targets the truth lacks"
    sel = np.asarray([order[s] for s in ptar])

    P = z["baseline"].astype(np.float64)
    T = t["effect"].astype(np.float64)[sel]
    b_build = z["bctrl"].astype(np.float64)
    b_ref = t["b_ctrl"].astype(np.float64)

    gpos = {g: i for i, g in enumerate(pg)}
    drop = {gpos[s] for s in ptar if s in gpos}
    mask = np.asarray([i for i in range(len(pg)) if i not in drop])
    return dict(P=P, T=T, b_build=b_build, b_ref=b_ref, mask=mask,
                targets=ptar, genes=pg)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(WORK, "screen.json"))
    args = ap.parse_args()

    results = {}
    for surface in ("hct116", "hek293t"):
        d = load(surface)
        P, T, mask = d["P"], d["T"], d["mask"]
        b_build, b_ref = d["b_build"], d["b_ref"]

        print(f"\n{'=' * 100}\n{surface}: {P.shape[0]} targets, "
              f"{len(mask)} scored genes ({P.shape[1] - len(mask)} target coords dropped)")
        print(f"  builder b_ctrl vs reference b_ctrl: max|diff| "
              f"{np.abs(b_build - b_ref).max():.3e}, "
              f"sum expm1(b_ref) = {np.expm1(b_ref).sum():.1f} (TS = {TS:.0f})")

        Pd = render(P, b_build)
        chk = invert_render(Pd, b_build)
        print(f"  render residual ||Pd-P||/||P|| median "
              f"{np.median(np.linalg.norm(Pd - P, axis=1) / np.linalg.norm(P, axis=1)):.4f}"
              f" | invert_render round-trip max err {np.abs(render(chk, b_build) - Pd).max():.2e}")

        rows = []
        print(f"\n  {'arm':<26} {'space':<10} {'pds':>8} {'d pds':>9} "
              f"{'cos':>8} {'row-row':>9} {'sse*':>8}")
        base = None
        for name, fn, kw in ARMS:
            for space in ("pre", "post"):
                if space == "pre":
                    Q = fn(P, **kw)
                    scored = render(Q, b_build)
                else:
                    Q = fn(Pd, **kw)
                    scored = render(invert_render(Q, b_build), b_build)
                e = endpoints(scored[:, mask], T[:, mask])
                if name == "identity":
                    if space == "post":
                        continue
                    base = e
                    label = "delivered"
                else:
                    label = space
                dp = e["pds"] - base["pds"]
                print(f"  {name:<26} {label:<10} {e['pds']:8.5f} {dp:+9.5f} "
                      f"{e['cos']:8.5f} {e['rr']:9.5f} {e['sse']:8.5f}")
                rows.append(dict(arm=name, space=label, **e, d_pds=dp))
            if name == "identity":
                e0 = endpoints(P[:, mask], T[:, mask])
                print(f"  {'(same effect, unrendered)':<26} {'INTENDED':<10} "
                      f"{e0['pds']:8.5f} {e0['pds'] - base['pds']:+9.5f} "
                      f"{e0['cos']:8.5f} {e0['rr']:9.5f} {e0['sse']:8.5f}"
                      f"   <-- what every prior screen scored")
                rows.append(dict(arm="identity", space="intended", **e0,
                                 d_pds=e0["pds"] - base["pds"]))
        results[surface] = rows

    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=1)
    print(f"\nwrote {args.out}")

    print(f"\n{'=' * 100}\nWORST-SURFACE RANKING (the repo's selection rule)")
    keys = {(r["arm"], r["space"]) for r in results["hct116"]}
    tab = []
    for k in keys:
        v = []
        for s in ("hct116", "hek293t"):
            m = [r for r in results[s] if (r["arm"], r["space"]) == k]
            v.append(m[0]["d_pds"] if m else np.nan)
        tab.append((min(v), v, k))
    tab.sort(reverse=True)
    print(f"  {'arm':<26} {'space':<10} {'worst d_pds':>12} {'hct116':>10} {'hek293t':>10}")
    for worst, v, k in tab:
        print(f"  {k[0]:<26} {k[1]:<10} {worst:+12.5f} {v[0]:+10.5f} {v[1]:+10.5f}")


if __name__ == "__main__":
    main()
