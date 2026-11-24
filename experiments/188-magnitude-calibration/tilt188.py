"""188 — the ONE extra covariate the mandate allows: CONTROL EXPRESSION.

WHY THIS COVARIATE AND WHY NOW
------------------------------
`identity188.py` decomposes 169's `nmae` cost exactly.  For the h-space arm at q = 1.15 the L1
BUDGET spent on the reference's DE pool rises by `dU/D = +0.0383`; for the lfc-space arm it
rises by only `+0.0022` — a 17x difference that accounts for essentially the whole 5.4x gap
between the two arms' `nmae`.  The mechanism is structural: the reference's DE pool is
EXPRESSION-BIASED (a Wilcoxon has more power on well-expressed genes) and the builder's effect
space `log1p(TS*p)` is EXPRESSION-WEIGHTED, so any transform that concentrates magnitude in `h`
concentrates it exactly onto the scored pool.  The `h -> lfc -> h` reparameterisation is one
particular expression tilt; this script asks whether the tilt itself, as a free parameter, buys
anything the power exponent does not.

    tilt_theta:   h_g  ->  h_g * (cpm_g / median cpm)^(-theta)   on the gate, row L2 restored

`theta > 0` moves magnitude DOWN the expression axis (off the pool), `theta < 0` up it.
`theta = 0` is the identity.  Composed with the power `q` this is a two-parameter family that
contains 169's two arms approximately and the incumbent exactly.

HONEST WARNING, STATED BEFORE THE NUMBERS: a tilt is NOT monotone in `|h|`, so unlike the power
law it does NOT preserve the `|effect|` ordering, and the shipped `keep_top 8000` block moves.
That is a second variable (AGENTS.md rule 5) and the overlap is reported for every row.

Also priced: the L1-norm-preserving variant of the power arm, as a control on whether the
budget term `dU` alone is what `nmae` charges.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "experiments/169-sharpening-six-member"))
from budget188 import FOLDS, N_DE, PDS_SPAN, deliver, load, pds_of   # noqa: E402
from mapfit import to_lfc                                            # noqa: E402
from sharp import sharpen, sharpen_lfc                                # noqa: E402


def tilt(P, gate, cpm, theta, q=1.0):
    """`sign(h)|h|^q * (cpm/med)^-theta` on the gate, each row's own L2 norm restored."""
    G = np.where(gate[None, :], P, 0.0)
    w = np.ones(gate.size)
    med = float(np.median(cpm[gate]))
    w[gate] = (cpm[gate] / med) ** (-theta)
    out = np.sign(G) * np.abs(G) ** q * w[None, :]
    n0, n1 = np.linalg.norm(G, axis=1), np.linalg.norm(out, axis=1)
    return out * (n0 / np.maximum(n1, 1e-300))[:, None]


def sharpen_l1(P, q, gate):
    """The power arm with the row's L1 preserved instead of its L2 — the `dU` control."""
    G = np.where(gate[None, :], P, 0.0)
    out = np.sign(G) * np.abs(G) ** q
    n0, n1 = np.abs(G).sum(1), np.abs(out).sum(1)
    return out * (n0 / np.maximum(n1, 1e-300))[:, None]


def umd(P, R, pool):
    U, M, D = [], [], []
    for i in range(P.shape[0]):
        s = np.flatnonzero(pool[i])
        h, r = P[i, s], R[i, s]
        d = np.abs(r).sum()
        if d == 0 or s.size < 10:
            continue
        m = np.where(np.sign(h) == np.sign(r), np.minimum(np.abs(h), np.abs(r)), 0.0)
        U.append(np.abs(h).sum()); M.append(m.sum()); D.append(d)
    return np.array(U), np.array(M), np.array(D)


def main() -> None:
    z = np.load(os.path.join(WORK, "chalframes.npz"))
    cg = z["gA"] & z["gB"] & z["gC"]
    bA, gA = z["bA"], z["gA"]
    expr = np.expm1(bA)
    cpmA = 1e6 * expr / float(expr.sum())

    ARMS = [("identity", None)]
    ARMS += [(f"tilt{t:+.2f}", ("tilt", t, 1.0)) for t in (-0.20, -0.10, 0.10, 0.20, 0.35)]
    ARMS += [(f"q1.15_tilt{t:+.2f}", ("tilt", t, 1.15)) for t in (-0.10, 0.10, 0.20, 0.35)]
    ARMS += [("q1.15", ("pow", 1.15)), ("q1.15L", ("powL", 1.15)),
             ("q1.25L", ("powL", 1.25)), ("q1.10L", ("powL", 1.10)),
             ("q1.20L", ("powL", 1.20)), ("q1.15_L1norm", ("l1", 1.15))]

    acc = {}
    for held in FOLDS:
        F = load(held)
        mk = np.zeros(cg.size, bool)
        mk[F["scored"]] = True
        c = np.flatnonzero(mk & cg)
        sc = F["scored"]
        Tr = F["truth"][:, sc]
        LR = to_lfc(F["truth"], F["b"], gA)[:, c]
        pool = np.zeros(LR.shape, bool)
        k = min(N_DE, c.size)
        np.put_along_axis(pool, np.argpartition(-np.abs(LR), k - 1, axis=1)[:, :k], True, 1)
        P0 = F["pred"]
        base_top = None
        for name, spec in ARMS:
            if spec is None:
                Q = np.where(gA[None, :], P0, 0.0)
            elif spec[0] == "tilt":
                Q = tilt(P0, gA, cpmA, spec[1], spec[2])
            elif spec[0] == "pow":
                Q = sharpen(P0, spec[1], gA, mode="row")
            elif spec[0] == "powL":
                Q = sharpen_lfc(P0, spec[1], bA, gA, mode="row")
            else:
                Q = sharpen_l1(P0, spec[1], gA)
            rA = float(np.median(np.linalg.norm(Q, axis=1)) /
                       np.median(np.linalg.norm(np.where(gA[None, :], P0, 0.0), axis=1)))
            kt = min(8000, int(gA.sum()))
            top = set(np.argsort(-np.abs(Q[0]))[:kt].tolist())
            if name == "identity":
                base_top = top
            ov = len(base_top & top) / kt
            Dl = deliver(Q, bA, gA)
            LP = to_lfc(Dl, bA, gA)[:, c]
            U, M, Dd = umd(LP, LR, pool)
            acc.setdefault(name, []).append(dict(
                pds=pds_of(Dl[:, sc], Tr), U=U, M=M, D=Dd, rA=rA, ov=ov))
        print(f"  [{held}] done", flush=True)

    print(f"\n{'=' * 112}\nEXPRESSION TILT on PLOLO — 4-fold mean, challenge frame A, "
          f"delivered.  `nmae*` pool = top-2000 |lfc| per row.")
    print(f"  {'arm':<18}{'d pds raw':>11}{'d pds scaled':>13}{'nmae*':>9}{'d nmae*':>10}"
          f"{'dU/D':>9}{'dM/D':>9}{'e':>9}{'rA':>8}{'keeptop ov':>12}")
    rows = {}
    for name, _ in ARMS:
        v = acc[name]
        b = acc["identity"]
        U = np.concatenate([a["U"] for a in v]); M = np.concatenate([a["M"] for a in v])
        Dd = np.concatenate([a["D"] for a in v])
        U0 = np.concatenate([a["U"] for a in b]); M0 = np.concatenate([a["M"] for a in b])
        nm = float(np.mean(1 + (U - 2 * M) / Dd))
        nm0 = float(np.mean(1 + (U0 - 2 * M0) / Dd))
        du = float(np.mean((U - U0) / Dd)); dm = float(np.mean((M - M0) / Dd))
        dp = float(np.mean([a["pds"] - c["pds"] for a, c in zip(v, b)]))
        rows[name] = dict(dpds=dp, dpds_scaled=dp / PDS_SPAN, nmae=nm, dnmae=nm - nm0,
                          dU=du, dM=dm, rA=float(np.mean([a["rA"] for a in v])),
                          ov=float(np.mean([a["ov"] for a in v])),
                          per_fold=[a["pds"] - c["pds"] for a, c in zip(v, b)])
        print(f"  {name:<18}{dp:>+11.6f}{dp/PDS_SPAN:>+13.6f}{nm:>9.5f}{nm-nm0:>+10.5f}"
              f"{du:>+9.5f}{dm:>+9.5f}{(dm/du if abs(du) > 1e-12 else float('nan')):>9.3f}"
              f"{rows[name]['rA']:>8.5f}{rows[name]['ov']:>12.5f}")

    with open(os.path.join(WORK, "tilt188.json"), "w") as fh:
        json.dump(rows, fh, indent=1, default=float)
    print(f"\nwrote {WORK}/tilt188.json")


if __name__ == "__main__":
    sys.exit(main())
