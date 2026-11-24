"""168d — break the HCT116 vs HEK293T tie on h1-test, the low-blob surface.

## The tie

168c fitted a non-parametric transform `w = (m/|D|)^tau` with `m(|D|) = E[t*sign(D) | |D|]`,
cross-surface (fit on one Orion line, score on the other). At `tau = 0.5` it is worth
+0.01845 / +0.01874 raw pds. The control that could fail half-failed:

| variant | HCT116 | HEK293T |
| --- | ---: | ---: |
| FITTED   | +0.01845 | +0.01874 |
| REVERSED | **-0.12838** | **+0.01712** |
| SHUFFLED | -0.02169 | -0.01434 |

`SHUFFLED` fails on both, so correct bin assignment matters and the arm is not a generic
reweighting artefact. But HEK293T is nearly indifferent to REVERSING the curve, and HCT116 is
not. I proposed HEK293T's 44.1% target-independent blob as the cause and **that is refuted**:
scoring against mean-removed truth leaves REVERSED at +0.01897 there (FITTED +0.02122) while
HCT116 stays at -0.13174.

So two surfaces disagree about whether the transform's DIRECTION matters, with no mechanism. A
transform whose optimal direction is surface-dependent will not transfer to D/E/F, which is the
stated priority. **h1-test breaks the tie**: 3.6% blob against 26.4% / 44.1%, it is the only
surface validated against a board `reach` point, and -- uniquely -- it can be scored through the
PRODUCTION operator without leakage, because it is not an Orion line and `union_pruned`'s basis
excludes it.

## Two design points

1. **The curve is indexed by |D| PERCENTILE, not absolute magnitude.** The Orion frozen
   artifacts ship alpha 0.982 and the h1-test artifact ships 1.08, so absolute bin edges would
   not align. Percentile indexing is scale-free, transfers across surfaces by construction, and
   is the form that would actually deploy.
2. **The h1-test base effect is RECOVERED from the emitted artifact**
   (`095/out-h1test/.../pred.h5ad`) as `log1p(TS*p) - b_ctrl`. That object is mass-consistent by
   construction, so re-rendering it after the mass fix is the identity and there is no
   double-render.

Leakage: `data/donors/h1_2025` covers 85/85 of h1-test's targets, but it is not used here at all
-- the prediction is the frozen 095 artifact built from `k562_gwps` + Orion, and the transform
curve comes from the Orion pair.
"""
from __future__ import annotations

import os
import sys

import h5py
import numpy as np
import scipy.sparse as sp

TS = 5.0e4
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
from massfix import GATE_CPM, mass_fix, production_effect      # noqa: E402
from screen import endpoints, render                           # noqa: E402

H1 = os.path.join(REPO, "data/harness/h1-test.h5ad")
PRED = os.path.join(REPO, "experiments/095-clean-harness/out-h1test/runs/"
                          "g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
NB = 24
CTRL = "non-targeting"


def labels(h, path):
    o = h[path]
    if hasattr(o, "keys") and "categories" in o:
        c = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                        for x in o["categories"][:]])
        return c[o["codes"][:]]
    v = o["values"][:] if hasattr(o, "keys") else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def group_profiles(path, chunk=4000):
    """Per-target summed counts -> (targets, genes, counts matrix)."""
    with h5py.File(path, "r") as h:
        n, g = (int(x) for x in h["X"].attrs["shape"])
        genes = labels(h, "var/_index")
        lab = labels(h, "obs/target_gene")
        uniq = np.asarray(sorted(set(lab)))
        gid = {u: i for i, u in enumerate(uniq)}
        codes = np.asarray([gid[x] for x in lab], dtype=np.int64)
        ip = h["X/indptr"][:]
        acc = np.zeros((len(uniq), g))
        for r0 in range(0, n, chunk):
            r1 = min(r0 + chunk, n)
            s, e = int(ip[r0]), int(ip[r1])
            X = sp.csr_matrix((h["X/data"][s:e].astype(np.float64), h["X/indices"][s:e],
                               ip[r0:r1 + 1] - ip[r0]), shape=(r1 - r0, g))
            oh = sp.csr_matrix((np.ones(r1 - r0), (codes[r0:r1], np.arange(r1 - r0))),
                               shape=(len(uniq), r1 - r0))
            acc += np.asarray((oh @ X).todense())
    return uniq, genes, acc


def pct_curve(D, T, gate, mask, nb=NB):
    """E[t*sign(D) | percentile of |D|], isotonic. Percentile index => scale-free."""
    Dm, Tm, gm = D[:, mask], T[:, mask], gate[mask]
    Dc = Dm - Dm.mean(0, keepdims=True)
    Tc = Tm - Tm.mean(0, keepdims=True)
    sel = np.zeros_like(Dm, dtype=bool)
    sel[:, gm] = True
    sel &= np.abs(Dm) > 0
    ap, al = np.abs(Dc[sel]), Tc[sel] * np.sign(Dc[sel])
    r = np.argsort(np.argsort(ap)) / max(len(ap) - 1, 1)          # percentile in [0,1]
    m = np.zeros(nb)
    for i in range(nb):
        s = (r >= i / nb) & (r < (i + 1) / nb) if i < nb - 1 else (r >= i / nb)
        m[i] = np.mean(al[s]) if s.sum() > 200 else np.nan
    ok = np.isfinite(m)
    m = np.interp(np.arange(nb), np.flatnonzero(ok), m[ok])
    m = np.maximum.accumulate(np.maximum(m, 0.0))
    xb = np.array([np.mean(ap[(r >= i / nb) & (r < (i + 1) / nb)]) if i < nb - 1
                   else np.mean(ap[r >= i / nb]) for i in range(nb)])
    mult = np.where(xb > 0, m / np.maximum(xb, 1e-300), 0.0)
    return mult / max(mult.max(), 1e-300)


def apply_pct(D, gate, mult, tau):
    """Per-row percentile-indexed reweighting, renormalised to the median row norm."""
    if tau == 0.0:
        return D.copy()
    G = np.where(gate[None, :], D, 0.0)
    out = np.zeros_like(G)
    nb = len(mult)
    for i in range(G.shape[0]):
        a = np.abs(G[i]); nz = np.flatnonzero(a > 0)
        if nz.size == 0:
            continue
        r = np.argsort(np.argsort(a[nz])) / max(nz.size - 1, 1)
        idx = np.clip((r * nb).astype(int), 0, nb - 1)
        out[i, nz] = G[i, nz] * mult[idx] ** tau
    bf = float(np.median(np.linalg.norm(G, axis=1)))
    af = float(np.median(np.linalg.norm(out, axis=1))) or 1.0
    return out * (bf / af)


def main() -> None:
    rng = np.random.default_rng(7)
    # ---- the curve, fitted on the Orion pair only (percentile-indexed, scale-free)
    curves = []
    for n in ("hct116", "hek293t"):
        D, T, b, bref, gate, mask, ptar = production_effect(n)
        curves.append(pct_curve(D, T, gate, mask))
    mult = np.mean(curves, axis=0)
    mult /= mult.max()
    print("curve fitted on the Orion pair, percentile-indexed (scale-free):")
    print("  " + "  ".join(f"{v:.3f}" for v in mult))
    print(f"  agreement between the two Orion fits: pearson "
          f"{np.corrcoef(curves[0], curves[1])[0, 1]:+.4f}\n")

    # ---- h1-test truth and the recovered production-shaped base effect
    tg, genes_t, acc_t = group_profiles(H1)
    ci = int(np.flatnonzero(tg == CTRL)[0])
    prof_t = np.log1p(TS * acc_t / acc_t.sum(1, keepdims=True))
    b_ctrl = prof_t[ci]
    keep_t = np.asarray([i for i in range(len(tg)) if i != ci])
    T = prof_t[keep_t] - b_ctrl[None, :]
    targets = tg[keep_t]

    tp, genes_p, acc_p = group_profiles(PRED)
    pos = {str(t): i for i, t in enumerate(tp)}
    assert np.array_equal(genes_t, genes_p), "h1-test and pred gene axes differ"
    sel = np.asarray([pos[str(t)] for t in targets if str(t) in pos])
    have = np.asarray([i for i, t in enumerate(targets) if str(t) in pos])
    prof_p = np.log1p(TS * acc_p[sel] / acc_p[sel].sum(1, keepdims=True))
    D = prof_p - b_ctrl[None, :]                     # mass-consistent by construction
    T = T[have]; targets = targets[have]

    expr = np.expm1(b_ctrl)
    cpm = 1e6 * expr / expr.sum()
    gate = (cpm > GATE_CPM) & (b_ctrl > 0)
    gpos = {g: i for i, g in enumerate(genes_t)}
    drop = {gpos[str(t)] for t in targets if str(t) in gpos}
    mask = np.asarray([i for i in range(len(genes_t)) if i not in drop])
    print(f"h1-test: {len(targets)} targets, {int(gate.sum())} gated genes, "
          f"{len(mask)} scored genes; median |D| row norm "
          f"{np.median(np.linalg.norm(D[:, mask], axis=1)):.4f}")

    variants = {"FITTED (Orion-fitted)": mult,
                "REVERSED": mult[::-1].copy(),
                "SHUFFLED": rng.permutation(mult)}
    for demean in (False, True):
        tag = "MEAN-REMOVED (diagnostic)" if demean else "RAW truth (as scored)"
        base_eff = mass_fix(D, b_ctrl, gate)
        P0 = render(base_eff, b_ctrl)[:, mask]; Q0 = T[:, mask]
        if demean:
            P0 = P0 - P0.mean(0, keepdims=True); Q0 = Q0 - Q0.mean(0, keepdims=True)
        base = endpoints(P0, Q0)
        print(f"\n########## h1-test, {tag}   baseline pds {base['pds']:.5f} "
              f"cos {base['cos']:.5f}")
        print(f"  {'variant':<24} {'tau':>5} {'pds':>9} {'d pds':>9} {'d cos %':>9}")
        for nm, mv in variants.items():
            for tau in (0.5, 1.0):
                eff = mass_fix(apply_pct(D, gate, mv, tau), b_ctrl, gate)
                P = render(eff, b_ctrl)[:, mask]; Q = T[:, mask]
                if demean:
                    P = P - P.mean(0, keepdims=True); Q = Q - Q.mean(0, keepdims=True)
                e = endpoints(P, Q)
                print(f"  {nm:<24} {tau:>5.2f} {e['pds']:9.5f} {e['pds'] - base['pds']:+9.5f} "
                      f"{100 * (e['cos'] / base['cos'] - 1):+8.1f}%")
    print("\nTIE-BREAK READING")
    print("  HCT116 said direction matters (fitted +0.018, reversed -0.128).")
    print("  HEK293T said it barely matters (fitted +0.019, reversed +0.017), and the blob")
    print("  explanation was refuted. Whichever h1-test agrees with is the answer, and if")
    print("  REVERSED wins here too the family is an artefact.")


if __name__ == "__main__":
    sys.exit(main())
