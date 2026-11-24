"""STEP 2/3 -- every donor-combination arm, measured as `W = X^2/A` and `sqrt(W/W_base)`.

One harness per process (the two Orion donor matrices are ~2.8 GB each, so holding both
harnesses' models at once is wasteful).  Results land in
`experiments/097-direction/work/<harness>-arms.json` for `w_report.py` to join.

Everything except the donor weighting is frozen at the shipped config: gamma 0.65 union
whitening from the CACHED whitener, q 1.0, alpha 1.08, min_cpm 5.0, gate_unexpressed, and the
emitter's per-perturbation `L` recovered from the shipped prediction.
"""
from __future__ import annotations

import json
import sys

import numpy as np
from wkit import BASE_DONORS, ROOT, Ctx, carriers, donor_vec, fixed, uniform


# ---------------------------------------------------------------------------- per-p key sets
def key_normrel(c) -> dict[str, np.ndarray]:
    """`||s_n(p)|| / median_p' ||s_n(p')||` -- how much this donor is carrying for `p`."""
    m = c.model
    out = {}
    for n in m.names:
        k = np.zeros(len(c.perts))
        for i, p in enumerate(c.perts):
            j = m.pert_index[n].get(p)
            if j is not None:
                k[i] = np.linalg.norm(m.s[n][j][m.dcols[n]]) / (m._scale[n] or 1.0)
        out[n] = k
    return out


def key_selfkd(c) -> dict[str, np.ndarray]:
    """On-target knockdown depth: `max(0, -s_n(p)[gene p])`, a pure QC signal, no reference."""
    m = c.model
    gpos = {g: i for i, g in enumerate(m.target_genes)}
    out = {}
    for n in m.names:
        k = np.zeros(len(c.perts))
        tmap = {int(t): int(d) for t, d in zip(m.tcols[n], m.dcols[n])}
        for i, p in enumerate(c.perts):
            j = m.pert_index[n].get(p)
            col = gpos.get(p)
            if j is None or col is None or col not in tmap:
                continue
            k[i] = max(0.0, -float(m.s[n][j][tmap[col]]))
        out[n] = k
    return out


def key_loo(c) -> dict[str, np.ndarray]:
    """`cos(s_n(p), mean of the OTHER carriers' s(p))`.  With two donors this is SYMMETRIC
    and therefore degenerate as a relative weight -- measured and reported, not assumed."""
    m = c.model
    out = {n: np.zeros(len(c.perts)) for n in m.names}
    for i, p in enumerate(c.perts):
        have = carriers(m, p)
        if len(have) < 2:
            continue
        V = {n: donor_vec(m, n, p) for n in have}
        for n in have:
            others = np.mean([V[o] for o in have if o != n], axis=0)
            a, b = V[n], others
            na, nb = np.linalg.norm(a), np.linalg.norm(b)
            out[n][i] = float(a @ b / (na * nb)) if na and nb else 0.0
    return out


def soft(c, keys, gam: float):
    """`w_n(p) propto key_n(p)^gam`, floored at the 10th percentile of the positive keys.

    The floor matters: a key that is exactly zero for a donor that DOES carry `p` (the
    on-target key hits zero whenever the donor shows no knockdown) would otherwise take all
    the weight at negative `gam`.  Flooring keeps the ordering and kills the blow-up.
    """
    idx = {p: i for i, p in enumerate(c.perts)}
    allk = np.concatenate([keys[n] for n in c.model.names])
    pos = allk[allk > 0]
    fl = float(np.percentile(pos, 10)) if pos.size else 1e-6

    def f(p):
        i = idx[p]
        return {n: float(max(keys[n][i], fl) ** gam) for n in c.model.names}
    f.floor = fl
    return f


def hard(c, keys):
    idx = {p: i for i, p in enumerate(c.perts)}

    def f(p):
        i = idx[p]
        best = max(c.model.names, key=lambda n: keys[n][i]
                   if p in c.model.pert_index[n] else -np.inf)
        return {n: (1.0 if n == best else 0.0) for n in c.model.names}
    return f


def dominant(c, wfun) -> list[str]:
    out = []
    for p in c.perts:
        w = wfun(p) or {n: 1.0 for n in c.model.names}
        have = carriers(c.model, p)
        out.append(max(have, key=lambda n: w.get(n, 0.0)) if have else "-")
    return out


# ------------------------------------------------------------------------------------ driver
def main():
    name = sys.argv[1]
    extra = tuple(a for a in sys.argv[2:] if a)
    c = Ctx(name, extra_donors=extra)
    dn = list(c.model.names)
    print(f"########## {name}   donors={dn}   n={len(c.perts)}  G={c.genes.size}")
    print(f"  wcols n={len(c.wcols)} sha256={c.wcols_sha[:16]}  (matched cached whitener)")
    print(f"  L per pert: median {int(np.median(c.Lvec))}  range {c.Lvec.min()}..{c.Lvec.max()}")
    cov = {n: int(sum(p in c.model.pert_index[n] for p in c.perts)) for n in dn}
    print(f"  donor coverage of the {len(c.perts)} targets: {cov}")

    base = c.run(uniform)
    Xr, Ar, Wr = c.XAW(c.dh_real)
    print(f"  base (uniform)  X={base['X']:.8g} A={base['A']:.8g} W={base['W']:.8g}")
    print(f"  scorer-realized X={Xr:.8g} A={Ar:.8g} W={Wr:.8g}  "
          f"(reconstruction sqrt-ratio {np.sqrt(base['W']/Wr):.6f})")
    Wb = base["W"]
    base_dom = dominant(c, uniform)
    rows = []

    def add(label, wfun, note=""):
        r = c.run(wfun)
        dom = dominant(c, wfun)
        nchg = sum(1 for a, b in zip(dom, base_dom) if a != b)
        cos = float(np.median(
            np.einsum("ij,ij->i", r["E"], base["E"])
            / np.maximum(np.linalg.norm(r["E"], axis=1)
                         * np.linalg.norm(base["E"], axis=1), 1e-30)))
        rec = dict(arm=label, X=r["X"], A=r["A"], W=r["W"],
                   mult=float(np.sqrt(r["W"] / Wb)), n_dom_changed=nchg,
                   med_cosE=cos, Enorm=r["Enorm"], dhnorm=r["dhnorm"], note=note)
        rows.append(rec)
        print(f"  {label:<38}X={r['X']:>11.6g} A={r['A']:>11.6g} W={r['W']:>11.6g} "
              f"mult={rec['mult']:>7.4f}  domchg={nchg:>4}  cos(E,Ebase)={cos:+.4f}")
        return rec

    print("\n--- (a1) uniform / sanity")
    add("a1 uniform (base)", uniform)

    print("\n--- (a2) single donor")
    for n in dn:
        w = {x: (1.0 if x == n else 0.0) for x in dn}
        miss = sum(1 for p in c.perts if p not in c.model.pert_index[n])
        add(f"a2 {n} alone (strict)", fixed(w), f"{miss} rows left at zero effect")

        def prefer(p, n=n, w=w):
            return dict(w) if p in c.model.pert_index[n] else None
        add(f"a2 {n} first, else others", prefer, "uncovered rows fall back to the rest")

    print("\n--- (a3) fixed non-uniform global weights (k562 share)")
    if set(BASE_DONORS) <= set(dn) and len(dn) == 2:
        for s in (0.0, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0):
            add(f"a3 k562 share {s:.2f}",
                fixed({"k562_gwps": s, "orion_hek293t": 1.0 - s}))

    print("\n--- (a4) basal weighting, exp(10*sim) as multi_donor.effect() implements it")
    sims = c.model.basal_weights(c.b_ctrl)
    ex = {n: float(np.exp(10.0 * sims[n])) for n in dn}
    tot = sum(ex.values())
    print(f"  basal cosines: { {k: round(v, 5) for k, v in sims.items()} }")
    print(f"  sharpened softmax weights: { {k: round(v/tot, 5) for k, v in ex.items()} }")
    add("a4 basal exp(10*sim)", fixed(ex),
        json.dumps({k: round(v / tot, 5) for k, v in ex.items()}))

    print("\n--- (a6) per-perturbation donor weighting")
    KEYS = {"normrel": key_normrel(c), "selfkd": key_selfkd(c), "loo": key_loo(c)}
    zero = {kn: {n: int(np.sum((keys[n] == 0) & np.array([p in c.model.pert_index[n] for p in c.perts])))
                 for n in dn} for kn, keys in KEYS.items()}
    print(f"  carriers whose key is exactly zero (the floor cases): {zero}")
    for kn, keys in KEYS.items():
        arr = np.array([keys[n] for n in dn])
        spread = float(np.mean(arr.max(0) > arr.min(0) + 1e-12))
        print(f"  key '{kn}': per-donor medians "
              f"{ {n: round(float(np.median(keys[n])), 4) for n in dn} }  "
              f"| rows where donors differ: {100*spread:.1f}%  "
              f"| argmax split { {n: int(np.sum(np.argmax(arr, 0) == i)) for i, n in enumerate(dn)} }")
    for kn, keys in KEYS.items():
        for gam in (-4.0, -3.0, -2.0, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25,
                    0.25, 0.5, 1.0, 2.0, 4.0):
            add(f"a6 soft {kn}^{gam:+.1f}", soft(c, keys, gam))
        add(f"a6 hard argmax {kn}", hard(c, keys))

    out = dict(harness=name, donors=dn, n=len(c.perts), G=int(c.genes.size),
               coverage=cov, sumD=c.sumD, claim=c.claim, Kc=c.Kc,
               W_base_reconstructed=Wb, W_base_realized=Wr,
               X_base_realized=Xr, A_base_realized=Ar,
               basal_sims=sims, arms=rows)
    tag = name if not extra else f"{name}+{'+'.join(extra)}"
    (ROOT / f"experiments/097-direction/work/{tag}-arms.json").write_text(
        json.dumps(out, indent=2))
    print(f"\nwrote work/{tag}-arms.json")


if __name__ == "__main__":
    main()
