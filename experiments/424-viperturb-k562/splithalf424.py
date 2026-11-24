"""424 diagnostic: VIPerturb's own split-half reliability vs its agreement with K562-GWPS (why the donor does not pay).

Per target (panel 284 / strongest-300 K562-GWPS / random 700 shared), in the donor table's space (b = log1p(5e4 P / sum P)
on the full probe axis, evaluated on the K562-GWPS gene axis minus the 300 panel genes, each row's own target zeroed,
centred over the target set within each half, as 159's `resid`):
  r_par    cos(half 1, half 2), random cell halves, DISJOINT NO-TARGET halves as controls
  r_guide  cos(top sgRNA, other sgRNAs) (disjoint control halves)
  c_half   cos(VIP half, K562-GWPS), mean over the two halves;  c_full = cos(VIP full table row, K562-GWPS)
Under a common-signal model: r_par = rho_h (a half's signal-energy fraction), rho_V = 2 rho_h / (1 + rho_h) (Spearman-
Brown), c_half^2 / rho_h = kappa^2 rho_G (what VIPerturb's reproducible signal shares with K562-GWPS).

    uv run python splithalf424.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
TS = 5e4


def main():
    ts = json.loads((WORK / "split/targets.json").read_text())
    Z = {b: np.load(WORK / f"split/{b}.npz") for b in ("binA", "binB", "binC")}
    feats = list(Z["binA"]["features"])
    labels = list(Z["binA"]["labels"])
    li = {l: k for k, l in enumerate(labels)}
    S = sum(Z[b]["sums"] for b in Z)            # each target's rows are nonzero in exactly one bin; controls in binA only
    N = sum(Z[b]["ncell"] for b in Z)
    panel = [l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
    gg = [l for l in (ROOT / "data/donors/k562_gwps/genes.txt").read_text().splitlines() if l]
    fi = {g: i for i, g in enumerate(feats)}
    ev = [g for g in gg if g in fi and g not in set(panel)]
    cols = np.asarray([fi[g] for g in ev])
    prof = lambda s: np.log1p(TS * s / s.sum())
    c0, c1 = prof(S[li["NO-TARGET|p0"]]), prof(S[li["NO-TARGET|p1"]])
    Dg = np.load(ROOT / "data/donors/k562_gwps/delta.npy", mmap_mode="r")
    gp = [l for l in (ROOT / "data/donors/k562_gwps/perts.txt").read_text().splitlines() if l]
    igp, igg = {p: i for i, p in enumerate(gp)}, {g: i for i, g in enumerate(gg)}
    gcols = np.asarray([igg[g] for g in ev])
    Dv = np.load(ROOT / "data/donors/viperturb_k562/delta.npy", mmap_mode="r")
    vp = [l for l in (ROOT / "data/donors/viperturb_k562/perts.txt").read_text().splitlines() if l]
    ivp = {p: i for i, p in enumerate(vp)}
    out = {}
    u = lambda M: M / np.maximum(np.linalg.norm(M, axis=1, keepdims=True), 1e-300)
    print(f"evaluation genes {len(ev)}; control halves {int(N[li['NO-TARGET|p0']])} / {int(N[li['NO-TARGET|p1']])} cells")
    for grp in ("panel", "strong", "random"):
        tg = [t for t in ts[grp] if N[li[f"{t}|p0"]] >= 10 and N[li[f"{t}|p1"]] >= 10]
        tgg = [t for t in tg if N[li[f"{t}|gA"]] >= 10 and N[li[f"{t}|gB"]] >= 10]
        epos = {g: i for i, g in enumerate(ev)}

        def rows(keys, ctrls):
            M = np.stack([(prof(S[li[k]]) - c)[cols] for k, c in zip(keys, ctrls)])
            for r, k in enumerate(keys):
                t = k.split("|")[0]
                if t in epos:
                    M[r, epos[t]] = 0.0
            return M - M.mean(0)

        H1, H2 = rows([f"{t}|p0" for t in tg], [c0] * len(tg)), rows([f"{t}|p1" for t in tg], [c1] * len(tg))
        GA, GB = rows([f"{t}|gA" for t in tgg], [c0] * len(tgg)), rows([f"{t}|gB" for t in tgg], [c1] * len(tgg))
        inG = [t for t in tg if t in igp]
        idx = [tg.index(t) for t in inG]
        Gm = np.stack([np.asarray(Dg[igp[t]], dtype=np.float64)[gcols] for t in inG])
        Vf = np.stack([np.asarray(Dv[ivp[t]], dtype=np.float64)[np.asarray([fi[g] for g in ev])] for t in inG])
        for r, t in enumerate(inG):
            if t in epos:
                Gm[r, epos[t]] = 0.0
                Vf[r, epos[t]] = 0.0
        Gm -= Gm.mean(0)
        Vf -= Vf.mean(0)
        r_par = (u(H1) * u(H2)).sum(1)
        r_gui = (u(GA) * u(GB)).sum(1)
        c_half = 0.5 * ((u(H1[idx]) * u(Gm)).sum(1) + (u(H2[idx]) * u(Gm)).sum(1))
        c_full = (u(Vf) * u(Gm)).sum(1)
        rh = float(r_par.mean())
        rho_v = 2 * rh / (1 + rh) if rh > 0 else float("nan")
        k2rg = float(c_half.mean() ** 2 / rh) if rh > 0 else float("nan")
        cells = np.asarray([N[li[f"{t}|p0"]] + N[li[f"{t}|p1"]] for t in tg])
        out[grp] = dict(n=len(tg), n_guide=len(tgg), cells_median=float(np.median(cells)), r_par_mean=rh,
                        r_par_median=float(np.median(r_par)), r_guide_mean=float(r_gui.mean()),
                        r_guide_median=float(np.median(r_gui)), c_half_mean=float(c_half.mean()),
                        c_full_mean=float(c_full.mean()), n_vs_gwps=len(inG), rho_V_spearman_brown=rho_v,
                        kappa2_rhoG_implied=k2rg)
        o = out[grp]
        print(f"{grp:7s} n {o['n']:4d} (cells median {o['cells_median']:.0f}): split-half r {o['r_par_mean']:+.4f} "
              f"(median {o['r_par_median']:+.4f}) -> rho_V {rho_v:.4f} | guide split r {o['r_guide_mean']:+.4f} (n {len(tgg)}) | "
              f"vs K562-GWPS (n {len(inG)}): half {o['c_half_mean']:+.4f}, full {o['c_full_mean']:+.4f} -> kappa^2 rho_G "
              f"implied {k2rg:.4f}")
    (WORK / "splithalf424.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
