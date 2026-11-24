"""204 Phase 1 — estimate each donor's attenuation `a_i` by routes that do not share a bias.

MODEL.  For a target t, donor i's measured effect is `d_i(t) = a_i s(t) + e_i(t)`, with
`s(t)` the cross-line-shared component, `a_i` the donor's attenuation (how much of the
true repression its guides actually achieved), and `e_i` line-specific biology plus
sampling.  `a_i` is not identifiable in isolation -- only ratios are -- so everything is
reported relative to `k562_gwps = 1`.

ROUTE A -- self-gene residual.  `a_i = 1 - r_i`, `r_i` the median achieved residual
expression of the target's own transcript (Phase 0, reproduced from 194 to 2e-7).
Measures the assay directly.  Its weakness is that it is a property of the target gene's
own row of the matrix, and 194 showed it inverts on Orion.

ROUTE B -- shared-target inner products.  `<d_i, d_j> = a_i a_j |s|^2` for i != j, because
the noise terms are independent ACROSS DONORS (different labs, different cells).  So
ratios of attenuation come out of the Gram matrix with `|s|^2` cancelled and with no
regression dilution:

    triple:      a_i / a_j = sqrt( G_ik G_jk^-1 * G_ij G_ij^-1 ) ... explicitly
                 a_i = sqrt(G_ij G_ik / G_jk)                      (3 donors, closed form)
    instrument:  a_i / a_j = G_iX / G_jX    for any FOURTH line X   (one product each)

The instrument form is the robust one: the triple form needs the Orion-Orion product
G(hct116, hek293t), and 194 measured that this exact product is inflated by a shared
X-Atlas technical axis.  An instrument outside both platforms never touches it.

ROUTE B' -- naive regression of a donor on the consensus of the others.  Reported because
it is the obvious thing to do and because its errors-in-variables bias (the regressor
carries noise) is the reason routes B and B' disagree in a predictable direction.

BLOB CONTROL.  Every route is computed twice: as-is, and after subtracting each donor's
own target-independent mean vector over the estimation set.  That vector is the
long-gene / nuclear-fraction technical axis 194 found knockdown weighting inflating
1.00% -> 10.45% on HEK293T; if it drives `a_i`, the number is measuring the platform.

    uv run python atten204.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))

from prep204 import DONORS, bank, prepared, shared_targets  # noqa: E402
from plolo import PANEL_SET  # noqa: E402

EFFECT_SCALE = 6.817
INSTRUMENTS = {
    "kolf2": ROOT / "data/external/kolf2",
    "rpe1": ROOT / "data/donors/rpe1",
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
    "h1_2025": ROOT / "data/donors/h1_2025",
    "k562_essential": ROOT / "data/donors/k562_essential",   # SAME LINE as k562_gwps
}


def gram(E: dict[str, np.ndarray], rows: np.ndarray) -> dict:
    """Sum of per-target inner products over `rows`, for every donor pair and self."""
    keys = list(E)
    out = {}
    for i, a in enumerate(keys):
        for b in keys[i:]:
            out[f"{a}|{b}"] = float((E[a][rows] * E[b][rows]).sum())
    return out


def triple_a(g: dict) -> dict[str, float]:
    """a_i from the 3-donor closed form, normalised to k562_gwps = 1."""
    k, h, e = DONORS

    def G(x, y):
        return g[f"{x}|{y}"] if f"{x}|{y}" in g else g[f"{y}|{x}"]
    u = {k: np.sqrt(G(k, h) * G(k, e) / G(h, e)),
         h: np.sqrt(G(h, k) * G(h, e) / G(k, e)),
         e: np.sqrt(G(e, k) * G(e, h) / G(k, h))}
    return {d: float(u[d] / u[k]) for d in DONORS}, u


def main() -> None:
    tg = shared_targets()
    panel_mask = np.array([t in PANEL_SET for t in tg])
    P = {d: prepared(d, tg, "shared") for d in DONORS}
    E = {d: P[d]["effects"] for d in DONORS}
    avail = {d: P[d]["target_available"] for d in DONORS}
    assert all(avail[d].all() for d in DONORS), "shared set must be complete"
    n_t = len(tg)

    out: dict = dict(n_targets=n_t, n_panel=int(panel_mask.sum()),
                     raw_median_norm={d: P[d]["raw_median_norm"] for d in DONORS},
                     mapped_gene_count={d: P[d]["mapped_gene_count"] for d in DONORS})

    # ---------------------------------------------------------------- route A: residual
    resid = {}
    for d in DONORS:
        r = np.load(WORK / f"resid-{d}.npy")
        names = (WORK / f"resid-{d}.txt").read_text().splitlines()
        pos = {p: i for i, p in enumerate(names)}
        idx = np.array([pos[t] for t in tg])
        resid[d] = r[idx]
    routeA = {}
    for scope, m in (("shared", np.ones(n_t, bool)), ("panel", panel_mask)):
        rec = {}
        for d in DONORS:
            v = resid[d][m]
            v = v[np.isfinite(v)]
            rec[d] = dict(n=int(v.size), median_residual=float(np.median(v)),
                          a=float(1.0 - np.median(v)),
                          frac_below_0p2=float(np.mean(v < 0.2)))
        base = rec[DONORS[0]]["a"]
        for d in DONORS:
            rec[d]["a_rel"] = rec[d]["a"] / base
        routeA[scope] = rec
    out["routeA_self_gene_residual"] = routeA

    # -------------------------------------------------------------------- blob geometry
    blob = {d: E[d].mean(0) for d in DONORS}
    out["blob"] = {
        d: dict(share_of_energy=float(n_t * float(blob[d] @ blob[d]) /
                                      float((E[d] * E[d]).sum())),
                norm=float(np.linalg.norm(blob[d]))) for d in DONORS}
    out["blob"]["cos_hct_hek"] = float(
        blob[DONORS[1]] @ blob[DONORS[2]] /
        (np.linalg.norm(blob[DONORS[1]]) * np.linalg.norm(blob[DONORS[2]])))
    out["blob"]["cos_k562_hct"] = float(
        blob[DONORS[0]] @ blob[DONORS[1]] /
        (np.linalg.norm(blob[DONORS[0]]) * np.linalg.norm(blob[DONORS[1]])))
    out["blob"]["cos_k562_hek"] = float(
        blob[DONORS[0]] @ blob[DONORS[2]] /
        (np.linalg.norm(blob[DONORS[0]]) * np.linalg.norm(blob[DONORS[2]])))

    Ed = {d: E[d] - blob[d][None, :] for d in DONORS}   # blob-deflated

    # ------------------------------------------------------------ route B: Gram / triple
    routeB: dict = {}
    for defl, EE in (("raw", E), ("blob_deflated", Ed)):
        rec: dict = {}
        for scope, m in (("shared", np.ones(n_t, bool)), ("panel", panel_mask),
                         ("extra", ~panel_mask)):
            rows = np.flatnonzero(m)
            g = gram(EE, rows)
            a_rel, u = triple_a(g)
            # signal / noise decomposition in the PRODUCTION frame
            sig = {d: float(u[d] ** 2) for d in DONORS}
            tot = {d: g[f"{d}|{d}"] for d in DONORS}
            noise = {d: float(tot[d] - sig[d]) for d in DONORS}
            snr = {d: sig[d] / noise[d] for d in DONORS}
            rec[scope] = dict(a_rel=a_rel, signal_energy=sig, total_energy=tot,
                              noise_energy=noise, snr=snr,
                              gram={k: v for k, v in g.items()},
                              # minimum-variance weight in the production frame:
                              # w_i ∝ (signal coefficient) / (noise variance) ∝ u_i / N_i
                              w_minvar={d: float(u[d] / noise[d]) for d in DONORS})
            base = rec[scope]["w_minvar"][DONORS[0]]
            rec[scope]["w_minvar_rel"] = {
                d: rec[scope]["w_minvar"][d] / base for d in DONORS}
        # split-half stability on the shared set
        rng = np.random.default_rng(20260903)
        perm = rng.permutation(n_t)
        halves = []
        for h in (perm[: n_t // 2], perm[n_t // 2:]):
            a_rel, _ = triple_a(gram(EE, np.sort(h)))
            halves.append(a_rel)
        rec["split_half"] = halves
        routeB[defl] = rec
    out["routeB_triple"] = routeB

    # -------------------------------------------------- route B: fourth-line instruments
    inst: dict = {}
    for iname, ipath in INSTRUMENTS.items():
        try:
            names = [p.strip() for p in (ipath / "perts.txt").read_text().splitlines()
                     if p.strip()]
        except FileNotFoundError:
            continue
        keep = [t for t in tg if t in set(names)]
        if len(keep) < 40:
            inst[iname] = dict(skipped=f"only {len(keep)} shared targets")
            continue
        from plolo import DONOR_PATHS
        DONOR_PATHS[iname] = ipath
        pi = prepared(iname, tuple(keep), f"inst-{iname}")
        rows = np.array([tg.index(t) for t in keep])
        X = pi["effects"]
        ok = pi["target_available"]
        rec = {}
        for defl, EE in (("raw", E), ("blob_deflated", Ed)):
            Xc = X - X[ok].mean(0)[None, :] if defl == "blob_deflated" else X
            gx = {d: float((EE[d][rows][ok] * Xc[ok]).sum()) for d in DONORS}
            base = gx[DONORS[0]]
            rec[defl] = dict(G=gx, a_rel={d: gx[d] / base for d in DONORS})
        rec["n_targets"] = int(ok.sum())
        rec["n_panel"] = int(sum(1 for t in keep if t in PANEL_SET))
        rec["same_line_as_k562"] = iname in ("k562_essential",)
        inst[iname] = rec
        print(f"[inst/{iname}] n={rec['n_targets']} a_rel(raw)="
              f"{ {d: round(v, 4) for d, v in rec['raw']['a_rel'].items()} }", flush=True)
    out["routeB_instrument"] = inst

    # ------------------------------------------------- route B': naive slope on consensus
    naive: dict = {}
    for defl, EE in (("raw", E), ("blob_deflated", Ed)):
        rec = {}
        for d in DONORS:
            others = [o for o in DONORS if o != d]
            M = 0.5 * (EE[others[0]] + EE[others[1]])
            num = float((EE[d] * M).sum())
            den = float((M * M).sum())
            rec[d] = float(num / den)
        base = rec[DONORS[0]]
        naive[defl] = dict(slope=rec, a_rel={d: rec[d] / base for d in DONORS})
    out["routeBp_naive_slope"] = naive

    # --------------------------------------------------- amplitude in the RAW donor frame
    # the production frame divides each donor by EFFECT_SCALE/raw_median_norm; undo it, so
    # `a_raw` is what the mandate's model calls a_i before any pipeline normalisation.
    conv = {d: P[d]["raw_median_norm"] / EFFECT_SCALE for d in DONORS}
    for defl in ("raw", "blob_deflated"):
        a = routeB[defl]["shared"]["a_rel"]
        out["routeB_triple"][defl]["shared"]["a_rel_donor_frame"] = {
            d: float(a[d] * conv[d] / conv[DONORS[0]]) for d in DONORS}
    out["prod_frame_conversion"] = conv

    (WORK / "atten204.json").write_text(json.dumps(out, indent=2))

    print("\n=== routeA (self-gene residual) ===")
    for scope in routeA:
        print(f"  {scope}: " + "  ".join(
            f"{d} r={routeA[scope][d]['median_residual']:.4f} a={routeA[scope][d]['a']:.4f} "
            f"rel={routeA[scope][d]['a_rel']:.4f}" for d in DONORS))
    print("\n=== routeB triple ===")
    for defl in routeB:
        for scope in ("shared", "panel", "extra"):
            r = routeB[defl][scope]
            print(f"  {defl:14s} {scope:6s} a_rel=" +
                  " ".join(f"{d.split('_')[-1]}:{r['a_rel'][d]:.4f}" for d in DONORS) +
                  "  snr=" + " ".join(f"{r['snr'][d]:.4f}" for d in DONORS) +
                  "  w_mv=" + " ".join(f"{r['w_minvar_rel'][d]:.4f}" for d in DONORS))
        print(f"  {defl:14s} split-half " +
              " | ".join(" ".join(f"{h[d]:.4f}" for d in DONORS)
                         for h in routeB[defl]["split_half"]))
    print("\n=== blob ===", json.dumps(out["blob"], indent=1))
    print("=== naive slope ===", json.dumps(naive, indent=1))
    print("wrote work/atten204.json")


if __name__ == "__main__":
    main()
