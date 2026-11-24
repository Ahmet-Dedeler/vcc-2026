"""518 — H1 as a lever surface: 446bj2's bulk recipe emulated on h1-test at the pseudobulk level, levers priced with the EXACT
mse identity (h1_518's legs; the prediction's jackknife held at the harness's, so the claim is fixed and only B moves).

Emulation (no operator, no emission):
  eff_p   = 323's native cap15 effect for target p (the builder's composition on H1's 85 targets; effects.npz).
  k_p     = the per-target scale ps_bulk applied (410's key x energy constant), recovered by a 1-D fit of the render of
            k * eff_p to the REALISED 435 row h_p (log-pseudobulk from the harness sums).  Reported: fit residual.
  render  = the builder's profile_from_effect at MIN_CPM (default 10) without mass fix (323's production MASS_FIX None):
            p = clip(expm1(b + gate*e), 0) / sum;  mu = log1p(TS p);  h = mu - b.
  base    = render(k_p eff_p) (gate 10, q 1)  -> must reproduce the realised mse within the emulation error (reported).
  recipe  = sharpen q 1.75 on the gate (row mode, H1-panel columns untouched, row norm restored), then a global level
            factor f on the effect, chosen at the exact mse optimum (the H1 cliff), as the re-level does on A/B/C.
Arms (each vs the levelled recipe, each re-levelled): render gate 7 / 5 / 2 CPM; across-target strength tilt; top-decile
x1.25; across-target centring; q 2. Also pds (cell-eval2's cosine discrimination: 85 decoys, all 85 H1 target genes dropped).

    uv run python experiments/518-bulk-mse-gap/h1arms518.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
TS = 50_000.0
EFF = ROOT / "experiments/323-current-family-validation/work/effects.npz"


def load():
    z = np.load(WORK / "h1_518_rows.npz")
    jj = json.loads((WORK / "h1_518.json").read_text())
    r = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz")
    genes = [str(g) for g in r["genes"]]
    targets = jj["arrays"]["targets"]
    e = np.load(EFF)
    eg = {g: i for i, g in enumerate(e["genes"].astype(str))}
    assert [str(t) for t in e["targets"]] == targets
    src = np.asarray(e["cap15"], np.float64)
    eff = np.zeros((src.shape[0], len(genes)))          # run323.map_axis: zeros where the H1 axis has no official gene
    for j, g in enumerate(genes):
        if g in eg:
            eff[:, j] = src[:, eg[g]]
    return dict(h=z["h"], dl=z["dl"], mask=z["mask"], b=z["b_ctrl"], wp=z["wp"], Chat=z["Chat"], Cp=z["Cp"], Cc=float(z["Cc"]),
                eff=eff, genes=genes, targets=targets, ident=jj["identity"])


def gate_at(b, m):
    ex = np.expm1(b)
    cpm = 1e6 * ex / ex.sum()
    return (cpm > m) & (b > 0)


def render(E, b, gate):
    G = np.where(gate[None, :], E, 0.0)
    P = np.clip(np.expm1(b[None, :] + G), 0.0, None)
    P = P / P.sum(1, keepdims=True)
    return np.log1p(TS * P) - b[None, :]


def sharpen(E, q, gate, protect):
    g = gate & ~protect
    H = E[:, g]
    S = np.sign(H) * np.abs(H) ** q
    out = E.copy()
    out[:, g] = S * (np.linalg.norm(H, axis=1) / np.maximum(np.linalg.norm(S, axis=1), 1e-300))[:, None]
    return out


class MSE:
    """exact legs with the harness's prediction jackknife (claim fixed), B from the arm's own rows."""
    def __init__(self, D):
        self.D = D
        m, wp = D["mask"], D["wp"]
        self.term = np.minimum(D["Chat"], D["Cp"])
        self.claim = float((wp * self.term).sum())
        self.Dd = ((D["dl"] ** 2) * m).sum(1)
        self.den = (self.Dd - D["Cp"] - D["Cc"]) / (1.0 / wp)
        self.W = m * wp[:, None]

    def raw(self, h):
        D, m, wp = self.D, self.D["mask"], self.D["wp"]
        mu = h + D["b"][None, :]
        wm = (self.W * mu).sum(0) / self.W.sum(0)
        B = float((self.W * (mu - wm[None, :]) ** 2).sum())
        rho = min(1.0, B / self.claim)
        dist = (((h - D["dl"]) ** 2) * m).sum(1)
        N = (dist - rho * self.term - D["Cp"]) * wp
        return float(N.sum() / self.den.sum()), rho


def pds(h, dl, panel_cols):
    keep = np.ones(h.shape[1], bool)
    keep[panel_cols] = False
    P = h[:, keep]
    T = dl[:, keep]
    Pu = P / np.linalg.norm(P, axis=1, keepdims=True)
    Tu = T / np.linalg.norm(T, axis=1, keepdims=True)
    S = Pu @ Tu.T
    own = np.diag(S).copy()
    n = S.shape[0]
    k = (S > own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)     # cosine DISTANCE smaller = similarity larger
    return float(np.mean(1 - k / (n - 1)))


def main():
    D = load()
    M = MSE(D)
    b, eff, h_real = D["b"], D["eff"], D["h"]
    g10 = gate_at(b, 10.0)
    gpos = {g: i for i, g in enumerate(D["genes"])}
    panel_cols = np.asarray([gpos[t] for t in D["targets"] if t in gpos])
    protect = np.zeros(b.size, bool)
    protect[panel_cols] = True
    raw_real, rho_real = M.raw(h_real)
    assert abs(raw_real - D["ident"]["raw"]) < 1e-10, "exact legs vs h1_518"
    pds_real = pds(h_real, D["dl"], panel_cols)
    # ---- recover the per-target scale k_p (1-D fit of render(k eff) to the realised row)
    ks = np.linspace(0.05, 4.0, 400)
    kfit = np.empty(eff.shape[0])
    resid = np.empty(eff.shape[0])
    for i in range(eff.shape[0]):
        errs = [np.sum((render(k * eff[i:i + 1], b, g10)[0] - h_real[i]) ** 2) for k in ks]
        j = int(np.argmin(errs))
        lo, hi = ks[max(j - 1, 0)], ks[min(j + 1, ks.size - 1)]
        for _ in range(40):                          # golden-section refine
            a1, a2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
            e1 = np.sum((render(a1 * eff[i:i + 1], b, g10)[0] - h_real[i]) ** 2)
            e2 = np.sum((render(a2 * eff[i:i + 1], b, g10)[0] - h_real[i]) ** 2)
            lo, hi = (lo, a2) if e1 < e2 else (a1, hi)
        kfit[i] = 0.5 * (lo + hi)
        hr = render(kfit[i] * eff[i:i + 1], b, g10)[0]
        resid[i] = np.linalg.norm(hr - h_real[i]) / np.linalg.norm(h_real[i])
    E0 = kfit[:, None] * eff
    h0 = render(E0, b, g10)
    raw0, rho0 = M.raw(h0)
    print(f"realised 435 row: raw {raw_real:.5f} (rho {rho_real:.4f}) pds {pds_real:.5f}")
    print(f"emulated base (fit k, gate 10, q 1): raw {raw0:.5f} (rho {rho0:.4f}) pds {pds(h0, D['dl'], panel_cols):.5f}; "
          f"row fit residual median {np.median(resid):.3f} max {resid.max():.3f}; k quantiles {np.quantile(kfit, [0, .1, .5, .9, 1]).round(3).tolist()}")

    def levelled(E, gate, q=1.75, fs=np.round(np.arange(0.5, 1.61, 0.01), 2)):
        Eq = sharpen(E, q, gate, protect) if q is not None else E
        best = None
        for f in fs:
            hh = render(f * Eq, b, gate)
            r_, rho_ = M.raw(hh)
            if best is None or r_ < best[0]:
                best = (r_, rho_, f, hh)
        return dict(raw=best[0], rho=best[1], f=float(best[2]), pds=pds(best[3], D["dl"], panel_cols), h=best[3])

    res = {}
    base_q1 = levelled(E0, g10, q=None)
    res["q1 gate10 (levelled)"] = base_q1
    rec = levelled(E0, g10, q=1.75)
    res["RECIPE q1.75 gate10"] = rec
    st = np.linalg.norm(E0[:, ~protect], axis=1)
    rel = st / np.median(st)
    top = st >= np.quantile(st, 0.9)
    arms = {
        "gate7": (E0, gate_at(b, 7.0), 1.75), "gate5": (E0, gate_at(b, 5.0), 1.75), "gate2": (E0, gate_at(b, 2.0), 1.75),
        "gate1": (E0, gate_at(b, 1.0), 1.75), "gate15": (E0, gate_at(b, 15.0), 1.75),
        "q2": (E0, g10, 2.0), "q1.5": (E0, g10, 1.5),
        "tilt^+0.25": (E0 * (rel ** 0.25)[:, None], g10, 1.75), "tilt^+0.5": (E0 * (rel ** 0.5)[:, None], g10, 1.75),
        "tilt^-0.25": (E0 * (rel ** -0.25)[:, None], g10, 1.75), "top10x1.25": (E0 * np.where(top, 1.25, 1.0)[:, None], g10, 1.75),
        "centre_full": (np.where(protect[None, :], E0, E0 - E0.mean(0, keepdims=True)), g10, 1.75),
    }
    for nm, (E, g, q) in arms.items():
        res[nm] = levelled(E, g, q=q)
    # oracle rows for the context: per-target magnitude oracle on the recipe (noise-corrected via h1_518's v)
    print(f"{'arm':26s}{'raw':>9s}{'d raw vs recipe':>17s}{'rho':>7s}{'level f':>9s}{'pds':>9s}{'d pds':>9s}")
    for nm, v in res.items():
        print(f"{nm:26s}{v['raw']:9.5f}{v['raw'] - rec['raw']:+17.5f}{v['rho']:7.3f}{v['f']:9.2f}{v['pds']:9.5f}{v['pds'] - rec['pds']:+9.5f}")
    out = {nm: {k: v for k, v in r_.items() if k != "h"} for nm, r_ in res.items()}
    out["_emulation"] = dict(realised_raw=raw_real, realised_pds=pds_real, emulated_raw=raw0, fit_resid_median=float(np.median(resid)),
                             fit_resid_max=float(resid.max()), k=kfit.tolist())
    (WORK / "h1arms518.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez_compressed(WORK / "h1arms518_recipe.npz", h=rec["h"], E0=E0, kfit=kfit)


if __name__ == "__main__":
    main()
