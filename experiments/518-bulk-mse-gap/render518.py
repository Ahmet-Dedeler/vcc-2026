"""518 — how the bulk effect is RENDERED on low/mid-expressed genes: the 10-CPM gate and the additive-log1p render, tested
where truth and prediction share one control frame (the board's situation).

The builder renders p = expm1(b + gate*h)/TS (normalised): h is an ADDITIVE change in log1p(TS p). For a gene near the
gate (TS p ~ 0.5) that turns h = -0.3 into a 78 % drop, where a fold-change reading of the same h (p' = p e^h) is a 26 %
drop; above ~100 CPM the two agree. PLOLO renders in the challenge frames while its truth lives in the held line's frame,
so it cannot see this. Two same-frame surfaces:
  (1) Orion HCT116 / HEK293T folds rendered in the HELD line's own control frame (fold cache b_ctrl), truth = the held
      line's raw delta rows (uncentred, its own log1p space), prediction = 446bj2's keyed stack (score518 dump, the fold's
      LODO key clip (0.5, 3)), sharpened q 1.75 on the frame's 10-CPM gate. Designed rows (no ballast). cos_agg + pds.
      (the K562 and kolf2 banks were pre-filtered to expressed genes: nothing below 10 CPM to score, so they are skipped.)
  (2) H1: h1arms518's emulated recipe, exact mse at each variant's own optimal level (the H1 cliff), and pds.
Variants: gate g in {10 (production), 7, 5, 2, 0 (= b > 0)} with the additive render; the fold-change render
p = p_ctrl exp(h) at gates 10 / 2 / 0; the HYBRID: additive at >= 10 CPM (production untouched) + fold-change render below.
Sharpening always on the 10-CPM gate (the production shape); genes below it carry their unsharpened h.

    uv run python experiments/518-bulk-mse-gap/render518.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit"))
import json  # noqa: E402

import numpy as np  # noqa: E402

TS = 50_000.0
WORK = HERE / "work"


def cpm_of(b):
    ex = np.expm1(b)
    return 1e6 * ex / ex.sum()


def render(E, b, mode, gate_add, gate_fc=None):
    """rows of delivered effect h = log1p(TS p) - b.  mode 'add': p ~ expm1(b + h) on gate_add; 'fc': p ~ expm1(b) e^h on
    gate_add; 'hyb': 'add' on gate_add plus 'fc' on gate_fc (disjoint)."""
    base = np.clip(np.expm1(b), 0.0, None)[None, :]
    if mode == "add":
        P = np.clip(np.expm1(b[None, :] + np.where(gate_add[None, :], E, 0.0)), 0.0, None)
    elif mode == "fc":
        P = base * np.exp(np.where(gate_add[None, :], E, 0.0))
    elif mode == "hyb":
        P = np.clip(np.expm1(b[None, :] + np.where(gate_add[None, :], E, 0.0)), 0.0, None)
        P = np.where(gate_fc[None, :], base * np.exp(np.where(gate_fc[None, :], E, 0.0)), P)
    P = P / P.sum(1, keepdims=True)
    return np.log1p(TS * P) - b[None, :]


def sharpen_on(E, q, gate):
    H = E[:, gate]
    S = np.sign(H) * np.abs(H) ** q
    out = E.copy()
    out[:, gate] = S * (np.linalg.norm(H, axis=1) / np.maximum(np.linalg.norm(S, axis=1), 1e-300))[:, None]
    return out


def cos_pds(P, T):
    cos = float((P * T).sum() / np.sqrt((P * P).sum() * (T * T).sum()))
    Pu = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    Tu = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = Pu @ Tu.T
    own = np.diag(S).copy()
    n = S.shape[0]
    pds = float(np.mean(((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)))
    return cos, pds


def variants(b):
    c = cpm_of(b)
    g = lambda m: (c > m) & (b > 0)
    g10 = g(10.0)
    V = {"add g10 (production)": ("add", g10, None)}
    for m in (7.0, 5.0, 2.0, 0.0):
        V[f"add g{m:g}"] = ("add", g(m), None)
    for m in (10.0, 2.0, 0.0):
        V[f"fc g{m:g}"] = ("fc", g(m), None)
    for m in (5.0, 2.0, 0.0):
        V[f"hyb add>=10 + fc>{m:g}"] = ("hyb", g10, g(m) & ~g10)
    return V, g10


def folds():
    from size410 import donor_rows
    from lib420 import key_scale
    out = {}
    for f in ("orion_hct116", "orion_hek293t"):
        dz = np.load(WORK / f"dump518-{f}.npz")
        fz = np.load(ROOT / "experiments/420-bulk-stack/work/plolo" / f"fold-{f}.npz", allow_pickle=True)
        targets = [str(t) for t in dz["targets"]]
        b = np.asarray(fz["b_ctrl"], np.float64)
        sc = dz["scored"]
        from lib420 import weights
        w = weights(dz["kv1"], dz["carried"], 0.5, 3.0)[0]
        E = key_scale(dz["fST"], w)[0]
        R, have = donor_rows(f, targets)
        T = R[:, sc]
        V, g10 = variants(b)
        Eq = sharpen_on(E, 1.75, g10 & np.isin(np.arange(b.size), dz["panel_positions"], invert=True))
        rec = {}
        cpm = cpm_of(b)
        low = (cpm[sc] <= 10) & (b[sc] > 0)
        tl = float((T[:, low] ** 2).sum() / (T ** 2).sum())
        for nm, (mode, ga, gf) in V.items():
            h = render(Eq, b, mode, ga, gf)[:, sc]
            cos, pds = cos_pds(h, T)
            rec[nm] = dict(cos=cos, pds=pds, low_share_of_our_energy=float((h[:, low] ** 2).sum() / (h ** 2).sum()),
                           low_cos=float((h[:, low] * T[:, low]).sum() / np.sqrt(max((h[:, low] ** 2).sum(), 1e-300) * (T[:, low] ** 2).sum())))
        b0 = rec["add g10 (production)"]
        print(f"== {f} (held frame; truth energy share on genes <= 10 CPM, noise-inflated: {100 * tl:.1f} %)")
        for nm, v in rec.items():
            print(f"   {nm:26s} cos {v['cos']:.5f} (r {100 * (v['cos'] / b0['cos'] - 1):+6.2f} %)  pds {v['pds']:.5f} ({v['pds'] - b0['pds']:+.5f})"
                  f"  our energy <10cpm {100 * v['low_share_of_our_energy']:5.1f} %  cos there {v['low_cos']:+.4f}")
        out[f] = dict(truth_low_share=tl, variants=rec)
    return out


def h1():
    import h1arms518 as H
    D = H.load()
    M = H.MSE(D)
    b = D["b"]
    z = np.load(WORK / "h1arms518_recipe.npz")
    E0 = z["E0"]
    gpos = {g: i for i, g in enumerate(D["genes"])}
    panel = np.zeros(b.size, bool)
    panel[[gpos[t] for t in D["targets"] if t in gpos]] = True
    V, g10 = variants(b)
    Eq = sharpen_on(E0, 1.75, g10 & ~panel)
    pcols = np.flatnonzero(panel)
    rec = {}
    for nm, (mode, ga, gf) in V.items():
        best = None
        for f in np.round(np.arange(0.6, 1.61, 0.01), 2):
            h = render(f * Eq, b, mode, ga, gf)
            r_, rho = M.raw(h)
            if best is None or r_ < best[0]:
                best = (r_, rho, float(f), h)
        rec[nm] = dict(raw=best[0], rho=best[1], level=best[2], pds=H.pds(best[3], D["dl"], pcols))
    b0 = rec["add g10 (production)"]
    print("== H1 (exact mse at each variant's optimal level)")
    for nm, v in rec.items():
        print(f"   {nm:26s} raw {v['raw']:.5f} ({v['raw'] - b0['raw']:+.5f})  pds {v['pds']:.5f} ({v['pds'] - b0['pds']:+.5f})  level {v['level']:.2f} rho {v['rho']:.3f}")
    return rec


def main():
    sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))
    out = dict(folds=folds(), h1=h1())
    (WORK / "render518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
