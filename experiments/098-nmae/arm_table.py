"""THE ranked table: every buildable route, on both harnesses, at three gate depths.

Gate depth matters because the challenge's reference calls ~2,000 genes per row and our
harnesses call 89 (orion) / 1,120 (h1) at BH 0.05.  `gateladder.py` shows the optimal global
scale moves from ~1.9 at BH 0.05 to ~0.5 at a challenge-deep gate, so an arm ranked at BH 0.05
alone is ranked on a cohort the board does not use.  Every arm is therefore read at
BH 0.05 / 0.35 / 0.60, whose mean |S_p| brackets the challenge's ~2,000.

Costs reported alongside: A/A0 (the quadratic term of `expr_mse`, exact) and `pds` (exact).
`nmae` is deterministic under re-seeding -- validate.py measures the same-policy twin gap at
1e-5 on both harnesses, four orders below every entry here -- so the error bars are the two
harnesses and the three gate depths, not the seed.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/098-nmae"); sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
import nmaekit as K
from linear import donor_keys
from msekit import TS
from scipy.stats import rankdata
from cell_eval2.metrics.delta import _across_pert_budget

DEPTHS = (0.05, 0.35, 0.60)


def desc_midrank(x):
    n = x.size; o = np.argsort(-x, kind="stable"); od = x[o]
    e = np.flatnonzero(np.concatenate(([True], od[1:] != od[:-1], [True])))
    st, en = e[:-1], e[1:]
    out = np.empty(n); out[o] = np.repeat((st + en - 1) / 2.0, en - st)
    return out


class Bench:
    def __init__(self, which):
        self.w = which
        self.H = H = K.Harness(which)
        self.eff0, _ = K.a12_eff(H)
        p = np.expm1(H.bctrl[None, :] + self.eff0); np.clip(p, 0, None, out=p)
        self.prof0 = p / p.sum(1, keepdims=True)
        g = H.gate
        self.lfc0 = np.zeros_like(p)
        self.lfc0[:, g] = np.log2((1e6 * self.prof0[:, g] + 1e-9) / (H.ref_lin[None, :] + 1e-9))
        # the control-normalisation offset: our zero-effect lfc is not zero
        z0 = np.expm1(H.bctrl); z0 = 1e6 * z0 / z0.sum()
        self.off = np.zeros(H.G); self.off[g] = np.log2((z0[g] + 1e-9) / (H.ref_lin + 1e-9))
        self.ag, self.un = donor_keys(H, which)
        z = np.load(f"{ROOT}/experiments/098-nmae/work/refde_{which}.npz", allow_pickle=True)
        self.q = z["q"]; self.rlfc = z["lfc"].astype(np.float64)
        self.base = None
        self.base = self.metrics(np.ones_like(self.eff0))

    def gate_rows(self, alpha):
        H = self.H; out = []
        for k, (i, j) in enumerate(zip(H.rows, H.rrows)):
            s = (self.q[j] < alpha) & np.isfinite(self.rlfc[j]) & (self.rlfc[j] != 0)
            gs = np.flatnonzero(s)
            if gs.size >= 10: out.append((k, j, gs))
        return out

    def dh_of(self, gain, shift=0.0, per_row=None):
        """gain: [P,G] or [G] multiplier on our own lfc; shift: constant lfc added on the gate."""
        H = self.H
        d = (np.asarray(gain) - 1.0) * self.lfc0
        if shift: d = d + shift * H.gate[None, :]
        p = self.prof0 * np.exp2(d)
        p /= p.sum(1, keepdims=True)
        dh = np.log1p(TS * p) - H.bctrl[None, :]
        if per_row is not None: dh = dh * np.asarray(per_row)[:, None]
        for k, col in H.cols.items(): dh[k, col] = 0.0
        return dh

    def metrics(self, gain, shift=0.0, per_row=None, hold_A=False):
        H = self.H
        dh = self.dh_of(gain, shift, per_row)
        if hold_A:
            A1, _ = K.AXR(H, dh)
            dh = dh * float(np.sqrt(self.base["A"] / A1))
        A, X = K.AXR(H, dh); B = float(_across_pert_budget(dh, H.cols, H.w))
        L = H.lfc_hat(dh[H.rows])
        out = dict(A=A, X=X, B=B, pds=K.pds(H, dh))
        for al in DEPTHS:
            v = []
            for k, j, gs in self.gate_rows(al):
                r = self.rlfc[j][gs]
                v.append(np.abs(L[k][gs] - r).sum() / np.abs(r).sum())
            out[f"n{al}"] = float(np.mean(v))
        return out


def arms(Bn):
    H = Bn.H; P, G = H.P, H.G; g = H.gate
    A = {}
    for c in (0.7, 0.85, 1.15, 1.5, 2.0):
        A[f"lfc-scale c={c}"] = dict(gain=np.full(G, c))
    # sharpen / flatten on our own |lfc| within each row
    for b in (-0.5, -0.25, 0.25, 0.5):
        gn = np.ones((P, G))
        a = np.abs(Bn.lfc0)
        med = np.where(a[:, g].mean(1) > 0, a[:, g].mean(1), 1.0)[:, None]
        gn[:, g] = np.power(np.maximum(a[:, g], 1e-6) / med, b)
        gn[:, g] /= gn[:, g].mean(1, keepdims=True)          # hold mean gain at 1
        A[f"sharpen |lfc|^{b:+.2f}"] = dict(gain=gn)
    # a13's key: the builder's exponential on the agreement midrank, at CONSTANT mean gain
    for gam in (0.25, 1.0):
        gn = np.ones((P, G))
        for i in range(P):
            gn[i] = 1.0 if np.unique(Bn.ag[i]).size < 2 else \
                np.exp(-gam * desc_midrank(Bn.ag[i]) / 2000.0)
            gn[i] /= gn[i].mean()
        A[f"agreement shrink g={gam}"] = dict(gain=gn)
        A[f"agreement AMPLIFY g={gam}"] = dict(gain=(1.0 / gn) / (1.0 / gn).mean(1, keepdims=True))
    # expression tilt
    x = np.zeros(G); x[g] = H.ref_lin / 20.0
    for b in (-0.5, 0.5):
        gn = np.ones(G)
        gn[g] = np.power(np.maximum(x[g], 1e-6), b); gn[g] /= gn[g].mean()
        A[f"cpm tilt ^{b:+.2f}"] = dict(gain=gn)
    A["de-offset control"] = dict(gain=np.ones(G), shift=None)   # handled specially
    for dl in (-0.02, +0.02):
        A[f"lfc shift {dl:+.3f}"] = dict(gain=np.ones(G), shift=dl)
    # per-perturbation: donor cosine rank
    return A


if __name__ == "__main__":
    for W_ in ("h1", "orion"):
        Bn = Bench(W_); b = Bn.base; H = Bn.H
        print(f"\n===== {W_}  a12 base: nmae {b['n0.05']:.4f}/{b['n0.35']:.4f}/{b['n0.6']:.4f}"
              f"  pds {b['pds']:.4f}  A {b['A']:.5f} =====")
        print(f"  {'arm':<26}{'A/A0':>7}{'pds':>8}{'Δpds':>8}"
              + "".join(f"{'nmae@'+str(d):>12}" for d in DEPTHS)
              + "".join(f"{'Δ@'+str(d):>10}" for d in DEPTHS))
        for nm, kw in arms(Bn).items():
            if nm == "de-offset control":
                r = Bn.metrics(np.ones(H.G), shift=0.0, per_row=None)
                dh = Bn.dh_of(np.ones(H.G))
                p = Bn.prof0 * np.exp2(-Bn.off[None, :]); p /= p.sum(1, keepdims=True)
                dh = np.log1p(TS * p) - H.bctrl[None, :]
                for k, col in H.cols.items(): dh[k, col] = 0.0
                A_, X_ = K.AXR(H, dh); L = H.lfc_hat(dh[H.rows])
                r = dict(A=A_, X=X_, pds=K.pds(H, dh))
                for al in DEPTHS:
                    v = [np.abs(L[k][gs] - Bn.rlfc[j][gs]).sum() / np.abs(Bn.rlfc[j][gs]).sum()
                         for k, j, gs in Bn.gate_rows(al)]
                    r[f"n{al}"] = float(np.mean(v))
            else:
                r = Bn.metrics(kw["gain"], kw.get("shift") or 0.0)
            print(f"  {nm:<26}{r['A']/b['A']:>7.3f}{r['pds']:>8.4f}{r['pds']-b['pds']:>+8.4f}"
                  + "".join(f"{r['n'+str(d)]:>12.4f}" for d in DEPTHS)
                  + "".join(f"{r['n'+str(d)]-b['n'+str(d)]:>+10.4f}" for d in DEPTHS))
