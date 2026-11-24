"""Arm sweep: per-gene gains on our own predicted log2 fold change.

    prof' = prof * 2^((c_g - 1) * lfc_hat_g),  renormalised    (then a global s on the log1p effect)

`c_g` is a function of the CONTROL's expression only -- buildable with no oracle.  Motivation:
the A-cost of one unit of lfc is not gene-independent.  With x_g = TS * (control composition)
= CPM_g / 20, a small lfc L costs

    dh = log((1 + x 2^L)/(1 + x))  ~  ln2 * L * x/(1+x),

so a unit of predicted fold change costs (x/(1+x))^2 in `A` -- 0.04 at the 5-CPM gate edge and
~1.0 above CPM 1000, a 25x spread -- while `nmae`'s own weight sits disproportionately on the
LOW-expression end (18.7 % / 21.1 % of the |lfc| mass in the bottom CPM decile).

Reported: nmae (measured), pds (exact), A/X/B (exact), and the challenge-mapped raw mse
(DERIVED, by ratio, on the per-context constants of A12-ARTIFACT-AUDIT.md).
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/098-nmae"); sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
import nmaekit as K
from msekit import TS

# ---- challenge constants (A12-ARTIFACT-AUDIT.md) -----------------------------------------
CH_A0, CH_B0, CH_X0 = 0.2829479, 0.2805034, 0.019087      # at s = 1 (a8's amplitude)
CH_KC, CH_SUMD, CH_S12 = 0.0054103, 0.230010, 0.9362727
CH_K = np.array([0.235434, 0.239126, 0.265026])            # per-context cap K_c
CH_BC = np.array([0.2496, 0.2496, 0.2496])                 # ctx budgets at a12 (approx, see note)


def ch_mse(rA, rB, rX, s):
    """raw challenge mse for an arm whose A/B/X are rA/rB/rX times a12's SHAPE, at scale s.
    DERIVED, not measured: the ratios come from a harness, the constants from the board."""
    A = s * s * CH_A0 * rA; B0c = CH_B0 * rB
    # per-context budget shares, back-solved so a12 reproduces its measured mean_c min = 0.240776
    bc = CH_BC * (s * s * rB / (CH_S12 ** 2))
    mn = float(np.mean(np.minimum(CH_K, bc)))
    return 1.0 + (A - 2 * s * CH_X0 * rX + CH_KC - mn) / CH_SUMD


class Arms:
    def __init__(self, which):
        self.H = H = K.Harness(which)
        self.eff0, _ = K.a12_eff(H)
        p = np.expm1(H.bctrl[None, :] + self.eff0); np.clip(p, 0, None, out=p)
        self.prof0 = p / p.sum(1, keepdims=True)
        self.lfc0 = np.zeros_like(p)
        g = H.gate
        self.lfc0[:, g] = np.log2((1e6 * self.prof0[:, g] + 1e-9) / (H.ref_lin[None, :] + 1e-9))
        self.x = np.zeros(H.G); self.x[g] = H.ref_lin / 20.0     # TS * composition
        self.base = self.eval(np.ones(H.G), 1.0)

    def dh(self, c, s=1.0):
        p = self.prof0 * np.exp2((c[None, :] - 1.0) * self.lfc0)
        p /= p.sum(1, keepdims=True)
        d = np.log1p(TS * p) - self.H.bctrl[None, :]
        d = s * d
        for k, col in self.H.cols.items():
            d[k, col] = 0.0
        return d

    def eval(self, c, s=1.0, hold_A=False):
        d = self.dh(c, 1.0)
        A1, X1 = K.AXR(self.H, d); B1 = K.budget(self.H, d)
        if hold_A:
            s = float(np.sqrt(self.base["A"] / A1))
        d = s * d
        A, X, B = A1 * s * s, X1 * s, B1 * s * s
        return dict(A=A, X=X, B=B, s=s, nmae=K.nmae_of(self.H, d), pds=K.pds(self.H, d), dh=d)


def shapes(a):
    return {
        "flat":      lambda x: np.full_like(x, 1.0),
        "invx":      lambda x: np.where(x > 0, 1.0 / np.maximum(x, 1e-9), 0.0),
        "lowband20": lambda x: ((x > 0) & (x < 1.0)).astype(float),      # CPM < 20
        "lowband50": lambda x: ((x > 0) & (x < 2.5)).astype(float),      # CPM < 50
        "lowband100":lambda x: ((x > 0) & (x < 5.0)).astype(float),      # CPM < 100
        "invsqrtx":  lambda x: np.where(x > 0, 1.0 / np.sqrt(np.maximum(x, 1e-9)), 0.0),
    }


if __name__ == "__main__":
    ARMS = []
    for c in (0.7, 0.85, 1.0, 1.15, 1.3, 1.6, 2.0):
        ARMS.append((f"flat c={c}", "flat", c - 1.0))
    for a in (0.25, 0.5, 1.0, 2.0):
        for sh in ("invx", "invsqrtx", "lowband20", "lowband50", "lowband100"):
            ARMS.append((f"{sh} a={a}", sh, a))
    for W_ in ("h1", "orion"):
        M = Arms(W_); H = M.H; g = H.gate
        b = M.base
        print(f"\n===== {W_}  base(a12): nmae {b['nmae']:.4f}  pds {b['pds']:.4f}  "
              f"A {b['A']:.5f}  X {b['X']:.5f}  B {b['B']:.5f} =====")
        print(f"  {'arm':<20}{'s(A-held)':>10}{'nmae':>9}{'Δnmae':>9}{'pds':>8}{'Δpds':>9}"
              f"{'A/A0':>7}{'X/X0':>7}{'B/B0':>7}{'ch mse':>9}{'Δmse':>8}")
        chb = ch_mse(1, 1, 1, CH_S12)
        print(f"  {'a12 base':<20}{1.0:>10.4f}{b['nmae']:>9.4f}{0:>9.4f}{b['pds']:>8.4f}"
              f"{0:>9.4f}{1:>7.3f}{1:>7.3f}{1:>7.3f}{chb:>9.4f}{0:>8.4f}")
        for nm, sh, a in ARMS:
            c = np.ones(H.G); c[g] = 1.0 + a * shapes(a)[sh](M.x)[g]
            r = M.eval(c, hold_A=True)
            rA, rX, rB = r["A"] / b["A"], r["X"] / b["X"], r["B"] / b["B"]
            m = ch_mse(rA, rB, rX, CH_S12)
            print(f"  {nm:<20}{r['s']:>10.4f}{r['nmae']:>9.4f}{r['nmae']-b['nmae']:>+9.4f}"
                  f"{r['pds']:>8.4f}{r['pds']-b['pds']:>+9.4f}{rA:>7.3f}{rX:>7.3f}{rB:>7.3f}"
                  f"{m:>9.4f}{m-chb:>+8.4f}")
