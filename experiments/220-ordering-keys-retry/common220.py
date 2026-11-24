"""220 — shared instrument for the ordering-key retry.

TWO INSTRUMENTS, AND THE WHOLE EXPERIMENT IS ABOUT THE DIFFERENCE BETWEEN THEM.

  OLD (189 / 190's)   `reach(key)` evaluated COUNTERFACTUALLY on a15's frozen h1-test DE
                      table `164/de164.a15.npz`.  Deterministic.  Its uncertainty is a
                      PAIRED BOOTSTRAP OVER THE 85 TARGETS (172 section 4.1), which is what
                      189 quoted as +-0.00424.

  NEW (214's)         the `rc_callset` emitter.  Under it the DELIVERED order IS the key
                      (`rc_emit.py`:983, delivered = prescribed U unsuppressible), so a key
                      can be shipped and read back.  214 measured a THREE-SEED paired sem of
                      0.00009 on this surface against 200's 0.00695 on the `quasi` emitter --
                      a 77x tightening of EMISSION-SEED noise.

Both are provided here so that every arm is measured on both and the two uncertainties are
never confused.  `harness172.check` gates every script (rule 10).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

for p in ("experiments/153-de-bench", "experiments/172-sign-confidence-key",
          "experiments/189-per-gene-reliability", "experiments/190-regulon-ordering"):
    sp = str(ROOT / p)
    if sp not in sys.path:
        sys.path.insert(0, sp)

from debench import Surface  # noqa: E402
import harness172 as H172  # noqa: E402

DE_A15 = ROOT / "experiments/164-callset-artefact/de164.a15.npz"
DE_RC = ROOT / "experiments/202-callset-rider/work/de202.rc_callset.npz"
PRIOR189 = ROOT / "experiments/189-per-gene-reliability/work/prior189.npz"

# 189 / 190 / 172's published rows on the OLD instrument.
EXPECT_OLD = dict(reach_production=0.18283, reach_absH=0.23533,
                  reach_sigfirst_absH=0.21914, reach_oracle=0.64553)
EXPECT_189 = 0.23941        # |H| * exp(+0.50 theta), 189 section 3
EXPECT_189_D = 0.00408
EXPECT_190 = 0.23670        # ct-prop-k25 at lambda = -0.25, 190 section 3
EXPECT_190_D = 0.00137
H_ORACLE_172 = 0.61862      # 172 section 4.2 (the INBOX quotes 0.62119; tie convention)


# ---------------------------------------------------------------------------- keys
def theta_on(genes) -> np.ndarray:
    """189's delivered per-gene prior mapped onto an arbitrary gene axis (0 where absent)."""
    z = np.load(PRIOR189, allow_pickle=True)
    pg = {str(g): i for i, g in enumerate(z["genes"])}
    idx = np.array([pg.get(str(g), -1) for g in genes])
    out = np.zeros(len(genes), dtype=np.float64)
    hit = idx >= 0
    out[hit] = np.asarray(z["theta"])[idx[hit]]
    return out


def key189(H: np.ndarray, theta: np.ndarray, lam: float = 0.5) -> np.ndarray:
    """`|H| * exp(lam * theta)` -- 189's best key.  Magnitude only; sign(H) is untouched."""
    return np.abs(H) * np.exp(np.clip(lam * theta, -50, 50))[None, :]


def signed_key189(H: np.ndarray, theta: np.ndarray, lam: float = 0.5) -> np.ndarray:
    """The same key carried with `H`'s sign -- what the emitter needs (it reads sign(key))."""
    return H * np.exp(np.clip(lam * theta, -50, 50))[None, :]


def agreement_R(N: np.ndarray, H: np.ndarray) -> np.ndarray:
    """190's `R(t, g)`, imported semantics, re-implemented here to avoid a load-order dep."""
    nz = N != 0
    if not nz.any():
        return np.zeros_like(N, dtype=np.float64)
    w = np.abs(N).astype(np.float64)
    med = float(np.median(w[nz]))
    w = w / (med if med > 0 else 1.0)
    return np.where(nz, w * np.sign(N) * np.sign(H), 0.0)


# ---------------------------------------------------------------------------- surfaces
class Instrument:
    """One DE table plus the design-space effect `H`, on the table's own gene axis."""

    def __init__(self, de_path: Path, label: str):
        self.label = label
        self.S = Surface(str(de_path))
        E = np.load(ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz",
                    allow_pickle=True)
        eg = {str(g): i for i, g in enumerate(E["genes"])}
        et = {str(t): i for i, t in enumerate(E["targets"])}
        self.gi = np.array([eg[g] for g in self.S.genes])
        self.pi = np.array([et[p] for p in self.S.perts])
        self.genes = self.S.genes
        self.perts = self.S.perts
        self.H = np.asarray(E["H"])[self.pi][:, self.gi]
        self.b_ctrl = np.asarray(E["b"])[self.gi]
        self.pool = self.S.real_sig & self.S.adj
        self.correct = (np.sign(self.S.pred_lfc) == np.sign(self.S.real_lfc))
        self.theta = theta_on(self.genes)
        self.n_real = self.S.real_sig.sum(1)

    # -- endpoints -------------------------------------------------------------------
    def reach(self, key=None, call=None) -> float:
        return self.S.reach(key=key, call=call)

    def reach_vec(self, key=None, call=None) -> np.ndarray:
        """Per-target `R_p = k*/n_real`, in debench's own order (scored targets only)."""
        S = self.S
        out = []
        for p in range(len(self.perts)):
            pool = np.flatnonzero(S.real_sig[p])
            if pool.size == 0:
                continue
            pool = pool[S._order(p, pool, key, call)]
            pool = pool[S.adj[p, pool]]
            if pool.size == 0:
                out.append(0.0)
                continue
            ok = np.sign(S.pred_lfc[p, pool]) == np.sign(S.real_lfc[p, pool])
            purity = np.cumsum(ok) / np.arange(1, pool.size + 1)
            deep = np.flatnonzero(purity >= 0.9)
            out.append((int(deep[-1] + 1) if deep.size else 0) / int(S.real_sig[p].sum()))
        return np.asarray(out, dtype=np.float64)

    def kstar_vec(self, key=None, call=None) -> np.ndarray:
        """The PURE PREFIX ITSELF, `k*` in genes -- 189's structural statistic."""
        S = self.S
        out = []
        for p in range(len(self.perts)):
            pool = np.flatnonzero(S.real_sig[p])
            if pool.size == 0:
                continue
            pool = pool[S._order(p, pool, key, call)]
            pool = pool[S.adj[p, pool]]
            if pool.size == 0:
                out.append(0)
                continue
            ok = np.sign(S.pred_lfc[p, pool]) == np.sign(S.real_lfc[p, pool])
            purity = np.cumsum(ok) / np.arange(1, pool.size + 1)
            deep = np.flatnonzero(purity >= 0.9)
            out.append(int(deep[-1] + 1) if deep.size else 0)
        return np.asarray(out, dtype=np.int64)

    def scored_mask(self) -> np.ndarray:
        return np.asarray([self.S.real_sig[p].sum() > 0 for p in range(len(self.perts))])

    def front_acc(self, key=None, call=None, depths=(5, 10, 20, 50)) -> dict:
        """Per-target sign accuracy in the first `d` adjudicable pool genes of the order."""
        S = self.S
        acc = {d: [] for d in depths}
        for p in range(len(self.perts)):
            pool = np.flatnonzero(S.real_sig[p])
            if pool.size == 0:
                continue
            pool = pool[S._order(p, pool, key, call)]
            pool = pool[S.adj[p, pool]]
            ok = (np.sign(S.pred_lfc[p, pool]) ==
                  np.sign(S.real_lfc[p, pool])).astype(float)
            for d in depths:
                acc[d].append(float(ok[:d].mean()) if ok.size >= d else np.nan)
        return {d: np.asarray(acc[d]) for d in depths}

    def key_absH(self):
        return np.abs(self.H)

    def key_oracle(self):
        return self.correct.astype(np.float64)

    def key_sigfirst_absH(self):
        return np.where(self.S.pred_sig, 1e6 + np.abs(self.H), np.abs(self.H))


def paired_bootstrap(a, b, n=20000, seed=20260904):
    """172 section 4.1's instrument: resample TARGETS, keep both keys on the same draw."""
    rng = np.random.default_rng(seed)
    d = np.asarray(a) - np.asarray(b)
    idx = rng.integers(0, len(d), size=(n, len(d)))
    boot = d[idx].mean(1)
    return dict(mean=float(d.mean()), se=float(boot.std(ddof=1)),
                lo=float(np.percentile(boot, 2.5)), hi=float(np.percentile(boot, 97.5)),
                p_better=float((boot > 0).mean()))


def check_old(T: Instrument) -> None:
    """rule 10 -- the four published rows of the OLD instrument."""
    got = dict(reach_production=T.reach(), reach_absH=T.reach(T.key_absH()),
               reach_sigfirst_absH=T.reach(T.key_sigfirst_absH()),
               reach_oracle=T.reach(T.key_oracle()))
    bad = []
    for k, v in got.items():
        ok = abs(v - EXPECT_OLD[k]) < 5e-6
        bad += [] if ok else [k]
        print(f"  {'OK ' if ok else 'FAIL'} {k:<24} {v:.5f}   expected {EXPECT_OLD[k]:.5f}")
    if bad:
        raise SystemExit(f"KNOWN-ANSWER ROWS FAILED: {bad} -- aborting (rule 10)")
    return got
