"""213 — THE TRANSPORT LAW, written down explicitly and then varied.

WHAT THE SHIPPED LAW IS
-----------------------
Every donor table on disk stores

    h_ij = log1p(TS * p^pert_ij) - log1p(TS * p^ctrl_j),      TS = 5e4

verified on this machine (`gate213.py` row 6): `sum_g expm1(b_ctrl + h_row) == 50000.00` on
every sampled row of all three production donors, so the definition is exact and `p^pert` is
exactly recoverable as `(p^ctrl + 1/TS) * exp(h) - 1/TS`.

`profile_from_effect` applies the SAME numbers to a different cell line by ADDITION in that
same space and renders `p^pert_rec = clip(expm1(b_rec + h), 0)/TS`, renormalised to sum 1.

Write `c0 := 1/TS = 2e-5` (= **20 CPM**).  Because `log1p(TS p) = log(p + c0) + log TS` and the
`log TS` cancels in a difference, the shipped law is EXACTLY

    f(p) = log(p + c0);   g = f(p^pert_don) - f(p^ctrl_don);   p^pert_rec = f^-1(f(p^ctrl_rec) + g)

i.e. `log(p + c)` at `c = c0`.  **The incumbent is an interior point of the `log(p+c)` family
by construction, not by fitting**, and `c0` = 20 CPM is exactly the crossover experiments 026
and 169 6.4 both identified as where `log1p` stops being linear and starts being logarithmic.

THE FAMILIES, BOTH THROUGH THE INCUMBENT
----------------------------------------
1. `logc(c)` — the OFFSET family.  `f_c(p) = log(p + c)`.
     c -> 0     : pure log-fold-change / multiplicative
     c = c0     : **SHIPPED**
     c -> inf   : additive in raw abundance
   NOTE the mandate's parametrisation said `c = 1/TS` was the log-fold-change end.  It is not:
   `c = 1/TS` is the SHIPPED point and log-fold-change is `c -> 0`.  Stated because that
   correction is what makes the incumbent interior rather than terminal.

2. `boxcox(q)` — the CURVATURE family.  `f_q(x) = (x^q - 1)/q` on `x = p + c0`.
     q -> 0 : **SHIPPED** (log) · q = 0.5 : sqrt VST · q = 1 : additive · q < 0 : more compressive
   The two families meet only at the incumbent and at the additive limit, so agreement between
   them is evidence about the law rather than about one parametrisation.

3. `rank(kind)` — non-parametric quantile transport on each frame's OWN control-profile CDF
   over genes.  The only member invariant to any monotone re-expression of abundance, and the
   only one whose link differs between donor and recipient.

NUMERICS — every direction is written in the cancellation-free form, and the incumbent is
reached EXACTLY (`D = c - c0` is exactly 0.0, `log1p(0.0)` is exactly 0.0):

    dp    = p^pert - p^ctrl = (p + c0) * expm1(h)                  [exact, from the definition]
    logc  : g = h + log1p( (c - c0) * expm1(-h) / (p + c) )        -> h BIT-EXACT at c = c0
            render p' = (p_rec + c) * exp(g) - c
    boxcox: g = (p + c0)^q * expm1(q h) / q                        -> h as q -> 0 (q == 0 exact)
            render p' = (p_rec + c0) * exp(log1p(q g / (p_rec+c0)^q)/q) - c0

At the incumbent's own parameter both `to_law` and `render` dispatch to the SHIPPED
implementation (`h` unchanged; `expm1(b + g)/TS`), so the identity row is bit-exact rather
than merely close; `gate213.py` additionally measures what the GENERAL path costs there.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

TS = 50_000.0
C0 = 1.0 / TS                      # 2e-5 in relative abundance == 20 CPM


@dataclass
class Law:
    key: str
    kind: str
    param: float
    note: str = ""
    _frames: dict = field(default_factory=dict, repr=False)

    # donor side: shipped log1p effect `h` + the DONOR's basal `p_ctrl` -> law-space effect
    def to_law(self, h: np.ndarray, p_ctrl: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    # recipient side: law-space effect `g` + the RECIPIENT's basal -> perturbed abundance
    def render(self, g: np.ndarray, p_ctrl: np.ndarray, b_ctrl: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def dp_dg(self, p_pert, p_ctrl, g):
        raise NotImplementedError

    @property
    def is_incumbent(self) -> bool:
        return ((self.kind == "logc" and self.param == C0)
                or (self.kind == "boxcox" and self.param == 0.0))

    def asdict(self):
        return dict(key=self.key, kind=self.kind, param=float(self.param), note=self.note,
                    is_incumbent=bool(self.is_incumbent))


class LogC(Law):
    """f(p) = log(p + c).  c == C0 is the shipped law, reached bit-exactly."""

    def to_law(self, h, p_ctrl):
        c = self.param
        d = c - C0                       # exactly 0.0 at the incumbent
        if d == 0.0:
            return np.asarray(h, dtype=np.float64)
        return h + np.log1p(d * np.expm1(-h) / (p_ctrl + c))

    def render(self, g, p_ctrl, b_ctrl):
        c = self.param
        if c == C0:
            return np.expm1(b_ctrl + g) / TS      # production's own expression
        return (p_ctrl + c) * np.exp(g) - c

    def dp_dg(self, p_pert, p_ctrl, g):
        return p_pert + self.param


class BoxCox(Law):
    """f(x) = (x^q - 1)/q on x = p + C0.  q == 0 is the shipped law."""

    def to_law(self, h, p_ctrl):
        q = self.param
        if q == 0.0:
            return np.asarray(h, dtype=np.float64)
        return np.power(p_ctrl + C0, q) * np.expm1(q * h) / q

    def render(self, g, p_ctrl, b_ctrl):
        q = self.param
        if q == 0.0:
            return np.expm1(b_ctrl + g) / TS
        x = p_ctrl + C0
        z = q * g / np.power(x, q)
        bad = z <= -1.0
        out = x * np.exp(np.log1p(np.where(bad, -0.5, z)) / q) - C0
        return np.where(bad, 0.0, out)

    def dp_dg(self, p_pert, p_ctrl, g):
        return np.power(np.maximum(p_pert + C0, 1e-300), 1.0 - self.param)


class Anscombe(Law):
    """f(p) = sqrt(p + 3/(8 TS)) — the count VST, in relative-abundance units."""

    A = 3.0 / (8.0 * TS)

    def to_law(self, h, p_ctrl):
        p_pert = np.maximum((p_ctrl + C0) * np.exp(h) - C0, 0.0)
        return np.sqrt(p_pert + self.A) - np.sqrt(p_ctrl + self.A)

    def render(self, g, p_ctrl, b_ctrl):
        r = np.sqrt(p_ctrl + self.A) + g
        return np.where(r > 0.0, r * r - self.A, 0.0)

    def dp_dg(self, p_pert, p_ctrl, g):
        return 2.0 * np.sqrt(np.maximum(p_pert + self.A, 0.0))


class RawAdd(Law):
    """p'_rec = p_rec + (p^pert_don - p^ctrl_don).  The c -> inf limit, written exactly."""

    def to_law(self, h, p_ctrl):
        return (p_ctrl + C0) * np.expm1(h)

    def render(self, g, p_ctrl, b_ctrl):
        return p_ctrl + g

    def dp_dg(self, p_pert, p_ctrl, g):
        return np.ones_like(np.asarray(p_ctrl))


class MultLFC(Law):
    """p'_rec = p_rec * (p^pert_don / p^ctrl_don).  The c -> 0 limit, written exactly.

    This is the law experiment 026 tested and rejected (-0.0041 macro pds under whitening),
    carried as a named member so 026 gets a known-answer row on the DELIVERED instrument.
    """

    LO = 1e-12

    def to_law(self, h, p_ctrl):
        p_pert = (p_ctrl + C0) * np.exp(h) - C0
        return np.log(np.maximum(p_pert, self.LO)) - np.log(np.maximum(p_ctrl, self.LO))

    def render(self, g, p_ctrl, b_ctrl):
        return p_ctrl * np.exp(g)

    def dp_dg(self, p_pert, p_ctrl, g):
        return p_pert


class Rank(Law):
    """Quantile transport on each frame's own control-profile CDF over genes.

    `param`: 0.0 = uniform scale (plain CDF difference), 1.0 = probit scale.
    """

    def _frame(self, p_ctrl):
        key = (p_ctrl.shape[0], float(p_ctrl.sum()), float(p_ctrl[0]), float(p_ctrl[-1]))
        fr = self._frames.get(key)
        if fr is None:
            xs = np.sort(np.asarray(p_ctrl, dtype=np.float64))
            us = (np.arange(xs.size) + 0.5) / xs.size
            fr = (xs, us)
            self._frames[key] = fr
        return fr

    def _to_u(self, p, p_ctrl):
        xs, us = self._frame(p_ctrl)
        return np.interp(p, xs, us, left=us[0], right=us[-1])

    def _from_u(self, u, p_ctrl):
        xs, us = self._frame(p_ctrl)
        return np.interp(u, us, xs, left=xs[0], right=xs[-1])

    @staticmethod
    def _z(u):
        from scipy.special import ndtri
        return ndtri(np.clip(u, 1e-9, 1.0 - 1e-9))

    @staticmethod
    def _iz(z):
        from scipy.special import ndtr
        return ndtr(z)

    def to_law(self, h, p_ctrl):
        p_pert = np.maximum((p_ctrl + C0) * np.exp(h) - C0, 0.0)
        u1 = self._to_u(p_pert, p_ctrl)
        u0 = self._to_u(p_ctrl, p_ctrl)
        if self.param == 1.0:
            return self._z(u1) - self._z(u0)
        return u1 - u0

    def render(self, g, p_ctrl, b_ctrl):
        u0 = self._to_u(p_ctrl, p_ctrl)
        if self.param == 1.0:
            u1 = self._iz(self._z(u0) + g)
        else:
            u1 = np.clip(u0 + g, 1e-9, 1.0 - 1e-9)
        return self._from_u(u1, p_ctrl)

    def dp_dg(self, p_pert, p_ctrl, g):
        raise NotImplementedError("rank transport has no analytic density; a15 geometry only")


# ------------------------------------------------------------------------------- the registry
def logc(c_cpm: float, note: str = "") -> Law:
    # NOTE `c_cpm / 1e6`, never `c_cpm * 1e-6`: `20.0 * 1e-6` is 1.9999999999999998e-05 and
    # differs from `1/TS` by one ulp, which puts the incumbent 3.4e-21 off its own parameter
    # and costs the identity row 1.3e-15 on every effect.  Division by the exactly
    # representable 1e6 lands on `1/50000.0` bit-for-bit.
    return LogC(key=f"logc{c_cpm:g}", kind="logc", param=c_cpm / 1e6,
                note=note or f"log(p + c), c = {c_cpm:g} CPM")


def boxcox(q: float, note: str = "") -> Law:
    return BoxCox(key=f"bc{q:g}", kind="boxcox", param=float(q),
                  note=note or f"Box-Cox on p + 20 CPM, q = {q:g}")


SHIPPED_KEY = "logc20"

# The offset grid is log-spaced over six decades so both structural limits are reached, and
# is REFINED between 24 and 80 CPM because the first pass put the cosine maximiser there.
C_GRID_CPM = [0.02, 0.2, 1.0, 2.0, 5.0, 10.0, 14.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0,
              45.0, 50.0, 60.0, 80.0, 100.0, 200.0, 500.0, 2000.0, 20000.0]
Q_GRID = [-1.0, -0.5, -0.25, -0.1, -0.05, -0.02, 0.0, 0.02, 0.05, 0.075, 0.1, 0.15,
          0.25, 0.5, 0.75, 1.0]


def all_laws() -> dict[str, Law]:
    laws: dict[str, Law] = {}
    for c in C_GRID_CPM:
        laws[f"logc{c:g}"] = logc(c)
    laws[SHIPPED_KEY].note = "SHIPPED — log1p(TS p); c = 1/TS = 20 CPM. IDENTITY ROW."
    for q in Q_GRID:
        laws[f"bc{q:g}"] = boxcox(q)
    laws["bc0"].note = "SHIPPED, reached from the curvature family. IDENTITY ROW."
    laws["mult"] = MultLFC(key="mult", kind="mult", param=0.0,
                           note="pure log-fold-change / multiplicative — 026's law, exact")
    laws["rawadd"] = RawAdd(key="rawadd", kind="rawadd", param=0.0,
                            note="additive in raw relative abundance (== CPM), exact")
    laws["anscombe"] = Anscombe(key="anscombe", kind="anscombe", param=0.0,
                                note="Anscombe count VST: f(p) = sqrt(p + 3/(8 TS))")
    laws["rank_u"] = Rank(key="rank_u", kind="rank", param=0.0,
                          note="quantile transport, each frame's own control CDF, uniform")
    laws["rank_z"] = Rank(key="rank_z", kind="rank", param=1.0,
                          note="quantile transport, each frame's own control CDF, probit")
    return laws


# ------------------------------------------------------------------------- delivery under a law
def deliver_law(effects: np.ndarray, b_ctrl: np.ndarray, gate: np.ndarray, law: Law,
                mass_fix: str | None = None) -> np.ndarray:
    """Render every row under `law` and read it back the way the scorer does.

    The read-back is ALWAYS `log1p(TS p) - b_ctrl`: that is the space `cell-eval2` measures
    in and it is not ours to choose.  Only the map (donor effect, recipient basal) -> p varies.
    """
    if law.is_incumbent:
        # At its own parameter the law IS the shipped law, so the shipped renderer runs.
        # `gate213.py` measures what the general path costs here (falsifier 1).
        import sys
        from pathlib import Path
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                               / "162-retrodicting-validator"))
        import plolo as _p
        return _p.deliver(effects, b_ctrl, gate, mass_fix)
    b = np.asarray(b_ctrl, dtype=np.float64)
    p_ctrl = np.expm1(b) / TS
    out = np.empty_like(effects)
    for i in range(effects.shape[0]):
        g = np.where(gate, effects[i], 0.0)
        if mass_fix is not None:
            g = mass_consistent_law(g, p_ctrl, b, gate, law, mode=mass_fix)
        p = np.clip(law.render(g, p_ctrl, b), 0.0, None)
        out[i] = np.log1p(TS * (p / float(p.sum()))) - b
    return out


def mass_consistent_law(g, p_ctrl, b_ctrl, on, law: Law, mode="minnorm", iters=64, rtol=1e-12):
    """156's mass fix, generalised to an arbitrary law.

    Production solves `sum_g clip(expm1(b + h), 0) == TS`, i.e. `sum_g clip(p', 0) == 1`, by
    Newton steps along `w = d(TS p')/dh = exp(b + h)`.  Under a general law the same constraint
    is solved along `w = TS * dp'/dg`.  At the incumbent `TS * dp'/dg = TS(p' + c0) = TS p' + 1
    = exp(b + h)` EXACTLY, so this reduces to production's iteration term for term.
    """
    h = np.asarray(g, dtype=np.float64).copy()
    for _ in range(iters):
        p = law.render(h, p_ctrl, b_ctrl)
        live = on & (p > 0.0)
        residual = float(np.clip(p, 0.0, None).sum() * TS - TS)
        if abs(residual) <= rtol * TS:
            break
        w = np.where(live, TS * law.dp_dg(p, p_ctrl, h), 0.0)
        denom = float(w @ w) if mode == "minnorm" else float(w.sum())
        if not np.isfinite(denom) or denom <= 0.0:
            break
        step = -residual / denom
        h = h + (step * w if mode == "minnorm" else np.where(live, step, 0.0))
    return h
