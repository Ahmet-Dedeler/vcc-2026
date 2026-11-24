"""Integer-total rounding rules for experiment 342. Pure NumPy, no emission RNG.

`x` is the intended per-gene count vector `n * L * profile_p` exactly as production
computes it at the one live rounding site (`src/vcc2026/donor_transfer.py`,
`emit_counts`, `mode="quasi"`). Three rules share `base = floor(x)` and the residual
count `k = round(sum(x)) - sum(base)`:

  floor       production, verbatim: `np.floor(x).astype(np.int64)`
  hamilton    `k` increments to the largest fractional parts, ties to the lower index
  random_up   `k` increments uniformly without replacement among fractional parts > 0,
              drawn from a DEDICATED generator seeded `3420909 + target_index`

Every coordinate of every rule differs from `floor(x)` by 0 or 1.

`BlockRounder` is the per-block adapter the seamed emitter calls. Production's
`rc_hybrid_block` calls `emit_counts(mode="quasi")` twice per block: pass 1 on the
full profile, pass 2 on the profile with the scorer-gated genes zeroed. The rule is
decided ONCE, on the first (full) call; any later call must be a coordinate-wise
restriction of the first (each entry either 0 or identical) and receives the
restriction of the decided totals. Floor is restriction-consistent by itself, so its
arm computes production's expression on every call and only checks the invariant.
"""
from __future__ import annotations

import numpy as np

CONTROL_SEED_BASE = 3420909
RULES = ("floor", "hamilton", "random_up")


def _validate(x) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 1:
        raise ValueError("intended totals must be one-dimensional")
    if not np.all(np.isfinite(x)):
        raise ValueError("intended totals carry non-finite entries")
    if np.any(x < 0):
        raise ValueError("intended totals carry negative entries")
    return x


def floor_totals(x) -> np.ndarray:
    """Production's rounding expression, verbatim."""
    return np.floor(x).astype(np.int64)


def fractional_parts(x) -> np.ndarray:
    x = _validate(x)
    return x - np.floor(x)


def residual_count(x) -> int:
    """`round(sum(x)) - sum(floor(x))`, bounded by the number of positive fractions."""
    x = _validate(x)
    base = floor_totals(x)
    k = int(round(float(x.sum()))) - int(base.sum())
    positive = int(np.count_nonzero(x - np.floor(x) > 0))
    if k < 0 or k > positive:
        raise ValueError(f"residual count {k} outside [0, {positive}]")
    return k


def hamilton_choice(x) -> np.ndarray:
    """Indices receiving +1 under largest remainder; stable lower-index tie breaking."""
    x = _validate(x)
    k = residual_count(x)
    frac = x - np.floor(x)
    order = np.argsort(-frac, kind="stable")
    chosen = np.sort(order[:k]).astype(np.int64)
    if np.any(frac[chosen] <= 0):
        raise AssertionError("largest-remainder selected a zero fractional part")
    return chosen


def hamilton_totals(x) -> tuple[np.ndarray, np.ndarray]:
    x = _validate(x)
    chosen = hamilton_choice(x)
    totals = floor_totals(x)
    totals[chosen] += 1
    return totals, chosen


def control_generator(target_index: int) -> np.random.Generator:
    """The rounding control's OWN generator; never the emission stream."""
    return np.random.default_rng(CONTROL_SEED_BASE + int(target_index))


def random_up_choice(x, rng: np.random.Generator) -> np.ndarray:
    x = _validate(x)
    k = residual_count(x)
    candidates = np.flatnonzero(x - np.floor(x) > 0)
    if k == 0:
        return np.zeros(0, dtype=np.int64)
    chosen = rng.choice(candidates, size=k, replace=False)
    return np.sort(np.asarray(chosen, dtype=np.int64))


def random_up_totals(x, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    x = _validate(x)
    chosen = random_up_choice(x, rng)
    totals = floor_totals(x)
    totals[chosen] += 1
    return totals, chosen


def decide(rule: str, x, rng: np.random.Generator | None = None):
    """`(totals, chosen)` for one rule on the full intended vector."""
    if rule == "floor":
        return floor_totals(x), np.zeros(0, dtype=np.int64)
    if rule == "hamilton":
        return hamilton_totals(x)
    if rule == "random_up":
        if rng is None:
            raise ValueError("random_up needs its dedicated generator")
        return random_up_totals(x, rng)
    raise ValueError(f"unknown rounding rule {rule!r}")


def summarize(x, totals, chosen, gate=None) -> dict:
    """Sufficient rounding statistics for one block (JSON-serialisable)."""
    x = _validate(x)
    base = floor_totals(x)
    delta = np.asarray(totals, dtype=np.int64) - base
    frac_positive = x - np.floor(x) > 0
    out = {
        "sum_x": float(x.sum()),
        "round_sum_x": int(round(float(x.sum()))),
        "floor_sum": int(base.sum()),
        "panel_total": int(np.asarray(totals, dtype=np.int64).sum()),
        "residual_k": int(round(float(x.sum())) - int(base.sum())),
        "n_fraction_positive": int(frac_positive.sum()),
        "moved": int(np.count_nonzero(delta)),
        "moved_up": int(np.count_nonzero(delta == 1)),
        "max_abs_delta": int(np.abs(delta).max(initial=0)),
        "min_delta": int(delta.min(initial=0)),
        "chosen_count": int(np.asarray(chosen).size),
        "chosen_all_fractional": bool(np.all(frac_positive[np.asarray(chosen, dtype=np.int64)])),
        "mass_conserved": bool(int(np.asarray(totals, dtype=np.int64).sum())
                               == int(round(float(x.sum())))),
        "zero_one_differences": bool(np.all((delta == 0) | (delta == 1))),
    }
    if gate is not None:
        gate = np.asarray(gate, dtype=bool)
        out["moved_in_scorer_gate"] = int(np.count_nonzero(delta[gate]))
        out["moved_outside_scorer_gate"] = int(np.count_nonzero(delta[~gate]))
    return out


class BlockRounder:
    """Per-block, per-arm totals rule bound into the seamed `emit_counts`.

    Receives ONLY the intended vector `n * L * profile_p`; it never sees, and can never
    consume, the emission generator. The first call decides; later calls must be
    restrictions of the first and get the restricted decision.
    """

    def __init__(self, rule: str, target_index: int):
        if rule not in RULES:
            raise ValueError(f"unknown rounding rule {rule!r}")
        self.rule = rule
        self.target_index = int(target_index)
        self.calls = 0
        self.x = None
        self.totals = None
        self.chosen = None
        self.control_rng = control_generator(target_index) if rule == "random_up" else None

    def __call__(self, x) -> np.ndarray:
        x = _validate(x)
        self.calls += 1
        if self.x is None:
            self.x = x.copy()
            self.totals, self.chosen = decide(self.rule, x, self.control_rng)
            if self.rule == "floor":
                return floor_totals(x)
            return self.totals.copy()
        if x.shape != self.x.shape:
            raise ValueError("a later call changed the gene axis")
        restricted = (x == 0.0) | (x == self.x)
        if not np.all(restricted):
            raise ValueError("a later call is not a restriction of the first")
        if self.rule == "floor":
            return floor_totals(x)
        return np.where(x == 0.0, 0, self.totals).astype(np.int64)

    def summary(self, gate=None) -> dict:
        if self.x is None:
            raise RuntimeError("the rounder was never called")
        out = summarize(self.x, self.totals, self.chosen, gate)
        out.update(rule=self.rule, target_index=self.target_index, calls=self.calls)
        return out
