"""319 -- the one statement of the weight rule, used by BOTH the fixture test and the real panel.

`_normrel_pre_whiten_effects` (MB:1500-1540) computes, for each panel target `t` and each donor `j`
that carries it,

    nv  = ||s_j[row(t)][dcols_j]||                     (the donor's own official-axis coverage)
    med = median over the panel targets j carries of that same norm
    r   = nv / med
    w_j = (r ** exp) * active_weights[j]               VCC_NORMREL_MODE=plain  (shipped)
    w_j = (max(1, r) ** exp) * active_weights[j]       VCC_NORMREL_MODE=cap    (319)

and falls back to `active_weights[j]` when `med <= 0` or `nv <= 0`.  `model.effect` then renormalises
`w` over the carrying donors, so only the SHARES matter downstream -- but the claim in 319's known-
answer row 3 ("a donor at or below its own panel median keeps weight 1") is about the RAW `w_j`
before that renormalisation, which is what `raw_weights` returns here.

Nothing in this module reads the builder; it is the independent statement the builder is checked
against.  Read-only outside `experiments/319-cap-variant/`.
"""
from __future__ import annotations

import numpy as np


def donor_panel_ratios(
    s: np.ndarray,
    dcols: np.ndarray,
    pert_index: dict[str, int],
    panel: list[str],
) -> tuple[dict[str, float], float]:
    """`r` for every panel target this donor carries, and the donor's own panel median norm.

    Exactly `MB:1519-1527`'s two quantities: the norm is taken over the donor's whole official-axis
    coverage `dcols` (the panel columns, i.e. the target's own knockdown, INCLUDED -- that is
    production's D5, not 308's), on the raw `s` rows.
    """
    vals: dict[str, float] = {}
    for t in panel:
        i = pert_index.get(t)
        if i is None:
            continue
        vals[t] = float(np.linalg.norm(s[i][dcols]))
    med = float(np.median(list(vals.values()))) if vals else 0.0
    if med <= 0.0:
        return {t: float("nan") for t in vals}, med
    return {t: v / med for t, v in vals.items()}, med


def donor_panel_ratios_rows(
    delta: np.ndarray,
    mean: np.ndarray,
    dcols: np.ndarray,
    pert_index: dict[str, int],
    panel: list[str],
) -> tuple[dict[str, float], float]:
    """The same quantities without materialising `s = delta - mean` for the whole bank.

    `(delta - mean)[i]` and `delta[i] - mean` are the same elementwise float32 operation on the
    same operands, so this is bitwise identical to `donor_panel_ratios` on the full `s` -- which is
    checked, on real donor data, by `capshare319.py`'s KA-C2 rather than asserted.  It exists only
    so that a 2.6 GB Orion bank does not need a 2.6 GB copy beside it.
    """
    vals: dict[str, float] = {}
    for t in panel:
        i = pert_index.get(t)
        if i is None:
            continue
        vals[t] = float(np.linalg.norm((delta[i] - mean)[dcols]))
    med = float(np.median(list(vals.values()))) if vals else 0.0
    if med <= 0.0:
        return {t: float("nan") for t in vals}, med
    return {t: v / med for t, v in vals.items()}, med


def raw_weight(ratio: float, exponent: float, mode: str, active: float) -> float:
    """One cell of the builder's `w` dict, before `model.effect` renormalises it."""
    r = ratio
    if mode == "cap":
        r = max(1.0, r)
    return float(r**exponent) * float(active)


def rule_weights(
    ratios: dict[str, dict[str, float]],
    active_weights: dict[str, float],
    panel: list[str],
    exponent: float,
    mode: str,
    carried: dict[str, set[str]],
) -> np.ndarray:
    """The RENORMALISED share matrix (n_targets, n_donors) the builder's `model.effect` sees.

    `ratios[donor][target]`; a donor missing a target, or with a non-positive norm or panel median,
    falls back to `active_weights[donor]` exactly as `MB:1528-1530` does.
    """
    names = list(active_weights)
    W = np.full((len(panel), len(names)), np.nan)
    for i, t in enumerate(panel):
        have = [n for n in names if t in carried[n]]
        if not have:
            continue
        w = np.zeros(len(names))
        for n in have:
            r = ratios[n].get(t, float("nan"))
            j = names.index(n)
            if np.isfinite(r) and r > 0.0:
                w[j] = raw_weight(r, exponent, mode, active_weights[n])
            else:
                w[j] = float(active_weights[n])
        if w.sum() <= 0.0:
            continue
        W[i] = w / w.sum()
    return W
