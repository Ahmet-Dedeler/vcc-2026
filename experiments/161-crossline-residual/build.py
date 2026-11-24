"""161 SS2 — build a mass-consistent multi-line surface.

A donor table stores `b = log1p(TS*p)` and `delta_p = b_p - b_ctrl`, so both `b` and
`b + delta_p` satisfy `sum(expm1(.)) == TS` **exactly**.  Experiment 156 showed that this
manifold is the whole story of the renderer: the truth lives on it and our effect does not.

Restricting to a shared gene subset breaks it, so the restriction is done in PROFILE space and
the frame is rebuilt:

    p      = expm1(b) / TS                        (sums to 1 on the full axis)
    p_t    = expm1(b + E_t) / TS                  (sums to 1, per target)
    q      = p[idx] / p[idx].sum() ;  b_sub = log1p(TS*q)
    q_t    = p_t[idx] / p_t[idx].sum() ;  E_sub = log1p(TS*q_t) - b_sub

after which `sum(expm1(b_sub + E_sub)) == TS` again, exactly.  159 §2 measured that this exact
renormalisation moves a donor table's within-line cosine by 0.0003, i.e. nothing.

The CPM gate is carried from the FULL axis, not recomputed in the sub-frame — a 5k-gene frame
concentrates the profile ~4x and a recomputed gate would keep everything, which is not what
production does.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from common import TS, WORK, load_table, panel_targets

GATE_CPM = 10.0                     # a15's shipped config.min_cpm

SURFACES = {
    # one table per distinct cell line; the whole target pool is OFF the challenge panel
    "core5": ["k562_gwps", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t"],
    # + the second K562 library, so a fold can see two views of one line
    "core6": ["k562_gwps", "k562_essential", "rpe1", "nadig_jurkat",
              "orion_hct116", "orion_hek293t"],
    # H1 as a held line; 400 cells/target, the challenge's own depth
    "wideh1": ["k562_gwps", "orion_hct116", "orion_hek293t", "h1_2025"],
    # rpe1 carries a target-independent blob (159 §4); this is the clean-truth control
    "core4": ["k562_gwps", "nadig_jurkat", "orion_hct116", "orion_hek293t"],
}


def subframe(b_full: np.ndarray, E_full: np.ndarray, idx: np.ndarray):
    """Restrict (b, E) to `idx` and rebuild the log1p(TS*p) frame exactly."""
    lin = np.expm1(b_full)                                   # == TS * p
    sub = lin[idx]
    b_sub = np.log1p(TS * sub / sub.sum())
    linp = np.clip(np.expm1(b_full[None, :] + E_full), 0.0, None)
    subp = linp[:, idx]
    E_sub = np.log1p(TS * subp / subp.sum(1, keepdims=True)) - b_sub[None, :]
    return b_sub, E_sub


def build(name: str, max_targets: int | None = None) -> None:
    members = SURFACES[name]
    panel = panel_targets()
    tabs = {m: load_table(m) for m in members}
    assert all(t is not None for t in tabs.values()), f"{name}: a member table is missing"

    T = sorted(set.intersection(*[set(t["perts"]) for t in tabs.values()]))
    G = set.intersection(*[set(t["genes"]) for t in tabs.values()]) - set(T)
    G = sorted(G)
    if max_targets and len(T) > max_targets:
        rng = np.random.default_rng(20260902)
        T = [T[i] for i in sorted(rng.choice(len(T), max_targets, replace=False))]
    n_panel = len(set(T) & panel)
    print(f"[{name}] {len(members)} tables, {len(T)} targets ({n_panel} on the challenge "
          f"panel), {len(G)} genes")

    out = {"targets": np.asarray(T), "genes": np.asarray(G),
           "members": np.asarray(members),
           "lines": np.asarray([tabs[m]["line"] for m in members])}
    for m in members:
        t = tabs[m]
        gi = {g: i for i, g in enumerate(t["genes"])}
        pi = {p: i for i, p in enumerate(t["perts"])}
        gidx = np.asarray([gi[g] for g in G])
        ridx = np.asarray([pi[p] for p in T])

        E_full = np.empty((len(T), len(t["genes"])), dtype=np.float64)
        for k in np.argsort(ridx):
            E_full[k] = np.nan_to_num(np.asarray(t["delta"][ridx[k]], dtype=np.float64))
        b_full = t["b"].astype(np.float64)

        # gate on the FULL axis, exactly as the builder does
        cpm_full = 1e6 * np.expm1(b_full) / TS
        gate_full = (cpm_full > GATE_CPM) & (b_full > 0)

        b_sub, E_sub = subframe(b_full, E_full, gidx)
        mass = np.expm1(b_sub[None, :] + E_sub).sum(1) / TS
        assert np.abs(mass - 1.0).max() < 1e-9, f"{m}: sub-frame mass broken {np.abs(mass-1).max()}"

        out[f"E_{m}"] = E_sub.astype(np.float32)
        out[f"b_{m}"] = b_sub
        out[f"gate_{m}"] = gate_full[gidx]
        out[f"btg_{m}"] = np.asarray(          # each target gene's own basal in this line
            [b_full[gi[p]] if p in gi else np.nan for p in T])
        print(f"   {m:<16} line={t['line']:<8} gated {out[f'gate_{m}'].sum():5d}/{len(G)}"
              f"  |E| median row norm {np.median(np.linalg.norm(E_sub, axis=1)):.3f}"
              f"  target-gene basal known for {np.isfinite(out[f'btg_{m}']).sum()}/{len(T)}")
        del E_full

    p = os.path.join(WORK, f"surface-{name}.npz")
    np.savez_compressed(p, **out)
    print(f"   wrote {p} ({os.path.getsize(p) / 1e6:.0f} MB)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("surfaces", nargs="*", default=list(SURFACES))
    ap.add_argument("--max-targets", type=int, default=None)
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    for s in (a.surfaces or list(SURFACES)):
        build(s, a.max_targets)


if __name__ == "__main__":
    main()
