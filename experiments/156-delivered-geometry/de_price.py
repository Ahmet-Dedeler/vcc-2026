"""156 — price the DE-side members of the mass fix (AGENTS.md 11.1 rule 12).

a19 lost -0.01352 on the board because its selecting screen priced `pds`, `reach`
and `nmae` and never computed `mse`. The rule that came out of it: an arm must be
priced on ALL SIX members, or must state which are unmeasured AND why they cannot
move. The mass fix changes the emitted profile, so nothing is exempt.

WHAT IS EXACT HERE AND WHAT IS A PROXY
--------------------------------------
Our emitter gives all 400 cells the same library L and pins per-gene totals to the
profile, so a gene's per-cell CPM MEAN is exactly `1e6 * p_g`. The DE members read
that mean. So the arm-to-arm DIFFERENCE in our predicted log2 fold change is
exactly `log2(p_new / p_old)` -- no proxy at all.

The ABSOLUTE lfc needs the reference control's per-cell CPM mean, which the cached
table does not carry; I use the control pseudobulk CPM in its place. That biases
both arms identically and by a per-gene factor, so signs of near-zero genes are the
one place it can matter -- which is exactly what `n_flip` below counts, so the
uncertainty is reported rather than hidden.

`reach` is ordered here by |predicted effect|, NOT by the emitted p-value order
(recomputing that needs the emitter). So the reach column prices the SIGN content
of the ordering, not the shipped order. Flagged, not laundered.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "153-de-bench"))
from debench import Surface, overall_delta  # noqa: E402

from massfix import GATE_CPM, mass_fix, production_effect  # noqa: E402
from screen import TS, render  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DE = os.path.join(REPO, "experiments/095-clean-harness/de_thinned.npz")
EPS = 1e-9


def lfc_from_profile(P_eff, b, ctrl_cpm_pseudobulk):
    """Our predicted log2 fold change, per the DE table's definition."""
    lin = np.clip(np.expm1(b[None, :] + P_eff), 0.0, None)
    p = lin / lin.sum(1, keepdims=True)
    return np.log2((1e6 * p + EPS) / (ctrl_cpm_pseudobulk[None, :] + EPS))


def main() -> None:
    surface = "hct116"
    D, T, b, b_ref, gate, mask, ptar = production_effect(surface)
    tz = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "work", f"truth-{surface}.npz"), allow_pickle=True)
    ctrl_counts = tz["ctrl_counts"].astype(np.float64)
    ctrl_cpm = 1e6 * ctrl_counts / ctrl_counts.sum()

    arms = {}
    arms["delivered (production)"] = render(D, b)
    arms["massfix-minnorm"] = render(mass_fix(D, b, gate, mode="minnorm"), b)
    arms["massfix-const"] = render(mass_fix(D, b, gate, mode="const"), b)

    s = Surface(DE)
    gpos = {g: i for i, g in enumerate(s.genes)}
    cols = np.asarray([gpos.get(g, -1) for g in
                       [str(x) for x in np.load(
                           os.path.join(REPO, "experiments/104-context-heldout-credit-adapter/"
                                        "work/hct116/frozen.npz"), allow_pickle=True)["genes"]]])
    ok = cols >= 0
    prow = {p: i for i, p in enumerate(s.perts)}
    rows = np.asarray([prow.get(p, -1) for p in ptar])
    keep = rows >= 0
    print(f"surface {surface}: {keep.sum()}/{len(ptar)} targets and {ok.sum()}/{len(ok)} genes "
          f"map onto the cached DE table")

    print(f"\n{'arm':<26} {'sign acc':>9} {'head pur':>9} {'reach|eff':>10} "
          f"{'nmae':>9} {'n_flip':>8} {'flip acc':>9}")
    base_lfc = None
    out = {}
    for name, P in arms.items():
        L = lfc_from_profile(P, b, ctrl_cpm)
        M = np.full((len(s.perts), len(s.genes)), np.nan)
        M[np.ix_(rows[keep], cols[ok])] = L[np.ix_(np.flatnonzero(keep), np.flatnonzero(ok))]
        M = np.where(np.isfinite(M), M, 0.0)
        s.pred_lfc = M
        key = [np.abs(M[p]) for p in range(len(s.perts))]
        sa, nm = s.sign_accuracy(), s.nmae()
        rc = s.reach(key=key)
        # head purity, top decile by |predicted effect| over the reference DE pool
        num = den = 0
        for p in range(len(s.perts)):
            i = np.flatnonzero(s.real_sig[p] & s.adj[p])
            if i.size < 10:
                continue
            m = i[np.abs(M[p, i]) >= np.percentile(np.abs(M[p, i]), 90.0)]
            num += int((np.sign(M[p, m]) == np.sign(s.real_lfc[p, m])).sum())
            den += m.size
        hp = num / max(den, 1)

        if base_lfc is None:
            base_lfc, nflip, facc = M, 0, float("nan")
        else:
            pool = s.real_sig & s.adj
            fl = pool & (np.sign(M) != np.sign(base_lfc))
            nflip = int(fl.sum())
            facc = (float((np.sign(M[fl]) == np.sign(s.real_lfc[fl])).mean())
                    if nflip else float("nan"))
        out[name] = dict(sign_acc=sa, head=hp, reach=rc, nmae=nm)
        print(f"{name:<26} {sa:9.4f} {hp:9.4f} {rc:10.5f} {nm:9.4f} "
              f"{nflip:8d} {facc:9.4f}")

    b0 = out["delivered (production)"]
    print(f"\n{'arm':<26} {'d sign':>9} {'d head':>9} {'d reach':>10} {'d nmae':>9} "
          f"{'reach->overall':>15}")
    for name, v in out.items():
        print(f"{name:<26} {v['sign_acc'] - b0['sign_acc']:+9.4f} "
              f"{v['head'] - b0['head']:+9.4f} {v['reach'] - b0['reach']:+10.5f} "
              f"{v['nmae'] - b0['nmae']:+9.4f} "
              f"{overall_delta(d_reach=v['reach'] - b0['reach']):+15.5f}")
    print("\nNOTE reach here is keyed on |predicted effect|, not the emitted p-order;"
          "\n     it prices the ordering's SIGN content only. `n_flip` is the count of"
          "\n     (pert, gene) pairs in the reference DE pool whose predicted SIGN moved,"
          "\n     and `flip acc` is how often the NEW sign is the correct one -- above 0.5"
          "\n     means the flips are net helpful.")


if __name__ == "__main__":
    main()
