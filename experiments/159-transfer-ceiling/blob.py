"""159 -- what is the Orion harness's shared blob made of?

`reconcile.py` shows the two thinned Orion harnesses carry a mean effect vector of norm
4.93 / 7.16 where the SAME two lines' donor tables, on the SAME 120 targets and the SAME
17,986 genes, carry 0.60 / 0.52.  That blob is 87% of the harness's true effect energy, it
sets the +0.1790 decoy level, and it is the reason the harness's same/diff ratio is 1.05
where the donor tables' is 54.

Candidate mechanism: a PSEUDOBULK-DEPTH MISMATCH.  b = log1p(TS * P / sum P) is not linear,
and a 150-cell pseudobulk (~2M counts) has a very different zero structure from an
18,400-cell one (~260M counts).  Subtracting the deep control from the shallow target then
leaves a systematic, target-independent vector.  If that is the mechanism, ||mean effect||
must be a function of the DEPTH RATIO between the two pseudobulks and must shrink as the
control is thinned toward the target's depth.

Both sweeps are run.  The scored gene mask is the challenge axis minus the 120 panel target
coordinates.
"""
from __future__ import annotations

import json
import os

import numpy as np

from reliability import load, prof, unit

WORK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "work")


def offdiag_mean(S):
    n = S.shape[0]
    return float(S[~np.eye(n, dtype=bool)].mean())


def main() -> None:
    out = {}
    for name in ("hct116", "hek293t", "h1test"):
        p = os.path.join(WORK, f"shards-{name}.npz")
        if not os.path.exists(p):
            continue
        d = load(name)
        mask, kt, kc = d["mask"], d["kt"], d["kc"]
        cells_t = float(d["tcell"].sum(1).mean())
        cell_per_cshard = float(d["ccell"].mean())
        print(f"\n=== {name}: {d['tsum'].shape[0]} targets x {cells_t:.0f} cells, "
              f"control {d['ccell'].sum():.0f} cells, {len(mask)} scored genes")

        print(f"  (a) control depth sweep, target held at {cells_t:.0f} cells")
        print(f"      {'ctrl cells':>11} {'depth ratio':>12} {'||mean eff||':>13} "
              f"{'med ||row||':>12} {'mean/row':>9} {'within-line cos':>16}")
        rows = []
        for nc in [1, 2, 5, 10, 20]:
            if nc > kc:
                continue
            b = prof(d["csum"][:nc].sum(0))[mask]
            E = prof(d["tsum"].sum(1))[:, mask] - b[None, :]
            me = float(np.linalg.norm(E.mean(0)))
            mr = float(np.median(np.linalg.norm(E, axis=1)))
            a = unit(E)
            w = offdiag_mean(a @ a.T)
            ncells = nc * cell_per_cshard
            rows.append(dict(ctrl_cells=ncells, ratio=ncells / cells_t, mean_norm=me,
                             med_row=mr, within=w))
            print(f"      {ncells:11.0f} {ncells / cells_t:12.2f} {me:13.3f} "
                  f"{mr:12.3f} {me / mr:9.3f} {w:16.4f}")

        print(f"  (b) target depth sweep, control held at "
              f"{d['ccell'].sum():.0f} cells")
        print(f"      {'targ cells':>11} {'depth ratio':>12} {'||mean eff||':>13} "
              f"{'med ||row||':>12} {'mean/row':>9} {'within-line cos':>16}")
        bfull = prof(d["csum"].sum(0))[mask]
        rows_t = []
        for j in (1, kt // 4 or 1, kt // 2, kt):
            E = prof(d["tsum"][:, :j].sum(1))[:, mask] - bfull[None, :]
            me = float(np.linalg.norm(E.mean(0)))
            mr = float(np.median(np.linalg.norm(E, axis=1)))
            w = offdiag_mean(unit(E) @ unit(E).T)
            nt = j * float(d["tcell"].mean())
            rows_t.append(dict(targ_cells=nt, ratio=d["ccell"].sum() / nt,
                               mean_norm=me, med_row=mr, within=w))
            print(f"      {nt:11.0f} {d['ccell'].sum() / nt:12.2f} {me:13.3f} "
                  f"{mr:12.3f} {me / mr:9.3f} {w:16.4f}")

        # (c) the matched-depth read: control thinned to the target's own cell count
        j = int(round(cells_t / cell_per_cshard))
        matched = None
        if 1 <= j <= kc:
            b = prof(d["csum"][:max(j, 1)].sum(0))[mask]
            E = prof(d["tsum"].sum(1))[:, mask] - b[None, :]
            matched = dict(ctrl_cells=max(j, 1) * cell_per_cshard,
                           mean_norm=float(np.linalg.norm(E.mean(0))),
                           within=offdiag_mean(unit(E) @ unit(E).T))
        out[name] = dict(ctrl_sweep=rows, targ_sweep=rows_t, matched=matched,
                         cells_per_target=cells_t)

    with open(os.path.join(WORK, "blob.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'blob.json')}")


if __name__ == "__main__":
    main()
