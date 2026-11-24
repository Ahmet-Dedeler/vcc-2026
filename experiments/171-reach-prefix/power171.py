"""171 — does the knee arm's reach gain survive challenge-like power?  Twin-averaged.

h1-test's native `n_real` is 1-5,484 (median 115) and `diag171.py` measures that 93.6% of its
reach comes from the perturbations with `n_real <= 13`.  The challenge sits near `n_real ~ 240`
with no such stratum.  157 section 3b's proportional instrument -- pool_p = the reference's top
round(f * n_real_p) by its own p-value -- is the only power match `reach` can be read through
(uniform `n_real_target` deletes the stratum that carries the score; 157 section 3a).

f = 1 is the known-answer row: it must reproduce the native BH set, i.e. a15 = 0.18283.

This is the check that decides the arm.  A gain that lives only at f = 1 lives only in a
stratum the challenge does not have.
"""
from __future__ import annotations

import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/171-reach-prefix")
from diag171 import Walk, kstar, spec_order  # noqa: E402

W = f"{ROOT}/experiments/171-reach-prefix/work"
D164 = f"{ROOT}/experiments/164-callset-artefact"
FACTORS = [0.5, 1.0, 1.5, 2.0, 3.5]
REACH_SPAN = 0.96414 - 0.07902

PAIRS = {
    "a15 (incumbent)":     (f"{W}/de171.a15.npz",          f"{D164}/de164.a15_b.npz"),
    "kt9000":              (f"{W}/de171.kt9000.npz",       f"{W}/de171.kt9000_b.npz"),
    "kt10000 = KNEE arm":  (f"{W}/de171.kt10000.npz",      f"{W}/de171.kt10000_b.npz"),
    "kt20000 (all gated)": (f"{W}/de171.kt20000.npz",      f"{W}/de171.kt20000_b.npz"),
    "armC a8 flat":        (f"{W}/de171.armC_a8_flat.npz", f"{W}/de171.armC_a8_flat_b.npz"),
}


def reach_at(path, f):
    Wk = Walk(path)
    z = np.load(path, allow_pickle=True)
    if f == 1.0:
        sig = Wk.real_sig
    else:
        rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
        sig = np.zeros(rp.shape, bool)
        for i in range(rp.shape[0]):
            k = int(round(f * int(Wk.real_sig[i].sum())))
            k = max(0, min(k, int(np.isfinite(z["real_p"][i]).sum())))
            if k:
                sig[i, np.argsort(rp[i], kind="stable")[:k]] = True
        sig &= np.isfinite(z["real_p"])
        for p, j in enumerate(Wk.self_idx):
            if j >= 0:
                sig[p, j] = False
    out, nr = [], []
    for p in range(Wk.perts.size):
        pool = np.flatnonzero(sig[p])
        if pool.size == 0:
            continue
        nr.append(pool.size)
        o = spec_order(pool, Wk.genes, Wk.pred_sig[p], Wk.pred_padj[p], Wk.pred_p[p],
                       Wk.pred_lfc[p])
        w = pool[o]
        w = w[Wk.adj[p, w]]
        if w.size == 0:
            out.append(0.0)
            continue
        ok = np.sign(Wk.pred_lfc[p, w]) == np.sign(Wk.real_lfc[p, w])
        out.append(kstar(ok)[0] / pool.size)
    return float(np.mean(out)), int(np.median(nr))


def main():
    print("TWIN-AVERAGED raw reach vs proportional power factor f")
    print(f"  {'arm':<22}" + "".join(f"{'f='+str(f):>10}" for f in FACTORS))
    med = {}
    vals = {}
    for nm, (a, b) in PAIRS.items():
        row, sp = [], []
        for f in FACTORS:
            ra, na = reach_at(a, f)
            rb, _ = reach_at(b, f)
            row.append((ra + rb) / 2)
            sp.append(abs(ra - rb))
            med[f] = na
        vals[nm] = (row, sp)
        print(f"  {nm:<22}" + "".join(f"{v:>10.5f}" for v in row))
    print(f"  {'median n_real':<22}" + "".join(f"{med[f]:>10d}" for f in FACTORS))
    print(f"  {'(the challenge is n_real ~ 240, i.e. between f = 2 and f = 3.5)'}")
    print()
    print("  SEED SPREAD |A - B| behind each cell")
    for nm, (row, sp) in vals.items():
        print(f"  {nm:<22}" + "".join(f"{v:>10.5f}" for v in sp))
    print()

    base = vals["a15 (incumbent)"][0]
    print("DELTA vs a15, twin-averaged, and whether it beats its own seed spread")
    print(f"  {'arm':<22}" + "".join(f"{'f='+str(f):>10}" for f in FACTORS))
    for nm, (row, sp) in vals.items():
        if nm == "a15 (incumbent)":
            continue
        cells = []
        for i, f in enumerate(FACTORS):
            d = row[i] - base[i]
            noise = max(sp[i], vals["a15 (incumbent)"][1][i])
            cells.append(f"{d:>+9.5f}{'*' if abs(d) > noise else ' '}")
        print(f"  {nm:<22}" + "".join(cells))
    print("  (* = the twin-averaged delta exceeds the larger of the two seed spreads)")
    print()

    shrink = 0.0079 / (vals["armC a8 flat"][0][1] - base[1])
    print(f"BOARD-SCALED at each f (shrink {shrink:.4f}, calibrated on armC at f = 1)")
    print(f"  {'arm':<22}" + "".join(f"{'f='+str(f):>10}" for f in FACTORS))
    for nm, (row, sp) in vals.items():
        if nm == "a15 (incumbent)":
            continue
        print(f"  {nm:<22}"
              + "".join(f"{(row[i]-base[i])*shrink/REACH_SPAN/6:>+10.6f}"
                        for i in range(len(FACTORS))))
    print()
    print("  Read the f = 2 / f = 3.5 columns: that is the challenge's n_real regime.")


if __name__ == "__main__":
    main()
