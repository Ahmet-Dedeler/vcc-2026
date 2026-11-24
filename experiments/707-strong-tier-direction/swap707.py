#!/usr/bin/env python3
"""707: the strong-tier direction swap on a stack dump (one construction, used both for the A/B/C arm and on the day).

On the top tertile of targets by z_src (the production donors' robust strength rank, 522's feature; computable for any
panel from the donor banks and the bundle's pert_counts.csv), each bulk row's NON-PANEL part is blended toward the
production donors' equal-weight unit mean U_p = Σ_j r_j/|r_j| (raw, unwhitened delta rows from data/donors/<bank>,
non-panel genes only), norm-preserving:
    b'_p = |b_p| · normalise((1 − β) b_p/|b_p| + β U_p/|U_p|)     (non-panel columns; panel columns untouched).
Targets without any donor row, and the lower two tertiles, keep their rows byte for byte. The carriers (the DE design),
perts, genes and env are copied unchanged, so this is a bulk-only change and every row keeps its norm: the P2g tercile
split (v = |bulk|/|carriers|) is unchanged by construction.

Held-out basis: experiments/706-strong-tier-budget (dirswap706: top tier, β 0.25).

usage:
  uv run python experiments/707-strong-tier-direction/swap707.py --dump-in D --dump-out O --bundle B [--beta 0.25]
  uv run python experiments/707-strong-tier-direction/swap707.py --ka --dump-in D --bundle data/challenge --dump-out O
      (KA-β0: β 0 reproduces the input bulk np.array_equal; KA-z: z on A/B/C == 518's abc518.zprod(targets, 'src'))
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/522-october-bulk-key"))
from def_strength_key522 import zsrc  # noqa: E402  (522's day-of z_src; its --ka checks it against 518 on A/B/C)
from strength522 import BANKS, panel_of  # noqa: E402

TOP = 2          # the top tertile
BANK_ORDER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def donor_unit_mean(perts, genes, nonpanel):
    U = np.zeros((len(perts), len(genes)))
    cover = {}
    for bk in BANK_ORDER:
        D = Path(BANKS[bk])
        dg = (D / "genes.txt").read_text().split()
        dp = (D / "perts.txt").read_text().split()
        delta = np.load(D / "delta.npy", mmap_mode="r")
        gpos = {g: i for i, g in enumerate(dg)}
        ppos = {p: i for i, p in enumerate(dp)}
        gc = np.array([gpos.get(g, -1) for g in genes])
        ok = (gc >= 0) & nonpanel
        n = 0
        for i, t in enumerate(perts):
            if t not in ppos:
                continue
            r = np.zeros(len(genes))
            r[ok] = np.asarray(delta[ppos[t]], np.float64)[gc[ok]]
            nr = np.linalg.norm(r)
            if nr > 0:
                U[i] += r / nr
                n += 1
        cover[bk] = n
    return U, cover


def swap(dump_in: Path, dump_out: Path, bundle: Path, beta: float) -> dict:
    D = np.load(dump_in, allow_pickle=False)
    perts = [str(x) for x in D["perts"]]
    genes = [str(x) for x in D["genes"]]
    panel = panel_of(bundle)
    assert sorted(perts) == sorted(panel) and len(perts) == 300, "dump perts != the bundle's panel"
    pset = set(panel)
    nonpanel = np.asarray([g not in pset for g in genes], dtype=bool)
    bulk = np.asarray(D["bulk"], dtype=np.float64)
    U, cover = donor_unit_mean(perts, genes, nonpanel)
    have = np.linalg.norm(U, axis=1) > 0
    z, zcover = zsrc(perts, panel)
    tier = np.full(len(perts), -1)
    r = np.argsort(np.argsort(z[have], kind="stable"), kind="stable")
    tier[np.flatnonzero(have)] = np.minimum((r * 3) // int(have.sum()), 2)
    sel = tier == TOP
    out = bulk.copy()
    nb = np.linalg.norm(bulk[:, nonpanel], axis=1)
    nu = np.linalg.norm(U[:, nonpanel], axis=1)
    moved = 0
    if beta != 0.0:
        for i in np.flatnonzero(sel & (nb > 0) & (nu > 0)):
            m = (1.0 - beta) * bulk[i, nonpanel] / nb[i] + beta * U[i, nonpanel] / nu[i]
            out[i, nonpanel] = nb[i] * m / np.linalg.norm(m)
            moved += 1
    dump_out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dump_out, bulk=out, carriers=np.asarray(D["carriers"]), perts=np.asarray(D["perts"]),
                        genes=np.asarray(D["genes"]), env=np.asarray(D["env"]))
    ch = out - bulk
    rows = np.flatnonzero(sel)
    cosr = [float(out[i, nonpanel] @ bulk[i, nonpanel] / max(nb[i] * np.linalg.norm(out[i, nonpanel]), 1e-300)) for i in rows
            if nb[i] > 0]
    side = {"out": str(dump_out), "sha256": sha256(dump_out), "dump_in": str(dump_in), "dump_in_sha256": sha256(dump_in),
            "bundle": str(bundle), "beta": beta, "rule": "top z_src tertile: non-panel row blended toward the production "
            "donors' unit mean, norm-preserving; panel columns, carriers, perts, genes, env unchanged",
            "n_have_donor_rows": int(have.sum()), "n_top": int(sel.sum()), "n_moved": moved, "donor_cover": cover,
            "z_cover": zcover, "z_top_min": float(z[sel].min()) if sel.any() else None,
            "row_norm_max_rel_change": float(np.max(np.abs(np.linalg.norm(out[:, nonpanel], axis=1) - nb) / np.maximum(nb, 1e-300))),
            "moved_row_cos_to_input_q": np.round(np.quantile(cosr, [0, .25, .5, .75, 1]), 4).tolist() if cosr else None,
            "change_energy_share_of_bulk": float((ch * ch).sum() / (bulk * bulk).sum()),
            "top_targets": [perts[i] for i in rows]}
    Path(str(dump_out) + ".swap707.json").write_text(json.dumps(side, indent=1))
    return side


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dump-in", type=Path, required=True)
    ap.add_argument("--dump-out", type=Path, required=True)
    ap.add_argument("--bundle", type=Path, required=True)
    ap.add_argument("--beta", type=float, default=0.25)
    ap.add_argument("--ka", action="store_true")
    a = ap.parse_args()
    if a.ka:
        s0 = swap(a.dump_in, a.dump_out, a.bundle, 0.0)
        ok0 = np.array_equal(np.load(a.dump_out)["bulk"], np.load(a.dump_in)["bulk"])
        print(f"KA-beta0 beta 0 reproduces the input bulk: {'PASS' if ok0 else 'FAIL'} ({s0['sha256'][:12]})")
        D = np.load(a.dump_in, allow_pickle=False)
        perts = [str(x) for x in D["perts"]]
        sys.path.insert(0, str(ROOT / "experiments/518-bulk-mse-gap"))
        from abc518 import zprod
        z, _ = zsrc(perts, panel_of(a.bundle))
        okz = np.array_equal(z, zprod(perts, "src"))
        print(f"KA-z     z on this panel == abc518.zprod(targets, 'src'): {'PASS' if okz else 'FAIL'}")
        a.dump_out.unlink()
        Path(str(a.dump_out) + ".swap707.json").unlink()
        raise SystemExit(0 if (ok0 and okz) else 1)
    side = swap(a.dump_in, a.dump_out, a.bundle, a.beta)
    print(json.dumps({k: v for k, v in side.items() if k != "top_targets"}, indent=1))


if __name__ == "__main__":
    main()
