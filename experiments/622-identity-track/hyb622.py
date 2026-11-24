#!/usr/bin/env python3
"""622c — HYBRID + identity on 612's folds (d4cd5f's (a)); nmae with the shipped MDE-law magnitude, so only signs differ.

Exactly 612's construction (kgh612.py: 480 surfaces, key proxy = same-target consensus of the OTHER donors, CLEAN4 minus
the held-out line's banks, per-target median-|key| gate), with the LOW-key half's sign changed from sign(CLEAN4) to
sign(z(P_C4 - .5) + lam z(I)) (identity alone where CLEAN4 has no sign). I = log1p ctrl CPM of the surface - mean of the
reference profiles' log1p ctrl CPM, with references FLEX = {A, B, C} (+ h1 when the surface is not h1; the D/E/F day-of
form) or MIXED (639's: every other profile). Read on 612's min-max-regret terms: delta vs mde_key and vs 612's hybrid,
target bootstrap (common480.boot_delta: mean(ref - X), positive = X better).
Run with cwd = experiments/480-nmae-instrument (fold480 reads relative paths). Analysis only, < 2 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument"))
sys.path.insert(0, str(ROOT / "experiments/537-arm-pipeline"))
import fold480 as F  # noqa: E402
from common480 import boot_delta, load_donor_rows, nmae_rows  # noqa: E402
import foldsign537 as FS  # noqa: E402

C4ALL = ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"]
C4BANKS = {"h1": C4ALL, "orion_hct116": C4ALL, "orion_hek293t": C4ALL, "k562_gwps": ["nadig_jurkat", "nadig_hepg2"]}
GATEF = {"h1": "gate-h1_full_ka.npz", "orion_hct116": "gate-orion_hct116.npz", "orion_hek293t": "gate-orion_hek293t.npz",
         "k562_gwps": "gate-k562_gwps.npz"}
G480 = ROOT / "experiments/480-nmae-instrument/work"


def profiles() -> dict[str, dict[str, float]]:
    prof = {}
    for s, f in GATEF.items():
        z = np.load(G480 / f, allow_pickle=True)
        prof[s] = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
    z = np.load(HERE / "work/abc_cpm622.npz", allow_pickle=True)
    g = z["genes"].astype(str)
    for k in ("ctxA", "ctxB", "ctxC"):
        prof[k] = dict(zip(g, z[k].astype(float)))
    return prof


def identity(prof, s, refs, genes):
    common = set(prof[s]).intersection(*(set(prof[r]) for r in refs))
    v = np.full(len(genes), np.nan)
    for j, g in enumerate(genes):
        if g in common:
            v[j] = np.log1p(prof[s][g]) - np.mean([np.log1p(prof[r][g]) for r in refs])
    return (v - np.nanmean(v)) / np.nanstd(v)


def main() -> None:
    prof = profiles()
    out = {}
    for name in ["h1", "orion_hct116", "orion_hek293t", "k562_gwps"]:
        S = F.load_surface(name)
        T, cnt, genes, targets = S["T"], S["cnt"], [str(g) for g in S["genes"]], S["targets"]
        R = [load_donor_rows(d, targets, np.array(genes), pc_cpm=1.0)[0] for d in S["donors"]]
        cons0 = np.nan_to_num(np.nanmean(np.stack(R), 0))
        P4 = FS.clean4_P(np.array(genes), C4BANKS[name])
        s4 = np.where(np.isfinite(P4), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
        zc = np.where(np.isfinite(P4), (P4 - 0.5) / np.nanstd(P4 - 0.5), np.nan)
        S4 = np.broadcast_to(s4, T.shape)
        sc = np.sign(cons0)
        q = np.full(T.shape, np.nan)
        for i in range(T.shape[0]):
            a = np.abs(cons0[i]); nz = a > 0
            if nz.sum() > 10:
                q[i, nz] = (np.argsort(np.argsort(a[nz])) + 0.5) / nz.sum()
        lowk = np.where(np.isfinite(q), q < 0.5, True)
        add = F.mde_add(cnt)[None, :]
        others = [k for k in prof if k != name]
        refsets = {"flex": [k for k in ("ctxA", "ctxB", "ctxC", "h1") if k != name], "mixed": others}
        P = {"mde_key": sc * add, "mde_clean4": S4 * add,
             "hyb612": np.where(lowk & (S4 != 0), S4, np.where(sc != 0, sc, S4)) * add}
        for rname, refs in refsets.items():
            zi = identity(prof, name, refs, genes)
            for lam in (0.5, 1.0):
                sc_low = np.where(np.isfinite(zc), zc + lam * np.nan_to_num(zi), np.nan_to_num(zi))   # positive -> down
                s_low = np.where(np.isfinite(sc_low) & (sc_low != 0), np.where(sc_low > 0, -1.0, 1.0), 0.0)
                S_low = np.broadcast_to(s_low, T.shape)
                P[f"hybI_{rname}_l{lam}"] = np.where(lowk & (S_low != 0), S_low, np.where(sc != 0, sc, S_low)) * add
                P[f"c4I_{rname}_l{lam}"] = S_low * add      # programme-only (no key) for reference
        rows = {k: nmae_rows(v, T) for k, v in P.items()}
        res = {"nmae": {k: float(np.nanmean(v)) for k, v in rows.items()},
               "vs_key": {k: boot_delta(rows[k], rows["mde_key"]) for k in P if k != "mde_key"},
               "vs_hyb612": {k: boot_delta(rows[k], rows["hyb612"]) for k in P if k != "hyb612"}}
        out[name] = res
        print(f"\n== {name} (donors {S['donors']})")
        for k in P:
            vk = res["vs_key"].get(k); vh = res["vs_hyb612"].get(k)
            s1 = f"vs key {vk[0]:+.5f} [{vk[1]:+.5f},{vk[2]:+.5f}]" if vk else ""
            s2 = f"vs hyb612 {vh[0]:+.5f} [{vh[1]:+.5f},{vh[2]:+.5f}]" if vh else ""
            print(f"   {k:22s} nmae {res['nmae'][k]:.5f}  {s1}  {s2}")
    (HERE / "work/hyb622.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
