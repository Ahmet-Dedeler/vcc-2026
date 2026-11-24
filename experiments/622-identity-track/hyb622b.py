#!/usr/bin/env python3
"""622b2 — 622's identity low half re-checked in 612b's DELIVERY-FAITHFUL shape (21dd4c's kgh612b.py), as asked.

The real emitter delivers base_lfc (key-signed design) + sign_hook x MDE add on uncalled genes; the hook only sets the
add's direction. 480's `small_base+mde` shape: P = BETA * cons + s * add. Signs s compared (same base, same add):
  nohook    key sign everywhere
  plain     CLEAN4 (|P - .5| >= .05) where covered, else key sign                       (496's plain hook)
  kg        key-gated: CLEAN4 on the weak-key half, key sign on the strong half       (612)
  plainI    S = z(C4 - .5) + lam z(I) sign wherever S is defined, else key sign
  kgI       key-gated with sign(S) on the weak half                                     (622's recommended delivery)
lam in {0.5, 1}; I with FLEX references (A/B/C [+ h1]), via hyb622's profiles/identity. Δ = mean(ref - X), + = X better.
Run with cwd = experiments/480-nmae-instrument. Analysis only, < 2 GB.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument")); sys.path.insert(0, str(ROOT / "experiments/537-arm-pipeline"))
sys.path.insert(0, str(HERE))
import fold480 as F  # noqa: E402
from common480 import boot_delta, load_donor_rows, nmae_rows, load_h1  # noqa: E402
import foldsign537 as FS  # noqa: E402
import hyb622 as H  # noqa: E402

C4ALL = ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"]
C4B = {"h1": C4ALL, "orion_hct116": C4ALL, "orion_hek293t": C4ALL, "k562_gwps": ["nadig_jurkat", "nadig_hepg2"]}
h = load_h1(); G = np.isfinite(h["T"])
L0 = np.nan_to_num(np.asarray(np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz")["lfc"], float))
c3 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, h["targets"], h["genes"], pc_cpm=1.0)[0] for d in F.DON]), 0))
BETA = float(np.median(np.abs(L0[G])) / np.median(np.abs(c3[G])))
print(f"BETA {BETA:.4f} (480's construction; must equal kgh612b's)")
prof = H.profiles()
out = {"BETA": BETA}
for name in C4B:
    S = F.load_surface(name)
    T, cnt, genes, targets = S["T"], S["cnt"], [str(g) for g in S["genes"]], S["targets"]
    cons0 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, np.array(genes), pc_cpm=1.0)[0] for d in S["donors"]]), 0))
    P4 = FS.clean4_P(np.array(genes), C4B[name])
    s4 = np.where(np.isfinite(P4) & (np.abs(P4 - 0.5) >= 0.05), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
    S4 = np.broadcast_to(s4, T.shape); sc = np.sign(cons0); add = F.mde_add(cnt)[None, :]
    weak = np.zeros(T.shape, bool)
    for i in range(T.shape[0]):
        a = np.abs(cons0[i]); nz = a > 0
        weak[i] = (~nz) | (a < (np.quantile(a[nz], 0.5) if nz.any() else 0.0))
    base = BETA * cons0
    zc = np.where(np.isfinite(P4), (P4 - 0.5) / np.nanstd(P4 - 0.5), np.nan)
    refs = [k for k in ("ctxA", "ctxB", "ctxC", "h1") if k != name]
    zi = H.identity(prof, name, refs, genes)
    P = {"nohook": base + sc * add,
         "plain": base + np.where(S4 != 0, S4, sc) * add,
         "kg": base + np.where(weak & (S4 != 0), S4, sc) * add}
    for lam in (0.5, 1.0):
        Sv = np.where(np.isfinite(zc), zc + lam * np.nan_to_num(zi), np.nan_to_num(zi))   # + -> down
        sI = np.where(np.isfinite(Sv) & (Sv != 0), np.where(Sv > 0, -1.0, 1.0), 0.0)
        SI = np.broadcast_to(sI, T.shape)
        P[f"plainI_l{lam}"] = base + np.where(SI != 0, SI, sc) * add
        P[f"kgI_l{lam}"] = base + np.where(weak & (SI != 0), SI, sc) * add
    rows = {k: nmae_rows(v, T) for k, v in P.items()}
    res = {"nmae": {k: float(np.nanmean(v)) for k, v in rows.items()}}
    comps = [("plain", "nohook"), ("kg", "nohook"), ("kg", "plain")]
    for lam in (0.5, 1.0):
        comps += [(f"plainI_l{lam}", "plain"), (f"kgI_l{lam}", "kg"), (f"kgI_l{lam}", "nohook"), (f"kgI_l{lam}", "plain")]
    for a, b in comps:
        d = boot_delta(rows[a], rows[b]); res[f"{a}-vs-{b}"] = [float(x) for x in d[:3]]
    out[name] = res
    print(f"\n== {name}")
    for a, b in comps:
        d = res[f"{a}-vs-{b}"]
        print(f"   {a:14s} vs {b:8s} {d[0]:+.5f} [{d[1]:+.5f},{d[2]:+.5f}]")
(HERE / "work/hyb622b.json").write_text(json.dumps(out, indent=1))
