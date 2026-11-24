#!/usr/bin/env python3
"""612 (council 600, [21dd4c], 2026-10-02) - does the KEY-GATED sign rule transfer to held-out cancer lines?

H1 finding (610/gatesign610, dosecurve610): on uncalled gate genes our target-specific key beats CLEAN4 where |key| is high
and loses where |key| is low; a hybrid (CLEAN4 where |key| is in the target's bottom half of nonzero |key|, key elsewhere)
prices ~3x 496's plain CLEAN4 hook on H1 nmae. The half-split was chosen on H1 after seeing the table (forking path), so
this checks it on 480's cross-study cancer folds before anything goes near a slot (AGENTS 11.1c).

Surfaces (480's gates; Orion truths carry a technical blob -> RELATIVE comparisons only): h1 (reference), orion_hct116,
orion_hek293t, k562_gwps. Key proxy = 480's same-target donor consensus from the OTHER donors (k562_gwps, orion_hct116,
orion_hek293t; log2 pseudobulk ratios, pc 1 CPM); CLEAN4 P_down from 371's banks, minus the held-out line's own banks
(K562 fold: Nadig Jurkat + HepG2 only), as 537's foldsign537 does. Magnitude = 480's shipped MDE law (min(0.1/sqrt(cnt),
0.1) for cnt >= 0.3) for every sign source, so only the SIGN differs between constructions. Analysis only, < 2 GB.
Writes work/kgh612.json.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument"))
sys.path.insert(0, str(ROOT / "experiments/537-arm-pipeline"))
import fold480 as F            # noqa: E402  (load_surface, mde_add, DON)
from common480 import boot_delta, load_donor_rows, nmae_rows  # noqa: E402
import foldsign537 as FS        # noqa: E402  (clean4_P)
C4ALL = ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"]
C4BANKS = {"h1": C4ALL, "orion_hct116": C4ALL, "orion_hek293t": C4ALL, "k562_gwps": ["nadig_jurkat", "nadig_hepg2"]}
out = {}
for name in ["h1", "orion_hct116", "orion_hek293t", "k562_gwps"]:
    S = F.load_surface(name)
    T, cnt, genes, targets = S["T"], S["cnt"], S["genes"], S["targets"]
    R = [load_donor_rows(d, targets, genes, pc_cpm=1.0)[0] for d in S["donors"]]
    cons = np.nanmean(np.stack(R), 0)
    cons0 = np.nan_to_num(cons)
    P4 = FS.clean4_P(genes, C4BANKS[name])
    s4 = np.where(np.isfinite(P4), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
    S4 = np.broadcast_to(s4, T.shape)
    sc = np.sign(cons0)
    # within-target |cons| rank among nonzero entries
    q = np.full(T.shape, np.nan)
    for i in range(T.shape[0]):
        a = np.abs(cons0[i]); nz = a > 0
        if nz.sum() > 10:
            q[i, nz] = (np.argsort(np.argsort(a[nz])) + 0.5) / nz.sum()
    gate = np.isfinite(T)
    res = dict(donors=S["donors"], clean4_banks=C4BANKS[name], n_targets=int(gate.any(1).sum()), gate_pairs=int(gate.sum()))
    # sign accuracy by |cons| quantile, on gate genes where BOTH sources carry a sign
    both = gate & (sc != 0) & (S4 != 0)
    tr = np.sign(np.where(gate, T, 0.0))
    bins = []
    for lo, hi in [(0, .5), (.5, .8), (.8, .95), (.95, 1.01)]:
        m = both & (q >= lo) & (q < hi)
        n = int(m.sum())
        bins.append(dict(lo=lo, hi=hi, n=n, key_acc=float((sc[m] == tr[m]).mean()) if n else None,
                         c4_acc=float((S4[m] == tr[m]).mean()) if n else None))
    res["by_key_quantile"] = bins
    # constructions: same MDE-law magnitude, different signs
    add = F.mde_add(cnt)[None, :]
    lowk = np.where(np.isfinite(q), q < 0.5, True)            # no key rank -> treat as low-key (CLEAN4 if it has a sign)
    hyb = np.where(lowk & (S4 != 0), S4, np.where(sc != 0, sc, S4))
    hyb_strict = np.where(lowk & (S4 != 0), S4, sc)            # CLEAN4 only on low-key genes, key elsewhere (0 if no key)
    P = {"mde_key": sc * add, "mde_clean4": S4 * add, "mde_hybrid": hyb * add, "mde_hybrid_strict": hyb_strict * add}
    rows = {k: nmae_rows(v, T) for k, v in P.items()}
    res["nmae"] = {k: float(np.nanmean(v)) for k, v in rows.items()}
    res["delta_vs_key"] = {k: boot_delta(rows[k], rows["mde_key"]) for k in P if k != "mde_key"}
    res["delta_vs_clean4"] = {k: boot_delta(rows[k], rows["mde_clean4"]) for k in P if k != "mde_clean4"}
    out[name] = res
    print(f"\n== {name}: targets {res['n_targets']}, gate pairs {res['gate_pairs']}, donors {S['donors']}, CLEAN4 banks {C4BANKS[name]}")
    for b in bins:
        if b['n']:
            print(f"   |key| quantile [{b['lo']:.2f},{b['hi']:.2f}): n={b['n']:6d}  key acc {b['key_acc']:.3f}  CLEAN4 acc {b['c4_acc']:.3f}")
    print("   nmae (MDE-law magnitude): " + "  ".join(f"{k} {v:.5f}" for k, v in res['nmae'].items()))
    for k, (m, lo, hi, pge0) in res["delta_vs_key"].items():
        print(f"   {k:18s} - mde_key    {m:+.5f} [{lo:+.5f}, {hi:+.5f}]")
    for k, (m, lo, hi, pge0) in res["delta_vs_clean4"].items():
        print(f"   {k:18s} - mde_clean4 {m:+.5f} [{lo:+.5f}, {hi:+.5f}]")
json.dump(out, open(ROOT / "experiments/612-key-gated-hook/work/kgh612.json", "w"), indent=1)
