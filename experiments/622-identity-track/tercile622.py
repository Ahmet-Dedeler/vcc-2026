#!/usr/bin/env python3
"""622l — where does the identity gain live: weak or strong targets? On 636's gates, gate-sign Δ of
S = z(C4 - .5) + 0.5 z(I) (Flex references) vs CLEAN4 by target-strength tercile (strength = the target's gate size on
that surface; reported, truth-derived, for interpretation only) and by the donor-consensus strength (621's mean |Zbal|,
truth-free). Analysis only, < 2 GB."""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402

I = Inst()
prof = {}
for s, f in SURF.items():
    z = np.load(G480 / f, allow_pickle=True); prof[s] = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
z = np.load(HERE / "work/abc_cpm622.npz", allow_pickle=True)
for k in ("ctxA", "ctxB", "ctxC"):
    prof[k] = dict(zip(z["genes"].astype(str), z[k].astype(float)))
out = {}
for s in I.S:
    D = I.S[s]
    refs = [k for k in ("ctxA", "ctxB", "ctxC", "h1") if k != s]
    common = sorted(set(prof[s]).intersection(*(set(prof[r]) for r in refs)))
    Iv = I.gene_vec(s, common, np.log1p([prof[s][g] for g in common]) - np.mean([np.log1p([prof[r][g] for g in common]) for r in refs], 0))
    zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
    C4, _ = clean4_lolo(I, s); zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
    Sv = zc + 0.5 * zi
    gate, dn = D["gate"], D["lfc"] < 0
    both = np.isfinite(zc) & np.isfinite(Sv)
    m = gate & both[None, :]
    okS = (((Sv > 0)[None, :] == dn) & m).sum(1).astype(float); okC = (((zc > 0)[None, :] == dn) & m).sum(1).astype(float)
    n = m.sum(1).astype(float)
    size = gate.sum(1)
    zb = np.load(ROOT / f"experiments/621-supervised-de/work/pred-{s}.npz", allow_pickle=True)["Zbal"]
    kstr = np.nanmean(np.abs(zb), 1)
    res = {}
    for name, v in (("gate_size", size), ("donor_strength", kstr)):
        ok = n > 0
        q = np.nanquantile(v[ok], [1 / 3, 2 / 3])
        for t, (lo, hi) in enumerate(((-np.inf, q[0]), (q[0], q[1]), (q[1], np.inf))):
            sel = ok & (v > lo) & (v <= hi)
            res[f"{name}_T{t + 1}"] = dict(n_targets=int(sel.sum()), pairs=int(n[sel].sum()),
                                           d=float((okS[sel].sum() - okC[sel].sum()) / max(n[sel].sum(), 1)),
                                           per_target_mean_d=float(np.nanmean((okS[sel] - okC[sel]) / np.maximum(n[sel], 1))))
    out[s] = res
    print(s, " | ".join(f"{k}: {v['d']:+.4f} (per-target {v['per_target_mean_d']:+.4f}, n {v['n_targets']})" for k, v in res.items()))
(HERE / "work/tercile622.json").write_text(json.dumps(out, indent=1))
