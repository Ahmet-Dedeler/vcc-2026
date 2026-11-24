"""812 step 1c — the truth's significant genes: UP/DOWN share and CLEAN4 / emitted sign accuracy by control CPM, called vs
uncalled (a Wilcoxon-sparsity selection effect would make low-CPM significant genes lean UP). usage: uv run python
experiments/812-de-signs/sigcpm812.py -> work/sigcpm812.json"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent)); import surf812 as SF
out = {}
for surf in ("h1test", "h1tv"):
    S = SF.load(surf, 0); adj, sig, call = SF.masks(S); rl, ol = S["rl"], S["ol"]; T, G = rl.shape
    c4 = S["c4"]; cov = np.isfinite(c4); c4s = np.where(c4 > 0.5, -1.0, 1.0); cpm = S["cpm"]
    mc = S["key"] != 0
    R = {}
    for lo, hi in ((5, 10), (10, 20), (20, 50), (50, 100), (100, 300), (300, 1e9)):
        g = (cpm >= lo) & (cpm < hi)
        for cn, cm in (("called", call), ("uncalled", ~call)):
            m = sig & adj & g[None, :] & cm
            n = int(m.sum())
            if not n: continue
            mm = m & cov[None, :]
            R[f"cpm{lo}-{hi:g}|{cn}"] = dict(n=n, up_share=round(float((rl[m] > 0).mean()), 4),
                acc_emit=round(float((np.sign(ol[m]) == np.sign(rl[m])).mean()), 4),
                acc_c4=round(float((c4s[None, :].repeat(T, 0)[mm] == np.sign(rl[mm])).mean()), 4) if mm.sum() else None,
                acc_key=round(float((np.sign(S["key"][m & mc]) == np.sign(rl[m & mc])).mean()), 4) if (m & mc).sum() else None,
                mean_abs_lfc=round(float(np.abs(rl[m]).mean()), 3), share_of_nmae_den=None)
    den = np.abs(np.where(sig & adj, rl, 0)).sum()
    for k in R:
        lo, hi = k.split("|")[0][3:].split("-"); lo = float(lo); hi = float(hi)
        g = (cpm >= lo) & (cpm < hi); cm = call if k.endswith("|called") else ~call
        R[k]["share_of_nmae_den"] = round(float(np.abs(np.where(sig & adj & g[None, :] & cm, rl, 0)).sum() / den), 4)
    out[surf] = R
    print(surf); [print("  ", k, v) for k, v in R.items()]
json.dump(out, open(SF.HERE / "work/sigcpm812.json", "w"), indent=1)
