"""516: Stouffer ablations on H1 (uncentred, cap 3, |Z| >= 0.5, all genes and uncalled-only)."""
import numpy as np
import importlib.util as u
s = u.spec_from_file_location("v", "vote516.py"); V = u.module_from_spec(s); s.loader.exec_module(V)
from common516 import price, summary_line

d = V.d; calls = d["calls"]
sets = {"K": ["K"], "H": ["H"], "E": ["E"], "R": ["R"], "S8": ["S8"], "S48": ["S48"], "V": ["V"],
        "KHE": ["K", "H", "E"], "KH": ["K", "H"], "KE": ["K", "E"], "cd4x3": ["R", "S8", "S48"],
        "K+cd4x3": ["K", "R", "S8", "S48"], "KH+cd4x3": ["K", "H", "R", "S8", "S48"],
        "KHE+S48": ["K", "H", "E", "S48"], "KHE+cd4x3": ["K", "H", "E", "R", "S8", "S48"],
        "KHEV+cd4x3": ["K", "H", "E", "V", "R", "S8", "S48"]}
for nm, keys in sets.items():
    Z = V.stouffer(keys, False, 3.0)
    for region, rm in (("all", np.ones_like(calls)), ("uncalled", ~calls)):
        C = np.where(rm, V.gate(Z, 0.5), 0.0)
        r = price(d, C, f"stz[{nm}] |Z|>=0.5 {region}", subsets=V.SUBS, verbose=False)
        for sn in V.SUBS:
            print(summary_line(r[sn]), flush=True)
# weighted: CD4 conditions as one source (mean z of the three) with weight w
for w in (1.0, 1.5, 2.0):
    Ss = [V.std_score(V.SRC[k], False, 3.0) for k in ("K", "H", "E")]
    c4 = np.nanmean(np.stack([V.std_score(V.SRC[k], False, 3.0) for k in ("R", "S8", "S48")]), 0)
    parts = Ss + [c4]; ws = [1, 1, 1, w]
    num = sum(wi * np.nan_to_num(p) for wi, p in zip(ws, parts)); den = np.sqrt(sum(wi**2 * np.isfinite(p) for wi, p in zip(ws, parts)))
    Z = np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)
    for region, rm in (("all", np.ones_like(calls)), ("uncalled", ~calls)):
        C = np.where(rm, V.gate(Z, 0.5), 0.0)
        r = price(d, C, f"stz[KHE + {w}*cd4mean] |Z|>=0.5 {region}", subsets=V.SUBS, verbose=False)
        for sn in V.SUBS:
            print(summary_line(r[sn]), flush=True)
