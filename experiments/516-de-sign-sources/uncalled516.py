"""516: the Stouffer-6 sign applied only where production does not call the gene (fid pinned by construction; the
reach walk's called head untouched), optionally only off the CPM block, by |Z| threshold."""
import numpy as np
import importlib.util as u
s = u.spec_from_file_location("v", "vote516.py"); V = u.module_from_spec(s); s.loader.exec_module(V)
from common516 import price, summary_line

d = V.d; G = d["genes"]; calls = d["calls"]
blk = np.zeros(len(G), bool); blk[np.argsort(-d["ctrl"])[:2200]] = True
for keys, nm in ((["K", "H", "E", "R", "S8", "S48"], "stz6"), (["K", "H", "E", "R", "S8", "S48", "V"], "stz7")):
    Z = V.stouffer(keys, False, 3.0)
    for thr in (0.0, 0.5, 1.0, 1.5):
        for rn, rm in (("uncalled-rest", ~calls & ~blk[None, :]), ("uncalled-all", ~calls)):
            C = np.where(rm, V.gate(Z, thr), 0.0)
            r = price(d, C, f"{nm} |Z|>={thr} {rn}", subsets=V.SUBS, verbose=False)
            for sn in V.SUBS:
                print(summary_line(r[sn]), flush=True)
