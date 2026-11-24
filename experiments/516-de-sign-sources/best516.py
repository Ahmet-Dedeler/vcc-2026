"""516: per-target table for the best H1 candidate (Stouffer-6 uncentred, cap 3, |Z| >= 0.5, genes production does not
call; fid pinned) and its CD4-reliability-gated variant, written to work/best516-per-target.csv + a summary json."""
import csv
import json
import numpy as np
import importlib.util as u
s = u.spec_from_file_location("v", "vote516.py"); V = u.module_from_spec(s); s.loader.exec_module(V)
import fold516 as FD
from common516 import WORK, hybrid, members, price

d = V.d; T = d["targets"]; calls = d["calls"]
six = ["K", "H", "E", "R", "S8", "S48"]
Z = V.stouffer(six, False, 3.0)
C = np.where(~calls & (np.abs(np.nan_to_num(Z)) >= 0.5), np.sign(np.nan_to_num(Z)), 0.0)
H = hybrid(d, C)
b = members(d); c = members(d, H)
m = d["pool"] & (d["prod"] != 0)
with open(WORK / "best516-per-target.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["target", "n_real", "band", "essential_covered", "pool_cells", "prod_acc", "cand_acc", "flips",
                "nmae_prod", "nmae_cand", "reach_prod", "reach_cand"])
    for k, t in enumerate(T):
        n = int(m[k].sum())
        w.writerow([t, int(d["n_real"][k]), ["<20", "20-100", ">=100"][d["band"][k]], bool(~V.SUBS["no_ess"][k]), n,
                    f"{(d['prod'][k][m[k]] == d['truth'][k][m[k]]).mean():.4f}" if n else "",
                    f"{(H[k][m[k]] == d['truth'][k][m[k]]).mean():.4f}" if n else "",
                    int((H[k] != d["prod"][k]).sum()),
                    f"{b['nmae'][k]:.5f}", f"{c['nmae'][k]:.5f}", f"{b['reach'][k]:.4f}", f"{c['reach'][k]:.4f}"])
r = price(d, C, "BEST stz6 |Z|>=0.5 uncalled", subsets=V.SUBS, verbose=True)
json.dump({sn: dict(acc=v["acc"], acc_w=v["acc_w"], acc_tmean=v["acc_tmean"], members=v["members"],
                    coverage=v["coverage_pool"], flips_pool=v["flips_pool"], flips_all=v["flips_all"])
           for sn, v in r.items()}, open(WORK / "best516.json", "w"), indent=1, default=float)
print("wrote work/best516-per-target.csv and work/best516.json")
