"""516: which targets lose reach under the uncalled-only Stouffer-6 candidate, and where their flipped genes sit."""
import numpy as np
import importlib.util as u
s = u.spec_from_file_location("v", "vote516.py"); V = u.module_from_spec(s); s.loader.exec_module(V)
from common516 import hybrid, members

d = V.d; G = d["genes"]; T = d["targets"]; calls = d["calls"]
Z = V.stouffer(["K", "H", "E", "R", "S8", "S48"], False, 3.0)
C = np.where(~calls, V.gate(Z, 0.5), 0.0)
H = hybrid(d, C)
b = d["_base_members"] if "_base_members" in d else members(d)
c = members(d, H)
dr = c["reach"] - b["reach"]
print("targets with reach change:")
for k in np.argsort(dr):
    if abs(np.nan_to_num(dr[k])) < 1e-12: continue
    rs = d["rs"][k]; order = d["orders"][k]
    seq = order[rs[order]]; seq = seq[np.isfinite(d["rl"][k][seq]) & (d["rl"][k][seq] != 0)]
    walkb = (d["prod"][k][seq] == d["truth"][k][seq]).astype(int)
    walkc = (H[k][seq] == d["truth"][k][seq]).astype(int)
    called = d["calls"][k][seq].astype(int)
    print(f"{T[k]:9s} n_real {d['n_real'][k]:5d} reach {b['reach'][k]:.3f} -> {c['reach'][k]:.3f}  "
          f"walk prod {''.join(map(str, walkb[:25]))}\n{'':44s}walk cand {''.join(map(str, walkc[:25]))}\n"
          f"{'':44s}called    {''.join(map(str, called[:25]))}  pos {V.pos[k][seq[:6]].tolist()}")
