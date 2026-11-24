"""516: where does the Stouffer-6 candidate gain on H1? Pool cells split by called/uncalled x block/rest, by |Z| bin."""
import numpy as np
from common516 import WORK, known_answers, load_h1
import vote516 as V

d = V.d; T, G = d["targets"], d["genes"]
pool, truth, prod, calls = d["pool"], d["truth"], d["prod"], d["calls"]
blk = np.zeros(len(G), bool); blk[np.argsort(-d["ctrl"])[:2200]] = True
Z = V.stouffer(["K", "H", "E", "R", "S8", "S48"], False, 3.0)
S = np.sign(np.nan_to_num(Z))
M0 = pool & (prod != 0) & (S != 0)
for subn, sub in V.SUBS.items():
    print(f"--- {subn}")
    for cn, cm in (("called", calls), ("uncalled", ~calls)):
        for rn, rm in (("block", blk[None, :]), ("rest", ~blk[None, :])):
            m = M0 & cm & rm & sub[:, None]
            print(f"  {cn:8s} {rn:5s} n {m.sum():6d}  prod {(prod[m]==truth[m]).mean():.4f}  stz6 {(S[m]==truth[m]).mean():.4f}"
                  f"  diff {(S[m]==truth[m]).mean()-(prod[m]==truth[m]).mean():+.4f}  agree {(S[m]==prod[m]).mean():.3f}")
    for lo, hi in ((0, 0.5), (0.5, 1), (1, 1.5), (1.5, 2), (2, 3), (3, 99)):
        m = M0 & (np.abs(Z) >= lo) & (np.abs(Z) < hi) & sub[:, None]
        print(f"  |Z| [{lo},{hi}) n {m.sum():6d}  prod {(prod[m]==truth[m]).mean():.4f}  stz6 {(S[m]==truth[m]).mean():.4f}"
              f"  diff {(S[m]==truth[m]).mean()-(prod[m]==truth[m]).mean():+.4f}")
