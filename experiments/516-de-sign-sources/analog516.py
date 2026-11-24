"""516: which simple donor construction best reproduces production's emitted signs on H1 (to serve as the production
analog on the cancer folds, where no emitter output exists)? Agreement with production's sign and pool accuracy,
on non-block genes (carrier signs) and block genes (CLEAN4 prior with q80 protection)."""
import numpy as np
import importlib.util as u
s = u.spec_from_file_location("v", "vote516.py"); V = u.module_from_spec(s); s.loader.exec_module(V)
from common516 import ROOT

d = V.d; z = V.z; G = d["genes"]; T = d["targets"]
pool, truth, prod = d["pool"], d["truth"], d["prod"]
blk = np.zeros(len(G), bool); blk[np.argsort(-d["ctrl"])[:2200]] = True
x = {k: z[f"{V.SRC[k]}__x"].astype(float) if not V.SRC[k].startswith("cd4_") else z[f"{V.SRC[k]}__log_fc"].astype(float)
     for k in V.SRC}
cd4m = np.nanmean(np.stack([x["R"], x["S8"], x["S48"]]), 0)
pan = np.stack([x["K"], x["H"], x["E"]])
A = {"mean KHE": np.nanmean(pan, 0),
     "mean KHE + 0.34 cd4": np.nan_to_num(np.nanmean(pan, 0)) + 0.34 * np.nan_to_num(cd4m),
     "mean KHE + 1 cd4": np.nan_to_num(np.nanmean(pan, 0)) + np.nan_to_num(cd4m),
     "stz KHE": V.stouffer(["K", "H", "E"], False, 3.0),
     "stz KHE+cd4x3": V.stouffer(["K", "H", "E", "R", "S8", "S48"], False, 3.0),
     "K only": x["K"]}
cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
bn = list(cp["bank_names"].astype(str)); print("P_banks:", bn)
gi = {g: i for i, g in enumerate(cp["genes"].astype(str))}
c4 = np.array([0.5 - cp["P_down"][gi[g]] if g in gi else np.nan for g in G])
for nm, R in A.items():
    S = np.sign(np.nan_to_num(R))
    for rn, rm in (("rest", ~blk), ("block", blk)):
        allm = (prod != 0) & (S != 0) & rm[None, :]
        pm = allm & pool
        print(f"{nm:22s} {rn:5s} agree w/ prod (all genes) {(S[allm]==prod[allm]).mean():.3f}  pool: agree "
              f"{(S[pm]==prod[pm]).mean():.3f} acc {(S[pm]==truth[pm]).mean():.4f} vs prod {(prod[pm]==truth[pm]).mean():.4f}")
# PPprior-like block emulation: CLEAN4 sign on block except top-20% |analog| per target
for nm in ("mean KHE + 0.34 cd4", "stz KHE+cd4x3"):
    R = np.nan_to_num(A[nm]); E = np.sign(R).copy()
    for k in range(len(T)):
        v = np.abs(R[k, blk]); thr = np.quantile(v, 0.8)
        flip = v <= thr
        cols = np.flatnonzero(blk)[flip]
        E[k, cols] = np.where(np.isfinite(c4[cols]), np.sign(c4[cols]), E[k, cols])
    m = (prod != 0) & (E != 0)
    for rn, rm in (("all", np.ones(len(G), bool)), ("block", blk)):
        mm = m & rm[None, :]; pm = mm & pool
        print(f"EMULATION [{nm} + CLEAN4 q80 block] {rn}: agree w/ prod {(E[mm]==prod[mm]).mean():.3f}; pool acc "
              f"{(E[pm]==truth[pm]).mean():.4f} vs prod {(prod[pm]==truth[pm]).mean():.4f}")
