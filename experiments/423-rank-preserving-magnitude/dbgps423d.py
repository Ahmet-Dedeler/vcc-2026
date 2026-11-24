import numpy as np, scipy.sparse as sp
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN, rp423
E = RN.E
r, cc, _ = E.load_reference()
z = np.load(HERE / "work/explore/block_s0_bps_58.npz")
X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
n_car = int(z["n_car"]); bal = np.arange(n_car, 400); sidx = np.flatnonzero(~cc.keep)
print("ballast libs", X[bal].sum(1), "sub-gate", X[np.ix_(bal, sidx)].sum(1))
s_, j_, p, q, l = E.block_stats(X, cc)
orig_sub = rp423._subgate
def sub(X_, car, bal_, sidx_, sub_delta, b_dep, L_bal):
    print("  before _subgate: dk", (X_[bal_].sum(1) - L_bal).tolist(), "carrier net sub need", int(sub_delta.sum()), "b_dep", b_dep)
    out = orig_sub(X_, car, bal_, sidx_, sub_delta, b_dep, L_bal)
    print("  after: row libs - target", (X_[bal_].sum(1) - L_bal).tolist(), "sub-gate left", X_[np.ix_(bal_, sidx_)].sum(1).tolist())
    return out
rp423._subgate = sub
try:
    rp423.postpass(X, n_car, cc, l, rp423.parse_variant("dn_k200"))
except AssertionError as e:
    print("FAIL", e)
