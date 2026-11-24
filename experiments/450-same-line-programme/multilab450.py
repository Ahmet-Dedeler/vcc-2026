"""450j: does a MULTI-LAB programme consensus beat CLEAN4 (all Weissman) on the Flex H1 truth? Consensus = mean of per-lab
z(P_down - 0.5) with the lab as the unit (Weissman = CLEAN4; Xaira = mean of HCT116, HEK293T; Satija = VIPerturb),
optionally weighted. Scored on h1-test gate-block genes (h1-test targets excluded from every bank's rows)."""
import runpy
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.argv = [sys.argv[0]]
g = runpy.run_path(str(HERE / "sameline450.py"), run_name="not_main")
gate, block, truth, T, G, c4p = (g[k] for k in ("gate", "block", "truth", "T", "G", "c4p"))
sys.path.insert(0, str(ROOT / "scripts"))
import def_sameline as S  # noqa: E402
S.BANK_DIRS.setdefault("orion_hek293t", ROOT / "data/donors/orion_hek293t")
h1 = set(map(str, T))


def pd(b):
    delta, perts, genes, _ = S.load_bank(b)
    return S.programme(delta, perts, genes, list(map(str, G)), h1)[1]


def z(v):
    m = np.isfinite(v)
    out = np.full_like(v, np.nan)
    out[m] = (v[m] - 0.5) / np.nanstd(v[m] - 0.5)
    return out


def acc(score):
    s = np.where(np.isfinite(score), np.where(score > 0, -1.0, 1.0), np.nan)
    m = gate & block & np.isfinite(s)[None, :]
    return float(np.mean((np.broadcast_to(s, truth.shape) == truth)[m])), int(m.sum())


W, X1, X2, V = z(c4p), z(pd("orion_hct116")), z(pd("orion_hek293t")), z(pd("viperturb_k562"))
Xa = np.nanmean(np.vstack([X1, X2]), axis=0)
both = np.isfinite(W) & np.isfinite(Xa) & np.isfinite(V)
print(f"genes with all three labs: {both.sum()}")
for name, sc in (("CLEAN4 (Weissman)", W), ("Xaira (HCT116+HEK293T)", Xa), ("Satija (VIPerturb)", V),
                 ("Weissman+Xaira", np.nanmean(np.vstack([W, Xa]), 0)), ("Weissman+Satija", np.nanmean(np.vstack([W, V]), 0)),
                 ("all three labs", np.nanmean(np.vstack([W, Xa, V]), 0))):
    a_all, n_all = acc(np.where(np.isfinite(W), sc, np.nan))
    a_c, n_c = acc(np.where(both, sc, np.nan))
    print(f"  {name:24s} on CLEAN4's genes {a_all:.4f} (n {n_all})   on the 3-lab common genes {a_c:.4f} (n {n_c})")
