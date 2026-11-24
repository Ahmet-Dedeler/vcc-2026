"""517 §2: a six-line consensus WEIGHTED by each line's measured cross-line predictive quality (xline517: mean sign accuracy
predicting the other five lines; weight = accuracy - 0.5), vs CLEAN4 and the equal-weight six-line consensus, on H1 (a line
none of them is): 462f-exact accuracy rows, and the analytic block swap and hook prices (c4 base, 3 seeds)."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "scripts"))
import price517 as P
import def_sameline as S
RC = P.RC
W = {"K562": 0.0792, "RPE1": 0.0757, "Jurkat": 0.0391, "HepG2": 0.0867, "HCT116": 0.1005, "HEK293T": 0.0212}
DIRS = {"K562": "k562_gwps", "RPE1": "rpe1", "Jurkat": "nadig_jurkat", "HepG2": None, "HCT116": "orion_hct116", "HEK293T": "orion_hek293t"}
S.BANK_DIRS.update({"rpe1": ROOT / "data/donors/rpe1", "orion_hct116": ROOT / "data/donors/orion_hct116",
                    "orion_hek293t": ROOT / "data/donors/orion_hek293t", "nadig_hepg2": ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2"})
DIRS["HepG2"] = "nadig_hepg2"
zr = np.load(RC.REF, allow_pickle=False); genes_full = zr["genes"]; keep = np.asarray(zr["keep"], bool); T = zr["targets"].astype(str)
lfc = np.asarray(zr["lfc"], float); q = np.asarray(zr["q"], float); bctrl = np.asarray(zr["b_ctrl"], float)
kg = np.asarray(genes_full)[keep].astype(str)
Z = {}
for L, b in DIRS.items():
    delta, perts, genes, _ = S.load_bank(b)
    pdn = S.programme(delta, perts, genes, list(kg), set(T))[1]
    Z[L] = np.where(np.isfinite(pdn), (0.5 - pdn) / np.nanstd(pdn), np.nan)
def wmean(ws):
    num = np.zeros(kg.size); den = np.zeros(kg.size)
    for L, w in ws.items():
        m = np.isfinite(Z[L]); num[m] += w * Z[L][m]; den[m] += w
    return np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)
SRC, _ = P.load_sources(genes_full, keep)
srcs = {"CLEAN4": 0.5 - SRC["c4"], "six-line equal": wmean({k: 1.0 for k in W}), "six-line quality-weighted": wmean(W)}
e_ = np.expm1(bctrl[keep]); cpm = 1e6 * e_ / e_.sum(); order = np.argsort(-cpm)
block = np.zeros(kg.size, bool); block[[i for i in order if kg[i] not in set(T)][:2200]] = True
for name, z in srcs.items():
    s = np.where(np.isfinite(z) & (z != 0), np.sign(z), np.nan); hit = {"fidlike": [0, 0], "gate_block": [0, 0], "gate_nonblock": [0, 0]}
    for t in range(len(T)):
        base = (kg != T[t]) & np.isfinite(lfc[t]) & (lfc[t] != 0) & np.isfinite(s)
        for k, m in (("fidlike", base & block), ("gate_block", base & block & (q[t] < .05)), ("gate_nonblock", base & ~block & (q[t] < .05))):
            hit[k][0] += int((np.sign(lfc[t][m]) == s[m]).sum()); hit[k][1] += int(m.sum())
    print(f"{name:26s} " + "  ".join(f"{k} {a / max(n, 1):.4f} (n {n})" for k, (a, n) in hit.items()), flush=True)
c4 = RC.prior_signs_full(genes_full)[keep]
n_c4 = int((np.isfinite(SRC["c4"]) & (np.abs(SRC["c4"] - .5) >= .05)).sum())
res = {}
for seed in (0, 1, 2):
    ctx = P.Ctx("c4", seed, lambda i, t: c4); rng = np.random.default_rng(5182 + seed)
    for name, z in srcs.items():
        if name == "CLEAN4":
            continue
        blk = np.where(np.isfinite(z) & (z != 0), np.sign(z), -1.0)
        a = np.abs(np.nan_to_num(z)); hs = np.where(a >= np.sort(a)[::-1][n_c4 - 1], np.sign(np.nan_to_num(z)), 0.0)
        for arm, fn in (("block", lambda: ctx.block(lambda i, t: blk)), ("hook@c4cov", lambda: ctx.hook(lambda i, t: hs))):
            lfc_new, _ = fn(); sc, _ = ctx.score(lfc_new, rng, n_boot=300)
            res.setdefault(f"{name} {arm}", []).append((sc["nmae"], sc["fid"], sc["reach"], sc["overall_equiv"]))
    del ctx
for k, v in res.items():
    a = np.array(v).mean(0)
    print(f"{k:38s} nmae {a[0]:+.5f} fid {a[1]:+.5f} reach {a[2]:+.5f} ov {a[3]:+.5f}")
