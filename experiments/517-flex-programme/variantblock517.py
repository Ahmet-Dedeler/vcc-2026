"""517: the Arc estimator variants as the BLOCK prior on H1 (analytic block path, validated in blockval517): gated (467's
strong 57) and ungated, vs the c4 base; pdown_all == 462l/467's construction (KA: equals grid517's r0 w1 q0.8 rows)."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import price517 as P
RC = P.RC
ESS = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}
zr = np.load(RC.REF, allow_pickle=False); genes_full = zr["genes"]; keep = np.asarray(zr["keep"], bool); T = zr["targets"].astype(str)
d = ROOT / "data/donors/h1_2025_full"; delta = np.load(d / "delta.npy"); perts = [l.strip() for l in open(d / "perts.txt")]
pos = {g: j for j, g in enumerate(l.strip() for l in open(d / "genes.txt"))}
col = np.array([pos[str(g)] for g in np.asarray(genes_full)[keep]])
rows = np.array([i for i, p in enumerate(perts) if p not in set(T)]); ne = np.array([perts[i] not in ESS for i in rows])
X = np.asarray(delta[rows][:, col], dtype=np.float64); sd = X.std(0); sd[sd == 0] = 1
def pdown(M):
    nz = (M != 0).sum(0); dn = (M < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
Z = {"pdown_all": 0.5 - pdown(X), "meanz_all": (X / sd).mean(0), "meanz_ne": (X[ne] / sd).mean(0), "pdown_ne": 0.5 - pdown(X[ne])}
c4 = RC.prior_signs_full(genes_full)[keep]; strong, _ = P.strong_set(T)
out = {}
for seed in (0, 1, 2):
    ctx = P.Ctx("c4", seed, lambda i, t: c4); rng = np.random.default_rng(5177 + seed)
    for name, z in Z.items():
        s = np.where(np.isfinite(z) & (z != 0), np.sign(z), -1.0)
        for gated in (True, False):
            pr = (lambda i, t, s=s: s if strong[i] else c4) if gated else (lambda i, t, s=s: s)
            lfc, _ = ctx.block(pr); sc, _ = ctx.score(lfc, rng, n_boot=300)
            out.setdefault(f"{name}_{'g' if gated else 'u'}", []).append({**{m: sc[m] for m in P.SP}, "ov": sc["overall_equiv"], "ne_ov": sc["ne_overall_equiv"]})
for k, v in out.items():
    m = lambda f: np.mean([r[f] for r in v])
    print(f"{k:12s} fid {m('fid'):+.5f} nmae {m('nmae'):+.5f} reach {m('reach'):+.5f} | ov {m('ov'):+.5f} ne {m('ne_ov'):+.5f}")
(HERE / "work/variantblock517.json").write_text(json.dumps(out, indent=1))
