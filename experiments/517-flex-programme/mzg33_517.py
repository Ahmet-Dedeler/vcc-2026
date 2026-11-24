"""517: the meanz_ne hook (variant517's construction) on the g33 (458g) and arc (458d) bases, analytic, 3 seeds."""
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
z = (X[ne] / sd).mean(0)
SRC, _ = P.load_sources(genes_full, keep)
n_cov = int((np.isfinite(SRC["arc"]) & (np.abs(SRC["arc"] - 0.5) >= 0.05)).sum())
a = np.abs(z); hs = np.where(a >= np.sort(a)[::-1][n_cov - 1], np.sign(z), 0.0)
hc = np.where(np.isfinite(SRC["c4"]) & (np.abs(SRC["c4"] - 0.5) >= 0.05), np.where(SRC["c4"] > 0.5, -1.0, 1.0), 0.0)
c4 = RC.prior_signs_full(genes_full)[keep]; arc = P.sign_from_P(SRC["arc"]); strong, _ = P.strong_set(T)
bases = {"g33": lambda i, t: arc if strong[i] else c4, "arc": lambda i, t: arc}
out = {}
for b, pr in bases.items():
    for seed in (0, 1, 2):
        ctx = P.Ctx(b, seed, pr); rng = np.random.default_rng(5178 + seed)
        for gated in (True, False):
            fn = (lambda i, t: hs if strong[i] else hc) if gated else (lambda i, t: hs)
            lfc, _ = ctx.hook(fn); sc, _ = ctx.score(lfc, rng, n_boot=500)
            out.setdefault(f"{b}_{'g' if gated else 'u'}", []).append({**{m: sc[m] for m in P.SP}, "ov": sc["overall_equiv"], "ci95": sc["ci95"], "ne_nmae": sc["ne_nmae"]})
        del ctx
for k, v in out.items():
    m = lambda f: np.mean([r[f] for r in v])
    print(f"{k:6s} meanz_ne hook: nmae {m('nmae'):+.5f} fid {m('fid'):+.5f} reach {m('reach'):+.5f} ov {m('ov'):+.5f} [{v[0]['ci95'][0]:+.4f},{v[0]['ci95'][1]:+.4f}] ne_nmae {m('ne_nmae'):+.5f}")
(HERE / "work/mzg33_517.json").write_text(json.dumps(out, indent=1))
