"""517: meanz_ne's attenuation curve (grid517's random-flip emulation, the same seed; flips invert z) for the hook alone on
the CLEAN4 block (gated) and for block + hook from one meanz_ne artifact (VCC_SAMELINE=$MZ with programme:0.05), gated and
ungated; analytic, 3 seeds. Feeds predict517's meanz_ne rows."""
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
z0 = (X[ne] / sd).mean(0)
SRC, _ = P.load_sources(genes_full, keep)
n_cov = int((np.isfinite(SRC["arc"]) & (np.abs(SRC["arc"] - 0.5) >= 0.05)).sum())
hc = np.where(np.isfinite(SRC["c4"]) & (np.abs(SRC["c4"] - 0.5) >= 0.05), np.where(SRC["c4"] > 0.5, -1.0, 1.0), 0.0)
c4 = RC.prior_signs_full(genes_full)[keep]; strong, _ = P.strong_set(T)
u = np.random.default_rng(51717).random(keep.sum())
rows_out = []
for seed in (0, 1, 2):
    ctx = P.Ctx("c4", seed, lambda i, t: c4); rng = np.random.default_rng(5179 + seed)
    for r in (0.0, 0.1, 0.2, 0.3, 0.4):
        z = np.where(u < r, -z0, z0)
        a = np.abs(z); hs = np.where(a >= np.sort(a)[::-1][n_cov - 1], np.sign(z), 0.0)
        blk = np.where(z != 0, np.sign(z), -1.0)
        for arm, kind in (("hook_g", "hook"), ("block_u+hook_u", "both_u"), ("block_g+hook_g", "both_g")):
            if kind == "hook":
                lfc, _ = ctx.hook(lambda i, t: hs if strong[i] else hc)
            elif kind == "both_u":
                lfc, _ = ctx.block(lambda i, t: blk, hook_for_target=lambda i, t: hs)
            else:
                lfc, _ = ctx.block(lambda i, t: blk if strong[i] else c4, hook_for_target=lambda i, t: hs if strong[i] else hc)
            sc, _ = ctx.score(lfc, rng, n_boot=300)
            rows_out.append({"seed": seed, "r": r, "arm": arm, **{m: sc[m] for m in P.SP}, "ov": sc["overall_equiv"], "ne_ov": sc["ne_overall_equiv"]})
    del ctx
from collections import defaultdict
agg = defaultdict(list)
for x in rows_out: agg[(x["arm"], x["r"])].append(x)
for (arm, r), v in sorted(agg.items()):
    m = lambda f: np.mean([x[f] for x in v])
    print(f"{arm:16s} r {r:.1f} | fid {m('fid'):+.5f} nmae {m('nmae'):+.5f} reach {m('reach'):+.5f} | ov {m('ov'):+.5f} ne {m('ne_ov'):+.5f}")
(HERE / "work/mzgrid517.json").write_text(json.dumps(rows_out, indent=1))
