"""462i — is 462h's K562/RPE1 agreement the lab or essential-gene biology? Same targets, other labs: Orion HCT116 and
HEK293T (Xaira, CRISPRi, genome-wide) and VIPerturb K562 (Satija) against the Weissman K562-essential and RPE1 programmes,
all over the targets every bank carries (Replogle essential set), on the genes all measure. Full-cell programmes (the
banks' deltas; Weissman from 462h's counts), so cross-bank r is honest (independent controls) but attenuated by each
bank's noise. usage: uv run python .../weissxaira462.py"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "scripts"))
import weissman462 as W  # noqa: E402
import def_sameline as S  # noqa: E402
S.BANK_DIRS.update({"orion_hek293t": ROOT / "data/donors/orion_hek293t"})
D = {}
for i, (line, path) in enumerate(W.FILES.items()):
    keys, genes, Sx, n = W.sums(path, 462 + i); D[line] = (keys, list(genes), Sx)
banks = {b: S.load_bank(b) for b in ("orion_hct116", "orion_hek293t", "viperturb_k562", "k562_gwps")}
T = set(D["K562"][0]) & set(D["RPE1"][0]) - {"control"}
for b, (delta, perts, genes, _) in banks.items():
    T &= set(perts)
T = sorted(T)
G = set(D["K562"][1]) & set(D["RPE1"][1])
for b, (delta, perts, genes, _) in banks.items():
    G &= set(genes)
G = sorted(G)
P = {}
for line, (keys, genes, Sx) in D.items():
    gi = {g: j for j, g in enumerate(genes)}; cols = np.array([gi[g] for g in G]); c = keys.index("control")
    C = Sx[c, 0, cols]; M = np.stack([W.lp(Sx[keys.index(t), 0, cols]) - W.lp(C) for t in T])
    P[f"weissman_{line}"] = W.pdn(M)
for b, (delta, perts, genes, _) in banks.items():
    gi = {g: j for j, g in enumerate(genes)}; cols = np.array([gi[g] for g in G])
    M = np.asarray(delta[[perts.index(t) for t in T]], dtype=np.float64)[:, cols]
    P[b] = W.pdn(M)
lab = {"weissman_K562": "W", "weissman_RPE1": "W", "k562_gwps": "W", "orion_hct116": "X", "orion_hek293t": "X", "viperturb_k562": "S"}
m = np.ones(len(G), bool)
rows = []
names = list(P)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]
        rows.append(dict(a=a, b=b, same_lab=lab[a] == lab[b], r=W.r(P[a], P[b], m)))
res = {"targets": len(T), "genes": len(G), "pairs": rows}
print(f"targets {len(T)}, genes {len(G)}")
for x in rows: print(f"  {x['a']:15s} {x['b']:15s} {'same lab ' if x['same_lab'] else 'cross-lab'} r {x['r']:+.3f}")
(HERE / "work/weissxaira462.json").write_text(json.dumps(res, indent=1))
