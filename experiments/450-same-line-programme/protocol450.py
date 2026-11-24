"""450i (2026-09-26): does a PROTOCOL-matched programme beat CLEAN4 for a Flex recipient? labline450 showed per-gene P_down
agrees within a lab (r ~0.5 across lines) and not across labs (r ~0.09 even within K562). The VCC H1 truth is 10x Flex;
VIPerturb-seq (Satija) is a Flex-style probe-based screen in K562. On h1-test (sameline450's truth, gate and block),
score gate-block sign accuracy of each bank's programme (P_down over its rows with the h1-test targets EXCLUDED):
CLEAN4 (Weissman 3'), Pan H1 (Whitehead 5', same line), VIPerturb K562 (Satija Flex-style, other line), Orion HCT116 /
HEK293T (Xaira), and the Weissman lines one by one.   usage: uv run python experiments/450-same-line-programme/protocol450.py"""
import runpy
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.argv = [sys.argv[0]]
g = runpy.run_path(str(HERE / "sameline450.py"), run_name="not_main")
P, gate, block, truth, T, G, adj = (g[k] for k in ("P", "gate", "block", "truth", "T", "G", "adj"))
sys.path.insert(0, str(ROOT / "scripts"))
import def_sameline as S  # noqa: E402
S.BANK_DIRS.setdefault("orion_hek293t", ROOT / "data/donors/orion_hek293t")
S.BANK_DIRS.setdefault("rpe1", ROOT / "data/donors/rpe1")
h1_targets = set(map(str, T))


def score(sign):
    out = {}
    for name, m in (("gate_block", gate & block), ("gate_all", gate), ("block_adj", adj & block)):
        mm = m & np.isfinite(sign)[None, :]
        out[name] = (float(np.mean((np.broadcast_to(sign, truth.shape) == truth)[mm])), int(mm.sum()))
    return out


res = {}
for name in ("clean4", "pan_pdown"):
    res[name] = score(P[name])
for b in ("viperturb_k562", "k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2", "rpe1", "orion_hct116", "orion_hek293t"):
    delta, perts, genes, _ = S.load_bank(b)
    psign, pdown, nrows, _ = S.programme(delta, perts, genes, list(map(str, G)), h1_targets)
    sign = np.where(np.isfinite(pdown), np.where(pdown > 0.5, -1.0, 1.0), np.nan)
    res[b] = score(sign)
print(f"{'predictor':16s} {'gate-block':>18s} {'gate-all':>18s} {'block (fid-like)':>18s}")
for k, v in res.items():
    print(f"{k:16s} " + " ".join(f"{v[m][0]:9.4f} n{v[m][1]:7d}" for m in ("gate_block", "gate_all", "block_adj")))
# common-coverage comparison: genes every predictor covers
cov = np.ones(len(G), bool)
signs = {}
for k in ("clean4", "pan_pdown"):
    signs[k] = P[k]
for b in ("viperturb_k562", "orion_hct116"):
    delta, perts, genes, _ = S.load_bank(b)
    _, pdown, _, _ = S.programme(delta, perts, genes, list(map(str, G)), h1_targets)
    signs[b] = np.where(np.isfinite(pdown), np.where(pdown > 0.5, -1.0, 1.0), np.nan)
for s in signs.values():
    cov &= np.isfinite(s)
print(f"\non the {int(cov.sum())} genes all of clean4 / pan / viperturb / hct116 cover:")
for k, s in signs.items():
    v = score(np.where(cov, s, np.nan))
    print(f"  {k:16s} gate-block {v['gate_block'][0]:.4f} (n {v['gate_block'][1]})  gate-all {v['gate_all'][0]:.4f}")
