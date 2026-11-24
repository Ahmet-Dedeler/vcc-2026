"""291 recon 2 — Broad-16 target availability COUNTS ONLY (no feature-outcome number)."""
import sys, json
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path("~/Code/virtual-cell-challenge-2026")
sys.path.insert(0, str(ROOT/'experiments/268-broad16-context-residuals/work'))
from load268 import *          # noqa
z = np.load(ROOT/'experiments/290-remainder-response-features/work/persrc290.npz', allow_pickle=True)
held = set(str(g) for g in z['genes'])
H, genes = load_halves(WORKING)
print('lines', len(WORKING), WORKING)
tg = ALL_TARGETS; ctrl = set(CONTROLS)
R, ok, avail = build_responses(H, genes, WORKING, tg, keys=("full",), min_cpm=5.0)
print('genes kept', int(ok.sum()))
allok = avail.all(1)
avail_t = [t for t,a in zip(tg, allok) if a]
print('targets available on all 12 working lines:', len(avail_t))
nonor = [t for t in avail_t if t not in ctrl]
print('  non-OR:', len(nonor), ' in held map:', len([t for t in nonor if t in held]))
print('  OR controls available:', len([t for t in avail_t if t in ctrl]),
      'in held map:', len([t for t in avail_t if t in ctrl and t in held]))
missing = [t for t in tg if t not in avail_t]
print('  unavailable:', missing)
print('  non-OR available but NOT in held map:', [t for t in nonor if t not in held])
# cells per KO
nc = {}
for t in avail_t:
    nc[t] = [int(H[l]['ncell'][H[l]['targets'].index(t),0]) for l in WORKING]
print('cells/KO median over lines: min', min(np.median(v) for v in nc.values()),
      'max', max(np.median(v) for v in nc.values()))
cq = pd.read_csv(ROOT/'data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv')
print(cq[cq['cell_line'].isin(WORKING)][['cell_line','avg_umis_per_ko','n_cells_per_ko']].to_string(index=False))
