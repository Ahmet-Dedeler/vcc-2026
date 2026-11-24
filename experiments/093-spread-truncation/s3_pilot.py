"""Step 3: pilot. Re-emit a handful of perturbations across a global-spread ladder and
measure (a) the wall cost, (b) the exactness of the panel totals, (c) how far `spread`
actually moves `n_pred`.

`spread` is clipped to [0, 1] inside `emit_counts`, so 1.0 -- pure Poisson scatter with the
gene total repaired -- is the hardest de-signification this channel can apply. If that does
not reach the 850 the counterfactual wants, no per-gene policy built on it can either, and
that is the answer to the open question.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib093 as L  # noqa: E402

gate = L.Gate()
tot = np.load(os.path.join(L.SCRATCH, "totals093.npy"))
d = np.load(os.path.join(L.E092, "de_tables.npz"), allow_pickle=True)
perts = d["perts"].astype(str)
pq_ref = d["pred_padj"]
rq = d["real_padj"]
G = gate.G
gpos = {g: i for i, g in enumerate(gate.genes)}

sel = np.arange(0, len(perts), 11)[:12]
ladder = [0.25, 0.4, 0.55, 0.7, 0.85, 1.0]
print(f"pilot on {len(sel)} perturbations, global spread ladder {ladder}")
gate.ctrl()
print("control block ready", flush=True)

rows = {s: [] for s in ladder}
ref_n = []
n_real = []
t0 = time.time()
for i in sel:
    ok = np.ones(G, bool)
    j = gpos.get(perts[i])
    if j is not None:
        ok[j] = False
    ref_n.append(int((ok & (pq_ref[i] < L.ALPHA)).sum()))
    n_real.append(int((ok & (rq[i] < L.ALPHA)).sum()))
    for s in ladder:
        rng = np.random.default_rng(1000 + i)
        te = time.time()
        block = L.emit_block(tot[i], np.full(gate.G_full, s), rng)
        assert np.array_equal(block.sum(axis=0).astype(np.int64), tot[i]), "totals moved!"
        td = time.time()
        pval, padj, lfc = gate.de(block)
        rows[s].append(int((ok & (padj < L.ALPHA)).sum()))
        if i == sel[0]:
            print(f"   spread {s}: emit {td - te:.2f}s  de {time.time() - td:.2f}s", flush=True)
        del block
    print(f"  {perts[i]:>10}  n_real {n_real[-1]:>5}  prod {ref_n[-1]:>5}  "
          + "  ".join(f"s{s}:{rows[s][-1]:>5}" for s in ladder), flush=True)

print(f"\ntotal {time.time() - t0:.0f}s for {len(sel)} perts x {len(ladder)} spreads "
      f"=> {(time.time() - t0) / (len(sel) * len(ladder)):.2f}s per (pert, policy)")
print(f"\n{'':>12}{'median n_pred':>15}{'min':>8}{'max':>8}")
print(f"{'n_real':>12}{np.median(n_real):>15.0f}{np.min(n_real):>8}{np.max(n_real):>8}")
print(f"{'production':>12}{np.median(ref_n):>15.0f}{np.min(ref_n):>8}{np.max(ref_n):>8}")
for s in ladder:
    v = np.array(rows[s])
    print(f"{f'spread {s}':>12}{np.median(v):>15.0f}{v.min():>8}{v.max():>8}")
