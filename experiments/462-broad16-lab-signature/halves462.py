"""462b — is the Broad lines' high programme reliability shared-control noise? Cell halves with DISJOINT controls.

For each working line, 268's seeded cell halves: half A of every target against half A of the 10 OR-KO controls, half B
against half B. P_down(A) and P_down(B) then share no cell on either side. Within-line r(P_A, P_B) is the honest
reliability (Spearman-Brown to full cells). A of line i against B of line j is the cross-line agreement with no shared
noise. As a control: the SHARED-control version, with both target halves measured against the full control, shows how
much a common control inflates the within-line number.
usage: uv run python experiments/462-broad16-lab-signature/halves462.py   -> work/halves462.json
"""
import itertools, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/268-broad16-context-residuals/work")); sys.path.insert(0, str(HERE))
import load268 as B  # noqa: E402
from lab462 import pdown  # noqa: E402

lines = list(B.WORKING)
H, genes = B.load_halves(lines)
T = [t for t in B.usable_targets(H, lines, min_cells=6) if t not in B.CONTROLS]
R, ok, avail = B.build_responses(H, genes, lines, T, keys=("full", "A", "B"), min_cpm=5)
# shared-control variant: target half A and half B each against the FULL control
Rs = {}
for k in ("A", "B"):
    M = np.zeros_like(R[k]); tl0 = None
    for j, l in enumerate(lines):
        tl = H[l]["targets"]; ci = [tl.index(c) for c in B.CONTROLS]
        ctrl = B.logprof(np.nansum(H[l]["full"][ci], 0), ok)
        for i, t in enumerate(T):
            if not R["avail"][k][i, j]: continue
            M[i, j] = B.logprof(H[l][k][tl.index(t)], ok) - ctrl
    Rs[k] = M
P = {}
for k in ("full", "A", "B"):
    X = R[k].astype(np.float64); X[~R["avail"][k]] = np.nan
    P[k] = [pdown(X[:, j, :])[0] for j in range(len(lines))]
for k in ("A", "B"):
    X = Rs[k].astype(np.float64); X[~R["avail"][k]] = np.nan
    P["shared" + k] = [pdown(X[:, j, :])[0] for j in range(len(lines))]
def r(a, b):
    m = np.isfinite(a) & np.isfinite(b); return float(np.corrcoef(a[m] - .5, b[m] - .5)[0, 1])
within_disjoint = [r(P["A"][j], P["B"][j]) for j in range(len(lines))]
within_shared = [r(P["sharedA"][j], P["sharedB"][j]) for j in range(len(lines))]
cross_disjoint = [r(P["A"][i], P["B"][j]) for i, j in itertools.permutations(range(len(lines)), 2)]
cross_full = [r(P["full"][i], P["full"][j]) for i, j in itertools.combinations(range(len(lines)), 2)]
sb = lambda x: 2 * x / (1 + x)
res = dict(lines=lines, n_targets=len(T), n_genes=int(ok.sum()),
           within_line_disjoint_controls=dict(zip(lines, within_disjoint)),
           within_line_shared_control=dict(zip(lines, within_shared)),
           mean_within_disjoint=float(np.mean(within_disjoint)), mean_within_shared=float(np.mean(within_shared)),
           mean_within_disjoint_SB=float(np.mean([sb(x) for x in within_disjoint])),
           mean_cross_line_disjoint=float(np.mean(cross_disjoint)), mean_cross_line_full=float(np.mean(cross_full)),
           cross_over_within_disjoint=float(np.mean(cross_disjoint) / np.mean(within_disjoint)))
(HERE / "work/halves462.json").write_text(json.dumps(res, indent=1))
print(json.dumps({k: v for k, v in res.items() if not isinstance(v, dict)}, indent=1))
print("within, disjoint controls:", {l: round(v, 3) for l, v in zip(lines, within_disjoint)})
print("within, shared control:   ", {l: round(v, 3) for l, v in zip(lines, within_shared)})
