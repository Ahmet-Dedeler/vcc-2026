"""382 item 3: even/odd split held-out accuracy of the Jurkat gene prior on A's block.

Own code throughout. Known answers first: A shape, median library, gene-order
head, Jurkat own-target negativity. Then: A mean-CPM via chunked h5py CSR scan;
block = top-2200 by A-CPM among Jurkat-mapped genes; P(down) fit on even target
rows, held-out eval on odd rows (and swapped); eval set = block minus the
target's head (top-400 |delta|); exact-zero deltas dropped. Unweighted accuracy
for the gene prior (P>0.5 -> down) and all-down. Expect ~0.5766 / ~0.5147.
"""
import json
import h5py
import numpy as np

CHUNK = 2000
out = {"known_answers": {}, "splits": {}}

# ---- context A mean CPM, chunked ----
f = h5py.File("data/challenge/context_A.h5ad", "r")
n_cells, n_genes = [int(v) for v in f["X"].attrs["shape"]]
out["known_answers"]["shape"] = [n_cells, n_genes]
data, indices, indptr = f["X/data"], f["X/indices"], f["X/indptr"]
# var gene names: find the string dataset under var/_index
vnames = None
def hunt(g, prefix=""):
    global vnames
    for k in g.keys():
        o = g[k]
        if isinstance(o, h5py.Dataset) and o.dtype.kind in ("O", "S", "U"):
            vnames = np.array([x.decode() if isinstance(x, bytes) else str(x)
                               for x in o[:]])
            return True
        if isinstance(o, h5py.Group) and hunt(o):
            return True
    return False
hunt(f["var"])
assert vnames is not None and len(vnames) == n_genes, (vnames is None, len(vnames) if vnames is not None else None)
csv = [l.strip() for l in open("data/challenge/gene_names.csv")][1:]
assert list(vnames[:3]) == ["TSPAN6", "TNMD", "DPM1"], list(vnames[:3])
assert csv[:3] == ["TSPAN6", "TNMD", "DPM1"], csv[:3]
out["known_answers"]["gene_head_match_csv_and_h5ad"] = True

sums = np.zeros(n_genes, dtype=np.float64)
libs = np.zeros(n_cells, dtype=np.float64)
for s in range(0, n_cells, CHUNK):
    e = min(s + CHUNK, n_cells)
    ip = indptr[s:e + 1]
    d = data[ip[0]:ip[-1]]
    ix = indices[ip[0]:ip[-1]]
    lib = np.add.reduceat(d, ip - ip[0]) if False else None
    # row sums via bincount on repeat counts
    counts = np.diff(ip)
    rows = np.repeat(np.arange(e - s), counts)
    lib = np.bincount(rows, weights=d, minlength=e - s)
    libs[s:e] = lib
    sums += np.bincount(ix, weights=d / lib[rows], minlength=n_genes)
mean_cpm = sums / n_cells * 1e6
out["known_answers"]["median_library"] = float(np.median(libs))
out["known_answers"]["median_library_match_20109"] = bool(abs(float(np.median(libs)) - 20109) < 1)

# ---- Jurkat bank ----
jg = [l.strip() for l in open("data/donors/nadig_jurkat/genes.txt")]
jp = [l.strip() for l in open("data/donors/nadig_jurkat/perts.txt")]
jd = np.load("data/donors/nadig_jurkat/delta.npy", mmap_mode="r")
assert jd.shape == (len(jp), len(jg)), (jd.shape, len(jp), len(jg))
J = {g: j for j, g in enumerate(jg)}
nT = len(jp)

# own-target negativity (sign-convention check)
jp_arr = np.array(jp)
tot, neg = 0, 0.0
for t in set(jp_arr):
    if t in J:
        r = np.where(jp_arr == t)[0]
        tot += len(r)
        neg += float((np.asarray(jd[r][:, J[t]]) < 0).sum())
out["known_answers"]["own_target_rows"] = tot
out["known_answers"]["own_target_neg_frac"] = neg / tot

# ---- block: top-2200 by A-CPM among mapped genes ----
A = {g: i for i, g in enumerate(csv)}
mapped = [g for g in csv if g in J]
cpm_m = np.array([mean_cpm[A[g]] for g in mapped])
order = np.argsort(-cpm_m, kind="stable")
block = [mapped[i] for i in order[:2200]]
block_j = np.array([J[g] for g in block], dtype=int)
out["block"] = {"size": len(block), "cutoff_cpm": float(cpm_m[order[2199]]),
                "n_mapped": len(mapped)}

def run_direction(fit_idx, eval_idx, label):
    P = (np.asarray(jd[fit_idx][:, block_j]) < 0).mean(axis=0)  # P(down)
    n_down = n_tot = 0
    n_down_base = n_tot_base = 0
    n_zero = 0
    for r in eval_idx:
        full = np.asarray(jd[r])                       # full 8881 donor axis
        head = set(np.argsort(-np.abs(full), kind="stable")[:400].tolist())
        row = full[block_j]                            # block-ordered values
        is_head = np.array([j in head for j in block_j])
        nz = row != 0
        keep = (~is_head) & nz
        n_zero += int((~nz).sum())
        s = row[keep] < 0                      # True = down
        p = P[keep] > 0.5
        n_down += int((p == s).sum())
        n_tot += int(keep.sum())
        n_down_base += int(s.sum())
        n_tot_base += int(keep.sum())
    out["splits"][label] = {
        "prior_acc": n_down / n_tot, "alldown_acc": n_down_base / n_tot_base,
        "n_eval": n_tot, "n_zero_dropped": n_zero}

idx = np.arange(nT)
run_direction(idx[0::2], idx[1::2], "fit_even_eval_odd")
run_direction(idx[1::2], idx[0::2], "fit_odd_eval_even")
a = out["splits"]["fit_even_eval_odd"]
b = out["splits"]["fit_odd_eval_even"]
out["mean"] = {"prior_acc": (a["prior_acc"] + b["prior_acc"]) / 2,
               "alldown_acc": (a["alldown_acc"] + b["alldown_acc"]) / 2}
out["expected"] = {"prior_acc": 0.5766, "alldown_acc": 0.5147, "tol": 0.003}
out["match"] = {"prior_acc": bool(abs(out["mean"]["prior_acc"] - 0.5766) <= 0.003),
                "alldown_acc": bool(abs(out["mean"]["alldown_acc"] - 0.5147) <= 0.003)}

print(json.dumps(out, indent=1))
json.dump(out, open("experiments/382-night-audit/audit3.json", "w"), indent=1,
          default=float)
