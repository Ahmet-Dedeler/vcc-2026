"""462c — how much of Arc's H1 programme (the prior 458 ships) is its control's noise? 462b's test on the H1 cells.

data/harness/h1-full.h5ad = the h1_2025 bank's source: 255 perturbations x 400 cells + 18,399 NTC cells (Flex). Seeded
random halves of every perturbation's cells and of the NTC cells.
  shared:   P_down(T1 vs full C) against P_down(T2 vs full C)   -> signal + control noise + half-cell noise
  disjoint: P_down(T1 vs C1)     against P_down(T2 vs C2)       -> signal only (no shared noise)
Also: the full programme (T vs C, what h1_2025's delta gives) against the bank's own P_down as a known answer, and each
half-programme's agreement with the A/B/C blocks' CLEAN4 prior.
Response = log1p(5e4 p / sum p) of the summed counts (the bank's frame), panel rows excluded as 458 does (A/B/C panel).
usage: uv run python experiments/462-broad16-lab-signature/h1halves462.py   -> work/h1halves462.json
"""
import csv, json, sys
from pathlib import Path
import h5py, numpy as np, scipy.sparse as sp
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts")); import def_sameline as S  # noqa: E402
S.BANK_DIRS.setdefault("h1_2025", ROOT / "data/donors/h1_2025")
TS = 5e4; rng = np.random.default_rng(462)
with h5py.File(ROOT / "data/harness/h1-full.h5ad", "r") as f:
    X = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=tuple(f["X"].attrs["shape"]))
    tg = f["obs/target_gene"]; tv = tg["categories"][:].astype(str)[tg["codes"][:]]
    vn = f["var/_index"]; genes = (vn[:] if not isinstance(vn, h5py.Group) else vn["values"][:]).astype(str)
panel = {r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))}
ctrl = np.where(tv == "non-targeting")[0]; perts = [p for p in np.unique(tv) if p != "non-targeting" and p not in panel]
print(f"cells {X.shape}, controls {ctrl.size}, non-panel perturbations {len(perts)}", flush=True)
def psum(idx): return np.asarray(X[idx].sum(0), dtype=np.float64).ravel()
def lp(c): return np.log1p(TS * c / c.sum())
cp = rng.permutation(ctrl); C1, C2 = cp[: ctrl.size // 2], cp[ctrl.size // 2:]
cF, c1, c2 = lp(psum(ctrl)), lp(psum(C1)), lp(psum(C2))
D = {k: np.zeros((len(perts), X.shape[1])) for k in ("full", "T1C", "T2C", "T1C1", "T2C2")}
for i, p in enumerate(perts):
    idx = rng.permutation(np.where(tv == p)[0]); h = idx.size // 2
    tF, t1, t2 = lp(psum(idx)), lp(psum(idx[:h])), lp(psum(idx[h:]))
    D["full"][i] = tF - cF; D["T1C"][i] = t1 - cF; D["T2C"][i] = t2 - cF; D["T1C1"][i] = t1 - c1; D["T2C2"][i] = t2 - c2
def pdn(M):
    nz = (M != 0).sum(0); dn = (M < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"): return np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
P = {k: pdn(v) for k, v in D.items()}
cpm = 1e6 * np.expm1(cF) / np.expm1(cF).sum(); expressed = cpm >= 5
def r(a, b, m):
    m = m & np.isfinite(a) & np.isfinite(b); return float(np.corrcoef(a[m] - .5, b[m] - .5)[0, 1])
res = {"n_perts": len(perts), "n_ctrl": int(ctrl.size), "n_expressed_genes": int(expressed.sum())}
# known answer: the full programme vs the bank's own P_down (same cells, same frame)
official = [x["gene_name"] for x in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
d, bp, bg, _ = S.load_bank("h1_2025"); _, bank_pd, nrows, _ = S.programme(d, bp, bg, official, panel)
pos = {g: i for i, g in enumerate(genes)}; col = np.array([pos.get(g, -1) for g in official])
mine = np.where(col >= 0, P["full"][np.maximum(col, 0)], np.nan)
m = np.isfinite(mine) & np.isfinite(bank_pd)
res["ka_full_vs_bank"] = dict(bank_rows=nrows, r=float(np.corrcoef(mine[m], bank_pd[m])[0, 1]),
                              max_abs=float(np.nanmax(np.abs(mine[m] - bank_pd[m]))))
for tag, mask in (("expressed_cpm5", expressed), ("all_measured", np.ones_like(expressed))):
    rs, rd = r(P["T1C"], P["T2C"], mask), r(P["T1C1"], P["T2C2"], mask)
    sb = lambda x: 2 * x / (1 + x)
    res[tag] = dict(shared_control_r=rs, disjoint_control_r=rd, shared_SB=sb(rs), disjoint_SB=sb(rd),
                    control_noise_share_of_reliability=float(1 - rd / rs) if rs > 0 else None)
# on the A/B/C blocks (2,200 top-CPM non-panel genes of each context's controls), in the official axis
blocks = {}
for c in "ABC":
    b = np.load(ROOT / f"experiments/438-october-stack-repricing/work/abc_frames/bctrl-{c}.npy"); e = np.expm1(b)
    order = np.argsort(-(e / e.sum())); blk = np.asarray([i for i in order if official[i] not in panel][:2200])
    to_h1 = col[blk]; ok = to_h1 >= 0
    bm = np.zeros(len(genes), bool); bm[to_h1[ok]] = True
    blocks[c] = dict(n=int(ok.sum()), shared_control_r=r(P["T1C"], P["T2C"], bm), disjoint_control_r=r(P["T1C1"], P["T2C2"], bm))
res["abc_blocks"] = blocks
print(json.dumps(res, indent=1)); (HERE / "work/h1halves462.json").write_text(json.dumps(res, indent=1))
