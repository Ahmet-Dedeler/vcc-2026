"""370 — Jurkat block-sign prior. See PLAN.md. Writes work/metrics370.json,
work/per-split370.csv, work/jurkat_block_prior.npz. Read-only on all inputs."""
import os, json, csv
import numpy as np

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

ROOT = "~/Code/virtual-cell-challenge-2026"
EXP = os.path.join(ROOT, "experiments/370-jurkat-block-sign-prior")
WORK = os.path.join(EXP, "work")
os.makedirs(WORK, exist_ok=True)

def read_lines(p):
    with open(p) as f:
        return [l.strip() for l in f if l.strip() != ""]

# ---------------- known-answer reproduction (own code) ----------------
print("== known answers ==", flush=True)
import h5py
fa = h5py.File(os.path.join(ROOT, "data/challenge/context_A.h5ad"), "r")
shape = tuple(fa["X"].attrs["shape"])
print("context_A X shape:", shape)
assert shape == (18400, 18533), shape
print("X datasets:", list(fa["X"].keys()))

chal_genes = read_lines(os.path.join(ROOT, "data/challenge/gene_names.csv"))
assert chal_genes[0] == "gene_name", chal_genes[0]
chal_genes = chal_genes[1:]
assert len(chal_genes) == 18533, len(chal_genes)
assert chal_genes[:3] == ["TSPAN6", "TNMD", "DPM1"], chal_genes[:3]
import anndata
ab = anndata.read_h5ad(os.path.join(ROOT, "data/challenge/context_A.h5ad"), backed="r")
assert ab.var_names[:3].tolist() == chal_genes[:3]
assert ab.shape == (18400, 18533)
print("gene order OK:", chal_genes[:3])

# 323 known-answer: mean of per-target nmae == metrics-cap15.json nmae
m15 = json.load(open(os.path.join(ROOT,
    "experiments/323-current-family-validation/work/metrics-cap15.json")))
per = {}
with open(os.path.join(ROOT,
    "experiments/323-current-family-validation/work/per-target-cap15.csv")) as f:
    for row in csv.DictReader(f):
        per.setdefault(row["metric"], []).append(float(row["value"]))
for k in ["de_wilcoxon_lfc_nmae", "de_wilcoxon_direction_fidelity_yield_raw",
           "de_wilcoxon_direction_reach_raw"]:
    v = np.array(per[k])
    print(f"per-target {k}: n={len(v)} nan={int(np.isnan(v).sum())} "
          f"nanmean={np.nanmean(v):.8f}")
print("recorded: nmae", m15["raw"]["nmae"], "fid", m15["raw"]["fid"],
      "reach", m15["raw"]["reach"])
assert abs(np.nanmean(np.array(per["de_wilcoxon_lfc_nmae"]))
           - m15["raw"]["nmae"]) < 1e-9
print("323 nmae known-answer MATCH (diff < 1e-9)")

# ---------------- step 1: mean CPM for context A (chunked CSR) ----------------
print("== context_A mean CPM ==", flush=True)
X = fa["X"]
indptr = X["indptr"][:]
indices = X["indices"]
data = X["data"]
n_cells, n_genes = shape
assert indptr.shape == (n_cells + 1,)
gene_sums = np.zeros(n_genes, dtype=np.float64)
lib = np.zeros(n_cells, dtype=np.float64)
CH = 2000
for s in range(0, n_cells, CH):
    e = min(s + CH, n_cells)
    ip0, ip1 = int(indptr[s]), int(indptr[e])
    d = data[ip0:ip1].astype(np.float64)
    ix = indices[ip0:ip1].astype(np.int64)
    ipl = indptr[s:e + 1].astype(np.int64) - ip0
    rowlen = np.diff(ipl)
    rep = np.repeat(np.arange(e - s), rowlen)
    lib[s:e] = np.bincount(rep, weights=d, minlength=e - s)
    w = d / lib[s:e][rep]
    gene_sums += np.bincount(ix, weights=w, minlength=n_genes) * 1e6
    if s == 0:
        print(f"chunk rows {s}:{e} data {ip0}:{ip1} peak-ish MB: "
              f"{(d.nbytes + ix.nbytes) / 1e6:.1f}", flush=True)
mean_cpm = gene_sums / n_cells
med_lib = float(np.median(lib))
print(f"median library size: {med_lib:.0f} (record ~20109)")
print(f"mean CPM: median {np.median(mean_cpm):.2f}, max {mean_cpm.max():.0f} "
      f"({chal_genes[int(np.argmax(mean_cpm))]}), frac>5: "
      f"{np.mean(mean_cpm > 5):.4f}")
np.save(os.path.join(WORK, "contextA_mean_cpm.npy"), mean_cpm)
np.save(os.path.join(WORK, "contextA_lib.npy"), lib)

# ---------------- donor tables ----------------
def load_donor(name):
    d = os.path.join(ROOT, "data/donors", name)
    D = np.load(os.path.join(d, "delta.npy"), mmap_mode="r")
    g = read_lines(os.path.join(d, "genes.txt"))
    p = read_lines(os.path.join(d, "perts.txt"))
    b = np.load(os.path.join(d, "b_ctrl.npy"), mmap_mode="r")
    assert D.shape == (len(p), len(g)), (name, D.shape, len(p), len(g))
    assert b.shape == (len(g),), (name, b.shape)
    return D, g, p, np.asarray(b)

Dj, gj, pj, bj = load_donor("nadig_jurkat")
print("jurkat:", Dj.shape, "frac<0:", float(np.mean(np.asarray(Dj) < 0)))

# own-target sign check (known answer: ~99% negative)
gidx = {g: i for i, g in enumerate(gj)}
n_cov = neg = 0
Djm = np.asarray(Dj)  # 80 MB, fine
for ti, t in enumerate(pj):
    j = gidx.get(t)
    if j is not None:
        n_cov += 1
        if Djm[ti, j] < 0:
            neg += 1
print(f"jurkat own-target coverage {n_cov}/{len(pj)}, frac<0 {neg / n_cov:.5f}")
assert n_cov > 2000 and neg / n_cov > 0.98

# ---------------- block for A ----------------
jset = set(gj)
mapped = np.array([g in jset for g in chal_genes])
print(f"challenge genes mapped to jurkat donor: {mapped.sum()}/18533; "
      f"donor genes absent from challenge: {len(gj) - len(set(gj) & set(chal_genes))}")
order = np.argsort(-mean_cpm, kind="stable")
rank = np.empty(n_genes, dtype=np.int64)
rank[order] = np.arange(n_genes)
mapped_ranks = rank[mapped]
cut = np.sort(mapped_ranks)[2199]
block_chal = mapped & (rank <= cut)
print(f"block cutoff rank {cut}, block size {block_chal.sum()}, "
      f"min block CPM {mean_cpm[block_chal].min():.1f}, "
      f"max excluded-mapped CPM {mean_cpm[mapped & ~block_chal].max():.1f}")
assert block_chal.sum() == 2200
chal_idx_of_j = {g: i for i, g in enumerate(chal_genes)}
block_j = np.array(sorted(chal_idx_of_j[g] for g in
                          [chal_genes[i] for i in np.where(block_chal)[0]])) \
    if False else None
# donor-axis block indices:
j_idx_of = {g: i for i, g in enumerate(gj)}
block_j = np.array([j_idx_of[chal_genes[i]] for i in np.where(block_chal)[0]])
print("block donor-axis size:", block_j.shape)
# challenge-wide top-2200 coverage context
top_all = set(order[:2200].tolist())
print(f"of challenge-wide top-2200 CPM genes, {len(top_all & set(np.where(block_chal)[0]))} "
      f"are in the mapped block (rest unmapped to donor)")

# per-target heads on jurkat axis
HEAD = 400
Ahead = np.argpartition(-np.abs(Djm), HEAD - 1, axis=1)[:, :HEAD]
is_head = np.zeros_like(Djm, dtype=bool)
is_head[np.arange(Djm.shape[0])[:, None], Ahead] = True
ov = is_head[:, block_j].sum(axis=1)
print(f"head-block overlap per target: mean {ov.mean():.1f} min {ov.min()} max {ov.max()}")
print(f"mean eval entries/target: {(2200 - ov).mean():.1f}")
print(f"bank-wide exact-zero delta entries: {int((Djm == 0).sum())}")

# ---------------- split-half engine ----------------
TAUS = [0.05, 0.10, 0.15, 0.20]

def fit_priors(Dh1):
    P = (Dh1 < 0).mean(axis=0)
    a = np.abs(Dh1)
    denom = a.sum(axis=0)
    num = (a * (Dh1 < 0)).sum(axis=0)
    Pw = np.divide(num, denom, out=np.full(Dh1.shape[1], 0.5),
                   where=denom > 0)
    return P, Pw, int((denom == 0).sum())

def split_eval(D, block_idx, is_head_full, seed, n_rep, label):
    n_t, n_g = D.shape
    B = len(block_idx)
    rng = np.random.default_rng(seed)
    out = {}  # key -> list over reps
    for r in range(n_rep):
        perm = rng.permutation(n_t)
        h1, h2 = perm[:n_t // 2], perm[n_t // 2:]
        P, Pw, n_zero_denom = fit_priors(D[h1])
        if n_zero_denom:
            print(f"[{label} rep{r}] NOTE {n_zero_denom} zero-denom genes -> Pw=0.5",
                  flush=True)
        D2 = D[h2][:, block_idx]          # (n2, B)
        H2 = is_head_full[h2][:, block_idx]
        Z = (D2 == 0)
        E = ~(H2 | Z)
        T = (D2 < 0)                       # true down
        W = np.abs(D2)
        Pb, Pwb = P[block_idx], Pw[block_idx]
        res = {}
        res["alldown"] = float((T & E).sum() / E.sum())
        res["alldown_w"] = float((W * (T & E)).sum() / (W * E).sum())
        for nm, Pv in (("prior", Pb), ("priorW", Pwb)):
            pred = (Pv > 0.5)[None, :]
            ok = ((pred == T) & E)
            res[nm] = float(ok.sum() / E.sum())
            res[nm + "_w"] = float((W * ok).sum() / (W * E).sum())
            for tau in TAUS:
                app = (np.abs(Pv - 0.5) > tau)[None, :]
                predt = np.where(app, pred, True)
                okt = ((predt == T) & E)
                k = f"{nm}_t{tau:.2f}"
                res[k] = float(okt.sum() / E.sum())
                res[k + "_w"] = float((W * okt).sum() / (W * E).sum())
                res[k + "_cov"] = float((app * np.ones_like(T)).sum()
                                        / T.size)
                res[k + "_gcov"] = float((np.abs(Pv - 0.5) > tau).mean())
        sub = np.random.default_rng(9000 + r)
        Psh = P[sub.permutation(n_g)][block_idx]
        preds = (Psh > 0.5)[None, :]
        oks = ((preds == T) & E)
        res["shuffled"] = float(oks.sum() / E.sum())
        res["shuffled_w"] = float((W * oks).sum() / (W * E).sum())
        res["n_eval"] = float(E.sum())
        res["n_zero_dropped"] = float((Z & ~H2).sum())
        for k, v in res.items():
            out.setdefault(k, []).append(v)
        print(f"[{label} rep{r}] alldown {res['alldown']:.4f} prior "
              f"{res['prior']:.4f} prior_w {res['prior_w']:.4f} shuf "
              f"{res['shuffled']:.4f} neval {res['n_eval']:.0f}", flush=True)
    return out

print("== jurkat split-half ==", flush=True)
jout = split_eval(Djm, block_j, is_head, 20260910, 5, "J")

# ---------------- H1 analogue ----------------
print("== H1 analogue ==", flush=True)
Dh, gh, ph, bh = load_donor("h1_2025")
Dhm = np.asarray(Dh)
ref = np.load(os.path.join(ROOT,
    "experiments/323-current-family-validation/work/reference.npz"))
assert gh == ref["genes"].tolist()
keep = ref["keep"]
ref_lin = ref["ref_lin"]
univ = ref["genes"][keep]
hb = np.argsort(-ref_lin, kind="stable")[:2200]
# universe-position -> donor-axis index
u_pos = np.where(keep)[0]
block_h = u_pos[hb]
print("H1 block size:", block_h.shape, "min ref_lin in block:",
      float(ref_lin[hb].min()))
Hhead = np.argpartition(-np.abs(Dhm), HEAD - 1, axis=1)[:, :HEAD]
is_head_h = np.zeros_like(Dhm, dtype=bool)
is_head_h[np.arange(Dhm.shape[0])[:, None], Hhead] = True
ovh = is_head_h[:, block_h].sum(axis=1)
print(f"H1 head-block overlap: mean {ovh.mean():.1f} min {ovh.min()} max {ovh.max()}")
gh_idx = {g: i for i, g in enumerate(gh)}
n_cov_h = neg_h = 0
for ti, t in enumerate(ph):
    j = gh_idx.get(t)
    if j is not None:
        n_cov_h += 1
        if Dhm[ti, j] < 0:
            neg_h += 1
print(f"H1 own-target coverage {n_cov_h}/{len(ph)}, frac<0 {neg_h / max(n_cov_h, 1):.4f}")
hout = split_eval(Dhm, block_h, is_head_h, 20260910, 5, "H")

# ---------------- cross-line Jurkat <-> K562 ----------------
print("== cross-line ==", flush=True)
Dk, gk, pk, bk = load_donor("k562_gwps")
print("k562 shape:", Dk.shape)
Dkm = np.asarray(Dk)
gk_idx = {g: i for i, g in enumerate(gk)}
n_cov_k = neg_k = 0
for ti, t in enumerate(pk):
    j = gk_idx.get(t)
    if j is not None:
        n_cov_k += 1
        if Dkm[ti, j] < 0:
            neg_k += 1
print(f"k562 own-target coverage {n_cov_k}/{len(pk)}, frac<0 {neg_k / max(n_cov_k, 1):.5f}")

shared = sorted(set(gj) & set(gk))
print("shared genes:", len(shared))
j_sh = np.array([j_idx_of[g] for g in shared])
k_idx_of = {g: i for i, g in enumerate(gk)}
k_sh = np.array([k_idx_of[g] for g in shared])
PJ_full, PJw_full, _ = fit_priors(Djm)
PK_full, PKw_full, _ = fit_priors(Dkm)
pj_s = PJ_full[j_sh]
pk_s = PK_full[k_sh]
print(f"P_down corr on shared: pearson {np.corrcoef(pj_s, pk_s)[0, 1]:.4f}, "
      f"spearman {np.corrcoef(np.argsort(np.argsort(pj_s)), np.argsort(np.argsort(pk_s)))[0, 1]:.4f}")
print(f"sign-call agreement: {np.mean((pj_s > 0.5) == (pk_s > 0.5)):.4f}, "
      f"frac J-down {np.mean(pj_s > 0.5):.4f} frac K-down {np.mean(pk_s > 0.5):.4f}")

def cross_eval(Dm, b, axis_sh, P_other_sh, label):
    # block = top-2200 b_ctrl among shared; head via own delta
    bsh = b[axis_sh]
    top = np.argsort(-bsh, kind="stable")[:2200]
    blk = axis_sh[top]
    H = np.argpartition(-np.abs(Dm), HEAD - 1, axis=1)[:, :HEAD]
    Hm = np.zeros_like(Dm, dtype=bool)
    Hm[np.arange(Dm.shape[0])[:, None], H] = True
    D2 = Dm[:, blk]
    E = ~(Hm[:, blk] | (D2 == 0))
    T = D2 < 0
    W = np.abs(D2)
    Pv = P_other_sh[top]
    pred = (Pv > 0.5)[None, :]
    ok = ((pred == T) & E)
    out = {
        "alldown": float((T & E).sum() / E.sum()),
        "alldown_w": float((W * (T & E)).sum() / (W * E).sum()),
        "cross_prior": float(ok.sum() / E.sum()),
        "cross_prior_w": float((W * ok).sum() / (W * E).sum()),
        "n_eval": float(E.sum()),
    }
    print(f"[{label}] " + " ".join(f"{k}={v:.4f}" for k, v in out.items()
                                   if k != "n_eval") + f" neval={out['n_eval']:.0f}",
          flush=True)
    return out

# Jurkat prior -> K562 entries ; K562 prior -> Jurkat entries
x_j2k = cross_eval(Dkm, bk, k_sh, pj_s, "Jprior_on_K562")
x_k2j = cross_eval(Djm, bj, j_sh, pk_s, "Kprior_on_Jurkat")

# ---------------- production prior file (ALL jurkat targets) ----------------
np.savez(os.path.join(WORK, "jurkat_block_prior.npz"),
         genes=np.array(gj), P_down=PJ_full, Pw_down=PJw_full,
         n_targets=np.int64(len(pj)), seed=np.int64(20260910),
         block_A_genes=np.array([chal_genes[i] for i in np.where(block_chal)[0]]),
         block_A_donor_idx=block_j,
         block_A_mean_cpm=mean_cpm[block_chal],
         head_k=np.int64(HEAD))
print("wrote jurkat_block_prior.npz")

# ---------------- machine-readable metrics ----------------
def summ(d):
    return {k: {"mean": float(np.mean(v)), "sd": float(np.std(v, ddof=1)),
                "n": len(v), "values": [float(x) for x in v]}
            for k, v in d.items()}

metrics = {
    "jurkat_Ablock": summ(jout),
    "h1_Hblock": summ(hout),
    "cross_Jprior_on_K562": x_j2k,
    "cross_Kprior_on_Jurkat": x_k2j,
    "coverage": {
        "challenge_mapped_to_jurkat": int(mapped.sum()),
        "donor_genes_absent_from_challenge": int(len(gj) - len(set(gj) & set(chal_genes))),
        "block_size": 2200,
        "jurkat_head_block_overlap_mean": float(ov.mean()),
        "h1_head_block_overlap_mean": float(ovh.mean()),
        "jurkat_own_target_cov": int(n_cov), "jurkat_own_target_neg": float(neg / n_cov),
        "h1_own_target_cov": int(n_cov_h), "h1_own_target_neg": float(neg_h / max(n_cov_h, 1)),
        "k562_own_target_cov": int(n_cov_k), "k562_own_target_neg": float(neg_k / max(n_cov_k, 1)),
        "contextA_median_lib": med_lib,
        "shared_JK_genes": len(shared),
    },
}
json.dump(metrics, open(os.path.join(WORK, "metrics370.json"), "w"), indent=1)
with open(os.path.join(WORK, "per-split370.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["experiment", "repeat", "estimator", "acc_unw", "acc_w", "coverage"])
    for exp, d in (("jurkat", jout), ("h1", hout)):
        reps = len(d["alldown"])
        base = ["alldown", "prior", "priorW", "shuffled"] + \
            [f"prior_t{t:.2f}" for t in TAUS] + [f"priorW_t{t:.2f}" for t in TAUS]
        keys = [k for k in base if k in d]
        for r in range(reps):
            for k in keys:
                cov = d.get(k + "_cov", [""] * reps)[r]
                accw = d.get(k + "_w", [""] * reps)[r]
                w.writerow([exp, r, k, d[k][r], accw, cov])
print("done")

# ---------------- block-definition addendum (cheap, full-data) ----------------
# Why is all-down 0.5147 on the A-CPM block but 0.5374 on the donor-b_ctrl block?
mapped_j = np.where(mapped)[0]
b_mapped = bj[np.array([j_idx_of[chal_genes[i]] for i in mapped_j])]
top_ctrl = top_b = None
ord_cpm = np.argsort(-mean_cpm[mapped_j], kind="stable")[:2200]
ord_b = np.argsort(-b_mapped, kind="stable")[:2200]
set_cpm = set(mapped_j[ord_cpm].tolist())
set_b = set(mapped_j[ord_b].tolist())
print(f"A-CPM block vs donor-b_ctrl block (both among {len(mapped_j)} mapped): "
      f"overlap {len(set_cpm & set_b)}/2200")
E2 = ~is_head[:, :]
for nm, S in (("A-CPM", set_cpm), ("b_ctrl", set_b)):
    cols = np.array(sorted(S))
    blk = np.array([j_idx_of[chal_genes[i]] for i in cols])
    D2 = Djm[:, blk]
    E = ~(is_head[:, blk] | (D2 == 0))
    print(f"  {nm}: alldown {float(((D2 < 0) & E).sum() / E.sum()):.4f}, "
          f"full-prior acc {float((((PJ_full[blk] > 0.5)[None, :] == (D2 < 0)) & E).sum() / E.sum()):.4f} "
          f"(in-sample, all 2270 targets; for comparison only)")
# full-data P_down distribution on the A block
Pb_full = PJ_full[block_j]
print(f"full-data P_down on A block: frac up-leaning (P<0.5) "
      f"{float(np.mean(Pb_full < 0.5)):.4f}, deciles "
      f"{np.round(np.quantile(Pb_full, np.arange(0, 1.01, 0.1)), 3).tolist()}, "
      f"frac |P-0.5|>0.10 {float(np.mean(np.abs(Pb_full - 0.5) > 0.10)):.4f}")
