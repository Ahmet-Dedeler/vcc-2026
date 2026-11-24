"""350: same-line priors — inventory + common vectors (chunked streaming).

Reads:
  data/external/nadig2025/{jurkat,hepg2}_raw_singlecell.h5ad  (dense float32 X)
  data/challenge/{pert_counts.csv,gene_names.csv,context_A.h5ad}
  experiments/323-current-family-validation/work/effects.npz (cap15 design)
  data/donors/k562_gwps/{delta.npy,genes.txt} (optional donor-mean comparison)
Writes (all under experiments/350-same-line-priors-data/work/):
  inventory.json, jurkat_common.npz, hepg2_common.npz, compare.json, summary.txt
Memory: one dense block (8192 x G x 4B <= ~320MB) + float64 accumulators (301 x G).
"""
import h5py
import json
import numpy as np
import anndata as ad
from pathlib import Path
from scipy.stats import spearmanr

ROOT = Path(".")
WORK = ROOT / "experiments/350-same-line-priors-data/work"
WORK.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(350)
BLOCK = 8192


def dec(a):
    return np.array([x.decode() if isinstance(x, (bytes, bytearray)) else str(x)
                     for x in a])


def load_headers(path):
    """h5py-only inventory; returns dict + open-file-light objects needed later."""
    f = h5py.File(path, "r")
    X = f["X"]
    obs_cols = list(f["obs"].keys())
    var_cols = list(f["var"].keys())
    obs_targets = dec(f["obs"]["__categories"]["gene"][:])
    gene_cats = dec(f["var"]["__categories"]["gene_name"][:])
    var_codes = np.asarray(f["var"]["gene_name"][:]).astype(np.int64)
    var_syms = gene_cats[var_codes]  # per-column symbols; may contain 1 duplicate
    n_cols = int(X.shape[1])
    var_ids = np.asarray(f["var"]["gene_id"][:])
    # gene_id may be categorical codes or strings
    try:
        gid_cats = dec(f["var"]["__categories"]["gene_id"][:])
        gid_head = list(gid_cats[:3])
        gid_type = f"categorical Ensembl ({len(gid_cats)} levels)"
    except Exception:
        gid_head = [str(x) for x in np.asarray(var_ids[:3]).ravel()]
        gid_type = str(var_ids.dtype)
    codes = np.asarray(f["obs"]["gene"][:]).astype(np.int64)
    ctrl_idx = list(obs_targets).index("non-targeting")
    n_ctrl = int((codes == ctrl_idx).sum())
    umi = np.asarray(f["obs"]["UMI_count"][:], dtype=np.float64)
    info = dict(
        path=str(path), n_cells=int(X.shape[0]), n_genes=n_cols,
        n_unique_symbols=int(len(set(var_syms))),
        X_dtype=str(X.dtype), X_chunks=None if X.chunks is None else list(X.chunks),
        obs_columns=obs_cols, var_columns=var_cols,
        pert_label_col="gene", n_pert_levels=int(len(obs_targets)),
        n_targets_excl_control=int(len(obs_targets) - 1),
        control_label="non-targeting", n_control_cells=n_ctrl,
        gene_id_type=gid_type, gene_id_head=gid_head,
        var_symbols_unique=int(len(set(var_syms))),
        umi_median=float(np.median(umi)), umi_mean=float(np.mean(umi)),
    )
    f.close()
    levels, level_counts = np.unique(codes, return_counts=True)
    count_by_code = dict(zip(levels.tolist(), level_counts.tolist()))
    return info, obs_targets, var_syms, codes, count_by_code


def stream_group_sums(path, want_codes, n_genes, block=BLOCK):
    """Sum X rows per wanted group code. Returns (sums[K,G] float64, counts[K])."""
    K = len(want_codes)
    code_to_row = {c: i for i, c in enumerate(want_codes)}
    sums = np.zeros((K, n_genes), dtype=np.float64)
    counts = np.zeros(K, dtype=np.int64)
    # need per-row codes; reload codes inside to keep function self-contained
    f = h5py.File(path, "r")
    codes = np.asarray(f["obs"]["gene"][:]).astype(np.int64)
    X = f["X"]
    n = X.shape[0]
    lut = np.full(int(codes.max()) + 1, -1, dtype=np.int64)
    for i, c in enumerate(want_codes):
        lut[c] = i
    for lo in range(0, n, block):
        hi = min(lo + block, n)
        blk = np.asarray(X[lo:hi], dtype=np.float64)
        rows = lut[codes[lo:hi]]
        ok = rows >= 0
        if ok.any():
            r = rows[ok]
            b = blk[ok]
            for g in np.unique(r):
                sums[g] += b[r == g].sum(axis=0)
            cc = np.bincount(r, minlength=K)
            counts += cc
        del blk
        if (lo // block) % 4 == 0:
            print(f"    {hi:,}/{n:,}", flush=True)
    f.close()
    return sums, counts


def common_vector(path, obs_targets, var_syms, panel_set, axis, count_by_code,
                  n_pick=300, min_cells=20):
    idx_of = {t: i for i, t in enumerate(obs_targets)}
    targets = [t for t in obs_targets
               if t != "non-targeting" and count_by_code.get(idx_of[t], 0) >= min_cells]
    rng.shuffle(targets)
    picked = sorted(targets[:n_pick])
    want = [list(obs_targets).index("non-targeting")] + \
           [list(obs_targets).index(t) for t in picked]
    print(f"  streaming {path} ...", flush=True)
    sums, counts = stream_group_sums(path, want, len(var_syms))
    # collapse duplicate symbols (Nadig carries two HSPA14 columns), as in
    # scripts/build_donor_table.py
    uniq, inv = np.unique(np.asarray(var_syms), return_inverse=True)
    if len(uniq) != len(var_syms):
        wide = np.zeros((sums.shape[0], len(uniq)), dtype=np.float64)
        for u in range(len(uniq)):
            wide[:, u] = sums[:, inv == u].sum(axis=1)
        sums = wide
        var_syms = uniq
        print(f"  collapsed {len(inv)-len(uniq)} duplicate symbol(s)", flush=True)
    ctrl_sum, pert_sums = sums[0], sums[1:]
    ctrl_n, pert_ns = counts[0], counts[1:]
    assert (pert_ns >= 20).all(), f"picked a low-cell target: min {pert_ns.min()}"
    # CPM of means, scorer's lfc formula
    cpm_ctrl = 1e6 * ctrl_sum / ctrl_sum.sum()
    cpm_pert = 1e6 * pert_sums / pert_sums.sum(axis=1, keepdims=True)
    lfc = np.log2((cpm_pert + 1e-9) / (cpm_ctrl + 1e-9))  # (300, G_bank)
    # map bank genes -> challenge axis
    axis_index = {g: i for i, g in enumerate(axis)}
    bank_to_axis = np.array([axis_index.get(g, -1) for g in var_syms])
    keep_bank = np.flatnonzero(bank_to_axis >= 0)
    ax_idx = bank_to_axis[keep_bank]
    order = np.argsort(ax_idx)
    lfc_ax = np.zeros((lfc.shape[0], len(axis)), dtype=np.float64)
    lfc_ax[:, ax_idx[order]] = lfc[:, keep_bank[order]]
    common_full = lfc.mean(axis=0)  # bank axis
    common_ax = lfc_ax.mean(axis=0)  # challenge axis (zeros where unmapped)
    mapped = np.zeros(len(axis), dtype=bool)
    mapped[ax_idx] = True
    # propensity on challenge axis
    propensity = np.zeros(len(axis))
    propensity[ax_idx] = (np.abs(lfc)[:, keep_bank] > 0.5).mean(axis=0)
    stats = dict(
        targets=picked, n_cells_ctrl=int(ctrl_n),
        cells_per_target=dict(min=int(pert_ns.min()),
                              median=float(np.median(pert_ns)),
                              max=int(pert_ns.max())),
        frac_entries_down=float((lfc < 0).mean()),
        frac_targets_mean_down=float((lfc.mean(axis=1) < 0).mean()),
        frac_common_down_bank=float((common_full < 0).mean()),
        n_bank_genes=int(len(var_syms)), n_mapped=int(mapped.sum()),
    )
    return lfc, common_full, common_ax, mapped, propensity, stats, var_syms


print("== headers ==")
panel = [l.strip() for l in open(ROOT / "data/challenge/pert_counts.csv")][1:]
axis = [l.strip() for l in open(ROOT / "data/challenge/gene_names.csv")][1:]
panel_set, axis_set = set(panel), set(axis)
inv = dict(panel_n=len(panel), axis_n=len(axis))
banks = {}
for name, path in [("jurkat", ROOT / "data/external/nadig2025/jurkat_raw_singlecell.h5ad"),
                   ("hepg2", ROOT / "data/external/nadig2025/hepg2_raw_singlecell.h5ad")]:
    info, obs_targets, var_syms, codes, cbc = load_headers(path)
    obs_set = set(obs_targets) - {"non-targeting"}
    info["panel_overlap_n"] = len(set(panel) & obs_set)
    info["panel_overlap"] = sorted(set(panel) & obs_set)
    info["var_axis_overlap_n"] = len(set(var_syms) & axis_set)
    info["targets_axis_overlap_n"] = len(obs_set & axis_set)
    info["n_targets_ge20"] = sum(1 for t in obs_set
                                if cbc.get(list(obs_targets).index(t), 0) >= 20)
    inv[name] = info
    banks[name] = (str(path), obs_targets, var_syms, cbc)
    print(name, info["n_cells"], "cells x", info["n_genes"], "genes;",
          info["n_control_cells"], "controls; panel overlap",
          info["panel_overlap_n"], "; var-axis overlap",
          info["var_axis_overlap_n"])
json.dump(inv, open(WORK / "inventory.json", "w"), indent=2)

print("== jurkat common ==")
lfc_j, cf_j, ca_j, map_j, prop_j, st_j, vs_j = common_vector(
    banks["jurkat"][0], banks["jurkat"][1], banks["jurkat"][2], panel_set, axis,
    banks["jurkat"][3])
print("== hepg2 common ==")
lfc_h, cf_h, ca_h, map_h, prop_h, st_h, vs_h = common_vector(
    banks["hepg2"][0], banks["hepg2"][1], banks["hepg2"][2], panel_set, axis,
    banks["hepg2"][3])
print("jurkat", st_j)
print("hepg2", st_h)

# top genes of common vectors (bank axis)
for nm, cf, vs in [("jurkat", cf_j, vs_j), ("hepg2", cf_h, vs_h)]:
    o = np.argsort(cf)
    print(f"-- {nm} most DOWN:", [(vs[i], round(float(cf[i]), 3)) for i in o[:25]])
    print(f"-- {nm} most UP:", [(vs[i], round(float(cf[i]), 3)) for i in o[-25:][::-1]])

np.savez_compressed(WORK / "jurkat_common.npz", genes=np.array(axis),
                    vector=ca_j, mapped=map_j, propensity=prop_j,
                    targets=np.array(st_j["targets"]), bank_genes=vs_j,
                    common_bank=cf_j, lfc_bank_mean_note="log2 CPM-ratio vs controls")
np.savez_compressed(WORK / "hepg2_common.npz", genes=np.array(axis),
                    vector=ca_h, mapped=map_h, propensity=prop_h,
                    targets=np.array(st_h["targets"]), bank_genes=vs_h,
                    common_bank=cf_h, lfc_bank_mean_note="log2 CPM-ratio vs controls")

# H1 design common vector
e = np.load(ROOT / "experiments/323-current-family-validation/work/effects.npz")
cap = e["cap15"]  # 85 x 18533, challenge axis already
h1_common = cap.mean(axis=0)
print("h1 cap15 mean norm:", float(np.linalg.norm(h1_common)),
      "frac down:", float((h1_common < 0).mean()))


def cos(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


def signagree(a, b):
    m = (a != 0) & (b != 0)
    return float((np.sign(a[m]) == np.sign(b[m])).mean()), int(m.sum())


shared_jh = map_j & map_h
shared_jh1 = shared_jh  # h1 covers full axis
out = dict(
    cos_jurkat_hepg2=cos(ca_j[shared_jh], ca_h[shared_jh]),
    sign_jurkat_hepg2=signagree(ca_j[shared_jh], ca_h[shared_jh]),
    cos_jurkat_h1=cos(ca_j[map_j], h1_common[map_j]),
    sign_jurkat_h1=signagree(ca_j[map_j], h1_common[map_j]),
    cos_hepg2_h1=cos(ca_h[map_h], h1_common[map_h]),
    sign_hepg2_h1=signagree(ca_h[map_h], h1_common[map_h]),
    n_shared_jh=int(shared_jh.sum()), n_mapped_j=int(map_j.sum()),
    n_mapped_h=int(map_h.sum()),
)
# donor-bank mean comparison (log1p pseudobulk space — caveated)
for dname in ["k562_gwps"]:
    gp = ROOT / f"data/donors/{dname}/genes.txt"
    dp = ROOT / f"data/donors/{dname}/delta.npy"
    if gp.exists() and dp.exists():
        dg = [l.strip() for l in open(gp)]
        dmean = np.load(dp, mmap_mode="r").mean(axis=0).astype(np.float64)
        aidx = {g: i for i, g in enumerate(axis)}
        rows = [(aidx[g], j) for j, g in enumerate(dg) if g in aidx]
        ai = np.array([r[0] for r in rows])
        dm_ax = np.zeros(len(axis))
        dm_ax[ai] = dmean[[r[1] for r in rows]]
        mm = np.zeros(len(axis), dtype=bool)
        mm[ai] = True
        mj = map_j & mm
        mh = map_h & mm
        out[f"cos_jurkat_{dname}"] = cos(ca_j[mj], dm_ax[mj])
        out[f"sign_jurkat_{dname}"] = signagree(ca_j[mj], dm_ax[mj])
        out[f"cos_hepg2_{dname}"] = cos(ca_h[mh], dm_ax[mh])
        out[f"cos_h1_{dname}"] = cos(h1_common[mm], dm_ax[mm])
        out[f"n_shared_donor"] = int(mm.sum())
        print(dname, "shared", mm.sum())
    else:
        print(dname, "missing — skipped")
print(json.dumps(out, indent=2))
json.dump(out, open(WORK / "compare.json", "w"), indent=2)

# context-A basal CPM (chunked, backed CSR) + Spearman with Jurkat propensity
print("== context A basal ==")
a = ad.read_h5ad(ROOT / "data/challenge/context_A.h5ad", backed="r")
n = a.n_obs
tot = np.zeros(n)
sums = None
CS = 2000
for lo in range(0, n, CS):
    blk = a.X[lo:min(lo + CS, n)].toarray().astype(np.float64)
    tot[lo:min(lo + CS, n)] = blk.sum(axis=1)
    sums = blk.sum(axis=0) if sums is None else sums + blk.sum(axis=0)
    if (lo // CS) % 4 == 0:
        print(f"    {min(lo+CS,n):,}/{n:,}", flush=True)
cpm = 1e6 * sums / sums.sum()
print("basal CPM median:", float(np.median(cpm)), "frac>5:", float((cpm > 5).mean()))
m = map_j  # jurkat-mapped genes
rho, pval = spearmanr(prop_j[m], cpm[m])
print("Spearman(propensity_jurkat, basalCPM_A):", float(rho), "p:", float(pval),
      "n:", int(m.sum()))
mh = map_h
rho_h, pval_h = spearmanr(prop_h[mh], cpm[mh])
print("Spearman(propensity_hepg2, basalCPM_A):", float(rho_h), "p:", float(pval_h),
      "n:", int(mh.sum()))
out["spearman_propJ_basalA"] = [float(rho), float(pval), int(m.sum())]
out["spearman_propH_basalA"] = [float(rho_h), float(pval_h), int(mh.sum())]
out["jurkat_stats"] = st_j
out["hepg2_stats"] = st_h
json.dump(out, open(WORK / "compare.json", "w"), indent=2)
with open(WORK / "summary.txt", "w") as fh:
    fh.write(json.dumps(dict(inventory=inv, compare=out), indent=2))
print("DONE")
