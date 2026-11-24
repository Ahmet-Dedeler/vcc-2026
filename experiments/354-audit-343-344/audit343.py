"""Independent re-derivation of 343's claims. numpy/pandas only; no 343 imports."""
import json
import numpy as np
import pandas as pd

W = '~/Code/virtual-cell-challenge-2026/experiments/323-current-family-validation/work/'
OUT = {}

ref = np.load(W + 'reference.npz', mmap_mode='r', allow_pickle=True)
st = np.load(W + 'stats-cap15.npz', allow_pickle=True)
ef = np.load(W + 'effects.npz', mmap_mode='r', allow_pickle=True)

genes80 = np.asarray(ref['genes'])
keep = np.asarray(ref['keep'])
uni = genes80[keep]                       # 10,779 universe
targets = np.asarray(ref['targets'])      # 85
rlfc = np.asarray(ref['lfc'])             # 85 x 10779
rq = np.asarray(ref['q'])
basal = np.asarray(ref['ref_lin'])        # 10,779 basal CPM
plfc = np.asarray(st['lfc'])              # ours, delivered
pq = np.asarray(st['q'])
assert (np.asarray(st['targets']) == targets).all()

# design matrix mapped by gene name
egenes = np.asarray(ef['genes'])
assert (np.asarray(ef['targets']) == targets).all()
cap = np.asarray(ef['cap15'])             # 85 x 18533
eidx = {g: i for i, g in enumerate(egenes)}
dsgn = np.zeros((85, len(uni)))
missing = 0
for j, g in enumerate(uni):
    if g in eidx:
        dsgn[:, j] = cap[:, eidx[g]]
    else:
        missing += 1
OUT['design_missing_universe_genes'] = missing

# parquet cross-check on one target
par = pd.read_parquet(W + 'reference_de.parquet')
p1 = par[par.target == targets[0]].set_index('feature')
OUT['parquet_vs_npz_max_lfc_diff'] = float(
    np.abs(p1.loc[uni, 'log2_fold_change'].to_numpy() - rlfc[0]).max())
OUT['parquet_vs_npz_max_q_diff'] = float(
    np.abs(p1.loc[uni, 'p_adj'].to_numpy() - rq[0]).max())

# ---- gates ----
gates = []   # boolean masks per target
elig = []
for i, t in enumerate(targets):
    g = (rq[i] < 0.05) & np.isfinite(rlfc[i]) & (rlfc[i] != 0) & (uni != t)
    gates.append(g)
    if g.sum() >= 10:
        elig.append(i)
gates = np.asarray(gates)
OUT['n_eligible'] = len(elig)
OUT['pooled_gate_genes'] = int(gates[elig].sum())
OUT['mean_gate_size'] = float(gates[elig].sum(axis=1).mean())

def nmae(a, b):
    return np.abs(a - b).sum() / np.abs(b).sum()

per_nmae = np.array([nmae(plfc[i][gates[i]], rlfc[i][gates[i]]) for i in elig])
OUT['mean_nmae'] = float(per_nmae.mean())

rec = json.load(open(W + 'metrics-cap15.json'))
OUT['recorded_nmae'] = float(rec['raw']['nmae'])
pt = pd.read_csv(W + 'per-target-cap15.csv')
pn = pt[pt.metric == 'de_wilcoxon_lfc_nmae']['value'].to_numpy()
OUT['per_target_csv_mean'] = float(pn.mean())
OUT['per_target_csv_n'] = int(len(pn))

# ---- claim 2: magnitude ratio + weighted sign accuracy (pooled and mean-of-targets) ----
num_a = num_b = den = 0.0
w_acc_del_pool = w_acc_dsg_pool = 0.0
mratios, wdel, wdsg = [], [], []
for i in elig:
    g = gates[i]
    l, h, d = rlfc[i][g], plfc[i][g], dsgn[i][g]
    al = np.abs(l)
    den += al.sum(); num_a += np.abs(h).sum(); num_b += 0  # placeholder
    w_acc_del_pool += (np.sign(h) == np.sign(l)).astype(float).dot(al)
    w_acc_dsg_pool += (np.sign(d) == np.sign(l)).astype(float).dot(al)
    mratios.append(np.abs(h).sum() / al.sum())
    wdel.append(((np.sign(h) == np.sign(l))).astype(float).dot(al) / al.sum())
    wdsg.append(((np.sign(d) == np.sign(l))).astype(float).dot(al) / al.sum())
OUT['mag_ratio_pooled'] = float(num_a / den)
OUT['mag_ratio_mean_of_targets'] = float(np.mean(mratios))
OUT['wacc_delivered_pooled'] = float(w_acc_del_pool / den)
OUT['wacc_delivered_mean'] = float(np.mean(wdel))
OUT['wacc_designed_pooled'] = float(w_acc_dsg_pool / den)
OUT['wacc_designed_mean'] = float(np.mean(wdsg))

# ---- claim 3: counterfactuals ----
cf_sign, cf_ls, cf_grid, best_a, ls_a = [], [], [], [], []
grid = np.arange(0.25, 6.0 + 1e-9, 0.25)
for i in elig:
    g = gates[i]
    l, h = rlfc[i][g], plfc[i][g]
    cf_sign.append(nmae(np.sign(l) * np.abs(h), l))
    a_ls = float(h.dot(l) / h.dot(h))
    ls_a.append(a_ls)
    cf_ls.append(nmae(a_ls * h, l))
    vals = [nmae(a * h, l) for a in grid]
    k = int(np.argmin(vals))
    cf_grid.append(vals[k]); best_a.append(float(grid[k]))
OUT['cf_all_signs_right'] = float(np.mean(cf_sign))
OUT['cf_ls_scalar'] = float(np.mean(cf_ls))
OUT['cf_ls_mean_a'] = float(np.mean(ls_a))
OUT['cf_grid_best'] = float(np.mean(cf_grid))
OUT['cf_grid_mean_best_a'] = float(np.mean(best_a))

# ---- claim 4: pool composition ----
Gl = np.concatenate([rlfc[i][gates[i]] for i in elig])
Gh = np.concatenate([plfc[i][gates[i]] for i in elig])
Gb = np.concatenate([basal[gates[i]] for i in elig])
Gq = np.concatenate([pq[i][gates[i]] for i in elig])
OUT['up_share_unweighted'] = float((Gl > 0).mean())
OUT['up_share_lfc_weighted'] = float(np.abs(Gl[Gl > 0]).sum() / np.abs(Gl).sum())
qs = np.quantile(Gb, [0, .2, .4, .6, .8, 1.0])
OUT['basal_quintile_edges_pooled_gate'] = [float(x) for x in qs]
called = Gq < 0.05
for b in range(5):
    m = (Gb >= qs[b]) & (Gb <= qs[b + 1] if b == 4 else Gb < qs[b + 1])
    OUT[f'called_frac_basalQ{b}'] = float(called[m].mean())
    OUT[f'mag_ratio_basalQ{b}'] = float(np.abs(Gh[m]).sum() / np.abs(Gl[m]).sum())
aq = np.quantile(np.abs(Gl), [0, .2, .4, .6, .8, 1.0])
OUT['abs_lfc_quintile_edges'] = [float(x) for x in aq]
m = np.abs(Gl) >= aq[4]
OUT['mag_ratio_top_absLFC_quintile'] = float(np.abs(Gh[m]).sum() / np.abs(Gl[m]).sum())

# ---- claim 5: UP-only accuracies ----
def wacc(l, x):
    al = np.abs(l)
    return float(((np.sign(x) == np.sign(l))).astype(float).dot(al) / al.sum())
up = Gl > 0
OUT['wacc_UP_delivered'] = wacc(Gl[up], Gh[up])
Gd = np.concatenate([dsgn[i][gates[i]] for i in elig])
OUT['wacc_UP_designed'] = wacc(Gl[up], Gd[up])
low = up & (Gb < 24)
OUT['wacc_UP_lowbasal_delivered'] = wacc(Gl[low], Gh[low])
OUT['wacc_UP_lowbasal_designed'] = wacc(Gl[low], Gd[low])
OUT['frac_gate_lowbasal_UP'] = float(low.sum() / len(Gl))

# ---- claim 6: emitter flip-rate bins over universe x 85 ----
uq = np.quantile(basal, [0, .2, .4, .6, .8, 1.0])
OUT['universe_basal_quintile_edges'] = [float(x) for x in uq]
Hall = plfc.ravel()
Dall = dsgn.ravel()
Ball = np.tile(basal, 85)
for b in range(5):
    lo, hi = uq[b], uq[b + 1]
    m = (Ball >= lo) & (Ball <= hi if b == 4 else Ball < hi)
    d = Dall[m] > 0
    OUT[f'P_delneg_given_despos_bin{b}'] = float((Hall[m][d] < 0).mean())
    OUT[f'bin{b}_edge'] = [float(lo), float(hi)]
OUT['frac_sig_calls_basal_lt10'] = float((Gb < 10).mean())
OUT['nmae_mass_basal_lt10_denom'] = float(np.abs(Gl[Gb < 10]).sum() / np.abs(Gl).sum())
OUT['nmae_mass_basal_lt10_numer'] = float(
    np.abs(Gh[Gb < 10] - Gl[Gb < 10]).sum() / np.abs(Gh - Gl).sum())

json.dump(OUT, open('~/Code/virtual-cell-challenge-2026/experiments/354-audit-343-344/work/audit343.json', 'w'), indent=1)
print(json.dumps(OUT, indent=1))
