"""347 -- sign-confidence strata on saved cap15 arrays. Descriptive only.

Known answers (343): mean NMAE 0.9508364630 over 58 eligible targets;
delivered w-acc 0.587968 (called 0.637666 / uncalled 0.551763);
designed w-acc 0.600185 (called 0.649299 / uncalled 0.558007);
agree design-vs-delivered 0.866750; mass_called 0.346388; m 0.205485.
"""
import json
import numpy as np, polars as pl
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
W = REPO / "experiments/323-current-family-validation/work"
OUT = Path(__file__).resolve().parent / "work"
OUT.mkdir(exist_ok=True)

s = np.load(W / "stats-cap15.npz")
r = np.load(W / "reference.npz")
e = np.load(W / "effects.npz")
genes_all = r["genes"].astype(str)
keep = r["keep"].astype(bool)
genes = genes_all[keep]
assert (s["targets"].astype(str) == r["targets"].astype(str)).all()
targets = s["targets"].astype(str)
assert (e["targets"].astype(str) == targets).all()
eg = {g: i for i, g in enumerate(e["genes"].astype(str))}
col = np.array([eg.get(g, -1) for g in genes])
miss = col < 0
print("universe genes absent from effects axis:", int(miss.sum()))
eff_full = np.asarray(e["cap15"], dtype=np.float64)
safe_col = np.where(miss, 0, col)
eff = eff_full[:, safe_col]
eff[:, miss] = 0.0
lhat = np.nan_to_num(np.asarray(s["lfc"], dtype=np.float64))
phat = np.asarray(s["p"], dtype=np.float64)
qhat = np.asarray(s["q"], dtype=np.float64)
ref_lin = np.asarray(r["ref_lin"], dtype=float)
ref = pl.read_parquet(W / "reference_de.parquet")
gi = {g: i for i, g in enumerate(genes)}

t_list, l_list, lh_list, p_list, q_list = [], [], [], [], []
ef_list, rlin_list, gene_list = [], [], []
elig = []
for ti, t in enumerate(targets):
    sub = ref.filter(pl.col("target") == t)
    g = sub["feature"].to_numpy().astype(str)
    q = sub["p_adj"].to_numpy().astype(float)
    l = sub["log2_fold_change"].to_numpy().astype(float)
    gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if gate.sum() < 10:
        continue
    elig.append(ti)
    idx = np.array([gi[x] for x in g[gate]])
    t_list.append(np.full(int(gate.sum()), ti))
    l_list.append(l[gate])
    lh_list.append(lhat[ti, idx])
    p_list.append(phat[ti, idx])
    q_list.append(qhat[ti, idx])
    ef_list.append(eff[ti, idx])
    rlin_list.append(ref_lin[idx])
    gene_list.append(g[gate])
Tt = np.concatenate(t_list)
Tl = np.concatenate(l_list)
Tlh = np.concatenate(lh_list)
Tp = np.concatenate(p_list)
Tq = np.concatenate(q_list)
Tef = np.concatenate(ef_list)
Trlin = np.concatenate(rlin_list)
Tgene = np.concatenate(gene_list)
N = len(Tl)
print("eligible targets: %d (expect 58); pooled gate pairs: %d" % (len(elig), N))
w = np.abs(Tl)
D = w.sum()
sd_des = (np.sign(Tef) == np.sign(Tl))
sd_del = (np.sign(Tlh) == np.sign(Tl))


def stats(mask, name):
    m = np.asarray(mask, dtype=bool)
    n = int(m.sum())
    wm = float(w[m].sum())
    if n == 0 or wm == 0:
        return dict(stratum=name, n=n, mass_share=0.0, design_wacc=float("nan"),
                    design_uacc=float("nan"), deliv_wacc=float("nan"),
                    deliv_uacc=float("nan"), mag_ratio=float("nan"))
    return dict(stratum=name, n=n, mass_share=wm / float(D),
                design_wacc=float((w[m] * sd_des[m]).sum() / wm),
                design_uacc=float(sd_des[m].mean()),
                deliv_wacc=float((w[m] * sd_del[m]).sum() / wm),
                deliv_uacc=float(sd_del[m].mean()),
                mag_ratio=float(np.abs(Tlh[m]).sum() / wm))


rows = []
# ---- known-answer reproduction (mean over targets, 343's convention) ----
nmaes, pw_d, pw_dd, pw_c, pw_u, dpw_c, dpw_u, mass_c = [], [], [], [], [], [], [], []
for ti in elig:
    m = Tt == ti
    l = Tl[m]
    Dp = np.abs(l).sum()
    nmaes.append(np.abs(Tlh[m] - l).sum() / Dp)
    pw_d.append((np.abs(l) * (np.sign(Tlh[m]) == np.sign(l))).sum() / Dp)
    pw_dd.append((np.abs(l) * (np.sign(Tef[m]) == np.sign(l))).sum() / Dp)
    c = Tq[m] < .05
    nc = np.abs(l[c]).sum()
    nu = np.abs(l[~c]).sum()
    pw_c.append((np.abs(l[c]) * (np.sign(Tlh[m][c]) == np.sign(l[c]))).sum() / max(nc, 1e-12))
    pw_u.append((np.abs(l[~c]) * (np.sign(Tlh[m][~c]) == np.sign(l[~c]))).sum() / max(nu, 1e-12))
    dpw_c.append((np.abs(l[c]) * (np.sign(Tef[m][c]) == np.sign(l[c]))).sum() / max(nc, 1e-12))
    dpw_u.append((np.abs(l[~c]) * (np.sign(Tef[m][~c]) == np.sign(l[~c]))).sum() / max(nu, 1e-12))
    mass_c.append(nc / Dp)
ka = dict(n_targets=len(elig), nmae=float(np.mean(nmaes)), deliv_wacc=float(np.mean(pw_d)),
          design_wacc=float(np.mean(pw_dd)), deliv_called=float(np.mean(pw_c)),
          deliv_uncalled=float(np.mean(pw_u)), design_called=float(np.mean(dpw_c)),
          design_uncalled=float(np.mean(dpw_u)), mass_called=float(np.mean(mass_c)),
          mag_ratio=float(np.abs(Tlh).sum() / float(D)))
print(json.dumps(ka, indent=1))

# ---- (a) design |effect| deciles: per-target-within-gate rank ----
ae = np.abs(Tef)
pt_rank = np.empty(N)
for ti in elig:
    m = np.where(Tt == ti)[0]
    pt_rank[m] = np.argsort(np.argsort(ae[m], kind="stable")) / max(len(m) - 1, 1)
for d in range(10):
    lo, hi = d / 10.0, (d + 1) / 10.0
    mm = (pt_rank >= lo) & (pt_rank <= hi + 1e-12) if d == 9 else (pt_rank >= lo) & (pt_rank < hi)
    rows.append(stats(mm, "a_pt_dec%d" % d))
# ---- (a) global deciles ----
qs = np.quantile(ae, np.linspace(0, 1, 11))
for d in range(10):
    lo, hi = qs[d], qs[d + 1]
    mm = (ae >= lo) & (ae <= hi) if d == 9 else (ae >= lo) & (ae < hi)
    rows.append(stats(mm, "a_glob_dec%d" % d))
print("global |ef| decile edges:", [str(round(float(x), 4)) for x in qs])

# ---- (b) delivered called vs not; delivered-p pooled deciles ----
rows.append(stats(Tq < .05, "b_called"))
rows.append(stats(~(Tq < .05), "b_uncalled"))
lpv = -np.log10(np.clip(Tp, 1e-300, 1.0))
qp = np.quantile(lpv, np.linspace(0, 1, 11))
for d in range(10):
    lo, hi = qp[d], qp[d + 1]
    mm = (lpv >= lo) & (lpv <= hi) if d == 9 else (lpv >= lo) & (lpv < hi)
    rows.append(stats(mm, "b_pdec%d" % d))
print("delivered -log10p decile edges:", [str(round(float(x), 3)) for x in qp])

# ---- (c) agreement with leave-one-out cross-target common vector ----
common = np.empty(N)
for ti in elig:
    m = np.where(Tt == ti)[0]
    mu = (eff_full.sum(axis=0) - eff_full[ti]) / (eff_full.shape[0] - 1)
    common[m] = mu[safe_col[np.array([gi[x] for x in Tgene[m]])]]
agree = np.sign(Tef) == np.sign(common)
zero_thr = float(np.quantile(np.abs(common), 0.10))
is_zero = np.abs(common) <= zero_thr
print("common~0 threshold (|common| p10):", round(zero_thr, 6),
      " agree-share-pairs:", round(float(agree.mean()), 4),
      " agree-mass-share:", round(float(w[agree].sum() / float(D)), 4))
rows.append(stats(agree & ~is_zero, "c_agree"))
rows.append(stats(~agree & ~is_zero, "c_disagree"))
rows.append(stats(is_zero, "c_common0"))
rows.append(stats(agree, "c_agree_incl0"))

# ---- (e) gene class prefixes ----
is_rpl = np.char.startswith(Tgene, "RPL")
is_rps = np.char.startswith(Tgene, "RPS")
is_mt = np.char.startswith(Tgene, "MT-")
is_hist = np.char.startswith(Tgene, "HIST")
rows.append(stats(is_rpl | is_rps, "e_ribo"))
rows.append(stats(is_mt, "e_mito"))
rows.append(stats(is_hist, "e_histone"))
rows.append(stats(is_rpl | is_rps | is_mt | is_hist, "e_special_all"))
rows.append(stats(~(is_rpl | is_rps | is_mt | is_hist), "e_rest"))
print("class counts: ribo=%d mito=%d hist=%d" % (int((is_rpl | is_rps).sum()),
      int(is_mt.sum()), int(is_hist.sum())))

# ---- (e) reconstructed emitter 2600-set: head400 by |design| + top2200 by ref_lin ----
in_call = np.zeros(N, dtype=bool)
in_head = np.zeros(N, dtype=bool)
for ti in elig:
    m = np.where(Tt == ti)[0]
    o = np.argsort(-ae[m], kind="stable")
    head = np.zeros(len(m), dtype=bool)
    head[o[:400]] = True
    rest_pool = np.where(~head)[0]
    rest_order = rest_pool[np.argsort(-Trlin[m][rest_pool], kind="stable")[:2200]]
    call = head.copy()
    call[rest_order] = True
    in_call[m] = call
    in_head[m] = head
rows.append(stats(in_head, "e_head400"))
rows.append(stats(in_call & ~in_head, "e_cpmrest2200"))
rows.append(stats(in_call, "e_call2600"))
rows.append(stats(~in_call, "e_outside2600"))

# ---- (f) top-3 per-target deciles AND common-agree ----
top3 = pt_rank >= 0.7 - 1e-12
rows.append(stats(top3 & agree & ~is_zero, "f_top3_and_agree"))
rows.append(stats(top3, "f_top3"))

pl.DataFrame(rows).write_csv(OUT / "strata-cap15.csv")
with open(OUT / "known-answer-cap15.json", "w") as f:
    json.dump(ka, f, indent=1)
with open(OUT / "notes-c.txt", "w") as f:
    f.write("common0_thr=%r\n" % zero_thr)
print("wrote", str(OUT / "strata-cap15.csv"))

# ---- exact per-target NMAE counterfactuals for candidate strata ----
base_nmae = float(np.mean(nmaes))


def counterfactual(mask_all, s_in, s_out):
    vals = []
    for ti in elig:
        m = np.where(Tt == ti)[0]
        l = Tl[m]
        Dp = np.abs(l).sum()
        a = np.where(mask_all[m], s_in, s_out)
        vals.append(np.abs(a * Tlh[m] - l).sum() / Dp)
    return float(np.mean(vals))


cands = {
    "c_agree": (agree & ~is_zero),
    "c_disagree": (~agree & ~is_zero),
    "e_call2600": in_call,
    "e_head400": in_head,
    "f_top3_and_agree": (top3 & agree & ~is_zero),
    "b_called": (Tq < .05),
    "e_ribo": (np.char.startswith(Tgene, "RPL") | np.char.startswith(Tgene, "RPS")),
}
cf_rows = [dict(base=base_nmae)]
for name, mk in cands.items():
    row = dict(stratum=name)
    for s_in, s_out in [(2.0, 1.0), (3.0, 1.0), (2.0, 0.5), (3.0, 0.5)]:
        row["x%.1f_o%.1f" % (s_in, s_out)] = counterfactual(mk, s_in, s_out)
    cf_rows.append(row)
pl.DataFrame(cf_rows).write_csv(OUT / "counterfactuals-cap15.csv")
for row in cf_rows:
    print(row)
