"""371 — Consensus block prior LOLO. See PLAN.md. Read-only on all inputs;
writes only into this experiment's work/. mmap + row-chunked streaming."""
import os, json, csv
import numpy as np

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

ROOT = "~/Code/virtual-cell-challenge-2026"
EXP = os.path.join(ROOT, "experiments/371-consensus-block-prior-lolo")
WORK = os.path.join(EXP, "work")
os.makedirs(WORK, exist_ok=True)

BANKS = [
    ("k562_gwps", os.path.join(ROOT, "data/donors/k562_gwps")),
    ("k562_essential", os.path.join(ROOT, "data/donors/k562_essential")),
    ("nadig_jurkat", os.path.join(ROOT, "data/donors/nadig_jurkat")),
    ("h1_2025", os.path.join(ROOT, "data/donors/h1_2025")),
    ("orion_hct116", os.path.join(ROOT, "data/donors/194-debias-orion_hct116")),
    ("orion_hek293t", os.path.join(ROOT, "data/donors/194-debias-orion_hek293t")),
    ("nadig_hepg2", os.path.join(ROOT,
     "experiments/166-donor-acquisition/donors/nadig_hepg2")),
]
HEAD = 400
BLOCK = 2200
CHUNK = 1000


def read_lines(p):
    with open(p) as f:
        return [l.strip() for l in f if l.strip() != ""]


# ---------------- 0. known-answer reproduction (own code) ----------------
print("== known answers ==", flush=True)
m15 = json.load(open(os.path.join(
    ROOT, "experiments/323-current-family-validation/work/metrics-cap15.json")))
per = {}
with open(os.path.join(ROOT, "experiments/323-current-family-validation/work/"
                         "per-target-cap15.csv")) as f:
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
assert abs(np.nanmean(np.array(
    per["de_wilcoxon_direction_fidelity_yield_raw"]))
    - m15["raw"]["fid"]) < 1e-9
assert abs(np.nanmean(np.array(
    per["de_wilcoxon_direction_reach_raw"]))
    - m15["raw"]["reach"]) < 1e-9
print("323 known-answer MATCH (nmae/fid/reach diff < 1e-9)")

import h5py
for ctx in ["A", "B", "C"]:
    fa = h5py.File(os.path.join(ROOT, f"data/challenge/context_{ctx}.h5ad"), "r")
    shape = tuple(fa["X"].attrs["shape"])
    print(f"context_{ctx} X shape:", shape)
    assert shape == (18400, 18533), shape
    fa.close()
chal_genes = read_lines(os.path.join(ROOT, "data/challenge/gene_names.csv"))
assert chal_genes[0] == "gene_name"
chal_genes = chal_genes[1:]
assert len(chal_genes) == 18533 and chal_genes[:3] == ["TSPAN6", "TNMD", "DPM1"]
print("challenge gene order OK:", chal_genes[:3])

# ---------------- 1. shared axis ----------------
bank_genes, bank_perts = {}, {}
for key, d in BANKS:
    bank_genes[key] = read_lines(os.path.join(d, "genes.txt"))
    bank_perts[key] = read_lines(os.path.join(d, "perts.txt"))
shared = sorted(set.intersection(*[set(v) for v in bank_genes.values()]))
print(f"7-way shared genes: {len(shared)} (sample {shared[:5]})")
assert len(shared) == 6346, len(shared)
SIDX = {g: i for i, g in enumerate(shared)}
chall_set = set(chal_genes)
print(f"shared genes in challenge axis: {sum(g in chall_set for g in shared)}"
      f"/{len(shared)}")

# ---------------- 2. streaming pass per bank ----------------
info = {}
for key, d in BANKS:
    D = np.load(os.path.join(d, "delta.npy"), mmap_mode="r")
    b = np.asarray(np.load(os.path.join(d, "b_ctrl.npy"), mmap_mode="r"))
    g = bank_genes[key]
    p = bank_perts[key]
    assert D.shape == (len(p), len(g)), (key, D.shape, len(p), len(g))
    assert b.shape == (len(g),)
    n_t, n_g = D.shape
    gpos = {gg: i for i, gg in enumerate(g)}
    scol = np.array([gpos[gg] for gg in shared], dtype=np.int64)
    b_sh = b[scol]
    block_pos = np.argsort(-b_sh, kind="stable")[:BLOCK]  # positions in shared
    block_full = scol[block_pos]
    in_block = np.zeros(n_g, dtype=bool)
    in_block[block_full] = True
    neg = np.zeros(len(shared), dtype=np.int64)
    nzero_sh = 0
    nzero_full = 0
    Db = np.empty((n_t, BLOCK), dtype=np.float32)
    Hb = np.empty((n_t, BLOCK), dtype=bool)
    for s in range(0, n_t, CHUNK):
        e = min(s + CHUNK, n_t)
        C = np.asarray(D[s:e])  # (c, G) transient
        Cs = C[:, scol]
        neg += (Cs < 0).sum(axis=0)
        nzero_sh += int((Cs == 0).sum())
        nzero_full += int((C == 0).sum())
        Db[s:e] = Cs[:, block_pos]
        Top = np.argpartition(-np.abs(C), HEAD - 1, axis=1)[:, :HEAD]
        # map full-axis head idx to block position via lookup table
        full2blk = np.full(n_g, -1, dtype=np.int64)
        full2blk[block_full] = np.arange(BLOCK)
        hb_idx = full2blk[Top]
        Hrows = np.zeros((e - s, BLOCK), dtype=bool)
        rr = np.broadcast_to(np.arange(e - s)[:, None], hb_idx.shape)
        valid = hb_idx >= 0
        Hrows[rr[valid], hb_idx[valid]] = True
        Hb[s:e] = Hrows
        del C, Cs, Top
    P = neg / n_t
    # own-target sign check (vectorised fancy read on mmap)
    rows, cols = [], []
    for ti, t in enumerate(p):
        j = gpos.get(t)
        if j is not None:
            rows.append(ti)
            cols.append(j)
    rows = np.array(rows, dtype=np.int64)
    cols = np.array(cols, dtype=np.int64)
    own = np.asarray(D[rows, cols]) if len(rows) else np.array([])
    cov, negfrac = len(rows), (float((own < 0).mean()) if len(rows) else np.nan)
    print(f"[{key}] {D.shape} P_down mean {P.mean():.4f} "
          f"zero_full {nzero_full/(n_t*n_g):.4f} zero_shared "
          f"{nzero_sh/(n_t*len(shared)):.4f} own-target {cov}/{n_t} "
          f"frac<0 {negfrac:.5f}", flush=True)
    info[key] = dict(P=P, Db=Db, Hb=Hb, n_t=n_t, n_g=n_g,
                     block_pos=block_pos, block_full=block_full,
                     bmin=float(b_sh[block_pos].min()),
                     zero_full=nzero_full / (n_t * n_g),
                     zero_shared=nzero_sh / (n_t * len(shared)),
                     own_cov=cov, own_neg=negfrac, n_targets=n_t)
    del D, b

# Jurkat own-target known answer from 370
assert info["nadig_jurkat"]["own_cov"] > 2000
assert abs(info["nadig_jurkat"]["own_neg"] - 0.99707) < 0.002, \
    info["nadig_jurkat"]["own_neg"]
print("jurkat own-target known-answer MATCH (~0.99707)")

# ---------------- 3. LOLO evaluation ----------------
keys = [k for k, _ in BANKS]
Pmat = np.stack([info[k]["P"] for k in keys])  # (7, S)
ntgt = np.array([info[k]["n_t"] for k in keys], dtype=np.float64)


def acc_pair(D, E, T, W, pred):
    ok = (pred == T) & E
    return float(ok.sum() / E.sum()), float((W * ok).sum() / (W * E).sum())


results = {}
for bi, L in enumerate(keys):
    DL = info[L]
    D2, H2 = DL["Db"], DL["Hb"]
    Z = (D2 == 0)
    E = ~(H2 | Z)
    T = (D2 < 0)
    W = np.abs(D2).astype(np.float64)
    n_eval = float(E.sum())
    others = [i for i in range(len(keys)) if i != bi]
    Po = Pmat[others]  # (6, S)
    Pv_eq = Po.mean(axis=0)
    w = ntgt[others] / ntgt[others].sum()
    Pv_w = (Po * w[:, None]).sum(axis=0)
    bpos = DL["block_pos"]
    ce = Pv_eq[bpos]
    cw = Pv_w[bpos]
    n_ties_eq = int((ce == 0.5).sum())
    out = {}
    base = float((T & E).sum() / E.sum())
    base_w = float((W * (T & E)).sum() / (W * E).sum())
    out["alldown"], out["alldown_w"] = base, base_w
    out["con_eq"], out["con_eq_w"] = acc_pair(D2, E, T, W, (ce > 0.5)[None, :])
    out["con_w"], out["con_w_w"] = acc_pair(D2, E, T, W, (cw > 0.5)[None, :])
    # (e) agree-only
    agree = ((Po > 0.5).all(axis=0) | (Po < 0.5).all(axis=0))
    ab = agree[bpos]
    pred_ag = np.where(ab[None, :], (ce > 0.5)[None, :], True)
    out["agree"], out["agree_w"] = acc_pair(D2, E, T, W, pred_ag)
    out["agree_entry_cov"] = float((ab[None, :] & np.ones_like(T)).sum()
                                   / T.size)
    out["agree_gene_cov"] = float(ab.mean())
    # (b) single-bank transfers
    for oi in others:
        Pv = Pmat[oi][bpos]
        a, aw = acc_pair(D2, E, T, W, (Pv > 0.5)[None, :])
        out[f"single_{keys[oi]}"] = a
        out[f"single_{keys[oi]}_w"] = aw
    # (f) shuffled consensus control, 5 draws
    rng_s = np.random.default_rng(5000 + bi)
    sh, sh_w = [], []
    for _ in range(5):
        Psh = Pv_eq[rng_s.permutation(len(shared))][bpos]
        a, aw = acc_pair(D2, E, T, W, (Psh > 0.5)[None, :])
        sh.append(a)
        sh_w.append(aw)
    out["shuf_mean"], out["shuf_sd"] = float(np.mean(sh)), \
        float(np.std(sh, ddof=1))
    out["shuf_w_mean"], out["shuf_w_sd"] = float(np.mean(sh_w)), \
        float(np.std(sh_w, ddof=1))
    out["n_eval"] = n_eval
    out["n_zero_dropped"] = float((Z & ~H2).sum())
    head_ov = H2.sum(axis=1)
    out["head_overlap_mean"] = float(head_ov.mean())
    out["head_overlap_min"] = float(head_ov.min())
    out["head_overlap_max"] = float(head_ov.max())
    out["zero_entry_frac"] = float((Z & ~H2).sum() / (BLOCK * DL["n_t"]
                                                      - H2.sum()))
    out["n_ties_eq"] = n_ties_eq
    # per-target vectors for the bootstrap (consensus-eq vs alldown)
    pred_c = ((ce > 0.5)[None, :] == T) & E
    pred_a = (T & E)
    out["_ct"] = pred_c.sum(axis=1).astype(np.float64)
    out["_at"] = pred_a.sum(axis=1).astype(np.float64)
    out["_nt"] = E.sum(axis=1).astype(np.float64)
    out["_wct"] = (W * pred_c).sum(axis=1)
    out["_wat"] = (W * pred_a).sum(axis=1)
    out["_wt"] = (W * E).sum(axis=1)
    # consensus-vs-donor sign agreement on shared axis (interpretation)
    out["_agree_with"] = {
        keys[oi]: float(np.mean((Pv_eq > 0.5) == (Pmat[oi] > 0.5)))
        for oi in others}
    results[L] = out
    sings = ", ".join(f"{keys[oi].split('_')[0][:4]}={out[f'single_{keys[oi]}']:.4f}"
                      for oi in others)
    print(f"[{L}] neval {n_eval:.0f} alldown {base:.4f}/{base_w:.4f} "
          f"con_eq {out['con_eq']:.4f}/{out['con_eq_w']:.4f} "
          f"con_w {out['con_w']:.4f} agree {out['agree']:.4f} "
          f"(cov {out['agree_entry_cov']:.3f}) shuf {out['shuf_mean']:.4f} "
          f"|| singles: {sings}", flush=True)

# ---------------- 4. bootstrap per bank (1000 draws) ----------------
for bi, L in enumerate(keys):
    r = results[L]
    n = info[L]["n_t"]
    rng = np.random.default_rng(20260910 + bi * 7919)
    IX = rng.integers(0, n, size=(1000, n))
    CT, AT, NT = r["_ct"], r["_at"], r["_nt"]
    WCT, WAT, WT = r["_wct"], r["_wat"], r["_wt"]
    du = CT[IX].sum(1) / NT[IX].sum(1) - AT[IX].sum(1) / NT[IX].sum(1)
    dw = WCT[IX].sum(1) / WT[IX].sum(1) - WAT[IX].sum(1) / WT[IX].sum(1)
    r["boot_du_mean"] = float(du.mean())
    r["boot_du_lo"], r["boot_du_hi"] = float(np.percentile(du, 2.5)), \
        float(np.percentile(du, 97.5))
    r["boot_du_frac_pos"] = float((du > 0).mean())
    r["boot_dw_mean"] = float(dw.mean())
    r["boot_dw_lo"], r["boot_dw_hi"] = float(np.percentile(dw, 2.5)), \
        float(np.percentile(dw, 97.5))
    r["boot_dw_frac_pos"] = float((dw > 0).mean())
    print(f"[boot {L}] d_unw {r['boot_du_mean']:+.4f} "
          f"[{r['boot_du_lo']:+.4f},{r['boot_du_hi']:+.4f}] "
          f"P>0 {r['boot_du_frac_pos']:.3f} | d_w {r['boot_dw_mean']:+.4f} "
          f"[{r['boot_dw_lo']:+.4f},{r['boot_dw_hi']:+.4f}]", flush=True)

# ---------------- 5. A/B/C CPM blocks + coverage + A sanity ----------------
print("== A/B/C mean CPM ==", flush=True)
cpm = {}
for ctx in ["A", "B", "C"]:
    fa = h5py.File(os.path.join(ROOT, f"data/challenge/context_{ctx}.h5ad"), "r")
    X, n_cells, n_genes = fa["X"], 18400, 18533
    indptr = X["indptr"][:]
    indices = X["indices"]
    data = X["data"]
    gs = np.zeros(n_genes, dtype=np.float64)
    lib = np.zeros(n_cells, dtype=np.float64)
    for s in range(0, n_cells, CHUNK * 2):
        e = min(s + CHUNK * 2, n_cells)
        ip0, ip1 = int(indptr[s]), int(indptr[e])
        dd = data[ip0:ip1].astype(np.float64)
        ix = indices[ip0:ip1].astype(np.int64)
        ipl = indptr[s:e + 1].astype(np.int64) - ip0
        rep = np.repeat(np.arange(e - s), np.diff(ipl))
        lib[s:e] = np.bincount(rep, weights=dd, minlength=e - s)
        gs += np.bincount(ix, weights=dd / lib[s:e][rep],
                          minlength=n_genes) * 1e6
    cpm[ctx] = gs / n_cells
    med_lib = float(np.median(lib))
    print(f"context_{ctx} median lib {med_lib:.0f} "
          f"CPM max {cpm[ctx].max():.0f} ({chal_genes[int(np.argmax(cpm[ctx]))]})",
          flush=True)
    if ctx == "A":
        assert med_lib == 20109.0, med_lib
        print("context_A median-lib known-answer MATCH (20109)")
    fa.close()
    np.save(os.path.join(WORK, f"context{ctx}_mean_cpm.npy"), cpm[ctx])
print("(median-lib known answer checked on the lines above; record 20109)")

con_genes_set = set(shared)
cover = {}
for ctx in ["A", "B", "C"]:
    order = np.argsort(-cpm[ctx], kind="stable")[:2200]
    blk = [chal_genes[i] for i in order]
    cover[ctx] = {"top2200_in_consensus": int(sum(g in con_genes_set
                                                  for g in blk)),
                  "min_cpm": float(cpm[ctx][order].min())}
    print(f"context_{ctx} top-2200: {cover[ctx]['top2200_in_consensus']}/2200 "
          f"in consensus axis; min CPM {cover[ctx]['min_cpm']:.1f}")

# A sanity on 370's exact A-CPM block (Jurkat-mapped), scored with Jurkat deltas
print("== A sanity (370 A-block, Jurkat deltas) ==", flush=True)
jg = bank_genes["nadig_jurkat"]
jpos = {g: i for i, g in enumerate(jg)}
mapped = np.array([g in jpos for g in chal_genes])
orderA = np.argsort(-cpm["A"], kind="stable")
rankA = np.empty(len(chal_genes), dtype=np.int64)
rankA[orderA] = np.arange(len(chal_genes))
cut = np.sort(rankA[mapped])[2199]
blockA_chal = np.where(mapped & (rankA <= cut))[0]
assert len(blockA_chal) == 2200
print(f"370 A-block cutoff reproduced: size {len(blockA_chal)}, "
      f"min CPM {cpm['A'][blockA_chal].min():.1f} (record 130.4)")
blk_genes = [chal_genes[i] for i in blockA_chal]
Peq_all = Pmat.mean(axis=0)
bi_j = keys.index("nadig_jurkat")
Peq_lolo = Pmat[[i for i in range(len(keys)) if i != bi_j]].mean(axis=0)
Pmap_all = np.array([Peq_all[SIDX[g]] if g in SIDX else np.nan
                     for g in blk_genes])
Pmap_lolo = np.array([Peq_lolo[SIDX[g]] if g in SIDX else np.nan
                      for g in blk_genes])
Dj = info["nadig_jurkat"]["Db"]  # NOTE: Db is on LOLO-block cols, not A-block
# rebuild Jurkat deltas on the A-block columns directly from mmap
Djm = np.load(os.path.join(dict(BANKS)["nadig_jurkat"], "delta.npy"),
              mmap_mode="r")
jcols = np.array([jpos[g] for g in blk_genes])
Da = np.asarray(Djm[:, jcols])
TopA = np.argpartition(-np.abs(np.asarray(Djm)), HEAD - 1, axis=1)[:, :HEAD]
full2blkA = np.full(len(jg), -1, dtype=np.int64)
full2blkA[jcols] = np.arange(len(blk_genes))
hbA = full2blkA[TopA]
Ha = np.zeros((Da.shape[0], len(blk_genes)), dtype=bool)
rrA = np.broadcast_to(np.arange(Da.shape[0])[:, None], hbA.shape)
validA = hbA >= 0
Ha[rrA[validA], hbA[validA]] = True
Za = (Da == 0)
Ea = ~(Ha | Za)
Ta = (Da < 0)
Wa = np.abs(Da).astype(np.float64)
sanity = {"ablock_size": 2200,
          "ablock_covered_by_consensus": int(np.sum(~np.isnan(Pmap_all)))}
for nm, Pv in (("full_consensus", Pmap_all), ("lolo_consensus", Pmap_lolo)):
    m = ~np.isnan(Pv)
    pred = np.where(m[None, :], (Pv > 0.5)[None, :], True)
    E2 = Ea & m[None, :]
    # covered_only: eval on covered entries; plus_alldown: uncovered->down
    a_cov = float((((pred == Ta) & E2).sum() / E2.sum()))
    aw_cov = float(((Wa * ((pred == Ta) & E2)).sum() / (Wa * E2).sum()))
    a_all = float((((pred == Ta) & Ea).sum() / Ea.sum()))
    base = float(((Ta & Ea).sum() / Ea.sum()))
    sanity[f"{nm}_covered_acc"] = a_cov
    sanity[f"{nm}_covered_acc_w"] = aw_cov
    sanity[f"{nm}_with_uncovered_as_down_acc"] = a_all
    sanity[f"{nm}_alldown_same_entries"] = base
    print(f"[A sanity {nm}] covered {m.sum()}/2200 acc {a_cov:.4f} "
          f"(w {aw_cov:.4f}); with-uncovered-as-down {a_all:.4f} vs "
          f"alldown {base:.4f} (370 same-line 0.5766)", flush=True)

# ---------------- 6. outputs ----------------
Peq = Pmat.mean(axis=0)
Pwt = (Pmat * (ntgt / ntgt.sum())[:, None]).sum(axis=0)
np.savez(os.path.join(WORK, "consensus_prior.npz"),
         genes=np.array(shared), P_down=Peq, P_down_weighted=Pwt,
         P_banks=Pmat, bank_names=np.array(keys),
         n_targets=np.array([info[k]["n_t"] for k in keys]),
         n_banks=np.int64(len(keys)), seed=np.int64(20260910))
print("wrote consensus_prior.npz")

def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items() if not k.startswith("_")}
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    return o

metrics = {
    "per_bank": clean(results),
    "pooled": {},
    "coverage": cover,
    "a_sanity": sanity,
    "bank_info": {k: {kk: v for kk, v in info[k].items()
                      if kk in ("n_t", "n_g", "bmin", "zero_full",
                                "zero_shared", "own_cov", "own_neg")}
                  for k in keys},
    "shared_axis": {"n": len(shared), "head": shared[:5], "tail": shared[-5:]},
}
# pooled (micro) over held-out evals
for m in ["alldown", "alldown_w", "con_eq", "con_eq_w", "con_w", "con_w_w",
          "agree", "agree_w"]:
    num = sum(results[L][m] * results[L]["n_eval"] for L in keys)
    den = sum(results[L]["n_eval"] for L in keys)
    metrics["pooled"][m] = num / den
metrics["pooled"]["macro_con_eq_minus_alldown"] = float(np.mean(
    [results[L]["con_eq"] - results[L]["alldown"] for L in keys]))
metrics["pooled"]["median_con_eq_minus_alldown"] = float(np.median(
    [results[L]["con_eq"] - results[L]["alldown"] for L in keys]))
json.dump(metrics, open(os.path.join(WORK, "metrics371.json"), "w"), indent=1)

with open(os.path.join(WORK, "per-bank371.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["heldout", "estimator", "acc_unw", "acc_w", "extra"])
    for L in keys:
        r = results[L]
        w.writerow([L, "alldown", r["alldown"], r["alldown_w"], ""])
        w.writerow([L, "con_eq", r["con_eq"], r["con_eq_w"],
                    f"ties={r['n_ties_eq']}"])
        w.writerow([L, "con_w", r["con_w"], r["con_w_w"], ""])
        w.writerow([L, "agree", r["agree"], r["agree_w"],
                    f"entrycov={r['agree_entry_cov']:.4f}"])
        w.writerow([L, "shuffled", r["shuf_mean"], r["shuf_w_mean"],
                    f"sd={r['shuf_sd']:.4f}"])
        for oi in range(len(keys)):
            if keys[oi] == L:
                continue
            w.writerow([L, f"single_{keys[oi]}", r[f"single_{keys[oi]}"],
                        r[f"single_{keys[oi]}_w"], ""])
        w.writerow([L, "boot_con_minus_alldown_unw", r["boot_du_mean"], "",
                    f"CI=[{r['boot_du_lo']:.4f},{r['boot_du_hi']:.4f}] "
                    f"P>0={r['boot_du_frac_pos']:.3f}"])
        w.writerow([L, "boot_con_minus_alldown_w", r["boot_dw_mean"], "",
                    f"CI=[{r['boot_dw_lo']:.4f},{r['boot_dw_hi']:.4f}] "
                    f"P>0={r['boot_dw_frac_pos']:.3f}"])
print("done")
