#!/usr/bin/env python3
"""406 -- where the reach walk breaks + neighbourhood-aware ORDER keys within the head (H1 analytic).

Frame + walk VERBATIM from experiments/397-knn-tilt/knn397.py (exec'd prefix); shipped-sign
construction VERBATIM from knn397d.py (six-bank family, K=25, shrink 400). Signs fixed at SG399
for every variant, so fid/nmae/jac are unchanged by construction; only the order moves.
Orion banks stay mmap'd and are streamed in chunks (bit-identical int64 vote counts)."""
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "~/Code/virtual-cell-challenge-2026/experiments/406-reach-front/work/numba-cache")
os.environ.setdefault("NUMBA_NUM_THREADS", "4"); os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np, json, time
t0 = time.time()

# ---- frame VERBATIM (knn397.py above "# ---- kNN tilts"; KA print suppressed) ----
exec(open("~/Code/virtual-cell-challenge-2026/experiments/397-knn-tilt/knn397.py").read().split("# ---- kNN tilts")[0].replace('print("KNOWN ANSWER', 'pass  # print("KNOWN ANSWER'))
OUT406 = REPO / "experiments/406-reach-front/work"; OUT406.mkdir(parents=True, exist_ok=True)
rng406 = np.random.default_rng(406)

# ---- 1. known answers ----
reach_E0 = float(np.nanmean(M0["reach"]))
n_scored = int((n_conf > 0).sum())
bp_E0 = block_prec(SG0)
print(f"KNOWN reach(E0,delivered)={reach_E0:.8f} (expect 0.25435; scored {n_scored}/85, expect 82)", flush=True)
print(f"KNOWN block_precision(E0)={bp_E0:.4f} (expect 0.5083)", flush=True)
assert abs(reach_E0 - 0.25435) < 5e-4, "reach known answer mismatch"
assert abs(bp_E0 - 0.5083) < 1e-4, "block precision known answer mismatch"

# ---- 2. shipped construction (knn397d.py VERBATIM logic) + chunked voter ----
for b in ("orion_hct116", "orion_hek293t", "nadig_hepg2"):
    d = (REPO / "data/donors" / b) if b != "nadig_hepg2" else (REPO / "experiments/166-donor-acquisition/donors/nadig_hepg2")
    D = np.load(d / "delta.npy", mmap_mode="r"); P = np.array([l.strip() for l in open(d / "perts.txt")]); G = np.array([l.strip() for l in open(d / "genes.txt")])
    pi = {p: i for i, p in enumerate(P)}; gi_b = {g: i for i, g in enumerate(G)}; extra[b] = (D, pi, np.array([gi_b.get(g, -1) for g in genes]))
FAM = ["k562_essential", "nadig_jurkat", "nadig_hepg2", "orion_hct116", "orion_hek293t"]

def votes_chunked(rows_idx, banks, chunk=500):
    dn = np.zeros(nG, dtype=np.int64); nn = np.zeros(nG, dtype=np.int64)
    sub = GW[rows_idx][:, np.where(covg, gmap, 0)]; nzm = sub != 0
    dn[covg] += ((sub < 0) & nzm).sum(0)[covg]; nn[covg] += nzm.sum(0)[covg]
    for b in banks:
        D, pi, gm = extra[b]; ridx = [pi[GWP[j]] for j in rows_idx if GWP[j] in pi]
        if not ridx: continue
        cvb = gm >= 0; col = np.where(cvb, gm, 0)
        for s in range(0, len(ridx), chunk):
            subb = D[ridx[s:s + chunk]][:, col]; nzb = subb != 0
            dn[cvb] += ((subb < 0) & nzb).sum(0)[cvb]; nn[cvb] += nzb.sum(0)[cvb]
    return np.where(nn > 0, dn / np.maximum(nn, 1), np.nan), nn

# control: chunked voter reproduces verbatim votes() EXACTLY
Pv, nnv = votes(np.arange(len(GWP)), ["k562_essential", "nadig_jurkat"])
Pc, nnc = votes_chunked(np.arange(len(GWP)), ["k562_essential", "nadig_jurkat"])
assert np.array_equal(nnv, nnc) and np.array_equal(Pv, Pc, equal_nan=True), "chunked voter mismatch (small banks)"
sub_idx = np.sort(rng406.choice(len(GWP), 200, replace=False))
for b in ("orion_hct116", "orion_hek293t"):
    Pv, nnv = votes(sub_idx, [b]); Pc, nnc = votes_chunked(sub_idx, [b])
    assert np.array_equal(nnv, nnc) and np.array_equal(Pv, Pc, equal_nan=True), f"chunked voter mismatch ({b})"
print("chunked voter: EXACT match on small banks (full) + Orion (200-row subset)", flush=True)

P_marg6, _ = votes_chunked(np.arange(len(GWP)), FAM); sg6 = signs_from(P_marg6, sg4)
K = 25; m = 400.0
has = np.array([t in gwp for t in targets])
SG399 = SG0.copy(); PS399 = np.full((nT, nG), np.nan); SGK399 = np.full((nT, nG), -1.0)
for i in range(nT):
    if not has[i]:
        SG399[i] = np.where(cb[i], sg6, SG0[i]); PS399[i] = P_marg6; SGK399[i] = sg6; continue
    v = GWn[gwp[t := targets[i]]]; s_ = GWn @ v; s_[gwp[t]] = -np.inf; nb = np.argsort(-s_)[:K]
    P, nn = votes_chunked(nb, FAM)
    Ps = np.where(np.isfinite(P) & np.isfinite(P_marg6), (nn * np.nan_to_num(P) + m * np.nan_to_num(P_marg6)) / (nn + m), np.where(np.isfinite(P), P, P_marg6))
    PS399[i] = Ps; sk = signs_from(Ps, sg6); SGK399[i] = sk
    SG399[i] = np.where(cb[i], sk, SG0[i])
bp_399 = block_prec(SG399)
R399_deliv = walk_reach(cur, adjud & (SG399 == np.sign(lrm)) & pool, adjud & pool, n_conf)
print(f"SHIPPED-399 block_precision={bp_399:.4f} (prereg ~0.5338); reach(SG399,delivered)={np.nanmean(R399_deliv):.5f} (prereg ~0.2519)", flush=True)

truth = np.sign(lrm)
conf399 = np.where(np.isfinite(PS399), np.abs(PS399 - 0.5), 0.0)
agree399 = (SGK399 == SG0)
match399 = adjud & (SG399 == truth) & pool
adjpool = adjud & pool

# ---- 3. descriptive: break rank per target (SG399 signs, delivered order) ----
scored = n_conf > 0
desc = []
for i in range(nT):
    if not scored[i]: continue
    o = cur[i]; am = adjpool[i, o]; adjord = o[am]; L = adjord.size
    corr = match399[i, adjord].astype(float)
    pur = np.cumsum(corr) / np.arange(1, L + 1)
    okk = np.flatnonzero(pur >= 0.9); ks = int(okk.max() + 1) if okk.size else 0
    R = ks / n_conf[i]
    assert abs(R - R399_deliv[i]) < 1e-12, f"walk mismatch t={targets[i]}"
    if ks < L:
        g = int(adjord[ks])
        cls = "head" if head[i, g] else ("block" if block[i, g] else "uncalled")
        desc.append(dict(t=targets[i], n_real=int(n_conf[i]), L=int(L), ks=ks, R=R, broke=True,
                         b=ks + 1, g=genes[g], cls=cls, E1=float(E1[i, g]), E0=float(E0[i, g]),
                         conf=float(conf399[i, g]), agree=bool(agree399[i, g]),
                         right=bool(SG399[i, g] == truth[i, g]),
                         P10=float(pur[9]) if L >= 10 else np.nan,
                         P20=float(pur[19]) if L >= 20 else np.nan,
                         P50=float(pur[49]) if L >= 50 else np.nan))
    else:
        desc.append(dict(t=targets[i], n_real=int(n_conf[i]), L=int(L), ks=ks, R=R, broke=False,
                         b=None, g=None, cls=None, E1=np.nan, E0=np.nan, conf=np.nan, agree=None,
                         right=None,
                         P10=float(pur[9]) if L >= 10 else np.nan,
                         P20=float(pur[19]) if L >= 20 else np.nan,
                         P50=float(pur[49]) if L >= 50 else np.nan))
D = desc
q1, q2 = np.quantile([d["n_real"] for d in D], [1 / 3, 2 / 3])
for d in D: d["terc"] = int(np.searchsorted([q1 + 1e-9, q2 + 1e-9], d["n_real"]))
print(f"n_real tercile edges: {q1:.0f}/{q2:.0f}; mean k*={np.mean([d['ks'] for d in D]):.1f} (388: ~14)", flush=True)

def Across(dd):
    br = [d for d in dd if d["broke"]]
    return dict(n=len(dd), n_real=np.mean([d["n_real"] for d in dd]),
                ks=np.mean([d["ks"] for d in dd]), R=np.mean([d["R"] for d in dd]),
                never=1 - len(br) / len(dd),
                b_mean=np.mean([d["b"] for d in br]) if br else np.nan,
                b_med=np.median([d["b"] for d in br]) if br else np.nan,
                head=np.mean([d["cls"] == "head" for d in br]) if br else np.nan,
                block=np.mean([d["cls"] == "block" for d in br]) if br else np.nan,
                uncalled=np.mean([d["cls"] == "uncalled" for d in br]) if br else np.nan,
                P10=np.nanmean([d["P10"] for d in dd]), P20=np.nanmean([d["P20"] for d in dd]),
                P50=np.nanmean([d["P50"] for d in dd]),
                conf=np.nanmean([d["conf"] for d in br]) if br else np.nan,
                agree=np.nanmean([float(d["agree"]) for d in br]) if br else np.nan,
                right=np.nanmean([float(d["right"]) for d in br]) if br else np.nan)
AGG = {"all": Across(D), **{f"terc{t}": Across([d for d in D if d["terc"] == t]) for t in (0, 1, 2)}}
for k, v in AGG.items(): print(k, json.dumps({a: (round(float(b), 4) if b == b else None) for a, b in v.items()}), flush=True)

# ---- 4. order variants (signs SG399 fixed) ----
def reorder_head(key=None, demote=None):
    out = cur.copy()
    for i in range(nT):
        o = cur[i]; hpos = np.flatnonzero(head[i, o]); hg = o[hpos]
        if demote is not None:
            dm = demote[i, hg]; new = np.concatenate([hg[~dm], hg[dm]])
        else:
            k = np.nan_to_num(-key[i, hg].astype(float), nan=np.inf)
            new = hg[np.lexsort((np.arange(hg.size), k))]
        out[i, hpos] = new
    return out

A = np.abs(E1)
AG = agree399.astype(float)
orders = {"a_delivered": cur}
for beta in (0.5, 1.0, 2.0):
    orders[f"b_b{beta}_agree"] = reorder_head(key=A * (1 + beta * AG))
    orders[f"c_b{beta}_conf"] = reorder_head(key=A * (1 + beta * conf399 * AG - beta * conf399 * (1 - AG)))
orders["d_demote_dis02"] = reorder_head(demote=(head & ~agree399 & (conf399 >= 0.2)))
dconf = conf399.copy(); dagree = agree399.copy()
for i in range(nT):
    hg = np.flatnonzero(head[i])
    dconf[i, hg] = dconf[i, rng406.permutation(hg)]
    dagree[i, hg] = dagree[i, rng406.permutation(hg)]
for beta in (0.5, 1.0, 2.0):
    orders[f"e_derconf_c_b{beta}"] = reorder_head(key=A * (1 + beta * dconf * AG - beta * dconf * (1 - AG)))
orders["e_derconf_d"] = reorder_head(demote=(head & ~agree399 & (dconf >= 0.2)))
orders["e_deragree_b1"] = reorder_head(key=A * (1 + 1.0 * dagree.astype(float)))

Ra = R399_deliv[scored]
res = {}
for name, co in orders.items():
    Rv = walk_reach(co, match399, adjpool, n_conf)[scored]
    dv = Rv - Ra
    B = rng406.integers(0, dv.size, (10000, dv.size)); bs = dv[B].mean(1)
    # break-rank distribution under this order
    bs_list = []
    for ii, i in enumerate(np.flatnonzero(scored)):
        o = co[i]; am = adjpool[i, o]; ao = o[am]; L = ao.size
        pur = np.cumsum(match399[i, ao].astype(float)) / np.arange(1, L + 1)
        okk = np.flatnonzero(pur >= 0.9); ks = int(okk.max() + 1) if okk.size else 0
        bs_list.append(None if ks >= L else ks + 1)
    brk = [b for b in bs_list if b is not None]
    res[name] = dict(reach=float(np.nanmean(Rv)), dmean=float(np.mean(dv)),
                     ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
                     b_mean=float(np.mean(brk)), b_med=float(np.median(brk)),
                     never=float(1 - len(brk) / len(bs_list)))
    print(f"{name}: reach {res[name]['reach']:.5f}  Δ {res[name]['dmean']:+.5f} [{res[name]['ci'][0]:+.5f},{res[name]['ci'][1]:+.5f}]  b_mean {res[name]['b_mean']:.1f} b_med {res[name]['b_med']:.0f} never {res[name]['never']:.2f}", flush=True)

(OUT406 / "results406.json").write_text(json.dumps(dict(
    known=dict(reach_E0=reach_E0, block_prec_E0=bp_E0, n_scored=n_scored),
    shipped=dict(block_prec_399=bp_399, reach_399_deliv=float(np.nanmean(R399_deliv))),
    descriptive_agg=AGG, variants=res), indent=1))
import pickle
pickle.dump(desc, open(OUT406 / "descriptive.pkl", "wb"))
print(f"done in {time.time()-t0:.0f}s", flush=True)
