#!/usr/bin/env python3
"""403 -- neighbourhood denoising of K562-GWPS donor rows: does shrinkage toward
the kNN mean help BULK-channel endpoints (retrieval / own cosine) cross-line?

Frame: source = K562-GWPS rows; recipients = nadig_jurkat, nadig_hepg2, rpe1
(cross-line), k562_essential (same-line upper bound). Truth = recipient's own
delta rows. Shared gene axis per pair; queries = targets in both banks; the
query target's own gene column is dropped from every vector (exact per-query
correction). Neighbours = top-K GWPS rows by Pearson on FULL GWPS rows, self
excluded (398's convention). Retrieval uses the pds form (brief S1) with the
full recipient bank as gallery. Paired 10k-target bootstrap, arm minus raw.
"""
import json, os, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np
from pathlib import Path
import resource, threading
_stop = threading.Event()
def _sampler():
    peak = 0
    while not _stop.wait(1.0):
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if rss > peak:
            peak = rss
            print(f"[rss] {time.strftime('%H:%M:%S')} maxrss {peak/2**30:.2f} GiB", flush=True)
threading.Thread(target=_sampler, daemon=True).start()

REPO = Path("~/Code/virtual-cell-challenge-2026")
OUT = REPO / "experiments/403-bulk-knn-denoise/work"
OUT.mkdir(parents=True, exist_ok=True)
DON = REPO / "data/donors"
HEPG2 = REPO / "experiments/166-donor-acquisition/donors/nadig_hepg2"
RECIPS = {"nadig_jurkat": DON / "nadig_jurkat", "nadig_hepg2": HEPG2,
          "rpe1": DON / "rpe1", "k562_essential": DON / "k562_essential"}
KS = (10, 25, 50)
WS = (0.25, 0.50, 0.75)
NBOOT = 10_000


def load_bank(d):
    d = Path(d)
    D = np.load(d / "delta.npy", mmap_mode="r")
    P = [l.strip() for l in open(d / "perts.txt")]
    G = [l.strip() for l in open(d / "genes.txt")]
    assert D.shape == (len(P), len(G)), (d, D.shape, len(P), len(G))
    return D, P, G


t0 = time.time()
GW, GWP, GWG = load_bank(DON / "k562_gwps")
gwp = {p: i for i, p in enumerate(GWP)}
gwg = {g: i for i, g in enumerate(GWG)}
nGW = len(GWP)
# row means / norms of full GWPS rows (chunked, for Pearson sims)
gw_mean = np.zeros(nGW, dtype=np.float64)
gw_norm = np.zeros(nGW, dtype=np.float64)
for a in range(0, nGW, 1000):
    C = np.asarray(GW[a:a + 1000], dtype=np.float32)
    gw_mean[a:a + 1000] = C.mean(1)
    Cc = C - gw_mean[a:a + 1000, None]
    gw_norm[a:a + 1000] = np.sqrt((Cc.astype(np.float64) ** 2).sum(1))
gw_norm = np.maximum(gw_norm, 1e-12)
grand_mean_full = None  # set per shared axis below
print(f"GWPS {GW.shape}, rows mapped {nGW}", flush=True)

rng_rand = np.random.default_rng(403)
rng_boot = np.random.default_rng(1403)
results, per_target = {}, {}


def score_arms(Rfull, T, J, q_r, E_dict):
    """Rfull (N,G) gallery, T (nT,G) truth, J (nT,) drop cols (-1 = absent),
    E_dict arm -> (nT,G) estimates. Returns dict arm -> dict of endpoints and
    per-target arrays."""
    N, G = Rfull.shape
    nT = T.shape[0]
    gn2 = (Rfull.astype(np.float64) ** 2).sum(1)  # (N,)
    RfullJ = Rfull[:, np.maximum(J, 0)]  # (N,nT); col 0 junk where J==-1
    hasJ = (J >= 0)
    RfullJ = RfullJ * hasJ[None, :]
    tn2 = (T.astype(np.float64) ** 2).sum(1)
    Tqj = np.where(hasJ, T[np.arange(nT), np.maximum(J, 0)], 0.0)
    out = {}
    for arm, E in E_dict.items():
        E = np.asarray(E, dtype=np.float32)
        DOTS = Rfull.astype(np.float32) @ E.T.astype(np.float32)  # (N,nT) f32
        DOTS = DOTS.astype(np.float64)
        Eqj = np.where(hasJ, E[np.arange(nT), np.maximum(J, 0)].astype(np.float64), 0.0)
        DOTS -= RfullJ.astype(np.float64) * Eqj[None, :]
        en2 = (E.astype(np.float64) ** 2).sum(1)
        en_d = np.sqrt(np.maximum(en2 - Eqj ** 2, 0.0))
        en_d = np.maximum(en_d, 1e-300)
        gn_d = np.sqrt(np.maximum(gn2[:, None] - RfullJ.astype(np.float64) ** 2, 0.0))
        with np.errstate(divide="ignore", invalid="ignore"):
            COS = np.where(gn_d > 0, DOTS / (en_d[None, :] * np.maximum(gn_d, 1e-300)), 0.0)
        own = (E.astype(np.float64) * T.astype(np.float64)).sum(1) - Eqj * Tqj
        tn_d = np.sqrt(np.maximum(tn2 - Tqj ** 2, 0.0))
        owncos = np.where(tn_d > 0, own / (en_d * np.maximum(tn_d, 1e-300)), 0.0)
        # retrieval, pds form with tie averaging; pin the self-match exact so
        # GEMM-vs-dot summation order cannot jitter the rank by half a step
        COS[q_r, np.arange(nT)] = owncos
        greater = (COS > owncos[None, :]).sum(0).astype(np.float64)
        equal = (COS == owncos[None, :]).sum(0).astype(np.float64)
        k = greater + 0.5 * (equal - 1.0)
        retr = 1.0 - k / (N - 1)
        num = float(own.sum())
        den = float(np.sqrt(max((en2 - Eqj ** 2).sum(), 0.0) * max((tn2 - Tqj ** 2).sum(), 0.0)))
        cosagg = num / den if den > 0 else 0.0
        out[arm] = {"retrieval": float(retr.mean()), "owncos": float(owncos.mean()),
                    "cosagg": cosagg, "retr_pt": retr, "cos_pt": owncos}
        print(f"    {arm:16s} retr {retr.mean():.5f} owncos {owncos.mean():.5f} cosagg {cosagg:.5f}", flush=True)
    return out


for rname, rdir in RECIPS.items():
    D, P, G = load_bank(rdir)
    gi = {g: i for i, g in enumerate(G)}
    shared = [g for g in G if g in gwg]  # recipient order
    idx_r = np.array([gi[g] for g in shared])
    idx_gw = np.array([gwg[g] for g in shared])
    qnames = [t for t in P if t in gwp]
    q_gw = np.array([gwp[t] for t in qnames])
    pi = {p: i for i, p in enumerate(P)}
    q_r = np.array([pi[t] for t in qnames])
    assert len(set(qnames)) == len(qnames) == len(set(q_gw)), "dup query targets"
    assert len(set(P)) == len(P), "dup recipient perts"
    nT = len(qnames)
    print(f"{rname}: {D.shape}, shared genes {len(shared)}, queries {nT}/{len(P)}", flush=True)
    # Pearson sims of query rows vs all GWPS rows (full axis), chunked
    Q = np.asarray(GW[q_gw], dtype=np.float32)  # (nT, 8248)
    qm = Q.mean(1, keepdims=True)
    Qc = Q - qm
    qn = np.sqrt((Qc.astype(np.float64) ** 2).sum(1))
    qn = np.maximum(qn, 1e-12)
    Qn = (Qc / qn[:, None]).astype(np.float32)
    sims = np.empty((nGW, nT), dtype=np.float32)
    for a in range(0, nGW, 1000):
        C = np.asarray(GW[a:a + 1000], dtype=np.float32)
        Cc = (C - gw_mean[a:a + 1000, None].astype(np.float32)).astype(np.float32)
        sims[a:a + 1000] = (Cc / gw_norm[a:a + 1000, None].astype(np.float32)) @ Qn.T
    # neighbour QC (self excluded: each query's own row reads exactly 1.0 in sims)
    s_no = sims.copy()
    s_no[q_gw, np.arange(nT)] = -np.inf
    s_sorted = np.sort(s_no, axis=0)
    med_top1 = float(np.median(s_sorted[-1]))
    med_top25 = float(np.median(s_sorted[-25]))
    print(f"  median top1 pearson {med_top1:.4f}, top25 {med_top25:.4f}", flush=True)
    nb = {}
    for K in KS:
        NB = np.empty((nT, K), dtype=np.int64)
        for j in range(nT):
            s_ = sims[:, j].copy()
            s_[q_gw[j]] = -np.inf
            NB[j] = np.argsort(-s_, kind="stable")[:K]
        nb[K] = NB
    # shared-axis source block (RAM) + recipient block
    GWs = np.asarray(GW[:, idx_gw], dtype=np.float32)  # (9787, Gs)
    Rfull = np.asarray(D[:, idx_r], dtype=np.float32)  # (N, Gs)
    T = Rfull[q_r]
    gpos = {g: i for i, g in enumerate(shared)}
    J = np.array([gpos[t] if t in gpos else -1 for t in qnames])
    print(f"  queries with target on shared axis: {(J >= 0).sum()}/{nT}", flush=True)
    grand = GWs.mean(0)

    def nb_mean(NB):
        # lean row-group mean: never materialise the (nT, K, Gs) gather
        # (2252 x 50 x 7352 f32 = 3.3 GiB transient -- blew the 8 GiB cap)
        M = np.empty((nT, GWs.shape[1]), dtype=np.float32)
        for j in range(nT):
            M[j] = GWs[NB[j]].mean(0)
        return M

    E = {"raw": GWs[q_gw]}
    for K in KS:
        M = nb_mean(nb[K])
        E[f"mean{K}"] = M
        for w in WS:
            E[f"shrink{K}_w{w:g}"] = (1 - w) * GWs[q_gw] + w * M
        del M
    for w in WS:
        E[f"grand_w{w:g}"] = (1 - w) * GWs[q_gw] + w * grand[None, :]
    for s in range(3):
        RR = np.empty((nT, 25), dtype=np.int64)
        for j in range(nT):
            pool = np.delete(np.arange(nGW), q_gw[j])
            RR[j] = rng_rand.choice(pool, 25, replace=False)
        E[f"random25_s{s}"] = nb_mean(RR)
        del RR
    scored = score_arms(Rfull, T, J, q_r, E)
    # bootstrap: every arm minus raw
    base_r, base_c = scored["raw"]["retr_pt"], scored["raw"]["cos_pt"]
    B = rng_boot.integers(0, nT, size=(NBOOT, nT))
    arms_out, pt = {}, {}
    for arm, v in scored.items():
        dr = v["retr_pt"] - base_r
        dc = v["cos_pt"] - base_c
        br = dr[B].mean(1, dtype=np.float64)
        bc = dc[B].mean(1, dtype=np.float64)
        arms_out[arm] = {"retrieval": round(v["retrieval"], 6), "owncos": round(v["owncos"], 6),
                         "cosagg": round(v["cosagg"], 6),
                         "d_retr": round(float(dr.mean()), 6),
                         "d_retr_ci": [round(float(np.quantile(br, .025)), 6), round(float(np.quantile(br, .975)), 6)],
                         "d_cos": round(float(dc.mean()), 6),
                         "d_cos_ci": [round(float(np.quantile(bc, .025)), 6), round(float(np.quantile(bc, .975)), 6)]}
        pt[arm + "_retr"] = v["retr_pt"].astype(np.float32)
        pt[arm + "_cos"] = v["cos_pt"].astype(np.float32)
    results[rname] = {"n_bank": int(D.shape[0]), "n_shared_genes": len(shared),
                      "n_queries": nT, "n_drop": int((J >= 0).sum()),
                      "med_top1_pearson": round(med_top1, 4),
                      "med_top25_pearson": round(med_top25, 4),
                      "arms": arms_out}
    per_target[rname] = pt
    del GWs, Rfull, T, E, sims, scored, Q, Qc, Qn, s_no, s_sorted, nb, B, grand
    print(f"  done {rname} ({time.time()-t0:.0f}s elapsed)", flush=True)

_stop.set()
(OUT / "results403.json").write_text(json.dumps(results, indent=1))
np.savez_compressed(OUT / "per_target403.npz",
                    **{f"{r}__{a}": v for r, d in per_target.items() for a, v in d.items()})
print(f"TOTAL {time.time()-t0:.0f}s", flush=True)
