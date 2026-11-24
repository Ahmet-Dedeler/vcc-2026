"""344 -- empirical-null key kill check on the challenge controls. See REGISTRATION.md. Truth-free."""
from __future__ import annotations
import json, sys, time, gc
import numpy as np, h5py, scipy.sparse as sp
from scipy.stats import spearmanr
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(REPO / "experiments/323-current-family-validation"))
import exact_de323 as E
OUT = Path(__file__).resolve().parent / "work"; OUT.mkdir(exist_ok=True)
K = 2200; rng = np.random.default_rng(20260909)

def load(ctx):
    with h5py.File(REPO / f"data/challenge/context_{ctx}.h5ad", "r") as h:
        X = h["X"]; M = sp.csr_matrix((X["data"][...], X["indices"][...], X["indptr"][...]), shape=tuple(X.attrs["shape"]))
        codes = h["obs/ntc_id/codes"][...]; cats = h["obs/ntc_id/categories"][...].astype(str)
    return M, codes, cats

res = {}
for ctx in "ABC":
    t0 = time.time(); M, codes, cats = load(ctx); n_c = len(cats)
    ctrl = E.control_block(M, min_cpm=5.)
    G = int(ctrl.keep.sum()); Z = np.zeros((n_c, G)); SIG = np.zeros((n_c, G), dtype=bool); nsig = []
    for k in range(n_c):
        rows = np.where(codes == k)[0]
        z, p, q, lfc = E.de_from_counts(M[rows], ctrl)
        Z[k] = np.abs(np.nan_to_num(z)); SIG[k] = q < .05; nsig.append(int(SIG[k].sum()))
    key = Z.mean(0); freq = SIG.mean(0); cpm = ctrl.ref_lin
    rho_key = float(spearmanr(key, cpm).correlation); rho_freq = float(spearmanr(freq, cpm).correlation)
    top_cpm = np.argsort(-cpm, kind="stable")[:K]
    rec_key, rec_cpm, rec_freq = [], [], []
    for k in range(n_c):
        own = SIG[k]
        if own.sum() == 0: continue
        others = np.ones(n_c, dtype=bool); others[k] = False
        kk = Z[others].mean(0); ff = SIG[others].mean(0)
        top_k = np.argsort(-kk, kind="stable")[:K]; top_f = np.argsort(-ff, kind="stable")[:K]
        rec_key.append(own[top_k].sum() / own.sum()); rec_cpm.append(own[top_cpm].sum() / own.sum()); rec_freq.append(own[top_f].sum() / own.sum())
    rec_key, rec_cpm, rec_freq = map(np.array, (rec_key, rec_cpm, rec_freq))
    d = rec_key - rec_cpm; df = rec_freq - rec_cpm
    boots = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(10000)]
    bootf = [df[rng.integers(0, len(df), len(df))].mean() for _ in range(10000)]
    res[ctx] = dict(n_constructs=n_c, universe=G, nsig_median=float(np.median(nsig)), nsig=nsig,
                    spearman_key_cpm=rho_key, spearman_freq_cpm=rho_freq,
                    n_constructs_with_sig=int(len(d)), recall_key=float(rec_key.mean()), recall_cpm=float(rec_cpm.mean()), recall_freq=float(rec_freq.mean()),
                    delta_key_minus_cpm=float(d.mean()), ci95=[float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                    delta_freq_minus_cpm=float(df.mean()), ci95_freq=[float(np.percentile(bootf, 2.5)), float(np.percentile(bootf, 97.5))],
                    seconds=round(time.time() - t0, 1))
    np.savez_compressed(OUT / f"nullkey-{ctx}.npz", key=key, freq=freq, cpm=cpm, keep=ctrl.keep, nsig=np.array(nsig))
    print(ctx, json.dumps({k: v for k, v in res[ctx].items() if k != "nsig"}), flush=True)
    del M, ctrl, Z, SIG; gc.collect()
r1 = all(res[c]["nsig_median"] < 10 for c in "ABC"); r2 = all(res[c]["spearman_key_cpm"] > .95 for c in "ABC")
wins = sum(res[c]["delta_key_minus_cpm"] > 0 and res[c]["ci95"][0] > 0 for c in "ABC")
verdict = "STOP_no_null" if r1 else "STOP_cpm_in_disguise" if r2 else ("ADVANCE_to_H1_test" if wins >= 2 else "STOP_recall_not_better")
res["verdict"] = verdict; res["wins"] = wins
(OUT / "results.json").write_text(json.dumps(res, indent=1)); print("VERDICT", verdict, "wins", wins)
