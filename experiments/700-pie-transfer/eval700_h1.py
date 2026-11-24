#!/usr/bin/env python3
"""700 — does PIE beat our own H1 predictions? A direct, label-level comparison on h1-test (85 targets, real truth).

Truth: lib521's H1 reference (rq, rp, rl on the kept gene axis; the same DE table the 537/521 instruments score against) and
the scorer's bulk effect computed here from data/harness/h1-test.h5ad (log1p(5e4 * p) pseudobulk minus control).
Ours: lib521's `key` (435's production cap15 bulk effect for H1, built without any H1 donor) and our production emitter's
DE output on H1 (op, ol: seed 0, 435's run).
PIE: predictions.parquet from infer700.py for dataset h1test / context ARC_H1 (p_de, lfc_pred, delta_p_pred).

Leakage note: PIE's xdataset checkpoint trained on VCC25 H1 train+val (190 perturbations, 1 % sampling weight). None of the 85
h1-test targets is among them. Same-context evidence is keyed by (dataset, context) and this dir's dataset is h1test, so
none is served. What PIE knows about H1 beyond the control mean is in its weights and H1's context text.

Reports (per target, then the mean over targets):
  DE ranking   AUPRC of a gene score against truth significance (q < 0.05, own target excluded), for PIE p_de, |PIE lfc|,
               |our key|, our emitted -log10 p;
  signs        accuracy of sign(score) against sign(truth lfc) on the truth's significant genes, for PIE lfc, PIE delta_p,
               our key, our emitted lfc;
  n_real       Spearman over targets of PIE sum(p_de) and of sum|our key| against n_real;
  bulk         cosine to the truth's bulk effect (panel genes removed, as pds does) for PIE delta_p, PIE lfc, our key, and
               blends key + w * z(PIE delta_p); plus a mini-pds: the rank of the right target among all 85 truths.
usage: uv run python experiments/700-pie-transfer/eval700_h1.py <predictions.parquet> [out.json]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402

TS = 5e4


def auprc(score: np.ndarray, y: np.ndarray) -> float:
    if y.sum() == 0:
        return float("nan")
    o = np.argsort(-score, kind="stable")
    yy = y[o]
    tp = np.cumsum(yy)
    prec = tp / np.arange(1, len(yy) + 1)
    return float((prec * yy).sum() / yy.sum())


def load_pie(path: str):
    t = pq.read_table(path)
    meta = json.loads(t.schema.metadata[b"pie"])
    df = t.to_pandas()
    df = df[df["dataset"] == "h1test"]
    genes = meta["genes"]["h1test"]
    out = {}
    for _, r in df.iterrows():
        out[r["perturbation"]] = {k: np.asarray(r[k], dtype=np.float64) for k in ("p_de", "lfc_pred", "delta_p_pred")}
    return genes, out


def truth_bulk(targets, genes):
    import anndata as ad
    a = ad.read_h5ad(ROOT / "data/harness/h1-test.h5ad")
    pos = {g: i for i, g in enumerate(a.var_names)}
    cols = np.array([pos[g] for g in genes])
    tg = a.obs["target_gene"].to_numpy()
    X = a.X.tocsr() if hasattr(a.X, "tocsr") else a.X

    def pb(mask):
        P = np.asarray(X[mask].sum(0)).ravel().astype(np.float64)
        return np.log1p(TS * P / P.sum())[cols]

    bc = pb(tg == "non-targeting")
    return np.stack([pb(tg == t) - bc for t in targets]), bc


def zrow(m):
    s = m.std(1, keepdims=True)
    return np.where(s > 0, (m - m.mean(1, keepdims=True)) / np.where(s > 0, s, 1), 0.0)


def main():
    pred_path = sys.argv[1]
    S = L.load_state(0)
    targets, genes = list(S["targets"]), list(S["genes"])
    pg, P = load_pie(pred_path)
    ppos = {g: i for i, g in enumerate(pg)}
    have = np.array([g in ppos for g in genes])
    idx = np.array([ppos.get(g, 0) for g in genes])
    T, G = len(targets), len(genes)
    pde = np.zeros((T, G)); plfc = np.zeros((T, G)); pdp = np.zeros((T, G))
    missing_t = [t for t in targets if t not in P]
    for i, t in enumerate(targets):
        if t in P:
            pde[i] = np.where(have, P[t]["p_de"][idx], 0.0)
            plfc[i] = np.where(have, P[t]["lfc_pred"][idx], 0.0)
            pdp[i] = np.where(have, P[t]["delta_p_pred"][idx], 0.0)
    own = S["own"]
    sig = (S["rq"] < L.ALPHA) & np.isfinite(S["rq"]) & ~own
    tested = np.isfinite(S["rq"]) & ~own
    rl = S["rl"]
    key = S["key"]
    op = np.where(np.isfinite(S["op"]), S["op"], 1.0)
    ol = np.where(np.isfinite(S["ol"]), S["ol"], 0.0)
    res = {"n_targets": T, "n_genes": G, "pie_gene_coverage": float(have.mean()), "pie_missing_targets": missing_t}
    # DE ranking
    scores = {"pie_p_de": pde, "pie_abs_lfc": np.abs(plfc), "our_abs_key": np.abs(key), "our_emit_neglog10p": -np.log10(np.clip(op, 1e-300, 1))}
    res["auprc"] = {}
    for k, sc in scores.items():
        v = [auprc(sc[i][tested[i]], sig[i][tested[i]]) for i in range(T)]
        res["auprc"][k] = float(np.nanmean(v))
    res["auprc_random"] = float(np.nanmean([sig[i][tested[i]].mean() if sig[i].any() else np.nan for i in range(T)]))
    # signs on truth DE genes
    res["sign_acc_on_truth_sig"] = {}
    for k, sc in {"pie_lfc": plfc, "pie_delta_p": pdp, "our_key": key, "our_emit_lfc": ol}.items():
        num = den = 0
        per = []
        for i in range(T):
            m = sig[i] & (sc[i] != 0) & np.isfinite(rl[i]) & (rl[i] != 0)
            if m.sum():
                a = (np.sign(sc[i][m]) == np.sign(rl[i][m]))
                num += a.sum(); den += m.sum(); per.append(a.mean())
        res["sign_acc_on_truth_sig"][k] = dict(pooled=float(num / max(den, 1)), per_target_mean=float(np.mean(per)),
                                               coverage=float(den / max(sig.sum(), 1)))
    # fid-like: sign agreement with the truth's lfc on OUR call set (head | CPM block, adjudicable genes, DE or not), the
    # quantity fid's k / n_pred reads when n_pred >= n_real; per source, and for PIE restricted to p_de >= 0.5
    c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=False)
    c4pos = {g: i for i, g in enumerate(map(str, c4["genes"]))}
    pdn = np.array([c4["P_down"][c4pos[g]] if g in c4pos else np.nan for g in genes], dtype=np.float64)
    c4sign = np.where(np.isfinite(pdn), np.where(pdn > 0.5, -1.0, 1.0), 0.0)
    adj = np.isfinite(rl) & (rl != 0) & ~own
    C = S["C"] & adj
    res["callset_sign_agreement"] = {}
    for k, sc in {"clean4": np.broadcast_to(c4sign, rl.shape), "pie_lfc": plfc, "our_key": key,
                  "pie_lfc_pde05_else_clean4": np.where(pde >= 0.5, np.sign(plfc), c4sign[None, :])}.items():
        m = C & (sc != 0)
        res["callset_sign_agreement"][k] = dict(acc=float((np.sign(sc[m]) == np.sign(rl[m])).mean()),
                                                coverage=float(m.sum() / max(C.sum(), 1)))
    res["callset_pie_pde05_frac"] = float(((pde >= 0.5) & C).sum() / max(C.sum(), 1))
    # n_real
    nreal = sig.sum(1)
    res["n_real_spearman"] = {"pie_sum_p_de": float(spearmanr(pde.sum(1), nreal).correlation),
                              "our_sum_abs_key": float(spearmanr(np.abs(key).sum(1), nreal).correlation),
                              "pie_sum_abs_delta_p": float(spearmanr(np.abs(pdp).sum(1), nreal).correlation)}
    # bulk direction
    tb, bc = truth_bulk(targets, genes)
    panel = np.isin(np.array(genes), np.array(targets))
    keepg = ~panel

    def cos_rows(a, b):
        a, b = a[:, keepg], b[:, keepg]
        na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
        return (a * b).sum(1) / np.where(na * nb > 0, na * nb, np.inf)

    def minipds(pred):
        a, b = pred[:, keepg], tb[:, keepg]
        an = a / np.maximum(np.linalg.norm(a, axis=1, keepdims=True), 1e-12)
        bn = b / np.maximum(np.linalg.norm(b, axis=1, keepdims=True), 1e-12)
        C = an @ bn.T
        r = [(C[i] > C[i, i]).sum() + 0.5 * ((C[i] == C[i, i]).sum() - 1) for i in range(T)]
        return float(1 - np.mean(r) / (T - 1))

    w_tb = (tb[:, keepg] ** 2).sum(1)
    res["bulk"] = {}
    cands = {"pie_delta_p": pdp, "pie_lfc": plfc, "our_key": key}
    # PIE lfc -> the scorer's effect space with H1's control pseudobulk (as make_pie_dump700 does for A/B/C)
    pc = np.expm1(bc); pc = pc / pc.sum()
    epie = np.zeros_like(plfc)
    tested_g = 1e6 * pc >= 5.0      # amendment 2: PIE's lfc only on genes >= 5 CPM in the controls
    for i in range(T):
        q = pc * np.exp2(np.where(tested_g, plfc[i], 0.0)); q = q / q.sum()
        epie[i] = np.log1p(5e4 * q) - np.log1p(5e4 * pc)
    cands["pie_effect"] = epie
    seen_set = set(json.load(open(HERE / "work/pie_train_perts.json")))
    seen_t = np.array([t in seen_set for t in targets])

    def ublend(base, pie, w, only):
        out = base.copy()
        for i in range(T):
            if only is not None and not only[i]:
                continue
            b, v = base[i, keepg], pie[i, keepg]
            nb, nv = np.linalg.norm(b), np.linalg.norm(v)
            if nb == 0 or nv == 0:
                continue
            m = (1 - w) * b / nb + w * v / nv
            out[i, keepg] = nb * m / np.linalg.norm(m)
        return out
    for w in (0.25, 0.5, 0.75):
        cands[f"blend_all_w{w}"] = ublend(key, epie, w, None)
        cands[f"blend_seen_w{w}"] = ublend(key, epie, w, seen_t)
    for k, v in cands.items():
        c = cos_rows(v, tb)
        res["bulk"][k] = dict(mean_cos=float(np.mean(c)), energy_weighted_cos=float((c * w_tb).sum() / w_tb.sum()),
                              minipds=minipds(v))
    # seen / unseen split: targets with same-pert evidence in PIE's training split (train.json perts of any dataset)
    seen = set(json.load(open(HERE / "work/pie_train_perts.json")))
    sm = np.array([t in seen for t in targets])
    res["n_seen"] = int(sm.sum())
    res["by_seen"] = {}
    for lab, msk in (("seen", sm), ("unseen", ~sm)):
        if msk.sum() < 3:
            continue
        r2 = {}
        for k, v in cands.items():
            c = cos_rows(v[msk], tb[msk])
            ww = w_tb[msk]
            r2[f"bulk_{k}"] = dict(mean_cos=float(np.mean(c)), energy_weighted_cos=float((c * ww).sum() / ww.sum()))
        for k, sc in {"pie_lfc": plfc, "our_key": key}.items():
            m = sig[msk] & (sc[msk] != 0) & np.isfinite(rl[msk]) & (rl[msk] != 0)
            r2[f"sign_{k}"] = float((np.sign(sc[msk][m]) == np.sign(rl[msk][m])).mean()) if m.sum() else None
        for k, sc in scores.items():
            r2[f"auprc_{k}"] = float(np.nanmean([auprc(sc[i][tested[i]], sig[i][tested[i]]) for i in np.flatnonzero(msk)]))
        res["by_seen"][lab] = dict(n=int(msk.sum()), **r2)
    out = sys.argv[2] if len(sys.argv) > 2 else str(HERE / "work/eval700_h1.json")
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
