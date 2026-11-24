#!/usr/bin/env python3
"""517 §2 — rank the target-independent sign sources that remain if Arc's programme does NOT carry (458d null / harm),
on H1, through the validated analytic emitter path (price517 / blockval517), plus accuracy rows.

Sources (P-like on H1's gated axis; z-sums mapped to a pseudo-P 0.5 - 0.5 tanh(z/2), 478b's mapping):
  c4        CLEAN4 (incumbent)                          arc   Arc's H1 programme (reference only: same line as h1-test)
  wx        Weissman (= CLEAN4) + Xaira (HCT116, HEK293T) lab-unit z consensus (450j), h1-test targets' rows excluded
  wxs       + Satija (VIPerturb K562, the only other Flex-style screen)
  ident     -I, the recipient's identity axis: log mean CPM of H1's controls minus the mean over A/B/C's (same chemistry)
  c4id1     z(CLEAN4) + z(-I)           c4id2  z(CLEAN4) + 2 z(-I)
  own       491c: z(CLEAN4) + z(proj_k2 of CLEAN4 on H1's control PCs) + z(-I)   (PCs by a streamed randomized SVD)
  down      all-down (470's weak-target tie)
Rows per source: accuracy (fid-like = all CPM-block genes vs the reference lfc sign; gate-block; gate-complement), the
block-prior swap vs CLEAN4 (ungated; and gated: source on 467's strong 57, CLEAN4 on the weak 28) and the complement hook
on the c4 base (at CLEAN4 tau 0.05's coverage, and on every covered gene). Known answers: price517's KA1-3 (in Ctx), the
none/c4 rows (Δ 0 for the c4 block swap), CLEAN4 fid-like 0.5244 and hook c4 tau .05 nmae -0.00147 (seed 0).
usage: uv run python experiments/517-flex-programme/rank517.py        (< 2 GB)
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts"))
import price517 as P  # noqa: E402
import def_sameline as S  # noqa: E402
RC = P.RC
S.BANK_DIRS.setdefault("orion_hek293t", ROOT / "data/donors/orion_hek293t")
S.BANK_DIRS.setdefault("orion_hct116", ROOT / "data/donors/orion_hct116")
S.BANK_DIRS.setdefault("viperturb_k562", ROOT / "data/donors/viperturb_k562")
H1T = ROOT / "data/harness/h1-test.h5ad"


def zc(P_):
    """z of (0.5 - P) over covered genes (positive = up); NaN where uncovered."""
    v = 0.5 - P_
    m = np.isfinite(v)
    out = np.full(v.shape, np.nan)
    out[m] = v[m] / np.std(v[m])
    return out


def zsum(*zs, w=None):
    w = w or [1.0] * len(zs)
    acc = np.zeros(zs[0].shape); cov = np.zeros(zs[0].shape, bool)
    for wi, z in zip(w, zs):
        acc += wi * np.nan_to_num(z); cov |= np.isfinite(z)
    return np.where(cov, 0.5 - 0.5 * np.tanh(acc / 2.0), np.nan)


def zmean(*zs):
    st = np.vstack(zs)
    with np.errstate(invalid="ignore"):
        m = np.nanmean(st, 0)
    return np.where(np.isfinite(m), 0.5 - 0.5 * np.tanh(m / 2.0), np.nan)


def h1_controls_stream(keep, fn, chunk_rows=4000):
    """Call fn(X) on chunks of H1-test's control cells: X = CPM on the kept genes (float32, cells x kept)."""
    with h5py.File(H1T, "r") as f:
        cats = f["obs"]["target_gene"]["categories"][:].astype(str)
        codes = f["obs"]["target_gene"]["codes"][:]
        rows = np.where(codes == int(np.where(cats == "non-targeting")[0][0]))[0]
        indptr = f["X"]["indptr"][:]
        kept_col = np.full(keep.size, -1, np.int64)
        kept_col[np.where(keep)[0]] = np.arange(int(keep.sum()))
        i = 0
        while i < rows.size:
            a = rows[i]; j = i
            while j + 1 < rows.size and rows[j + 1] < a + chunk_rows:
                j += 1
            lo, hi = indptr[a], indptr[rows[j] + 1]
            dat = f["X"]["data"][lo:hi].astype(np.float64); ind = f["X"]["indices"][lo:hi]
            X = np.zeros((j - i + 1, int(keep.sum())), np.float32)
            for k in range(i, j + 1):
                rr = rows[k]; s_, t_ = indptr[rr] - lo, indptr[rr + 1] - lo
                d, ix = dat[s_:t_], ind[s_:t_]; c = kept_col[ix]; m = c >= 0
                X[k - i, c[m]] = (1e6 * d[m] / d.sum()).astype(np.float32)
            fn(X)
            i = j + 1
    return int(rows.size)


def own_axes(keep, genes_keep, Pc4, k=2, n_iter=2, seed=0):
    """-I and proj_k(CLEAN4) on H1's control PCs (log1p CPM, centred), streamed (491c / 478b's construction)."""
    G = int(keep.sum())
    s = {"n": 0, "sum": np.zeros(G), "sumlog": np.zeros(G)}

    def p1(X):
        s["n"] += X.shape[0]; s["sum"] += X.sum(0, dtype=np.float64); s["sumlog"] += np.log1p(X).sum(0, dtype=np.float64)
    n = h1_controls_stream(keep, p1)
    m_cpm = s["sum"] / n; mu = (s["sumlog"] / n).astype(np.float32)
    rg = np.random.default_rng(seed)
    Om = rg.standard_normal((G, k + 20)).astype(np.float32)
    Ys = []
    h1_controls_stream(keep, lambda X: Ys.append((np.log1p(X) - mu) @ Om))
    Y = np.vstack(Ys)
    for _ in range(n_iter):
        Q, _ = np.linalg.qr(Y)
        Z = np.zeros((G, Q.shape[1]), np.float64); off = [0]

        def p2(X, Q=Q, Z=Z, off=off):
            r0 = off[0]; r1 = r0 + X.shape[0]
            Z[:] += ((np.log1p(X) - mu).T @ Q[r0:r1]).astype(np.float64); off[0] = r1
        h1_controls_stream(keep, p2)
        Ys = []
        h1_controls_stream(keep, lambda X, Z=Z: Ys.append((np.log1p(X) - mu) @ Z.astype(np.float32)))
        Y = np.vstack(Ys)
    Q, _ = np.linalg.qr(Y)
    B = np.zeros((Q.shape[1], G), np.float64); off = [0]

    def p3(X):
        r0 = off[0]; r1 = r0 + X.shape[0]
        B[:] += (Q[r0:r1].T @ (np.log1p(X) - mu)).astype(np.float64); off[0] = r1
    h1_controls_stream(keep, p3)
    _, sv, Vt = np.linalg.svd(B, full_matrices=False)
    V = Vt[:k].T
    mC = np.where(np.isfinite(Pc4), 0.5 - Pc4, 0.0)
    proj = V @ (V.T @ mC)
    offi = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    opos = {g: i for i, g in enumerate(offi)}
    abc = []
    for ctx in "ABC":
        mc = np.load(ROOT / f"experiments/491-own-control-programme/work/ctrl_meancpm_{ctx}.npy")
        abc.append(np.array([mc[opos[g]] if g in opos else np.nan for g in genes_keep]))
    with np.errstate(invalid="ignore"):
        abc = np.nanmean(np.vstack([np.log(np.maximum(v, 1e-3)) for v in abc]), 0)
    ident = np.log(np.maximum(m_cpm, 1e-3)) - abc
    return dict(ident=np.where(np.isfinite(ident), ident, np.nan), proj=proj, sv=sv[:5], n_ctrl=n, m_cpm=m_cpm)


def main():
    t0 = time.monotonic()
    zr = np.load(RC.REF, allow_pickle=False)
    genes_full = zr["genes"]; keep = np.asarray(zr["keep"], dtype=bool); T = zr["targets"].astype(str)
    zr.close()
    gk = np.asarray(genes_full)[keep].astype(str)
    with h5py.File(H1T, "r") as f:
        vi = f["var"]["_index"]
        vg = (vi[:] if not isinstance(vi, h5py.Group) else vi["values"][:]).astype(str)
    assert np.array_equal(vg, np.asarray(genes_full).astype(str)), "h1-test var order != reference genes"
    SRC, _ = P.load_sources(genes_full, keep)
    Pc4, Parc = SRC["c4"], SRC["arc"]
    h1 = set(T)

    def bankP(b):
        delta, perts, genes, _ = S.load_bank(b)
        return S.programme(delta, perts, genes, list(map(str, gk)), h1)[1]
    Pxa = np.nanmean(np.vstack([bankP("orion_hct116"), bankP("orion_hek293t")]), 0)
    Psa = bankP("viperturb_k562")
    ax = own_axes(keep, gk, Pc4)
    print(f"[517r] own axes: {ax['n_ctrl']} H1 controls; top sv {np.round(ax['sv'], 1)}; ({time.monotonic() - t0:.0f}s)", flush=True)
    zI = -ax["ident"] / np.nanstd(ax["ident"])
    zP = ax["proj"] / np.std(ax["proj"][ax["proj"] != 0])
    Z4 = zc(Pc4)
    sources = {
        "c4": Pc4, "arc": Parc,
        "wx": zmean(Z4, zc(Pxa)), "wxs": zmean(Z4, zc(Pxa), zc(Psa)),
        "ident": np.where(np.isfinite(zI), 0.5 - 0.5 * np.tanh(zI / 2), np.nan),
        "c4id1": zsum(Z4, zI), "c4id2": zsum(Z4, zI, w=[1, 2]),
        "own": zsum(Z4, zP, zI),
        "down": np.full(Pc4.shape, 1.0),
    }
    c4_block = RC.prior_signs_full(genes_full)[keep]
    strong, _ = P.strong_set(T)
    n_c4_hook = int((np.isfinite(Pc4) & (np.abs(Pc4 - 0.5) >= 0.05)).sum())
    out = {"n_c4_hook_genes": n_c4_hook, "sources": {}}
    for seed in (0, 1, 2):
        ctx = P.Ctx("c4", seed, lambda i, t: c4_block)
        rng = np.random.default_rng(5172 + seed)
        for name, Pv in sources.items():
            rec = out["sources"].setdefault(name, {"coverage": float(np.isfinite(Pv).mean()),
                                                   "frac_down": float((P.sign_from_P(Pv) < 0).mean())})
            sgn = P.sign_from_P(Pv)                    # block rule: uncovered -> down (run400)
            if seed == 0:
                # accuracy rows (seed-independent: the truth and the block membership are fixed)
                hit = {"fidlike": [0, 0], "gate_block": [0, 0], "gate_comp": [0, 0]}
                for B in ctx.blocks:
                    i = B["i"]; lr, qr = ctx.lfc_real[i], ctx.q_real[i]
                    own_gene = ctx.G == B["t"]
                    adj = np.isfinite(lr) & (lr != 0) & ~own_gene
                    gate = adj & (qr < 0.05)
                    cov = np.isfinite(Pv)
                    for key, m in (("fidlike", adj & B["blk"] & cov), ("gate_block", gate & B["blk"] & cov),
                                   ("gate_comp", gate & B["comp"] & cov)):
                        hit[key][0] += int((sgn[m] == np.sign(lr[m])).sum()); hit[key][1] += int(m.sum())
                rec["acc"] = {k: {"acc": v[0] / max(v[1], 1), "n": v[1]} for k, v in hit.items()}
            # block swap, ungated and gated
            for tagb, pr in (("block_all", lambda i, t, s=sgn: s), ("block_g33", lambda i, t, s=sgn: s if strong[i] else c4_block)):
                lfc_new, _ = ctx.block(pr)
                sc, _ = ctx.score(lfc_new, rng, n_boot=500)
                rec.setdefault(tagb, []).append(sc)
            # hooks on the c4 base: coverage matched to CLEAN4 tau .05, and all covered genes
            conf = np.where(np.isfinite(Pv), np.abs(Pv - 0.5), -1.0)
            thr = np.sort(conf)[::-1][min(n_c4_hook, conf.size) - 1] if name != "down" else 0.0
            for tagh, th in (("hook_matched", thr), ("hook_all", 0.0)):
                hs = np.where((conf >= th) & (conf >= 0), P.sign_from_P(Pv), 0.0)
                lfc_new, hrec = ctx.hook(lambda i, t, hs=hs: hs)
                sc, _ = ctx.score(lfc_new, rng, n_boot=500)
                sc["flips"] = float(np.median(hrec["flips"])); sc["acc"] = ctx.gate_accuracy(lambda i, t, hs=hs: hs, "comp")
                rec.setdefault(tagh, []).append(sc)
            print(f"s{seed} {name:6s} " + " | ".join(
                f"{k} ov {rec[k][-1]['overall_equiv']:+.5f} nmae {rec[k][-1]['nmae']:+.5f} fid {rec[k][-1]['fid']:+.5f} "
                f"reach {rec[k][-1]['reach']:+.5f}" for k in ("block_all", "block_g33", "hook_matched", "hook_all"))
                + (f" | acc {json.dumps({k: round(v['acc'], 4) for k, v in rec['acc'].items()})}" if seed == 0 else ""),
                flush=True)
        del ctx
    (HERE / "work/rank517.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"done {time.monotonic() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
