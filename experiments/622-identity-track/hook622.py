#!/usr/bin/env python3
"""622o: how often is a new line "K562-like" for hooks? The weight of idall3's stress tail (asked by [21dd4c] 2026-10-02
~14:00Z for idall3's board prereg).

613d's add rule (block613d.pred_dnmae, copied verbatim below and checked against the import) on the 12 working Broad-16
lines. The 4 reserved lines are never parsed: matrices via broad622.load_matrix (usecols), covariate rows skipped before
any field is read.
Per line L, on UNCALLED gate genes:
  gate     = authors' fdr < 0.05, finite lfc != 0, gene != KO; truth T = lfc on gate genes.
  call set = head 400 by the other 11 lines' mean |lfc| for the same KO + block 2,200 by L's OR-control CPM
             (bank613f's construction); uncalled = tested & ~head & ~block.
  key sign = sign of the other lines' mean lfc for the same KO (0 where no other line has it).
  add      = fold480.mde_add(cnt), cnt = control CPM × library / 1e6. The library is L's own mean UMIs per cell (Broad's
             avg_umis_per_ko), with a beside-run at 20,000.
  hooks    c4h    = CLEAN4 sign where |P - 1/2| >= 0.05 (496a's form; no bank shares a Broad line);
           idall3 = sign of z(P_C4 - 1/2) + 0.5 z(I) on CLEAN4-covered genes, ungated (all uncalled genes). I = same-lab
                    identity (622g): log1p CPM_L - mean over the other 11 lines, genes > 5 CPM in L. zc is centred at 1/2
                    (604 v3 / 613d).
  pools    = 521's regimes f 1 / 0.5 / 0.35 via ds612b.pool_T. Broad gives no raw p, so the authors' fdr orders the pool
             (BH is monotone in p).
Reads: (a) c4h - no hook and (b) idall3 - c4h, predicted Δnmae (raw; + = worse). Counted: lines with (a) > 0 (the hook
hurts) and (b) > 0 (idall3 worse than plain).
Caveat stated in advance: KO gates lean UP (622g), so read the ORDER of hook forms, not the levels.
KEY variants (added after the first run; the spec's key is the first):
  `samelab` = the spec: head and key sign from the other 11 Broad lines' SAME-KO response. That is a near-oracle same-lab
              per-target programme, far stronger than any D/E/F key.
  `donor`   = production-like: head 400 and key sign from K562 GWPS's CRISPRi response to the same target
              (data/donors/k562_gwps delta, another lab and modality, as our donor key is for D/E/F). Targets absent from
              GWPS get no head and key 0.
usage: uv run python experiments/622-identity-track/hook622.py [samelab|donor]   -> work/hook622_<mode>.json
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import broad622 as B  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument"))
import fold480 as F  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments/612-key-gated-hook"))
from ds612b import pool_T  # noqa: E402

REG = (1.0, 0.5, 0.35)


def pred_dnmae(hook_sign, key_sign, add, Tpool, uncalled):
    """block613d.pred_dnmae, verbatim: the mean over scored targets of Σ Δ|err| / Σ|T| on the pool (426's add rule)."""
    G = np.isfinite(Tpool)
    den = np.where(G, np.abs(Tpool), 0.0).sum(1)
    ok = (G.sum(1) >= 10) & (den > 0)
    tdn = np.sign(np.nan_to_num(Tpool))
    hs = np.broadcast_to(hook_sign, Tpool.shape)
    act = G & uncalled & (hs != 0)
    flip = act & (key_sign != 0) & (hs != key_sign)
    new = act & (key_sign == 0)
    corr = hs == tdn
    d = np.where(flip, np.where(corr, -2.0, 2.0) * add, 0.0) + np.where(new, np.where(corr, -1.0, 1.0) * add, 0.0)
    rows = d.sum(1) / np.where(den > 0, den, 1.0)
    return float(np.mean(rows[ok])), int(ok.sum()), int(flip.sum()), int(new.sum())


def pred_dnmae_bal(hook_sign, key_sign, add, Tpool, uncalled):
    """pred_dnmae with the base rate BALANCED per target. Each gate gene's Δ|err| is weighted so that UP and DOWN truth
    carry equal total weight on the target's acting cells. This removes the KO gates' UP lean (622g), which otherwise
    favours whichever form signs fewer genes DOWN. Not 613d's rule; a diagnostic beside it."""
    G = np.isfinite(Tpool)
    den = np.where(G, np.abs(Tpool), 0.0).sum(1)
    ok = (G.sum(1) >= 10) & (den > 0)
    tdn = np.sign(np.nan_to_num(Tpool))
    hs = np.broadcast_to(hook_sign, Tpool.shape)
    act = G & uncalled & (hs != 0)
    flip = act & (key_sign != 0) & (hs != key_sign)
    new = act & (key_sign == 0)
    corr = hs == tdn
    nu = (G & (tdn > 0)).sum(1, keepdims=True); nd = (G & (tdn < 0)).sum(1, keepdims=True)
    w = np.where(tdn > 0, 0.5 * (nu + nd) / np.maximum(nu, 1), 0.5 * (nu + nd) / np.maximum(nd, 1))
    d = w * (np.where(flip, np.where(corr, -2.0, 2.0) * add, 0.0) + np.where(new, np.where(corr, -1.0, 1.0) * add, 0.0))
    rows = d.sum(1) / np.where(den > 0, den, 1.0)
    return float(np.mean(rows[ok])), int(ok.sum())


def check_verbatim():
    """The copy must equal 613d's function (imported with its own path setup) on random inputs."""
    try:
        sys.path.insert(0, str(ROOT / "experiments/613-block-prior-validation"))
        sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
        sys.path.insert(0, str(ROOT / "experiments/639-identity-axis-lolo"))
        import block613d as BD  # noqa: PLC0415
    except Exception as e:                       # noqa: BLE001
        return f"not imported ({type(e).__name__}); copy unchecked"
    rng = np.random.default_rng(6220)
    for _ in range(20):
        T = np.where(rng.random((30, 200)) < 0.3, rng.normal(0, 1, (30, 200)), np.nan)
        hs = rng.choice([-1.0, 0.0, 1.0], 200); ks = rng.choice([-1.0, 0.0, 1.0], (30, 200))
        add = rng.uniform(0, 0.1, (1, 200)); unc = rng.random((30, 200)) < 0.6
        assert pred_dnmae(hs, ks, add, T, unc) == BD.pred_dnmae(hs, ks, add, T, unc)
    return "equal to block613d.pred_dnmae on 20 random cases"


def line_umis(lines):
    out = {}
    with open(B.P / "cell_quality_covariates.csv") as fh:
        r = csv.reader(fh)
        hdr = next(r)
        for row in r:
            if row[0] in B.RESERVED:             # skipped before any field is parsed
                continue
            if row[0] in lines:
                out[row[0]] = float(row[hdr.index("avg_umis_per_ko")])
    return out


KEYMODE = sys.argv[1] if len(sys.argv) > 1 else "samelab"
assert KEYMODE in ("samelab", "donor"), KEYMODE


def donor_key(gl):
    D = np.load(ROOT / "data/donors/k562_gwps/delta.npy", mmap_mode="r")
    perts = [ln.strip() for ln in open(ROOT / "data/donors/k562_gwps/perts.txt") if ln.strip()]
    genes = [ln.strip() for ln in open(ROOT / "data/donors/k562_gwps/genes.txt") if ln.strip()]
    ppos = {p: i for i, p in enumerate(perts)}
    gpos = {g: i for i, g in enumerate(genes)}
    col = np.array([gpos.get(g, -1) for g in gl])
    def row(t):
        i = ppos.get(t)
        if i is None:
            return None
        r = np.asarray(D[i], float)
        out = np.full(len(gl), np.nan)
        out[col >= 0] = r[col[col >= 0]]
        return out
    return row


def main() -> None:
    print("pred_dnmae copy:", check_verbatim(), f"| key mode {KEYMODE}", flush=True)
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
    assert list(gl) == list(gf) and cl == cf
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl})
    assert len(lines) == 12 and not (set(lines) & B.RESERVED)
    umis = line_umis(set(lines))
    all_t = set(k for _, k in cl)
    ctrl_cols = {ln: [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t] for ln in lines}
    gpos = {g: i for i, g in enumerate(gp)}
    rowmap = np.array([gpos.get(g, -1) for g in gl])
    cpm = {}
    for ln in lines:
        s = np.nansum(PB[:, ctrl_cols[ln]], axis=1)
        absent = np.all(~np.isfinite(PB[:, ctrl_cols[ln]]), axis=1)
        v = 1e6 * s / np.nansum(s); v[absent] = np.nan
        cpm[ln] = np.where(rowmap >= 0, v[np.maximum(rowmap, 0)], np.nan)
    c4 = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    nm = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[nm.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    cmap = dict(zip(c4["genes"].astype(str), Pc))
    C4 = np.array([cmap.get(g, np.nan) for g in gl])
    c4h = np.where(np.isfinite(C4) & (np.abs(C4 - 0.5) >= 0.05), np.where(C4 > 0.5, -1.0, 1.0), 0.0)
    zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
    colline = np.array([ln for ln, _ in cl]); coltgt = np.array([k for _, k in cl])
    drow = donor_key(gl) if KEYMODE == "donor" else None
    out = {}
    for ln in lines:
        others = [o for o in lines if o != ln]
        I = np.where(cpm[ln] > 5, np.log1p(cpm[ln]) - np.nanmean([np.log1p(cpm[o]) for o in others], 0), np.nan)
        zi = (I - np.nanmean(I)) / np.nanstd(I)
        sc = np.where(np.isfinite(zc), zc + 0.5 * np.nan_to_num(zi), np.nan)
        idall3 = np.where(np.isfinite(sc) & (sc != 0), np.where(sc > 0, -1.0, 1.0), 0.0)
        c4all = np.where(np.isfinite(C4) & (C4 != 0.5), np.where(C4 > 0.5, -1.0, 1.0), 0.0)     # idall3's coverage, no identity
        id_on_c4h = np.where(c4h != 0, idall3, 0.0)                                            # identity inside c4h's coverage
        # c4h + an identity-TAIL extension: where c4h leaves a covered gene unsigned (|P - 1/2| < 0.05), sign it by the
        # combined score only if the line's identity is extreme (|z(I)| >= 1.5); the low-I tail goes UP (622n)
        tailx = np.where((c4h == 0) & np.isfinite(C4) & (np.abs(np.nan_to_num(zi)) >= 1.5), idall3, c4h)
        jl = np.flatnonzero(colline == ln)
        tl = coltgt[jl]
        L, Fq = LFC[:, jl].T, FDR[:, jl].T                     # targets x genes
        own = tl[:, None] == gl[None, :]
        tested = np.isfinite(L) & ~own
        gate = tested & (L != 0) & (Fq < 0.05)
        T = np.where(gate, L, np.nan)
        head = np.zeros(L.shape, bool)
        key = np.zeros(L.shape)
        for a, t in enumerate(tl):
            if KEYMODE == "donor":
                r = drow(t)
                if r is not None:
                    head[a, np.argsort(-np.nan_to_num(np.abs(r), nan=-1.0), kind="stable")[:400]] = True
                    key[a] = np.sign(np.nan_to_num(r))
                continue
            jo = [j for j, (l2, k) in enumerate(cl) if k == t and l2 != ln]
            if jo:
                k = np.nanmean(np.abs(LFC[:, jo]), 1)
                head[a, np.argsort(-np.nan_to_num(k, nan=-1.0), kind="stable")[:400]] = True
                key[a] = np.sign(np.nan_to_num(np.nanmean(LFC[:, jo], 1)))
        rank_c = np.nan_to_num(cpm[ln], nan=-1.0)
        block = np.zeros(L.shape, bool)
        for a in range(L.shape[0]):
            block[a, np.argsort(-np.where(head[a], -np.inf, rank_c), kind="stable")[:2200]] = True
        block &= ~head
        unc = tested & ~head & ~block
        res = {"umis": umis[ln]}
        for libname, lib in (("own", umis[ln]), ("20k", 20000.0)):
            add = F.mde_add(cpm[ln] * lib / 1e6)
            add = np.nan_to_num(add)[None, :]
            for f in REG:
                Tp = pool_T(T, Fq, f)
                a = pred_dnmae(c4h, key, add, Tp, unc)
                b = pred_dnmae(idall3, key, add, Tp, unc)
                c = pred_dnmae(c4all, key, add, Tp, unc)
                e = pred_dnmae(id_on_c4h, key, add, Tp, unc)
                ab, bb = pred_dnmae_bal(c4h, key, add, Tp, unc), pred_dnmae_bal(idall3, key, add, Tp, unc)
                cb, eb = pred_dnmae_bal(c4all, key, add, Tp, unc), pred_dnmae_bal(id_on_c4h, key, add, Tp, unc)
                tx, txb = pred_dnmae(tailx, key, add, Tp, unc), pred_dnmae_bal(tailx, key, add, Tp, unc)
                res[f"{libname}/{f}"] = dict(c4h_minus_nohook=a[0], idall3_minus_c4h=b[0] - a[0], idall3_minus_nohook=b[0],
                                             coverage_c4all_minus_c4h=c[0] - a[0], identity_idall3_minus_c4all=b[0] - c[0],
                                             identity_within_c4h=e[0] - a[0], tailx_minus_c4h=tx[0] - a[0], bal_tailx_minus_c4h=txb[0] - ab[0],
                                             bal_c4h_minus_nohook=ab[0], bal_idall3_minus_c4h=bb[0] - ab[0],
                                             bal_coverage=cb[0] - ab[0], bal_identity=bb[0] - cb[0], bal_identity_within_c4h=eb[0] - ab[0],
                                             n_scored=a[1], c4h_flips=a[2], c4h_new=a[3], idall3_flips=b[2], idall3_new=b[3])
        out[ln] = res
        r1, r35 = res["own/1.0"], res["own/0.35"]
        print(f"{ln:9s} UMIs {umis[ln]:7.0f} | own depth f1: c4h-none {r1['c4h_minus_nohook']:+.5f}  idall3-c4h {r1['idall3_minus_c4h']:+.5f} "
              f"| f0.35: {r35['c4h_minus_nohook']:+.5f} / {r35['idall3_minus_c4h']:+.5f} | 20k f1: {res['20k/1.0']['c4h_minus_nohook']:+.5f} / "
              f"{res['20k/1.0']['idall3_minus_c4h']:+.5f} (n {r1['n_scored']})", flush=True)
    print("\ndecomposition, own depth f 1 (raw predicted nmae; + = worse): coverage = c4all - c4h, identity = idall3 - c4all, "
          "identity within c4h's coverage; then the same BASE-RATE-BALANCED")
    for ln in lines:
        r = out[ln]["own/1.0"]
        print(f"  {ln:9s} coverage {r['coverage_c4all_minus_c4h']:+.5f} identity {r['identity_idall3_minus_c4all']:+.5f} "
              f"within-c4h {r['identity_within_c4h']:+.5f} | balanced: c4h-none {r['bal_c4h_minus_nohook']:+.5f} idall3-c4h "
              f"{r['bal_idall3_minus_c4h']:+.5f} coverage {r['bal_coverage']:+.5f} identity {r['bal_identity']:+.5f} within-c4h {r['bal_identity_within_c4h']:+.5f}")
    summ = {}
    for libname in ("own", "20k"):
        for f in REG:
            k = f"{libname}/{f}"
            a = [out[ln][k]["c4h_minus_nohook"] for ln in lines]; b = [out[ln][k]["idall3_minus_c4h"] for ln in lines]
            extra = {m: dict(n_worse=int(sum(out[ln][k][m] > 0 for ln in lines)), mean=float(np.mean([out[ln][k][m] for ln in lines])))
                     for m in ("coverage_c4all_minus_c4h", "identity_idall3_minus_c4all", "identity_within_c4h", "bal_c4h_minus_nohook",
                               "bal_idall3_minus_c4h", "bal_coverage", "bal_identity", "bal_identity_within_c4h",
                               "tailx_minus_c4h", "bal_tailx_minus_c4h")}
            summ[k] = dict(extra, hook_hurts=int(sum(v > 0 for v in a)), idall3_worse=int(sum(v > 0 for v in b)),
                           mean_c4h_minus_nohook=float(np.mean(a)), mean_idall3_minus_c4h=float(np.mean(b)),
                           max_idall3_minus_c4h=float(np.max(b)), argmax=lines[int(np.argmax(b))])
            print(f"{k:9s}: (a) c4h hurts on {summ[k]['hook_hurts']}/12 (mean {summ[k]['mean_c4h_minus_nohook']:+.5f}); "
                  f"(b) idall3 worse than c4h on {summ[k]['idall3_worse']}/12 (mean {summ[k]['mean_idall3_minus_c4h']:+.5f}, "
                  f"worst {summ[k]['max_idall3_minus_c4h']:+.5f} {summ[k]['argmax']})")
            print("            " + "; ".join(f"{m} worse on {summ[k][m]['n_worse']}/12 (mean {summ[k][m]['mean']:+.5f})"
                                           for m in ("coverage_c4all_minus_c4h", "identity_idall3_minus_c4all", "identity_within_c4h",
                                                     "bal_c4h_minus_nohook", "bal_idall3_minus_c4h", "bal_coverage", "bal_identity",
                                                     "bal_identity_within_c4h", "tailx_minus_c4h", "bal_tailx_minus_c4h")))
    out["_summary"] = summ
    out["_key_mode"] = KEYMODE
    (HERE / f"work/hook622_{KEYMODE}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
