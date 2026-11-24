"""185 Part B/D — realise emission arms on h1-test and write their DE tables.

Every arm delivers a15's per-gene panel totals BIT-IDENTICALLY (sha256 of the delivered
column totals per target must equal a15's, which is re-derived here from a fresh emission
through production's `emit_counts`, exactly `realize171.py`'s code path and seeds).

ARMS
  a15, a15_b        known-answer gates: median n_pred 6035, reach 0.18283 (157/164/171/177)
  armC_a8_flat      known-answer gate: median n_pred 3581, reach 0.20291
  rc_hist_a         EME, |z| multiset = a15's own per-perturbation |z| histogram, RE-ASSIGNED
                    to genes by descending |H| (sign free -> "abs")
  rc_hist_a_b       twin (a15_b's histogram, other placement / sub-gate seed)
  rc_hist_s         same targets, sign forced to sign(H) and clipped to the reachable side
  rc_perm           CONTROL: a15's histogram assigned by a RANDOM permutation of the genes.
                    What would this read if delivering an order did nothing?  The same as
                    rc_hist_a.  If the mechanism is real it reads the constant-key floor.
  rc_k1000/2000/4000 step designs for Part C: top-K by |H| take a15's top-K |z|, the rest 0
  a21_like, a22_like Part D retrodiction rows: the resample law at a21's dose (keep the top
                    28.4% of the gated axis parametric, resample the rest) and a22's dose
                    (resample the top 20% of the gated axis, keep the rest on a15's law),
                    through production's `emit_counts(mode="resample_mix")` as 177 ran it.

Donors: 095's h1-test artifact (k562_gwps + orion_hek293t).  h1_2025 never enters.
No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

import h5py
import numpy as np
from scipy.stats import false_discovery_control, norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import harness185 as H  # noqa: E402
from eme185 import EME  # noqa: E402
import vcc2026.donor_transfer as _dt  # noqa: E402


def _fast_repair_totals(block, totals, rng):
    """177's vectorised equivalent of `donor_transfer._repair_totals` (same contract)."""
    out = np.asarray(block, dtype=np.int64).copy()
    n = out.shape[0]
    resid = np.asarray(totals, dtype=np.int64) - out.sum(0)
    for j in np.flatnonzero(resid > 0):
        np.add.at(out[:, j], rng.integers(0, n, int(resid[j])), 1)
    for j in np.flatnonzero(resid < 0):
        d = int(-resid[j])
        col = out[:, j]
        while d > 0:
            nz = np.flatnonzero(col > 0)
            if nz.size == 0:
                break
            take = rng.choice(nz, size=min(d, nz.size), replace=False)
            col[take] -= 1
            d -= take.size
        out[:, j] = col
    return out


_dt._repair_totals = _fast_repair_totals
from vcc2026.donor_transfer import emit_counts  # noqa: E402

n, EMIT_L, PHI = 400, H.EMIT_L, 4.0
OUT = os.path.join(H.WORK, os.environ.get("OUTTAG", "de185"))
# a21 kept round(0.14*18533)=2595 genes parametric out of a ~9,100-gene challenge gate (28.4%);
# a22 resampled the top 10% = 1,853 genes = 20.4% of the gate.  Fraction-matched to h1-test's
# 10,779-gene DE universe (181 §1 addendum's convention).
K_A21 = int(round(0.284 * 10779))
K_A22 = int(round(0.204 * 10779))

ARMS = [
    ("a15",          "quasi", dict(tail=1.0, ktop=8000, ks=0.25), 90_000),
    ("a15_b",        "quasi", dict(tail=1.0, ktop=8000, ks=0.25), 7_000_000),
    ("armC_a8_flat", "quasi", dict(tail=None, ktop=0, ks=0.25), 90_000),
    ("rc_hist_a",    "eme",   dict(mode="abs", template="a15", K=None, perm=False), 90_000),
    ("rc_hist_a_b",  "eme",   dict(mode="abs", template="a15_b", K=None, perm=False), 7_000_000),
    ("rc_hist_s",    "eme",   dict(mode="signed", template="a15", K=None, perm=False), 90_000),
    ("rc_perm",      "eme",   dict(mode="abs", template="a15", K=None, perm=True), 90_000),
    ("rc_k1000",     "eme",   dict(mode="abs", template="a15", K=1000, perm=False), 90_000),
    ("rc_k2000",     "eme",   dict(mode="abs", template="a15", K=2000, perm=False), 90_000),
    ("rc_k4000",     "eme",   dict(mode="abs", template="a15", K=4000, perm=False), 90_000),
    # hybrid: a15's CALL SET (a p-threshold set, so BH on the same |z| multiset reproduces it
    # exactly) with the |H| order INSIDE the set and inside the complement -- INBOX-4's
    # "sig-first, then |H|" counterfactual (0.21914), realised
    ("rc_hybrid",    "eme",   dict(mode="abs", template="a15", K=None, perm=False, set="a15"), 90_000),
    # the key is a free choice: a15's own |pred lfc| (171's `lfconly` counterfactual 0.22173)
    ("rc_lfc",       "eme",   dict(mode="abs", template="a15", K=None, perm=False, key="lfc"), 90_000),
    ("a21_like",     "rs",    dict(head_param=True, K=K_A21), 90_000),
    ("a22_like",     "rs",    dict(head_param=False, K=K_A22), 90_000),
]
if os.environ.get("ONLY"):
    _only = set(os.environ["ONLY"].split(","))
    ARMS = [a for a in ARMS if a[0] in _only]


class RowReader:
    """Dense rows out of a CSR h5ad without loading it."""

    def __init__(self, path):
        self.f = h5py.File(path, "r")
        self.indptr = self.f["X"]["indptr"][:]
        self.data = self.f["X"]["data"]
        self.indices = self.f["X"]["indices"]

    def rows(self, rows, G):
        rows = np.sort(np.asarray(rows))
        out = np.zeros((rows.size, G), np.float64)
        for k, r in enumerate(rows):
            a, b = int(self.indptr[r]), int(self.indptr[r + 1])
            out[k, self.indices[a:b]] = self.data[a:b]
        return out


def to_block(vals, cols, ptr, G):
    B = np.zeros((n, G), np.int64)
    for r in range(n):
        sl = slice(ptr[r], ptr[r + 1])
        B[r, cols[sl]] = vals[sl]
    return B


def main() -> None:
    t0 = time.time()
    C = H.load_controls()
    cs, tie, keep, cpm, ref_lin = C["cs"], C["tie"], C["keep"], C["cpm"], C["ref_lin"]
    genes_full = C["genes"]
    G_full, G, m = genes_full.size, cs.shape[0], cs.shape[1]
    b_cell = np.log1p(cpm)
    gated = (cpm > H.MIN_CPM) & (b_cell > 0)
    assert gated.sum() == G and np.array_equal(np.flatnonzero(gated), np.flatnonzero(keep))
    gidx = np.flatnonzero(keep)
    TOT = H.load_totals()
    assert np.array_equal(TOT["genes"], genes_full)
    tmask = TOT["targets"] != "non-targeting"
    targets = list(TOT["targets"][tmask])
    totals_all = TOT["totals"][tmask]
    if os.environ.get("MAXT"):                      # smoke-test path only
        targets = targets[:int(os.environ["MAXT"])]
        totals_all = totals_all[:len(targets)]
    Hd = H.load_H()
    assert np.array_equal(Hd["genes"], genes_full)
    hpos = {t: i for i, t in enumerate(Hd["targets"])}
    T_n = len(targets)

    import anndata as ad
    R = ad.read_h5ad(H.REAL, backed="r")
    rlab = np.asarray(R.obs[H.COL].astype(str).values, dtype="U64")
    assert np.array_equal(np.asarray([str(g) for g in R.var_names]), genes_full)
    R.file.close()
    rd = RowReader(H.REAL)
    ci = np.flatnonzero(rlab == H.CTRL)
    # pooled control profile for the resample law (177 §7 contract b)
    pooled_p = os.path.join(H.WORK, "ctrl_pooled.npy")
    if os.path.exists(pooled_p):
        ctot = np.load(pooled_p)
    else:
        ctot = np.zeros(G_full)
        for i in range(0, ci.size, 2000):
            ctot += rd.rows(ci[i:i + 2000], G_full).sum(0)
        np.save(pooled_p, ctot)
    ctrl_rel = ctot / ctot.sum()
    print(f"harness: {ci.size} controls, {G}/{G_full} gated; {T_n} perturbations; arms "
          f"{[a[0] for a in ARMS]} [{time.time()-t0:.0f}s]", flush=True)

    eme = EME(cs, tie, 20000, n)
    real = {k: np.full((T_n, G), np.nan) for k in ("p", "q", "l")}
    out = {a[0]: {k: np.full((T_n, G), np.nan) for k in ("p", "q", "l", "z")} for a in ARMS}
    zdes = {a[0]: np.full((T_n, G), np.nan) for a in ARMS if a[1] == "eme"}
    ztgt = {a[0]: np.full((T_n, G), np.nan) for a in ARMS if a[1] == "eme"}
    svec = {a[0]: np.zeros((T_n, G), np.float32) for a in ARMS}
    npred = {a[0]: [] for a in ARMS}
    digest = {a[0]: [] for a in ARMS}
    diag = {a[0]: [] for a in ARMS}
    absH_all = np.zeros((T_n, G), np.float32)

    for i, t in enumerate(targets):
        # ---- reference DE table (the real perturbation cells)
        Br = rd.rows(np.flatnonzero(rlab == t), G_full)
        s = Br.sum(1)
        s[s == 0] = 1.0
        Br = np.log1p(Br[:, keep] * (1e6 / s[:, None]))
        z = H.wilcoxon_z(cs, tie, np.ascontiguousarray(Br.T))
        real["p"][i] = 2 * norm.sf(np.abs(z))
        real["q"][i] = false_discovery_control(real["p"][i], method="bh")
        real["l"][i] = np.log2((np.expm1(Br).mean(0) + 1e-9) / (ref_lin + 1e-9))
        del Br

        Mfull = totals_all[i].astype(np.int64)
        S = int(Mfull.sum())
        T = S // n
        Tc = np.full(n, T, np.int64)
        Tc[: S % n] += 1
        profile = (Mfull + 0.5) / (n * EMIT_L)          # realize171: floor(n L profile) == Mfull
        pcpm = 1e6 * Mfull / max(S, 1)
        eff = np.abs(np.log2((pcpm + 1e-9) / (cpm + 1e-9)))
        score = np.where(gated, eff, -np.inf)
        order = np.lexsort((np.arange(G_full), -score))
        absH = np.abs(Hd["H"][hpos[t]])
        absH_all[i] = absH[keep]

        blocks = {}
        for name, law, cfg, seed in ARMS:
            rng = np.random.default_rng(seed + i)
            if law == "quasi":
                if cfg["tail"] is None:
                    sv = np.full(G_full, float(cfg["ks"]))
                else:
                    sv = np.full(G_full, float(cfg["tail"]))
                    k = min(int(cfg["ktop"]), int(gated.sum()))
                    sv[order[:k]] = float(cfg["ks"])
                vals, cols, ptr = emit_counts(profile, np.full(n, EMIT_L), 0.1, rng,
                                              mode="quasi", spread=sv, overdispersion=PHI)
                B = to_block(vals, cols, ptr, G_full)
                svec[name][i] = sv[keep]
                dg = {}
            elif law == "eme":
                tmpl = out[cfg["template"]]["z"][i]
                assert np.isfinite(tmpl).all(), "template arm must run before its rc arm"
                zabs_sorted = np.sort(np.abs(tmpl))[::-1]
                if cfg.get("key") == "lfc":
                    key = np.abs(out[cfg["template"]]["l"][i]).astype(np.float64)
                else:
                    key = absH[keep].astype(np.float64)
                jit = np.random.default_rng(1000 + i).random(G) * 1e-12
                if cfg["perm"]:
                    rank_order = np.random.default_rng(5000 + i).permutation(G)
                elif cfg.get("set") == "a15":
                    called = out[cfg["template"]]["q"][i] < 0.05
                    kk = key + jit + np.where(called, 1e6, 0.0)     # called block first
                    rank_order = np.argsort(-kk, kind="stable")
                else:
                    rank_order = np.argsort(-(key + jit), kind="stable")
                pos = np.empty(G, np.int64)
                pos[rank_order] = np.arange(G)
                zt = zabs_sorted[pos]
                if cfg["K"] is not None:
                    zt = np.where(pos < int(cfg["K"]), zt, 0.0)
                sgn = np.sign(Hd["H"][hpos[t]][keep])
                sgn[sgn == 0] = 1.0
                zt_signed = sgn * zt
                dsg = eme.design(Mfull[keep], zt_signed, mode=cfg["mode"], T=T)
                # sub-gate block on a15's tail law, then place + repair
                prof_sub = np.zeros(G_full)
                prof_sub[~keep] = profile[~keep]
                v2, c2, p2 = emit_counts(prof_sub, np.full(n, EMIT_L), 0.1, rng, mode="quasi",
                                         spread=1.0, overdispersion=PHI)
                B_sub = to_block(v2, c2, p2, G_full)
                assert B_sub[:, keep].sum() == 0
                B, dg = eme.place(dsg, gidx, B_sub, Tc, rng)
                assert np.array_equal(B.sum(1), Tc), name
                zdes[name][i] = dsg["z"]
                ztgt[name][i] = zt_signed
                dg.update(clip_frac=float(np.mean(np.abs(np.abs(dsg["z"]) - zt) > 0.5)),
                          sign_flip_frac=float(np.mean((np.sign(dsg["z"]) != sgn) & (zt > 0.5))))
            else:  # resample law, 177's realisation
                cidx = rng.choice(ci, n, replace=False)
                Cc = rd.rows(cidx, G_full)
                dlib = Cc.sum(1)
                dlib[dlib == 0] = 1.0
                Ld = float(int(min(np.median(dlib), 1_000_000)))
                pp = (Mfull + 0.5) / (n * Ld)             # floor(n Ld pp) == Mfull exactly
                sv = np.full(G_full, 1.0)
                k8 = min(8000, int(gated.sum()))
                sv[order[:k8]] = 0.25
                pmask = np.zeros(G_full, dtype=bool)
                kk = min(int(cfg["K"]), int(gated.sum()))
                if cfg["head_param"]:
                    pmask[order[:kk]] = True               # a21: head parametric, tail resampled
                else:
                    pmask[:] = True
                    pmask[order[:kk]] = False              # a22: head resampled, rest parametric
                vals, cols, ptr = emit_counts(pp, dlib, 0.1, rng, mode="resample_mix", spread=sv,
                                              overdispersion=PHI, control_counts=Cc,
                                              control_profile=ctrl_rel, param_mask=pmask)
                B = to_block(vals, cols, ptr, G_full)
                svec[name][i] = sv[keep]
                dg = dict(param_frac=float(pmask[keep].mean()))
                del Cc
            delivered = B.sum(0).astype(np.int64)
            assert np.array_equal(delivered, Mfull), (name, int(np.abs(delivered - Mfull).sum()))
            digest[name].append(hashlib.sha256(delivered.tobytes()).hexdigest()[:16])
            zz, p, q, lfc = H.de_from_block(B, keep, cs, tie, ref_lin)
            out[name]["z"][i], out[name]["p"][i], out[name]["q"][i], out[name]["l"][i] = zz, p, q, lfc
            npred[name].append(int((q < 0.05).sum()))
            if law == "eme":
                dg["z_realised_vs_designed_max"] = float(np.abs(zz - dsg["z"]).max())
            dg["cell_total_sd"] = float(B.sum(1).std())
            diag[name].append(dg)
            del B
        if (i + 1) % 5 == 0 or i == 0:
            print(f"  {i+1}/{T_n} [{time.time()-t0:.0f}s]  n_pred " +
                  " ".join(f"{a[0]}={npred[a[0]][-1]}" for a in ARMS), flush=True)

    print("\nPSEUDOBULK BIT-IDENTITY (sha256 of the delivered per-gene panel total, per target)")
    ref = digest[ARMS[0][0]]
    for name, *_ in ARMS:
        same = sum(a == b for a, b in zip(digest[name], ref))
        print(f"  {name:<14} identical to {ARMS[0][0]} on {same}/{T_n} targets   first digest {digest[name][0]}")
    print()
    for name, law, cfg, seed in ARMS:
        path = f"{OUT}.{name}.npz"
        extra = {}
        if law == "eme":
            extra = dict(z_design=zdes[name], z_target=ztgt[name])
        np.savez_compressed(
            path, genes=genes_full[keep], perts=np.array(targets),
            real_p=real["p"], real_padj=real["q"], real_lfc=real["l"],
            pred_p=out[name]["p"], pred_padj=out[name]["q"], pred_lfc=out[name]["l"],
            pred_z=out[name]["z"], spread=svec[name], absH=absH_all, digest=np.array(digest[name]),
            **extra,
        )
        d = diag[name]
        summ = {k: float(np.median([x[k] for x in d])) for k in d[0]} if d else {}
        print(f"wrote {path}   median n_pred {np.median(npred[name]):.0f}   diag(median) {summ}",
              flush=True)
    with open(os.path.join(H.WORK, "realize185_diag.json"), "w") as fh:
        json.dump(dict(npred=npred, digest=digest, diag=diag), fh)
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
