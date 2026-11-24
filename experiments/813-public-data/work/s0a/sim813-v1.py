"""813 stage 0a: KA-SIM. Sets the KA-SG bounds, decides A vs B, and gives the power verdict, BEFORE any real TeloHAEC
expression row is read (PROPOSAL-2 section 2; VERDICT-2 G1/G5; C3). Runs on EC2 through labcloud.sh only.

Reads from the real IGVF file (IGVFFI2522NCFR, resting rep 1) ONLY metadata, by byte-capped HTTP range reads:
uns/guide_names, obsm/guide_assignment (indptr, indices), obs n_guides_total, total_counts, pool codes, and the cell
barcodes (for the seeded half split). X is never opened.

Model (fixed in PROPOSAL-2):
  base p  = the TeloHAEC PILOT's control profile (data/external/igvf-telohaec/b_ctrl.npy, TS-5e4 log1p space), on a
            3,000-gene axis: pilot CPM > 5, in K562-GWPS's genes, on the official axis, no panel genes (seed 813);
            the rest of the transcriptome is one remainder column that no perturbation touches.
  truth   = K562-GWPS donor rows (data/donors/k562_gwps) for the library targets K562 carries, zero for the rest; each row
            becomes a multiplicative effect m = expm1(b + delta) / expm1(b) (floored at 1e-6) so that a single-guide
            experiment measures delta* = log1p(TS p m / Z_t) - log1p(TS p), Z_t = sum(p m) + remainder (the analytic
            target every estimator is compared with).
  cells   = the real assignment and real total_counts; a cell's expected counts are L_c * p * prod(m) / Z_c.
  designs = (1) the real multi-guide design; (2) its single-guide twin (one random guide kept per cell, seed 813);
            (3) design 1 with an injected guide-load offset: exp(k_c * tau), tau a fixed random direction (seed 813) of
            norm 0.5 x median ||delta*||, k_c = n_guides_total.
  noise   = Poisson draws of every group sum (sums of independent Poissons are Poisson), R replicates.

Estimators (rows in TS-5e4 log1p space, kasg813.row):
  A      stratum-matched background, strata s = (min(n_guides_total, 12), pool); P_bg = sum_s w_s (S_s - C_ts).
  B      first-order co-occurrence deconvolution of A: B_t = A_t - sum_{u != t} (n_tu / n_t) (A~_u - mean A~), with A~ the
         stratum-matched rows of all library targets (stratum totals not carrier-subtracted). DEVIATION from PROPOSAL-2's
         full conjugate-gradient B: a full CG solve per replicate does not fit stage 0a's r7i.xlarge / $0.5 budget; the
         first Neumann term is the part of B that removes the co-occurrence fluctuation. Recorded in the freeze.
  S      KA-SG's single-guide rows: carriers with n_guides == 1 vs non-carriers with n_guides == 1, pool-matched.
  naive  carriers with n_guides >= 2 vs ALL zero-guide cells (the power row's biased estimator).

Usage (EC2): uv run python experiments/813-public-data/sim813.py --out experiments/813-public-data/work/s0a [--smoke]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import kasg813 as K  # noqa: E402

URL = "https://api.data.igvf.org/matrix-files/IGVFFI2522NCFR/@@download/IGVFFI2522NCFR.h5ad"
TS = 5e4
SEED = 813
KCAP = 12
G_AXIS = 3000
BLOCK = 40000
t0 = time.time()


def log(*a):
    print(f"[sim813 {time.time() - t0:7.1f}s]", *a, flush=True)


def sha_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def sha_arr(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ metadata (range reads)
def read_metadata(cap_mb: float):
    src = (HERE / "work" / "h5peek813.py").read_text().split("rf = RangeFile(URL)")[0]
    ns = {"__name__": "x"}
    argv = sys.argv
    sys.argv = ["x", "u", "o", str(cap_mb)]
    exec(src, ns)
    sys.argv = argv
    rf = ns["RangeFile"](URL)
    with h5py.File(rf, "r") as f:
        names = [x.decode() if isinstance(x, bytes) else str(x) for x in f["uns/guide_names"][:]]
        ga = f["obsm/guide_assignment"]
        indptr = ga["indptr"][:].astype(np.int64)
        indices = ga["indices"][:].astype(np.int64)
        n_g = f["obs/n_guides_total"][:].astype(np.int64)
        L = f["obs/total_counts"][:].astype(np.float64)
        pool = f["obs/pool/codes"][:].astype(np.int64)
        n_pool = int(f["obs/pool/categories"].shape[0])
        halves_rule = "sha1(b'813:' + barcode)[0] & 1"
        try:
            bc = f["obs/_index"][:]
            half = np.array([hashlib.sha1(b"813:" + (b if isinstance(b, bytes) else str(b).encode())).digest()[0] & 1
                             for b in bc], np.int64)
        except Exception:                          # byte cap (h5py may wrap it): a seeded hash of the row index instead
            halves_rule = "sha1(b'813:' + str(row_index))[0] & 1 (barcodes over the byte cap)"
            half = np.array([hashlib.sha1(b"813:" + str(i).encode()).digest()[0] & 1 for i in range(len(n_g))], np.int64)
    return dict(names=names, indptr=indptr, indices=indices, n_g=n_g, L=L, pool=pool, n_pool=n_pool, half=half,
                halves_rule=halves_rule, bytes=rf.fetched, size=rf.size)


def target_of(n: str) -> str:
    if n.startswith("non-targeting"):
        return "NTC"
    p = n.split("_")
    return p[1] if len(p) >= 3 else n


# ------------------------------------------------------------------------------------------------ the model
class Model:
    def __init__(self, md: dict, smoke: bool):
        rng = np.random.default_rng(SEED)
        gt = np.array([target_of(n) for n in md["names"]])
        lib = sorted(set(gt) - {"NTC"})
        tix = {t: i for i, t in enumerate(lib)}
        g2t = np.array([tix.get(t, -1) for t in gt])
        n_cells = len(md["n_g"])
        rows = np.repeat(np.arange(n_cells), np.diff(md["indptr"]))
        cols = g2t[md["indices"]]
        keep = cols >= 0
        M = sp.csr_matrix((np.ones(keep.sum(), np.float64), (rows[keep], cols[keep])), shape=(n_cells, len(lib)))
        M.sum_duplicates()
        M.data[:] = 1.0                                         # two guides of one target in a cell count once
        self.lib, self.tix = lib, tix
        self.M = M
        self.n_g, self.L, self.pool, self.half, self.n_pool = md["n_g"], md["L"], md["pool"], md["half"], md["n_pool"]
        if smoke:                                               # a 60k-cell subsample, to fail fast before the full run
            sel = np.sort(rng.choice(n_cells, 60000, replace=False))
            self.M, self.n_g, self.L, self.pool, self.half = M[sel], self.n_g[sel], self.L[sel], self.pool[sel], self.half[sel]
        self.n = self.M.shape[0]
        self.inc = np.asarray(self.M.sum(0)).ravel()            # incidences per library target

        # panel and stage-0 targets
        panel = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))]
        self.focal = sorted(t for t in panel if t in tix)
        elig = sorted(t for t in self.focal if self.inc[tix[t]] >= (100 if not smoke else 1))
        self.stage0 = sorted(np.random.default_rng(SEED).choice(elig, 100, replace=False).tolist())
        self.fidx = np.array([tix[t] for t in self.focal])

        # gene axis and base profile
        off = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
        pil_g = (ROOT / "data/external/igvf-telohaec/genes.txt").read_text().split()
        pil_b = np.load(ROOT / "data/external/igvf-telohaec/b_ctrl.npy").astype(np.float64)
        k_g = (ROOT / "data/donors/k562_gwps/genes.txt").read_text().split()
        k_p = (ROOT / "data/donors/k562_gwps/perts.txt").read_text().split()
        p_full = np.expm1(pil_b) / TS
        p_full = p_full / p_full.sum()
        cpm = p_full * 1e6
        kset, offset, pset = set(k_g), set(off), set(panel)
        cand = sorted(g for i, g in enumerate(pil_g) if cpm[i] > 5 and g in kset and g in offset and g not in pset)
        self.axis = sorted(np.random.default_rng(SEED).choice(cand, G_AXIS, replace=False).tolist())
        pi = {g: i for i, g in enumerate(pil_g)}
        ki = {g: i for i, g in enumerate(k_g)}
        p = p_full[[pi[g] for g in self.axis]]
        self.p, self.R = p, 1.0 - p.sum()
        b = np.log1p(TS * p)
        delta = np.load(ROOT / "data/donors/k562_gwps/delta.npy", mmap_mode="r")
        kp = {t: i for i, t in enumerate(k_p)}
        cols = np.array([ki[g] for g in self.axis])
        self.has_truth = np.array([t in kp for t in lib])
        logm = np.zeros((len(lib), G_AXIS), np.float64)
        for i, t in enumerate(lib):
            if t in kp:
                d = np.asarray(delta[kp[t], cols], np.float64)
                m = np.expm1(b + d) / np.expm1(b)
                logm[i] = np.log(np.maximum(m, 1e-6))
        self.logm = logm
        # the analytic single-guide target delta*
        w = p[None, :] * np.exp(logm[self.fidx])
        Zt = w.sum(1) + self.R
        self.dstar = np.log1p(TS * w / Zt[:, None]) - np.log1p(TS * p)[None, :]
        self.truth_f = self.has_truth[self.fidx]
        nrm = np.linalg.norm(self.dstar[self.truth_f], axis=1)
        self.med_norm = float(np.median(nrm))
        tau = np.random.default_rng(SEED + 1).normal(size=G_AXIS)
        self.tau = tau / np.linalg.norm(tau) * 0.5 * self.med_norm
        # strata
        self.n_s = (KCAP + 1) * self.n_pool
        self.kbin = np.minimum(self.n_g, KCAP)
        log(f"model: {self.n} cells, {len(lib)} library targets ({int(self.has_truth.sum())} with K562 truth), "
            f"{len(self.focal)} focal (panel) targets ({int(self.truth_f.sum())} with truth), axis {G_AXIS} genes, "
            f"remainder mass {self.R:.3f}, median ||delta*|| {self.med_norm:.3f}")

    # one design's expected group sums ------------------------------------------------------------
    def expected(self, design: str) -> dict:
        n, G1 = self.n, G_AXIS + 1
        if design == "twin":
            rng = np.random.default_rng(SEED + 2)
            Mc = self.M.tocsr()
            ip, ind = Mc.indptr, Mc.indices
            keep_rows, keep_cols = [], []
            for c in range(n):
                a, b = ip[c], ip[c + 1]
                if b > a:
                    keep_rows.append(c)
                    keep_cols.append(ind[a + rng.integers(b - a)])
            M = sp.csr_matrix((np.ones(len(keep_rows)), (keep_rows, keep_cols)), shape=Mc.shape)
            kb = np.asarray(M.sum(1)).ravel().astype(np.int64)              # 0 or 1
            kg = kb
        else:
            M, kb, kg = self.M, self.kbin, self.n_g
        s = kb * self.n_pool + self.pool
        sh = s * 2 + self.half
        n_sh = self.n_s * 2
        # focal (target, stratum, half) pairs that exist -> compact index
        Mf = M[:, self.fidx].tocoo()
        key = Mf.col.astype(np.int64) * n_sh + sh[Mf.row]
        ukey, inv = np.unique(key, return_inverse=True)
        Fm = sp.csr_matrix((np.ones(len(key)), (Mf.row, inv)), shape=(n, len(ukey)))
        H = sp.csr_matrix((np.ones(n), (np.arange(n), sh)), shape=(n, n_sh))
        E_f = np.zeros((len(ukey), G1))
        E_sh = np.zeros((n_sh, G1))
        E_all = np.zeros((len(self.lib), G1)) if design == "multi" else None
        for a in range(0, n, BLOCK):
            b = min(n, a + BLOCK)
            Mb = M[a:b]
            ls = np.asarray(Mb @ self.logm)
            if design == "inject":
                ls = ls + kg[a:b, None] * self.tau[None, :]
            w = self.p[None, :] * np.exp(ls)
            Z = w.sum(1) + self.R
            rate = np.empty((b - a, G1))
            rate[:, :G_AXIS] = self.L[a:b, None] * w / Z[:, None]
            rate[:, G_AXIS] = self.L[a:b] * self.R / Z
            E_f += np.asarray(Fm[a:b].T @ rate)
            E_sh += np.asarray(H[a:b].T @ rate)
            if E_all is not None:
                E_all += np.asarray(Mb.T @ rate)
        # cell counts per (focal, stratum, half) and per stratum-half; all-target stratum composition (for B)
        n_f = np.bincount(inv, minlength=len(ukey)).astype(np.float64)
        N_sh = np.bincount(sh, minlength=n_sh).astype(np.float64)
        out = dict(design=design, ukey=ukey, E_f=E_f, n_f=n_f, E_sh=E_sh, N_sh=N_sh, n_sh=n_sh, kb=kb)
        if E_all is not None:
            Sx = sp.csr_matrix((np.ones(n), (np.arange(n), s)), shape=(n, self.n_s))
            out["E_all"] = E_all
            out["comp_all"] = np.asarray((M.T @ Sx).todense())                 # lib x strata carrier counts
            Mfc = M[:, self.fidx].tocsc()
            co = (Mfc.T @ M).tocsr().astype(np.float64)                         # focal x lib co-occurrence counts
            co[np.arange(len(self.fidx)), self.fidx] = 0
            co.eliminate_zeros()
            out["co"] = co
            out["n_t"] = self.inc[self.fidx].astype(np.float64)
        tk = ukey // n_sh
        order = np.argsort(tk, kind="stable")
        bounds_ = np.searchsorted(tk[order], np.arange(len(self.fidx) + 1))
        out["pairs"] = [order[bounds_[i]:bounds_[i + 1]] for i in range(len(self.fidx))]   # pair rows per focal target
        return out


# ------------------------------------------------------------------------------------------------ rows from group sums
def rows_for(model: Model, ex: dict, Cf: np.ndarray, Ssh: np.ndarray, targets: list, kset: str, halves: bool):
    """Estimator-A-type rows: carriers in strata of k-class kset ('all' k>=1, 'ge2', 'eq1') vs stratum-matched
    non-carriers (matched per (stratum, half) cell group; P_bg = sum w (S - C), w = n_carriers / n_noncarriers).
    Returns {'full': (n_t, G), 'h0': ..., 'h1': ...} (halves only if asked); NaN rows where a half has no carriers."""
    n_sh, npool = ex["n_sh"], model.n_pool
    kk = (np.arange(model.n_s) // npool)
    ok_s = {"all": kk >= 1, "ge2": kk >= 2, "eq1": kk == 1}[kset]
    pos = {t: i for i, t in enumerate(model.focal)}
    names = ["full", "h0", "h1"] if halves else ["full"]
    res = {k: np.full((len(targets), G_AXIS), np.nan) for k in names}
    for j, t in enumerate(targets):
        sel = ex["pairs"][pos[t]]
        shk = ex["ukey"][sel] % n_sh
        s_of, h_of = shk // 2, shk % 2
        m = ok_s[s_of]
        sel, shk, h_of = sel[m], shk[m], h_of[m]
        if len(sel) == 0:
            continue
        nonc = ex["N_sh"][shk] - ex["n_f"][sel]
        w = np.where(nonc > 0, ex["n_f"][sel] / np.maximum(nonc, 1.0), 0.0)
        Cs, Ss = Cf[sel], Ssh[shk]
        for name in names:
            mm = np.ones(len(sel), bool) if name == "full" else (h_of == int(name[1]))
            if not mm.any():
                continue
            P = Cs[mm].sum(0)
            bg = w[mm] @ (Ss[mm] - Cs[mm])
            res[name][j] = K.row(P, bg)[:G_AXIS]
    return res


def naive_rows(model: Model, ex: dict, Cf, Ssh, targets):
    """carriers with k >= 2 vs ALL zero-guide cells (no pool match), per half and full."""
    n_sh, npool = ex["n_sh"], model.n_pool
    kk = (np.arange(model.n_s) // npool)
    pos = {t: i for i, t in enumerate(model.focal)}
    zero = {h: Ssh[[s * 2 + h for s in range(model.n_s) if kk[s] == 0]].sum(0) for h in (0, 1)}
    res = {k: np.full((len(targets), G_AXIS), np.nan) for k in ("full", "h0", "h1")}
    for j, t in enumerate(targets):
        sel = ex["pairs"][pos[t]]
        shk = ex["ukey"][sel] % n_sh
        s_of, h_of = shk // 2, shk % 2
        m = kk[s_of] >= 2
        sel, h_of = sel[m], h_of[m]
        for name, hs in (("full", (0, 1)), ("h0", (0,)), ("h1", (1,))):
            mm = np.isin(h_of, hs)
            if not mm.any():
                continue
            P = Cf[sel[mm]].sum(0)
            bg = sum(zero[h] for h in hs)
            res[name][j] = K.row(P, bg)[:G_AXIS]
    return res


def cosines(X, Y):
    nx, ny = np.linalg.norm(X, axis=1), np.linalg.norm(Y, axis=1)
    return (X * Y).sum(1) / np.maximum(nx * ny, 1e-300)


# ------------------------------------------------------------------------------------------------ replicates
_G: dict = {}


def rep_multi(r: int) -> dict:
    m, ex = _G["model"], _G["multi"]
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 1, r]))
    Cf, Ssh = rng.poisson(ex["E_f"]).astype(np.float64), rng.poisson(ex["E_sh"]).astype(np.float64)
    A = rows_for(m, ex, Cf, Ssh, m.focal, "all", False)["full"]
    tr = m.truth_f
    out = dict(cosA=float(np.nanmedian(cosines(A[tr], m.dstar[tr]))),
               offA=float(np.linalg.norm(np.nanmean(A - m.dstar, 0)) / m.med_norm))
    # KA-SG statistics on the stage-0 targets
    Ag = rows_for(m, ex, Cf, Ssh, m.stage0, "ge2", True)
    Sg = rows_for(m, ex, Cf, Ssh, m.stage0, "eq1", True)
    ok = ~(np.isnan(Ag["h0"]).any(1) | np.isnan(Ag["h1"]).any(1) | np.isnan(Sg["h0"]).any(1) | np.isnan(Sg["h1"]).any(1))
    out["beta"] = K.beta(Ag["full"][ok], Sg["full"][ok], Sg["h0"][ok], Sg["h1"][ok])
    out["dsh"] = K.dsh(Ag["h0"][ok], Ag["h1"][ok], Sg["h0"][ok], Sg["h1"][ok])
    out["n_sg_targets"] = int(ok.sum())
    if r < _G["R_B"]:                                              # estimator B (first-order deconvolution)
        Call = rng.poisson(ex["E_all"]).astype(np.float64)
        S_full = Ssh.reshape(m.n_s, 2, -1).sum(1)
        N_full = ex["N_sh"].reshape(m.n_s, 2).sum(1)
        comp = ex["comp_all"]
        W = comp / np.maximum(N_full[None, :] - comp, 1.0)
        Pbg = W @ S_full
        At = (np.log1p(TS * Call / np.maximum(Call.sum(1, keepdims=True), 1.0))
              - np.log1p(TS * Pbg / np.maximum(Pbg.sum(1, keepdims=True), 1e-300)))
        At = At[:, :G_AXIS]
        good = (ex["n_t_all"] > 0) & (Call.sum(1) > 0) & np.isfinite(At).all(1)
        mean_t = (At[good] * ex["n_t_all"][good, None]).sum(0) / ex["n_t_all"][good].sum()
        corr = ex["co"] @ (np.where(good[:, None], At - mean_t[None, :], 0.0))
        B = A - corr / ex["n_t"][:, None]
        out["cosB"] = float(np.nanmedian(cosines(B[tr], m.dstar[tr])))
        out["offB"] = float(np.linalg.norm(np.nanmean(B - m.dstar, 0)) / m.med_norm)
    return out


def rep_inject(r: int) -> dict:
    m, ex = _G["model"], _G["inject"]
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 3, r]))
    Cf, Ssh = rng.poisson(ex["E_f"]).astype(np.float64), rng.poisson(ex["E_sh"]).astype(np.float64)
    A = rows_for(m, ex, Cf, Ssh, m.focal, "all", False)["full"]
    out = dict(offA=float(np.linalg.norm(np.nanmean(A - m.dstar, 0)) / m.med_norm))
    Ag = rows_for(m, ex, Cf, Ssh, m.stage0, "ge2", True)
    Sg = rows_for(m, ex, Cf, Ssh, m.stage0, "eq1", True)
    Nv = naive_rows(m, ex, Cf, Ssh, m.stage0)
    ok = ~(np.isnan(Ag["h0"]).any(1) | np.isnan(Ag["h1"]).any(1) | np.isnan(Sg["h0"]).any(1) | np.isnan(Sg["h1"]).any(1))
    out["beta_A"] = K.beta(Ag["full"][ok], Sg["full"][ok], Sg["h0"][ok], Sg["h1"][ok])
    out["dsh_A"] = K.dsh(Ag["h0"][ok], Ag["h1"][ok], Sg["h0"][ok], Sg["h1"][ok])
    out["beta_naive"] = K.beta(Nv["full"][ok], Sg["full"][ok], Sg["h0"][ok], Sg["h1"][ok])
    out["dsh_naive"] = K.dsh(Nv["h0"][ok], Nv["h1"][ok], Sg["h0"][ok], Sg["h1"][ok])
    return out


def rep_twin(r: int) -> dict:
    m, ex = _G["model"], _G["twin"]
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 2, r]))
    Cf, Ssh = rng.poisson(ex["E_f"]).astype(np.float64), rng.poisson(ex["E_sh"]).astype(np.float64)
    A = rows_for(m, ex, Cf, Ssh, m.focal, "all", False)["full"]
    tr = m.truth_f
    return dict(cosT=float(np.nanmedian(cosines(A[tr], m.dstar[tr]))))


def run_pool(fn, R, procs):
    with mp.get_context("fork").Pool(procs) as P:
        return P.map(fn, range(R), chunksize=max(1, R // (procs * 4)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--procs", type=int, default=os.cpu_count() or 4)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    R, R_T, R_B = (200, 50, 50) if not a.smoke else (4, 2, 2)
    md = read_metadata(cap_mb=400)
    log(f"metadata: {len(md['n_g'])} cells, {len(md['names'])} guides, {md['bytes'] / 2**20:.1f} MB fetched of "
        f"{md['size'] / 1e9:.1f} GB; halves rule: {md['halves_rule']}")
    m = Model(md, a.smoke)
    _G.update(model=m, R_B=R_B)

    # noise-free code check: twin design, single-guide carriers vs zero-guide cells, expected sums -> delta*
    tw = m.expected("twin")
    _G["twin"] = tw
    kk = np.arange(m.n_s) // m.n_pool
    zero = tw["E_sh"][[s * 2 + h for s in range(m.n_s) for h in (0, 1) if kk[s] == 0]].sum(0)
    nf = np.full((len(m.focal), G_AXIS), np.nan)
    for j in range(len(m.focal)):
        sel = tw["pairs"][j]
        if len(sel):
            nf[j] = K.row(tw["E_f"][sel].sum(0), zero)[:G_AXIS]
    tr = m.truth_f
    cos_nf = cosines(nf[tr], m.dstar[tr])
    log(f"noise-free twin check: min cos {cos_nf.min():.8f}, median {np.median(cos_nf):.8f}")

    mu = m.expected("multi")
    mu["n_t_all"] = m.inc.astype(np.float64)
    _G["multi"] = mu
    log("expected sums: twin + multi done")
    inj = m.expected("inject")
    _G["inject"] = inj
    log("expected sums: inject done")

    res_m = run_pool(rep_multi, R, a.procs)
    log("multi replicates done")
    res_i = run_pool(rep_inject, R, a.procs)
    log("inject replicates done")
    res_t = run_pool(rep_twin, R_T, a.procs)
    log("twin replicates done")

    col = lambda rs, k: np.array([x[k] for x in rs if k in x], float)   # noqa: E731
    beta_lo, beta_hi = K.bounds(col(res_m, "beta"))
    dsh_lo, dsh_hi = K.bounds(col(res_m, "dsh"))
    bn, dn = col(res_i, "beta_naive"), col(res_i, "dsh_naive")
    out_n = ((bn < beta_lo) | (bn > beta_hi) | (dn < dsh_lo) | (dn > dsh_hi))
    power = float(out_n.mean())
    bA, dA = col(res_i, "beta_A"), col(res_i, "dsh_A")
    in_A_inj = float((~((bA < beta_lo) | (bA > beta_hi) | (dA < dsh_lo) | (dA > dsh_hi))).mean())
    cosA, cosT, cosB = np.median(col(res_m, "cosA")), np.median(col(res_t, "cosT")), np.median(col(res_m, "cosB"))
    offA1, offA3 = np.median(col(res_m, "offA")), np.median(col(res_i, "offA"))
    ka = {"i_cos_ratio_A_vs_twin": float(cosA / cosT), "i_pass": bool(cosA / cosT >= 0.9),
          "ii_offset_design1": float(offA1), "ii_offset_design3": float(offA3),
          "ii_pass": bool(offA1 <= 0.10 and offA3 <= 0.10),
          "iii_noise_free_min_cos": float(cos_nf.min()), "iii_pass": bool(cos_nf.min() > 0.9999)}
    ka["pass"] = ka["i_pass"] and ka["ii_pass"] and ka["iii_pass"]
    choice = "B" if (cosB - cosA) >= 0.02 else "A"
    freeze = {
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "smoke": a.smoke,
        "estimator": choice, "estimator_B_form": "first-order co-occurrence deconvolution (deviation from CG; see sim813 doc)",
        "kasg": {"beta_bounds": [beta_lo, beta_hi], "dsh_bounds": [dsh_lo, dsh_hi], "R": R,
                 "power_naive_outside": power, "informative": bool(power >= 0.8),
                 "A_inside_bounds_under_injection": in_A_inj},
        "ka_sim": ka,
        "medians": {"cosA": float(cosA), "cosB": float(cosB), "cos_twin": float(cosT), "offB": float(np.median(col(res_m, "offB")))},
        "R": {"multi": R, "inject": R, "twin": R_T, "B": R_B},
        "stage0_targets": m.stage0, "stage0_targets_sha256": hashlib.sha256("\n".join(m.stage0).encode()).hexdigest(),
        "focal_targets_n": len(m.focal), "focal_with_truth": int(m.truth_f.sum()),
        "axis_sha256": hashlib.sha256("\n".join(m.axis).encode()).hexdigest(), "axis_n": G_AXIS,
        "strata": {"kcap": KCAP, "n_pool": m.n_pool, "rule": "(min(n_guides_total,12), pool)"},
        "halves_rule": md["halves_rule"], "seed": SEED,
        "inputs": {"igvf_file": "IGVFFI2522NCFR (IGVFDS2068OODO, resting rep 1)", "igvf_size_bytes": md["size"],
                   "metadata_bytes_fetched": md["bytes"],
                   "guide_names_sha256": hashlib.sha256("\n".join(md["names"]).encode()).hexdigest(),
                   "assignment_indices_sha256": sha_arr(md["indices"]),
                   "pilot_b_ctrl_sha256": sha_file(ROOT / "data/external/igvf-telohaec/b_ctrl.npy"),
                   "k562_delta_sha256": sha_file(ROOT / "data/donors/k562_gwps/delta.npy")},
        "code": {"sim813.py": sha_file(Path(__file__)), "kasg813.py": sha_file(HERE / "kasg813.py"),
                 "h5peek813.py": sha_file(HERE / "work/h5peek813.py")},
        "model": {"cells": m.n, "library_targets": len(m.lib), "with_truth": int(m.has_truth.sum()),
                  "remainder_mass": m.R, "median_dstar_norm": m.med_norm, "tau_norm": float(np.linalg.norm(m.tau))},
    }
    name = "freeze813.json" if not a.smoke else "freeze813-smoke.json"
    (out / name).write_text(json.dumps(freeze, indent=1))
    (out / ("detail813.json" if not a.smoke else "detail813-smoke.json")).write_text(
        json.dumps({"multi": res_m, "inject": res_i, "twin": res_t, "noise_free_cos": cos_nf.tolist()}, indent=0))
    log(f"wrote {out / name} sha256 {sha_file(out / name)}")
    log(json.dumps({k: freeze[k] for k in ("estimator", "kasg", "ka_sim", "medians")}, indent=1))


if __name__ == "__main__":
    main()
