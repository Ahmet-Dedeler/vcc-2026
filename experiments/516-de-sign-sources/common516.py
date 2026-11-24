"""516 shared helpers: H1 surface, production signs, pricing a candidate sign source, donor loaders.

Surfaces
  H1     h1-test truth (323's reference.npz, keep axis) and the production emitter's H1 stats (435 coarse MDE) — the
         four DE members reproduce walk513's known answers exactly.
  folds  480's cancer-fold gates (Orion HCT116 / HEK293T thinned harnesses, K562-GWPS panel cells; 150 cells per
         target; relative direction only, 480 section 4).

Pricing a candidate sign source S (targets x genes; +1 / -1 where it has an opinion, 0 or NaN where it abstains):
  the hybrid sign H = S where S != 0, else production's sign. Every member keeps production's magnitudes, calls and order:
  (i)   sign accuracy on the recipient's DE pool (reference q < 0.05, finite non-zero lfc, target excluded), pooled per
        pool band (< 20, 20-100, >= 100) with a target-bootstrap 95 % CI of H - production;
  (ii)  nmae with pred = H * |our lfc| on the scorer's gate (>= 10 gate genes);
  (iii) reach with H substituted in our order (walk513's order and walk, BRIEF section 5);
  (iv)  fid with H on our calls (calls unchanged).
Nothing here writes outside experiments/516-de-sign-sources/. Run with PYTHONDONTWRITEBYTECODE=1 (imports 513/480 code).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/513-reach-walk-anatomy"))
sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument"))
import walk513  # noqa: E402
from common480 import gate_from, load_donor_rows, nmae_rows  # noqa: E402

ALPHA, P0 = 0.05, 0.9
MET435 = ROOT / "experiments/435-zs-on-coarse/work/h1/metrics-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.json"
KA = dict(fid=0.5089291277653674, jac=0.07469853420125332, nmae=0.9427040517722616, reach=0.26925697759272765)
BANDS = [(0, 20, "<20"), (20, 100, "20-100"), (100, 10 ** 9, ">=100")]
# panel-available donors (cover A/B/C targets); essential-only screens cover 0/300 A/B/C targets (430)
PANEL_DONORS = ["k562_gwps", "orion_hct116", "orion_hek293t"]
ESS_DONORS = ["k562_essential", "rpe1", "nadig_jurkat"]
CD4_H5AD = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"


# ----------------------------------------------------------------------------------------------------------------- H1
def load_h1():
    """walk513's loader plus the pool mask, n_real, band and nmae gate."""
    d = walk513.load()
    G, T = d["genes"], d["targets"]
    own = G[None, :] == T[:, None]
    d["own"] = own
    d["rs"] = (d["rq"] < ALPHA) & ~own                                        # reference significant (reach pool)
    d["pool"] = d["rs"] & np.isfinite(d["rl"]) & (d["rl"] != 0)              # adjudicable pool genes
    d["n_real"] = d["rs"].sum(1)
    d["band"] = np.select([d["n_real"] < 20, d["n_real"] < 100], [0, 1], 2)
    d["gateT"] = gate_from(d["rq"], d["rl"], G, T)                             # nmae gate (>= 10 genes)
    d["truth"] = np.sign(np.nan_to_num(d["rl"]))
    ol = d["ol"]
    d["prod"] = np.where(np.isnan(ol), 0.0, np.sign(ol))                     # +-inf keep their sign (walk513)
    d["nz"] = ~np.isnan(ol) & (ol != 0)                                        # a sign the walk/fid can read
    d["absl"] = np.abs(np.nan_to_num(ol, nan=0.0, posinf=0.0, neginf=0.0))    # nmae magnitude (walk513: inf -> 0)
    d["calls"] = (d["oq"] < ALPHA) & ~own
    d["orders"] = [walk513.order_keys(d["oq"][k], d["op"][k], d["ol"][k], G) for k in range(len(T))]
    return d


def members(d, H=None):
    """fid, jac, nmae, reach per target (arrays) with signs H (default production) at production magnitudes/order."""
    if H is None:
        H = d["prod"]
    T = d["targets"]
    fid = np.full(len(T), np.nan); jac = np.full(len(T), np.nan); reach = np.full(len(T), np.nan)
    for k in range(len(T)):
        rs, os_ = d["rs"][k], d["calls"][k]
        adj = np.isfinite(d["rl"][k]) & (d["rl"][k] != 0) & ~d["own"][k]
        n_real = int(rs.sum()); calls = os_ & adj; n_pred = int(calls.sum())
        kk = int(((H[k][calls] == d["truth"][k][calls]) & d["nz"][k][calls]).sum()) if n_pred else 0
        if n_real > 0 or n_pred > 0:
            fid[k] = kk / max(n_pred, n_real)
        u = (rs | os_).sum(); jac[k] = (rs & os_).sum() / u if u else 1.0
        if n_real:
            order = d["orders"][k]
            seq = order[rs[order]]
            seq = seq[np.isfinite(d["rl"][k][seq]) & (d["rl"][k][seq] != 0)]
            walk = ((H[k][seq] == d["truth"][k][seq]) & d["nz"][k][seq]).astype(float)
            if len(walk):
                cum = np.cumsum(walk) / np.arange(1, len(walk) + 1)
                ok = np.nonzero(cum >= P0)[0]
                reach[k] = (ok[-1] + 1 if len(ok) else 0) / n_real
            else:
                reach[k] = 0.0
    nm = nmae_rows(H * d["absl"], d["gateT"])
    return dict(fid=fid, jac=jac, nmae=nm, reach=reach)


def known_answers(d):
    m = members(d)
    got = {k: float(np.nanmean(v)) for k, v in m.items()}
    ok = all(abs(got[k] - KA[k]) < 1e-9 for k in KA)
    want = json.load(open(MET435))["raw"]
    ok &= all(abs(got[k] - want[k]) < 1e-9 for k in KA)
    return got, ok, m


# ----------------------------------------------------------------------------------------------------------- pricing
def hybrid(d, S):
    """Candidate sign where it has an opinion (+1/-1), production elsewhere."""
    S = np.nan_to_num(np.asarray(S, float))
    return np.where(S != 0, np.sign(S), d["prod"])


def boot_pooled(hits_a, hits_b, n, groups, B=4000, seed=516):
    """Target bootstrap of pooled accuracy difference. hits_*: per-target hit counts; n: per-target cell counts.
    groups: boolean mask of targets to include. Returns (acc_a, acc_b, delta, lo, hi)."""
    a, b, n = hits_a[groups], hits_b[groups], n[groups]
    if n.sum() == 0:
        return (np.nan,) * 5
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(n), size=(B, len(n)))
    na = n[idx].sum(1)
    dd = (b[idx].sum(1) - a[idx].sum(1)) / np.maximum(na, 1)
    return (a.sum() / n.sum(), b.sum() / n.sum(), (b.sum() - a.sum()) / n.sum(),
            float(np.quantile(dd, 0.025)), float(np.quantile(dd, 0.975)))


def boot_mean_delta(ra, rb, B=4000, seed=516):
    ok = np.isfinite(ra) & np.isfinite(rb)
    dlt = (rb - ra)[ok]
    rng = np.random.default_rng(seed)
    bs = dlt[rng.integers(0, len(dlt), size=(B, len(dlt)))].mean(1)
    return float(dlt.mean()), float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


def sign_table(truth, pool, band, A, S, cover=None, weight=None):
    """Pooled sign accuracy of A (baseline) and S (candidate) on pool cells by band, with target-bootstrap CI of S - A.
    Cells where A == 0 are dropped (production's zero lfc cannot score). cover: optional mask (targets x genes) to
    restrict to. weight: optional per-cell weights (e.g. |l|). Returns dict band -> tuple."""
    m = pool & (A != 0)
    if cover is not None:
        m &= cover
    w = np.ones_like(truth) if weight is None else weight
    ha = ((A == truth) & m) * w
    hb = ((S == truth) & m) * w                                    # S is a hybrid: never 0 where A != 0
    n = (m * w).sum(1)
    out = {}
    for lo, hi, name in [(0, 10 ** 9, "all")] + BANDS:
        g = (band >= 0) if name == "all" else np.isin(band, [i for i, b in enumerate(BANDS) if b[2] == name])
        out[name] = boot_pooled(ha.sum(1), hb.sum(1), n, g)
        out[name] = out[name] + (int((m & g[:, None]).sum()),)
    return out


def price(d, S, label, subsets=None, verbose=True, members_too=True):
    """Price candidate signs S on H1 against production, on each target subset (dict name -> bool mask over targets;
    default all 85). Returns a dict of results per subset."""
    H = hybrid(d, S)
    Sn = np.nan_to_num(np.asarray(S, float))
    cov = (Sn != 0)
    if subsets is None:
        subsets = {"all85": np.ones(len(d["targets"]), bool)}
    wl = np.abs(np.nan_to_num(d["rl"]))
    if members_too:
        base = d.get("_base_members")
        if base is None:
            base = d["_base_members"] = members(d)
        cand = members(d, H)
    out = {}
    for sn, sub in subsets.items():
        pool = d["pool"] & sub[:, None]
        res = dict(label=label, subset=sn, coverage_pool=float((cov & pool).sum() / pool.sum()),
                   flips_pool=int(((H != d["prod"]) & pool).sum()),
                   flips_all=int(((H != d["prod"]) & (d["prod"] != 0) & sub[:, None]).sum()))
        res["acc"] = sign_table(d["truth"], pool, d["band"], d["prod"], H)
        res["acc_cov"] = sign_table(d["truth"], pool, d["band"], d["prod"], H, cover=cov)
        res["acc_w"] = sign_table(d["truth"], pool, d["band"], d["prod"], H, weight=wl)
        # target-mean accuracy over targets with >= 10 pool cells (nmae-like weighting of targets)
        m = pool & (d["prod"] != 0)
        n = m.sum(1); ok = n >= 10
        ha = ((d["prod"] == d["truth"]) & m).sum(1) / np.maximum(n, 1)
        hb = ((H == d["truth"]) & m).sum(1) / np.maximum(n, 1)
        res["acc_tmean"] = (float(ha[ok].mean()), float(hb[ok].mean()),
                            *boot_mean_delta(np.where(ok, ha, np.nan), np.where(ok, hb, np.nan)), int(ok.sum()))
        if members_too:
            res["members"] = {}
            for k in ("fid", "nmae", "reach"):
                a = np.where(sub, base[k], np.nan); b = np.where(sub, cand[k], np.nan)
                res["members"][k] = (float(np.nanmean(a)), float(np.nanmean(b)), *boot_mean_delta(a, b))
        out[sn] = res
        if verbose:
            print_price(res)
    return out


def fmt_acc(t):
    a, b, dl, lo, hi, n = t
    return f"{a:.4f}->{b:.4f} {dl:+.4f} [{lo:+.4f},{hi:+.4f}] n{n}"


def print_price(res):
    print(f"== {res['label']} [{res['subset']}]: pool coverage {res['coverage_pool']:.3f}, pool flips "
          f"{res['flips_pool']}, all flips {res['flips_all']}")
    for key, nm in (("acc", "pool sign acc (hybrid)"), ("acc_cov", "pool sign acc (covered cells)"),
                    ("acc_w", "|l|-weighted pool acc")):
        print(f"   {nm:30s} " + " | ".join(f"{b}: {fmt_acc(res[key][b])}" for b in ("all", "<20", "20-100", ">=100")))
    t = res["acc_tmean"]
    print(f"   target-mean pool acc (>=10 cells): {t[0]:.4f}->{t[1]:.4f} {t[2]:+.4f} [{t[3]:+.4f},{t[4]:+.4f}] "
          f"n{t[5]}")
    if "members" in res:
        print("   members (prod -> cand, delta [95% target bootstrap]): " +
              " | ".join(f"{k} {v[0]:.5f}->{v[1]:.5f} {v[2]:+.5f} [{v[3]:+.5f},{v[4]:+.5f}]"
                         for k, v in res["members"].items()), flush=True)


def summary_line(res):
    """One compact line: pooled acc delta [CI], band deltas, target-mean delta, nmae/reach/fid deltas."""
    a = res["acc"]; t = res["acc_tmean"]; m = res.get("members", {})
    s = (f"{res['label'][:44]:44s} {res['subset']:6s} cov {res['coverage_pool']:.2f} "
         f"pool {a['all'][1]:.4f} {a['all'][2]:+.4f} [{a['all'][3]:+.4f},{a['all'][4]:+.4f}] "
         f"b<20 {a['<20'][2]:+.3f} b20 {a['20-100'][2]:+.3f} b100 {a['>=100'][2]:+.4f} | tmean {t[2]:+.4f} "
         f"[{t[3]:+.4f},{t[4]:+.4f}]")
    if m:
        s += (f" | nmae {m['nmae'][2]:+.5f} [{m['nmae'][3]:+.5f},{m['nmae'][4]:+.5f}] reach {m['reach'][2]:+.4f} "
              f"[{m['reach'][3]:+.4f},{m['reach'][4]:+.4f}] fid {m['fid'][2]:+.5f}")
    return s


# ------------------------------------------------------------------------------------------------------------ donors
def donor_rows(name, targets, genes, pc_cpm=1.0):
    """log2 pseudobulk ratio rows for a donor table in data/donors (NaN where absent), plus donor ctrl CPM."""
    L2, has, okc, cpm = load_donor_rows(name, targets, genes, pc_cpm=pc_cpm)
    return L2, has, cpm


def cd4_rows(targets, genes, layers=("log_fc", "zscore", "lfcSE", "baseMean")):
    """GWCD4i DESeq2 stats per condition on (targets x genes). Returns dict cond -> dict layer -> array (NaN absent),
    plus per (cond, target) obs metadata (n_cells_target, donor_correlation_all_mean, ontarget_significant)."""
    import h5py
    f = h5py.File(CD4_H5AD, "r")
    obs = f["obs"]

    def cat(k):
        o = obs[k]
        return o["categories"][:].astype(str)[o["codes"][:]]
    tname = cat("target_contrast_gene_name"); cond = cat("culture_condition")
    gname = f["var"]["gene_name"][:].astype(str)
    gcol = {}
    for i, g in enumerate(gname):
        gcol.setdefault(g, i)
    cols = np.array([gcol.get(g, -1) for g in genes]); okc = cols >= 0
    meta_keys = ["n_cells_target", "donor_correlation_all_mean", "ontarget_significant", "ontarget_effect_size",
                 "n_total_de_genes", "guide_correlation_all"]
    out = {}
    for c in sorted(set(cond)):
        rows = {}
        for k, t in enumerate(targets):
            r = np.flatnonzero((tname == t) & (cond == c))
            if len(r):
                rows[k] = int(r[0])
        ks = sorted(rows, key=lambda k: rows[k])
        ridx = [rows[k] for k in ks]
        res = {}
        for L in layers:
            A = np.full((len(targets), len(genes)), np.nan)
            if ridx:
                blk = f["layers"][L][ridx, :]
                A[np.array(ks)[:, None], np.flatnonzero(okc)[None, :]] = blk[:, cols[okc]]
            res[L] = A
        meta = {}
        for mk in meta_keys:
            v = np.full(len(targets), np.nan)
            if ridx:
                v[np.array(ks)] = obs[mk][ridx].astype(float)
            meta[mk] = v
        res["meta"] = meta
        out[c] = res
    f.close()
    return out
