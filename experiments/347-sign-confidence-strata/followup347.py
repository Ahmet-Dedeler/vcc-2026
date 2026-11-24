"""347 follow-up: counterfactuals for the gate-passing strata + mechanism (overshoot)."""
import numpy as np, polars as pl
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
W = REPO / "experiments/323-current-family-validation/work"
OUT = Path(__file__).resolve().parent / "work"

s = np.load(W / "stats-cap15.npz")
r = np.load(W / "reference.npz")
e = np.load(W / "effects.npz")
genes = r["genes"].astype(str)[r["keep"].astype(bool)]
targets = s["targets"].astype(str)
eg = {g: i for i, g in enumerate(e["genes"].astype(str))}
col = np.array([eg.get(g, -1) for g in genes])
miss = col < 0
eff = np.asarray(e["cap15"], dtype=np.float64)[:, np.where(miss, 0, col)]
eff[:, miss] = 0.0
lhat = np.nan_to_num(np.asarray(s["lfc"], dtype=np.float64))
qhat = np.asarray(s["q"], dtype=np.float64)
ref = pl.read_parquet(W / "reference_de.parquet")
gi = {g: i for i, g in enumerate(genes)}

per_t = {}
elig = []
for ti, t in enumerate(targets):
    sub = ref.filter(pl.col("target") == t)
    g = sub["feature"].to_numpy().astype(str)
    q = sub["p_adj"].to_numpy().astype(float)
    l = sub["log2_fold_change"].to_numpy().astype(float)
    gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if gate.sum() < 10:
        continue
    elig.append(ti)
    idx = np.array([gi[x] for x in g[gate]])
    per_t[ti] = dict(l=l[gate], lh=lhat[ti, idx], ef=eff[ti, idx], q=qhat[ti, idx])

# rebuild masks per target
common_full = eff  # (85, G); LOO mean per target computed on the fly for gate genes only
strata = {}
for ti, d in per_t.items():
    ae = np.abs(d["ef"])
    rk = np.argsort(np.argsort(ae, kind="stable")) / max(len(ae) - 1, 1)
    top1 = rk >= 0.9 - 1e-12
    # LOO common on this target's gate genes: need gene positions; recompute via universe idx
    sub = ref.filter(pl.col("target") == targets[ti])
    g = sub["feature"].to_numpy().astype(str)
    q = sub["p_adj"].to_numpy().astype(float)
    l = sub["log2_fold_change"].to_numpy().astype(float)
    gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != targets[ti])
    idx = np.array([gi[x] for x in g[gate]])
    mu = (eff.sum(axis=0) - eff[ti]) / (eff.shape[0] - 1)
    cm = mu[idx]
    thr = 0.000527  # from main script (pooled p10 of |common|)
    agree = (np.sign(d["ef"]) == np.sign(cm)) & (np.abs(cm) > thr)
    head = np.zeros(len(d["l"]), dtype=bool)
    head[np.argsort(-ae, kind="stable")[:400]] = True
    strata[ti] = dict(top1=top1, agree=agree, top1_agree=top1 & agree, head=head,
                      called=d["q"] < .05)


def cf(mask_fn, s_in, s_out):
    vals = []
    for ti, d in per_t.items():
        Dp = np.abs(d["l"]).sum()
        a = np.where(mask_fn(ti), s_in, s_out)
        vals.append(np.abs(a * d["lh"] - d["l"]).sum() / Dp)
    return float(np.mean(vals))


base = cf(lambda ti: np.zeros(len(per_t[ti]["l"]), dtype=bool), 2.0, 1.0)
print("base NMAE:", round(base, 6))
for name, fn in [("a_pt_dec9", lambda ti: strata[ti]["top1"]),
                 ("f_top3_and_agree", lambda ti: strata[ti]["top1_agree"]),
                 ("e_head400", lambda ti: strata[ti]["head"]),
                 ("b_called", lambda ti: strata[ti]["called"])]:
    r2 = cf(fn, 2.0, 1.0)
    r3 = cf(fn, 3.0, 1.0)
    r2h = cf(fn, 2.0, 0.5)
    print("%-16s x2=%0.4f (d=%+.4f) x3=%0.4f x2/o.5=%0.4f" % (name, r2, r2 - base, r3, r2h))

# mechanism: for each stratum, pooled overshoot share + e = dM/dU at x2 (188's identity)
print("\nmechanism at x2 (pooled over gate pairs):")
print("stratum | P(sign ok, w) | overshoot share of stratum |lh| | e=dM/dU | dNMAE_sign(scale-invariant dir)")
for name, fn in [("a_pt_dec9", lambda ti: strata[ti]["top1"]),
                 ("f_top3_and_agree", lambda ti: strata[ti]["top1_agree"]),
                 ("e_head400", lambda ti: strata[ti]["head"]),
                 ("b_called", lambda ti: strata[ti]["called"])]:
    Uh = Ml = Ul = Dl = 0.0
    corr_mass = tot_mass = 0.0
    over = 0.0
    htot = 0.0
    for ti, d in per_t.items():
        m = fn(ti)
        h = np.abs(d["lh"][m])
        rr = np.abs(d["l"][m])
        ok = np.sign(d["lh"][m]) == np.sign(d["l"][m])
        htot += h.sum()
        over += h[~ok].sum() + np.maximum(h[ok] - rr[ok], 0).sum()
        tot_mass += rr.sum()
        corr_mass += rr[ok].sum()
        # e at x2: dU = sum h; dM = sum_ok (min(2h,r)-min(h,r))
        Ul += h.sum()
        Ml += (np.minimum(2 * h[ok], rr[ok]) - np.minimum(h[ok], rr[ok])).sum()
        Dl += rr.sum()
    print("%-16s qw(deliv)=%.3f overshoot=%.3f e=%.3f" % (
        name, corr_mass / tot_mass, over / htot if htot else 0, Ml / Ul if Ul else 0))

# per-target distribution of top-decile design accuracy (is it a few targets?)
print("\nper-target a_pt_dec9: design w-acc and mass share")
accs, masses = [], []
for ti, d in per_t.items():
    m = strata[ti]["top1"]
    Dp = np.abs(d["l"]).sum()
    accs.append((np.abs(d["l"][m]) * (np.sign(d["ef"][m]) == np.sign(d["l"][m]))).sum() / np.abs(d["l"][m]).sum())
    masses.append(np.abs(d["l"][m]).sum() / Dp)
accs = np.array(accs)
print("targets with top-dec accuracy>=0.65: %d/58; median %.3f; q10 %.3f" % (
    int((accs >= 0.65).sum()), float(np.median(accs)), float(np.quantile(accs, 0.1))))
print("mass share median %.3f min %.3f max %.3f" % (float(np.median(masses)), float(np.min(masses)), float(np.max(masses))))
