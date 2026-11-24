"""426d: reach lives in the weak tercile (406: n_real ~6, reach 0.67). What are the weak targets' reference DE genes -- cis
neighbours (TSS <= 5 kb), and where do we rank them / with which sign? H1, 426 base stats (delivered q/p/lfc)."""
from pathlib import Path
import numpy as np, pandas as pd, polars as pl
ROOT = Path(__file__).resolve().parents[2]
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; gi = {g: i for i, g in enumerate(genes)}
s = np.load(Path(__file__).parent / "work/h1/stats-base_s0_bps.npz"); targets = s["targets"].astype(str)
L = s["lfc"]; Q = s["q"]; P = s["p"]
tss = pd.read_csv(ROOT / "experiments/425-crispri-cis-neighbours/work/tss_v44.tsv", sep="\t").set_index("gene_name")
ref = pl.read_parquet(W323 / "reference_de.parquet")
rows = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    nreal = int(ok.sum())
    if nreal == 0 or nreal > 13:
        continue
    # our order: called first, then ascending q, p, descending |lfc|
    lk = np.nan_to_num(L[k]); order_key = np.lexsort((-np.abs(lk), P[k], Q[k], ~(Q[k] < .05)))
    rank = np.empty(order_key.size, int); rank[order_key] = np.arange(order_key.size)
    tt = tss.loc[t] if t in tss.index else None
    for gg, ll in zip(g[ok], l[ok]):
        j = gi[gg]
        d = np.inf
        if tt is not None and gg in tss.index and tss.loc[gg, "chrom"] == tt["chrom"]:
            d = abs(int(tss.loc[gg, "tss"]) - int(tt["tss"]))
        rows.append(dict(target=t, n_real=nreal, gene=gg, ref_lfc=ll, cis=d <= 5000, dist=d, our_rank=int(rank[j]),
                         our_called=bool(Q[k, j] < .05), sign_ok=bool(np.sign(lk[j]) == np.sign(ll))))
df = pd.DataFrame(rows)
print(f"weak targets (1 <= n_real <= 13): {df.target.nunique()}, DE genes {len(df)}")
print(f"cis (<=5 kb) share of weak-target DE genes: {df.cis.mean():.3f}; targets with a cis DE gene: {df[df.cis].target.nunique()}")
print(f"sign ok: cis {df[df.cis].sign_ok.mean():.3f} (n {df.cis.sum()}), non-cis {df[~df.cis].sign_ok.mean():.3f}")
print(f"our rank of the DE genes (0 = first): cis median {df[df.cis].our_rank.median():.0f}, non-cis median {df[~df.cis].our_rank.median():.0f}; called share cis {df[df.cis].our_called.mean():.2f} non-cis {df[~df.cis].our_called.mean():.2f}")
print(f"ref lfc: cis median {df[df.cis].ref_lfc.median():+.2f}, non-cis median |lfc| {df[~df.cis].ref_lfc.abs().median():.2f}, non-cis down share {(df[~df.cis].ref_lfc < 0).mean():.2f}")
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 200)
print(df.sort_values(["n_real", "target", "our_rank"]).to_string(index=False, float_format=lambda v: f"{v:.2f}"))
