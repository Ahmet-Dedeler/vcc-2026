"""814 -- build the frozen inputs of predict814.py (the 10-22 reader of PROPOSAL-2's P-A / P-B). Run ONCE, before 10-22.

--table (EC2: reads the three production banks' delta.npy, 2.7 GB each):
    the six raw prop10 features for EVERY gene of the official axis u the production banks' perts u Pan-H1 u VIPerturb
    perts, computed by the same code that built work/run1/features814.csv.gz (feat814.bank_features with its fixed
    2,000-row bank-wide scale; flags814's weak-KD rule; flags814.prop10's H1-control and CD4-guide definitions).
    KNOWN ANSWER (abort): on the 4,825 genes of features814.csv.gz the table equals what flags814.prop10 fed the fit,
    feature by feature, bit-equal (NaN positions included).
    -> work/predict/feat6_allgenes.csv.gz
--constants (local, seconds, < 10 MB): the frozen scoring constants, derived from the frozen 10-16 files only:
    medians used for the fit's fillna (over the 4,355 fit-pool rows), mean / sd (prop10.json), the full-precision
    coefficients and intercept recovered by least squares from propc.tsv's scores (the 10-16 fit's own decision values;
    prop10.json stores them rounded), and the nine decile edges = np.quantile(s over the 4,355 fit-pool rows, 0.1 .. 0.9),
    exactly as flags814.prop10 cut them.
    KNOWN ANSWERS (abort): the fit-pool sha == e3b6b04c... (4,355); the recovered coefficients round to prop10.json's;
    max |s_rebuilt - propc| < 1e-9; every one of the 4,825 decile labels == prop10.tsv; A/B/C in deciles 9-10 == 164.
    -> work/predict/constants814.json (its numbers are then written into predict814.py as literals)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "work/predict"
FEAT = HERE / "work/run1/features814.csv.gz"
FLAGS = HERE / "work/flags"
FEATS = ("prod_logcnt_mean", "weakkd", "in_vip", "pan_cov", "cd4_guide_anysig", "h1_ctl_logcpm")
FIT_POOL_SHA = "e3b6b04c110b6c9a893ccf0ef4e31bd07b210abba99de47ebea1a7af803f183d"


def lines(p):
    return [ln.strip() for ln in Path(p).read_text().splitlines() if ln.strip()]


def sha_list(xs):
    return hashlib.sha256("\n".join(sorted(xs)).encode()).hexdigest()


def six_from(F: pd.DataFrame, weak: pd.Series, h1: pd.Series) -> pd.DataFrame:
    """flags814.prop10's X before its fillna (the same expressions)."""
    return pd.DataFrame(dict(prod_logcnt_mean=F.prod_logcnt_mean, weakkd=weak.astype(float), in_vip=F.in_vip, pan_cov=F.pan_cov,
                             cd4_guide_anysig=F.cd4_guide_anysig.fillna(0.0), h1_ctl_logcpm=np.log10(h1 + 1.0)), index=F.index)


def h1_cpm(index) -> pd.Series:
    g = lines(ROOT / "data/donors/h1_2025/genes.txt")
    ex = np.expm1(np.load(ROOT / "data/donors/h1_2025/b_ctrl.npy").astype(np.float64))
    h1 = pd.Series(1e6 * ex / ex.sum(), index=g)
    return h1[~h1.index.duplicated()].reindex(index)


def table():
    sys.path.insert(0, str(HERE))
    import feat814 as FT
    import flags814 as FL
    genes = [str(g) for g in pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"]]
    univ = set(genes)
    for p in ("data/donors/k562_gwps", "data/donors/orion_hct116", "data/donors/orion_hek293t",
              "experiments/115-gse295214-h1/work/compact-donor", "data/donors/viperturb_k562"):
        univ |= set(lines(ROOT / p / "perts.txt"))
    out = sorted(univ)
    F = pd.DataFrame(index=pd.Index(out, name="gene"))
    for name in ("k562", "hct", "hek"):
        F = F.join(FT.bank_features(name, FT.BANKS[name], out), how="left")
        F[f"{name}_cov"] = F[f"{name}_cov"].fillna(0.0)
        F[f"{name}_logcnt"] = np.log1p(F[f"{name}_cnt"])
    F["prod_logcnt_mean"] = F[["k562_logcnt", "hct_logcnt", "hek_logcnt"]].mean(axis=1)
    F["pan_cov"] = [float(t in set(lines(FT.BANKS["pan"] / "perts.txt"))) for t in out]
    F["in_vip"] = [float(t in set(lines(FT.LISTS["vip"]))) for t in out]
    kd = pd.read_csv(FT.CD4_KD)
    kd = kd[kd["culture_condition"] == "Rest"].copy()
    kd["gene"] = kd.iloc[:, 0].astype(str).str.rsplit("-", n=1).str[0]
    F["cd4_guide_anysig"] = kd.groupby("gene")["signif_knockdown"].any().astype(float).reindex(out).values
    F["in_essscreen"] = [float(t in set(lines(FT.LISTS["rpe1"])) | set(lines(FT.LISTS["jurkat"]))) for t in out]
    weak = FL.flags(F)["weakkd"]
    X = six_from(F, weak, h1_cpm(F.index))
    # KA: on features814's genes, the same six inputs as the 10-16 fit (before fillna), bit-equal
    R = pd.read_csv(FEAT, index_col=0, float_precision="round_trip")
    Xr = six_from(R, FL.flags(R)["weakkd"], h1_cpm(R.index))
    a, b = X.reindex(R.index).to_numpy(), Xr.to_numpy()
    same = np.array_equal(np.nan_to_num(a, nan=-9e99), np.nan_to_num(b, nan=-9e99))
    assert same, f"KA FAIL: the all-gene table differs from features814 on its genes ({np.nanmax(np.abs(a - b))})"
    print(f"KA PASS: six features bit-equal to features814 on its {len(R)} genes; table {len(X)} genes", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    X.to_csv(OUT / "feat6_allgenes.csv.gz")


def constants():
    sys.path.insert(0, str(HERE))
    import flags814 as FL
    R = pd.read_csv(FEAT, index_col=0, float_precision="round_trip")
    P = (R.in_null_pool > 0) & (R.in_essscreen == 0)
    assert int(P.sum()) == 4355 and sha_list(R.index[P]) == FIT_POOL_SHA, "KA FAIL: fit pool"
    X = six_from(R, FL.flags(R)["weakkd"], h1_cpm(R.index))
    med = X[P].median()
    Xf = X.fillna(med)
    pj = json.loads((FLAGS / "prop10.json").read_text())
    mu, sd = pd.Series(pj["mean"]), pd.Series(pj["sd"])
    Z = ((Xf - mu[list(FEATS)]) / sd[list(FEATS)]).to_numpy()
    s = pd.Series({g: float(v) for g, v in (ln.split("\t") for ln in (FLAGS / "propc.tsv").read_text().splitlines())}).reindex(R.index)
    A1 = np.hstack([Z, np.ones((Z.shape[0], 1))])
    sol, *_ = np.linalg.lstsq(A1, s.to_numpy(), rcond=None)
    coef, b0 = sol[:6], sol[6]
    assert np.allclose(np.round(coef, 4), [pj["coef"][f] for f in FEATS], atol=0) and abs(b0 - pj["intercept"]) < 1e-9, \
        f"KA FAIL: recovered coefficients {coef} {b0}"
    s2 = A1 @ sol
    err = float(np.max(np.abs(s2 - s.to_numpy())))
    assert err < 1e-9, f"KA FAIL: rebuilt scores differ from propc by {err}"
    edges = np.quantile(s[P], np.arange(1, 10) / 10)
    # labels from the frozen scores themselves (np.quantile at q 0.5 over 4,355 rows lands exactly on one pool score, so a
    # last-bit difference in a rebuilt score can flip that tie: predict814 uses propc's frozen scores wherever they exist)
    lab = pd.Series([f"p{v}" for v in np.digitize(s.to_numpy(), edges)], index=R.index)
    ref = pd.read_csv(FLAGS / "prop10.tsv", sep="\t", header=None, index_col=0)[1].reindex(R.index)
    assert (lab == ref).all(), f"KA FAIL: {int((lab != ref).sum())} decile labels differ from prop10.tsv"
    lab2 = pd.Series([f"p{v}" for v in np.digitize(s2, edges)], index=R.index)
    flips = sorted(lab2.index[lab2 != ref])
    on_edge = [g for g in flips if float(s[g]) in set(edges.tolist())]
    assert flips == on_edge, f"KA FAIL: rebuilt-score label flips not explained by an exact edge tie: {flips}"
    abc = int(lab[R.in_abc > 0].isin(["p8", "p9"]).sum())
    assert abc == 164, f"KA FAIL: A/B/C in deciles 9-10 = {abc}"
    rec = dict(features=list(FEATS), fill_median=med.to_dict(), mean=mu[list(FEATS)].to_dict(), sd=sd[list(FEATS)].to_dict(),
               coef=dict(zip(FEATS, coef.tolist())), intercept=float(b0), edges=edges.tolist(), max_abs_err_vs_propc=err,
               abc_deciles_9_10=abc, fit_pool_sha=FIT_POOL_SHA,
               rebuilt_score_flips_at_exact_edge=flips)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "constants814.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))
    print("KA PASS: fit pool, coefficients, scores, 4,825 decile labels, A/B/C 164")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", action="store_true")
    ap.add_argument("--constants", action="store_true")
    a = ap.parse_args()
    if a.constants:
        constants()
    if a.table:
        table()
