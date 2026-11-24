"""814 -- target flag lists for gains814's restricted panelrank (light: reads features814.csv.gz, ~1 MB).
Writes work/flags/<name>.txt (one target per line, all flagged genes of the feature table) and prints each flag's share in
A/B/C's panel and in L1-L6 (the matched keep share = ABC share / L share).
  essscreen  in the RPE1 or Jurkat essential-gene screen (data/donors/{rpe1,nadig_jurkat}/perts.txt)
  weakkd     estimated knockdown < 0.5 in K562-GWPS or Orion HCT116 where the bank measures the target's own gene
             (kd = 1 - expm1(log1p(m) + own_lfc) / m, m = the bank's control mean of the target at TS 50,000)
  silent     zero robust DE genes (bank-wide 518 statistic) in at least one production bank that carries the target
  esskd      essscreen or weakkd
    uv run python experiments/814-panel-rule/flags814.py            # flag lists
    uv run python experiments/814-panel-rule/flags814.py --classes  # class files for gains814 --match
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FEAT = HERE / "work/run1/features814.csv.gz"
BANKS = {"k562": ROOT / "data/donors/k562_gwps", "hct": ROOT / "data/donors/orion_hct116"}


def flags(F: pd.DataFrame) -> dict[str, pd.Series]:
    kd = {}
    for b, p in BANKS.items():
        S = np.expm1(np.load(p / "b_ctrl.npy").astype(np.float64)).sum()
        m = F[f"{b}_basal_cpm"] * S / 1e6
        kd[b] = 1 - np.expm1(np.log1p(m) + F[f"{b}_own_lfc"]) / m
    weak = (kd["k562"] < 0.5).fillna(False) | (kd["hct"] < 0.5).fillna(False)
    silent = (F.k562_cnt == 0) | (F.hct_cnt == 0) | (F.hek_cnt == 0)
    ess = F.in_essscreen > 0
    return dict(essscreen=ess, weakkd=weak, silent=silent, esskd=ess | weak)


def main():
    F = pd.read_csv(FEAT, index_col=0)
    A = F.in_abc > 0
    L = F[[f"in_L{i}" for i in range(1, 7)]].sum(axis=1) > 0
    out = HERE / "work/flags"
    out.mkdir(parents=True, exist_ok=True)
    for name, s in flags(F).items():
        (out / f"{name}.txt").write_text("\n".join(sorted(F.index[s])) + "\n")
        a, l = float(s[A].mean()), float(s[L].mean())
        print(f"{name:10s} ABC {a:.3f}  L {l:.3f}  matched keep share {a / l if l else float('nan'):.3f}  ({int(s.sum())} genes)")


if __name__ == "__main__" and "--classes" not in sys.argv:
    main()


def write_classes():
    """class files for gains814 --match: resp5 = pool quintile of prod_logcnt_mean (production-bank responsiveness, edges
    from the null pool minus the essential screens); kd2 = weakkd yes/no; resp5kd = their product."""
    F = pd.read_csv(FEAT, index_col=0)
    P = (F.in_null_pool > 0) & (F.in_essscreen == 0)
    x = F.prod_logcnt_mean
    e = np.nanquantile(x[P], [0.2, 0.4, 0.6, 0.8])
    q = np.digitize(x, e)
    wk = flags(F)["weakkd"]
    out = HERE / "work/flags"
    for name, lab in (("resp5", [f"r{v}" for v in q]), ("kd2", ["kd0" if w else "kd1" for w in wk]),
                      ("resp5kd", [f"r{v}_{'w' if w else 'k'}" for v, w in zip(q, wk)])):
        (out / f"{name}.tsv").write_text("".join(f"{g}\t{c}\n" for g, c in zip(F.index, lab)))
    print("resp5 edges", np.round(e, 4).tolist())
    prop10(F, P, wk)


PROP_FEATURES = ("prod_logcnt_mean", "weakkd", "in_vip", "pan_cov", "cd4_guide_anysig", "h1_ctl_logcpm")


def prop10(F: pd.DataFrame, P: pd.Series, wk: pd.Series):
    """prop10 = pool decile of a propensity score P(A/B/C panel | features) vs the arc800 null pool minus the essential
    screens. Six features, frozen before the fit (RESULTS.md section 4): production-bank responsiveness, weak knockdown,
    VIPerturb-K562 coverage, Pan-2025 H1 coverage, CD4 guide-level knockdown, H1 control expression (VCC-2025 H1 controls,
    data/donors/h1_2025/b_ctrl.npy, CPM over its gene axis). L2 logistic regression (C 1, standardised); 5-fold
    out-of-fold AUC reported. A/B/C targets sit outside the pool (it excludes them), so the fit is panel vs pool."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import roc_auc_score
    g = [ln.strip() for ln in (ROOT / "data/donors/h1_2025/genes.txt").read_text().splitlines() if ln.strip()]
    ex = np.expm1(np.load(ROOT / "data/donors/h1_2025/b_ctrl.npy").astype(np.float64))
    h1 = pd.Series(1e6 * ex / ex.sum(), index=g)
    h1 = h1[~h1.index.duplicated()].reindex(F.index)
    X = pd.DataFrame(dict(prod_logcnt_mean=F.prod_logcnt_mean, weakkd=wk.astype(float), in_vip=F.in_vip, pan_cov=F.pan_cov,
                          cd4_guide_anysig=F.cd4_guide_anysig.fillna(0.0), h1_ctl_logcpm=np.log10(h1 + 1.0)), index=F.index)
    X = X.fillna(X[P].median())
    A = F.in_abc > 0
    fit = A | P
    mu, sd = X[fit].mean(), X[fit].std()
    Z = (X - mu) / sd
    y = A[fit].astype(int).to_numpy()
    lr = LogisticRegression(C=1.0, max_iter=1000)
    oof = cross_val_predict(lr, Z[fit].to_numpy(), y, cv=StratifiedKFold(5, shuffle=True, random_state=814), method="predict_proba")[:, 1]
    lr.fit(Z[fit].to_numpy(), y)
    s = pd.Series(lr.decision_function(Z.to_numpy()), index=F.index)
    e = np.quantile(s[P], np.arange(1, 10) / 10)
    d = np.digitize(s, e)
    (HERE / "work/flags/prop10.tsv").write_text("".join(f"{gg}\tp{v}\n" for gg, v in zip(F.index, d)))
    (HERE / "work/flags/propc.tsv").write_text("".join(f"{gg}\t{float(v)!r}\n" for gg, v in zip(F.index, s.to_numpy())))
    rec = dict(features=list(PROP_FEATURES), coef=dict(zip(PROP_FEATURES, lr.coef_[0].round(4).tolist())),
               intercept=float(lr.intercept_[0]), oof_auc=float(roc_auc_score(y, oof)), mean=mu.to_dict(), sd=sd.to_dict(),
               abc_decile_share=np.bincount(d[A.to_numpy()], minlength=10).tolist())
    (HERE / "work/flags/prop10.json").write_text(json.dumps(rec, indent=1))
    print("prop10:", json.dumps({k: rec[k] for k in ("coef", "oof_auc", "abc_decile_share")}))


if __name__ == "__main__" and "--classes" in sys.argv:
    write_classes()
