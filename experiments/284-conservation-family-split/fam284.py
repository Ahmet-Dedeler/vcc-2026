"""284 — the family axis under S1 (PREREG §6): grouped out-of-fold predictions in which a test gene
shares NO family with its training rows (folds = whole components), observed vs predicted family
medians on the balanced core, and the zinc-finger deficit.

  PYTHONDONTWRITEBYTECODE=1 uv run python fam284.py

Outputs: work/fam284.json, work/fam284.log, work/oof284.npz
"""
from __future__ import annotations

import hashlib
import json
import sys
import time

import numpy as np
import pandas as pd

import lib284 as L

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def oof(X, y, fold, gbm: bool = False):
    p = np.full(len(y), np.nan)
    ok = np.isfinite(y)
    for k in range(int(fold.max()) + 1):
        te, tr = ok & (fold == k), ok & (fold != k)
        if te.sum() == 0 or tr.sum() < 50:
            continue
        if gbm:
            p[te], _, _ = L.fit_gbm(X[tr], y[tr], X[te])
        else:
            p[te], _ = L.fit_ridge(X[tr], y[tr], X[te])
    return p


def family_table(mem, y, pred, mask, min_n: int, first_only: np.ndarray | None = None):
    """observed vs predicted family medians over families with >= min_n genes."""
    acc: dict[str, list[int]] = {}
    for i in np.flatnonzero(mask):
        fams = [first_only[i]] if first_only is not None else mem[i]
        for f in fams:
            if f:
                acc.setdefault(f, []).append(i)
    rows = []
    for f, idx in acc.items():
        if len(idx) >= min_n:
            rows.append(dict(family=f, n=len(idx), observed=float(np.median(y[idx])),
                             predicted=float(np.median(pred[idx]))))
    return pd.DataFrame(rows).sort_values("observed", ascending=False).reset_index(drop=True)


def main() -> None:
    h = hashlib.sha256((L.HERE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    if h != (L.WORK / "PREREGISTRATION.sha256").read_text().split()[0]:
        sys.exit("preregistration digest mismatch — aborting")
    P(f"prereg {h} MATCH")

    held = L.load_held()
    mem = L.memberships(held)
    sp = np.load(L.WORK / "splits284.npz", allow_pickle=True)
    comp = sp["comp"]
    block = np.asarray(sp["block"]).astype(str)
    fold_c = L.component_folds(comp, 5)
    X = held[L.FEATS].to_numpy(float)
    y_raw = held["y_raw"].to_numpy(float)
    y_core = held["c3_balanced_core"].to_numpy(float)
    isz = held["znf_class"].to_numpy() == "ZNF"
    isn = held["znf_class"].to_numpy() == "NONZNF"
    P(f"component folds (STRICT: a test gene shares no family with its training rows): "
      f"sizes {np.bincount(fold_c).tolist()}")

    # 281's own grouped folds (whole first-family blocks) for comparison
    ub = np.array(sorted(set(block.tolist())))
    bidx = {b: i for i, b in enumerate(ub)}
    bcode = np.array([bidx[b] for b in block])
    fold_b = L.component_folds(bcode, 5)
    P(f"281-style first-family folds: sizes {np.bincount(fold_b).tolist()}")

    out: dict = {"prereg_sha256": h}
    res = {}
    for tag, fold in (("component_folds", fold_c), ("first_family_folds", fold_b)):
        for tname, y in (("core", y_core), ("raw", y_raw)):
            p = oof(X, y, fold)
            ok = np.isfinite(y) & np.isfinite(p)
            rho = L.spearman(y[ok], p[ok])
            res[f"{tag}::{tname}"] = dict(n=int(ok.sum()), oof_rho=rho, oof_r2=L.r2(y[ok], p[ok]),
                                          oof_auc10=L.auc_bottom_decile(y[ok], p[ok]))
            P(f"  OOF {tag:20s} {tname:5s}: n {int(ok.sum())}, rho {rho:+.4f}, "
              f"R2 {L.r2(y[ok], p[ok]):+.4f}, AUC10 {L.auc_bottom_decile(y[ok], p[ok]):.4f}")
            if tag == "component_folds" and tname == "core":
                pred_core = p
            if tag == "component_folds" and tname == "raw":
                pred_raw = p
            if tag == "first_family_folds" and tname == "core":
                pred_core_ff = p
    out["oof"] = res

    # 281's own out-of-fold column, read off its deliverable (the reference for §6)
    p281 = pd.read_csv(L.PRED281, usecols=["gene", "set", "pred_oof_core_ridge"])
    p281 = p281[p281["set"] == "held"].reset_index(drop=True)
    assert (p281["gene"].to_numpy() == held["gene"].to_numpy()).all()
    pred_core_281 = p281["pred_oof_core_ridge"].to_numpy(float)
    ok281 = np.isfinite(y_core) & np.isfinite(pred_core_281)
    P(f"  281's own OOF core column: n {int(ok281.sum())}, rho "
      f"{L.spearman(y_core[ok281], pred_core_281[ok281]):+.4f} (R281 §5 reports +0.3398)")
    out["oof"]["pred281_column::core"] = dict(
        n=int(ok281.sum()), oof_rho=L.spearman(y_core[ok281], pred_core_281[ok281]))

    # ---- the family axis on the balanced core, component folds ------------------------------
    okc = np.isfinite(y_core) & np.isfinite(pred_core)
    P(f"\nfamily axis on the balanced core, {int(okc.sum())} genes")
    for label, first in (("ALL memberships (primary)", None),
                         ("first HGNC group (281's version)", held["gene_group1"].fillna("")
                          .to_numpy().astype(str))):
        tb = family_table(mem, y_core, pred_core, okc, 20, first)
        rho = L.spearman(tb["observed"].to_numpy(), tb["predicted"].to_numpy())
        rb = np.random.default_rng(L.SEED + 601)
        bs = []
        for _ in range(2000):
            s = rb.integers(0, len(tb), len(tb))
            bs.append(L.spearman(tb["observed"].to_numpy()[s], tb["predicted"].to_numpy()[s]))
        lo, hi = float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))
        P(f"  {label}: {len(tb)} families with >= 20 genes -> Spearman(observed median, "
          f"predicted median) {rho:+.4f} [{lo:+.4f}, {hi:+.4f}]")
        key = "family_axis_all" if first is None else "family_axis_first"
        out[key] = dict(n_families=int(len(tb)), rho=rho, lo=lo, hi=hi,
                        table=tb.to_dict(orient="records"))
        if first is None:
            P("    top 8 / bottom 8 by observed median:")
            for _, r in pd.concat([tb.head(8), tb.tail(8)]).iterrows():
                P(f"      {r['family'][:44]:46s} n {int(r['n']):4d}  obs {r['observed']:+.4f}  "
                  f"pred {r['predicted']:+.4f}")

    # the same family axis through 281's own OOF column, and through the 281-style folds
    for lab, pp in (("281's own OOF column", pred_core_281),
                    ("first-family folds, my code", pred_core_ff)):
        okx = np.isfinite(y_core) & np.isfinite(pp)
        tb = family_table(mem, y_core, pp, okx, 20,
                          held["gene_group1"].fillna("").to_numpy().astype(str))
        P(f"  [comparison] {lab}: {len(tb)} first-group families -> "
          f"{L.spearman(tb['observed'].to_numpy(), tb['predicted'].to_numpy()):+.4f}")
        out[f"family_axis_cmp::{lab}"] = dict(
            n_families=int(len(tb)),
            rho=L.spearman(tb["observed"].to_numpy(), tb["predicted"].to_numpy()))

    # ---- the zinc-finger deficit, predicted vs observed --------------------------------------
    for tname, y, p in (("core", y_core, pred_core), ("raw", y_raw, pred_raw),
                        ("core_first_family_folds", y_core, pred_core_ff),
                        ("core_pred281_column", y_core, pred_core_281)):
        ok = np.isfinite(y) & np.isfinite(p)
        oz, on_ = float(np.median(y[isz & ok])), float(np.median(y[isn & ok]))
        pz, pn = float(np.median(p[isz & ok])), float(np.median(p[isn & ok]))
        out[f"znf_{tname}"] = dict(n_znf=int((isz & ok).sum()), n_nonznf=int((isn & ok).sum()),
                                   observed_d=oz - on_, predicted_d=pz - pn,
                                   frac_predicted=float((pz - pn) / (oz - on_))
                                   if (oz - on_) != 0 else float("nan"))
        P(f"  ZNF deficit ({tname}): observed Δmedian {oz - on_:+.5f}, predicted {pz - pn:+.5f} "
          f"({100 * (pz - pn) / (oz - on_):.0f} % of it); n {int((isz & ok).sum())} ZNF")

    np.savez_compressed(L.WORK / "oof284.npz", pred_core=pred_core, pred_raw=pred_raw,
                        fold_component=fold_c, fold_first_family=fold_b,
                        gene=held["gene"].to_numpy().astype(str))
    out["log"] = LOG
    (L.WORK / "fam284.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "fam284.log").write_text("\n".join(LOG))
    P("wrote work/fam284.json, work/oof284.npz")


if __name__ == "__main__":
    main()
