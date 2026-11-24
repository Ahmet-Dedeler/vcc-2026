"""284 — POST-HOC (PREREG §8) and the deliverable.

(a) the remainder test repeated with STRONGER predictors subtracted — features + measured
    reliability `r_bar`, and + response magnitude — i.e. is the reproducible remainder still there
    once the confounder 281 named is removed? (Both are measured on ALL 28 sources, so subtracting
    them removes information shared by the two halves: the test is conservative.)
(b) the remainder outside the giant component, and the attenuation-corrected ceiling.
(c) the descriptive tail: 20 genes each way whose remainder agrees in both halves.
(d) work/pred284.csv.

  PYTHONDONTWRITEBYTECODE=1 uv run python extra284.py

Outputs: work/extra284.json, work/extra284.log, work/pred284.csv
"""
from __future__ import annotations

import hashlib
import json
import sys
import time

import numpy as np
import pandas as pd

import lib284 as L
from rem284 import oof_ridge, spearman_brown

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def main() -> None:
    h = hashlib.sha256((L.HERE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    if h != (L.WORK / "PREREGISTRATION.sha256").read_text().split()[0]:
        sys.exit("preregistration digest mismatch — aborting")
    P(f"prereg {h} MATCH")

    held = L.load_held()
    mem = L.memberships(held)
    sp = np.load(L.WORK / "splits284.npz", allow_pickle=True)
    comp = sp["comp"]
    fold = L.component_folds(comp, 5)
    rem = np.load(L.WORK / "rem284.npz", allow_pickle=True)
    C1, C2, r1, r2, ok = rem["C1"], rem["C2"], rem["r1"], rem["r2"], rem["ok"]
    n1, n2 = rem["n1"], rem["n2"]
    z = np.load(L.HELD279, allow_pickle=True)
    mag = np.asarray(z["mag"], dtype=float)
    rbar = held["r_bar"].to_numpy(float)
    X = held[L.FEATS].to_numpy(float)
    out: dict = {"prereg_sha256": h}

    # ---------- (a) stronger predictors -------------------------------------------------------
    P("\n(a) POST-HOC: the same remainder test with stronger predictors subtracted")
    sets = {
        "features (frozen, PREREG §7)": X,
        "features + r_bar": np.column_stack([X, rbar]),
        "features + r_bar + magnitude": np.column_stack([X, rbar, mag]),
        "r_bar alone": rbar.reshape(-1, 1),
        "r_bar + magnitude": np.column_stack([rbar, mag]),
    }
    rows = {}
    for nm, M in sets.items():
        p1, p2 = oof_ridge(M, C1.astype(float), fold), oof_ridge(M, C2.astype(float), fold)
        q1, q2 = C1 - p1, C2 - p2
        m = ok & np.isfinite(q1) & np.isfinite(q2)
        b = L.boot_rho(q1[m], q2[m], comp[m], np.random.default_rng(L.SEED + 701))
        rows[nm] = dict(**b,
                        oof_rho_half1=L.spearman(C1[np.isfinite(p1)], p1[np.isfinite(p1)]),
                        oof_rho_half2=L.spearman(C2[np.isfinite(p2)], p2[np.isfinite(p2)]))
        P(f"  {nm:32s} remainder rho {b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}]  "
          f"(the predictor's own OOF rho {rows[nm]['oof_rho_half1']:+.3f} / "
          f"{rows[nm]['oof_rho_half2']:+.3f}; n {int(m.sum())})")
    out["stronger_predictors"] = rows

    # ---------- (b) outside the giant component, and the ceiling's own attenuation -------------
    big = int(np.argmax(np.bincount(comp)))
    m = ok & (comp != big)
    b = L.boot_rho(r1[m], r2[m], comp[m], np.random.default_rng(L.SEED + 703))
    bc = L.boot_rho(C1[m], C2[m], comp[m], np.random.default_rng(L.SEED + 704))
    P(f"\n(b) outside the giant component (n {int(m.sum())} of {int(ok.sum())}): remainder "
      f"{b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}], ceiling {bc['point']:+.4f} "
      f"[{bc['lo']:+.4f}, {bc['hi']:+.4f}]")
    m2 = ok & (comp == big)
    b2 = L.boot_rho(r1[m2], r2[m2], comp[m2], np.random.default_rng(L.SEED + 705), n=200)
    P(f"    inside the giant component (n {int(m2.sum())}): remainder {b2['point']:+.4f} "
      f"(one cluster — no interval is defined by the component bootstrap)")
    out["outside_giant"] = dict(n=int(m.sum()), remainder=b, ceiling=bc,
                                inside_giant=dict(n=int(m2.sum()), remainder=b2["point"]))

    # ceiling attenuation, same construction as the remainder's
    j = json.loads((L.WORK / "rem284.json").read_text())
    zz = np.load(L.HELD279, allow_pickle=True)
    Ccross, pairs = np.asarray(zz["Ccross"]), np.asarray(zz["pairs_cross"])
    keys, lines = np.asarray(zz["keys"]).astype(str), np.asarray(zz["lines"]).astype(str)
    h1l = set(j["primary"]["half1_lines"])
    s1 = np.flatnonzero(np.isin(lines, list(h1l)))
    s2 = np.flatnonzero(~np.isin(lines, list(h1l)))
    p1i = np.array([q for q, (a, bb) in enumerate(pairs) if a in set(s1) and bb in set(s1)])
    p2i = np.array([q for q, (a, bb) in enumerate(pairs) if a in set(s2) and bb in set(s2)])
    rr = np.random.default_rng(L.SEED + 401)          # the same stream as rem284's attenuation
    att = {}
    for nm, pp in (("half1", p1i), ("half2", p2i)):
        q = rr.permutation(len(pp))
        a, bb = pp[q[: len(pp) // 2]], pp[q[len(pp) // 2:]]
        Ca, Cb = L.cmap(Ccross, a).astype(float), L.cmap(Ccross, bb).astype(float)
        mm = np.isfinite(Ca) & np.isfinite(Cb) & ok
        rho = L.spearman(Ca[mm], Cb[mm])
        att[nm] = dict(rho_subhalf=rho, rho_full=spearman_brown(rho), n=int(mm.sum()))
    den = float(np.sqrt(att["half1"]["rho_full"] * att["half2"]["rho_full"]))
    ceil_pt = j["primary"]["boot_ceiling"]["point"]
    out["ceiling_attenuation"] = dict(**att, denominator=den, corrected=float(ceil_pt / den))
    P(f"    ceiling attenuation: half1 {att['half1']['rho_full']:+.4f}, half2 "
      f"{att['half2']['rho_full']:+.4f} -> corrected ceiling {ceil_pt / den:+.4f}")

    # ---------- (c) the descriptive tail ------------------------------------------------------
    P("\n(c) POST-HOC descriptive: genes whose remainder agrees in BOTH halves")
    gene = held["gene"].to_numpy().astype(str)
    agree_pos = ok & (r1 > 0) & (r2 > 0)
    agree_neg = ok & (r1 < 0) & (r2 < 0)
    score_pos = np.where(agree_pos, np.minimum(r1, r2), -np.inf)
    score_neg = np.where(agree_neg, np.maximum(r1, r2), np.inf)
    top_pos = np.argsort(-score_pos)[:20]
    top_neg = np.argsort(score_neg)[:20]

    def describe(idx):
        return [dict(gene=gene[i], r1=float(r1[i]), r2=float(r2[i]), C1=float(C1[i]),
                     C2=float(C2[i]), n_pairs_1=int(n1[i]), n_pairs_2=int(n2[i]),
                     families="|".join(mem[i]) if mem[i] else "",
                     humap=float(held["humap_n_complexes"].to_numpy()[i]),
                     reactome=float(held["reactome_n_complexes"].to_numpy()[i]),
                     dm_mean_effect=float(held["dm_mean_effect"].to_numpy()[i]),
                     replogle_essential=float(held["replogle_essential"].to_numpy()[i]),
                     r_bar=float(rbar[i])) for i in idx]

    out["tail_above"] = describe(top_pos)
    out["tail_below"] = describe(top_neg)
    P("  20 genes ABOVE prediction in both halves (gene, r1, r2, pairs, families):")
    for d in out["tail_above"]:
        P(f"    {d['gene']:12s} {d['r1']:+.3f} / {d['r2']:+.3f}  pairs {d['n_pairs_1']:2d}/"
          f"{d['n_pairs_2']:2d}  ess {d['replogle_essential']:.0f}  {d['families'][:48]}")
    P("  20 genes BELOW prediction in both halves:")
    for d in out["tail_below"]:
        P(f"    {d['gene']:12s} {d['r1']:+.3f} / {d['r2']:+.3f}  pairs {d['n_pairs_1']:2d}/"
          f"{d['n_pairs_2']:2d}  ess {d['replogle_essential']:.0f}  {d['families'][:48]}")
    for nm, idx in (("above", top_pos), ("below", top_neg)):
        ess = float(np.nanmean(held["replogle_essential"].to_numpy()[idx]))
        cx = float(np.nanmean(held["humap_n_complexes"].to_numpy()[idx]))
        P(f"    {nm}: mean replogle_essential {ess:.2f}, mean hu.MAP complexes {cx:.2f}, "
          f"{sum(1 for i in idx if mem[i])}/20 with an HGNC family")

    # ---------- (d) the deliverable -----------------------------------------------------------
    pr = np.load(L.WORK / "preds284.npz", allow_pickle=True)
    oo = np.load(L.WORK / "oof284.npz", allow_pickle=True)
    teS1, teS2 = sp["teS1"], sp["teS2"]
    d = pd.DataFrame(dict(
        gene=gene, observed=held["y_raw"].to_numpy(float),
        observed_rel=held["y_rel"].to_numpy(float),
        observed_core=held["c3_balanced_core"].to_numpy(float),
        r_bar=rbar, n_pairs=held["n_pairs"].to_numpy(float),
        gene_group_all=held["gene_group_all"].fillna("").to_numpy(),
        component=comp, component_fold=fold,
        split_s1=np.where(teS1, "held", "train"), split_s2=np.where(teS2, "held", "train"),
        split_281=np.where(sp["te281"], "held", "train"),
        split_random=np.where(sp["teRND"], "held", "train")))
    for tag, te in (("s1", teS1), ("s2", teS2)):
        for arm in ("features_ridge", "features_gbm", "reliability_only",
                    "features_plus_reliability"):
            col = np.full(len(d), np.nan)
            key = "S1_all_membership" if tag == "s1" else "S2_component"
            col[pr[f"{key}::idx"]] = pr[f"{key}::{arm}"]
            d[f"pred_{tag}_{arm}"] = col
    d["pred_oof_raw"] = oo["pred_raw"]
    d["pred_oof_core"] = oo["pred_core"]
    d["C_half1"] = C1
    d["C_half2"] = C2
    d["n_pairs_half1"] = n1
    d["n_pairs_half2"] = n2
    d["remainder_half1"] = np.where(ok, r1, np.nan)
    d["remainder_half2"] = np.where(ok, r2, np.nan)
    d.to_csv(L.WORK / "pred284.csv", index=False)
    P(f"\nwrote work/pred284.csv ({len(d)} rows x {d.shape[1]} columns)")

    out["log"] = LOG
    (L.WORK / "extra284.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "extra284.log").write_text("\n".join(LOG))
    P("wrote work/extra284.json")


if __name__ == "__main__":
    main()
