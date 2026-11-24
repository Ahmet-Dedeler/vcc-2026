"""296 -- the frozen stratification, P1-P5, the failing-construct phenomenon, the failure rates and
`efficacy296.csv`.  Reads only the per-library outputs written by k562_296 / nadig_296 / broad_296 / cd4_296.

    PYTHONDONTWRITEBYTECODE=1 uv run python strat296.py
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib296 as X  # noqa: E402
from lib296 import L  # noqa: E402

P = L.Logger(X.WORK / "strat296.log")
W = X.WORK
CRISPRi = ("K562-GWPS", "Nadig-jurkat", "Nadig-hepg2", "CD4-GWCD4i")
CUT_ABS = 1.0
CUT_FRAC = 0.5


# ---------------------------------------------------------------------------- loaders
def load_k562_like(tag: str, name: str, n: int):
    j = json.load(open(W / f"{tag}.json"))
    z = np.load(W / f"{tag}.npz", allow_pickle=True)
    if str(n) not in j["rungs"]:
        return None
    kdc = j["kd_cons"]
    tgts = z[f"{n}|targets"].astype(str)
    big = z[f"{n}|c_big"].astype(str); small = z[f"{n}|c_small"].astype(str)
    d = dict(library=name, unit=f"n={n} cells/side", targets=tgts,
             a=z[f"{n}|a"], a2=z[f"{n}|a2"], b=z[f"{n}|b"], x=z[f"{n}|x"],
             pds_a=z[f"{n}|pds_a"], pds_b=z[f"{n}|pds_b"], pds_x=z[f"{n}|pds_x"],
             c_big=big, c_small=small,
             kd_big=np.asarray([kdc[c]["lfc"] for c in big], float),
             kd_small=np.asarray([kdc[c]["lfc"] for c in small], float),
             kf_big=np.asarray([kdc[c]["kfrac"] for c in big], float),
             kf_small=np.asarray([kdc[c]["kfrac"] for c in small], float),
             cells_big=np.asarray([kdc[c]["cells"] for c in big], float),
             cells_small=np.asarray([kdc[c]["cells"] for c in small], float),
             n_for_sb=float(n), json=j)
    d["constructs"] = z[f"{n}|constructs"].astype(str)
    d["cos_to_mean"] = z[f"{n}|cos_to_mean"]
    d["cross"] = [(d["constructs"], z[f"{n}|cross_resid"], z[f"{n}|cross_uncent"])]
    d["cx_target"] = np.asarray([c.split("||")[0] if "||" in c else t
                                 for c, t in zip(d["constructs"],
                                                 np.repeat(tgts, 2))], str)
    d["cx_lfc"] = np.asarray([kdc[c]["lfc"] for c in d["constructs"]], float)
    d["cx_kfrac"] = np.asarray([kdc[c]["kfrac"] for c in d["constructs"]], float)
    # each construct's OWN split: big -> (a), small -> (a2)
    d["cx_own"] = np.empty(len(d["constructs"]))
    d["cx_own"][0::2] = d["a"]; d["cx_own"][1::2] = d["a2"]
    d["cx_partner"] = np.repeat(d["b"], 2)
    return d


def load_broad(key: str):
    j = json.load(open(W / "broad_296.json"))
    if key not in j["budgets"]:
        return None
    z = np.load(W / "broad_296.npz", allow_pickle=True)
    rows = j["budgets"][key]["rows"]
    lab = np.asarray([f"{r['line']}|{r['target']}" for r in rows])
    d = dict(library="Broad-16", unit=f"T={key} UMIs/side", targets=lab,
             a=z[f"{key}|a"], a2=np.nanmean(np.vstack([z[f"{key}|a1"], z[f"{key}|a2"]]), 0),
             b=z[f"{key}|b"], x=z[f"{key}|x"],
             pds_a=z[f"{key}|pds_a"], pds_b=z[f"{key}|pds_b"], pds_x=z[f"{key}|pds_x"],
             c_big=np.asarray([f"{r['line']}|{r['target']}|v1" for r in rows]),
             c_small=np.asarray([f"{r['line']}|{r['target']}|v2" for r in rows]),
             kd_big=z[f"{key}|kd_v1"], kd_small=z[f"{key}|kd_v2"],
             kf_big=z[f"{key}|kf_v1"], kf_small=z[f"{key}|kf_v2"],
             cells_big=np.asarray([r["ncell"] for r in rows], float),
             cells_small=np.asarray([r["ncell"] for r in rows], float),
             n_for_sb=None, json=j,
             a1=z[f"{key}|a1"], a2_raw=z[f"{key}|a2"])
    cz = np.load(W / f"broad_296_cross_{key}.npz", allow_pickle=True)
    blocks, i = [], 0
    while f"{i}" in cz:
        n = cz[f"{i}"].shape[0]
        blocks.append((None, cz[f"{i}"], cz[f"u{i}"], n))
        i += 1
    labels = cz["labels"].astype(str)
    out, off = [], 0
    for _, R, U, n in blocks:
        out.append((labels[off:off + n], R, U)); off += n
    d["cross"] = out
    d["constructs"] = labels
    d["cos_to_mean"] = z[f"{key}|cos_to_mean"]
    kdj = j["kd"]
    d["cx_lfc"] = np.asarray([kdj[f"{c.split('|')[0]}|{c.split('|')[1]}|{c.split('|')[2]}"]["lfc"]
                              for c in labels], float)
    d["cx_kfrac"] = np.asarray([kdj[f"{c.split('|')[0]}|{c.split('|')[1]}|{c.split('|')[2]}"]["kfrac"]
                                for c in labels], float)
    d["cx_target"] = np.asarray(["|".join(c.split("|")[:2]) for c in labels])
    own = {f"{r['line']}|{r['target']}|v1": r["a1"] for r in rows}
    own.update({f"{r['line']}|{r['target']}|v2": r["a2"] for r in rows})
    part = {f"{r['line']}|{r['target']}|v1": r["b"] for r in rows}
    part.update({f"{r['line']}|{r['target']}|v2": r["b"] for r in rows})
    d["cx_own"] = np.asarray([own.get(c, np.nan) for c in labels], float)
    d["cx_partner"] = np.asarray([part.get(c, np.nan) for c in labels], float)
    return d


def load_cd4(key: str, cond: str):
    j = json.load(open(W / "cd4_296.json"))
    if key not in j["budgets"] or cond not in j["budgets"][key]:
        return None
    z = np.load(W / "cd4_296.npz", allow_pickle=True)
    r = j["budgets"][key][cond]
    tg = z[f"{key}|{cond}|targets"].astype(str)
    kdp = z[f"{key}|{cond}|kd_pair"]; kfp = z[f"{key}|{cond}|kf_pair"]
    d = dict(library="CD4-GWCD4i", unit=f"T={key} UMIs/side, {cond}", targets=tg,
             a=z[f"{key}|{cond}|a"], a2=z[f"{key}|{cond}|a2"], b=z[f"{key}|{cond}|b"],
             x=z[f"{key}|{cond}|x"],
             pds_a=z[f"{key}|{cond}|pds_a"], pds_b=z[f"{key}|{cond}|pds_b"],
             pds_x=z[f"{key}|{cond}|pds_x"],
             c_big=np.asarray([r["guides"][t][0] for t in tg]),
             c_small=np.asarray([r["guides"][t][1] for t in tg]),
             kd_big=kdp[:, 0], kd_small=kdp[:, 1], kf_big=kfp[:, 0], kf_small=kfp[:, 1],
             cells_big=np.full(len(tg), np.nan), cells_small=np.full(len(tg), np.nan),
             n_for_sb=None, json=j)
    cz = np.load(W / f"cd4_296_cross_{key}_{cond}.npz", allow_pickle=True)
    blocks, i = [], 0
    while f"lab{i}" in cz:
        blocks.append((cz[f"lab{i}"].astype(str), cz[f"r{i}"], cz[f"u{i}"]))
        i += 1
    d["cross"] = blocks
    labs = np.concatenate([b[0] for b in blocks]) if blocks else np.asarray([])
    d["constructs"] = labs
    d["cx_target"] = np.asarray([s.split("|")[0] for s in labs])
    slot = np.asarray([int(s.split("|s")[1]) for s in labs])
    kdmap = {t: (kdp[i, 0], kdp[i, 1]) for i, t in enumerate(tg)}
    kfmap = {t: (kfp[i, 0], kfp[i, 1]) for i, t in enumerate(tg)}
    d["cx_lfc"] = np.asarray([kdmap.get(t, (np.nan, np.nan))[s] for t, s in zip(d["cx_target"], slot)])
    d["cx_kfrac"] = np.asarray([kfmap.get(t, (np.nan, np.nan))[s] for t, s in zip(d["cx_target"], slot)])
    own0 = dict(zip(tg, d["a"])); own1 = dict(zip(tg, d["a2"])); pb = dict(zip(tg, d["b"]))
    d["cx_own"] = np.asarray([(own0 if s == 0 else own1).get(t, np.nan)
                              for t, s in zip(d["cx_target"], slot)], float)
    d["cx_partner"] = np.asarray([pb.get(t, np.nan) for t in d["cx_target"]], float)
    d["cos_to_mean"] = np.asarray([np.nan] * len(labs))
    ctm_lab = z[f"{key}|{cond}|cos_to_mean_lab"].astype(str)
    ctm = z[f"{key}|{cond}|cos_to_mean"]
    m = dict(zip(ctm_lab, ctm))
    d["cos_to_mean"] = np.asarray([m.get(s, np.nan) for s in labs], float)
    return d


# ---------------------------------------------------------------------------- stratification
def strat_block(d: dict, cut_name: str, sel_both: np.ndarray, sel_weak: np.ndarray) -> dict:
    rows = {}
    for nm, sel in (("ALL", np.ones(len(d["a"]), bool)), ("both >= cut", sel_both),
                    ("weaker < cut", sel_weak)):
        if sel.sum() == 0:
            continue
        rows[nm] = X.summarise(nm, sel, d, d.get("n_for_sb"))
    if "both >= cut" in rows and "weaker < cut" in rows:
        rows["diff"] = X.two_sample_ratio_diff(d["b"][sel_both], d["a"][sel_both],
                                               d["b"][sel_weak], d["a"][sel_weak])
    return dict(cut=cut_name, rows=rows)


def analyse(d: dict, out: dict) -> None:
    lib, unit = d["library"], d["unit"]
    kdw = np.minimum(np.abs(d["kd_big"]), np.abs(d["kd_small"]))
    kfw = np.minimum(d["kf_big"], d["kf_small"])
    cons_lfc = np.abs(np.concatenate([d["kd_big"], d["kd_small"]]))
    med_cut = float(np.nanmedian(cons_lfc))
    cuts = {"C1 |lfc| >= 1.0": (kdw >= CUT_ABS, kdw < CUT_ABS),
            f"C2 |lfc| >= library median {med_cut:.3f}": (kdw >= med_cut, kdw < med_cut),
            "C3 knockdown fraction >= 0.5": (kfw >= CUT_FRAC, kfw < CUT_FRAC)}
    key = f"{lib} [{unit}]"
    out[key] = dict(library=lib, unit=unit, k=int(len(d["a"])), median_cut=med_cut,
                    kd_weak_median=float(np.nanmedian(kdw)),
                    cons_lfc_median=float(np.nanmedian(cons_lfc)), strata={})
    P(f"\n=== {key}  ({len(d['a'])} pairs; weaker-construct |lfc| median {np.nanmedian(kdw):.3f}; "
      f"library median single-construct |lfc| {med_cut:.3f}) ===")
    for cname, (sb, sw) in cuts.items():
        sb = sb & np.isfinite(d["a"]) & np.isfinite(d["b"])
        sw = sw & np.isfinite(d["a"]) & np.isfinite(d["b"])
        blk = strat_block(d, cname, sb, sw)
        out[key]["strata"][cname] = blk
        P(f"  -- {cname}")
        for nm, r in blk["rows"].items():
            if nm == "diff":
                continue
            rd = ("   n/a" if r["rho_dis"] is None else f"{r['rho_dis']:+.3f}")
            P(f"     {nm:14s} k {r['k']:4d}  (a) {r['a']['point']:+.4f}  (a2) {r['a2']['point']:+.4f}  "
              f"(b) {r['b']['point']:+.4f}  (x) {r['x']['point']:+.4f}  rho_guide "
              f"{X.fmt_rho(r['rho_guide'])}  rho_dis {rd}  pds {r['pds_a']:.3f}/{r['pds_b']:.3f}/"
              f"{r['pds_x']:.3f}")
        if "diff" in blk["rows"]:
            dd = blk["rows"]["diff"]
            P(f"     STRATUM BOOTSTRAP  rho(both) - rho(weak) = {dd['point']:+.3f} "
              f"[{dd['lo']:+.3f}, {dd['hi']:+.3f}]  P(<=0) {dd['p_le0']:.4f}  "
              f"-> P1 rule {'PASS' if dd['lo'] > 0 else 'FAIL'}")
    # P4: knockdown of the weaker construct vs the per-target disattenuated ratio
    with np.errstate(invalid="ignore", divide="ignore"):
        rdis = np.where((d["a"] > 0) & (d["a2"] > 0), d["b"] / np.sqrt(d["a"] * d["a2"]), np.nan)
    s1, p1 = X.spearman(kdw, rdis)
    s2, p2 = X.spearman(kdw, d["b"])
    s3, p3 = X.spearman(kdw, d["a"])
    out[key]["P4"] = dict(spearman_rho_dis=s1, p_rho_dis=p1, n_rho_dis=int(np.isfinite(rdis).sum()),
                          spearman_b=s2, p_b=p2, spearman_a=s3, p_a=p3, rule="> +0.30",
                          outcome="PASS" if s1 > 0.30 else "FAIL")
    P(f"  P4  Spearman(|lfc| of the weaker construct, per-target rho_dis) {s1:+.3f} (p {p1:.3g}, "
      f"n {int(np.isfinite(rdis).sum())}) -> {out[key]['P4']['outcome']};  with (b) {s2:+.3f} "
      f"(p {p2:.3g});  with (a), the cell split, {s3:+.3f} (p {p3:.3g})")
    out[key]["rho_dis_per_target"] = rdis.tolist()
    out[key]["kd_weak"] = kdw.tolist()


# ---------------------------------------------------------------------------- failing constructs (P5)
def failing(d: dict, out: dict) -> None:
    key = f"{d['library']} [{d['unit']}]"
    lfc = np.abs(d["cx_lfc"])
    cls = np.where(lfc < 0.5, "fail", np.where(lfc >= 1.0, "work", "mid"))
    own, part, ctm = d["cx_own"], d["cx_partner"], d["cos_to_mean"]
    r = dict(n_fail=int((cls == "fail").sum()), n_work=int((cls == "work").sum()),
             n_mid=int((cls == "mid").sum()))
    for c in ("fail", "work"):
        m = cls == c
        r[f"own_split_{c}"] = L.boot_median(own[m])
        r[f"cos_partner_{c}"] = L.boot_median(part[m])
        r[f"cos_to_libmean_{c}"] = L.boot_median(ctm[m])
    r["own_split_diff"] = X.two_sample_median_diff(own[cls == "work"], own[cls == "fail"])
    r["cos_to_libmean_diff"] = X.two_sample_median_diff(ctm[cls == "fail"], ctm[cls == "work"])
    # cross-construct agreement between DIFFERENT genes, by class pair, in both spaces
    tgt = d["cx_target"]
    off = 0
    acc = {sp: {"ff": [], "ww": [], "fw": []} for sp in ("resid", "uncent")}
    for labs, R, U in d["cross"]:
        n = len(labs)
        sl = slice(off, off + n)
        cl, tg = cls[sl], tgt[sl]
        diff = tg[:, None] != tg[None, :]
        iu = np.triu(np.ones((n, n), bool), 1) & diff
        for sp, Mx in (("resid", R), ("uncent", U)):
            for tag, msk in (("ff", (cl == "fail")[:, None] & (cl == "fail")[None, :]),
                             ("ww", (cl == "work")[:, None] & (cl == "work")[None, :]),
                             ("fw", ((cl == "fail")[:, None] & (cl == "work")[None, :])
                                    | ((cl == "work")[:, None] & (cl == "fail")[None, :]))):
                v = np.asarray(Mx)[iu & msk]
                if v.size:
                    acc[sp][tag].append(v)
        off += n
    def cap(v, m=40000, seed=L.SEED):
        """Cross-GENE construct pairs run to >1e6 per library; the bootstrap is over a seeded
        subsample of at most `m` pairs (the median itself is computed on ALL of them)."""
        v = np.asarray(v, float); v = v[np.isfinite(v)]
        if v.size <= m:
            return v, v.size
        return np.random.default_rng(seed).choice(v, m, replace=False), int(v.size)

    for sp in acc:
        raw = {}
        for tag in acc[sp]:
            v = np.concatenate(acc[sp][tag]) if acc[sp][tag] else np.asarray([])
            raw[tag] = v
            if v.size:
                s, n_all = cap(v)
                bm = L.boot_median(s, n=4000)
                bm["point"] = float(np.median(v[np.isfinite(v)]))
                bm["n"] = n_all
                r[f"cross_{sp}_{tag}"] = bm
            else:
                r[f"cross_{sp}_{tag}"] = dict(point=np.nan, lo=np.nan, hi=np.nan, n=0)
        if raw["ff"].size and raw["ww"].size:
            f_, _ = cap(raw["ff"]); w_, _ = cap(raw["ww"])
            d_ = X.two_sample_median_diff(f_, w_, n=4000)
            d_["point"] = float(np.median(raw["ff"]) - np.median(raw["ww"]))
            r[f"cross_{sp}_ff_minus_ww"] = d_
    out[key] = r
    P(f"\n--- failing constructs, {key}: {r['n_fail']} failing (|lfc| < 0.5), {r['n_work']} working "
      f"(>= 1.0), {r['n_mid']} between")
    P(f"      own split-half   failing {L.fmt(r['own_split_fail'])}   working {L.fmt(r['own_split_work'])}"
      f"   (working - failing {L.fmt(r['own_split_diff'])})")
    P(f"      cos to its working PARTNER      failing {L.fmt(r['cos_partner_fail'])}   "
      f"working {L.fmt(r['cos_partner_work'])}")
    P(f"      cos to the LIBRARY MEAN (uncentred)  failing {L.fmt(r['cos_to_libmean_fail'])}   "
      f"working {L.fmt(r['cos_to_libmean_work'])}   (failing - working {L.fmt(r['cos_to_libmean_diff'])})")
    for sp in ("resid", "uncent"):
        P(f"      cross-GENE agreement [{sp:7s}]  fail x fail {L.fmt(r[f'cross_{sp}_ff'])} "
          f"(n {r[f'cross_{sp}_ff']['n']})   work x work {L.fmt(r[f'cross_{sp}_ww'])} "
          f"(n {r[f'cross_{sp}_ww']['n']})   diff {L.fmt(r.get(f'cross_{sp}_ff_minus_ww', dict(point=np.nan, lo=np.nan, hi=np.nan)))}")


def main() -> None:
    t0 = time.time()
    P(f"296 strat296.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    libs: list[dict] = []
    for n in (25, 50):
        d = load_k562_like("k562_296", "K562-GWPS", n)
        if d:
            libs.append(d)
    for wh in ("jurkat", "hepg2"):
        for n in (25, 50):
            d = load_k562_like(f"nadig_{wh}_296", f"Nadig-{wh}", n)
            if d:
                libs.append(d)
    for key in ("2e+05", "5e+05"):
        d = load_broad(key)
        if d:
            libs.append(d)
    for key in ("1e+05", "2e+05"):
        for c in ("Rest", "Stim8hr", "Stim48hr"):
            d = load_cd4(key, c)
            if d:
                libs.append(d)
    P(f"  {len(libs)} library x budget cells loaded")

    out: dict = {"strata": {}, "failing": {}}
    for d in libs:
        analyse(d, out["strata"])
    for d in libs:
        failing(d, out["failing"])

    # ------------------------------------------------------------------ efficacy296.csv
    with open(W / "efficacy296.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["library", "unit", "target", "construct", "role", "cells", "lfc", "abs_lfc",
                    "knockdown_frac", "a", "a2", "b", "x", "rho_guide_pair", "rho_dis_pair",
                    "pds_a", "pds_b", "pds_x", "own_split", "cos_to_library_mean"])
        for d in libs:
            with np.errstate(invalid="ignore", divide="ignore"):
                rp = np.where(np.abs(d["a"]) > 1e-12, d["b"] / d["a"], np.nan)
                rd = np.where((d["a"] > 0) & (d["a2"] > 0), d["b"] / np.sqrt(d["a"] * d["a2"]), np.nan)
            ctm = dict(zip(d["constructs"], d["cos_to_mean"]))
            own = dict(zip(d["constructs"], d["cx_own"]))
            for i, t in enumerate(d["targets"]):
                for role, c, kd, kf, ncl in (("big/v1/slot0", d["c_big"][i], d["kd_big"][i],
                                              d["kf_big"][i], d["cells_big"][i]),
                                             ("small/v2/slot1", d["c_small"][i], d["kd_small"][i],
                                              d["kf_small"][i], d["cells_small"][i])):
                    w.writerow([d["library"], d["unit"], t, c, role,
                                "" if not np.isfinite(ncl) else int(ncl),
                                f"{kd:.6f}" if np.isfinite(kd) else "",
                                f"{abs(kd):.6f}" if np.isfinite(kd) else "",
                                f"{kf:.6f}" if np.isfinite(kf) else "",
                                f"{d['a'][i]:.6f}", f"{d['a2'][i]:.6f}" if np.isfinite(d["a2"][i]) else "",
                                f"{d['b'][i]:.6f}", f"{d['x'][i]:.6f}" if np.isfinite(d["x"][i]) else "",
                                f"{rp[i]:.6f}" if np.isfinite(rp[i]) else "",
                                f"{rd[i]:.6f}" if np.isfinite(rd[i]) else "",
                                f"{d['pds_a'][i]:.4f}", f"{d['pds_b'][i]:.4f}", f"{d['pds_x'][i]:.4f}",
                                f"{own.get(c, np.nan):.6f}" if np.isfinite(own.get(c, np.nan)) else "",
                                f"{ctm.get(c, np.nan):.6f}" if np.isfinite(ctm.get(c, np.nan)) else ""])
    P(f"\n  wrote work/efficacy296.csv "
      f"({sum(2 * len(d['targets']) for d in libs)} (library, unit, target, construct) rows)")

    out["seconds"] = round(time.time() - t0, 1)
    L.Logger.dump(P, out, W / "strat296.json")
    P(f"\ndone {out['seconds']}s")
    P.save()


if __name__ == "__main__":
    main()
