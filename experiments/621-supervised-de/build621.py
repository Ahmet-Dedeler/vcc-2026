#!/usr/bin/env python3
"""621 — supervised DE direction + membership model (chill-style), leave-one-surface-out, scored by 636.

Council portfolio item (2), 2026-10-02 (owner 7772dc; spec 21dd4c; instrument d4cd5f's 636).

Per (target t, gene g) on a held-out surface s (h1 / hct116 / hek293t / k562; gates from 480's gates480.py), the
features come ONLY from banks 636 allows on s (no own line, no same-study sibling, no H1 bank on h1):
  per-target source evidence  z_d(t,g) = delta_d(t,g) / sigma_d(g)   (sigma = robust SD of the bank's rows for g)
  per study family f          Z_f = Stouffer over the family's banks that carry t
  family-balanced consensus   Zbal = mean_f Z_f (equal weight per study family, chill's "family-balanced")
  agreement, max |Z_f|, number of families
  gene-level programme        logit P_down_f(g) per family, family-balanced mean (gene aggregate; context only)
  gene responsiveness         fraction of bank rows with |z| > 3 (family-balanced)
  recipient expression        log10 control CPM of g in the held-out surface's OWN controls (available on the day)
  target strength             mean over carried families of mean_g Z_f(t,g)^2
Heads (HistGradientBoosting): direction P(down | gate) trained on gate pairs; membership P(gate | tested).
Training surfaces for held-out s exclude s's targets (no per-target label leakage) and, for the Orion pair, the sibling
surface (same study). Each training surface gets equal total weight.
Outputs work/pred-<s>.npz (targets, genes, P_down, member) + json sidecars, plus baseline rows on identical axes
(consensus sign, programme). Analysis only: no build, no vcc, no cell-eval2.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import BANK_FORBIDDEN, Inst  # noqa: E402

BANKS = {  # name -> (dir, study family)
    "k562_gwps": (ROOT / "data/donors/k562_gwps", "weissman"),
    "k562_essential": (ROOT / "data/donors/k562_essential", "weissman"),
    "rpe1": (ROOT / "data/donors/rpe1", "weissman"),
    "nadig_jurkat": (ROOT / "data/donors/nadig_jurkat", "weissman"),
    "nadig_hepg2": (ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2", "weissman"),
    "viperturb_k562": (ROOT / "data/donors/viperturb_k562", "satija"),
    "orion_hct116": (ROOT / "data/donors/orion_hct116", "xaira"),
    "orion_hek293t": (ROOT / "data/donors/orion_hek293t", "xaira"),
    "h1_2025_full": (ROOT / "data/donors/h1_2025_full", "arc"),
}
FAMILIES = ["weissman", "satija", "xaira", "arc"]
TRAIN_FOR = {"h1": ["hct116", "hek293t", "k562"], "hct116": ["h1", "k562"], "hek293t": ["h1", "k562"],
             "k562": ["h1", "hct116", "hek293t"]}
N_SUB = 3000
RNG = np.random.default_rng(6210)


def read_list(p: Path) -> list[str]:
    return [x.strip() for x in p.read_text().split("\n") if x.strip()]


def bank_stats(name: str) -> dict:
    """Gene axis, perts, memmapped delta, robust per-gene sigma, programme P_down and responsiveness."""
    d, fam = BANKS[name]
    genes = read_list(d / "genes.txt")
    perts = read_list(d / "perts.txt")
    D = np.load(d / "delta.npy", mmap_mode="r")
    assert D.shape == (len(perts), len(genes)), (name, D.shape, len(perts), len(genes))
    rows = np.sort(RNG.choice(len(perts), size=min(N_SUB, len(perts)), replace=False))
    S = np.asarray(D[rows], dtype=np.float32)
    med = np.median(S, axis=0)
    mad = np.median(np.abs(S - med), axis=0) * 1.4826
    sigma = np.where(mad > 1e-6, mad, np.nan).astype(np.float32)
    z = S / sigma
    pdown = ((S < 0).sum(0) + 0.5) / ((S != 0).sum(0) + 1.0)
    resp = (np.abs(z) > 3).mean(0)
    del S, z
    return dict(name=name, fam=fam, genes=genes, gpos={g: i for i, g in enumerate(genes)},
                ppos={p: i for i, p in enumerate(perts)}, D=D, sigma=sigma, pdown=pdown, resp=resp)


def surface_features(inst: Inst, s: str, B: dict) -> dict:
    T, G = inst.axes(s)
    nT, nG = len(T), len(G)
    allowed = [b for b in B if s not in BANK_FORBIDDEN.get(b, set())]
    fam_sum = {f: np.zeros((nT, nG), np.float32) for f in FAMILIES}
    fam_n = {f: np.zeros((nT, nG), np.float32) for f in FAMILIES}
    prog = {f: [] for f in FAMILIES}
    resp = {f: [] for f in FAMILIES}
    for b in allowed:
        st = B[b]
        gidx = np.array([st["gpos"].get(g, -1) for g in G])
        okg = gidx >= 0
        pv = np.full(nG, np.nan, np.float32); pv[okg] = st["pdown"][gidx[okg]]
        rv = np.full(nG, np.nan, np.float32); rv[okg] = st["resp"][gidx[okg]]
        prog[st["fam"]].append(pv); resp[st["fam"]].append(rv)
        tidx = [(i, st["ppos"][t]) for i, t in enumerate(T) if t in st["ppos"]]
        if not tidx:
            continue
        ti = np.array([a for a, _ in tidx]); ri = np.array([r for _, r in tidx])
        order = np.argsort(ri)
        rows = np.asarray(st["D"][ri[order]], dtype=np.float32)
        Z = np.full((len(ti), nG), np.nan, np.float32)
        Z[:, okg] = rows[:, gidx[okg]] / st["sigma"][gidx[okg]]
        ti_sorted = ti[order]
        fin = np.isfinite(Z)
        fam_sum[st["fam"]][ti_sorted] += np.where(fin, Z, 0)
        fam_n[st["fam"]][ti_sorted] += fin
    Zf = {f: np.where(fam_n[f] > 0, fam_sum[f] / np.sqrt(np.maximum(fam_n[f], 1)), np.nan) for f in FAMILIES}
    stackZ = np.stack([Zf[f] for f in FAMILIES])  # F x T x G
    nfam = np.isfinite(stackZ).sum(0).astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        Zbal = np.nanmean(stackZ, axis=0)
        agree = np.nanmean(np.sign(stackZ) == np.sign(Zbal)[None], axis=0)
        agree = np.where(nfam > 0, agree, np.nan)
        maxabs = np.nanmax(np.abs(stackZ), axis=0)
        famprog = np.stack([np.nanmean(np.stack(prog[f]), 0) if prog[f] else np.full(nG, np.nan) for f in FAMILIES])
        famresp = np.stack([np.nanmean(np.stack(resp[f]), 0) if resp[f] else np.full(nG, np.nan) for f in FAMILIES])
        logit = np.log(famprog / (1 - famprog))
        prog_bal = np.nanmean(logit, axis=0)
        resp_bal = np.nanmean(famresp, axis=0)
        tstrength = np.nanmean(np.where(np.isfinite(Zbal), Zbal ** 2, np.nan), axis=1)
    z480 = np.load(inst_path(s), allow_pickle=True)
    lcpm = np.log10(np.asarray(z480["ctrl_cpm"], dtype=np.float64) + 1e-3).astype(np.float32)
    F = {
        "Zbal": Zbal, "absZbal": np.abs(Zbal), "agree": agree, "maxabs": maxabs, "nfam": nfam,
        "prog_bal": np.broadcast_to(prog_bal, (nT, nG)), "resp_bal": np.broadcast_to(resp_bal, (nT, nG)),
        "lcpm": np.broadcast_to(lcpm, (nT, nG)), "tstrength": np.broadcast_to(tstrength[:, None], (nT, nG)),
        "Zxprog": Zbal * np.broadcast_to(prog_bal, (nT, nG)),
    }
    for f in FAMILIES:
        F[f"Z_{f}"] = Zf[f]
    return dict(T=T, G=G, F=F, allowed=allowed)


def inst_path(s: str) -> Path:
    from inst636 import G480, SURF
    return G480 / SURF[s]


FEATS = ["Zbal", "absZbal", "agree", "maxabs", "nfam", "prog_bal", "resp_bal", "lcpm", "tstrength", "Zxprog",
         "Z_weissman", "Z_satija", "Z_xaira", "Z_arc"]


def design(SF: dict, inst: Inst, s: str, exclude_targets: set[str], which: str, max_neg: int = 200_000):
    D = inst.S[s]
    keepT = np.array([t not in exclude_targets for t in SF["T"]])
    mask = D["gate"] if which == "dir" else D["tested"]
    mask = mask & keepT[:, None]
    ii, jj = np.nonzero(mask)
    if which == "mem":
        lab = D["gate"][ii, jj]
        neg = np.flatnonzero(~lab)
        if len(neg) > max_neg:
            drop = RNG.choice(neg, size=len(neg) - max_neg, replace=False)
            keep = np.ones(len(ii), bool); keep[drop] = False
            ii, jj = ii[keep], jj[keep]
    X = np.stack([np.asarray(SF["F"][f])[ii, jj] for f in FEATS], axis=1).astype(np.float32)
    if which == "dir":
        y = (D["lfc"][ii, jj] < 0).astype(np.int8)
        w = np.ones(len(y), np.float64)
    else:
        y = D["gate"][ii, jj].astype(np.int8)
        pos = y.sum(); negn = len(y) - pos
        w = np.where(y == 1, 0.5 / max(pos, 1), 0.5 / max(negn, 1)) * len(y)
    return X, y, w


def main() -> None:
    from sklearn.ensemble import HistGradientBoostingClassifier
    t0 = time.time()
    inst = Inst()
    B = {b: bank_stats(b) for b in BANKS}
    print(f"[621] bank stats {time.time() - t0:.0f}s", flush=True)
    SF = {s: surface_features(inst, s, B) for s in ["h1", "hct116", "hek293t", "k562"]}
    print(f"[621] surface features {time.time() - t0:.0f}s", flush=True)
    code_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:12]
    for s in SF:
        heads = {}
        for which in ("dir", "mem"):
            Xs, ys, ws = [], [], []
            for tr in TRAIN_FOR[s]:
                X, y, w = design(SF[tr], inst, tr, set(SF[s]["T"]), which)
                if len(y) == 0:
                    continue
                Xs.append(X); ys.append(y); ws.append(w / w.sum())  # equal total weight per training surface
            X = np.concatenate(Xs); y = np.concatenate(ys); w = np.concatenate(ws) * len(np.concatenate(ys))
            clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
                                                 min_samples_leaf=200, l2_regularization=1.0, random_state=6210)
            clf.fit(X, y, sample_weight=w)
            heads[which] = clf
        T, G = SF[s]["T"], SF[s]["G"]
        Xall = np.stack([np.asarray(SF[s]["F"][f]).reshape(-1) for f in FEATS], axis=1).astype(np.float32)
        P_down = heads["dir"].predict_proba(Xall)[:, 1].reshape(len(T), len(G)).astype(np.float32)
        member = heads["mem"].predict_proba(Xall)[:, 1].reshape(len(T), len(G)).astype(np.float32)
        Zb = np.asarray(SF[s]["F"]["Zbal"])
        cons = np.where(np.isfinite(Zb), (Zb < 0).astype(np.float32), np.nan).astype(np.float32)
        prog = (1 / (1 + np.exp(-np.asarray(SF[s]["F"]["prog_bal"])))).astype(np.float32)
        np.savez_compressed(WORK / f"pred-{s}.npz", targets=T, genes=G, P_down=P_down, member=member,
                            P_down_consensus=cons, P_down_programme=prog,
                            Zbal=Zb.astype(np.float32), agree=np.asarray(SF[s]["F"]["agree"], np.float32),
                            nfam=np.asarray(SF[s]["F"]["nfam"], np.float32))
        side = {"surface": s, "model": "HGB dir+mem v1", "code_sha12": code_sha, "train_surfaces": TRAIN_FOR[s],
                "excluded_targets": "all of the held-out surface's targets", "banks_used": SF[s]["allowed"],
                "features": FEATS, "n_sub_rows_per_bank": N_SUB}
        (WORK / f"pred-{s}.json").write_text(json.dumps(side, indent=1))
        print(f"[621] {s}: wrote pred ({len(T)}x{len(G)}), banks {SF[s]['allowed']}, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
