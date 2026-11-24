#!/usr/bin/env python3
"""708 Part B (and Part A's confirmation): price the late switches on a panel's PLOLO folds from dump708's dumps.

Per fold (work/dump708-<tag>-<fold>.npz), every arm is a per-target weight on the composition, then 522's render:
key_scale (energy-preserving), q 1.75 row sharpening on each frame's gate, minnorm render, ballast KMAX 16, scored on the
fold's truth (fold522.Oct's arithmetic). The incumbent is the percentile cap P 96.647 on fin = h2c1.
Arms vs the incumbent:
  r    522's θ 0.25 floor-preserving rider: w' = max(0.5, w exp(0.25 z)) above the floor, w at the floor. z = mean over the
       fold's two SOURCE donors of crank(robust strength over the panel's targets) (522's src feature, fold form).
  s    707's swap: on the top z tertile (targets with source rows), fin's non-panel row blended toward the sources'
       unit mean U = Σ_j r_j/|r_j| at β 0.25, norm-preserving (the composition analogue of swap707.py).
  rs   both.
  s|r  the swap's MARGINAL price on top of the rider: (rs) priced against (r) (708 AMENDMENT 3: when the rider ships, the
       swap is judged on top of it), with its own leave-top-3-out.
  A:<design>  fin_<design> (dump708 --designs) vs fin, under the incumbent key (708 AMENDMENT 1).
Price per fold, in 522 board units: 0.0710 · (X_arm/X_inc − 1) · ret + 0.1245 · Δpds_lit, with k 1 and ret 1. X is the
frame-averaged cos_agg. Leave-top-3-out (236 / robust438.top3): the union over folds of each fold's 3 targets with the
largest Σ_frames Δx is removed, and every arm is re-scored on the rest.
Rules (PREREGISTRATION): 438's move rule = mean ≥ +0.0010, 3/3 folds > 0, and mean after leave-top-3-out > 0. The veto
rule = mean ≤ −0.0005 with ≥ 2/3 folds < 0.
usage: uv run python experiments/708-def-panel-repricing/score708.py <tag> [--designs h2.5c1.5]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/522-october-bulk-key"))
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
import numpy as np  # noqa: E402

import fold522 as FD  # noqa: E402
import lib438 as L  # noqa: E402
from lib438 import SLOTS  # noqa: E402
from plolo import TS, profile_from_effect  # noqa: E402
from strength522 import BANKS, robust_strength  # noqa: E402

FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
GENES = None


class Fold:
    """fold522.Oct's arithmetic, loaded from a dump708 file (no 546 known answer: dump708 checked its own)."""

    def __init__(self, path: Path):
        z = np.load(path)
        self.targets = [str(t) for t in z["targets"]]
        self.sources = [str(s) for s in z["sources"]]
        self.fin, self.fPP = z["fin"], z["fPP"]
        self.extra = {k[4:]: z[k] for k in z.files if k.startswith("fin_")}
        self.kv, self.carried = z["kv"], z["carried"]
        self.scored, self.pp_pos = z["scored"], z["panel_positions"]
        self.truth_full = np.asarray(z["truth"], np.float64)
        self.truth = self.truth_full[:, self.scored]
        self.betas = {c: z[f"beta_{c}"] for c in SLOTS}
        self.frames = {c: (z[f"b_{c}"], z[f"gate_{c}"]) for c in SLOTS}
        self.nonpanel = np.ones(self.fin.shape[1], bool)
        self.nonpanel[self.pp_pos] = False
        self.n = len(self.targets)
        self.P_DE = {c: [profile_from_effect(self.fPP[i], b, g, "minnorm") for i in range(self.n)]
                     for c, (b, g) in self.frames.items()}

    def score(self, w, fin=None, mask=None, q: float = 1.75):
        F0 = self.fin if fin is None else fin
        idx = np.arange(self.n) if mask is None else np.flatnonzero(mask)
        bulk = L.key_scale(F0[idx], np.asarray(w, float)[idx])[0]
        rea = {}
        for c, (b, gate) in self.frames.items():
            Fs = L.sharpen(bulk, float(q), gate, mode="row")
            R = np.empty_like(Fs)
            for r, i in enumerate(idx):
                p_bk = profile_from_effect(Fs[r], b, gate, "minnorm")
                p_re, _, _, _, _ = L.realise_k(self.P_DE[c][i], p_bk, 16)
                R[r] = np.log1p(TS * p_re) - b
            rea[c] = R
        s, arr = L.row_score(rea, self.scored, self.truth[idx], self.betas, self.nonpanel)
        return s, arr, idx


def zsrc_fold(targets, panel, sources):
    cols = []
    for k in sources:
        s, _, _ = robust_strength(BANKS[k], list(targets), set(panel))
        r = FD.crank(s)
        r[~np.isfinite(s)] = np.nan
        cols.append(r)
    with np.errstate(all="ignore"):
        m = np.nanmean(np.vstack(cols), 0)
    return np.where(np.isfinite(m), m, 0.0)


def unit_mean(targets, sources, nonpanel):
    global GENES
    if GENES is None:
        import csv
        GENES = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    U = np.zeros((len(targets), len(GENES)))
    for k in sources:
        D = Path(BANKS[k])
        dg = (D / "genes.txt").read_text().split()
        dp = (D / "perts.txt").read_text().split()
        delta = np.load(D / "delta.npy", mmap_mode="r")
        gpos = {g: i for i, g in enumerate(dg)}
        ppos = {p: i for i, p in enumerate(dp)}
        gc = np.array([gpos.get(g, -1) for g in GENES])
        ok = (gc >= 0) & nonpanel
        for i, t in enumerate(targets):
            if t in ppos:
                r = np.zeros(len(GENES))
                r[ok] = np.asarray(delta[ppos[t]], np.float64)[gc[ok]]
                nr = np.linalg.norm(r)
                if nr > 0:
                    U[i] += r / nr
    return U


def swap_fin(fin, U, sel, nonpanel, beta=0.25):
    F = fin.copy()
    nb = np.linalg.norm(fin[:, nonpanel], axis=1)
    nu = np.linalg.norm(U[:, nonpanel], axis=1)
    for i in np.flatnonzero(sel & (nb > 0) & (nu > 0)):
        m = (1 - beta) * fin[i, nonpanel] / nb[i] + beta * U[i, nonpanel] / nu[i]
        F[i, nonpanel] = nb[i] * m / np.linalg.norm(m)
    return F


def X_of(s):
    return float(s["cos_agg"])


def price(s_arm, s_inc):
    return 0.0710 * (X_of(s_arm) / X_of(s_inc) - 1.0) + 0.1245 * (s_arm["pds"] - s_inc["pds"])


def price_all(F_all, arms_all, inc_r, bm, label=""):
    """Every arm against its incumbent, then leave-top-3-out (union over folds, per arm), on the targets bm[f] keeps
    (None = all of them, the preregistered primary; 708 AMENDMENT 4: the floor-kept targets). With bm all None this is
    the pre-AMENDMENT-4 arithmetic call for call (the known answer: the primary json equals the old code's)."""
    out = {}
    base = {f: F_all[f].score(arms_all[f][0], mask=bm[f]) for f in FOLDS}

    def top3(S, sa, si):
        dx = (np.asarray(sa[1]["x"]) - np.asarray(si[1]["x"])).sum(0)
        return {S.targets[sa[2][i]] for i in np.argsort(-dx)[:3]}

    def lt3_keep(S, f, dx_top):
        k = np.array([t not in dx_top for t in S.targets])
        return k if bm[f] is None else (k & bm[f])

    def record(name, per, dx_top, lt3, armsc, incsc, gap):
        vals = np.array([per[f] for f in FOLDS])
        mean, npos, nneg = float(vals.mean()), int((vals > 0).sum()), int((vals < 0).sum())
        lt3m = float(np.mean([lt3[f] for f in FOLDS]))
        move = mean >= 0.0010 and npos == 3 and lt3m > 0
        veto = mean <= -0.0005 and nneg >= 2
        out[name] = dict(per_fold=per, mean=mean, n_pos=npos, lt3_mean=lt3m, top3_union=sorted(dx_top),
                         move_rule_438=bool(move), veto_rule=bool(veto),
                         x_rel={f: X_of(armsc[f][0]) / X_of(incsc[f][0]) - 1 for f in FOLDS},
                         dpds={f: armsc[f][0]["pds"] - incsc[f][0]["pds"] for f in FOLDS})
        print(f"{label}{name:12s} mean {mean:+.5f} ({npos}/3 > 0) per fold " + " / ".join(f"{per[f]:+.5f}" for f in FOLDS)
              + f" | lt3 {lt3m:+.5f}" + ("" if gap else f" | MOVE(438) {move} | VETO {veto}"), flush=True)

    for name in arms_all[FOLDS[0]][1]:
        per, dx_top, armsc = {}, set(), {}
        for f in FOLDS:
            S = F_all[f]
            w, fin = arms_all[f][1][name]
            sa = S.score(w, fin=fin, mask=bm[f])
            armsc[f] = sa
            per[f] = price(sa[0], base[f][0])
            dx_top |= top3(S, sa, base[f])
        lt3 = {}
        for f in FOLDS:
            S = F_all[f]
            keep = lt3_keep(S, f, dx_top)
            w, fin = arms_all[f][1][name]
            sa = S.score(w, fin=fin, mask=keep)
            si = S.score(arms_all[f][0], mask=keep)
            lt3[f] = price(sa[0], si[0])
        record(name, per, dx_top, lt3, armsc, base, False)
    # s|r: the swap on top of the rider ((rs) vs (r)), its own top-3 union and leave-top-3-out
    per, dx_top, armsc, incsc = {}, set(), {}, {}
    for f in FOLDS:
        S = F_all[f]
        sa = S.score(arms_all[f][1]["rs"][0], fin=arms_all[f][1]["rs"][1], mask=bm[f])
        si = S.score(inc_r[f], mask=bm[f])
        armsc[f], incsc[f] = sa, si
        per[f] = price(sa[0], si[0])
        dx_top |= top3(S, sa, si)
    lt3 = {}
    for f in FOLDS:
        S = F_all[f]
        keep = lt3_keep(S, f, dx_top)
        sa = S.score(arms_all[f][1]["rs"][0], fin=arms_all[f][1]["rs"][1], mask=keep)
        si = S.score(inc_r[f], mask=keep)
        lt3[f] = price(sa[0], si[0])
    record("s|r", per, dx_top, lt3, armsc, incsc, True)
    return out, base


def board_sized(out, cls):
    """708 AMENDMENT 4's information column: each arm's mean divided by the lab's same-panel class ratio (522 units)."""
    rows = {}
    for name, a in out.items():
        c = cls["arms"].get(name) or (cls["arms"].get("A:*") if name.startswith("A:") else None)
        if c is None:
            continue
        m = a["mean"]
        if "rho_522" in c:
            txt = f"{m / c['rho_522']:+.5f} (rho {c['rho_522']}{', provisional' if c.get('provisional') else ''})"
            if "rho_522_heavy_tail" in c:
                txt += f"; if D/E/F's tail is heavy {m / c['rho_522_heavy_tail']:+.5f} (rho {c['rho_522_heavy_tail']})"
        else:
            lo, hi = c["rho_522_range"]
            txt = f"{m / hi:+.5f} .. {m / lo:+.5f} (rho {lo}-{hi})"
        rows[name] = dict(cls=c["class"], text=txt)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--panel-bundle", type=Path, default=None, help="bundle whose pert_counts.csv is the panel (default: from dump708's json spec)")
    ap.add_argument("--floor-json", type=Path, default=None, help="708 AMENDMENT 4: floor708.py's day record; prices every arm on the floor-kept targets too")
    ap.add_argument("--dumps-dir", type=Path, default=HERE / "work", help="where dump708-<tag>-<fold>.npz/.json live (tests)")
    ap.add_argument("--out", type=Path, default=None, help="output json (default work/score708-<tag>.json)")
    ap.add_argument("--class-json", type=Path, default=HERE / "work" / "class708.json")
    a = ap.parse_args()
    t0 = time.time()
    res = {"tag": a.tag, "folds": {}}
    F_all, arms_all, inc_r = {}, {}, {}
    for f in FOLDS:
        p = a.dumps_dir / f"dump708-{a.tag}-{f}.npz"
        meta = json.loads(p.with_suffix(".json").read_text())
        bundle = Path(a.panel_bundle or meta["spec"]["bundle"])
        import csv
        panel = [r["target_gene"] for r in csv.DictReader(open(bundle / "pert_counts.csv")) if r["target_gene"] != "non-targeting"]
        S = Fold(p)
        wp = FD.wclip(S.kv, S.carried, 0.5, FD.pct_hi(S.kv, S.carried, FD.P_PROD))
        z = zsrc_fold(S.targets, panel, S.sources)
        wr = np.where(wp <= 0.5 + 1e-12, wp, np.maximum(0.5, wp * np.exp(0.25 * z)))
        U = unit_mean(S.targets, S.sources, S.nonpanel)
        have = np.linalg.norm(U, axis=1) > 0
        tier = np.full(S.n, -1)
        r_ = np.argsort(np.argsort(z[have], kind="stable"), kind="stable")
        tier[np.flatnonzero(have)] = np.minimum((r_ * 3) // max(int(have.sum()), 1), 2)
        Fs = swap_fin(S.fin, U, tier == 2, S.nonpanel)
        arms = {"r": (wr, None), "s": (wp, Fs), "rs": (wr, Fs)}
        arms.update({f"A:{d}": (wp, S.extra[d]) for d in S.extra})
        F_all[f] = S
        arms_all[f] = (wp, arms)
        inc_r[f] = wr
    out, base = price_all(F_all, arms_all, inc_r, {f: None for f in FOLDS})
    res["arms"] = out
    res["base"] = {f: {"cos_agg": base[f][0]["cos_agg"], "pds": base[f][0]["pds"]} for f in FOLDS}
    # 708 AMENDMENT 4: the floor-kept version, the tail statistic and the board-sized information column
    if a.floor_json is not None:
        fl = json.loads(a.floor_json.read_text())
        bmf = {}
        for f in FOLDS:
            r = fl["folds"][f]
            S = F_all[f]
            assert sorted(r["kept"] + r["below"] + r["unmeasured"]) == sorted(S.targets), \
                f"{f}: the floor record's targets differ from the dump's"
            kept = set(r["kept"])
            bmf[f] = np.array([t in kept for t in S.targets])
        res["floor"] = dict(json=str(a.floor_json), threshold_cpm=fl["threshold_cpm"],
                            kept={f: int(bmf[f].sum()) for f in FOLDS}, n={f: F_all[f].n for f in FOLDS})
        print(f"floor F {fl['threshold_cpm']:.4f} CPM; floor-kept targets " + " / ".join(
            f"{f} {int(bmf[f].sum())} of {F_all[f].n}" for f in FOLDS), flush=True)
        res["arms_floor"] = price_all(F_all, arms_all, inc_r, bmf, "floor ")[0]
    res["fold_cap_hi"] = {f: FD.pct_hi(F_all[f].kv, F_all[f].carried, FD.P_PROD) for f in FOLDS}
    print("fold cap hi (P 96.647 of key / median, carried): " + " / ".join(f"{f} {res['fold_cap_hi'][f]:.2f}" for f in FOLDS)
          + "  [A/B/C 5.2-7.0; the lab's L panels 7.8-19.9; 199's real 9.3-16.0]", flush=True)
    if a.class_json.exists():
        cls = json.loads(a.class_json.read_text())
        res["board_sized"] = {"source": cls.get("source"), "all": board_sized(res["arms"], cls)}
        if "arms_floor" in res:
            res["board_sized"]["floor"] = board_sized(res["arms_floor"], cls)
        for ver in ("all", "floor"):
            for name, row in res["board_sized"].get(ver, {}).items():
                print(f"board-sized (information only; {ver}) {name} [{row['cls']}]: {row['text']}", flush=True)
    res["seconds"] = round(time.time() - t0, 1)
    (a.out or HERE / "work" / f"score708-{a.tag}.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
