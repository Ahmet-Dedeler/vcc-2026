"""518 — per-target features the LODO key does not use, per PLOLO fold, and how well each predicts the residual ideal scale.

Response (from arms540's saved per-target arrays of the base `V1_0.5_3@R16` = 446bj2's bulk on the fold):
    u_p = (sum_f x_pf / sum_f a_pf) / (sum x / sum a)   the magnitude oracle's multiplier on target p (can be < 0).
Features (operator-free: raw donor delta rows via size410.donor_rows, controls via plolo._b_ctrl_on_axis):
  per source s of the fold (the lines the fold's prediction is built from):
    kd_s    own-gene delta (the target's own knockdown, log units; more negative = deeper)
    ex_s    the target's control expression (log1p CPM-like, TS 5e4)
    nrm_s   non-panel row norm over genes > 10 CPM in s
    de_s    number of genes > 10 CPM in s with |delta| > 0.25 (a DE-count proxy)
    pan_s   share of the row's energy on OTHER panel genes' columns (the mse scores them; the bulk zeroes them)
  cross-source: mean pairwise cosine of the rows (non-panel, genes > 10 CPM in both); sd of log norms; sd of kd.
  context: the target's expression in the challenge contexts A/B/C (what a production per-context key could read).
  held line (diagnostic ONLY; unavailable for A/B/C): kd_held, ex_held, the held row's panel share.
Screen: Spearman(u, f) per fold, and the analytic value of a tilt  w_p -> w_p (f_p / median f)^beta  (first order on the
saved x, a: cos_agg(s) = mean_f sum s x / sqrt(sum s^2 a sum d); headroom540's validated route), beta in [-1, 1].

    uv run python experiments/518-bulk-mse-gap/feat518.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit"))
import plolo  # noqa: E402

plolo.WORK = ROOT / "experiments/420-bulk-stack/work/plolo"      # read-only: challenge b_ctrl caches
from plolo import E062, DONOR_PATHS, FOLDS, OFFICIAL, PANEL, TS, _b_ctrl_on_axis  # noqa: E402
from size410 import donor_rows  # noqa: E402

W540 = ROOT / "experiments/540-lodo-clip-scan/work"
FOLD_NAMES = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLD_NAMES[:3]
BASE = "V1_0.5_3@R16"
OFF = {g: i for i, g in enumerate(OFFICIAL)}
PANEL_POS = np.asarray([OFF[g] for g in PANEL])
BETAS = (-1.0, -0.75, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, 1.0)


def spearman(u, v):
    ok = np.isfinite(u) & np.isfinite(v)
    ru, rv = np.argsort(np.argsort(u[ok])), np.argsort(np.argsort(v[ok]))
    return float(np.corrcoef(ru, rv)[0, 1])


def cosagg(x, a, d, s):
    return float(np.mean([(s * x[j]).sum() / np.sqrt((s * s * a[j]).sum() * d.sum()) for j in range(3)]))


_CTRL = {}


def ctrl(name):
    if name not in _CTRL:
        b, _ = _b_ctrl_on_axis(E062.Bank.load(name, DONOR_PATHS[name]))
        cpm = 1e6 * np.expm1(b) / np.expm1(b).sum()
        _CTRL[name] = (b, cpm > 10.0)
    return _CTRL[name]


def line_features(name, targets):
    R, have = donor_rows(name, targets)
    b, expd = ctrl(name)
    own = np.asarray([OFF[t] for t in targets])
    nonpanel = np.ones(len(OFFICIAL), bool)
    nonpanel[PANEL_POS] = False
    gate = expd & nonpanel
    kd = R[np.arange(len(targets)), own]
    ex = b[own]
    nrm = np.linalg.norm(R[:, gate], axis=1)
    de = (np.abs(R[:, gate]) > 0.25).sum(1).astype(float)
    pmask = np.zeros(len(OFFICIAL), bool)
    pmask[PANEL_POS] = True
    tot = (R * R).sum(1) - kd ** 2
    pan = ((R[:, pmask] ** 2).sum(1) - kd ** 2) / np.maximum(tot, 1e-300)
    for arr in (kd, ex, nrm, de, pan):
        arr[~have] = np.nan
    return dict(R=R, have=have, gate=gate, kd=kd, ex=ex, nrm=nrm, de=de, pan=pan)


def main():
    out = dict(folds={})
    chal = {c: plolo.challenge_bctrl(c) for c in "ABC"}
    allrows = {}
    for f in FOLD_NAMES:
        z = np.load(W540 / f"arrays540-{f}.npz")
        x, a, d = z[f"{BASE}__x"], z[f"{BASE}__a"], z["d"]
        dump = np.load(WORK / f"dump518-{f}.npz")
        targets = [str(t) for t in dump["targets"]]
        sources = [str(s) for s in dump["sources"]]
        assert x.shape[1] == len(targets)
        u = (x.sum(0) / a.sum(0)) / (x.sum() / a.sum())
        kv1, carried = dump["kv1"], dump["carried"]
        med = np.median(kv1[carried])
        w = np.where(carried, np.clip(kv1 / med, 0.5, 3.0), 1.0)
        feats = {}
        src = {s: line_features(s, targets) for s in sources}
        held = line_features(f, targets)
        for s in sources:
            for k in ("kd", "ex", "nrm", "de", "pan"):
                feats[f"{k}_{s}"] = src[s][k]
        stack = lambda k: np.vstack([src[s][k] for s in sources])
        with np.errstate(all="ignore"):
            feats["kd_mean"] = np.nanmean(stack("kd"), 0)
            feats["kd_min"] = np.nanmin(stack("kd"), 0)
            feats["kd_sd"] = np.nanstd(stack("kd"), 0)
            feats["ex_mean"] = np.nanmean(stack("ex"), 0)
            feats["lognrm_mean"] = np.nanmean(np.log(stack("nrm")), 0)
            feats["lognrm_sd"] = np.nanstd(np.log(stack("nrm")), 0)
            feats["logde_mean"] = np.nanmean(np.log1p(stack("de")), 0)
            feats["logde_sd"] = np.nanstd(np.log1p(stack("de")), 0)
            feats["pan_mean"] = np.nanmean(stack("pan"), 0)
        # cross-source agreement: mean pairwise cosine over genes expressed in both, non-panel
        cs = []
        for i in range(len(sources)):
            for j in range(i + 1, len(sources)):
                si, sj = src[sources[i]], src[sources[j]]
                g = si["gate"] & sj["gate"]
                num = (si["R"][:, g] * sj["R"][:, g]).sum(1)
                den = np.linalg.norm(si["R"][:, g], axis=1) * np.linalg.norm(sj["R"][:, g], axis=1)
                c = num / np.maximum(den, 1e-300)
                c[~(si["have"] & sj["have"])] = np.nan
                cs.append(c)
        with np.errstate(all="ignore"):
            feats["src_cos"] = np.nanmean(np.vstack(cs), 0)
        own = np.asarray([OFF[t] for t in targets])
        for c in "ABC":
            feats[f"ex_ctx{c}"] = chal[c][own]
        feats["ex_ctx_mean"] = np.mean([chal[c][own] for c in "ABC"], 0)
        feats["ex_ctx_minus_src"] = feats["ex_ctx_mean"] - feats["ex_mean"]
        feats["key_w"] = w
        feats["log_kv1"] = np.log(np.maximum(kv1, 1e-6))
        feats["our_energy"] = a.mean(0)
        feats["pre_key_energy"] = a.mean(0) / w ** 2
        # held-line diagnostics (not usable in production)
        feats["HELD_kd"] = held["kd"]
        feats["HELD_ex"] = held["ex"]
        feats["HELD_ex_minus_src"] = held["ex"] - feats["ex_mean"]
        feats["HELD_kd_minus_src"] = held["kd"] - feats["kd_mean"]
        feats["HELD_pan"] = held["pan"]
        rec = dict(n=len(targets), sources=sources, spearman={}, tilt={})
        base_cos = cosagg(x, a, d, np.ones(len(targets)))
        for k, v in feats.items():
            rec["spearman"][k] = spearman(u, v)
            vv = np.asarray(v, float)
            ok = np.isfinite(vv) & (vv > 0) if not k.startswith(("kd", "HELD_kd", "log", "ex_ctx_minus", "HELD_ex_minus", "src_cos", "pan", "HELD_pan")) else np.isfinite(vv)
            if ok.sum() < 0.8 * len(vv):
                continue
            # tilt: positive features use (f/median)^beta; signed ones use exp(beta * (f - median)/sd)
            med_f = np.nanmedian(vv[ok])
            sd_f = np.nanstd(vv[ok]) or 1.0
            positive = np.all(vv[ok] > 0) and not k.startswith(("kd", "HELD_kd", "log", "ex_", "HELD_ex", "src_cos", "pan", "HELD_pan"))
            tl = {}
            for be in BETAS:
                if positive:
                    m = np.where(ok, (np.where(ok, vv, med_f) / med_f) ** be, 1.0)
                else:
                    m = np.where(ok, np.exp(be * (np.where(ok, vv, med_f) - med_f) / sd_f), 1.0)
                tl[f"{be:+g}"] = cosagg(x, a, d, m) / base_cos - 1
            rec["tilt"][k] = tl
        out["folds"][f] = rec
        allrows[f] = dict(u=u, feats=feats)
        print(f"== {f} (sources {sources}); base cos {base_cos:.5f}")
    # summary: features ranked by mean cancer Spearman and best consistent tilt
    keys = sorted(out["folds"]["k562_gwps"]["spearman"])
    summ = []
    for k in keys:
        sp = {f: out["folds"][f]["spearman"].get(k, np.nan) for f in FOLD_NAMES}
        best = None
        if all(k in out["folds"][f]["tilt"] for f in FOLD_NAMES):
            for be in BETAS:
                r = {f: out["folds"][f]["tilt"][k][f"{be:+g}"] for f in FOLD_NAMES}
                cm = float(np.mean([r[f] for f in CANCER]))
                if best is None or cm > best[1]:
                    best = (be, cm, r)
        summ.append(dict(feature=k, spearman=sp, best=best))
    summ.sort(key=lambda r: -abs(np.nanmean([r["spearman"][f] for f in CANCER])))
    print(f"{'feature':24s} " + " ".join(f"{f[:7]:>8s}" for f in FOLD_NAMES) + "   best tilt beta: cancer mean (per fold), kolf2")
    for r in summ:
        b = r["best"]
        bt = (f"b {b[0]:+g}: {100 * b[1]:+.2f} % (" + " ".join(f"{100 * b[2][f]:+.2f}" for f in CANCER) + f") kolf2 {100 * b[2]['kolf2']:+.2f}") if b else ""
        print(f"{r['feature']:24s} " + " ".join(f"{r['spearman'][f]:+8.3f}" for f in FOLD_NAMES) + "   " + bt)
    out["summary"] = summ
    (WORK / "feat518.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez_compressed(WORK / "feat518_arrays.npz", **{f"{f}__{k}": np.asarray(v, float) for f, r in allrows.items()
                                                         for k, v in [("u", r["u"])] + list(r["feats"].items())})


if __name__ == "__main__":
    main()
