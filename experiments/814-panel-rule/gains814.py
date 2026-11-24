"""814 -- which targets carry the key / direction gains on the L folds, and does A/B/C's panel lack them?

Reads only arrays already on disk (the standard-format npz that panelrank800 reads; MBs) and small target lists.
For every (arm pair, held line, fold) it computes the fold overall (stats800's bulk_fold_raw -> overall, pds + mse) and the
exact leave-one-target-out contribution of every target, then reports
  (a) the share of each fold's gain carried by a flagged target set (default: targets of the RPE1 / Jurkat essential-gene
      screens, data/donors/{rpe1,nadig_jurkat}/perts.txt, which A/B/C's panel has none of),
  (b) panelrank800's statistic (ABCp's rank among L1-L6) after the flagged targets are removed from the L folds,
  (c) truth-strength profiles per fold: quantiles of sqrt(d) = ||delta_truth|| on the held line, top-30 share of sum d.
KNOWN ANSWER (abort): the all-target fold values equal panelrank800.json's stored values (|diff| < 1e-12).
    uv run python experiments/814-panel-rule/gains814.py [--flag rpe1jurkat|none|file:<list>] [--keep-frac p]
ABCp is never cut (a v2 rule changes the lab's draws, never Arc's panel); --keep-frac p keeps a random share p of the
flagged L targets (20 fixed seeds 814100+, fold values averaged) to emulate a MATCHED share rather than an exclusion.
Light: < 50 MB, ~1 s.
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import argparse  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/800-checker"))
import panelrank800 as PR  # noqa: E402
import stats800 as S  # noqa: E402

HELD = PR.HELD


def lines(p):
    return {ln.strip() for ln in Path(p).read_text().splitlines() if ln.strip()}


def flagset(name: str) -> set[str]:
    if name == "none":
        return set()
    if name == "rpe1jurkat":
        return lines(ROOT / "data/donors/rpe1/perts.txt") | lines(ROOT / "data/donors/nadig_jurkat/perts.txt")
    if name.startswith("file:"):
        return lines(name[5:])
    raise SystemExit(f"unknown flag {name}")


def fold_overall(zb, zc, idx):
    raw = S.bulk_fold_raw(zb, zc, idx)
    return S.overall_of({m: raw[m] for m in ("pds", "mse")})


def loo(zb, zc):
    """exact leave-one-target-out fold overall, vectorised over targets (frame-averaged X, as stats800._bulk_X)."""
    n = zb["d"].size
    D = zb["d"].sum()

    def X_minus(z):
        x, a, d = z["x"], z["a"], z["d"]
        return np.mean((x.sum(1)[:, None] - x) / np.sqrt((a.sum(1)[:, None] - a) * (D - d)[None, :]), 0)

    r = X_minus(zc) / X_minus(zb) - 1.0
    dp = ((zc["pds"] - zb["pds"]).sum(1)[:, None] - (zc["pds"] - zb["pds"])).mean(0) / (n - 1)
    return (S.SIGN["pds"] * dp * S.INV["pds"] + S.SIGN["mse"] * S.MSE_PER_PCT_X * 100.0 * r * S.INV["mse"]) / 6.0


REF = "ABCp"      # the A/B/C-panel comparator: ABCp (A/B/C fold rules) or LABC (801 follow-up 4: the L-fold recipe)


def load_folds():
    """{pair: {held: {fold: (zb, zc)}}} with panelrank800's root order (first root that has the fold wins); the comparator
    fold is stored under the key "ABCp" whichever REF is chosen."""
    out = {}
    for name, (cand, base, roots, _) in PR.PAIRS.items():
        out[name] = {}
        for h in HELD:
            v = {}
            for r in list(roots) + (PR.ROOTS_LABC if REF == "LABC" else []):
                for d in sorted(r.glob(f"*__{h}")):
                    b, c = d / f"{base}.npz", d / f"{cand}.npz"
                    key = d.name.split("__")[0]
                    if REF == "LABC" and key == "ABCp":
                        continue
                    key = "ABCp" if key == REF else key
                    if b.exists() and c.exists() and key not in v and (key == "ABCp" or (key[:1] == "L" and key[1:].isdigit())):
                        v[key] = (dict(np.load(b, allow_pickle=True)), dict(np.load(c, allow_pickle=True)))
            out[name][h] = v
    return out


def rankstat(vals: dict):
    L = np.array([vals[k] for k in sorted(vals) if k != "ABCp"])
    rank = int(1 + np.sum(L < vals["ABCp"]))
    return rank, len(L) + 1


REPEATS = 20


def read_classes(p: str) -> dict[str, str]:
    out = {}
    for ln in Path(p).read_text().splitlines():
        if ln.strip():
            g, c = ln.split("\t")
            out[g] = c
    return out


def main(flag: str, keep_frac: float = 0.0, excl: str = "none", control: str = "none", match: str = "none"):
    F = flagset(flag)
    E = flagset(excl)
    CL = read_classes(match) if match != "none" else None
    abc_share = None
    if CL is not None:                       # A/B/C's own class shares (its 300 targets), the template every L fold is thinned to
        import pandas as pd
        abc = [str(t) for t in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"]]
        cnt = {}
        for t in abc:
            cnt[CL.get(t, "NA")] = cnt.get(CL.get(t, "NA"), 0) + 1
        abc_share = {c: v / len(abc) for c, v in cnt.items()}
    folds = load_folds()
    stored = json.loads((ROOT / "experiments/800-checker/work/panelrank800.json").read_text())
    rep, ranks_all, ranks_cut = {}, [], []
    ka_max = 0.0
    for name, byh in folds.items():
        rep[name] = {}
        for h, fs in byh.items():
            if "ABCp" not in fs or len(fs) < 4:
                continue
            row = {}
            v_all, v_cut = {}, {}
            for k, (zb, zc) in sorted(fs.items()):
                t = [str(x) for x in zb["targets"]]
                n = len(t)
                v_all[k] = fold_overall(zb, zc, None)
                if k != "ABCp" and k in stored[name][h]["L"]:
                    ka_max = max(ka_max, abs(v_all[k] - stored[name][h]["L"][k]))
                if k == "ABCp" and REF == "ABCp":
                    ka_max = max(ka_max, abs(v_all[k] - stored[name][h]["ABCp"]["value"]))
                fl = np.array([x in F for x in t])
                ex = np.array([x in E for x in t])
                fl = fl & ~ex
                if CL is not None and k != "ABCp":
                    # MATCH: thin the fold (after the full exclusion) so its class mix equals A/B/C's: keep prob per class
                    # proportional to abc_share / fold_share, scaled to max 1; 20 seeded repeats, fold values averaged
                    base_keep = ~ex
                    cls_t = np.array([CL.get(x, "NA") for x in t])
                    fs_ = {c: np.mean(cls_t[base_keep] == c) for c in set(cls_t[base_keep])}
                    ratio = {c: abc_share.get(c, 0.0) / fs_[c] for c in fs_}
                    mx = max(ratio.values())
                    pk = np.array([ratio.get(c, 0.0) / mx for c in cls_t])
                    vs = []
                    for rr in range(REPEATS):
                        g = np.random.default_rng(814300 + rr)
                        kp = base_keep & (g.random(n) < pk)
                        vs.append(fold_overall(zb, zc, np.where(kp)[0]))
                    v_cut[k] = float(np.mean(vs))
                    row_keep = float(np.mean(pk[base_keep]))
                elif control != "none" and k != "ABCp":
                    # CONTROL: remove as many targets as the flag holds in this fold, chosen by truth strength d (top) or
                    # at random (20 repeats), to tell the flag's identity from its count / its strength
                    m = int(fl.sum())
                    if control == "topd":
                        drop = np.argsort(zb["d"])[::-1][:m]
                        v_cut[k] = fold_overall(zb, zc, np.setdiff1d(np.arange(n), drop))
                    elif control == "topa":
                        drop = np.argsort(zb["a"].mean(0))[::-1][:m]
                        v_cut[k] = fold_overall(zb, zc, np.setdiff1d(np.arange(n), drop))
                    else:
                        vs = [fold_overall(zb, zc, np.setdiff1d(np.arange(n), np.random.default_rng(814200 + rr).choice(n, m, replace=False)))
                              for rr in range(REPEATS)]
                        v_cut[k] = float(np.mean(vs))
                elif keep_frac > 0 and k != "ABCp" and fl.any():       # matched share: keep a random keep_frac of the flagged
                    vs = []
                    for rr in range(REPEATS):
                        g = np.random.default_rng(814100 + rr)
                        f2 = ex | (fl & (g.random(fl.size) >= keep_frac))
                        vs.append(fold_overall(zb, zc, np.where(~f2)[0]))
                    v_cut[k] = float(np.mean(vs))
                elif k == "ABCp":
                    v_cut[k] = v_all[k]                               # Arc's panel is never cut
                else:
                    v_cut[k] = fold_overall(zb, zc, np.where(~(fl | ex))[0]) if (fl | ex).any() else v_all[k]
                row[k] = dict(n=n, n_flag=int(fl.sum()), v_all=v_all[k], v_cut=v_cut[k],
                              v_flag_effect=v_all[k] - v_cut[k])
            ra, na = rankstat(v_all)
            rc, nc = rankstat(v_cut)
            ranks_all.append((ra, na))
            ranks_cut.append((rc, nc))
            rep[name][h] = dict(folds=row, rank_all=f"{ra}/{na}", rank_cut=f"{rc}/{nc}")
            Ls = [row[k] for k in sorted(row) if k != "ABCp"]
            print(f"{name:30s} {h:14s} L mean {np.mean([x['v_all'] for x in Ls]):+.5f} -> cut {np.mean([x['v_cut'] for x in Ls]):+.5f} "
                  f"(n_flag {np.mean([x['n_flag'] for x in Ls]):.1f}) | ABCp {row['ABCp']['v_all']:+.5f} | rank {ra}/{na} -> {rc}/{nc}")
    assert ka_max < 1e-12, f"KA FAIL: fold values differ from panelrank800.json by {ka_max}"
    for tag, rr in (("all targets", ranks_all), (f"L folds without flag '{flag}'", ranks_cut)):
        r = np.array([a for a, _ in rr]); n = np.array([b for _, b in rr])
        print(f"ABCp vs L ({tag}): mean percentile {np.mean((r - 0.5) / n):.3f}; bottom {int(np.sum(r == 1))}/{len(r)}")
    rep["_summary"] = {tag: dict(mean_percentile=float(np.mean((np.array([a for a, _ in rr]) - 0.5) / np.array([b for _, b in rr]))),
                                 bottom=int(sum(a == 1 for a, _ in rr)), cells=len(rr))
                       for tag, rr in (("all", ranks_all), ("cut", ranks_cut))}
    rep["_ka_max_abs_diff_vs_panelrank800"] = ka_max
    rep["_flag"] = dict(name=flag, size=len(F))
    # truth-strength profile per fold (any arm's base has d; use the 495 pair's base arrays)
    prof = {}
    for h, fs in folds["495 cap6-cap3 (key)"].items():
        prof[h] = {}
        for k, (zb, _) in sorted(fs.items()):
            d = zb["d"]; s = np.sqrt(d)
            t = [str(x) for x in zb["targets"]]
            fl = np.array([x in F for x in t])
            prof[h][k] = dict(n=int(d.size), q50=float(np.median(s)), q90=float(np.quantile(s, .9)), q99=float(np.quantile(s, .99)),
                              top30_share=float(np.sort(d)[::-1][:30].sum() / d.sum()),
                              flag_q50=float(np.median(s[fl])) if fl.any() else None,
                              a_q90=float(np.quantile(np.sqrt(zb["a"].mean(0)), .9)),
                              a_top30_share=float(np.sort(zb["a"].mean(0))[::-1][:30].sum() / zb["a"].mean(0).sum()))
            p = prof[h][k]
            print(f"  profile {h:14s} {k:5s} n {p['n']:3d} |dT| q50 {p['q50']:.3f} q90 {p['q90']:.3f} q99 {p['q99']:.3f} "
                  f"top30 {p['top30_share']:.3f} | flag q50 {p['flag_q50'] if p['flag_q50'] is None else round(p['flag_q50'], 3)} "
                  f"| pred |P| q90 {p['a_q90']:.3f} top30 {p['a_top30_share']:.3f}")
    rep["_profile"] = prof
    tag = Path(flag[5:]).stem if flag.startswith("file:") else flag
    etag = "" if excl == "none" else "-excl" + (Path(excl[5:]).stem if excl.startswith("file:") else excl)
    ctag = ("" if control == "none" else f"-control{control}") + ("" if match == "none" else f"-match{Path(match).stem}")
    rep["_flag"]["match"] = match
    out = HERE / "work" / f"gains814-{tag}{f'-keep{keep_frac:g}' if keep_frac else ''}{etag}{ctag}.json"
    rep["_flag"]["control"] = control
    rep["_flag"]["ref"] = REF
    if REF != "ABCp":
        out = out.with_name(out.stem + f"-ref{REF}.json")
    rep["_flag"]["excl"] = excl
    rep["_flag"]["keep_frac"] = keep_frac
    out.write_text(json.dumps(rep, indent=1, default=float))
    print(f"KA PASS (max |diff| {ka_max:.2e}) -> {out.relative_to(ROOT)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--flag", default="rpe1jurkat")
    ap.add_argument("--excl", default="none", help="a second flag set, always removed in full from the L folds")
    ap.add_argument("--ref", default="ABCp", choices=("ABCp", "LABC"), help="the A/B/C-panel comparator folds")
    ap.add_argument("--match", default="none", help="gene<TAB>class file: thin L folds to A/B/C's class mix (20 repeats)")
    ap.add_argument("--control", default="none", choices=("none", "topd", "topa", "random"),
                    help="remove the flag's per-fold COUNT of targets by top truth d / top predicted a / at random instead")
    ap.add_argument("--keep-frac", type=float, default=0.0, help="keep this random share of flagged L targets (20 repeats)")
    a = ap.parse_args()
    REF = a.ref
    main(a.flag, a.keep_frac, a.excl, a.control, a.match)
