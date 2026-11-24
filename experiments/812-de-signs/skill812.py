"""812 — skill812: a DE sign-skill scorer for ANY model's predictions on H1's 211 targets (h1-test 85 + h1-trainval 126).

WHY. The DE gap to the field is sign skill on responding genes (812 RESULTS §4). This scores a source's per-target signs
by the classes that pay, against CLEAN4 (the best public sign source) and against the incumbent's delivered sign, and
prices the difference in overall-equivalent units through rates812's per-class H1 rates. The model must never have seen
H1: VCC25 H1 (h1_2025 / h1_vcc2025) and Pan's H1 must be out of every bank it uses (the caller's responsibility; declare it
with --h1-free).

INPUT. One npz per surface: `targets` (str), `genes` (str), `pred` (T x G float; its sign is the prediction; 0 or NaN = no
prediction). Targets and genes are matched by name onto the surface's targets and 5-CPM gene axis (10,779); anything
missing counts as no prediction.

CLASSES (from the incumbent's realized call set on the surface: the production emitter's zs arm, seed 0; own gene out;
truth = the surface's reference DE table at BH 0.05):
  a0      called, truth non-significant (fid's filler; 88 % of calls)
  a1      called, truth significant
  a1head  called, truth significant, in the head (top-400 |key|): the front of the order (reach)
  u1      uncalled, truth significant (nmae's bulk; the hook's genes)
each also split by control CPM (5-20, 20-50, 50-100, 100-300, >= 300).

COMPARATORS on the cells the source signs:
  CLEAN4 (391's P_down sign; paired cells = both sign) and the incumbent's delivered sign (the hook arm's realized lfc:
  h1-test 612 h1_c4kgall s0, h1-trainval 618 hk_c4h s0; CLEAN4 tau 0.05 on uncalled genes, emitted signs on calls).

STATISTICS. Paired target bootstrap (4,000 draws, seed 812; a target's cells move together); ratio of sums. Pooled = the
equal-weight mean of the two surfaces, each resampled independently.

PRICE (checker units: overall-eq). If the source's sign replaced the comparator on the cells it signs, the class
accuracy moves by D = sum(right_src - right_cmp) / n_class, and overall-eq = D / 0.01 x rates812's "+0.01" price for that
class and surface (linear; rates812 measured +0.01 and +0.05 as near-linear). CAVEAT: the rates are H1's (one line, the
production emitter at f 1); board transfer differs by class (H1-BOARD-CALIBRATION: sign priors 0.07-0.67x on fid).

KNOWN ANSWERS (--ka): CLEAN4 through skill812 gives a0 0.5171 (h1-test; 700 signsplit row
`callset_head|block|all|b_nonsig` clean4 and a0class812 `a0_c4_on_c4cov`) and 0.5225 (h1-trainval, a0class812); the zs arm's
own realized signs give rates812's class accuracies (a0 0.5153 / 0.5174, a1 0.6005 / 0.6216, a1head 0.7222 / 0.7383,
u1 0.5366 / 0.5444) to 5e-5.

usage:
  uv run python experiments/812-de-signs/skill812.py --name NAME --pred h1test=A.npz --pred h1tv=B.npz [--h1-free yes]
  uv run python experiments/812-de-signs/skill812.py --ka
writes work/skill/<NAME>.json and prints the table.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import surf812 as SF  # noqa: E402

HOOK = {"h1test": SF.ROOT / "experiments/612-key-gated-hook/work/h1_c4kgall/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz",
        "h1tv": SF.ROOT / "experiments/618-h1-trainval-surface/work/hk_c4h/stats-base_s0_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"}
RATES = json.loads((SF.HERE / "work/rates812.json").read_text())
CLASSES = ("a0", "a1", "a1head", "u1")
BINS = ((5, 20), (20, 50), (50, 100), (100, 300), (300, 1e12))
B = 4000


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_pred(path, targets, genes):
    z = np.load(path, allow_pickle=False)
    pt, pg, P = z["targets"].astype(str), z["genes"].astype(str), np.asarray(z["pred"], np.float64)
    ti = {t: i for i, t in enumerate(pt)}
    gi = {g: i for i, g in enumerate(pg)}
    r = np.array([ti.get(t, -1) for t in targets])
    c = np.array([gi.get(g, -1) for g in genes])
    out = np.zeros((len(targets), len(genes)))
    sub = P[np.ix_(np.maximum(r, 0), np.maximum(c, 0))]
    sub[r < 0, :] = 0
    sub[:, c < 0] = 0
    out[:] = np.sign(np.nan_to_num(sub, nan=0.0))
    return out, int((r >= 0).sum()), int((c >= 0).sum())


def surface_state(surf):
    S = SF.load(surf, 0)
    adj, sig, call = SF.masks(S)
    h = np.load(HOOK[surf])
    assert list(h["targets"].astype(str)) == list(S["targets"]), "hook-arm targets differ"
    c4 = np.where(np.isfinite(S["c4"]), np.where(S["c4"] > 0.5, -1.0, 1.0), 0.0)
    cls = {"a0": call & adj & ~sig, "a1": call & sig & adj, "a1head": call & sig & adj & S["head"], "u1": ~call & sig & adj}
    return dict(S=S, cls=cls, c4=np.broadcast_to(c4[None, :], S["rl"].shape), inc=np.sign(np.nan_to_num(h["lfc"])),
                truth=np.sign(S["rl"]), cpm=S["cpm"])


def per_target(mask, src, cmp_, truth):
    """per-target counts on `mask`: n, right_src, right_cmp (cells where src signs; cmp may be 0 = wrong)."""
    n = mask.sum(1)
    rs = (mask & (src == truth)).sum(1)
    rc = (mask & (cmp_ == truth) & (cmp_ != 0)).sum(1)
    return n.astype(float), rs.astype(float), rc.astype(float)


def boot_idx(T, seed):
    return np.random.default_rng(seed).integers(0, T, size=(B, T))


def score_surface(st, src, seed):
    truth, c4, inc, cpm = st["truth"], st["c4"], st["inc"], st["cpm"]
    T = truth.shape[0]
    idx = boot_idx(T, seed)
    out = {}
    signed = src != 0
    for cn in CLASSES:
        cm = st["cls"][cn]
        n_class_t = cm.sum(1).astype(float)
        rows = {}
        for bn, gm in [("all", np.ones(cpm.size, bool))] + [(f"cpm{lo}-{hi:g}" if hi < 1e12 else f"cpm{lo}+", (cpm >= lo) & (cpm < hi)) for lo, hi in BINS]:
            m = cm & gm[None, :]
            mP = m & signed
            nc = float(m.sum())
            if nc == 0:
                continue
            r = dict(n_class=int(nc), coverage=round(float(mP.sum()) / nc, 4))
            # accuracy over the class (unsigned = wrong) and over the signed cells
            r["acc_all"] = round(float((m & (src == truth) & signed).sum()) / nc, 4)
            r["acc_signed"] = round(float((mP & (src == truth)).sum()) / max(float(mP.sum()), 1), 4) if mP.any() else None
            for cname, cmp_, paired in (("clean4", c4, mP & (c4 != 0)), ("incumbent", inc, mP)):
                n_t, rs_t, rc_t = per_target(paired, src, cmp_, truth)
                if n_t.sum() == 0:
                    continue
                d = (rs_t.sum() - rc_t.sum()) / n_t.sum()
                bd = (rs_t[idx].sum(1) - rc_t[idx].sum(1)) / np.maximum(n_t[idx].sum(1), 1)
                # price: class-accuracy change if the source replaced the comparator on these cells
                ncl_t = n_class_t
                D = (rs_t.sum() - rc_t.sum()) / ncl_t.sum()
                bD = (rs_t[idx].sum(1) - rc_t[idx].sum(1)) / np.maximum(ncl_t[idx].sum(1), 1)
                r[f"vs_{cname}"] = dict(n_paired=int(n_t.sum()), acc_src=round(rs_t.sum() / n_t.sum(), 4), acc_cmp=round(rc_t.sum() / n_t.sum(), 4),
                                        delta=round(d, 5), ci95=[round(float(np.quantile(bd, 0.025)), 5), round(float(np.quantile(bd, 0.975)), 5)],
                                        p_le0=round(float((bd <= 0).mean()), 4), dclass=D, boot_dclass=bD)
            rows[bn] = r
        out[cn] = rows
    return out


def price(res, surf):
    for cn in CLASSES:
        rate = RATES[surf]["rates"][f"{cn}+0.01"]["d_overall"]
        for bn, r in res[cn].items():
            for cname in ("clean4", "incumbent"):
                if f"vs_{cname}" not in r:
                    continue
                v = r[f"vs_{cname}"]
                D, bD = v.pop("dclass"), v.pop("boot_dclass")
                v["price_overall_eq"] = round(D / 0.01 * rate, 6)
                v["price_ci95"] = [round(float(np.quantile(bD, 0.025)) / 0.01 * rate, 6), round(float(np.quantile(bD, 0.975)) / 0.01 * rate, 6)]
                v["_boot_price"] = bD / 0.01 * rate
    return res


def run(name, preds, h1_free="undeclared", quiet=False):
    res = {"name": name, "h1_free": h1_free, "inputs": {}, "surfaces": {}}
    boots = {}
    for k, surf in enumerate(("h1test", "h1tv")):
        if surf not in preds:
            continue
        st = surface_state(surf)
        src, nt, ng = load_pred(preds[surf], st["S"]["targets"], st["S"]["genes"])
        res["inputs"][surf] = dict(path=str(preds[surf]), sha256=sha(preds[surf]), targets_matched=nt, genes_matched=ng)
        r = price(score_surface(st, src, 812 + k), surf)
        boots[surf] = {(cn, bn, c): r[cn][bn][f"vs_{c}"].pop("_boot_price") for cn in CLASSES for bn in r[cn]
                       for c in ("clean4", "incumbent") if f"vs_{c}" in r[cn][bn]}
        res["surfaces"][surf] = r
    if len(boots) == 2:   # pooled: equal-weight mean of the two surfaces' prices, independent resamples
        pooled = {}
        for key in set(boots["h1test"]) & set(boots["h1tv"]):
            cn, bn, c = key
            pt = 0.5 * (res["surfaces"]["h1test"][cn][bn][f"vs_{c}"]["price_overall_eq"] + res["surfaces"]["h1tv"][cn][bn][f"vs_{c}"]["price_overall_eq"])
            bp = 0.5 * (boots["h1test"][key] + boots["h1tv"][key])
            pooled.setdefault(cn, {}).setdefault(bn, {})[f"vs_{c}"] = dict(
                price_overall_eq=round(pt, 6), ci95=[round(float(np.quantile(bp, 0.025)), 6), round(float(np.quantile(bp, 0.975)), 6)],
                p_le0=round(float((bp <= 0).mean()), 4))
        res["pooled"] = pooled
        tot = {c: sum(pooled[cn]["all"][f"vs_{c}"]["price_overall_eq"] for cn in ("a0", "a1", "u1") if f"vs_{c}" in pooled[cn]["all"]) for c in ("clean4", "incumbent")}
        res["pooled_total_a0_a1_u1"] = {c: round(v, 6) for c, v in tot.items()}
    res["caveat"] = "rates are H1's (production emitter, f 1); board transfer differs by class (sign priors 0.07-0.67x on fid)"
    out = SF.HERE / "work/skill" / f"{name}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1))
    if not quiet:
        print(table(res))
    return res


def table(res):
    L = [f"skill812 [{res['name']}] h1_free={res['h1_free']}  (overall-eq at H1's per-class rates; {res['caveat']})",
         "| surface | class | coverage | acc (signed) | vs CLEAN4: acc src / C4 (paired n) | Δ [95 %] | price | vs incumbent: Δ [95 %] | price |",
         "|---|---|---|---|---|---|---|---|---|"]
    for surf, r in res["surfaces"].items():
        for cn in CLASSES:
            a = r[cn]["all"]
            c = a.get("vs_clean4", {})
            i = a.get("vs_incumbent", {})
            L.append(f"| {surf} | {cn} | {a['coverage']} | {a['acc_signed']} | {c.get('acc_src')} / {c.get('acc_cmp')} ({c.get('n_paired')}) | "
                     f"{c.get('delta')} {c.get('ci95')} | {c.get('price_overall_eq')} | {i.get('delta')} {i.get('ci95')} | {i.get('price_overall_eq')} |")
    if "pooled" in res:
        for cn in CLASSES:
            p = res["pooled"][cn]["all"]
            L.append(f"| pooled | {cn} | | | | | {p['vs_clean4']['price_overall_eq']} {p['vs_clean4']['ci95']} | | "
                     f"{p['vs_incumbent']['price_overall_eq']} {p['vs_incumbent']['ci95']} |")
        L.append(f"pooled total (a0 + a1 + u1): vs CLEAN4 {res['pooled_total_a0_a1_u1']['clean4']}, vs incumbent {res['pooled_total_a0_a1_u1']['incumbent']}")
    return "\n".join(L)


def ka():
    src = SF.HERE / "work/sources"
    r = run("ka_clean4", {s: src / f"clean4_{s}.npz" for s in ("h1test", "h1tv")}, "yes (CLEAN4 has no H1 bank)", quiet=True)
    a = {s: r["surfaces"][s]["a0"]["all"]["acc_signed"] for s in ("h1test", "h1tv")}
    assert a["h1test"] == 0.5171 and a["h1tv"] == 0.5225, f"KA FAIL clean4 a0 {a}"
    e = run("ka_emitted_zs", {s: src / f"emitted_zs_{s}.npz" for s in ("h1test", "h1tv")}, "n/a (our own emission)", quiet=True)
    for s in ("h1test", "h1tv"):
        for cn in CLASSES:
            got = e["surfaces"][s][cn]["all"]["acc_all"]
            want = RATES[s]["class_acc"][cn]["acc"]
            assert abs(got - want) < 5e-5, f"KA FAIL emitted {s} {cn}: {got} vs rates812 {want}"
    print(f"KA PASS: clean4 a0 {a}; emitted_zs class accuracies == rates812 (" +
          ", ".join(f"{s} {cn} {e['surfaces'][s][cn]['all']['acc_all']}" for s in ("h1test", "h1tv") for cn in CLASSES) + ")")
    print(table(r))
    print(table(e))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name")
    ap.add_argument("--pred", action="append", default=[], help="surface=path.npz (h1test, h1tv)")
    ap.add_argument("--h1-free", default="undeclared")
    ap.add_argument("--ka", action="store_true")
    a = ap.parse_args()
    if a.ka:
        ka()
    else:
        run(a.name, dict(p.split("=", 1) for p in a.pred), a.h1_free)
