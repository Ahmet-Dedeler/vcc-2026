"""522 (4a, continued) — is the October-like panels' heavy key tail a property of how they were DRAWN? Re-draw them to
A/B/C's response-strength mix and re-price the key forms there (and the reverse: A/B/C drawn to the October-like mix).

Arc built the A/B/C panel from an H1 pilot binned by DE count (negligible / subtle / strong) and will build D/E/F the same
way; 438's October-like panels are uniform random draws from donor coverage. The strength matching variable is leak-free:
the mean robust strength (strength522) over the FOLD'S OWN two sources only (never the held-out line).
  matched   quintile bins of the reference panel's distribution of that variable (A/B/C's, same two sources); keep
            n = 5 x (smallest bin count) targets, n/5 per bin, sampled without replacement; 20 draws
  random    the same n drawn uniformly (20 draws): separates "matched mix" from "smaller panel"
  reverse   A/B/C folds matched to the POOLED October-like panels' distribution (efp + efr + real, same two sources)
  essential-like removed (October-like panels): every target with DepMap mean Chronos < -0.5 (17 lines) or on the K562
            essential list dropped (one deterministic sub-panel); and the same pool re-drawn to A/B/C's mix (20 draws):
            the "Arc-like" emulation (Arc's panels are essential-free and built from a DE-count-binned pilot)
Arms per draw (everything re-derived on the sub-panel: key median, the percentile, the strength ranks, key_scale, pds):
cap 3 (base), cap 6, the percentile cap, cap 8, uncapped, S1, pct x S0.25, pct x S0.5, cap 6 x S0.5.

    uv run python experiments/522-october-bulk-key/subpanel522.py [panel ...]     (~1-2 GB)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import fold522 as FD  # noqa: E402

WORK = HERE / "work"
NDRAW = 20
OCTP = ("efp", "efr", "real")


def src_strength(st, panel, targets, sources):
    pt = [str(t) for t in st[f"{panel}__targets"]]
    idx = np.asarray([pt.index(t) for t in targets])
    M = np.vstack([np.asarray(st[f"{panel}__{s}"], float)[idx] for s in sources])
    with np.errstate(all="ignore"):
        return np.nanmean(M, 0)


def zsub(st, S, mask):
    """zsrc on the sub-panel: crank over the kept targets only (0 elsewhere; unused there)."""
    pt = [str(t) for t in st[f"{S.panel}__targets"]]
    keep = np.flatnonzero(mask)
    idx = np.asarray([pt.index(S.targets[i]) for i in keep])
    cols = []
    for s in S.sources:
        v = np.asarray(st[f"{S.panel}__{s}"], float)[idx]
        r = FD.crank(v)
        r[~np.isfinite(v)] = np.nan
        cols.append(r)
    with np.errstate(all="ignore"):
        m = np.nanmean(np.vstack(cols), 0)
    z = np.zeros(S.n)
    z[keep] = np.where(np.isfinite(m), m, 0.0)
    return z


def arms(S, z, mask):
    kv, car = S.kv, S.carried
    w3 = FD.wclip(kv, car, 0.5, 3.0, mask)
    w6 = FD.wclip(kv, car, 0.5, 6.0, mask)
    hp = FD.pct_hi(kv, car, FD.P_PROD, mask)
    wp = FD.wclip(kv, car, 0.5, hp, mask)
    A = {"cap3": w3, "cap6": w6, "pctPROD": wp, "cap8": FD.wclip(kv, car, 0.5, 8.0, mask),
         "uncapped": FD.wclip(kv, car, 0.5, 1e9, mask), "S1": w3 * np.exp(z), "pctPRODxS0.25": wp * np.exp(0.25 * z),
         "pctPRODxS0.5": wp * np.exp(0.5 * z), "cap6xS0.5": w6 * np.exp(0.5 * z)}
    return A, hp


def draws(values, ref, n_draw, seed0, matched=True, pool=None):
    ok = np.isfinite(values) if pool is None else (np.isfinite(values) & pool)
    edges = np.nanquantile(ref[np.isfinite(ref)], [0.2, 0.4, 0.6, 0.8])
    b = np.digitize(values, edges)
    counts = [int(np.sum(ok & (b == k))) for k in range(5)]
    per = min(counts)
    n = 5 * per if matched != "all" else int(ok.sum())
    out = []
    for s in range(n_draw if matched != "all" else 1):
        rng = np.random.default_rng(seed0 + s)
        m = np.zeros(values.size, bool)
        if matched == "all":
            m = ok.copy()
        elif matched:
            for k in range(5):
                pool_k = np.flatnonzero(ok & (b == k))
                m[rng.choice(pool_k, per, replace=False)] = True
        else:
            m[rng.choice(np.flatnonzero(ok), n, replace=False)] = True
        out.append(m)
    return out, dict(bin_counts=counts, n_sub=n, edges=edges.tolist())


def run(panel, held, st):
    t0 = time.time()
    S = FD.load(panel, held, check=True)
    val = src_strength(st, panel, S.targets, S.sources)
    pt = [str(t) for t in st[f"{panel}__targets"]]
    ii = np.asarray([pt.index(t) for t in S.targets])
    essl = (np.asarray(st[f"{panel}__dep_mean"], float)[ii] < -0.5) | np.asarray(st[f"{panel}__k562_essential"])[ii]
    if panel == "abc":   # reverse: A/B/C drawn to the pooled October-like mix of the same two sources
        ref = np.concatenate([src_strength(st, p, [str(t) for t in st[f"{p}__targets"]], S.sources) for p in OCTP])
        modes = (("to October-like mix", True, None), ("random", False, None))
    else:
        ref = src_strength(st, "abc", [str(t) for t in st["abc__targets"]], S.sources)
        modes = (("to A/B/C mix", True, None), ("random", False, None), ("essential-like removed", "all", ~essl),
                 ("essential-like removed, to A/B/C mix", True, ~essl))
    rec = dict(panel=panel, fold=held, n=S.n, n_essential_like=int(essl.sum()), modes={})
    for label, matched, pool in modes:
        ms, info = draws(val, ref, NDRAW, 5220 if matched else 5320, matched, pool)
        res = {}
        his = []
        for m in ms:
            z = zsub(st, S, m)
            A, hp = arms(S, z, m)
            his.append(hp)
            base = None
            for nm, w in A.items():
                s, _ = S.score(w, mask=m)
                if base is None:
                    base = s
                res.setdefault(nm, dict(r=[], dpds=[]))
                res[nm]["r"].append(s["cos_agg"] / base["cos_agg"] - 1)
                res[nm]["dpds"].append(s["pds"] - base["pds"])
        summ = {nm: dict(r_mean=float(np.mean(v["r"])), r_sd=float(np.std(v["r"], ddof=1)) if len(v["r"]) > 1 else 0.0,
                         frac_pos=float(np.mean(np.asarray(v["r"]) > 0)), dpds_mean=float(np.mean(v["dpds"])),
                         r=v["r"]) for nm, v in res.items()}
        for nm in summ:   # paired: the form vs cap 6 on the same draws
            summ[nm]["vs_cap6_mean"] = float(np.mean(np.asarray(res[nm]["r"]) - np.asarray(res["cap6"]["r"])))
            summ[nm]["vs_cap6_frac_pos"] = float(np.mean(np.asarray(res[nm]["r"]) > np.asarray(res["cap6"]["r"])))
        rec["modes"][label] = dict(info=info, hi_mean=float(np.mean(his)), hi_q=np.quantile(his, [0, .5, 1]).tolist(),
                                   arms=summ)
        print(f"[{panel}/{held}] {label:20s} n_sub {info['n_sub']} (bins {info['bin_counts']}) hi {np.mean(his):.2f} "
              + " ".join(f"{nm} {100 * summ[nm]['r_mean']:+.2f}" for nm in summ if nm != "cap3")
              + f" ({time.time() - t0:.0f}s)", flush=True)
    full = src_strength(st, panel, S.targets, S.sources)
    rec["full_panel"] = dict(hi=FD.pct_hi(S.kv, S.carried, FD.P_PROD),
                             strength_quantiles=np.nanquantile(full, [.1, .25, .5, .75, .9]).tolist(),
                             ref_quantiles=np.nanquantile(ref, [.1, .25, .5, .75, .9]).tolist())
    return rec


def main():
    st = np.load(WORK / "strength522.npz")
    which = sys.argv[1:] or ("efp", "efr", "real", "abc")
    # one file per invocation (concurrent runs must not overwrite each other); analyse merges subpanel522*.json
    path = WORK / f"subpanel522-{'-'.join(which)}.json"
    out = json.loads(path.read_text()) if path.exists() else {}
    for p, f in FD.surfaces(which):
        if p == "abc" and f == "kolf2":
            continue
        if p != "abc" and not (WORK / f"dump522-{p}-{f}.npz").exists():
            continue
        out[f"{p}/{f}"] = run(p, f, st)
        path.write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
