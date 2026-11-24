"""513: anatomy of our reach walk on H1 (production emitter, 435's coarse-MDE stats), and what orderings could buy.

Reach (BRIEF section 5): pool = reference-significant genes of target p (target excluded); order = our significant calls
first, then ascending our adjusted p, then ascending our p, then descending |our lfc|, then gene name; walk counting only
adjudicable genes (reference lfc finite and non-zero); k* = deepest prefix with sign purity >= 0.9; R_p = k*/n_real.

Known answers first: fid, jac, nmae and reach recomputed from the saved stats must equal 435's metrics json.
    uv run python experiments/513-reach-walk-anatomy/walk513.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OURS = ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
MET = ROOT / "experiments/435-zs-on-coarse/work/h1/metrics-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.json"
REF = ROOT / "experiments/323-current-family-validation/work/reference.npz"
ALPHA, P0 = 0.05, 0.9


def load():
    o = np.load(OURS, allow_pickle=True)
    r = np.load(REF, allow_pickle=True)
    keep = r["keep"].astype(bool)
    genes = r["genes"].astype(str)[keep]
    assert list(o["targets"]) == list(r["targets"]), "target order differs"
    return dict(targets=r["targets"].astype(str), genes=genes, rq=np.asarray(r["q"], float),
                rl=np.asarray(r["lfc"], float), oq=np.asarray(o["q"], float), op=np.asarray(o["p"], float),
                ol=np.asarray(o["lfc"], float), ctrl=np.asarray(r["ref_lin"], float)[keep]
                if r["ref_lin"].shape[0] == keep.shape[0] else np.asarray(r["ref_lin"], float))


def order_keys(oq, op, ol, genes):
    """Our confidence order (BRIEF 5): sig calls first, adj p asc, p asc, |lfc| desc, name."""
    sig = (oq < ALPHA)
    return np.lexsort((genes, -np.abs(np.nan_to_num(ol)), np.nan_to_num(op, nan=1.0),
                       np.nan_to_num(oq, nan=1.0), ~sig))


def reach_row(pool_idx, order, rl, ol):
    """Returns (k*, n_real, walk array of +1 match / 0 mismatch for adjudicable pool genes in our order)."""
    n_real = len(pool_idx)
    if n_real == 0:
        return None
    inpool = np.zeros(len(order), bool)
    inpool[pool_idx] = True
    seq = [g for g in order if inpool[g]]
    walk = []
    for g in seq:
        t = rl[g]
        if not np.isfinite(t) or t == 0:
            continue                      # not adjudicable: skipped
        walk.append(1 if np.sign(ol[g]) == np.sign(t) and ol[g] != 0 else 0)
    walk = np.asarray(walk, int)
    if len(walk) == 0:
        return 0, n_real, walk
    cum = np.cumsum(walk) / np.arange(1, len(walk) + 1)
    ok = np.nonzero(cum >= P0)[0]
    kstar = int(ok[-1] + 1) if len(ok) else 0
    return kstar, n_real, walk


def main():
    d = load()
    T, G = d["targets"], d["genes"]
    res = {"per_target": []}
    fid, jac, reach, nm = [], [], [], []
    for k, t in enumerate(T):
        own = G == t
        rs = (d["rq"][k] < ALPHA) & ~own
        adj = np.isfinite(d["rl"][k]) & (d["rl"][k] != 0) & ~own
        os_ = (d["oq"][k] < ALPHA) & ~own
        # fid (BRIEF 4)
        n_real = int(rs.sum()); calls = os_ & adj; n_pred = int(calls.sum())
        kk = int((np.sign(d["ol"][k][calls]) == np.sign(d["rl"][k][calls])).sum())
        if n_real > 0 or n_pred > 0:
            fid.append(kk / max(n_pred, n_real) if max(n_pred, n_real) > 0 else np.nan)
        # jac (BRIEF 6)
        u = (rs | os_).sum(); jac.append((rs & os_).sum() / u if u else 1.0)
        # nmae (BRIEF 7): gate = ref sig, finite lfc, non-own; >= 10
        gate = rs & np.isfinite(d["rl"][k])
        if gate.sum() >= 10:
            lh = np.nan_to_num(d["ol"][k][gate], nan=0.0, posinf=0.0, neginf=0.0)
            nm.append(np.abs(lh - d["rl"][k][gate]).sum() / np.abs(d["rl"][k][gate]).sum())
        # reach
        order = order_keys(d["oq"][k], d["op"][k], d["ol"][k], G)
        rr = reach_row(np.nonzero(rs)[0], order, d["rl"][k], d["ol"][k])
        if rr is not None:
            ks, nr, walk = rr
            reach.append(ks / nr)
            first_miss = int(np.argmin(walk)) if (walk == 0).any() else -1
            res["per_target"].append(dict(target=t, n_real=nr, kstar=ks, reach=ks / nr, n_pred=n_pred,
                                          walk_len=len(walk), head10_acc=float(walk[:10].mean()) if len(walk) else None,
                                          first_miss=first_miss, walk_acc=float(walk.mean()) if len(walk) else None))
    got = dict(fid=float(np.nanmean(fid)), jac=float(np.mean(jac)), nmae=float(np.mean(nm)), reach=float(np.mean(reach)))
    want = json.load(open(MET))["raw"]
    res["known_answers"] = {m: dict(got=got[m], want=want[m], ok=abs(got[m] - want[m]) < 1e-9) for m in got}
    print(json.dumps(res["known_answers"], indent=1))
    pt = res["per_target"]
    nr = np.array([r["n_real"] for r in pt]); ks = np.array([r["kstar"] for r in pt]); rc = np.array([r["reach"] for r in pt])
    wa = np.array([r["walk_acc"] for r in pt], float); h10 = np.array([r["head10_acc"] for r in pt], float)
    fm = np.array([r["first_miss"] for r in pt])
    print(f"targets {len(pt)}; n_real median {np.median(nr):.0f} (IQR {np.percentile(nr,25):.0f}-{np.percentile(nr,75):.0f}); "
          f"k* median {np.median(ks):.0f}; reach mean {rc.mean():.4f}")
    print(f"walk sign accuracy (whole pool, our signs) median {np.nanmedian(wa):.3f} mean {np.nanmean(wa):.3f}; "
          f"head-10 accuracy mean {np.nanmean(h10):.3f}; first miss position median {np.median(fm[fm>=0]):.0f}")
    for lo, hi in [(1, 20), (20, 100), (100, 1000), (1000, 20000)]:
        m = (nr >= lo) & (nr < hi)
        if m.any():
            print(f"  n_real [{lo},{hi}): {m.sum():3d} targets, reach {rc[m].mean():.3f}, walk acc {np.nanmean(wa[m]):.3f}, "
                  f"head10 {np.nanmean(h10[m]):.3f}")
    (HERE / "work").mkdir(exist_ok=True)
    json.dump(res, open(HERE / "work/walk513.json", "w"), indent=1)


if __name__ == "__main__":
    main()
