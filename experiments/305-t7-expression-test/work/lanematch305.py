"""305 step 6 -- POST-HOC technical control (declared as a deviation added AFTER the outcomes):
LANE-MATCHED controls.

The frozen test references a construct to the whole non-targeting pool.  If a construct's cells sit in
lanes whose control composition differs, a purely technical lane effect can look like repression.  This
rebuilds every effect against a control pseudobulk carrying the *construct's own lane composition*:

    S_ctrl_matched(c) = sum_lane ( n_{c,lane} / N_{ctrl,lane} ) * S_{ctrl,lane}

and recomputes e_pool, both halves and the empirical rank (every null construct also lane-matched to
itself) for the 15 candidate pairs.  The 31 non-targeting constructs are byte-identical across the three
2025 files, so one stream of the train file supplies the per-lane control sums for every split.

    PYTHONDONTWRITEBYTECODE=1 uv run python lanematch305.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib305 as X  # noqa: E402
from expr305 import MAG, P_RAW, Split  # noqa: E402


def ctrl_lane_sums(P):
    """Per (lane, NTC parity) exact count sums, from ONE stream of the train file."""
    cache = X.HERE / "ctrl_lane_sums.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        return z["S"], z["lanes"].astype(str)
    path = X.SPLITS["train"]
    gid, tgt, bat, genes, shape = X.read_obs(path)
    cons = np.asarray([X.construct_of(g) for g in gid])
    is_ctrl = tgt == X.CTRL
    ntc = sorted(set(cons[is_ctrl].tolist()))
    par = {c: i % 2 for i, c in enumerate(ntc)}
    lanes = sorted(set(bat.tolist()))
    lidx = {b: i for i, b in enumerate(lanes)}
    groups = np.full(shape[0], -1, dtype=np.int64)
    sel = np.flatnonzero(is_ctrl)
    groups[sel] = np.asarray([2 * lidx[bat[i]] + par[cons[i]] for i in sel], dtype=np.int64)
    P(f"  streaming per-lane control sums: {len(lanes)} lanes x 2 parities = {2 * len(lanes)} groups")
    acc, cnt, _ = X.stream_group_sums(path, groups, 2 * len(lanes), log=P)
    S = np.rint(acc).astype(np.int64)
    np.savez_compressed(cache, S=S, lanes=np.asarray(lanes), counts=cnt)
    return S, np.asarray(lanes)


def matched(sp: Split, Sl: np.ndarray, lanes: np.ndarray, members: list[str], parities: tuple[int, ...]):
    """Lane-matched control sums for each construct in `members`, using the given control parities."""
    li = {b: i for i, b in enumerate(lanes)}
    rows = [Sl[2 * i + p] for i in range(len(lanes)) for p in parities]
    Sc = np.zeros((len(lanes), Sl.shape[1]), dtype=np.float64)
    Nc = np.zeros(len(lanes))
    for i in range(len(lanes)):
        for p in parities:
            Sc[i] += Sl[2 * i + p]
            Nc[i] += Sl[2 * i + p].sum()
    W = np.zeros((len(members), len(lanes)))
    cov = np.zeros(len(members))
    for k, c in enumerate(members):
        comp = {}
        for h in (parities if False else (0, 1)):
            for b, n in sp.lanes[f"{c}|{h}"].items():
                comp[b] = comp.get(b, 0) + n
        tot = sum(comp.values())
        for b, n in comp.items():
            if Nc[li[b]] > 0:
                W[k, li[b]] = n / Nc[li[b]]
                cov[k] += n
        cov[k] /= max(tot, 1)
    return W @ Sc, cov


def main() -> None:
    P = X.L.Logger(X.HERE / "lanematch305.log")
    Sl, lanes = ctrl_lane_sums(P)
    SP = {s: Split(s) for s in X.SPLITS}
    for s, sp in SP.items():
        assert set(lanes.tolist()) == set(sp.meta["lane_composition"][f"{sp.constructs[0]}|0"].keys()) | \
            {b for c in sp.constructs for h in (0, 1) for b in sp.meta["lane_composition"][f"{c}|{h}"]}, s

    pairs = list(csv.DictReader(open(X.HERE / "candidates305.csv")))
    out = []
    P(f"\n  {'group':8s} {'construct':26s} {'gene':9s} {'e_pool':>9s} {'lane-matched':>13s} "
      f"{'h0':>8s} {'h1':>8s} {'rank':>9s} {'p_raw':>7s}  verdict")
    for r in pairs:
        c, g, grp = r["construct"], r["candidate_gene"], r["group"]
        s = r["split"] if grp != "ntc" else "train"
        sp = SP[s]
        j = sp.gpos[g]
        members = [x for x in sp.targeting if x != c] + [c]
        Cm, cov = matched(sp, Sl, lanes, members, (0, 1))
        Pm = X.L.prof(np.stack([sp.tot[x] for x in members]))
        E = Pm[:, j] - X.L.prof(Cm)[:, j]
        e_pool = float(E[-1])
        # halves against lane-matched halves of their own control parity
        Ch0, _ = matched(sp, Sl, lanes, [c], (0,))
        Ch1, _ = matched(sp, Sl, lanes, [c], (1,))
        # (halves keep their own lane composition; rebuild the weights per half)
        e_h = []
        for h, Ch in ((0, Ch0), (1, Ch1)):
            li = {b: i for i, b in enumerate(lanes)}
            par = (h,)
            Sc = np.zeros((len(lanes), Sl.shape[1]))
            Nc = np.zeros(len(lanes))
            for i in range(len(lanes)):
                Sc[i] = Sl[2 * i + h]
                Nc[i] = Sl[2 * i + h].sum()
            w = np.zeros(len(lanes))
            for b, n in sp.lanes[f"{c}|{h}"].items():
                if Nc[li[b]] > 0:
                    w[li[b]] = n / Nc[li[b]]
            e_h.append(float(X.L.prof(sp.half[(c, h)])[j] - X.L.prof(w @ Sc)[j]))
        null = E[:-1]
        rank = int((null <= e_pool).sum()) + 1
        p_raw = (1 + int((null <= e_pool).sum())) / (1 + null.size)
        conj = (e_pool < 0) and (e_pool <= MAG) and (e_h[0] < 0) and (e_h[1] < 0) and (p_raw <= P_RAW)
        orig = None
        out.append(dict(pair_id=r["pair_id"], group=grp, construct=c, candidate_gene=g, split=s,
                        e_pool_lanematched=e_pool, e_half0_lanematched=e_h[0], e_half1_lanematched=e_h[1],
                        rank=rank, n_null=int(null.size), p_raw=p_raw, meets_conjunction=bool(conj),
                        lane_coverage=float(cov[-1])))
        P(f"  {grp:8s} {c[:26]:26s} {g:9s} {'':9s} {e_pool:+13.4f} {e_h[0]:+8.4f} {e_h[1]:+8.4f} "
          f"{rank:4d}/{null.size + 1:<4d} {p_raw:7.4f}  {'MEETS C1-C4' if conj else '-'}")

    # compare with the frozen (pooled-control) reading
    frozen = {r["pair_id"]: r for r in csv.DictReader(open(X.HERE / "expression305.csv"))
              if r["group"] != "ntc" or r["split"] == "train"}
    P("\n  frozen pooled-control vs lane-matched, per pair:")
    flips = []
    for r in out:
        f = frozen[r["pair_id"]]
        was = (f["C1_negative"] == "True" and f["C2_magnitude"] == "True"
               and f["C3_both_halves"] == "True" and f["C4_unusual"] == "True")
        if was != r["meets_conjunction"]:
            flips.append(r["pair_id"])
        P(f"    {r['pair_id'][:44]:44s} e {float(f['e_pool']):+.4f} -> {r['e_pool_lanematched']:+.4f}  "
          f"conjunction {was} -> {r['meets_conjunction']}"
          + ("   <== FLIP" if was != r["meets_conjunction"] else ""))
    P(f"\n  conjunction verdict flips under lane matching: {len(flips)} ({flips})")
    X.dump(dict(rows=out, flips=flips, lanes=int(len(lanes))), X.HERE / "lanematch305.json")
    P.save()


if __name__ == "__main__":
    main()
