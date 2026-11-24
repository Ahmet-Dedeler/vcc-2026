"""300 item 4 -- the FLEX REFERENCE with an interval: the best shape-only reliability tilt (arm A) and the best
direct spread manipulation (arm C) on 273's union-OFF Flex-400 instrument (85 targets, the scorer's own per-target
statistic), against their MATCHED un-pruned refit, with a paired target bootstrap over the 85 targets.

The production operator `union_pruned` is fitted on a PRUNED, per-donor-capped row set whose exact seed-0
subsample IS reproducible (300's `prune_rows`, KA-B) but which is a three-bank object; the comparison here is
still a matched pair of refits on the same three banks and the same axis, with `flex:prod` (= 273's INCUMBENT)
carried as the known-answer row.  Everything else is 273's incumbent verbatim: base + H1-CD4 deflated k = 5
re-solved at f = 0.34 (all-rows rule), union OFF.

    PYTHONDONTWRITEBYTECODE=1 uv run python flex300.py tilt:d2.5 spread:d1.0
"""
from __future__ import annotations

import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E265 = ROOT / "experiments/265-majority-permutation-test"
E273 = ROOT / "experiments/273-flex-reference-instrument"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(E273))
sys.path.insert(0, str(E265))
import lib273 as X  # noqa: E402
import run300 as RU  # noqa: E402
import op272 as O  # noqa: E402
import boot265 as BT  # noqa: E402

plolo, L, E062 = X.plolo, X.L, X.E062
GAMMA, LAM, ALPHA, CD4_F = X.GAMMA, X.LAM, X.ALPHA, X.CD4_F
OFFICIAL, PANEL = X.OFFICIAL, list(plolo.PANEL)
CD4_H1 = X.WORK / "cd4/cd4_compact_h1"
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def main(argv: list[str]) -> None:
    t0 = time.time()
    want = argv or ["tilt:d1.0"]
    st0 = X.staged170()
    tr = X.load_truth()
    targets = st0["h1_targets"]
    cols, cand = st0["cols"], st0["cand_h1"]
    tset = set(targets.tolist())
    n = targets.size
    P(f"[flex] {n} targets, {cols.size} operator columns; arms {want}")

    op = plolo.Operator(plolo.PROD_OPS / "whitener_union_pruned.npz", name="union_pruned")
    base_prod = op.apply(cand, cols, GAMMA, LAM)
    dop0 = O.DampedOperator.from_operator(op, label="union_pruned")
    b0 = O.preserve(cand, dop0.transform(cand[:, cols], GAMMA, LAM), cols)
    ka_f0 = float(np.abs(b0 - base_prod).max())
    P(f"  KA-F0 DampedOperator+preserve == Operator.apply: max|d| {ka_f0:.3e}")
    del dop0, b0, op
    gc.collect()

    banks_h1 = [E062.Bank.load(nm, CD4_H1 / nm, require_cd4_scale=True) for nm in X.CD4_NAMES]
    prep = E062.prepare_cd4_mean(banks_h1, list(targets), OFFICIAL, tset, "gene_norm", chunk_rows=256, max_rss_gib=20.0)
    cd4h1 = X.complete_mean3(prep)
    d5, aud5 = L.deflate(cd4h1, "svd", 5)
    P(f"  H1 CD4 deflated k=5: energy removed {aud5['energy_removed_frac']:.6f} on {aud5['n_carried']} carried rows")
    del prep, banks_h1
    gc.collect()

    arms: dict[str, dict] = {}
    per_target: dict[str, np.ndarray] = {}
    prev = WORK / "flex300.json"
    if prev.exists():
        old = json.loads(prev.read_text())
        arms.update({k: v for k, v in old.get("arms", {}).items()})
        z = np.load(WORK / "pertarget300-flex.npz")
        per_target.update({k: z[k] for k in z.files if k != "targets"})
        P(f"  resuming: arms already measured {sorted(arms)}")

    def add(name: str, base: np.ndarray, **meta) -> dict:
        fr, pu = X.solve_fraction(base, d5, CD4_F)
        e = X.score_h1(X.compose(base, [(d5, fr)], cols), tr)
        per_target[f"{name}|pds"] = np.asarray(e["pds_a23"], dtype=np.float64)
        per_target[f"{name}|cos"] = np.asarray(e["cos_a23"], dtype=np.float64)
        rec = dict(pds=float(e["pds"]), cos=float(e["cos"]), a23_cos_agg=float(e["a23_cos_agg"]),
                   delivered=float(e["a23_delivered_norm"]), fraction=float(fr), per_unit=float(pu),
                   base_cos_to_prod=float(np.median(X.row_cos(base[:, cols], base_prod[:, cols]))), **meta)
        arms[name] = rec
        P(f"  [{name:18s}] pds {rec['pds']:.6f}  own cos {rec['cos']:.5f}  cos_agg {rec['a23_cos_agg']:.5f}  "
          f"delivered {rec['delivered']:.4f}  base~prod {rec['base_cos_to_prod']:.4f}  ({time.time()-t0:.0f}s)")
        return rec

    add("flex:prod", base_prod, note="273's INCUMBENT: the production operator, union OFF")
    ref = json.loads((E273 / "work/instrument273.json").read_text())
    inc = ref.get("arms", {}).get("INCUMBENT")
    ka_f1 = float("nan")
    if inc is not None:
        ka_f1 = max(abs(float(np.mean(inc[k])) - float(arms["flex:prod"][kk])) for k, kk in
                    (("pds", "pds"), ("cos_a23", "cos"), ("a23_cos_agg", "a23_cos_agg")) if k in inc)
        P(f"  KA-F1 flex:prod == 273's INCUMBENT (pds, own cos, cos_agg): max|d| {ka_f1:.3e}")

    banks = [E062.Bank.load(nm, plolo.DONOR_PATHS[nm]) for nm in X.DONORS]
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    pp = np.asarray([pos[g] for g in PANEL if g in pos], dtype=np.int64)
    raw_full = {}
    for b in banks:
        pr = E062.prepare_main_source(b, PANEL, OFFICIAL, pp, chunk_rows=256, max_rss_gib=20.0)
        raw_full[b.name] = float(pr.raw_median_norm)
        del pr
    gc.collect()
    P(f"  raw medians {raw_full}")
    rel = {b.name: RU.load_rel(b.name) for b in banks}
    stx = dict(raw_full=raw_full, fit_banks=banks)
    full_rows = {b.name: np.arange(len(b.perts), dtype=np.int64) for b in banks}
    U = {b.name: RU.tilt_u(b, rel[b.name], full_rows[b.name]) for b in banks}
    V: dict = {}

    def flush() -> None:
        o = dict(n_targets=int(n), n_cols=int(cols.size), ka=dict(KA_F0=ka_f0, KA_F1=ka_f1),
                 arms={k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")} for k, v in arms.items()},
                 seconds=round(time.time() - t0, 1), peak_rss_gib=RU.rss())
        if "flex:full" in arms:
            o["deltas"] = {}
            for a in arms:
                if a == "flex:full":
                    continue
                d = per_target[f"{a}|pds"] - per_target["flex:full|pds"]
                b = BT.boot(d, np.random.default_rng(BT.SEED + abs(hash(a)) % 99991))
                o["deltas"][a] = dict(d_pds=arms[a]["pds"] - arms["flex:full"]["pds"],
                                      d_cos=arms[a]["cos"] - arms["flex:full"]["cos"],
                                      d_cos_agg=arms[a]["a23_cos_agg"] - arms["flex:full"]["a23_cos_agg"],
                                      amp_ratio=arms[a]["delivered"] / arms["flex:full"]["delivered"],
                                      boot=b)
        (WORK / "flex300.json").write_text(json.dumps(o, indent=1, default=float))
        (WORK / "flex300.log").write_text("\n".join(LOG) + "\n")
        np.savez_compressed(WORK / "pertarget300-flex.npz", targets=np.asarray(targets), **per_target)

    flush()
    for arm in ["flex:full"] + [f"flex:{a}" for a in want]:
        wts = None
        if arm != "flex:full":
            kind, dose = arm.split(":")[1], float(arm.split(":")[2][1:])
            if kind == "spread" and not V:
                P("    arm C's basis is not available (flex:full's eigenvectors were not kept) -- skipping")
                continue
            q = U if kind == "tilt" else V
            wts = {b.name: RU.dose_weights(q[b.name], dose) for b in banks}
            for b in banks:
                w = wts[b.name]
                P(f"    {arm} {b.name}: {w.size} rows, sum w {w.sum():.3f}, n_eff "
                  f"{w.sum()**2/(w*w).sum():.0f} ({100*(w.sum()**2/(w*w).sum())/w.size:.1f} %)")
        if arm in arms and arm != "flex:full":
            P(f"  [{arm:18s}] already measured, kept")
            continue
        evals, evecs, fmeta = RU.fit300(banks, full_rows, raw_full, cols, weights_by=wts, log=P)
        if arm == "flex:full" and any(a.split(":")[0] == "spread" for a in want):
            top = np.argsort(evals)[-RU.N_PC:][::-1]
            V50 = np.ascontiguousarray(evecs[:, top])
            for b in banks:
                V[b.name], aud = RU.spread_v(stx, b, full_rows[b.name], cols, V50)
                P(f"    arm C {b.name}: top-{RU.N_PC} energy share median {aud['share_median']:.4f}")
            del V50
            gc.collect()
        dop = O.DampedOperator(evals, evecs, None, 0.0, arm)
        base = O.preserve(cand, dop.transform(cand[:, cols], GAMMA, LAM), cols)
        del dop, evecs
        gc.collect()
        add(arm, base, fit={k: v for k, v in fmeta.items() if k != "raw_medians"})
        del base
        gc.collect()
        flush()

    P("\n  deltas vs flex:full (the matched un-pruned refit), paired target bootstrap over the 85 targets:")
    o = json.loads((WORK / "flex300.json").read_text())
    for a, d in o.get("deltas", {}).items():
        P(f"    {a:20s} d pds {d['d_pds']:+.6f} [{d['boot']['lo']:+.6f}, {d['boot']['hi']:+.6f}]  "
          f"d own cos {d['d_cos']:+.6f}  d cos_agg {d['d_cos_agg']:+.6f}  amplitude x{d['amp_ratio']:.4f}")
    (WORK / "flex300.log").write_text("\n".join(LOG) + "\n")
    P(f"[flex] done {(time.time()-t0)/60:.1f} min, peak {RU.rss():.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1:])
