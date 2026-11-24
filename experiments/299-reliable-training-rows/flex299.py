"""299 item 4 -- the FLEX REFERENCE re-read: the reliability-weighted operator on 273's union-OFF Flex-400
instrument (85 targets, the scorer's own per-target statistic).

The production operator `union_pruned` is itself fitted on a PRUNED, per-donor-capped row set (12,000 rows per
donor, then the essential-gene rows dropped: 7,784 / 10,807 / 10,803 of 9,787 / 17,767 / 17,977), whose exact
seed-0 subsample is not reproducible here without guessing 062's RNG call order.  The comparison is therefore a
MATCHED PAIR of refits on the same three banks and the same axis:

    flex:full        every row of the three donor banks, 062's recipe (this is the un-pruned reference)
    flex:wgt:mild    the same rows, each weighted by 0.5 + u (u = the rank quantile of the row's split-half
                     reliability) -- run299.row_weights, the same object PLOLO reads
    flex:wgt:mildrev the reversed tilt

with `flex:prod` (the production operator, = 273's INCUMBENT) carried as the known-answer row against 273's saved
`instrument273.json`.  Everything else is 273's incumbent verbatim: base + H1-CD4 deflated k = 5 re-solved at
f = 0.34 by the all-rows rule, union OFF.

    PYTHONDONTWRITEBYTECODE=1 uv run python flex299.py
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
E273 = ROOT / "experiments/273-flex-reference-instrument"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(E273))
import lib273 as X  # noqa: E402
import run299 as RU  # noqa: E402
import op272 as O  # noqa: E402

plolo, L, E062 = X.plolo, X.L, X.E062
GAMMA, LAM, ALPHA, CD4_F = X.GAMMA, X.LAM, X.ALPHA, X.CD4_F
OFFICIAL, PANEL = X.OFFICIAL, list(plolo.PANEL)
CD4_H1 = X.WORK / "cd4/cd4_compact_h1"
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def main() -> None:
    t0 = time.time()
    st = X.staged170()
    tr = X.load_truth()
    targets = st["h1_targets"]
    cols, cand = st["cols"], st["cand_h1"]
    tset = set(targets.tolist())
    n = targets.size
    P(f"[flex] {n} targets, {cols.size} operator columns; 170 pre {st['sha']['pre'][:12]}")

    op = plolo.Operator(plolo.PROD_OPS / "whitener_union_pruned.npz", name="union_pruned")
    base_prod = op.apply(cand, cols, GAMMA, LAM)
    # KA-F0: DampedOperator + preserve == Operator.apply, bit for bit (273's KA-3, re-run on this frame)
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
    P(f"  H1 CD4 deflated k=5: energy removed {aud5['energy_removed_frac']:.6f} on {aud5['n_carried']} carried rows; "
      f"{time.time()-t0:.0f}s, peak {RU.rss():.2f} GiB")
    del prep, banks_h1
    gc.collect()

    arms: dict[str, dict] = {}

    def add(name: str, base: np.ndarray, **meta) -> dict:
        fr, pu = X.solve_fraction(base, d5, CD4_F)
        e = X.score_h1(X.compose(base, [(d5, fr)], cols), tr)
        rec = dict(pds=float(e["pds"]), cos=float(e["cos"]), a23_cos_agg=float(e["a23_cos_agg"]),
                   delivered=float(e["a23_delivered_norm"]), fraction=float(fr), per_unit=float(pu),
                   base_cos_to_prod=float(np.median(X.row_cos(base[:, cols], base_prod[:, cols]))), **meta)
        arms[name] = rec
        P(f"  [{name:16s}] pds {rec['pds']:.6f}  own cos {rec['cos']:.5f}  cos_agg {rec['a23_cos_agg']:.5f}  "
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
    else:
        P("  KA-F1 273's instrument273.json does not expose INCUMBENT in the expected shape -- reported, not asserted")

    # ---------------------------------------------------------------- the refits on the three donor banks
    banks = [E062.Bank.load(nm, plolo.DONOR_PATHS[nm]) for nm in X.DONORS]
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    pp = np.asarray([pos[g] for g in PANEL if g in pos], dtype=np.int64)
    raw_full = {}
    for b in banks:
        pr = E062.prepare_main_source(b, PANEL, OFFICIAL, pp, chunk_rows=256, max_rss_gib=20.0)
        raw_full[b.name] = float(pr.raw_median_norm)
        del pr
    gc.collect()
    P(f"  raw medians {raw_full}; {time.time()-t0:.0f}s, peak {RU.rss():.2f} GiB")
    rel = {b.name: RU.load_rel(b.name) for b in banks}

    def flush(extra: dict | None = None) -> None:
        o = dict(n_targets=int(n), n_cols=int(cols.size), ka=dict(KA_F0=ka_f0, KA_F1=ka_f1),
                 arms={k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")} for k, v in arms.items()},
                 seconds=round(time.time() - t0, 1), peak_rss_gib=RU.rss(), **(extra or {}))
        (WORK / "flex299.json").write_text(json.dumps(o, indent=1, default=float))
        (WORK / "flex299.log").write_text("\n".join(LOG) + "\n")

    flush()
    for arm in ("flex:full", "flex:wgt:mildrev"):
        rows_by, wts, auds = {}, (None if arm == "flex:full" else {}), {}
        for b in banks:
            if arm == "flex:full":
                rows_by[b.name] = np.arange(b.shape[0], dtype=np.int64)
            else:
                kind = "wgtmild" if arm.endswith("mild") else "wgtmildrev"
                idx, w, a = RU.row_weights(b, rel[b.name], kind, 299777)
                rows_by[b.name], wts[b.name], auds[b.name] = idx, w, a
                P(f"    {arm} {b.name}: {a['n']} rows, n_eff {a['n_eff']:.0f} ({100*a['n_eff_over_n']:.1f} %)")
        evals, evecs, fmeta = RU.fit299(banks, rows_by, raw_full, cols, weights_by=wts, log=P)
        dop = O.DampedOperator(evals, evecs, None, 0.0, arm)
        base = O.preserve(cand, dop.transform(cand[:, cols], GAMMA, LAM), cols)
        del dop, evecs
        gc.collect()
        add(arm, base, fit={k: v for k, v in fmeta.items() if k != "raw_medians"}, selection=auds)
        del base
        gc.collect()
        flush()

    out = dict(n_targets=int(n), n_cols=int(cols.size), ka=dict(KA_F0=ka_f0, KA_F1=ka_f1),
               arms={k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")} for k, v in arms.items()},
               deltas={a: dict(d_pds=arms[a]["pds"] - arms["flex:full"]["pds"],
                               d_cos=arms[a]["cos"] - arms["flex:full"]["cos"],
                               d_cos_agg=arms[a]["a23_cos_agg"] - arms["flex:full"]["a23_cos_agg"],
                               amp_ratio=arms[a]["delivered"] / arms["flex:full"]["delivered"])
                       for a in arms if a != "flex:full"},
               seconds=round(time.time() - t0, 1), peak_rss_gib=RU.rss())
    (WORK / "flex299.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "flex299.log").write_text("\n".join(LOG) + "\n")
    P("\n  deltas vs flex:full (the matched un-pruned refit):")
    for a, d in out["deltas"].items():
        P(f"    {a:16s} d pds {d['d_pds']:+.6f}  d own cos {d['d_cos']:+.6f}  d cos_agg {d['d_cos_agg']:+.6f}  "
          f"amplitude x{d['amp_ratio']:.4f}")
    P(f"[flex] done {(time.time()-t0)/60:.1f} min, peak {RU.rss():.2f} GiB")


if __name__ == "__main__":
    main()
