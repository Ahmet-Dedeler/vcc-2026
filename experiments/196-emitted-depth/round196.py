"""196 Phase 0(b) — the ANATOMY of the integer rounding, measured on the shipped totals.

The quasi emitter delivers `T_g = floor(x_g)`, `x_g = n * L * p_g`.  Conditional on the
delivered total, the fractional part `u = x - T` is EXACTLY uniform on [0,1) for any smooth
profile (equidistribution; it fails only where `x < 1`, which is far below the scorer's
5-CPM gate at `x = 8*cpm`).  So every rounding statistic below is a closed form in `T`:

    E[ u/(T+u) ]        = 1 - T ln((T+1)/T)                 relative shortfall of the total
    E[ log2(T/(T+u)) ]  = -[(T+1) ln((T+1)/T) - 1] / ln 2   log2 bias of the delivered lfc

and NOTHING here needs the model's effect vector.  What the model does supply -- and what
this file also reports -- is the common-mode: `q_g = T_g / sum(T)` is renormalised by the
scorer, so the panel-wide `sum(u)/sum(T)` comes back as a POSITIVE shift on every gene.  The
net per-gene tilt is therefore `-u/(T+u) + sum(u)/sum(T)`, zero at one crossover expression
and signed either side of it.  That crossover is the number this file exists to produce.

    uv run python round196.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
ROOT = HERE.parent.parent
LN2 = float(np.log(2.0))


def e_rel(T):
    """E[u/(T+u)] -- expected relative shortfall of the delivered total. T=0 -> 1."""
    T = np.asarray(T, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = 1.0 - T * np.log1p(1.0 / np.maximum(T, 1e-300))
    return np.where(T <= 0, 1.0, v)


def e_log2(T):
    """E[log2(T/(T+u))] -- expected log2 bias of the delivered per-gene mean. T=0 -> -inf."""
    T = np.asarray(T, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = -((T + 1.0) * np.log1p(1.0 / np.maximum(T, 1e-300)) - 1.0) / LN2
    return np.where(T <= 0, -np.inf, v)


def deciles(cpm, gate):
    """Ten equal-count bins of the GATED genes by control CPM, low to high."""
    idx = np.flatnonzero(gate)
    order = idx[np.argsort(cpm[idx], kind="stable")]
    return np.array_split(order, 10)


def table(name, T, cpm, gate, n, L, fh):
    """One surface's decile table.  `T` is (blocks, genes) delivered totals."""
    Tg = T[:, gate].astype(np.float64)
    common = float(np.mean(T.sum(1) / (n * L)))          # sum(T)/(nL) = 1 - sum(u)/(nL)
    shift = 1.0 / common - 1.0                           # the positive renormalisation
    bins = deciles(cpm, gate)
    rows = []
    print(f"\n### {name}   n={n} L={L}  n*L={n*L:,}  gated={int(gate.sum())}  "
          f"blocks={T.shape[0]}", file=fh)
    print(f"    sum(T)/(nL) = {common:.6f}  ->  common-mode profile shift "
          f"{shift*1e2:+.4f}%  ({shift/LN2*1e3:+.3f} millilog2)", file=fh)
    print(f"\n| decile | ctrl CPM | median T | share T<100 | E[rel shortfall] | "
          f"E[log2 bias] | net tilt (log2) |", file=fh)
    print("|---|---|---:|---:|---:|---:|---:|", file=fh)
    for k, b in enumerate(bins):
        sub = T[:, b].astype(np.float64)
        lo, hi = cpm[b].min(), cpm[b].max()
        med = float(np.median(sub))
        frac100 = float((sub < 100).mean())
        rel = float(np.mean(e_rel(sub)))
        lg = float(np.mean(e_log2(np.maximum(sub, 1e-12))[np.isfinite(e_log2(np.maximum(sub, 1e-12)))]))
        net = lg + shift / LN2
        rows.append(dict(decile=k + 1, cpm_lo=float(lo), cpm_hi=float(hi), median_T=med,
                         share_below_100=frac100, e_rel=rel, e_log2=lg, net_log2=net))
        print(f"| {k+1} | {lo:.1f}–{hi:.1f} | {med:.0f} | {frac100*100:.1f}% | "
              f"{rel*100:.3f}% | {lg:+.5f} | {net:+.5f} |", file=fh)
    allg = T[:, gate].astype(np.float64)
    over = dict(
        share_below_100=float((allg < 100).mean()),
        share_below_40=float((allg < 40).mean()),
        share_below_10=float((allg < 10).mean()),
        share_zero=float((allg == 0).mean()),
        median_T=float(np.median(allg)),
        mean_rel=float(np.mean(e_rel(allg))),
        common_shift=shift,
        crossover_T=float(1.0 / shift * 0.5) if shift > 0 else float("nan"),
    )
    over["crossover_cpm"] = over["crossover_T"] / (n * L / 1e6)
    print(f"\n    gated pool: median T {over['median_T']:.0f}; share below 100 counts "
          f"**{over['share_below_100']*100:.1f}%**, below 40 {over['share_below_40']*100:.1f}%, "
          f"below 10 {over['share_below_10']*100:.2f}%, exactly zero "
          f"{over['share_zero']*100:.3f}%", file=fh)
    print(f"    mean relative shortfall over the gated pool "
          f"**{over['mean_rel']*100:.3f}%**", file=fh)
    print(f"    tilt crossover: T = {over['crossover_T']:.0f} counts = "
          f"{over['crossover_cpm']:.1f} control CPM  (genes below it are delivered LOW, "
          f"above it HIGH)", file=fh)
    return dict(rows=rows, overall=over)


def main() -> None:
    out = {}
    lines = []

    class Tee:
        def write(self, s):
            lines.append(s)
        def flush(self):
            pass
    fh = Tee()

    # ---------------- h1-test (the surface with truth) ----------------
    E = np.load(ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz",
                allow_pickle=True)
    Th = E["sums"].astype(np.int64)
    cpm_h = E["cpm"].astype(np.float64)
    meta = np.load(ROOT / "experiments/185-rank-controlled-emission/work/ctrl_meta.npz",
                   allow_pickle=True)
    gate_h = meta["keep"]                       # the SCORER's gate: mean-per-cell CPM > 5
    Lh = int(np.median(meta["ctrl_lib"]))
    out["h1_test"] = table("h1-test (a15-analogue artifact, 85 perturbations)",
                           Th, cpm_h, gate_h, 400, Lh, fh)

    # ---------------- challenge A/B/C (the real package) ----------------
    z = np.load(WORK / "chal-totals-a15.npz", allow_pickle=True)
    blocks = np.asarray(z["blocks"]).astype(str)
    T = z["totals"]
    ctxs = np.array([b.split("|")[0] for b in blocks])
    summ = json.load(open(WORK / "ctrl-summary.json"))
    for c in "ABC":
        cz = np.load(WORK / f"ctrl-{c}.npz")
        mcpm = cz["mean_cpm"].astype(np.float64)
        gate = mcpm > 5.0
        out[f"chal_{c}"] = table(f"challenge context {c} (a15 as shipped, 300 perturbations)",
                                 T[ctxs == c], mcpm, gate, 400, summ[c]["L"], fh)
    (WORK / "round196.md").write_text("".join(lines))
    json.dump(out, open(WORK / "round196.json", "w"), indent=1)
    print("".join(lines))


if __name__ == "__main__":
    main()
