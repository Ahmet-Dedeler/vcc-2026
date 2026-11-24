"""612 ds: the coordinator's 'disagreement shrink' (shrink |lfc| where sign sources disagree, never flip), priced
delivery-faithfully on the folds in kgh612c's shape (base = BETA * donor consensus, plus a hook-signed MDE add), with the
H1 power ladder (h1_c100 / h1_c150 / h1_c250 / h1_thin / h1) as analogs of 521's regimes.

Why: on the H1 real emitter (612/diag612b), the plain hook's opposite add on STRONG-key genes (a shrink where CLEAN4
disagrees with the key) helped at f 0.35 (+0.0011 raw nmae vs idkg) but not at f 1. On the full-power folds the key gate,
which removes that shrink, gained on all four lines. This asks whether the sign of the strong-half shrink turns with power,
and whether a never-flip shrink keeps the gain without the tail.

Forms (all base + s*add unless noted; s = the key sign sc where a source has no sign):
  nohook      s = sc
  plain       s = CLEAN4 (tau 0.05)
  kgC4        weak half: CLEAN4 (tau 0); strong half: sc
  kgID        weak half: identity S (622 form, tau 0); strong half: sc          <- the current recommendation
  IDw_C4s     weak half: identity S; strong half: CLEAN4 (tau 0.05)             <- kgID + the plain hook's strong shrink
  plainID     identity S everywhere (ungated)
  *_nf        the same with a never-flip clamp: pred -> 0 where it would take the opposite sign of the base
  msh_<D>_<r> magnitude shrink: rho * base (no add) where source D disagrees with sc; nohook elsewhere
The decomposition splits each form's gain over nohook EXACTLY (nmae rows are sums over cells) into cells where the form
flips the base sign, shrinks it, grows it, or signs a gene the key left at 0.
Analysis only; < 3 GB; no build. Output: work/ds612.out, work/ds612.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
for p in ("experiments/480-nmae-instrument", "experiments/537-arm-pipeline", "experiments/622-identity-track"):
    sys.path.insert(0, str(ROOT / p))
import fold480 as F  # noqa: E402
from common480 import boot_delta, load_donor_rows, load_h1  # noqa: E402
import foldsign537 as FS  # noqa: E402
import hyb622 as H  # noqa: E402

C4ALL = ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"]
SURFS = ["h1_c100", "h1_c150", "h1_c250", "h1_thin", "h1", "orion_hct116", "orion_hek293t", "k562_gwps"]
BANKS = {s: (["nadig_jurkat", "nadig_hepg2"] if s == "k562_gwps" else C4ALL) for s in SURFS}
PROFKEY = {s: ("h1" if s.startswith("h1") else s) for s in SURFS}


def nmae_cells(pred, T):
    """Per-target nmae and the per-cell numerators (|P - T| / den on scored cells, 0 elsewhere)."""
    P = np.nan_to_num(np.asarray(pred, float))
    G = np.isfinite(T)
    Tz = np.where(G, T, 0.0)
    den = np.abs(Tz).sum(1)
    ok = G.any(1) & (den > 0)
    cell = np.where(G, np.abs(P - Tz), 0.0) / np.where(den > 0, den, 1.0)[:, None]
    rows = np.where(ok, cell.sum(1), np.nan)
    return rows, cell, ok


def setup():
    """BETA (kgh612c's magnitude map, from H1) and the basal profiles for the identity axis."""
    h = load_h1()
    G = np.isfinite(h["T"])
    L0 = np.nan_to_num(np.asarray(np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz")["lfc"],
                                  float))
    c3 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, h["targets"], h["genes"], pc_cpm=1.0)[0]
                                            for d in F.DON]), 0))
    BETA = float(np.median(np.abs(L0[G])) / np.median(np.abs(c3[G])))
    print(f"BETA {BETA:.4f} (kgh612c: 0.5025)")
    return BETA, H.profiles()


def surface_forms(name, BETA, prof):
    """Every form's prediction on one fold surface. Returns dict(T, P, weak, sc, S, cons0)."""
    S = F.load_surface(name)
    T, cnt, genes, targets = S["T"], S["cnt"], [str(g) for g in S["genes"]], S["targets"]
    cons0 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, np.array(genes), pc_cpm=1.0)[0]
                                               for d in S["donors"]]), 0))
    P4 = FS.clean4_P(np.array(genes), BANKS[name])
    s4 = np.where(np.isfinite(P4) & (np.abs(P4 - 0.5) >= 0.05), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
    s40 = np.where(np.isfinite(P4) & (P4 != 0.5), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
    zc = np.where(np.isfinite(P4), (P4 - 0.5) / np.nanstd(P4 - 0.5), np.nan)
    pk = PROFKEY[name]
    zi = H.identity(prof, pk, [k for k in ("ctxA", "ctxB", "ctxC", "h1") if k != pk], genes)
    sc_low = np.where(np.isfinite(zc), zc + 0.5 * np.nan_to_num(zi), np.nan_to_num(zi))
    sI = np.where(np.isfinite(sc_low) & (sc_low != 0), np.where(sc_low > 0, -1.0, 1.0), 0.0)
    S4, S40, SI = (np.broadcast_to(v, T.shape) for v in (s4, s40, sI))
    sc = np.sign(cons0)
    add = F.mde_add(cnt)[None, :]
    weak = np.zeros(T.shape, bool)
    for i in range(T.shape[0]):
        a = np.abs(cons0[i])
        nz = a > 0
        weak[i] = (~nz) | (a < (np.quantile(a[nz], 0.5) if nz.any() else 0.0))
    base = BETA * cons0
    hk = lambda s: base + np.where(s != 0, s, sc) * add  # noqa: E731
    nf = lambda p: np.where((sc != 0) & (np.sign(p) * sc < 0), 0.0, p)  # noqa: E731
    P = {"nohook": hk(sc), "plain": hk(S4), "kgC4": hk(np.where(weak, S40, 0)),
         "kgID": hk(np.where(weak, SI, 0)), "IDw_C4s": hk(np.where(weak, SI, S4)), "plainID": hk(SI)}
    for k in ("plain", "kgID", "IDw_C4s", "plainID"):
        P[k + "_nf"] = nf(P[k])
    dis = {"C4": (S4 != 0) & (S4 != sc) & (sc != 0), "ID": (SI != 0) & (SI != sc) & (sc != 0)}
    dis["C4andID"] = dis["C4"] & dis["ID"]
    for d, m in dis.items():
        for r in (0.5, 0.0):
            P[f"msh_{d}_{r:g}"] = np.where(m, r * base, P["nohook"])
    # kgID on the weak half + a magnitude shrink on the strong half where CLEAN4 disagrees (the 'both' candidate)
    for r in (0.5, 0.0):
        P[f"kgID+mshC4s_{r:g}"] = np.where(~weak & dis["C4"], r * base, P["kgID"])
    return dict(T=T, P=P, weak=weak, sc=sc, S=S, cons0=cons0)


def main() -> None:
    BETA, prof = setup()
    out, lines = {}, []

    def say(s=""):
        print(s)
        lines.append(s)

    for name in SURFS:
        D = surface_forms(name, BETA, prof)
        T, P, weak, sc = D["T"], D["P"], D["weak"], D["sc"]
        rows, cells = {}, {}
        for k, v in P.items():
            rows[k], cells[k], ok = nmae_cells(v, T)
        n_ok = int(ok.sum())
        res = {k: boot_delta(rows[k], rows["nohook"]) for k in P if k != "nohook"}
        vs = {k: boot_delta(rows[k], rows["kgID"]) for k in P if k not in ("nohook", "kgID")}
        # exact decomposition of the gain over nohook by cell class (mean over scored targets; + = better)
        dec = {}
        p0 = P["nohook"]
        for k in ("plain", "kgID", "IDw_C4s", "plainID", "kgID+mshC4s_0.5"):
            p1 = P[k]
            d = (cells["nohook"] - cells[k])[ok]
            cls = {"flip": (sc != 0) & (np.sign(p1) * sc < 0),
                   "shrink": (sc != 0) & (np.sign(p1) * sc >= 0) & (np.abs(p1) < np.abs(p0) - 1e-15),
                   "grow": (sc != 0) & (np.sign(p1) * sc >= 0) & (np.abs(p1) > np.abs(p0) + 1e-15),
                   "new": (sc == 0) & (p1 != p0)}
            dec[k] = {c: float(np.where(m[ok], d, 0.0).sum(1).mean()) for c, m in cls.items()}
            dec[k]["total"] = float(d.sum(1).mean())
            dec[k]["strong"] = float(np.where(~weak[ok], d, 0.0).sum(1).mean())
            dec[k]["weak"] = float(np.where(weak[ok], d, 0.0).sum(1).mean())
        out[name] = {"n_scored": n_ok, "vs_nohook": res, "vs_kgID": vs, "decomp": dec}
        say(f"\n== {name}: {n_ok} scored targets; Δnmae vs no hook (+ = better) [95% target bootstrap]")
        for k, v in res.items():
            w = vs.get(k)
            say(f"  {k:18s} {v[0]:+.5f} [{v[1]:+.5f},{v[2]:+.5f}]" + (f"   vs kgID {w[0]:+.5f} [{w[1]:+.5f},{w[2]:+.5f}]"
                                                                        if w else ""))
        say("  decomposition of the gain over nohook (flip / shrink / grow / new | weak / strong = total):")
        for k, v in dec.items():
            say(f"    {k:16s} " + " ".join(f"{c} {v[c]:+.5f}" for c in ("flip", "shrink", "grow", "new"))
                + f" | weak {v['weak']:+.5f} strong {v['strong']:+.5f} = {v['total']:+.5f}")
    (HERE / "work/ds612.out").write_text("\n".join(lines) + "\n")
    (HERE / "work/ds612.json").write_text(json.dumps({"BETA": BETA, "surfaces": out}, indent=1))
    say("\nwrote work/ds612.out, work/ds612.json")


if __name__ == "__main__":
    main()
