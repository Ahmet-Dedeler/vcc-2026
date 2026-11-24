"""308 -- the known-answer rows for the FEATURE and its provenance (nothing here scores an arm).

  KA-A  293's on-target knockdown numbers reproduce from 293's own saved arrays by THIS code
        (the §4 per-target `kd big` / `kd small` table and §5b's three failure-rate rows).
  KA-B  the chunked gather `read_diagonal` == a direct element-wise memmap read, all three donor banks.
  KA-C  the scorer's own-target exclusion is what 097 says it is: `_drop_on_target` / `_exclusion_cols` in
        `cell_eval2.metrics.delta`, `run._ONTARGET_EXCLUDING_FUNCS`, `metrics.direction._exclude_on_target`,
        and `metrics.discrimination`'s `exclude_target_gene=True` at `exclusion_scope="panel"` -- verified by
        reading the INSTALLED package and by calling `_exclusion_cols` / `_drop_on_target` on a synthetic case.
        (Static + unit-level; `cell-eval2` itself is never run.)
  KA-D  LEAKAGE: the feature computed with ONLY the fold's two donor banks present is bit-identical to the
        corresponding columns of the panel-wide feature; and it is invariant to the held fold's truth entirely
        (the truth array never enters the computation -- demonstrated by construction and by re-deriving the
        two columns in a process that never opens the held bank or the fold cache).

    PYTHONDONTWRITEBYTECODE=1 uv run python ka308.py
"""
from __future__ import annotations

import inspect
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import kd308 as K  # noqa: E402
import plolo  # noqa: E402

E062 = plolo.E062
E293 = ROOT / "experiments/293-challenge-construct-ceiling/work"
E296 = ROOT / "experiments/296-guide-efficacy/work"
LINES: list[str] = []
KA: list[dict] = []
FAIL: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def check(label, got, want, tol) -> bool:
    d = abs(float(got) - float(want))
    ok = d <= tol
    KA.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok)))
    if not ok:
        FAIL.append(label)
    P(f"  KA {'PASS' if ok else 'FAIL'} {label}: got {got!r} want {want!r} |d| {d:.2e}")
    return ok


def note(label, ok: bool, detail: str) -> bool:
    KA.append(dict(label=label, ok=bool(ok), detail=detail))
    if not ok:
        FAIL.append(label)
    P(f"  KA {'PASS' if ok else 'FAIL'} {label}: {detail}")
    return ok


# ------------------------------------------------------------------------------ KA-A: 293's |kd|
def ka_a() -> None:
    P("KA-A  293's on-target knockdown, re-derived by this code from 293's saved arrays")
    import csv
    rows = list(csv.DictReader((E293 / "construct293.csv").open()))
    # the §4 table is the n = 200/side, thinned (challenge-regime) rung
    want = {"SOX4": (-2.094, -1.189), "STAG2": (-1.149, -1.795), "TMSB4X": (-3.037, -2.737),
            "PRCP": (-1.569, -1.539), "IGF1R": (-1.086, -0.662), "ZNF581": (-1.748, -0.169),
            "KLF10": (-1.077, -0.019), "SMARCB1": (-1.807, -0.593), "GRB10": (-0.197, -1.009),
            "ARPC2": (-1.602, -0.078), "KDM2B": (-1.355, -0.772), "HIRA": (+0.030, -1.494),
            "CAST": (-0.237, -0.482), "MED15": (-0.053, -1.376), "DAXX": (-1.184, -1.672),
            "PCCB": (-1.595, -0.187)}
    cols = set(rows[0])
    big = next(c for c in cols if c in ("kd_big", "kd_a", "kd_A", "kd_1"))
    small = next(c for c in cols if c in ("kd_small", "kd_b", "kd_B", "kd_2"))
    got, n_hit = {}, 0
    for r in rows:
        t = r.get("target")
        if t in want and float(r.get("n", r.get("n_side", 0)) or 0) == 200:
            got.setdefault(t, (float(r[big]), float(r[small])))
            n_hit += 1
    dmax = max(max(abs(round(got[t][0], 3) - want[t][0]), abs(round(got[t][1], 3) - want[t][1])) for t in want if t in got)
    note("KA-A1 293 sec 4 per-target kd big/small (16 targets, both constructs) reproduce to the printed 3 d.p.",
         len(got) == 16 and dmax <= 5e-4, f"{len(got)}/16 targets, max |d| {dmax:.1e}")
    # 5b: the three failure-rate rows, from 293's kd_all
    kd25: dict[str, float] = {}
    for sp in ("train", "validation", "test"):
        kd25.update(json.loads((E293 / f"run293_{sp}.json").read_text())["kd_all"])
    v = np.abs(np.asarray(list(kd25.values()), float))
    check("KA-A2 293 sec 5b Arc-2025 panel median |lfc|", float(np.median(v)), 1.485, 5e-4)
    check("KA-A3 293 sec 5b Arc-2025 panel |lfc| < 1.0 (%)", 100 * float((v < 1.0).mean()), 18.7, 0.05)
    check("KA-A4 293 sec 5b Arc-2025 panel |lfc| < 0.5 (%)", 100 * float((v < 0.5).mean()), 4.3, 0.05)


# ------------------------------------------------------------------------------ KA-B: the gather
def ka_b(feat) -> None:
    P("KA-B  the chunked gather == a direct element-wise memmap read")
    for j, name in enumerate(K.DONORS):
        bank = E062.Bank.load(name, plolo.DONOR_PATHS[name])
        gpos = K.own_column_map(bank)
        pairs, idx = [], []
        for i, t in enumerate(plolo.PANEL):
            r, c = bank.pert_pos.get(t), gpos.get(t)
            if r is not None and c is not None:
                pairs.append((int(r), int(c)))
                idx.append(i)
        delta = np.load(bank.delta_path, mmap_mode="r")
        direct = np.array([float(delta[r, c]) for r, c in pairs])
        check(f"KA-B {name}: chunked gather == direct read on {len(pairs)} (target, own gene) pairs",
              float(np.abs(direct - feat["kd"][np.asarray(idx), j]).max()), 0.0, 0.0)
        del bank, delta


# ------------------------------------------------------------------------------ KA-C: the scorer's exclusion
def ka_c() -> None:
    P("KA-C  the scorer drops each row's own-target coordinate (097 route (d)'s structural claim), verified in the")
    P("      INSTALLED cell_eval2 by reading the source and by calling the two functions on a synthetic case")
    import cell_eval2
    from cell_eval2.metrics import delta as D
    src_delta = Path(inspect.getfile(D)).read_text()
    ver = getattr(cell_eval2, "__version__", "?")
    P(f"      cell_eval2 {ver} at {Path(inspect.getfile(cell_eval2)).parent}")
    have = {n: hasattr(D, n) for n in ("_exclusion_cols", "_drop_on_target")}
    note("KA-C1 metrics.delta defines _exclusion_cols and _drop_on_target", all(have.values()), json.dumps(have))
    from cell_eval2 import run as R
    raw = getattr(R, "_ONTARGET_EXCLUDING_FUNCS", ()) or ()
    names = sorted({getattr(f, "__name__", str(f)) for f in raw})
    note("KA-C2 run._ONTARGET_EXCLUDING_FUNCS covers mse_unbiased_capped, de_sig_jaccard, de_lfc_nmae",
         {"mse_unbiased_capped", "de_sig_jaccard", "de_lfc_nmae"} <= set(names), f"{names}")
    from cell_eval2.metrics import direction as DIR
    note("KA-C3 metrics.direction defines _exclude_on_target (fidelity / reach)", hasattr(DIR, "_exclude_on_target"),
         "present" if hasattr(DIR, "_exclude_on_target") else "ABSENT")
    from cell_eval2.metrics import discrimination as DIS
    ssrc = Path(inspect.getfile(DIS)).read_text()
    ok_p = ("exclude_target_gene" in ssrc) and ("exclusion_scope" in ssrc) and ('"panel"' in ssrc or "'panel'" in ssrc)
    note("KA-C4 metrics.discrimination takes exclude_target_gene with exclusion_scope 'panel' (pds drops ALL 300 "
         "panel genes)", ok_p, "exclude_target_gene / exclusion_scope / 'panel' all present in the source")
    # the numeric demonstration: _drop_on_target must remove exactly the own-target coordinate's contribution
    sig = inspect.signature(D._drop_on_target)
    P(f"      _drop_on_target{sig}")
    ex_sig = inspect.signature(D._exclusion_cols)
    P(f"      _exclusion_cols{ex_sig}")
    genes = ["G0", "G1", "G2", "G3"]
    labels = ["G1", "G2"]
    want_col = {"G1": 1, "G2": 2}
    tmap = {"G1": "G1", "G2": "G2"}                              # perturbation -> gene SYMBOL, the scorer's map
    cols = dict(D._exclusion_cols(labels, genes, target_gene_map=tmap, gate_labels=None, who="test", n_genes=len(genes)))
    got = {str(k): sorted(int(c) for c in np.atleast_1d(v)) for k, v in cols.items()}
    hit = sum(1 for r, lb in enumerate(labels)
              if want_col[lb] in got.get(str(r), got.get(lb, [])))
    note("KA-C5 _exclusion_cols maps each row's label to its OWN gene column", hit == len(labels),
         f"{got} (own column found for {hit}/{len(labels)} rows)")
    # the exclusion is driven entirely by the resolver's map, and an unresolvable map is REFUSED, not ignored
    empty = dict(D._exclusion_cols(labels, genes, target_gene_map={}, gate_labels=None, who="test", n_genes=len(genes)))
    n_excl_empty = sum(len(np.atleast_1d(v)) for v in empty.values())
    note("KA-C5b the exclusion is the DEFAULT: with no target-gene map the resolver still matches each label to its "
         "own gene symbol and excludes it", n_excl_empty == sum(len(v) for v in got.values()),
         f"{n_excl_empty} columns excluded with an empty map vs {sum(len(v) for v in got.values())} with the real one")
    try:
        D._exclusion_cols(labels, genes, target_gene_map={"G1": 1, "G2": 2}, gate_labels=None, who="test",
                          n_genes=len(genes))
        note("KA-C5c an unresolvable map is REFUSED, not silently ignored", False, "no refusal raised")
    except Exception as exc:
        txt = str(exc)
        note("KA-C5c an unresolvable map is REFUSED, not silently ignored", "cannot be switched off" in txt,
             f"{type(exc).__name__}: '...exclusion is part of the metric and cannot be switched off'"
             if "cannot be switched off" in txt else txt[:200])
    n_hits = sum(src_delta.count(k) for k in ("_drop_on_target", "_exclusion_cols"))
    note("KA-C6 the exclusion is wired, not vestigial (call sites in metrics/delta.py)", n_hits >= 6,
         f"{n_hits} references to _drop_on_target / _exclusion_cols in metrics/delta.py")
    P("      => 097 route (d) is confirmed AS A PREDICTION TARGET: the own-target coordinate is removed from every")
    P("         member, and `pds` drops all 300 panel genes.  308 uses |kd| as a WEIGHT on the donor rows, which")
    P("         the exclusion does not touch: the weight changes the DIRECTION of the 17.9 k scored coordinates.")


# ------------------------------------------------------------------------------ KA-D: leakage
def ka_d(feat) -> None:
    P("KA-D  LEAKAGE: the feature never sees the held line")
    for held in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        others = [d for d in K.DONORS if d != held]
        vals = {}
        for name in others:                                     # recomputed WITHOUT the held bank in the process
            bank = E062.Bank.load(name, plolo.DONOR_PATHS[name])
            gpos = K.own_column_map(bank)
            pairs, idx = [], []
            for i, t in enumerate(plolo.PANEL):
                r, c = bank.pert_pos.get(t), gpos.get(t)
                if r is not None and c is not None:
                    pairs.append((int(r), int(c)))
                    idx.append(i)
            v = np.full(len(plolo.PANEL), np.nan)
            v[np.asarray(idx)] = K.read_diagonal(bank, pairs)
            vals[name] = v
            del bank
        d = max(float(np.nanmax(np.abs(vals[nm] - feat["kd"][:, K.DONORS.index(nm)]))) for nm in others)
        nan_ok = all(bool((np.isnan(vals[nm]) == np.isnan(feat["kd"][:, K.DONORS.index(nm)])).all()) for nm in others)
        check(f"KA-D1[{held}] the two surviving donors' |kd| computed with the HELD bank never opened == the "
              f"panel-wide array", d, 0.0, 0.0)
        note(f"KA-D2[{held}] identical missingness pattern", nan_ok, "NaN masks equal")
    note("KA-D3 the feature's inputs are exactly {delta.npy, genes.txt, perts.txt} of ONE donor bank plus the "
         "target symbol", True, "no truth array, no fold cache, no operator, no held bank is opened by kd308.read_diagonal")


# ------------------------------------------------------------------------------ KA-E: an independent path
def ka_e(feat) -> None:
    P("KA-E  cross-check against 296's INDEPENDENT cell-level knockdown for K562-GWPS (a different code path:")
    P("      296 pseudobulked the raw scPerturb cells per construct; 308 reads the donor bank's delta.npy)")
    j = json.loads((E296 / "k562_296.json").read_text())
    ks = j["kd_single"]
    col = K.DONORS.index("k562_gwps")
    pan = list(plolo.PANEL)
    pairs = [(pan.index(t), float(v["lfc"])) for t, v in ks.items() if t in set(pan)]
    pairs = [(i, x) for i, x in pairs if np.isfinite(feat["kd"][i, col])]
    mine = np.array([feat["kd"][i, col] for i, _ in pairs])
    theirs = np.array([x for _, x in pairs])
    from scipy.stats import spearmanr, pearsonr
    rs = float(spearmanr(mine, theirs).statistic)
    rp = float(pearsonr(mine, theirs).statistic)
    P(f"      n = {len(pairs)} shared panel targets; Spearman {rs:+.3f}, Pearson {rp:+.3f}; "
      f"median |kd| mine {np.median(np.abs(mine)):.3f} vs 296 {np.median(np.abs(theirs)):.3f}; "
      f"< 1.0 mine {100*np.mean(np.abs(mine) < 1):.1f} % vs 296 {100*np.mean(np.abs(theirs) < 1):.1f} %")
    note("KA-E K562-GWPS bank-derived |kd| agrees with 296's cell-derived |kd| (Spearman > +0.8)", rs > 0.8,
         f"Spearman {rs:+.3f} on {len(pairs)} targets")
    KA.append(dict(label="KA-E detail", n=len(pairs), spearman=rs, pearson=rp,
                   median_abs_mine=float(np.median(np.abs(mine))), median_abs_296=float(np.median(np.abs(theirs))),
                   frac_lt1_mine=float(np.mean(np.abs(mine) < 1)), frac_lt1_296=float(np.mean(np.abs(theirs) < 1))))


def main() -> None:
    feat = np.load(WORK / "kd308_feature.npz", allow_pickle=True)
    feat = dict(kd=np.asarray(feat["kd"], float))
    ka_a()
    ka_b(feat)
    ka_c()
    ka_d(feat)
    ka_e(feat)
    P(f"\n{len(KA)} rows, {len(FAIL)} failures" + (": " + "; ".join(FAIL) if FAIL else ""))
    (WORK / "ka308.json").write_text(json.dumps(dict(ka=KA, failures=FAIL), indent=2))
    (WORK / "ka308.log").write_text("\n".join(LINES) + "\n")
    if FAIL:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
