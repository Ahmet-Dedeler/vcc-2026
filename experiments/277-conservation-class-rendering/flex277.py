"""277 -- the SAME class-conditional renderings on 273's Flex-400 reference (h1-test, union OFF).

    PYTHONDONTWRITEBYTECODE=1 uv run python flex277.py

273's instrument is re-used verbatim (`lib273`: 170's staged `cand_h1`, the production `union_pruned`
operator, the CD4 bank re-extracted for h1-test's 85 targets, `score_h1` = the scorer's own per-target
`pds_cosine`).  The class is `h1_c_lolo` from `experiments/276-conservation-map/work/predict276.npz` --
276's conservation score rebuilt on h1-test's own targets with BOTH H1 banks removed
(`predict276.py:263`, `lib276.HELD_LINE_SOURCES["H1"]`), i.e. already the H1-excluded score; it is used
as-is, not recomputed.  Bottom 30 % of the targets with a finite score (DEFINITION-277 section 1).

KNOWN-ANSWER ROWS (exit 1 on failure): the re-composed `BASE`, `cd4h1:off` and `INCUMBENT` per-target
`pds` (a23, a15, effect level) and own cosines reproduce `work/instrument273.npz` at max|d| == 0.0.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib277 as Z  # noqa: E402

sys.path.insert(0, str(Z.E273))
import lib273 as X3  # noqa: E402  (re-pins plolo.WORK to 273's cache; restored below for anything we stage)

L, plolo, E062 = X3.L, X3.plolo, X3.E062
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def main() -> None:
    t0 = time.time()
    fail: list[str] = []
    st = X3.staged170()
    tr = X3.load_truth()
    targets = [str(t) for t in st["h1_targets"]]
    n = len(targets)
    cols, cand = st["cols"], st["cand_h1"]
    P(f"277 flex277.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}; {n} h1-test targets, "
      f"{cols.size} operator columns; 170 pre.npz {st['sha']['pre'][:12]}")

    op = plolo.Operator(plolo.PROD_OPS / "whitener_union_pruned.npz", name="union_pruned")
    base = op.apply(cand, cols, X3.GAMMA, X3.LAM)
    del op
    P(f"  base whitened ({time.time()-t0:.0f}s, peak {Z.peak_gib():.2f} GiB)")

    # --------------------------------------------------------------------------- the CD4 object (273's)
    banks = [E062.Bank.load(nm, X3.WORK / "cd4/cd4_compact_h1" / nm, require_cd4_scale=True)
             for nm in X3.CD4_NAMES]
    prep = E062.prepare_cd4_mean(banks, targets, X3.OFFICIAL, set(targets), "gene_norm",
                                 chunk_rows=256, max_rss_gib=40.0)
    cd4 = X3.complete_mean3(prep)
    carried = np.asarray(cd4.target_condition_count) == 3
    d5, aud5 = L.deflate(cd4, "svd", 5)
    fr_off, _ = X3.solve_fraction(base, cd4, X3.CD4_F)
    fr5, _ = X3.solve_fraction(base, d5, X3.CD4_F)
    P(f"  CD4 (re-extracted for the 85 targets): {int(carried.sum())}/{n} complete rows; k = 5 removes "
      f"{aud5['energy_removed_frac']:.6f}; fractions off {fr_off:.6f} / svd5 {fr5:.6f}")

    # --------------------------------------------------------------------------- the class
    zp = np.load(Z.E276 / "work/predict276.npz", allow_pickle=True)
    hmap = {str(t): float(v) for t, v in zip(zp["h1_targets"], zp["h1_c_lolo"])}
    cl = np.asarray([hmap.get(t, np.nan) for t in targets], dtype=np.float64)
    C30 = Z.bottom_class(cl, targets, 0.30)
    C20 = Z.bottom_class(cl, targets, 0.20)
    C40 = Z.bottom_class(cl, targets, 0.40)
    TOP30 = Z.bottom_class(cl, targets, 0.30, top=True)
    dd_h = np.asarray([{str(t): float(v) for t, v in zip(zp["h1_targets"], zp["h1_dd"])}.get(t, np.nan)
                       for t in targets], dtype=np.float64)
    AGREE = Z.bottom_class(dd_h, targets, 0.30)
    rands = {f"C_rand{i+1}": Z.random_class(cl, int(C30.sum()), s)
             for i, s in enumerate((20260906, 20260907, 20260908))}
    P(f"  CLASS: {int(np.isfinite(cl).sum())}/{n} eligible; CLASS30 n = {int(C30.sum())} "
      f"(cut {np.nanmax(cl[C30]):+.5f}); of those CD4-carried {int((C30 & carried).sum())}; "
      f"AGREE30 {int(AGREE.sum())} (overlap {int((AGREE & C30).sum())}); TOPCONS30 {int(TOP30.sum())}")

    arms: dict[str, dict] = {}
    arrays: dict[str, np.ndarray] = {}
    audits: dict[str, dict] = {}

    def add(name: str, final: np.ndarray, **meta) -> dict:
        e = X3.score_h1(final, tr)
        for k in ("pds_a23", "cos_a23", "pds_a15", "cos_a15", "pds_eff", "cos_eff"):
            arrays[f"{name}|{k}"] = np.asarray(e[k], dtype=np.float64)
        row = {k: float(np.mean(v)) for k, v in e.items()}
        arms[name] = row
        audits[name] = dict(meta)
        P(f"  [{name:14s}] pds a23 {row['pds_a23']:.6f} eff {row['pds_eff']:.6f} | own cos "
          f"{row['cos_a23']:.5f} cos_agg {row['a23_cos_agg']:.5f} | delivered {row['a23_delivered_norm']:.4f}")
        return row

    add("BASE", base)
    add("cd4h1:off", X3.compose(base, [(cd4, fr_off)], cols), fraction=fr_off)
    add("INCUMBENT", X3.compose(base, [(d5, fr5)], cols), fraction=fr5, **aud5)

    # --------------------------------------------------------------------------- KA vs 273
    z3 = np.load(Z.E273 / "work/instrument273.npz", allow_pickle=True)
    assert [str(t) for t in z3["targets"]] == targets
    ka: list[dict] = []
    worst = 0.0
    for nm in ("BASE", "cd4h1:off", "INCUMBENT"):
        for k in ("pds_a23", "cos_a23", "pds_a15", "cos_a15", "pds_eff", "cos_eff"):
            dd_ = float(np.abs(np.asarray(z3[f"{nm}|{k}"]) - arrays[f"{nm}|{k}"]).max())
            worst = max(worst, dd_)
    ok = worst == 0.0
    ka.append(dict(row="KA-2 273's BASE / cd4h1:off / INCUMBENT per-target arrays", value=worst, ok=bool(ok)))
    P(f"  KA-2 re-composed 273 rows vs instrument273.npz: max|d| over 18 arrays {worst:.1e} "
      f"-> {'PASS' if ok else 'FAIL'}")
    if not ok:
        fail.append("KA-2 273 rows")

    # --------------------------------------------------------------------------- the renderings
    def cd4_arm(name: str, m, *, resolve: bool = True):
        fr = X3.solve_fraction(base, m, X3.CD4_F)[0] if resolve else fr5
        return add(name, X3.compose(base, [(m, fr)], cols), fraction=fr, resolved=resolve)

    cd4_arm("V1_undefl", Z.mix_rows(cd4, d5, C30, "V1"))
    cd4_arm("V1f_undefl", Z.mix_rows(cd4, d5, C30, "V1"), resolve=False)
    cd4_arm("V2_mirror", Z.mix_rows(d5, cd4, C30, "V2"))
    cd4_arm("V3_noresid", Z.zero_rows(d5, C30, "V3"))
    cd4_arm("V5_20", Z.mix_rows(cd4, d5, C20, "V5_20"))
    cd4_arm("V5_40", Z.mix_rows(cd4, d5, C40, "V5_40"))
    cd4_arm("C_top30", Z.mix_rows(cd4, d5, TOP30, "C_top30"))
    cd4_arm("C_agree30", Z.mix_rows(cd4, d5, AGREE, "C_agree30"))
    for nm, m in rands.items():
        cd4_arm(nm, Z.mix_rows(cd4, d5, m, nm))

    # --------------------------------------------------------------------------- per-class reads
    rng = np.random.default_rng(Z.SEED)
    base_cos = np.asarray(z3["base_cos_eff"], dtype=np.float64)
    grp = X3.groups_262(base_cos)
    tab: dict[str, dict] = {}
    inc = arrays["INCUMBENT|pds_a23"]
    for nm in arms:
        if nm in ("BASE", "cd4h1:off", "INCUMBENT"):
            continue
        d = arrays[f"{nm}|pds_a23"] - inc
        tab[nm] = dict(d_pds=float(d.mean()), boot=Z.boot(d, np.random.default_rng(Z.SEED), Z.B_BOOT),
                       d_pds_class=float(d[C30].mean()), d_pds_rest=float(d[~C30].mean()),
                       rel_cos_agg=float(arms[nm]["a23_cos_agg"] / arms["INCUMBENT"]["a23_cos_agg"] - 1.0),
                       lt3o=Z.lt3o(d),
                       by_group={g: float(d[m].mean()) for g, m in grp.items()})
    P()
    P(f"  {'arm':14s} {'d pds':>9s} {'[95% CI]':>22s} {'class':>9s} {'rest':>9s} {'cos_agg':>9s} {'LT3O':>7s}")
    for nm, r in tab.items():
        P(f"  {nm:14s} {r['d_pds']:+9.5f} [{r['boot']['lo']:+.5f}, {r['boot']['hi']:+.5f}] "
          f"{r['d_pds_class']:+9.5f} {r['d_pds_rest']:+9.5f} {100*r['rel_cos_agg']:+8.2f}% "
          f"{r['lt3o']['retention']:+7.2f}")

    arrays["targets"] = np.asarray(targets)
    arrays["c_lolo_h1"] = cl
    arrays["mask|CLASS30"] = C30
    arrays["mask|CLASS20"] = C20
    arrays["mask|CLASS40"] = C40
    arrays["mask|TOPCONS30"] = TOP30
    arrays["mask|AGREE30"] = AGREE
    arrays["cd4_carried"] = carried
    arrays["base_cos_eff"] = base_cos
    np.savez_compressed(Z.WORK / "flex277.npz", **arrays)
    (Z.WORK / "flex277.json").write_text(json.dumps(
        dict(meta=dict(n=n, n_class=int(C30.sum()), class_cut=float(np.nanmax(cl[C30])),
                       fr_off=fr_off, fr5=fr5, deflation=aud5, ka=ka, failures=fail,
                       seconds=round(time.time() - t0, 1), peak_rss_gib=Z.peak_gib()),
             arms=arms, audits=audits, table=tab, log=LOG), indent=1, default=float))
    P(f"done: {len(fail)} KA failures, {time.time()-t0:.0f}s, peak RSS {Z.peak_gib():.2f} GiB")
    if fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
