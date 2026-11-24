"""308 -- the six-member price of every |kd| arm OVER cd4svd5, on the three COHERENT scenarios (READ-FIRST section 3;
261's recipe imported verbatim through 271's path: `price261.scenarios` / `price261.retention` with 257's constants).
The RETRODICTION row runs first on THIS experiment's own `cd4svd5` and `STACK` rows: it must reproduce BR249's board
step (+0.0043) and 257's R1 to 1e-6, or nothing is priced.

Gate: AGENTS section 8.4 / READ-FIRST section 3 -- amplitude-changing arms +0.0089, shape-only +0.0019.
    PYTHONDONTWRITEBYTECODE=1 uv run python price308.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E257 = ROOT / "experiments/257-soft-deflation-price"
E261 = ROOT / "experiments/261-union-deflation"
sys.path.insert(0, str(E261))
sys.path.insert(0, str(E257))
import price261 as PR  # noqa: E402

P = PR.P
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
BASE = "cd4svd5"
LINES: list[str] = []


def Pr(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def load():
    rows, pt, idn, meta = {}, {}, {}, {}
    for held in CANCER:
        d = json.loads((WORK / f"sweep308-{held}.json").read_text())
        meta[held] = d["meta"]
        for r in d["rows"]:
            rows[(held, r["arm"])] = r
        s = json.loads((WORK / f"stack308-{held}.json").read_text())
        for r in s["rows"]:
            if r["arm"] == "STACK":
                rows[(held, "STACK")] = r
        z = np.load(WORK / f"pertarget308-{held}.npz", allow_pickle=True)
        for k in z.files:
            if k.endswith("|pds") or k.endswith("|cos"):
                pt[(held, k)] = np.asarray(z[k], dtype=np.float64)
        z = np.load(WORK / f"identity308-{held}.npz")
        for k in z.files:
            idn[(held, k)] = z[k]
        zs = np.load(WORK / f"stackid308-{held}-STACK.npz")
        for k in zs.files:
            idn[(held, f"STACK|{k}")] = zs[k]
    return rows, pt, idn, meta


def nmae_concat(idn, arm, folds) -> float:
    U = np.concatenate([idn[(f, f"{arm}|U")] for f in folds])
    M = np.concatenate([idn[(f, f"{arm}|M")] for f in folds])
    D = np.concatenate([idn[(f, f"{arm}|D")] for f in folds])
    return float(np.mean(1 + (U - 2 * M) / D))


def compare(rows, idn, arm, base, folds) -> dict:
    """271's `compare` verbatim (mean over folds of per-fold ratios / deltas)."""
    ra, rm, dp, dq, dn_a, dn_b, fr = [], [], [], [], [], [], []
    for f in folds:
        a, b = rows[(f, arm)], rows[(f, base)]
        ra.append(a["cos_agg"] / b["cos_agg"] - 1.0)
        rm.append(a["a23_cos"] / b["a23_cos"] - 1.0)
        dp.append(a["pds"] - b["pds"])
        dq.append(a["a15_sgn_w"] - b["a15_sgn_w"])
        dn_a.append(a["delivered_norm"])
        dn_b.append(b["delivered_norm"])
        fr.append(a["audit"]["cd4_fraction"] / b["audit"]["cd4_fraction"])
    nm_a, nm_b = nmae_concat(idn, arm, folds), nmae_concat(idn, base, folds)
    return dict(arm=arm, base=base, n=len(folds), rel_agg=float(np.mean(ra)), rel_macro=float(np.mean(rm)),
                d_pds=float(np.mean(dp)), per_fold_d_pds=[float(x) for x in dp],
                pds_pos=int(sum(x > 0 for x in dp)), sigma=float(np.mean(dp)) / P.SIGMA_PDS,
                d_qw=float(np.mean(dq)), delivered_ratio=float(np.mean(dn_a) / np.mean(dn_b)),
                fraction_ratio=float(np.mean(fr)), nmae_arm=nm_a, nmae_base=nm_b, d_nmae=nm_a - nm_b,
                d_prec89=float(np.mean([rows[(f, arm)]["top89_precision"] - rows[(f, base)]["top89_precision"]
                                        for f in folds])), energy_removed=0.0)


def main() -> None:
    rows, pt, idn, meta = load()
    an = json.loads((WORK / "analyse308.json").read_text())
    targets = {f: [str(t) for t in np.load(WORK / f"pertarget308-{f}.npz", allow_pickle=True)["targets"]]
               for f in CANCER}
    arms = [r["arm"] for r in json.loads((WORK / "sweep308-k562_gwps.json").read_text())["rows"] if r["arm"] != BASE]

    Pr("=" * 132)
    Pr("308 -- the |kd|-weighted-donor arms priced OVER cd4svd5 (261's recipe, three COHERENT scenarios); retrodiction first")
    Pr("=" * 132)
    for f in CANCER:
        m = meta[f]
        Pr(f"  {f:14s} {m['n_targets']} targets, donors {m['donors']}, KA failures {m['failures']}, "
           f"peak {m['peak_gib']:.2f} GiB")

    # ---- the retrodiction row (built inside THIS experiment: stack308.py's STACK + the sweep's cd4svd5)
    c_r = compare(rows, idn, BASE, "STACK", CANCER)
    s_r = PR.p257_sets(c_r, P.COS_H1UNION, P.PER_H1UNION)
    r2 = s_r["mine"]["point"]
    p257 = json.loads((E257 / "work/price257.json").read_text())
    ref_r1 = p257["retrodiction"]["R1"]["price"]["point"]["total"]
    c257 = p257["retrodiction"]["R1"]["compare"]
    d_cmp = max(abs(c_r["d_pds"] - c257["d_pds"]), abs(c_r["rel_agg"] - c257["rel_agg"]),
                abs(c_r["d_nmae"] - c257["d_nmae"]))
    board = P.BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_avg"] - P.BOARD["iUziqE1gewfcPKzYUDaB"]["score_avg"]
    ok_r = abs(r2["total"] - ref_r1) < 1e-6 and d_cmp < 1e-6 and s_r["max_abs_diff"] < 1e-12
    Pr(f"\nRETRODICTION  cd4svd5 vs STACK: d pds {c_r['d_pds']:+.6f} (257 {c257['d_pds']:+.6f}); agg "
       f"{100*c_r['rel_agg']:+.4f}% (257 {100*c257['rel_agg']:+.4f}%); d nmae89 {c_r['d_nmae']:+.5f} "
       f"(257 {c257['d_nmae']:+.5f}); max|d| {d_cmp:.1e}")
    Pr(f"              priced POINT: pds {r2['pds']:+.5f} mse {r2['mse']:+.5f} DE3 {r2['de3']:+.5f} "
       f"nmae {r2['nmae']:+.5f} = {r2['total']:+.5f}  (257 R1 {ref_r1:+.5f}; BR249 board {board:+.7f})"
       f"  -> {'PASS' if ok_r else 'FAIL'}")
    out: dict = dict(retrodiction=dict(compare=c_r, price=r2, ref257=ref_r1, board=board, ok=ok_r),
                     gate_amplitude=P.GATE_AMPL, gate_shape=P.GATE_SHAPE)
    if not ok_r:
        Pr("*** RETRODICTION FAILED -- nothing is priced")
        (WORK / "price308.json").write_text(json.dumps(out, indent=1, default=float))
        (WORK / "price308.log").write_text("\n".join(LINES) + "\n")
        raise SystemExit(1)

    Pr("\n" + "=" * 132)
    Pr("EVERY CELL as an increment over cd4svd5 (3 cancer folds); legs at POINT; PESSIMISTIC / POINT / OPTIMISTIC totals")
    Pr(f"gate: amplitude-changing >= {P.GATE_AMPL:+.4f}; shape-only >= {P.GATE_SHAPE:+.4f} (an arm is shape-only if the "
       f"delivered-norm ratio is within {100*P.AMP_SHAPE_BOUND:.1f} %)")
    Pr("=" * 132)
    Pr(f"{'cell':26s}{'d pds':>10s}{'f':>4s}{'boot 95% CI':>24s}{'agg%':>8s}{'dlv%':>7s}{'d nmae':>9s}"
       f"{'pds':>9s}{'mse':>9s}{'DE3':>9s}{'nmae':>9s}{'PESS':>9s}{'POINT':>9s}{'OPT':>9s}{'ret3%':>7s}{'shape':>6s}")
    table = {}
    for name in arms:
        c = compare(rows, idn, name, BASE, CANCER)
        sc = PR.scenarios(c)
        ret = PR.retention(pt, name, BASE, CANCER, targets)
        b = an["arms"][name]["pooled"]
        shape = abs(c["delivered_ratio"] - 1) <= P.AMP_SHAPE_BOUND
        gate = P.GATE_SHAPE if shape else P.GATE_AMPL
        table[name] = dict(compare=c, scenarios=sc, retention=ret, boot=b, shape_only=bool(shape), gate=gate,
                           cleared=bool(sc["PESSIMISTIC"]["total"] >= gate and c["pds_pos"] == 3 and b["lo"] > 0
                                        and ret["retained"] >= 0.5))
        pe, po, op = sc["PESSIMISTIC"], sc["POINT"], sc["OPTIMISTIC"]
        Pr(f"{name:26s}{c['d_pds']:+10.6f}{c['pds_pos']:>3d}/3 [{b['lo']:+.4f}, {b['hi']:+.4f}]"
           f"{100*c['rel_agg']:+8.3f}{100*(c['delivered_ratio']-1):+7.3f}{c['d_nmae']:+9.5f}"
           f"{po['pds']:+9.5f}{po['mse']:+9.5f}{po['de3']:+9.5f}{po['nmae']:+9.5f}"
           f"{pe['total']:+9.5f}{po['total']:+9.5f}{op['total']:+9.5f}{100*ret['retained']:+7.0f}"
           f"{('yes' if shape else 'no'):>6s}")
    out["cells"] = {k: dict(d_pds=v["compare"]["d_pds"], pds_pos=v["compare"]["pds_pos"],
                            d_nmae=v["compare"]["d_nmae"], rel_agg=v["compare"]["rel_agg"],
                            delivered_ratio=v["compare"]["delivered_ratio"], boot=v["boot"],
                            pessimistic=v["scenarios"]["PESSIMISTIC"]["total"],
                            point=v["scenarios"]["POINT"]["total"], optimistic=v["scenarios"]["OPTIMISTIC"]["total"],
                            retention=v["retention"]["retained"], shape_only=v["shape_only"], gate=v["gate"],
                            cleared=v["cleared"]) for k, v in table.items()}
    ARMS = [a for a in arms if (a.startswith("gate:") or a.startswith("soft:"))
            and ":rev" not in a and ":shuf" not in a and ":rand" not in a and ":lit" not in a and a != "soft:rel:p0"]
    best = max(ARMS, key=lambda n: table[n]["scenarios"]["PESSIMISTIC"]["total"])
    cleared = [n for n in ARMS if table[n]["cleared"]]
    Pr("\n" + "=" * 132)
    Pr(f"ARMS screened: {len(ARMS)}; CLEARED: {len(cleared)} {cleared}")
    t = table[best]
    Pr(f"best arm by pessimistic total: {best}  PESS {t['scenarios']['PESSIMISTIC']['total']:+.5f} "
       f"= {t['scenarios']['PESSIMISTIC']['total']/t['gate']:.2f} x its gate ({t['gate']:+.4f}, "
       f"{'shape-only' if t['shape_only'] else 'amplitude-changing'}); POINT "
       f"{t['scenarios']['POINT']['total']:+.5f}; OPT {t['scenarios']['OPTIMISTIC']['total']:+.5f}; "
       f"d pds {t['compare']['d_pds']:+.5f} [{t['boot']['lo']:+.5f}, {t['boot']['hi']:+.5f}] on "
       f"{t['compare']['pds_pos']}/3 folds; d nmae {t['compare']['d_nmae']:+.5f}")
    out["best_arm"] = best
    out["cleared"] = cleared
    out["n_arms_screened"] = len(ARMS)
    out["n_cells_screened"] = len(arms)
    (WORK / "price308.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "price308.log").write_text("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
