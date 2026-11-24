"""811 -- aggregate a run's STANDARD-FORMAT arrays with the checker's own statistics (stats800: paired target bootstrap,
equal weight per fold, best-of-N, class calibration, ship rule), WITHOUT check800.py's side effects (it appends to 800's wear
ledger and calls git; this worker writes only in its own directory). The would-be wear rows go to work/wear811.jsonl for the
lab lead to append.

Surfaces: W = the 18 L-folds, ABC = the three A/B/C cancer folds, ABCk = ABC/kolf2 (reported apart, PHASEB).
Every report carries the literal-pds overall (the checker's standard) and the board-emulated variant (pdsl, lambda 0.5).

    uv run python experiments/811-lolo-model/agg811.py <run-dir> --family F --tries N --pairs cand:base[,cand:base...]
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/800-checker"))
import stats800 as S  # noqa: E402

KA = HERE / "work" / "ka800.copy.json"     # byte copy (cp -c) of experiments/800-checker/ka800.json; sha256 in RESULTS.md
SURF = {"W": lambda f: f.startswith("L"), "ABC": lambda f: f.startswith("ABC") and not f.endswith("kolf2"),
        "ABCk": lambda f: f == "ABC/kolf2"}


def load(run: Path, base: str, cand: str, pick, pdsl=False):
    out = []
    for d in sorted((run / "bulk").glob("*")):
        fid = d.name.replace("__", "/")
        if not pick(fid) or not (d / f"{base}.npz").exists() or not (d / f"{cand}.npz").exists():
            continue
        zb, zc = dict(np.load(d / f"{base}.npz", allow_pickle=True)), dict(np.load(d / f"{cand}.npz", allow_pickle=True))
        assert list(zb["targets"]) == list(zc["targets"])
        if pdsl:
            zb, zc = {**zb, "pds": zb["pdsl"]}, {**zc, "pds": zc["pdsl"]}
        out.append(S.Fold(fid, "bulk", zb, zc))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--family", required=True)
    ap.add_argument("--tries", type=int, required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--calib", default="bulk-external")
    ap.add_argument("--boot", type=int, default=4000)
    a = ap.parse_args()
    run = Path(a.run)
    cal = json.loads(KA.read_text()).get("calibration", {})
    res, wear = {}, []
    for pair in a.pairs.split(","):
        cand, base = pair.split(":")
        for sname, pick in SURF.items():
            for variant in ("pds", "pdsl"):
                F = load(run, base, cand, pick, pdsl=(variant == "pdsl"))
                if not F:
                    continue
                r = S.price(F, members=("pds", "mse"), n_boot=a.boot, seed=800)
                bon = S.best_of_n(r, a.tries)
                c = cal.get(f"{a.calib}|{'W' if sname == 'W' else 'ABC'}")
                proj = S.project(dict(point=bon["point"], ci95=bon["ci95"]), c)
                rec = dict(cand=cand, base=base, surface=sname, pds_variant=variant, n_folds=r["n_folds"],
                           point=r["point"], ci95=r["ci95"], sd=r["sd"], p_le0=r["p_le0"], folds_pos=r["folds_pos"],
                           leave_top3_out=r["leave_top_out_mean"], members=r["members"],
                           best_of_n=bon, bin_checker=S.decide(bon["ci95"][0], bon["point"], bon["ci95"][1]),
                           projection=proj, bin_board=(S.decide(proj["ci95"][0], proj["point"], proj["ci95"][1])
                                                       if proj.get("ci95") else None),
                           per_fold=[dict(fold=p["fold"], n=p["n"], overall=p["overall"], **{k: float(v) for k, v in p["raw"].items()})
                                     for p in r["per_fold"]])
                res[f"{cand}|{base}|{sname}|{variant}"] = rec
                if variant == "pds":
                    wear.append(dict(utc=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), candidate=cand,
                                     family=a.family, cls="transfer", base=base, cand=cand, run=str(run),
                                     surfaces=[f.fid for f in F], config_key=f"811:{cand}-vs-{base}",
                                     summary=dict(bulk=dict(point=r["point"], ci95=r["ci95"], bin=rec["bin_checker"]))))
                m = r["members"]
                print(f"{cand:11s} vs {base:8s} {sname:4s} {variant:4s} n={r['n_folds']:2d} overall {r['point']:+.5f} "
                      f"[{r['ci95'][0]:+.5f}, {r['ci95'][1]:+.5f}] P<=0 {r['p_le0']:.3f} folds+ {r['folds_pos']}/{r['n_folds']} | "
                      f"best-of-{a.tries} {bon['point']:+.5f} [{bon['ci95'][0]:+.5f}, {bon['ci95'][1]:+.5f}] {rec['bin_checker']:4s} | "
                      f"r {100 * m['r']['raw']:+.2f}% pds {m['pds']['raw']:+.5f} mse {m['mse']['raw']:+.5f} | board proj "
                      + (f"{proj['point']:+.5f} [{proj['ci95'][0]:+.5f}, {proj['ci95'][1]:+.5f}] {rec['bin_board']}" if proj.get('ci95') else "n/a"),
                      flush=True)
    tag = hashlib.sha256(a.pairs.encode()).hexdigest()[:8]
    (run / f"agg811-{a.family}-{tag}.json").write_text(json.dumps(dict(args=vars(a), results=res), indent=1, default=float))
    with open(HERE / "work" / "wear811.jsonl", "a") as fh:
        for w in wear:
            fh.write(json.dumps(w, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
