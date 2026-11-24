"""810 -- PROPOSAL-2.md's decision rule, as code, written before any price exists. Light (MBs, local).

Reads price810's STANDARD-FORMAT arrays (<run>/bulk/<fold>/{inc,707,sgu1,sgu1all}.npz) and, with stats800's own per-fold
arithmetic (bulk_fold_raw: r, literal pds, pdsl, mse via 518's identity; overall_of: board units), computes on each surface
set (BULK-ABC cancer = ABC/k562_gwps, ABC/orion_hct116, ABC/orion_hek293t; BULK-W = the 18 L folds; ABC/kolf2 reported):
  * per arm vs inc: overall = mse leg + pds leg (literal pds DECIDES; the pdsl version is reported beside), point, paired
    target-bootstrap 95 % interval (4,000 draws, seed 810; every fold resampled independently, equal weight per fold, the
    SAME resampled targets for every arm in a draw), folds positive;
  * the exchange E = -(pds leg) / (mse leg) (the share of the mse gain the pds leg eats), per arm, point and P(E_sgu1 <
    E_sgu1all) over the paired draws;
  * paired differences sgu1 - sgu1all and sgu1 - 707 (overall);
  * the best-of-N shift on sgu1: N = 6 (PROPOSAL-2.md section 3), e(6) = stats800.e_max_normal(6) x the bootstrap sd.
Then it prints the verdict by PROPOSAL-2.md section 6's rule (STOP / CLOSE / HOLD-707 / DOSSIER), and writes
<run>/decide810.json. It never edits the wear ledger (the run's check800 aggregate does that).

    uv run python experiments/810-strong-direction/decide810.py <run-dir> [--max-dk 0.05]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/800-checker"))
import stats800 as S  # noqa: E402

ARMS = ("707", "sgu1", "sgu1all")
N_TRIES = 6
SETS = {"ABC": ("ABC__k562_gwps", "ABC__orion_hct116", "ABC__orion_hek293t"), "kolf2": ("ABC__kolf2",)}
B, SEED = 4000, 810


def legs(raw, pds_key="pds"):
    mse = S.overall_of({"mse": raw["mse"]})
    pds = S.overall_of({"pds": raw[pds_key]})
    return mse, pds


def load(run: Path, fold: str):
    d = run / "bulk" / fold
    return {a: dict(np.load(d / f"{a}.npz", allow_pickle=True)) for a in ("inc",) + ARMS}


def evaluate(run: Path, folds: list[str]) -> dict:
    Z = {f: load(run, f) for f in folds}
    n = {f: Z[f]["inc"]["d"].size for f in folds}

    def stat(idx_of):
        out = {}
        for a in ARMS:
            m, p, pl, r = [], [], [], []
            for f in folds:
                raw = S.bulk_fold_raw(Z[f]["inc"], Z[f][a], idx_of(f))
                mm, pp = legs(raw)
                _, ppl = legs(raw, "pdsl")
                m.append(mm); p.append(pp); pl.append(ppl); r.append(raw["r"])
            out[a] = dict(mse=float(np.mean(m)), pds=float(np.mean(p)), pdsl=float(np.mean(pl)), r=float(np.mean(r)),
                          overall=float(np.mean(m) + np.mean(p)), overall_pdsl=float(np.mean(m) + np.mean(pl)),
                          per_fold=[mi + pi for mi, pi in zip(m, p)])
        return out

    point = stat(lambda f: None)
    rng = np.random.default_rng(SEED)
    draws = {a: {k: np.empty(B) for k in ("overall", "overall_pdsl", "E")} for a in ARMS}
    diff = {k: np.empty(B) for k in ("sgu1-sgu1all", "sgu1-707")}
    lessE = np.empty(B, bool)
    for b in range(B):
        idx = {f: rng.integers(0, n[f], n[f]) for f in folds}
        st = stat(lambda f: idx[f])
        for a in ARMS:
            draws[a]["overall"][b] = st[a]["overall"]
            draws[a]["overall_pdsl"][b] = st[a]["overall_pdsl"]
            draws[a]["E"][b] = -st[a]["pds"] / st[a]["mse"] if st[a]["mse"] > 0 else np.nan
        diff["sgu1-sgu1all"][b] = st["sgu1"]["overall"] - st["sgu1all"]["overall"]
        diff["sgu1-707"][b] = st["sgu1"]["overall"] - st["707"]["overall"]
        lessE[b] = draws["sgu1"]["E"][b] < draws["sgu1all"]["E"][b]
    q = lambda v: [float(np.nanquantile(v, .025)), float(np.nanquantile(v, .975))]
    res = {}
    for a in ARMS:
        P = point[a]
        res[a] = dict(**{k: P[k] for k in ("mse", "pds", "pdsl", "r", "overall", "overall_pdsl")},
                      ci95=q(draws[a]["overall"]), sd=float(np.std(draws[a]["overall"])), ci95_pdsl=q(draws[a]["overall_pdsl"]),
                      E=(-P["pds"] / P["mse"]) if P["mse"] > 0 else None, folds_pos=int(sum(v > 0 for v in P["per_fold"])),
                      n_folds=len(folds))
    s = S.e_max_normal(N_TRIES) * res["sgu1"]["sd"]
    res["sgu1"]["best_of_n"] = dict(n=N_TRIES, shift=s, point=res["sgu1"]["overall"] - s, lo=res["sgu1"]["ci95"][0] - s)
    res["paired"] = {k: dict(point=(point["sgu1"]["overall"] - point[k.split("-")[1]]["overall"]), ci95=q(v))
                     for k, v in diff.items()}
    res["P(E_sgu1 < E_sgu1all)"] = float(np.mean(lessE))
    return res


def dk_ok(run: Path, max_dk: float) -> tuple[bool, float]:
    worst = 0.0
    for p in sorted((run / "meta").glob("*.json")):
        m = json.loads(p.read_text())
        for a in ARMS:
            v = (m["arms"][a].get("vs_inc") or {}).get("dK_rel")
            if v is not None:
                worst = max(worst, abs(v))
    return worst <= max_dk, worst


def verdict(R: dict, dk: tuple) -> tuple[str, list[str]]:
    """PROPOSAL-2.md section 6, in order."""
    A, W = R["ABC"], R["W"]
    why = []
    if A["sgu1all"]["ci95"][0] > 0:
        return "STOP", ["sgu1all (the closed global form) reads clearly net-positive on BULK-ABC (lower bound > 0): the checker "
                        "does not see the pds mechanism that closed gamma; the gated read is unpriceable"]
    if not dk[0]:
        return "DE-H1 FIRST", [f"|dK_rel| {dk[1]:.3f} > 5 %: the bulk engine cannot price it"]
    E1, Eall = A["sgu1"]["E"], A["sgu1all"]["E"]
    close = []
    if A["sgu1"]["overall"] <= 0:
        close.append("C1: sgu1 BULK-ABC point <= 0")
    if E1 is None or Eall is None or E1 >= Eall:
        close.append(f"C2 (falsifier 2): E(sgu1) {E1} >= E(sgu1all) {Eall}: the gate does not concentrate mse relative to pds")
    if A["paired"]["sgu1-sgu1all"]["point"] <= 0:
        close.append("C3: sgu1 - sgu1all paired point <= 0 on BULK-ABC")
    if A["paired"]["sgu1-707"]["point"] <= 0:
        close.append("C4: sgu1 - 707 paired point <= 0 on BULK-ABC (adds nothing over 707)")
    if W["sgu1"]["overall"] <= 0:
        close.append("C5: sgu1 BULK-W point <= 0")
    if close:
        return "CLOSE", close
    ship = (A["sgu1"]["best_of_n"]["lo"] > 0 and A["sgu1"]["overall_pdsl"] > 0 and W["sgu1"]["folds_pos"] >= 14)
    if ship:
        return "DOSSIER", ["all of D1-D4 hold; the board projection stays a range including 0 until 707's board read"]
    why.append(f"sgu1 BULK-ABC lower bound after e(6) {A['sgu1']['best_of_n']['lo']:+.5f} <= 0, or pdsl net "
               f"{A['sgu1']['overall_pdsl']:+.5f} <= 0, or BULK-W folds+ {W['sgu1']['folds_pos']}/18 < 14")
    return "HOLD-707", why


def main():
    run = Path(sys.argv[1])
    max_dk = float(sys.argv[sys.argv.index("--max-dk") + 1]) if "--max-dk" in sys.argv else 0.05
    folds = sorted(p.name for p in (run / "bulk").iterdir() if p.is_dir())
    sets = dict(SETS, W=tuple(f for f in folds if f.startswith("L")))
    assert len(sets["W"]) == 18 and all(f in folds for f in SETS["ABC"]), "need all 18 L folds and the 3 ABC cancer folds"
    R = {k: evaluate(run, list(v)) for k, v in sets.items() if all(f in folds for f in v)}
    dk = dk_ok(run, max_dk)
    v, why = verdict(R, dk)
    out = dict(results=R, dk_worst=dk[1], verdict=v, why=why, n_tries=N_TRIES)
    (run / "decide810.json").write_text(json.dumps(out, indent=1, default=float))
    for k in ("ABC", "W", "kolf2"):
        if k in R:
            for a in ARMS:
                x = R[k][a]
                print(f"{k:5s} {a:8s} overall {x['overall']:+.5f} [{x['ci95'][0]:+.5f}, {x['ci95'][1]:+.5f}] (pdsl {x['overall_pdsl']:+.5f}) "
                      f"mse {x['mse']:+.5f} pds {x['pds']:+.5f} r {100 * x['r']:+.2f}% E {x['E']} folds+ {x['folds_pos']}/{x['n_folds']}")
            print(f"{k:5s} paired {json.dumps(R[k]['paired'])}  P(E_sgu1 < E_sgu1all) {R[k]['P(E_sgu1 < E_sgu1all)']:.3f}")
    print(f"dK worst {dk[1]:.4f}\nVERDICT: {v}\n  " + "\n  ".join(why))


if __name__ == "__main__":
    main()
