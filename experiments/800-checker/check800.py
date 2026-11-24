"""800 — the lab's checker: one command from a candidate file to a board-unit price with intervals, tries and a track record.

    uv run python experiments/800-checker/check800.py aggregate <run-dir> --cand <arm> [--base <arm>] [--candidate cand.json]
    uv run python experiments/800-checker/check800.py wear                  # configurations seen per surface, N per family
    uv run python experiments/800-checker/check800.py ka                    # the known-answer table (ka800.json -> markdown)

A run directory holds the engines' standard-format arrays (stats800's docstring):
    <run-dir>/bulk/<fold-id with / as __>/<arm>.npz            e.g. bulk/L3__orion_hct116/inc.npz
    <run-dir>/de/<surface>/<arm>/s<seed>.npz                    e.g. de/h1tv/inc/s0.npz
`aggregate` pairs cand with base on every fold that has both, prices bulk and DE folds separately (they measure different
members) and together, applies best-of-N from the wear ledger and the class calibration from ka800.json, decides the bin,
appends the evaluation to wear/ledger.jsonl, and writes <run-dir>/report-<cand>.json + .md.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import stats800 as S  # noqa: E402

WEAR = HERE / "wear" / "ledger.jsonl"
KA = HERE / "ka800.json"
ENGINE_FILES = ("experiments/800-checker/stats800.py", "experiments/800-checker/check800.py")


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def git_blob(rel: str) -> str:
    r = subprocess.run(["git", "hash-object", rel], cwd=ROOT, capture_output=True, text=True)
    return r.stdout.strip()


# ------------------------------------------------------------------------------------------------- candidate files
def load_candidate(path: str | None) -> dict:
    if not path:
        return {}
    p = ROOT / path if not Path(path).is_absolute() else Path(path)
    c = json.loads(p.read_text())
    for k in ("name", "family", "class"):
        assert k in c, f"candidate.json needs '{k}'"
    assert c["class"] in ("transfer", "scorer", "abc"), "class must be transfer | scorer | abc (AGENTS 11.1c.1)"
    shas = {"candidate.json": sha256(p)}
    if c.get("prereg"):
        got = sha256(ROOT / c["prereg"])
        assert got == c.get("prereg_sha256"), f"prereg sha256 mismatch: file {got} vs candidate {c.get('prereg_sha256')}"
        shas[c["prereg"]] = got
    for sec in ("bulk", "de"):
        mod = (c.get(sec) or {}).get("module")
        if mod:
            got = sha256(ROOT / mod)
            assert got == c[sec].get("sha256"), f"{sec} module sha256 mismatch: {mod}"
            shas[mod] = got
    c["_shas"] = shas
    return c


# ------------------------------------------------------------------------------------------------- wear ledger
def wear_append(rec: dict) -> None:
    WEAR.parent.mkdir(parents=True, exist_ok=True)
    with open(WEAR, "a") as fh:
        fh.write(json.dumps(rec, sort_keys=True) + "\n")


def wear_rows() -> list[dict]:
    return [json.loads(ln) for ln in WEAR.read_text().splitlines() if ln.strip()] if WEAR.exists() else []


def tries(family: str, surfaces: list[str], config_key: str | None = None) -> int:
    """N = distinct configurations of `family` priced on any of `surfaces` (this evaluation included)."""
    rows = [r for r in wear_rows() if r.get("family") == family]
    keys = {r["config_key"] for r in rows if set(r.get("surfaces", [])) & set(surfaces)}
    if config_key:
        keys.add(config_key)
    declared = max([int(r.get("declared_n", 0)) for r in rows] + [0])   # a family's declared prior reads (e.g. 810: N = 6)
    return max(1, len(keys), declared)


def cmd_wear(_):
    rows = wear_rows()
    surf, fam, decl = {}, {}, {}
    for r in rows:
        for s in r.get("surfaces", []):
            surf.setdefault(s, set()).add(r["config_key"])
        if r.get("declared_n"):
            decl[r.get("family", "?")] = max(decl.get(r.get("family", "?"), 0), int(r["declared_n"]))
            continue
        fam.setdefault(r.get("family", "?"), set()).add(r["config_key"])
    print(f"{len(rows)} evaluations in {WEAR.relative_to(ROOT)}")
    print("\nconfigurations seen per surface (wear):")
    for s in sorted(surf):
        print(f"  {s:28s} {len(surf[s]):4d}")
    print("\ntries per family (N for best-of-N):")
    for f in sorted(set(fam) | set(decl)):
        n = max(len(fam.get(f, ())), decl.get(f, 0))
        print(f"  {f:28s} {n:4d}  e(N) = {S.e_max_normal(n):.2f} sd" + (f"  (declared prior reads {decl[f]})" if f in decl else ""))


# ------------------------------------------------------------------------------------------------- loading a run
def load_folds(run: Path, base: str, cand: str, regime: str) -> list[S.Fold]:
    folds = []
    for d in sorted((run / "bulk").glob("*")) if (run / "bulk").exists() else []:
        b, c = d / f"{base}.npz", d / f"{cand}.npz"
        if b.exists() and c.exists():
            zb, zc = dict(np.load(b, allow_pickle=True)), dict(np.load(c, allow_pickle=True))
            assert list(zb["targets"]) == list(zc["targets"]), f"{d.name}: target lists differ"
            folds.append(S.Fold(d.name.replace("__", "/"), "bulk", zb, zc, meta=dict(change=zc.get("change_energy"))))
    for d in sorted((run / "de").glob("*")) if (run / "de").exists() else []:
        bs, cs = sorted((d / base).glob("s*.npz")), sorted((d / cand).glob("s*.npz"))
        seeds = sorted({p.stem for p in bs} & {p.stem for p in cs})
        if not seeds:
            continue
        B = [dict(np.load(d / base / f"{s}.npz", allow_pickle=True)) for s in seeds]
        C = [dict(np.load(d / cand / f"{s}.npz", allow_pickle=True)) for s in seeds]
        assert all(list(x["targets"]) == list(B[0]["targets"]) for x in B + C), f"{d.name}: target lists differ"
        folds.append(S.Fold(d.name, "de", B, C, regime=regime, meta=dict(seeds=seeds)))
    return folds


def config_key(cand: dict, cand_arm: str, base_arm: str) -> str:
    if cand.get("_shas"):
        return hashlib.sha256(json.dumps(cand["_shas"], sort_keys=True).encode()).hexdigest()[:16]
    return f"{cand_arm}-vs-{base_arm}"


def cmd_aggregate(a):
    run = Path(a.run)
    cand = load_candidate(a.candidate)
    folds = load_folds(run, a.base, a.cand, a.regime)
    assert folds, f"no fold in {run} carries both {a.base} and {a.cand}"
    family = cand.get("family", a.family or a.cand)
    cls = cand.get("class", a.cls or "?")
    calib_class = a.calib_class or cand.get("calib_class")
    out = dict(utc=utc(), run=str(run.resolve().relative_to(ROOT)) if run.resolve().is_relative_to(ROOT) else str(run), cand=a.cand, base=a.base,
               family=family, cls=cls, regime=a.regime, candidate=cand.get("name"), engine={f: git_blob(f) for f in ENGINE_FILES})
    blocks = {}
    # surfaces are read SEPARATELY (804 PROTOCOL §5): A/B/C-panel folds (ABC/*, ABCp/*) and fresh folds (L*/...) are different
    # instruments with different calibrations; a pooled bulk number is reported but flagged as not a decision number
    # (810, 10-16: a pool of 18 fresh + 3 A/B/C folds read SHIP while A/B/C alone read −0.0020)
    def surface_set(fid):
        # LABC (A/B/C's targets on the fresh-fold recipe) is an Arc-drawn panel, not a fresh draw: its own set (10-16 11:10Z)
        if fid.startswith("LABC/"):
            return "BULK-LABC"
        return "BULK-ABC" if fid.startswith("ABC/") else ("BULK-ABCp" if fid.startswith("ABCp/") else "BULK-W")
    bulk = [f for f in folds if f.kind == "bulk"]
    for ss in sorted({surface_set(f.fid) for f in bulk}):
        F = [f for f in bulk if surface_set(f.fid) == ss]
        blocks[ss] = S.price(F, members=("pds", "mse"), n_boot=a.boot, seed=800)
    if len({surface_set(f.fid) for f in bulk}) > 1:
        blocks["bulk-pooled (NOT a decision number)"] = S.price(bulk, members=("pds", "mse"), n_boot=a.boot, seed=800)
    de_members = tuple(m for m in ("fid", "jac", "nmae", "reach") if m in a.de_members)
    F = [f for f in folds if f.kind == "de"]
    if F:
        blocks["DE-H1"] = S.price(F, members=de_members, n_boot=a.boot, seed=800)
    if "DE-H1" in blocks:
        for ss in [k for k in blocks if k.startswith("BULK-")]:
            c = S.combine([blocks[ss], blocks["DE-H1"]])
            c.update(per_fold=blocks[ss]["per_fold"] + blocks["DE-H1"]["per_fold"], members={},
                     folds_pos=blocks[ss]["folds_pos"] + blocks["DE-H1"]["folds_pos"],
                     n_folds=blocks[ss]["n_folds"] + blocks["DE-H1"]["n_folds"], leave_top_out_mean=None,
                     overall_members=blocks[ss]["overall_members"] + blocks["DE-H1"]["overall_members"])
            blocks[f"combined {ss} + DE-H1"] = c
    surfaces = [f.fid for f in folds]
    ck = config_key(cand, a.cand, a.base)
    n = tries(family, surfaces, ck)
    cal = json.loads(KA.read_text()).get("calibration", {}) if KA.exists() else {}
    def calib_for(kind):
        # each block is projected with ITS surface's calibration (10-16: ABCp has its own rho; a W read must not be projected
        # with an A/B/C rho). calib_class may carry a surface suffix ("bulk-key|ABC"); only its class part is used for bulk
        # blocks. DE-H1 takes a "de-…" class as given. Pooled and combined blocks get no projection (no validated rho).
        if not calib_class:
            return None
        if kind.startswith("BULK-"):
            return cal.get(f"{calib_class.split('|')[0]}|{kind[len('BULK-'):]}")
        if kind == "DE-H1" and calib_class.startswith("de-"):
            return cal.get(calib_class)
        return None
    for kind, res in blocks.items():
        res["best_of_n"] = S.best_of_n(res, n)
        res["projection"] = S.project(dict(point=res["best_of_n"]["point"], ci95=res["best_of_n"]["ci95"]), calib_for(kind))
        lo, hi = res["best_of_n"]["ci95"]
        res["bin_checker"] = S.decide(lo, res["best_of_n"]["point"], hi)
        if res["projection"].get("ci95"):
            pl, ph = res["projection"]["ci95"]
            res["bin_board"] = S.decide(pl, res["projection"]["point"], ph)
    out["blocks"] = blocks
    out["n_tries"] = n
    out["calib_class"] = calib_class
    rpt = run / f"report-{a.cand}-vs-{a.base}.json"
    rpt.write_text(json.dumps(out, indent=1, default=lambda v: v.tolist() if isinstance(v, np.ndarray) else float(v)))
    if not a.no_wear:
        wear_append(dict(utc=out["utc"], candidate=cand.get("name", a.cand), family=family, cls=cls, config_key=ck,
                         surfaces=surfaces, base=a.base, cand=a.cand, run=out["run"],
                         summary={k: dict(point=v["point"], ci95=v["ci95"], bin=v["bin_checker"]) for k, v in blocks.items()},
                         engine=out["engine"], shas=cand.get("_shas", {})))
    print(render(out))


def render(o: dict) -> str:
    L = [f"## {o['cand']} vs {o['base']}  (family {o['family']}, class {o['cls']}, tries N={o['n_tries']}, regime {o['regime']})"]
    for kind, r in o["blocks"].items():
        L.append(f"\n**{kind}**: overall {r['point']:+.5f} [{r['ci95'][0]:+.5f}, {r['ci95'][1]:+.5f}] board units, "
                 f"P(<=0) {r['p_le0']:.3f}, folds + {r['folds_pos']}/{r['n_folds']}, leave-top-3-out "
                 f"{r['leave_top_out_mean'] if r['leave_top_out_mean'] is None else format(r['leave_top_out_mean'], '+.5f')}; "
                 f"members {r['overall_members']}")
        b = r["best_of_n"]
        L.append(f"- best-of-{b['n_tries']}: shift {b['shift']:.5f} -> {b['point']:+.5f} [{b['ci95'][0]:+.5f}, {b['ci95'][1]:+.5f}]"
                 f" -> checker bin **{r['bin_checker']}**")
        p = r["projection"]
        L.append(f"- board projection: " + (f"{p['point']:+.5f} [{p['ci95'][0]:+.5f}, {p['ci95'][1]:+.5f}] (rho {p['rho']:.2f})"
                                            f" -> **{r.get('bin_board')}**" if p.get("ci95") else p["note"]))
        if not r["members"]:
            continue
        L.append("| member | raw Δ | 95 % CI | overall |")
        L.append("|---|---:|---|---:|")
        for m, v in r["members"].items():
            ov = "" if v["overall"] is None else f"{v['overall']:+.5f}"
            L.append(f"| {m} | {v['raw']:+.6f} | [{v['ci95'][0]:+.6f}, {v['ci95'][1]:+.6f}] | {ov} |")
        L.append("| fold | n | overall | " + " | ".join(k for k in r["per_fold"][0]["raw"]) + " |")
        L.append("|---|---:|---:|" + "---:|" * len(r["per_fold"][0]["raw"]))
        for pf in r["per_fold"]:
            L.append(f"| {pf['fold']} | {pf['n']} | {pf['overall']:+.5f} | " + " | ".join(f"{v:+.5f}" for v in pf["raw"].values()) + " |")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    g = sp.add_parser("aggregate")
    g.add_argument("run")
    g.add_argument("--cand", required=True)
    g.add_argument("--base", default="inc")
    g.add_argument("--candidate", default=None, help="candidate.json (sha-checked); omit for named known-answer arms")
    g.add_argument("--family", default=None)
    g.add_argument("--cls", default=None)
    g.add_argument("--calib-class", default=None)
    g.add_argument("--regime", default="f035")
    g.add_argument("--de-members", default="fid,jac,nmae", help="reach is unpriced for delivery levers (619)")
    g.add_argument("--boot", type=int, default=4000)
    g.add_argument("--no-wear", action="store_true", help="known-answer reruns of the same config only")
    g.set_defaults(fn=cmd_aggregate)
    w = sp.add_parser("wear")
    w.set_defaults(fn=cmd_wear)
    a = ap.parse_args()
    if a.cmd == "aggregate":
        a.de_members = tuple(a.de_members.split(","))
    a.fn(a)


if __name__ == "__main__":
    main()
