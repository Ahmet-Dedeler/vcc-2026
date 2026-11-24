#!/usr/bin/env python3
"""803 — the lab's DE engine: the production emitter on the H1 real-emitter surfaces, one command per (surface, arms, seed).

  uv run python experiments/803-de-engine/de803.py run --surface {h1test,h1tv} --arms A[,B,...] --seed S
      Every arm named must share ONE emitter configuration ("emit" in arms/<arm>.json); the emitter runs once (emit803.py,
      which imports 618/run618e.py or 612/run612.py and 496/run496.py verbatim), then 630/rnup630.py once for the union of
      the arms' rnup variants, then the members step for every arm. Per arm it writes
          work/out/<surface>/<arm>/s<S>/members.npz   per-target arrays (below)
          work/out/<surface>/<arm>/s<S>/manifest.json  argv, arm file + sha, every step's manifest (git blobs of imported
                                                       modules, input sha256s, wall, peak RSS), output sha256s
      Scratch (stats, blocks, diags) stays in work/scratch/ and is not collected.
  uv run python experiments/803-de-engine/de803.py ka --name N (--new <stats.npz> | --arm surface/arm/seed) --ref <stats.npz>
      Known answer: arrays sums / jk / p / q / lfc bit-equal (NaN pattern included) -> PASS; otherwise max |Δ| per array.
      Writes work/out/ka/<N>.json.
  uv run python experiments/803-de-engine/de803.py members --surface S --out-root R --items arm=stats=diag[,...] --seed S
      (the members step alone; used by `run`.)

members.npz (T = the surface's targets, in the reference's order):
  targets; for every regime r in f1 / f0.61 / f0.5 / f0.35 / cmp: fid_<r>, jac_<r>, reach_<r>, nmae_<r> (528's rows_at on 521's
  truth pools, the realized order: exactly score619's functions) and nreal_<r> (pool size per target);
  ce_nmae, ce_fid, ce_reach, ce_jac, ce_pds: cell-eval2's own per-target values at f 1 (evaluate323.score: the scorer's
  dispatch_de_metrics and discrimination_score, panel genes out) and mse_num, mse_den: the scorer's per-target legs of
  expr_mse_unbiased_capped_norm (ratio of sums over targets = the member); ce_raw (json): the six raw members;
  K: ballast K per target; ncall: our calls per target (q < 0.05, own gene out); sums (T x G full axis, int64) and jk:
  the per-target pseudobulk sums and jackknife term the bulk members read.
Surfaces: h1test = 323's reference (h1-test, 85 targets, 435's harness); h1tv = 618's ref_tv (h1-trainval, 126 targets).
data/donors/h1_2025 never enters (the key is 323's / 618's cap15 on k562_gwps + orion_hct116 + orion_hek293t).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NS = os.environ.get("DE803_NS", "")          # optional namespace: parallel variants of the same arms on one instance
OUT = HERE / "work" / "out" / NS if NS else HERE / "work" / "out"
SCR = HERE / "work" / "scratch" / NS if NS else HERE / "work" / "scratch"
W618 = ROOT / "experiments/618-h1-trainval-surface/work"
SURF = {
    "h1test": dict(refdir=ROOT / "experiments/323-current-family-validation/work",
                   eff=ROOT / "experiments/323-current-family-validation/work/effects.npz",
                   st435=ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz",
                   n=85),
    "h1tv": dict(refdir=W618 / "ref_tv", eff=W618 / "tv/effects.npz",
                 st435=W618 / "g2mde/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz", n=126),
}
REGS = {"f1": "1.0", "f0.61": "0.61", "f0.5": "0.5", "f0.35": "0.35", "cmp": "cmp"}
MEM = ("fid", "jac", "reach", "nmae")
MDE = ["--mde-c", "0.1", "--mde-cap", "0.1", "--mde-min-count", "0.3", "--mde-coarse"]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def load_arm(name: str) -> dict:
    p = HERE / "arms" / f"{name}.json"
    a = json.loads(p.read_text())
    assert a["name"] == name, (a["name"], name)
    a["_path"], a["_sha256"] = str(p.relative_to(ROOT)), sha256(p)
    return a


def emit_args(emit: dict, seed: int, n_targets: int) -> list[str]:
    """run496's argv for one emitter configuration (its own tag logic names the outputs)."""
    a = ["base", str(seed), "--bulk", "ps"]
    if emit.get("mde", "coarse") == "coarse":
        a += MDE
    else:
        assert emit["mde"] == "none", emit
    if float(emit.get("zero_repair", -5)) < 0:
        a += ["--zero-repair", str(emit.get("zero_repair", -5))]
    if emit.get("zs", True):
        a += ["--zs", "dn_k250_zs"]
    if emit.get("dump", False):
        a += ["--dump", str(n_targets)]
    return a


def emit_key(emit: dict) -> str:
    return hashlib.sha256(json.dumps(emit, sort_keys=True).encode()).hexdigest()[:10]


def run_step(cmd: list[str], log: Path) -> dict:
    t0 = time.time()
    with open(log, "w") as f:
        rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT).returncode
    rec = dict(cmd=cmd, rc=rc, wall_s=round(time.time() - t0, 1),
               children_peak_rss_kb=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    if rc != 0:
        tail = log.read_text()[-3000:]
        raise SystemExit(f"step failed rc={rc}: {' '.join(cmd)}\n{tail}")
    return rec


def cmd_run(a):
    surf = SURF[a.surface]
    arms = [load_arm(n) for n in a.arms.split(",")]
    emit = arms[0]["emit"]
    for x in arms:
        assert x["emit"] == emit, f"arms do not share one emitter configuration: {x['name']}"
        assert a.surface in x.get("surfaces", list(SURF)), f"{x['name']} not defined on {a.surface}"
    variants = sorted({x["stage"].split(":", 1)[1] for x in arms if x["stage"].startswith("rnup:")})
    if variants:
        assert emit.get("dump") and emit.get("zs", True), "rnup stages need a zs run with --dump"
    key = emit_key(emit)
    rundir = SCR / a.surface / f"e{key}" / f"s{a.seed}"
    if rundir.exists():
        shutil.rmtree(rundir)
    rundir.mkdir(parents=True)
    (rundir / "emit.json").write_text(json.dumps(emit, indent=1, sort_keys=True))
    steps = []
    py = sys.executable
    ea = emit_args(emit, a.seed, surf["n"])
    if emit.get("sharpen"):          # follow-up 1: per-target bulk sharpening of the ballast render (emit803.install_sharpen)
        os.environ["DE803_SHARPEN"] = str(emit["sharpen"])
    if emit.get("bulk_transform"):   # §8: a module transforming the composition rows before key / sharpening / ballast
        os.environ["DE803_BULK_TRANSFORM"] = json.dumps(emit["bulk_transform"], sort_keys=True)
    try:
        steps.append(run_step([py, str(HERE / "emit803.py"), "emit", a.surface, emit.get("hook", "none"),
                               str(int(emit.get("callset_cpm", 2200))), str(rundir)] + ea, rundir / "emit.log"))
    finally:
        os.environ.pop("DE803_SHARPEN", None)
        os.environ.pop("DE803_BULK_TRANSFORM", None)
    mets = sorted(rundir.glob("metrics-*.json"))
    pre = [m for m in mets if "pre_tag" not in json.loads(m.read_text())]
    assert len(pre) == 1, [m.name for m in mets]
    tag = json.loads(pre[0].read_text())["tag"]
    stage_files = {"pre": rundir / f"stats-{tag}.npz", "zs": rundir / f"stats-{tag}_dnk250zs.npz"}
    diag = rundir / f"diag-{tag}.json"
    if variants:
        bl = rundir / "work" / "blocks"
        assert len(list(bl.glob(f"{tag}_*.npz"))) == surf["n"], "dumped blocks missing"
        bl.rename(rundir / "blocks")
        tmpl = tag.replace(f"base_s{a.seed}_", "base_s{seed}_", 1)
        assert tmpl != tag, tag
        steps.append(run_step([py, str(HERE / "emit803.py"), "rnup", a.surface, str(rundir), str(a.seed)] + variants
                              + ["--src", str(rundir), "--name", "e", "--tag", tmpl], rundir / "rnup.log"))
        for v in variants:
            stage_files[f"rnup:{v}"] = rundir / f"stats-rnup_e_s{a.seed}_{v}.npz"
        shutil.rmtree(rundir / "blocks")
    items = []
    for x in arms:
        sf = stage_files[x["stage"]]
        assert sf.exists(), sf
        items.append(f"{x['name']}={sf}={diag}")
    mlog = rundir / "members.log"
    steps.append(run_step([py, str(HERE / "de803.py"), "members", "--surface", a.surface, "--seed", str(a.seed),
                           "--items", ",".join(items)], mlog))
    for x in arms:
        od = OUT / a.surface / x["name"] / f"s{a.seed}"
        man = dict(argv=sys.argv, utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), surface=a.surface, seed=a.seed,
                   arm={k: v for k, v in x.items()}, emit_key=key, emit_args=ea, steps=steps,
                   step_manifests={p.stem: json.loads(p.read_text()) for p in sorted(rundir.glob("manifest-*.json"))},
                   stats=str(stage_files[x["stage"]].relative_to(ROOT)), stats_sha256=sha256(stage_files[x["stage"]]),
                   members_sha256=sha256(od / "members.npz"),
                   members_meta=json.loads((od / "members-meta.json").read_text()))
        (od / "manifest.json").write_text(json.dumps(man, indent=1, default=str))
    # the stats the KA step compares are kept; blocks were removed above
    print(f"[803] {a.surface} s{a.seed} {[x['name'] for x in arms]} done; scratch {rundir}")


def cmd_members(a):
    import numpy as np  # noqa: PLC0415
    surf = SURF[a.surface]
    sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
    sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
    sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))
    import lib521 as L  # noqa: PLC0415  (read-only; module constants re-pointed before the first load_state, as 618/619 do)
    L.REF = surf["refdir"] / "reference.npz"
    L.EFF = surf["eff"]
    L.ST435 = Path(str(surf["st435"]).replace("_s0_", "_s{seed}_"))
    L.load_state.cache_clear()
    import price528 as P5  # noqa: PLC0415
    from table521 import regime_pool  # noqa: PLC0415
    import run400 as H  # noqa: PLC0415  (R = run323, E = evaluate323)
    R, E = H.R, H.E
    R.WORK = surf["refdir"]
    S = L.load_state(0)
    pools = {k: regime_pool(S, r)[0] for k, r in REGS.items()}
    r, cc, _ = E.load_reference()
    targets = [str(t) for t in r["targets"]]
    assert targets == list(S["targets"]), "reference / lib521 target order differ"
    own = S["own"]
    for it in a.items.split(","):
        arm, stats, diag = it.split("=")
        z = np.load(stats, allow_pickle=False)
        assert [str(t) for t in z["targets"]] == targets
        A = dict(p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float))
        out = {"targets": np.asarray(targets)}
        for k, pool in pools.items():
            rows = P5.rows_at(S, A, pool)
            for j, m in enumerate(MEM):
                out[f"{m}_{k}"] = rows[j]
            out[f"nreal_{k}"] = pool.sum(1)
        sums = np.asarray(z["sums"])
        jk = np.asarray(z["jk"], float)
        raw, crow, legs = E.score(r, cc, list(sums), list(jk), list(A["p"]), list(A["q"]), list(A["lfc"]))
        full = {"nmae": "de_wilcoxon_lfc_nmae", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
                "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard", "pds": "pds_cosine"}
        for m, fn in full.items():
            d = {str(x["perturbation"]): x["value"] for x in crow if x["metric"] == fn}
            out[f"ce_{m}"] = np.array([d.get(t, np.nan) if d.get(t) is not None else np.nan for t in targets], float)
        out["mse_num"] = np.array([legs["numerator"].get(t, np.nan) for t in targets], float)
        out["mse_den"] = np.array([legs["denominator"].get(t, np.nan) for t in targets], float)
        out["ce_raw"] = np.asarray(json.dumps(raw))
        dg = json.loads(Path(diag).read_text())
        assert [d["target"] for d in dg] == targets
        out["K"] = np.array([int(d["ballast_K"]) for d in dg])
        out["ncall"] = ((A["q"] < 0.05) & ~own).sum(1)
        out["sums"] = sums.astype(np.int64)
        out["jk"] = jk
        # internal KA: 521's f 1 means == cell-eval2's own members (619's KA, here on every run)
        chk = {m: float(np.nanmean(out[f"{m}_f1"]) - raw[m]) for m in MEM}
        chk["mse_ratio_of_sums"] = float(np.nansum(out["mse_num"]) / np.nansum(out["mse_den"]) - raw["mse"])
        bad = {k: v for k, v in chk.items() if not abs(v) <= 1e-9}
        od = OUT / a.surface / arm / f"s{a.seed}"
        od.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(od / "members.npz", **out)
        meta = dict(arm=arm, stats=str(Path(stats).resolve().relative_to(ROOT)), raw=raw, ka_f1_vs_cell_eval2=chk,
                    ka_pass=not bad, mean_K=float(out["K"].mean()), mean_ncall=float(out["ncall"].mean()),
                    peak_rss_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (od / "members-meta.json").write_text(json.dumps(meta, indent=1))
        print(f"[members] {arm}: raw {json.dumps(raw)}; f1-vs-cell-eval2 {json.dumps(chk)}", flush=True)
        if bad:
            raise SystemExit(f"internal KA FAILED for {arm}: {bad}")


def cmd_ka(a):
    import numpy as np  # noqa: PLC0415
    if a.arm:   # surface/arm/seed -> the stats file the arm's manifest records
        sf, arm, sd = a.arm.split("/")
        a.new = str(ROOT / json.loads((OUT / sf / arm / f"s{sd}" / "manifest.json").read_text())["stats"])
    zn, zr = np.load(a.new, allow_pickle=False), np.load(a.ref, allow_pickle=False)
    rec = dict(name=a.name, new=a.new, ref=a.ref, new_sha256=sha256(Path(a.new)), ref_sha256=sha256(Path(a.ref)), arrays={})
    ok = True
    for k in ("targets", "sums", "jk", "p", "q", "lfc"):
        x, y = zn[k], zr[k]
        if k == "targets":
            eq = [str(t) for t in x] == [str(t) for t in y]
            rec["arrays"][k] = dict(equal=eq)
        else:
            xf, yf = np.asarray(x, float), np.asarray(y, float)
            eq = x.shape == y.shape and bool(np.array_equal(xf, yf, equal_nan=True))
            nanpat = x.shape == y.shape and bool(np.array_equal(np.isnan(xf), np.isnan(yf)))
            fin = np.isfinite(xf) & np.isfinite(yf) if x.shape == y.shape else None
            md = float(np.max(np.abs(xf[fin] - yf[fin]), initial=0.0)) if fin is not None else float("nan")
            rec["arrays"][k] = dict(equal=eq, nan_pattern_equal=nanpat, max_abs_diff=md,
                                    n_diff=int((xf[fin] != yf[fin]).sum()) if fin is not None else -1)
        ok &= eq
    rec["bit_equal"] = ok
    (OUT / "ka").mkdir(parents=True, exist_ok=True)
    (OUT / "ka" / f"{a.name}.json").write_text(json.dumps(rec, indent=1))
    print(f"[ka] {a.name}: {'PASS (bit-equal)' if ok else 'NOT bit-equal'} {json.dumps(rec['arrays'])}")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--surface", choices=list(SURF), required=True)
    r.add_argument("--arms", required=True)
    r.add_argument("--seed", type=int, required=True)
    m = sub.add_parser("members")
    m.add_argument("--surface", choices=list(SURF), required=True)
    m.add_argument("--seed", type=int, required=True)
    m.add_argument("--items", required=True)
    k = sub.add_parser("ka")
    k.add_argument("--name", required=True)
    k.add_argument("--new")
    k.add_argument("--arm", help="surface/arm/seed: compare the stats file that arm's manifest records")
    k.add_argument("--ref", required=True)
    a = ap.parse_args()
    {"run": cmd_run, "members": cmd_members, "ka": cmd_ka}[a.cmd](a)


if __name__ == "__main__":
    main()
