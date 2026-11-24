#!/usr/bin/env python3
"""554 -- real-data identity: with every new knob unset, the patched def_preflight.py / def_postbuild.py print the same text
and return the same exit code as today's (538-patched) scripts, on every bundle and meta the coordinator listed. Plus the
guards' behaviour through the real CLIs: a wide-library bundle (R7b's controls as built, tripled) fails preflight and
passes with DEF_DAY_ALLOW_WIDE_LIBS=1; postbuild with --rc-audit on 433a (healthy) and 429a (the known failure).

Trees are built in a temp dir by test_guard554.make_tree (copies + 538 [+ 554]; data/ and submissions/ are read-only
symlinks). The only lines masked when comparing are the run-specific ones: preflight's elapsed-seconds line and its
"wrote <path>" line (the --out paths differ). Writes work/identity554.json.

    uv run python experiments/554-dose-guard/identity554.py
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import test_guard554 as G  # noqa: E402

PY = G.PY
BUNDLES = ["data/challenge", "experiments/180-def-readiness/work/def-bundle", "experiments/238-essential-free-panel/work/ef-bundle",
           "experiments/238-essential-free-panel/work/ef-bundle-r", "experiments/534-generalization-audit/work/partB/S1/bundle",
           "experiments/552-r7-reserve/work/T1/bundle"]
INC = "submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
DEFB = "experiments/180-def-readiness/work/def-bundle"
EFP, EFR = "experiments/238-essential-free-panel/work/ef-bundle", "experiments/238-essential-free-panel/work/ef-bundle-r"
METAS = [  # (meta, bundle, cd4 root or None for the default)
    ("submissions/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json", "data/challenge", None),
    ("submissions/446bj2zcp-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cp.h5ad.meta.json", "data/challenge", None),
    ("submissions/446bj2zcph-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cph.h5ad.meta.json", "data/challenge", None),
    ("submissions/453dr-october-path-abc-q1.75.h5ad.meta.json", DEFB, None),
    ("submissions/457dr-october-path-abc-cd4x15-q1.75.h5ad.meta.json", DEFB, None),
    ("submissions/460dr-october-path-lodo-abc-q1.75.h5ad.meta.json", DEFB, None),
    ("submissions/497dr-october-path-lodopct-spt-abc.h5ad.meta.json", DEFB, None),
    ("experiments/526-cloud-october-review2/work/efp_inert/474efp_final.h5ad.meta.json", EFP, "experiments/526-cloud-october-review2/work/efp_inert/cd4_compact"),
    ("experiments/497-october-switches/work/efr497/497efr_final.h5ad.meta.json", EFR, "experiments/497-october-switches/work/efr497/cd4_compact"),
    ("experiments/510-october-recipe-rehearsal/work/efr510/510efr_final.h5ad.meta.json", EFR, "experiments/510-october-recipe-rehearsal/work/efr510/cd4_compact"),
    ("experiments/452-def-day-rehearsal/work/efr/def452_final.h5ad.meta.json", EFR, "experiments/452-def-day-rehearsal/work/efr/cd4_compact"),
    ("experiments/474-october-hardening/work/efr_retrace/465efr_final.h5ad.meta.json", EFR, "experiments/474-october-hardening/work/efr_retrace/cd4_compact"),
]
MASK = [re.compile(r"^\d+(\.\d+)?s$"), re.compile(r"^wrote ")]


def clean_env(**extra) -> dict:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    env.update(extra)
    return env


def run(tree: Path, script: str, args: list[str], env: dict) -> tuple[int, str, str]:
    p = subprocess.run([str(PY), str(tree / "scripts" / script), *args], cwd=tree, env=env, capture_output=True, text=True)
    out = "\n".join(ln for ln in p.stdout.splitlines() if not any(m.search(ln.strip()) for m in MASK))
    return p.returncode, out, p.stderr


def main() -> None:
    t0 = time.time()
    tmp = Path(tempfile.mkdtemp(prefix="identity554-"))
    live = G.live_hashes()
    res = {"preflight": {}, "postbuild": {}, "behaviour": {}}
    try:
        base = G.make_tree(tmp / "base", with554=False)
        tree = G.make_tree(tmp / "patched", with554=True)
        print("== def_preflight.py, base vs patched (knobs unset) ==", flush=True)
        for b in BUNDLES:
            outs = {}
            for name, t in (("base", base), ("patched", tree)):
                o = tmp / name / "pre.json"
                outs[name] = run(t, "def_preflight.py", ["--bundle", str(ROOT / b), "--out", str(o)], clean_env())
                outs[name] = outs[name] + (json.loads(o.read_text()) if o.exists() else None,)
            same = outs["base"][:2] == outs["patched"][:2]
            spreads = {c: r.get("library_spread_q75_q25") for c, r in (outs["patched"][3] or {}).get("control_files", {}).items()}
            res["preflight"][b] = dict(rc=outs["patched"][0], identical_text_and_rc=same, spreads=spreads)
            print(f"  {b:58s} rc {outs['base'][0]}/{outs['patched'][0]}  identical {same}  q75/q25 {spreads}", flush=True)
        print("== def_postbuild.py, base vs patched (no --rc-audit; knobs unset) ==", flush=True)
        for meta, bundle, cd4 in METAS:
            args = ["--meta", str(ROOT / meta), "--bundle", str(ROOT / bundle), "--incumbent", str(ROOT / INC)]
            if cd4 and (ROOT / cd4).exists():
                args += ["--cd4-root", str(ROOT / cd4)]
            outs = {name: run(t, "def_postbuild.py", args + ["--out", str(tmp / name / "post.json")], clean_env())
                    for name, t in (("base", base), ("patched", tree))}
            same = outs["base"][:2] == outs["patched"][:2]
            dg = json.loads((tmp / "patched" / "post.json").read_text()).get("dose_guard") if (tmp / "patched" / "post.json").exists() else None
            res["postbuild"][meta] = dict(rc=outs["patched"][0], identical_text_and_rc=same,
                                          dose_guard_source=(dg or {}).get("source"),
                                          pooled=(dg or {}).get("per_context", {}).get("pooled"))
            pl = (dg or {}).get("per_context", {}).get("pooled") or {}
            print(f"  {Path(meta).name[:64]:64s} rc {outs['base'][0]}/{outs['patched'][0]}  identical {same}  "
                  f"pooled prior {pl.get('prior_lambda')} mde {pl.get('mde_lambda')} missed {pl.get('call_set_missed')}", flush=True)
            if not same:
                res["postbuild"][meta]["diff"] = {"base": outs["base"][1][-800:], "patched": outs["patched"][1][-800:]}
        print("== behaviour through the CLIs ==", flush=True)
        # (a) a wide-library bundle: data/challenge's manifest/panel/genes with R7b's as-built controls as every context
        wide = tmp / "wide"
        wide.mkdir()
        man = json.loads((ROOT / "data/challenge/manifest.json").read_text())
        src = ROOT / "experiments/552-r7-reserve/work/r7b/bundle/context_R.h5ad"
        import anndata as ad
        a = ad.read_h5ad(src)
        for c in man["contexts"]:
            b = a.copy()
            b.obs["context"] = c
            b.obs["target_gene"] = "non-targeting"
            b.write_h5ad(wide / f"context_{c}.h5ad")
            man["per_context"][c]["control_cells"] = int(a.n_obs)
        (wide / "manifest.json").write_text(json.dumps(man))
        for f in ("gene_names.csv", "pert_counts.csv"):
            shutil.copy(ROOT / "data/challenge" / f, wide / f)
        r0 = run(tree, "def_preflight.py", ["--bundle", str(wide), "--out", str(tmp / "w0.json")], clean_env())
        r1 = run(tree, "def_preflight.py", ["--bundle", str(wide), "--out", str(tmp / "w1.json")], clean_env(DEF_DAY_ALLOW_WIDE_LIBS="1"))
        rb = run(base, "def_preflight.py", ["--bundle", str(wide), "--out", str(tmp / "wb.json")], clean_env())
        w0, w1 = json.loads((tmp / "w0.json").read_text()), json.loads((tmp / "w1.json").read_text())
        row = lambda js: [c for c in js["checks"] if "q75/q25" in c["check"]]
        res["behaviour"]["wide_bundle"] = dict(
            base_spread_rows=len(row(json.loads((tmp / "wb.json").read_text()))),
            patched_rows=row(w0), allowed_rows=row(w1),
            allowed_flag={c: r.get("library_spread_allowed") for c, r in w1["control_files"].items()})
        print(f"  wide bundle (R7b controls as built, x3): base spread rows {res['behaviour']['wide_bundle']['base_spread_rows']}; "
              f"patched {[(r['pass'], r['value']) for r in row(w0)][:1]} x{len(row(w0))}; ALLOW {[(r['pass'], r['value'][-40:]) for r in row(w1)][:1]}", flush=True)
        # (b) postbuild --rc-audit, real dumps
        for tag, meta, dump in (
                ("433a", "submissions/433a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse.h5ad.meta.json",
                 "experiments/433-coarse-mde-arm/work/rc_audit_433a.jsonl"),
                ("429a", "submissions/429a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010-zr5-ctxlevel-zs250.h5ad.meta.json",
                 "experiments/429-zs-production/work/rc_audit_429a.jsonl")):
            for allow in ("", "1"):
                o = tmp / f"post_{tag}_{allow or 0}.json"
                rc, _, _ = run(tree, "def_postbuild.py", ["--meta", str(ROOT / meta), "--bundle", str(ROOT / "data/challenge"),
                                                          "--incumbent", str(ROOT / INC), "--rc-audit", str(ROOT / dump),
                                                          "--out", str(o)], clean_env(**({"DEF_DAY_ALLOW_LOW_DOSE": "1"} if allow else {})))
                js = json.loads(o.read_text())
                rows = [c for c in js["checks"] if "(554)" in c["check"]]
                res["behaviour"][f"postbuild_{tag}_allow{allow or 0}"] = dict(
                    rc=rc, rows=rows, per_context={c: {k: v.get(k) for k in ("prior_lambda", "mde_lambda", "call_set_missed")}
                                                   for c, v in js["dose_guard"]["per_context"].items()})
                print(f"  postbuild {tag} --rc-audit (ALLOW {allow or 0}): rc {rc}; 554 rows {[(r['pass'], r['value'][:110]) for r in rows]}", flush=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
        assert G.live_hashes() == live, "a live file changed"
    res["seconds"] = round(time.time() - t0, 1)
    res["all_identical"] = all(v["identical_text_and_rc"] for v in res["preflight"].values()) and \
        all(v["identical_text_and_rc"] for v in res["postbuild"].values())
    out = HERE / "work" / "identity554.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, default=str))
    print(f"\nall identical: {res['all_identical']}  ({res['seconds']} s) -> {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
