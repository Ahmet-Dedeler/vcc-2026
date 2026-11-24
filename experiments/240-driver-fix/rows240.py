"""240 — the cheap known-answer rows: env resolution before and after, and the guard firing.

Every row here is a MEASUREMENT of the real code path, not an argument about it: each case
builds a `Runner` from the driver under test (the pre-240 byte copy from git HEAD, or the
current `scripts/def_pipeline.py`), hands it the stage env exactly as the `build` stage does,
and runs `work/armprobe.py` as the child — which imports the real builder and prints the values
the builder itself parsed out of the environment it was handed.

    uv run python experiments/240-driver-fix/rows240.py

Writes `work/rows240.json` and one child log per case under `work/probe-<case>/logs/`.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
PROBE = WORK / "armprobe.py"
PY = ["uv", "run", "python"]
# The config the pre-240 driver produces on A/B/C: `rc_callset`'s, which is the incumbent it was
# pinned to. It is the artifact 239 §2.1 says the Oct-22 path would have built.
RC_CALLSET_META = ROOT / "submissions/202-rc_callset.h5ad.meta.json"
H1UNION_META = ROOT / "submissions/222-h1union-g028.h5ad.meta.json"


def load(path: Path, name: str, asif: Path | None = None):
    """Load a driver. `asif` runs the byte copy under the path it would occupy in `scripts/`,
    because both scripts derive `ROOT` from `__file__` -- 235's byte-copy trick, same reason."""
    spec = importlib.util.spec_from_file_location(name, asif or path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    if asif is None:
        spec.loader.exec_module(mod)
    else:
        exec(compile(path.read_text(), str(asif), "exec"), mod.__dict__)
    return mod


def probe(driver, case: str, operator_env: dict[str, str]) -> dict:
    """Run the driver's own Runner with the driver's own stage env, under `operator_env`."""
    work = WORK / f"probe-{case}"
    work.mkdir(parents=True, exist_ok=True)
    stage_env = getattr(driver, "SHIPPED_ENV", None) or driver.INCUMBENT_ENV
    saved = {k: os.environ.get(k) for k in ("VCC_ALPHA", "VCC_CD4_TARGET_F", "VCC_H1_UNION_G",
                                            "VCC_UNION_SOURCES")}
    for k in saved:
        os.environ.pop(k, None)
    os.environ.update(operator_env)
    try:
        runner = driver.Runner(work, dry=False)
        runner.run("probe", [*PY, str(PROBE)], stage_env)
    finally:
        for k, v in saved.items():
            os.environ.pop(k, None)
            if v is not None:
                os.environ[k] = v
    text = (work / "logs" / "probe.log").read_text()
    payload = json.loads(text[text.index("{"):text.rindex("}") + 1])
    return {"case": case, "driver": driver.__file__, "operator_env": operator_env,
            "driver_stage_env": stage_env, **payload}


def postbuild(script: Path, meta: Path, incumbent: Path, extra: list[str]) -> dict:
    args = ["--meta", str(meta), "--bundle", str(ROOT / "data/challenge"),
            "--incumbent", str(incumbent), "--cd4-root",
            str(ROOT / "data/external/cd4t/cloud-compact/cd4_compact"), *extra]
    if script.name.endswith("PRE240.py"):
        # The pre-240 bytes, run under the path they had in `scripts/` (see `load`).
        cmd = [*PY, str(WORK / "run_pre240.py"), str(script), *args]
    else:
        cmd = [*PY, str(script), *args]
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    return {"script": str(script), "meta": meta.name, "incumbent": incumbent.name,
            "extra": extra, "returncode": p.returncode, "stdout": p.stdout + p.stderr}


def main() -> None:
    pre = load(WORK / "def_pipeline.PRE240.py", "pre240",
               asif=ROOT / "scripts/def_pipeline.PRE240.py")
    post = load(ROOT / "scripts/def_pipeline.py", "post240")
    out: dict[str, object] = {"pre240_copies_match_git_head": {
        name: subprocess.run(["git", "show", f"HEAD:scripts/{name}.py"], cwd=ROOT,
                             capture_output=True, text=True).stdout
        == (WORK / f"{name}.PRE240.py").read_text()
        for name in ("def_pipeline", "def_postbuild")}}

    # --- ROW 2 / ROW 3: what the builder actually receives ---------------------------------
    cases = [
        ("pre240-asdocumented", pre, {}),
        ("pre240-operator-tries-both-arms", pre,
         {"VCC_CD4_TARGET_F": "0.34", "VCC_H1_UNION_G": "0.28"}),
        ("post240-defaults", post, {}),
        ("post240-operator-f010", post, {"VCC_CD4_TARGET_F": "0.10"}),
        ("post240-operator-arm-off", post, {"VCC_H1_UNION_G": ""}),
    ]
    out["env_resolution"] = [probe(d, name, env) for name, d, env in cases]

    # --- ROW 4: the guard fires ------------------------------------------------------------
    resolved, _ = post.resolve_env(post.SHIPPED_ENV)
    fired: dict[str, object] = {"resolved_intent": resolved}
    try:
        post.assert_build_carries_arms(RC_CALLSET_META, resolved)
        fired["raised"] = False
        fired["message"] = "DID NOT RAISE — the guard is not a guard"
    except SystemExit as exc:
        fired["raised"] = True
        fired["message"] = str(exc)
    # and the control: the same guard on the artifact that DOES carry the arms
    try:
        post.assert_build_carries_arms(H1UNION_META, resolved)
        fired["control_h1union_raised"] = False
    except SystemExit as exc:
        fired["control_h1union_raised"] = True
        fired["control_message"] = str(exc)
    out["guard"] = fired

    # --- ROW 3 second half + ROW 4 second half: the two postbuilds -------------------------
    out["postbuild"] = [
        # ROW 3: the pre-240 check, on the pre-240 driver's own artifact, with the pre-240
        # driver's own arguments. This is the run that PASSES while both arms are missing.
        {"label": "PRE-240 postbuild on the PRE-240 driver's config (rc_callset)",
         **postbuild(WORK / "def_postbuild.PRE240.py", RC_CALLSET_META, RC_CALLSET_META,
                     ["--allow-config-diff", "whitener_id"])},
        # ROW 4: the same artifact through the 240 check, with 240's default reference.
        {"label": "240 postbuild on the same config, default reference",
         **postbuild(ROOT / "scripts/def_postbuild.py", RC_CALLSET_META, H1UNION_META,
                     ["--allow-config-diff", "whitener_id"])},
        # and re-pointing the reference at the stale incumbent must NOT rescue it.
        {"label": "240 postbuild on the same config, reference re-pointed at rc_callset",
         **postbuild(ROOT / "scripts/def_postbuild.py", RC_CALLSET_META, RC_CALLSET_META,
                     ["--allow-config-diff", "whitener_id"])},
        # control: the shipped package through the 240 check must still pass.
        {"label": "240 postbuild on the SHIPPED h1union config (control)",
         **postbuild(ROOT / "scripts/def_postbuild.py", H1UNION_META, H1UNION_META, [])},
    ]

    (WORK / "rows240.json").write_text(json.dumps(out, indent=1) + "\n")
    print("\n\n================ SUMMARY ================")
    for row in out["env_resolution"]:
        b = row["builder_parsed"]
        print(f"{row['case']:<34} operator={row['operator_env']}  ->  builder read "
              f"alpha={b['alpha']} f={b['cd4_target_f']} g={b['h1_union_g']} "
              f"sources={b['h1_union_sources'] if b['h1_union_g'] is not None else '(inert)'}")
    print(f"\nguard on rc_callset's config raised: {out['guard']['raised']}; "
          f"on the shipped h1union: {out['guard']['control_h1union_raised']} (must be False)")
    for row in out["postbuild"]:
        tail = [ln for ln in row["stdout"].splitlines() if "checks pass" in ln]
        print(f"{row['label']:<62} rc={row['returncode']}  {tail[-1] if tail else ''}")
    print(f"\nwrote {WORK / 'rows240.json'}")


if __name__ == "__main__":
    main()
