"""247 item A — `def_postbuild.py`'s operator block and row-count row, demonstrated FIRING.

Four metas, each derived from the shipped `submissions/222-h1union-g028.h5ad.meta.json`:

    KA-1  the shipped meta, untouched                       -> every row PASSES (the control)
    KA-2  operator block replaced by `whitener_union.npz`'s  -> the registry rows FAIL
          (sha 5608e807..., same n_cols, same axis by construction -- 242 showed 5/5 PASS)
    KA-3  operator.axis_sha256 set to 235's `def-bundle-badaxis` digest f2d6cc07...
                                                            -> the SYMBOL-axis row FAILS while
                                                               columns_sha256 still passes
    KA-4  output.n_rows 360000 -> 359600, context C 119600  -> the row-count rows FAIL
    KA-5  cd4 fully_supported_targets 276 -> 90 (+ pin)      -> the coverage floor FAILS

Each case is ALSO run through the byte copy of the pre-247 script (`work/pre247/def_postbuild.py`,
git HEAD before this experiment) -- the BLIND rows, which show the old guard passing KA-2 and
KA-4.  The pre-247 module is loaded with its `__file__` pinned to `scripts/` so its `ROOT` is
the repository root (stage222's trick).

    uv run python experiments/247-oct22-guards/postbuild247.py
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work" / "postbuild"
SHIPPED = ROOT / "submissions/222-h1union-g028.h5ad.meta.json"
UNION_JSON = ROOT / "data/donors/whitener_union.json"
BADAXIS = "f2d6cc07"   # 235 work/axis235.json: def-bundle-badaxis, prefix; full value read below


def load_script(path: Path, name: str, pinned_file: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(pinned_file)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def run(module, meta: Path) -> tuple[int, str]:
    argv = sys.argv
    sys.argv = ["def_postbuild.py", "--meta", str(meta), "--bundle", str(ROOT / "data/challenge")]
    buf = io.StringIO()
    rc = 0
    try:
        with contextlib.redirect_stdout(buf):
            module.main()
    except SystemExit as exc:
        rc = int(exc.code or 0)
    finally:
        sys.argv = argv
    return rc, buf.getvalue()


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    shipped = json.loads(SHIPPED.read_text())
    axis235 = json.loads((ROOT / "experiments/235-october-hardening/work/axis235.json").read_text())
    badaxis = axis235["def199_badaxis"]["axis_sha256"]
    assert badaxis.startswith(BADAXIS), badaxis
    union_meta = json.loads(UNION_JSON.read_text())
    import hashlib
    union_json_sha = hashlib.sha256(UNION_JSON.read_bytes()).hexdigest()

    cases: dict[str, dict] = {}
    cases["KA-1 shipped meta, untouched"] = json.loads(json.dumps(shipped))
    m = json.loads(json.dumps(shipped))
    m["operator"] = {
        **m["operator"],
        "path": str(ROOT / "data/donors/whitener_union.npz"),
        "metadata_path": str(UNION_JSON),
        "sha256": "5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3",
        "metadata_sha256": union_json_sha,
        "metadata": union_meta,
    }
    cases["KA-2 operator block = whitener_union.npz (5608e807..., same n_cols)"] = m
    m = json.loads(json.dumps(shipped))
    m["operator"]["axis_sha256"] = badaxis
    m["operator"]["axis_sha256_checked_against"] = "metadata:whitener_union_pruned.json"
    cases["KA-3 operator.axis_sha256 = def-bundle-badaxis f2d6cc07..."] = m
    m = json.loads(json.dumps(shipped))
    m["output"]["n_rows"] = 359600
    m["output"]["contexts"]["C"]["emitted_rows"] = 119600
    cases["KA-4 output.n_rows 359600, context C 119600 rows"] = m
    m = json.loads(json.dumps(shipped))
    m["effects"]["cd4"]["fully_supported_targets"] = 90
    m["frozen_bundle"]["panel_support"]["cd4_complete_targets"] = 90
    cases["KA-5 cd4 fully_supported_targets = 90 (30% coverage)"] = m

    new = load_script(ROOT / "scripts/def_postbuild.py", "def_postbuild_247", ROOT / "scripts/def_postbuild.py")
    old = load_script(HERE / "work/pre247/def_postbuild.py", "def_postbuild_pre247", ROOT / "scripts/def_postbuild.py")

    result = {}
    for name, meta in cases.items():
        path = WORK / (name.split()[0].replace("-", "").lower() + ".meta.json")
        path.write_text(json.dumps(meta, indent=1))
        rc_new, out_new = run(new, path)
        rc_old, out_old = run(old, path)
        fails_new = [ln.strip() for ln in out_new.splitlines() if ln.strip().startswith("FAIL")]
        fails_old = [ln.strip() for ln in out_old.splitlines() if ln.strip().startswith("FAIL")]
        tally_new = next((ln for ln in out_new.splitlines() if ln.endswith("checks pass")), "")
        tally_old = next((ln for ln in out_old.splitlines() if ln.endswith("checks pass")), "")
        result[name] = {"meta": str(path.relative_to(ROOT)),
                        "new": {"rc": rc_new, "tally": tally_new, "fail_lines": fails_new},
                        "pre247": {"rc": rc_old, "tally": tally_old, "fail_lines": fails_old}}
        (WORK / (path.stem + ".new.log")).write_text(out_new)
        (WORK / (path.stem + ".pre247.log")).write_text(out_old)
        print(f"\n=== {name}")
        print(f"  247:     rc={rc_new}  {tally_new}")
        for ln in fails_new:
            print(f"           {ln}")
        print(f"  pre-247: rc={rc_old}  {tally_old}")
        for ln in fails_old:
            print(f"           {ln}")

    # the verdict: the control passes both; KA-2 and KA-4 are BLIND spots of the old script
    verdict = {
        "KA-1 passes (247)": result["KA-1 shipped meta, untouched"]["new"]["rc"] == 0,
        "KA-2 FAILS (247)": result["KA-2 operator block = whitener_union.npz (5608e807..., same n_cols)"]["new"]["rc"] == 1,
        "KA-2 was BLIND (pre-247 passes it)": result["KA-2 operator block = whitener_union.npz (5608e807..., same n_cols)"]["pre247"]["rc"] == 0,
        "KA-3 FAILS on the symbol-axis row (247)": any("axis_sha256 ==" in ln for ln in result["KA-3 operator.axis_sha256 = def-bundle-badaxis f2d6cc07..."]["new"]["fail_lines"]),
        "KA-3 columns_sha256 row still passes (247)": not any("columns_sha256" in ln for ln in result["KA-3 operator.axis_sha256 = def-bundle-badaxis f2d6cc07..."]["new"]["fail_lines"]),
        "KA-3 was BLIND (pre-247 passes it)": result["KA-3 operator.axis_sha256 = def-bundle-badaxis f2d6cc07..."]["pre247"]["rc"] == 0,
        "KA-4 FAILS on the row-count rows (247)": any("n_rows" in ln for ln in result["KA-4 output.n_rows 359600, context C 119600 rows"]["new"]["fail_lines"]),
        "KA-4 was BLIND (pre-247 passes it)": result["KA-4 output.n_rows 359600, context C 119600 rows"]["pre247"]["rc"] == 0,
        "KA-5 FAILS on the coverage floor (247)": any(">= 100" in ln for ln in result["KA-5 cd4 fully_supported_targets = 90 (30% coverage)"]["new"]["fail_lines"]),
    }
    result["_verdict"] = verdict
    (WORK / "postbuild247.json").write_text(json.dumps(result, indent=1))
    print()
    for k, v in verdict.items():
        print(f"[{'PASS' if v else 'FAIL'}] {k}")
    ok = all(verdict.values())
    print(f"\nITEM A {'PASS' if ok else 'FAIL'}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
