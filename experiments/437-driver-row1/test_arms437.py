"""432: known-answer rows for the PROPOSED driver arms guard (def_pipeline.proposed.py), and the trap in the CURRENT one.

Rows (each prints PASS/FAIL; exit 1 if any FAIL):
  T1  current guard, A/B/C defaults, on cd4svd5's meta (no September knob at all)      -> ACCEPTS  (the trap, measured)
  T2  current guard, A/B/C defaults, on 433h's meta                                     -> accepts
  P1  proposed guard, A/B/C defaults + 433h's exported LOAD/CTX, on 433h's meta         -> accepts
  P2  proposed guard, same env, on 417a's meta (no MDE / zero repair / LOAD / CTX)       -> REFUSES, naming them
  P3  proposed guard, same env on cd4svd5's meta                                         -> REFUSES
  P4  proposed guard, + VCC_DE_MAG=dn_k250_zs exported, on 433h's meta (built without)   -> REFUSES on de_mag
  P5  proposed guard, D/E/F defaults on each 432 build meta found (b0: no LOAD/CTX; b1: + LOAD; b2: + LOAD + CTX + zs),
      with that build's own day-of exports                                                -> accepts
usage: uv run python proposals/test_arms432.py
"""
from __future__ import annotations

import glob
import importlib.util
import io
import json
import os
import sys
from contextlib import redirect_stdout
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W = ROOT / "experiments/432-def-rehearsal-new-knobs/work"
sys.path.insert(0, str(ROOT / "scripts"))


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CUR = load(ROOT / "scripts/def_pipeline.py", "def_pipeline_current")
NEW = load(ROOT / "scripts/def_pipeline.py", "def_pipeline_installed")
META = {k: glob.glob(str(ROOT / f"submissions/{k}-*.h5ad.meta.json"))[0] for k in ("249", "417a", "433h")}
LOAD_433H = "experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz:bb19f0e58bd00aa70e7aeb28b749fb3eee26d69dd62668e036da6370b8b792f2"
CTX_433H = "A=1.0088,B=1.05,C=1.0736"
DAYOF = ("VCC_BULK_LOAD", "VCC_BULK_CTX_SCALE", "VCC_DE_MAG")


def verdict(mod, meta: str, is_abc: bool, exports: dict[str, str]) -> tuple[bool, str]:
    saved = {k: os.environ.pop(k) for k in list(os.environ) if k.startswith("VCC_")}
    os.environ.update(exports)
    try:
        resolved, _ = mod.resolve_env(mod.day_env(is_abc))
        buf = io.StringIO()
        try:
            with redirect_stdout(buf):
                mod.assert_build_carries_arms(Path(meta), resolved)
            return True, buf.getvalue()
        except SystemExit as e:
            return False, buf.getvalue() + str(e)
    finally:
        for k in [k for k in os.environ if k.startswith("VCC_")]:
            del os.environ[k]
        os.environ.update(saved)


def main() -> None:
    rows = []

    def row(name, mod, meta, is_abc, exports, want_accept, must_name=()):
        ok, text = verdict(mod, meta, is_abc, exports)
        no_lines = [ln.strip() for ln in text.splitlines() if ln.strip().startswith("NO ")]
        named = all(any(f"config.{k}" in ln for ln in no_lines) for k in must_name)
        good = (ok == want_accept) and named
        rows.append(dict(row=name, meta=Path(meta).name[:40], accepted=ok, want_accept=want_accept, no_lines=no_lines,
                         result="PASS" if good else "FAIL"))
        print(f"{'PASS' if good else 'FAIL'}  {name}: accepted={ok} (want {want_accept}); NO rows: "
              f"{[ln.split(':')[0] for ln in no_lines]}", flush=True)

    row("T1 current guard accepts cd4svd5 (no September knob)", CUR, META["249"], True, {}, True)
    row("T2 current guard accepts 433h", CUR, META["433h"], True, {}, True)
    ex = {"VCC_BULK_LOAD": LOAD_433H, "VCC_BULK_CTX_SCALE": CTX_433H}
    row("P1 proposed guard accepts 433h", NEW, META["433h"], True, ex, True)
    row("P2 proposed guard refuses 417a", NEW, META["417a"], True, ex, False,
        ("de_mde", "de_zero_repair.floor_log2", "bulk_load_sha256", "bulk_ctx_scale"))
    row("P3 proposed guard refuses cd4svd5", NEW, META["249"], True, ex, False,
        ("normrel_exp", "cd4_partial_g", "de_sign_prior", "pert_scale", "pert_scale_ballast", "bulk_sharpen_q"))
    row("P4 proposed guard refuses 433h when zs was asked for", NEW, META["433h"], True, {**ex, "VCC_DE_MAG": "dn_k250_zs"},
        False, ("de_mag.variant",))
    for b in ("b0", "b1", "b2", "ef_b0", "ef_b1", "ef_b2"):
        meta = W / b / f"{b}.h5ad.meta.json"
        if not meta.exists():
            continue
        cfg = json.loads(meta.read_text())["config"]
        exb = {}
        if "bulk_load_sha256" in cfg:
            exb["VCC_BULK_LOAD"] = f"x:{cfg['bulk_load_sha256']}"
        if "bulk_ctx_scale" in cfg:
            exb["VCC_BULK_CTX_SCALE"] = ",".join(f"{k}={v}" for k, v in cfg["bulk_ctx_scale"].items())
        if "de_mag" in cfg:
            exb["VCC_DE_MAG"] = cfg["de_mag"]["variant"]
        row(f"P5 proposed guard accepts 432's {b} (D/E/F defaults + its day-of exports {sorted(exb)})",
            NEW, str(meta), False, exb, True)
    (ROOT / "experiments/437-driver-row1/work/test_arms437.json").write_text(json.dumps(rows, indent=1))
    bad = [r for r in rows if r["result"] != "PASS"]
    print(f"{len(rows) - len(bad)}/{len(rows)} rows PASS")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
