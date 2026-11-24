#!/usr/bin/env python3
"""gate524.py passes a wrong artifact. Exits 0 while that is still true.

Two finished-looking trees, both scored by the real gate524.py:

  dropped-sameline  — recipe asked for w 0.5 + noprog; the meta has no sameline record.
                      453dr's record differed, so the key is in the diff, and BASE_ALLOWED
                      contains it. The else branch logs the absence and does not fail.
  wrong-values      — sameline is present and satisfies the flag check, but de_mde_sign.tau
                      is 0.40 against recipe clean4:0.05 and de_prior_gate.le is 99 against
                      recipe le1. The gate only tests bool(meta key). The def_day log does
                      not name the programme bank; the gate never looks.

A third tree (right) is the same gate on a meta whose values match the recipe. It must PASS,
so the two false passes are not a broken harness. riders_match, which def_day runs, rejects
the wrong tau; the gate does not.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
GATE = ROOT / "experiments/524-october-recipe-dr/gate524.py"
REF = ROOT / "submissions/453dr-october-path-abc-q1.75.h5ad.meta.json"
HELPERS = ROOT / "scripts/def_day_helpers.py"
PIN = "1658ae2c400f75ee1659a8ffcaebca569bf3c8afd55ffd1c3006414344ef2791"

BASE_CFG = json.loads(REF.read_text())["config"]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


def tree(tmp: Path, name: str, cfg: dict, env: dict, flags: str, log_extra: str) -> tuple[Path, Path, Path]:
    w = tmp / name
    run = w / "run"
    run.mkdir(parents=True)
    (w / "recipe524.json").write_text(json.dumps({
        "env": env,
        "args": ["bundle", "work", "524dr", "1.75", "0.5", flags, "4"],
        "problems": [],
        "pending_treated_off": [],
    }))
    (run / "postbuild.json").write_text(json.dumps({"checks": [{"check": "synthetic", "pass": True}]}))
    (run / "def_day.log").write_text(
        "LODO key ON for the stack dump\n"
        "CONVERGED at D=1,E=1,F=1 after 1 levelled build(s)\n"
        "final: foo factors D=1,E=1,F=1 -> /tmp/foo.h5ad\n"
        + log_extra
    )
    (run / "foo.ntclevel.json").write_text(json.dumps(
        {"per_context": {"D": {"rho": 0.90}, "E": {"rho": 0.91}, "F": {"rho": 0.92}}}
    ))
    (run / "audit428.txt").write_text(
        "blocks 3; gated genes with ZERO block total (target excluded): 0 in 0 blocks; lfc < -5: 0\n"
    )
    stem = tmp / name / "pkg"
    (Path(str(stem) + ".h5ad.meta.json")).write_text(json.dumps({"config": cfg, "relabel453": {"categories_after": ["A", "B", "C"], "cells": 1}}))
    pkg = Path(str(stem) + ".vcc")
    pkg.write_bytes(b"synthetic-vcc-bytes")
    digest = sha256(pkg)
    chain = tmp / name / "chain.log"
    chain.write_text(f"novelty rc=0\npackage {pkg} sha256 {digest}\n")
    return w, stem, chain


def run_gate(w: Path, stem: Path, chain: Path) -> str:
    r = subprocess.run(
        [sys.executable, str(GATE), "--work", str(w), "--stem", str(stem), "--chainlog", str(chain)],
        cwd=ROOT, capture_output=True, text=True,
    )
    text = (r.stdout or "") + (r.stderr or "")
    if r.returncode not in (0, 1):
        raise SystemExit(f"gate524 rc={r.returncode}\n{text}")
    return text


def sl(prog: bool) -> dict:
    return {
        "path": "work/sameline.npz", "sha256": "ab" * 32, "w": 0.5, "pt": False, "prog": prog,
        "contexts": {"D": {"w": 0.5}, "E": {"w": 0.5}, "F": {"w": 0.5}},
    }


def base_moved() -> dict:
    cfg = json.loads(json.dumps(BASE_CFG))
    cfg["bulk_load"] = "dump:synthetic"
    cfg["bulk_load_sha256"] = "cd" * 32
    cfg["de_sign_emitter_sha256"] = PIN
    cfg["de_mag"] = {"variant": "dn_k250_zs"}
    return cfg


def main() -> int:
    ref = json.loads(REF.read_text())["config"]
    assert ref["de_sign_emitter_sha256"].startswith("edb5b8dc")
    env_plain = {"DEF_DAY_LODO": "1", "DEF_DAY_DE_MAG": "dn_k250_zs"}
    env_full = {
        **env_plain,
        "DEF_DAY_PROGRAMME_SOURCE": "h1_2025_full",
        "DEF_DAY_PRIOR_GATE": "le1",
        "DEF_DAY_MDE_SIGN": "clean4:0.05",
    }
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        # right: values match the plain recipe. Must PASS.
        cfg = base_moved()
        cfg["sameline"] = sl(False)
        text = run_gate(*tree(tmp, "right", cfg, env_plain, "noprog", ""))
        print("--- right artifact ---")
        print(text)
        assert "ALL GATES PASS" in text, text

        # wrong 1: same-line record dropped. Must PASS today (the defect).
        cfg = base_moved()
        cfg.pop("sameline", None)
        text = run_gate(*tree(tmp, "dropped", cfg, env_plain, "noprog", ""))
        print("--- dropped sameline ---")
        print(text)
        assert "ALL GATES PASS" in text, text
        assert "no record" in text

        # wrong 2: tau and le disagree with the recipe; programme bank is not in the log.
        cfg = base_moved()
        cfg["sameline"] = sl(True)
        cfg["de_mde_sign"] = {"source": "clean4", "tau": 0.40}
        cfg["de_prior_gate"] = {"le": 99, "tau": 0.2, "weak": "clean4"}
        w, stem, chain = tree(tmp, "values", cfg, env_full, "", "same-line mode ON\n")
        text = run_gate(w, stem, chain)
        print("--- wrong tau / le / unnamed programme ---")
        print(text)
        assert "ALL GATES PASS" in text, text
        assert "0.4" in text or "0.40" in text

        meta = Path(str(stem) + ".h5ad.meta.json")
        rm = subprocess.run(
            [sys.executable, str(HELPERS), "riders_match", str(meta), "dn_k250_zs", "clean4:0.05"],
            capture_output=True, text=True,
        )
        print("--- riders_match against the recipe tau (def_day's own check) ---")
        print(rm.stderr.strip() or rm.stdout.strip(), "rc", rm.returncode)
        assert rm.returncode == 1

    # The fixed predicates. Not applied; they are what patches/gate524_values.patch installs.
    ms = {"source": "clean4", "tau": 0.40}
    tau_ok = isinstance(ms, dict) and ms.get("source") == "clean4" and abs(float(ms["tau"]) - 0.05) < 1e-12
    pg = {"le": 99}
    le_ok = isinstance(pg, dict) and pg.get("le") == 1
    log = "same-line mode ON\n"
    bank_ok = "programme-source h1_2025_full" in log
    print(f"fixed predicates: tau_ok={tau_ok} le_ok={le_ok} bank_ok={bank_ok}")
    assert not tau_ok and not le_ok and not bank_ok
    print("DEFECT CONFIRMED: gate524 ALL GATES PASS on a dropped sameline record and on tau/le that disagree with the recipe")
    return 0


if __name__ == "__main__":
    sys.exit(main())
