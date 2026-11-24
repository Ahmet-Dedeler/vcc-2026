#!/usr/bin/env python3
"""G: 549's EXPECT against the draft falsifiers, using the tracked 446bj2z meta and the cp snapshot.
No h5ad. Does not score pds or mse.
"""
from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
Z = ROOT / "submissions/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json"
CP = ROOT / "experiments/550-cloud-1005-review/snapshots/446bj2zcp.h5ad.meta.json"
KEYS = ROOT / "experiments/550-cloud-1005-review/repro/config_keys550.py"
EMITTER = ROOT / "src/vcc2026/rc_emit_k.py"
DRAFT = {"rc_callset_block", "de_mde_sign", "de_sign_emitter_sha256"}


def run(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(KEYS), *args], capture_output=True, text=True)


def main() -> int:
    z = json.loads(Z.read_text())["config"]
    cp = json.loads(CP.read_text())["config"]
    emitter = hashlib.sha256(EMITTER.read_bytes()).hexdigest()
    inc = z.get("de_sign_emitter_sha256", "")
    expect = {"rc_callset_block", "de_mde_sign"}
    if emitter != inc:
        expect.add("de_sign_emitter_sha256")
    print(f"emitter file {emitter[:16]}  446bj2z meta {inc[:16]}  EXPECT {sorted(expect)}")
    print(f"draft set {sorted(DRAFT)}  equal {expect == DRAFT}")
    if expect != DRAFT:
        return 1
    # the cp snapshot against 446bj2z: the block and the emitter pin, no hook
    a = run([str(CP), str(Z), "de_sign_emitter_sha256,rc_callset_block"])
    print("cp vs z", a.stdout.strip(), "rc", a.returncode)
    if a.returncode != 0:
        return 1
    # a built meta that is the cp config plus the hook, and nothing else
    built = {"config": copy.deepcopy(cp)}
    built["config"]["de_mde_sign"] = {"source": "clean4", "tau": 0.05}
    tmp = Path(tempfile.mkdtemp(prefix="550cfg-"))
    built_path = tmp / "built.meta.json"
    built_path.write_text(json.dumps(built))
    b = run([str(built_path), str(CP), "de_mde_sign"])
    print("built vs cp", b.stdout.strip(), "rc", b.returncode)
    if b.returncode != 0:
        return 1
    # a drifted swap is invisible against 446bj2z (the key is absent there either way)
    # and visible against the cp snapshot
    drifted = copy.deepcopy(built)
    drifted["config"]["rc_callset_block"]["swap"] = 0.5
    drift_path = tmp / "drift.meta.json"
    drift_path.write_text(json.dumps(drifted))
    c = run([str(drift_path), str(Z), "de_mde_sign,de_sign_emitter_sha256,rc_callset_block"])
    d = run([str(drift_path), str(CP), "de_mde_sign"])
    print("drifted swap vs z (key set only)", c.stdout.strip(), "rc", c.returncode)
    print("drifted swap vs cp", d.stdout.strip(), "rc", d.returncode)
    if c.returncode != 0 or d.returncode == 0:
        print("FAIL: the z key-set gate passes a drifted swap; the cp gate must not")
        return 1
    block = cp["rc_callset_block"]
    print("cp block", {k: block[k] for k in ("source", "swap", "steer", "sha256", "block_keep")})
    print("cp has de_mde_sign", "de_mde_sign" in cp)
    print("CONFIG KEYS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
