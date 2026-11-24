#!/usr/bin/env python3
"""Freeze Experiment 127 without loading any biological NumPy value array."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL = HERE / "preregistered_protocol.json"
FREEZE = HERE / "IMPLEMENTATION_FREEZE.json"
IMPLEMENTATION = ("run.py", "test_run.py", "freeze.py", "README.md")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def iter_bindings(protocol: dict) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for surface in protocol["surfaces"].values():
        for key in ("reference_de", "incumbent_prediction", "exact_prediction_de_table"):
            record = surface.get(key)
            if record:
                out.append((record["path"], record["sha256"]))
        calibration = surface.get("scaled_nmae_calibration")
        if calibration:
            out.extend((record["path"], record["sha256"]) for record in calibration["sources"])
    out.extend((record["path"], record["sha256"]) for record in protocol["shared_inputs"].values())
    return out


def main() -> None:
    if FREEZE.exists():
        raise FileExistsError(f"refusing to overwrite {FREEZE}")
    protocol = json.loads(PROTOCOL.read_text())
    verified = {}
    for relative, expected in iter_bindings(protocol):
        path = ROOT / relative
        actual = sha256_file(path)
        if actual != expected:
            raise RuntimeError(f"input hash mismatch: {relative}: {actual} != {expected}")
        verified[relative] = {"bytes": path.stat().st_size, "sha256": actual}
    implementation = {}
    for name in IMPLEMENTATION:
        path = HERE / name
        implementation[name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    payload = {
        "schema_version": 1,
        "experiment": "127-panel-metric-ceiling",
        "status": "frozen_before_biological_value_array_access",
        "protocol": {"bytes": PROTOCOL.stat().st_size, "sha256": sha256_file(PROTOCOL)},
        "inputs_verified_as_opaque_bytes": verified,
        "implementation": implementation,
        "biological_value_arrays_loaded": False,
    }
    FREEZE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

