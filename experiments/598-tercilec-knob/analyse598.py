#!/usr/bin/env python3
"""598 -- tercilec on three real DAY-OF stack dumps (def_day.sh rehearsals), through the patched builder's own functions.

Each run's registry is rebuilt from its pin records (pin-artifacts.json, pin-lodo.json: id, path, sha256), which is what the
builder's PERT_SCALE_ARTIFACTS held during that run's level builds. dump_key_carried then resolves the dump env's
VCC_PERT_SCALE exactly as a level build would; tercile / stercile / tercilec are computed on the dump's own arrays.

  497dr  A/B/C relabelled as a new panel, LODO percentile key (lodo_lo0.5_pct_def)    experiments/497-october-switches/work/run
  efp    526's efp_inert: an October-like panel with many non-carried targets (410 key) experiments/526-cloud-october-review2/work/efp_inert
  efr    510's efr: realistic coverage, LODO percentile key                            experiments/510-october-recipe-rehearsal/work/efr510

Light: one 85 MB dump at a time. Writes work/analyse598.json.
  uv run python experiments/598-tercilec-knob/analyse598.py
"""
from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
X = ROOT / "experiments/598-tercilec-knob"
BUILDER = X / "work/builder.post598.py"     # test598's patched builder (the applied text)
RUNS = {"497dr": ROOT / "experiments/497-october-switches/work/run",
        "efp": ROOT / "experiments/526-cloud-october-review2/work/efp_inert",
        "efr": ROOT / "experiments/510-october-recipe-rehearsal/work/efr510"}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def ns_from_builder(registry: dict) -> dict:
    text = BUILDER.read_text()
    ns: dict = {"np": np, "Path": Path, "hashlib": hashlib, "ROOT": ROOT}
    for node in ast.parse(text).body:
        if isinstance(node, ast.FunctionDef) and node.name in ("sharpen_q_by_tercile", "dump_key_carried", "file_sha256"):
            exec(compile(ast.get_source_segment(text, node), "<builder>", "exec"), ns)
    ns["PERT_SCALE_ARTIFACTS"] = registry
    return ns


def main() -> int:
    out = {"builder_sha256": sha(BUILDER), "runs": {}}
    for name, d in RUNS.items():
        reg = {}
        for rec in sorted(d.glob("pin*.json")):
            ps = json.loads(rec.read_text()).get("pert_scale")
            if isinstance(ps, dict):
                reg[ps["id"]] = (ROOT / ps["path"], ps["sha256"])
        ns = ns_from_builder(reg)
        dump = d / "dump_h1x2_cd4x1.npz"
        with np.load(dump, allow_pickle=False) as z:
            bulk = np.asarray(z["bulk"], np.float64)
            car = np.asarray(z["carriers"], np.float64)
            perts = [str(x) for x in z["perts"]]
            genes = [str(x) for x in z["genes"]]
            env = json.loads(str(z["env"]))
        carried, key = ns["dump_key_carried"](env, perts, bulk, car, genes)
        res = {"dump": str(dump.relative_to(ROOT)), "dump_sha256": sha(dump), "dump_env_pert_scale": env.get("VCC_PERT_SCALE"),
               "registry": {k: [str(v[0].relative_to(ROOT)), v[1]] for k, v in reg.items()}, "key": key}
        for mode in ("tercile", "stercile", "tercilec"):
            kw = {"carried": carried} if mode == "tercilec" else {}
            q, a = ns["sharpen_q_by_tercile"](bulk, car, genes, perts, {"mode": mode, "q": [2.5, 2.0, 1.75]}, **kw)
            t = np.searchsorted(-np.asarray([2.5, 2.0, 1.75]), -q)   # q -> tier (q strictly decreasing in tier)
            res[mode] = {**a, "carried_counts": [int(((t == i) & carried).sum()) for i in range(3)],
                         "noncarried_counts": [int(((t == i) & ~carried).sum()) for i in range(3)]}
            if mode == "tercilec":
                res[mode]["n_moved"] = len(a["moved_vs_tercile"])
                res[mode].pop("moved_vs_tercile")
                res["moved_vs_tercile"] = a["moved_vs_tercile"]
        with np.load(reg[key["carried_key"]][0], allow_pickle=False) as z:
            w = np.asarray(z["w"], np.float64)
        res["w_noncarried_unique"] = sorted({round(float(x), 12) for x in w[~carried]})[:5]
        res["w_floor_carried"] = int((w[carried] <= w.min() + 1e-12).sum())
        out["runs"][name] = res
        tc, tp = res["tercilec"], res["tercile"]
        print(f"{name}: key {key['carried_key']} ({key['carried_key_sha256'][:12]}) maxdev {key['key_ratio_maxdev']:.1e}; "
              f"carried {int(carried.sum())}/300; tercile {tp['counts']} (non-carried {tp['noncarried_counts']}); "
              f"tercilec {tc['counts']} (carried {tc['carried_counts']}); moved {tc['n_moved']}; stercile {res['stercile']['counts']}",
              flush=True)
    (X / "work/analyse598.json").write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
