"""801 -- write work/registry.json: every fold id with its dump path + sha256, its operator sha256, target count, coverage,
and status (working / abc / ka / sealed-drawn). Reads only small files and hashes the dumps (light, seconds).

  L<k>/<held>   WORKING (L1-L6): work/dumps/dump708-L<k>-<held>.npz (EC2, run G), operator s3://<bucket>/output/lab/ops801/L<k>/
  ABC/<held>    A/B/C's own panel through 518's dumps (k562_gwps / orion_hct116 / orion_hek293t / kolf2); operators: 062's
                held-out fits (plolo.OPERATOR_REGISTRY) and, for kolf2, the production operator (518's path)
  efp/<held>    the known-answer folds (238's efp inputs), when their EC2 dump exists
  L7, L8        sealed-drawn: target-list sha256 only; no fold exists

    uv run python experiments/801-lab-folds/registry801.py <run-A dir with out/repo/...>
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BUCKET = "<your-s3-bucket>"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
ABC_OPS = {  # 518's dumps were built on these operators (plolo.OPERATOR_REGISTRY; kolf2 = lib438.PROD_OP)
    "k562_gwps": ("experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz",
                  "13d988fcc3b129d25d254524db73e8f8253143fa493ee64d9654f0f61bc3ac55"),
    "orion_hct116": ("experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hct116.npz",
                     "7cecbf7f36e40d22a6bea1ff3adc86a8788b329def397044fb5b9ba675247cbc"),
    "orion_hek293t": ("experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hek293t.npz",
                      "97e486e0d1a805e1f24a622f1cefbfb0917e9282c5d0491c4ddb41c8475f2af1"),
    "kolf2": ("data/donors/whitener_union_pruned.npz", "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"),
}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def dump_entry(path: Path) -> dict:
    z = np.load(path)
    t = [str(x) for x in z["targets"]]
    return dict(dump=str(path.relative_to(ROOT)), dump_sha256=sha256(path), n_targets=len(t),
                n_carried=int(z["carried"].sum()), sources=[str(s) for s in z["sources"]], coverage=f"{len(t)}/300")


def main(run_dir: Path):
    draws = json.loads((WORK / "draws.json").read_text())["draws"]
    out = dict(written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), draws_sha256=sha256(WORK / "draws.json"),
               engine_sha256=sha256(HERE / "engine801.py"), dump_script_sha256=sha256(HERE / "dump801.py"), folds={}, sealed={})
    cloud = run_dir / "out/repo/experiments/801-lab-folds/work"
    for k in ("L1", "L2", "L3", "L4", "L5", "L6"):
        shas = {}
        for h in FOLDS:
            for ln in (cloud / "cloud" / f"{k}-{h}-outputs.sha256").read_text().splitlines():
                s, p = ln.split(None, 1)
                shas[p.strip()] = s
        for h in FOLDS:
            p = WORK / "dumps" / f"dump708-{k}-{h}.npz"
            e = dump_entry(p)
            assert e["dump_sha256"] == shas[f"dumps/dump708-{k}-{h}.npz"], f"{p}: sha differs from the EC2 record"
            opj = json.loads((cloud / k / "ops" / f"held_{h}.json").read_text())
            assert opj["sha256"] == shas[f"{k}/ops/held_{h}.npz"]
            e.update(status="working", panel=k, held=h, panel_sha256=draws[k]["sha256"],
                     operator=f"s3://{BUCKET}/output/lab/ops801/{k}/held_{h}.npz", operator_sha256=opj["sha256"],
                     operator_axis_sha256=opj["axis_sha256"], operator_donors=opj["donors"],
                     dump_record=json.loads((WORK / "dumps" / f"dump708-{k}-{h}.json").read_text()).get("ka"))
            out["folds"][f"{k}/{h}"] = e
    for h in FOLDS + ("kolf2",):
        p = ROOT / f"experiments/518-bulk-mse-gap/work/dump518-{h}.npz"
        op, ops = ABC_OPS[h]
        out["folds"][f"ABC/{h}"] = dict(status="abc", panel="ABC", held=h, dump=str(p.relative_to(ROOT)), dump_sha256=sha256(p),
                                        operator=op, operator_sha256=ops or "see plolo.OPERATOR_REGISTRY / lib438.PROD_OP",
                                        loader="fold522.Abc (score518.Fold)")
    for h in FOLDS:
        p = WORK / "dumps" / f"dump708-efp-{h}.npz"
        if p.exists():
            e = dump_entry(p)
            e.update(status="ka", panel="efp", held=h, operator=f"experiments/238-essential-free-panel/work/ops/efp__held_{h}.npz")
            out["folds"][f"efp/{h}"] = e
    for k in ("L7", "L8"):
        out["sealed"][k] = dict(status="sealed-drawn", seed=draws[k]["seed"], sha256=draws[k]["sha256"],
                                target_list=draws[k]["target_list"])
    (WORK / "registry.json").write_text(json.dumps(out, indent=1))
    print(f"registry: {len(out['folds'])} folds, {len(out['sealed'])} sealed -> work/registry.json sha256 {sha256(WORK / 'registry.json')}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
