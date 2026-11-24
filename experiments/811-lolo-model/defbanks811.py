"""811 -- turn kit811ext.py's extended banks (work/banks_ext, its own manifest schema) into hook811's bank format (work/banks_def,
ext811.py's manifest schema), so the C0 dump hook can blend D/E/F's real panel with EXACTLY the rows 811 priced on DEFo / DEF.
Coordinator's X2 ruling (SIXA 344c4c31), item 1. Written 2026-10-22 by the lab lead.

kit811ext.py (801's copy of ext811a.py, unchanged; the code that built the priced DEFo / DEF rows on EC2) writes
banks_ext/<bank>.npz = the pinned banks' rows byte-for-byte + the new panels' rows appended, and a manifest WITHOUT the "out" key that
hook811's load811() checks. This adapter changes no row and no array:
  - checks the old manifest is the pinned one (37e34335...), every banks_ext file against banks_ext's manifest, and that each
    banks_ext bank was merged from the pinned bank file (old_sha256 == the pinned out_sha256);
  - copies every banks_ext/<bank>.npz byte-for-byte to banks_def/<bank>.npz (sha256 asserted equal);
  - writes banks_def/manifest.json in ext811's schema: panels = the pinned panels + the kit panels (KIT811_PANELS), per bank
    sources / n_targets / n_zero_rows / n_measured / coverage / out (= work/banks/<bank>.npz, where it will live) / out_sha256,
    plus provenance (extended_from, kit811ext's manifest sha256, its L1 CD4 KA line, the CD4 compacts' sha256).
The caller then moves work/banks -> work/banks.pre-def and work/banks_def -> work/banks, and re-pins hook811.BANKS_MANIFEST_SHA.
    KIT811_PANELS=DEF uv run python experiments/811-lolo-model/defbanks811.py
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD, EXT, NEW = HERE / "work" / "banks", HERE / "work" / "banks_ext", HERE / "work" / "banks_def"
PINNED = "37e343351000b9b65ceb4adc6a18004563219908cf03b805ed8b27215e57d83b"
KIT811EXT = ROOT / "experiments/801-lab-folds/kit/kit811ext.py"
PANELS = os.environ["KIT811_PANELS"].split()


def sha(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def targets(p: Path) -> list[str]:
    with open(p) as f:
        return [r["target_gene"] for r in csv.DictReader(f)]


def main():
    assert sha(OLD / "manifest.json") == PINNED, "work/banks is not the pinned bank set"
    mo, me = json.loads((OLD / "manifest.json").read_text()), json.loads((EXT / "manifest.json").read_text())
    assert set(mo["banks"]) == set(me["banks"]), "bank lists differ"
    panels = dict(mo["panels"])
    for p in PANELS:
        assert p not in panels, p
        panels[p] = f"experiments/801-lab-folds/work/{p}/bundle/pert_counts.csv"
    pans = {k: targets(ROOT / v) for k, v in panels.items()}
    NEW.mkdir(parents=True, exist_ok=False)
    banks = {}
    for b in sorted(mo["banks"]):
        o, e = mo["banks"][b], me["banks"][b]
        src = EXT / f"{b}.npz"
        assert sha(src) == e["out_sha256"], f"{b}: banks_ext file differs from its manifest"
        assert e["old_sha256"] == o["out_sha256"], f"{b}: banks_ext was not merged from the pinned bank"
        dst = NEW / f"{b}.npz"
        shutil.copyfile(src, dst)
        s = sha(dst)
        assert s == e["out_sha256"], f"{b}: copy differs"
        z = np.load(dst)
        T, R, M = [str(t) for t in z["targets"]], np.asarray(z["rows"]), np.asarray(z["measured"], bool)
        assert T[:o["n_targets"]] == [str(t) for t in np.load(OLD / f"{b}.npz")["targets"]], f"{b}: old target order changed"
        banks[b] = dict(sources={**o["sources"], **e["sources"]}, n_targets=len(T), n_zero_rows=int((np.abs(R).sum(1) == 0).sum()),
                        n_measured=int(M.sum()), coverage={k: len(set(v) & set(T)) for k, v in pans.items()},
                        out=f"experiments/811-lolo-model/work/banks/{b}.npz", out_sha256=s,
                        extended_from_sha256=o["out_sha256"], n_old=e["n_old"], n_added=e["n_added"],
                        measured_changed=bool(o["n_measured"] != int(M.sum())))
        print(f"{b:16s} {o['n_targets']:5d} + {e['n_added']:4d} targets, measured {o['n_measured']} -> {int(M.sum())}, "
              + " ".join(f"{k}:{v}" for k, v in banks[b]["coverage"].items() if k in PANELS + ["ABC"]), flush=True)
    man = dict(panels=panels, n_wanted=len(set().union(*map(set, pans.values()))), banks=banks,
               extended_from_manifest_sha256=PINNED, kit811ext_sha256=sha(KIT811EXT),
               ext_manifest_sha256=sha(EXT / "manifest.json"), ka_L1_cd4=me.get("ka_L1_cd4"), cd4=me.get("cd4"),
               adapter_sha256=sha(Path(__file__)))
    (NEW / "manifest.json").write_text(json.dumps(man, indent=1))
    print(f"banks_def/manifest.json sha256 {sha(NEW / 'manifest.json')}", flush=True)


if __name__ == "__main__":
    main()
