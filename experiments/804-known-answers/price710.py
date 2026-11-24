"""804 — the prospective row 710 (H1 dose: VCC_H1_UNION_G 0.70 = H1 x2.5 vs ALLSW's x2) on every bulk fold, through 801's
engine, in 800's STANDARD FORMAT (PROTOCOL.md amendment 2).

Arms per fold: `inc` (the fold's incumbent: L folds pct + P2g on fin = h2c1; A/B/C folds cap6 + P2g(abc) on fin = 518's fST =
h2c1.5, ALLSW's stack) and `710` (the same key weights and sharpening, composition rows swapped for 801's H1 x2.5 design:
fin_h2.5c1 on L folds, fin_h2.5c1.5 on A/B/C folds). Key weights stay: the LODO key is read from the h1c1 carriers, not from
the H1 dose; P2g's cut reads the key ratio only. KA: the design file's base composition equals the fold's fin (|d| < 1e-12).
    KA804_OUT=/work/ka804 uv run python experiments/804-known-answers/price710.py <fold-id>
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/801-lab-folds"))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import engine801 as E  # noqa: E402

OUT = Path(os.environ.get("KA804_OUT", str(HERE / "work" / "run")))
DES = ROOT / "experiments/801-lab-folds/work/designs"


def main(fold_id: str):
    t0 = time.time()
    F = E.load(fold_id, check=True)
    abc = F.panel == "ABC"
    z = np.load(DES / f"{fold_id.replace('/', '__')}.npz")
    base_key, cand_key = ("fin_h2c1.5", "fin_h2.5c1.5") if abc else ("fin_h2c1", "fin_h2.5c1")
    assert [str(t) for t in z["targets"]] == list(F.targets), "design targets != fold targets"
    ka = float(np.max(np.abs(np.asarray(z[base_key], np.float64) - np.asarray(F.fin, np.float64))))
    assert ka < 1e-12, f"KA: {base_key} != fold fin ({ka})"
    inc = E.named(F, "cap6+P2g(abc)" if abc else "pct+P2g")
    arms = {"inc": E.Arm(w=inc.w, q=inc.q, name="inc"),
            "710": E.Arm(w=inc.w, q=inc.q, fin=np.asarray(z[cand_key], np.float64), name="710")}
    rec = dict(fold=fold_id, ka_base_composition=ka, base=base_key, cand=cand_key, arms={})
    d = OUT / "bulk" / fold_id.replace("/", "__")
    d.mkdir(parents=True, exist_ok=True)
    for nm, arm in arms.items():
        s, arr = E.price(F, arm, base=arms["inc"] if nm == "710" else None)
        keep = {k: v for k, v in arr.items() if not k.startswith("base_")}
        keep["targets"] = np.asarray(F.targets)
        np.savez_compressed(d / f"{'inc710' if nm == 'inc' else nm}.npz", **keep)
        rec["arms"][nm] = dict(cos_agg=s["cos_agg"], pds=s["pds"], K_mean=s["K_mean"], vs_base=s.get("vs_base"))
        print(f"[{fold_id}] {nm:4s} cos_agg {s['cos_agg']:.6f} pds {s['pds']:.5f}" +
              (f" | r {100 * s['vs_base']['r']:+.3f}% dpds {s['vs_base']['dpds']:+.5f}" if nm == "710" else ""), flush=True)
    (OUT / "meta").mkdir(parents=True, exist_ok=True)
    (OUT / "meta" / f"710-{fold_id.replace('/', '__')}.json").write_text(json.dumps(rec, indent=1, default=float))
    print(f"[{fold_id}] done ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
