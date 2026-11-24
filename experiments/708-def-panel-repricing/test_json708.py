#!/usr/bin/env python3
"""Test for the 10-16 fix: dump708's json record must serialize with --fold-ops held (Path objects in spec["fold_ops"]),
numpy scalars and arrays, and score708 / day708 must be able to read it back. It builds a record shaped exactly like
dump708's rec_ (same keys and value types, from dump708's __main__ and KA block), writes it with dump708's encoder, then
reads it as score708 does (meta["spec"]["bundle"]). The old encoder (default=float) must FAIL on the same record.
usage: uv run python experiments/708-def-panel-repricing/test_json708.py
"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json708  # noqa: E402

FOLDS = ["k562_gwps", "orion_hct116", "orion_hek293t"]
spec = dict(bundle=str(Path("data/challenge-def").resolve()), rules="def", cd4_root=str(Path("x/cd4_compact").resolve()),
            full_op=str(Path("x/whitener_union_pruned_def.npz").resolve()), folds=FOLDS, synthetic=False, meta=None,
            fold_ops={h: Path(f"/tmp/heldops/held_{h}.npz").resolve() for h in FOLDS})   # Paths, as dump708's __main__
rec = dict(panel="def1022", fold="k562_gwps", n=300, sources=["orion_hct116", "orion_hek293t"],
           ka={"438 h2c1|V1_5_3+sh1.5": dict(realised=np.float64(0.0), design=np.float64(0.0), cos_agg=np.float64(0.07),
                                            pds=np.float64(0.89))},
           designs=["h2.5c1.5"], spec=spec, n_carried=np.int64(294), peak_rss_gib=np.float64(5.6), seconds=123.4,
           extra_array=np.arange(3), flag=np.bool_(True))
ok = True
try:
    json.dumps(rec, indent=1, default=float)
    print("FAIL: the old encoder (default=float) did not raise on this record; the test does not reproduce the bug")
    ok = False
except TypeError as e:
    print(f"PASS old encoder reproduces the bug: TypeError ({str(e)[:80]})")
txt = json.dumps(rec, indent=1, default=json708.default)
with tempfile.TemporaryDirectory() as d:
    f = Path(d) / "dump708-def1022-k562_gwps.json"
    f.write_text(txt)
    meta = json.loads(f.read_text())
    b = Path(meta["spec"]["bundle"])                     # score708's read
    fo = meta["spec"]["fold_ops"]["k562_gwps"]
    good = (b == Path(spec["bundle"]) and isinstance(fo, str) and fo.endswith("held_k562_gwps.npz")
            and meta["n_carried"] == 294 and meta["flag"] is True and meta["extra_array"] == [0, 1, 2])
    print(f"{'PASS' if good else 'FAIL'} new encoder: json written and read back (bundle {b.name}, fold_op {Path(fo).name}, "
          f"n_carried {meta['n_carried']})")
    ok &= good
sys.exit(0 if ok else 1)
