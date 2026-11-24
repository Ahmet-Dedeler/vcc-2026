#!/usr/bin/env python3
"""614 KA for ref614.py: its H1 build must equal 323's reference.npz array for array (and the DE parquet by sha)."""
import json
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
a = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
b = np.load(Path(__file__).resolve().parent / "work/ref_h1ka/reference.npz", allow_pickle=False)
bad = [k for k in a.files if k not in b.files or not np.array_equal(a[k], b[k], equal_nan=a[k].dtype.kind == "f")]
ja = json.loads((ROOT / "experiments/323-current-family-validation/work/reference.json").read_text())
jb = json.loads((Path(__file__).resolve().parent / "work/ref_h1ka/reference.json").read_text())
print(f"KA ref614 on H1: arrays {'IDENTICAL' if not bad else 'DIFFER: ' + str(bad)}; DE parquet sha "
      f"{'equal' if ja['de_sha256'] == jb['de_sha256'] else 'DIFFERS'}; npz sha {'equal' if ja['statistics_sha256'] == jb['statistics_sha256'] else 'differs (key order / extra fields)'}")
sys.exit(0 if not bad and ja["de_sha256"] == jb["de_sha256"] else 1)
