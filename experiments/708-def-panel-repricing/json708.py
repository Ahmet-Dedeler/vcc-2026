"""708: one JSON encoder for every record the 708 tools write (dump708's per-fold json, decisions).

dump708's record holds the day's spec, whose fold_ops values are pathlib.Path objects (--fold-ops held), plus numpy
scalars and arrays from the KA block. `json.dumps(..., default=float)` raised TypeError on the Paths AFTER the npz was
saved (rc 1, no json): the lab's EC2 run of its copy (dump801.py) found it; the coordinator verified it (10-16).
"""
from __future__ import annotations

from pathlib import PurePath


def default(o):
    """Paths -> str; numpy scalars -> Python scalars; numpy arrays -> lists; anything else -> str (never raises)."""
    if isinstance(o, PurePath):
        return str(o)
    if hasattr(o, "tolist") and hasattr(o, "shape"):
        return o.tolist()
    if hasattr(o, "item"):
        try:
            return o.item()
        except (TypeError, ValueError):
            pass
    return str(o)
