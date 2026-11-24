"""455: the October driver's SHIPPED q follows the board (446b, 2026-09-26: q 1.75 levelled, +0.0022834 vs q 1.5).
Two literal edits in scripts/def_pipeline.py; anchors must each match once.  usage: patch455.py [--check]"""
import sys
from pathlib import Path

p = Path(__file__).resolve().parents[2] / "scripts/def_pipeline.py"
s = p.read_text()
edits = [
    ('"VCC_BULK_SHARPEN_Q": "1.25",', '"VCC_BULK_SHARPEN_Q": "1.75",'),
    ('"VCC_BULK_SHARPEN_Q": "417a bulk sharpening q 1.25: +0.0069612",',
     '"VCC_BULK_SHARPEN_Q": "bulk sharpening, levelled: q 1.25 +0.0069612 (417a), q 1.5 +0.0046475 (433k), '
     'q 1.75 +0.0022834 (446b, UDgjOKL9gI1udAexq6ue)",'),
]
for a, _ in edits:
    if s.count(a) != 1:
        raise SystemExit(f"anchor count {s.count(a)} != 1: {a[:70]!r}")
if "--check" in sys.argv:
    print("anchors OK; not writing")
else:
    for a, b in edits:
        s = s.replace(a, b)
    p.write_text(s)
    print("patched", p)
