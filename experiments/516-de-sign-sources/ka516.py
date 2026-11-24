"""516 known answer: the four H1 DE members of the production emitter, recomputed by common516 (walk513's loaders,
order and walk) from 435's coarse-MDE stats and 323's reference. Must equal walk513's / 435's numbers to 1e-9.
    PYTHONDONTWRITEBYTECODE=1 uv run python ka516.py
"""
import json

from common516 import WORK, known_answers, load_h1

d = load_h1()
got, ok, _ = known_answers(d)
line = (f"KNOWN ANSWER {'PASS' if ok else 'FAIL'}: fid {got['fid']!r} jac {got['jac']!r} nmae {got['nmae']!r} "
        f"reach {got['reach']!r} (want fid 0.5089291277653674 jac 0.0746985342 nmae 0.9427040518 reach 0.2692569776)")
print(line)
WORK.mkdir(exist_ok=True)
json.dump(dict(ok=ok, got=got), open(WORK / "ka516.json", "w"), indent=1)
if not ok:
    raise SystemExit(1)
