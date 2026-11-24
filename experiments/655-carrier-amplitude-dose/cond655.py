"""655's run condition (PREREGISTRATION): 654's production-faithful x1.25 - x1 overall-eq > 0 at f 0.35 on >= 2 of 3 seeds.
Reads 654's work/analyze654zr.json (re-run first by seq655.sh). Prints GO / WITHDRAW and exits 0 / 3."""
import json, sys
d = json.load(open("experiments/654-carrier-amplitude/work/analyze654zr.json"))["deltas"]
vals = [d[f"x125|0.35|{s}"]["overall"] for s in (0, 1, 2) if f"x125|0.35|{s}" in d]
pos = sum(v > 0 for v in vals)
print(f"655 condition: x125 zr - x1 at f 0.35 per seed {[round(v, 5) for v in vals]} -> {pos}/{len(vals)} positive")
ok = len(vals) == 3 and pos >= 2
print("GO" if ok else "WITHDRAW")
sys.exit(0 if ok else 3)
