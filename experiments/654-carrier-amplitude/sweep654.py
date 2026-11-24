"""654 sweep (coordinator's request): zero-total gated cells (lfc < -30, own target excluded) in every H1-harness stats file
(shape 85 x 10779, lib521's gated axis) written since 09-29, inside lib521's f 1 / f 0.35 pools. Zero-repair runs (tag zr-5)
are listed separately. Non-H1 shapes are counted and skipped. Output: work/sweep654.{out,json}."""
import glob, json, os, sys, time
from collections import defaultdict
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
S = L.load_state(0)
OWN = S["own"]
P1, P35 = L.pools(S, 1.0) & ~OWN, L.pools(S, 0.35) & ~OWN
cut = 0.0 if "--all" in sys.argv else time.mktime(time.strptime("2026-09-29", "%Y-%m-%d"))
files = [f for f in glob.glob(str(ROOT / "experiments/**/stats-*.npz"), recursive=True) if os.path.getmtime(f) >= cut]
rows, other = [], defaultdict(int)
for f in sorted(files):
    try:
        z = np.load(f)
        if "lfc" not in z.files:
            other["no lfc"] += 1; continue
        l = z["lfc"]
    except Exception as e:  # noqa: BLE001
        other[f"unreadable {type(e).__name__}"] += 1; continue
    if l.shape != OWN.shape:
        other[f"shape {l.shape}"] += 1; continue
    zz = (np.asarray(l, float) < -30)
    nonown = zz & ~OWN
    rel = os.path.relpath(f, ROOT)
    rows.append(dict(file=rel, exp=rel.split("/")[1], zr="zr-5" in os.path.basename(f), own=int((zz & OWN).sum()),
                     nonown=int(nonown.sum()), in_f1=int((nonown & P1).sum()), in_f35=int((nonown & P35).sum())))
out = []
say = lambda s="": (print(s), out.append(s))  # noqa: E731
say(f"H1-shaped stats files ({'all dates' if '--all' in sys.argv else 'since 09-29'}): {len(rows)}; skipped: {dict(other)}")
flag = [r for r in rows if r["nonown"] > 0 and not r["zr"]]
say(f"without zero repair and with a non-own zero total: {len(flag)} files")
by = defaultdict(list)
for r in flag:
    by[r["exp"]].append(r)
for e, rs in sorted(by.items()):
    say(f"\n[{e}] {len(rs)} files; non-own zero-total cells per file: min {min(r['nonown'] for r in rs)} max {max(r['nonown'] for r in rs)}; "
        f"inside f 1 pools: max {max(r['in_f1'] for r in rs)}; inside f 0.35 pools: max {max(r['in_f35'] for r in rs)}")
    for r in sorted(rs, key=lambda r: -r["in_f1"])[:12]:
        say(f"    {r['file'].split('/work/')[1]}: non-own {r['nonown']}, in f1 {r['in_f1']}, in f0.35 {r['in_f35']}")
clean = defaultdict(int)
for r in rows:
    if not (r["nonown"] > 0 and not r["zr"]):
        clean[r["exp"]] += 1
say(f"\nclean files (no non-own zero total, or zero repair on) by experiment: {dict(sorted(clean.items()))}")
(Path(__file__).resolve().parent / ("work/sweep654all.json" if "--all" in sys.argv else "work/sweep654.json")).write_text(json.dumps(rows, indent=0))
(Path(__file__).resolve().parent / ("work/sweep654all.out" if "--all" in sys.argv else "work/sweep654.out")).write_text("\n".join(out) + "\n")
