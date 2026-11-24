"""511 (501 section 6): the per-target sharpening tie rules in `sharpen_q_by_tercile`.
- tercile: the tie-to-MID loop runs only for a cut at the floor (min r); a cap block containing the upper cut stays HIGH
  (493h: capped targets belong in HIGH).
- stercile: a lower-edge tie block of >= 2 goes LOW; equal cuts send the tied mass LOW and everything above it HIGH
  (before: LOW emptied, or every row got the mildest q).
- both: quantiles over finite values only; a non-finite key/strength lands in MID without moving the cuts.
On A/B/C (tie-free, all finite) the assignment is unchanged; `test_ties511.py` proves the q vectors equal the old ones on
446bj2's dump for tercile:2.5,2.0,1.75 and stercile:2.5,2.0,1.75 before any build.
usage: uv run python experiments/511-tercile-ties/patch_builder511.py <builder.py> [--check]
"""
import sys
from pathlib import Path

path = Path(sys.argv[1]); src = path.read_text()
old = '''    fill = float(np.nanmedian(v))
    vv = np.where(np.isfinite(v), v, fill)
    edges = np.quantile(vv, [1.0 / 3.0, 2.0 / 3.0])
    t = np.digitize(vv, edges)
    tie_mid = np.zeros(len(vv), dtype=bool)
    if pt["mode"] == "tercile":
        for e in edges:   # a block of identical ratios containing a cut goes to MID (non-carried targets sit at w = 1)
            blk = np.abs(vv - e) <= abs(e) * 1e-9
            if int(blk.sum()) >= 2:
                tie_mid |= blk
        t = np.where(tie_mid, 1, t)
        t = np.where(floor, 0, t)
    q = np.asarray(pt["q"], dtype=np.float64)[t]
'''
new = '''    # 511 (501 section 6): cuts over finite values only; tie blocks (within 1e-9 relative, which also absorbs the ~1e-14
    # float noise of r = c*w) are assigned explicitly instead of by digitize's edge semantics.
    finite = np.isfinite(v)
    if not finite.any():
        raise SystemExit("VCC_BULK_SHARPEN_Q: no finite value in the loaded dump")
    sample = v[finite]
    fill = float(np.median(sample))
    vv = np.where(finite, v, fill)
    edges = np.quantile(sample, [1.0 / 3.0, 2.0 / 3.0])
    t = np.digitize(vv, edges)   # right=False; the tie blocks below override it
    tie_mid = np.zeros(len(vv), dtype=bool)

    def _near(e):
        return np.abs(vv - e) <= abs(float(e)) * 1e-9
    if pt["mode"] == "tercile":
        vmin_ = float(np.min(sample))
        for e in edges:   # a cut AT THE FLOOR: its tie block -> MID, then the floor guard takes the floor rows LOW
            if e <= vmin_ * (1.0 + 1e-9):
                blk = _near(e)
                if int(blk.sum()) >= 2:
                    tie_mid |= blk
        t = np.where(tie_mid, 1, t)
        up = _near(edges[1])
        if int(up.sum()) >= 2 and not (edges[1] <= vmin_ * (1.0 + 1e-9)):
            t = np.where(up, 2, t)   # a cap block containing the upper cut stays HIGH (493h), float noise included
        t = np.where(floor, 0, t)
    else:
        lo = _near(edges[0])
        if abs(float(edges[1]) - float(edges[0])) <= abs(float(edges[0])) * 1e-9:
            t = np.where(vv <= edges[0] * (1.0 + 1e-9), 0, 2)   # collapsed cuts: the tied (weakest) mass LOW, the rest HIGH
        elif int(lo.sum()) >= 2:
            t = np.where(lo, 0, t)   # a tie block on the lower cut belongs to LOW, not MID
    t = np.where(finite, t, 1)       # a non-finite value does not enter the cuts; it sits in MID
    q = np.asarray(pt["q"], dtype=np.float64)[t]
'''
if src.count(old) != 1:
    raise SystemExit(f"anchor count {src.count(old)} != 1")
src = src.replace(old, new)
compile(src, str(path), "exec")
if "--check" in sys.argv:
    print("patch applies cleanly (not written)")
else:
    path.write_text(src)
    print(f"patched {path}")
