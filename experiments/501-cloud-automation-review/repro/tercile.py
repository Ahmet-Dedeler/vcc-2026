#!/usr/bin/env python3
"""sharpen_q_by_tercile, copied from patch_builder479.py, on synthetic panels.

No stack dump: rows are synthetic key ratios / strengths. The A/B/C priced case
(no exact tie on a cut) is the control. An October-like panel has a large exact
tie at the cap and another at the floor.
"""
import numpy as np

def sharpen_q_by_tercile(bulk, carriers, genes, perts, pt):
    panel = set(map(str, perts))
    nonpanel = np.fromiter((str(g) not in panel for g in genes), dtype=bool, count=len(genes))
    cn = np.linalg.norm(np.asarray(carriers, dtype=np.float64)[:, nonpanel], axis=1)
    floor = np.zeros(len(cn), dtype=bool)
    if pt["mode"] == "stercile":
        v = cn
    else:
        bn = np.linalg.norm(np.asarray(bulk, dtype=np.float64)[:, nonpanel], axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.where(cn > 0, bn / np.where(cn > 0, cn, 1.0), np.nan)
        if not np.isfinite(v).any():
            raise SystemExit("no finite key ratio")
        vmin = float(np.nanmin(v))
        floor = np.isfinite(v) & (v <= vmin * (1.0 + 1e-9))
    fill = float(np.nanmedian(v))
    vv = np.where(np.isfinite(v), v, fill)
    edges = np.quantile(vv, [1.0 / 3.0, 2.0 / 3.0])
    t = np.digitize(vv, edges)
    tie_mid = np.zeros(len(vv), dtype=bool)
    if pt["mode"] == "tercile":
        for e in edges:
            blk = np.abs(vv - e) <= abs(e) * 1e-9
            if int(blk.sum()) >= 2:
                tie_mid |= blk
        t = np.where(tie_mid, 1, t)
        t = np.where(floor, 0, t)
    q = np.asarray(pt["q"], dtype=np.float64)[t]
    audit = {"mode": pt["mode"], "edges": [float(e) for e in edges],
             "counts": [int((t == i).sum()) for i in range(3)],
             "floor_to_low": int(floor.sum()), "tie_to_mid": int(tie_mid.sum()),
             "nonfinite": int((~np.isfinite(v)).sum())}
    return q, t, audit

def panel(ratios, zero_carrier=()):
    """One non-panel gene. bulk row norm == ratio, carrier norm == 1, except zero_carrier rows."""
    n = len(ratios)
    bulk = np.zeros((n, 1))
    carriers = np.ones((n, 1))
    for i, r in enumerate(ratios):
        bulk[i, 0] = r
    for i in zero_carrier:
        carriers[i, 0] = 0.0
        bulk[i, 0] = 0.0
    genes, perts = ["G0"], ["NOT_A_ROW"]
    return bulk, carriers, genes, perts

def show(title, ratios, mode, zero=()):
    b, c, g, p = panel(ratios, zero)
    q, t, a = sharpen_q_by_tercile(b, c, g, p, {"mode": mode, "q": [2.5, 2.0, 1.75]})
    # how many rows at the maximum landed in each bin
    r = np.asarray(ratios, float)
    rmax = np.nanmax(r)
    at_cap = np.where(np.isfinite(r) & (np.abs(r - rmax) <= 1e-12))[0]
    bins, counts = np.unique(t[at_cap], return_counts=True) if len(at_cap) else ([], [])
    cap_bins = dict(zip(map(int, bins), map(int, counts)))
    print(f"{title}")
    print(f"  mode={mode} edges={a['edges']} counts(L/M/H)={a['counts']} floor={a['floor_to_low']} "
          f"tie_to_mid={a['tie_to_mid']} nonfinite={a['nonfinite']}")
    print(f"  rows at the max value {rmax}: bin counts {cap_bins} (0=LOW,1=MID,2=HIGH)")

# control: 300 distinct ratios, no tie on a cut. Terciles should be 100/100/100, tie_to_mid 0.
show("CONTROL distinct ratios", np.linspace(0.5, 3.0, 300), "tercile")

# Cap of 49/300 (540: cap 3 binds 49 of 294 carried). The 2/3 cut falls below the cap.
ratios = np.concatenate([np.full(40, 1.0), np.linspace(1.05, 2.9, 211), np.full(49, 3.0)])
show("TERCILE cap block of 49 (does not contain the 2/3 cut)", ratios, "tercile")
# Cap block large enough to contain the 2/3 quantile.
ratios = np.concatenate([np.full(40, 1.0), np.linspace(1.2, 5.9, 110), np.full(150, 6.0)])
show("TERCILE cap block of 150 (contains the 2/3 cut)", ratios, "tercile")

# same panel, stercile, feeding the ratio in as the carrier norm directly
def show_strength(title, strength):
    n = len(strength)
    bulk = np.ones((n, 1))
    carriers = np.asarray(strength, float).reshape(n, 1)
    q, t, a = sharpen_q_by_tercile(bulk, carriers, ["G0"], ["X"], {"mode": "stercile", "q": [2.5, 2.0, 1.75]})
    print(f"{title}")
    print(f"  edges={a['edges']} counts(L/M/H)={a['counts']} tie_to_mid={a['tie_to_mid']} nonfinite={a['nonfinite']}")
    print(f"  digitize of a value exactly equal to the lower edge -> {int(np.digitize(np.array([a['edges'][0]]), a['edges'])[0])}"
          f"  (0 would be LOW; numpy default right=False returns the next bin)")
    print(f"  digitize of a value exactly equal to the upper edge -> {int(np.digitize(np.array([a['edges'][1]]), a['edges'])[0])}")

show_strength("STERCILE heavy tie at both cuts", np.concatenate([np.full(140, 1.0), np.linspace(1.2, 5.9, 80), np.full(80, 6.0)]))

# >2/3 of the rows identical: both quantiles collapse
show("TERCILE two-thirds tied at the floor", np.concatenate([np.full(220, 1.0), np.linspace(2.0, 6.0, 80)]), "tercile")
show_strength("STERCILE two-thirds tied", np.concatenate([np.full(220, 1.0), np.linspace(2.0, 6.0, 80)]))
show_strength("STERCILE 140 exact zeros (no carrier) then a spread", np.concatenate([np.zeros(140), np.linspace(0.2, 5.0, 160)]))

# non-finite / zero-carrier rows: print the bins those rows actually received
base = np.linspace(0.5, 3.0, 300)
b, c, g, p = panel(base, range(10))
q, t, a = sharpen_q_by_tercile(b, c, g, p, {"mode": "tercile", "q": [2.5, 2.0, 1.75]})
print("TERCILE 10 zero-carrier rows")
print(f"  counts={a['counts']} nonfinite={a['nonfinite']} floor={a['floor_to_low']} "
      f"bins of the 10 zero-carrier rows={t[:10].tolist()} (median fill, not LOW)")

# one huge tail
tail = np.concatenate([np.linspace(0.5, 2.0, 299), [1e6]])
show("TERCILE one 1e6 tail among 299 ordinary rows", tail, "tercile")

def sharpen_fixed(bulk, carriers, genes, perts, pt):
    """Same as the builder function, with the §6 assignment fixes."""
    panel = set(map(str, perts))
    nonpanel = np.fromiter((str(g) not in panel for g in genes), dtype=bool, count=len(genes))
    cn = np.linalg.norm(np.asarray(carriers, dtype=np.float64)[:, nonpanel], axis=1)
    floor = np.zeros(len(cn), dtype=bool)
    if pt["mode"] == "stercile":
        v = cn
    else:
        bn = np.linalg.norm(np.asarray(bulk, dtype=np.float64)[:, nonpanel], axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.where(cn > 0, bn / np.where(cn > 0, cn, 1.0), np.nan)
        if not np.isfinite(v).any():
            raise SystemExit("no finite key ratio")
        vmin = float(np.nanmin(v))
        floor = np.isfinite(v) & (v <= vmin * (1.0 + 1e-9))
    finite = np.isfinite(v)
    if not finite.any():
        raise SystemExit("no finite value")
    sample = v[finite]
    fill = float(np.median(sample))
    vv = np.where(finite, v, fill)
    edges = np.quantile(sample, [1.0 / 3.0, 2.0 / 3.0])
    t = np.digitize(vv, edges)
    tie_mid = np.zeros(len(vv), dtype=bool)
    if pt["mode"] == "tercile":
        for e in edges:
            if not (e <= vmin * (1.0 + 1e-9)):
                continue
            blk = np.abs(vv - e) <= abs(e) * 1e-9
            if int(blk.sum()) >= 2:
                tie_mid |= blk
        t = np.where(tie_mid, 1, t)
        t = np.where(floor, 0, t)
    elif float(edges[0]) == float(edges[1]):
        t = np.where(vv <= edges[0], 0, 2)
    else:
        e0 = float(edges[0])
        blk = np.abs(vv - e0) <= abs(e0) * 1e-9
        if int(blk.sum()) >= 2:
            t = np.where(blk, 0, t)
    t = np.where(finite, t, 1)
    q = np.asarray(pt["q"], dtype=np.float64)[t]
    audit = {"mode": pt["mode"], "edges": [float(e) for e in edges],
             "counts": [int((t == i).sum()) for i in range(3)],
             "floor_to_low": int(floor.sum()), "tie_to_mid": int(tie_mid.sum()),
             "nonfinite": int((~finite).sum())}
    return q, t, audit

def show2(title, ratios, mode, zero=()):
    b, c, g, p = panel(ratios, zero)
    q, t, a = sharpen_fixed(b, c, g, p, {"mode": mode, "q": [2.5, 2.0, 1.75]})
    print(f"FIXED {title}")
    print(f"  counts(L/M/H)={a['counts']} floor={a['floor_to_low']} tie_to_mid={a['tie_to_mid']} "
          f"nonfinite={a['nonfinite']} edges={a['edges']}")

def show2s(title, strength):
    n = len(strength)
    bulk = np.ones((n, 1))
    carriers = np.asarray(strength, float).reshape(n, 1)
    q, t, a = sharpen_fixed(bulk, carriers, ["G0"], ["X"], {"mode": "stercile", "q": [2.5, 2.0, 1.75]})
    print(f"FIXED {title}")
    print(f"  counts(L/M/H)={a['counts']} tie_to_mid={a['tie_to_mid']} edges={a['edges']}")

print("--- fixed assignment (control must stay 100/100/100, tie_to_mid 0) ---")
show2("CONTROL", np.linspace(0.5, 3.0, 300), "tercile")
show2s("CONTROL stercile", np.linspace(0.5, 3.0, 300))
show2("cap 49", np.concatenate([np.full(40, 1.0), np.linspace(1.05, 2.9, 211), np.full(49, 3.0)]), "tercile")
show2("cap 150", np.concatenate([np.full(40, 1.0), np.linspace(1.2, 5.9, 110), np.full(150, 6.0)]), "tercile")
show2s("140 at 1.0", np.concatenate([np.full(140, 1.0), np.linspace(1.2, 5.9, 80), np.full(80, 6.0)]))
show2("floor collapse", np.concatenate([np.full(220, 1.0), np.linspace(2.0, 6.0, 80)]), "tercile")
show2s("stercile collapse", np.concatenate([np.full(220, 1.0), np.linspace(2.0, 6.0, 80)]))
show2s("140 zeros", np.concatenate([np.zeros(140), np.linspace(0.2, 5.0, 160)]))
b, c, g, p = panel(np.linspace(0.5, 3.0, 300), range(10))
q, t, a = sharpen_fixed(b, c, g, p, {"mode": "tercile", "q": [2.5, 2.0, 1.75]})
print(f"FIXED 10 zero-carrier: counts={a['counts']} nonfinite={a['nonfinite']} bins={t[:10].tolist()}")
show2("one 1e6 tail", np.concatenate([np.linspace(0.5, 2.0, 299), [1e6]]), "tercile")
