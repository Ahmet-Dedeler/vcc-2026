"""317 — the six-member sweep table: raw members per operating point, and their board price."""
import json
from pathlib import Path
import numpy as np
from read_agg import read, MEMBERS
from price317 import board_slopes

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CO = HERE / "cloud-output"
L221 = REPO / "experiments/221-cloud-scorer/cloud-output"

# `base` = 221's `rc_callset`, rebuilt byte-identically on every worker (gate 5). Its recorded
# row is the anchor; KA-1 re-scores the same file to prove the instrument reproduces it.
BASE_221 = L221 / "ladder8-20260904/rc_callset/run/agg_results.csv"
TWIN = (L221 / "ladder6-20260904/a15/run/agg_results.csv",
        L221 / "ladder7-20260904/a15_seed1/run/agg_results.csv")

ARMS = [("od1", "--overdispersion 4.0 -> 1.0", "emission"),
        ("od8", "--overdispersion 4.0 -> 8.0", "emission"),
        ("ks010", "--keep-spread 0.25 -> 0.10", "emission"),
        ("ts070", "--tail-spread 1.0 -> 0.70", "emission"),
        ("kt16000", "--keep-top 8000 -> 16000", "emission"),
        ("amp070", "amplitude x0.70 (alpha 0.982 -> 0.687)", "amplitude"),
        ("amp085", "amplitude x0.85 (alpha 0.982 -> 0.835)", "amplitude"),
        ("amp120", "amplitude x1.20 (alpha 0.982 -> 1.178)", "amplitude"),
        ("amp145", "amplitude x1.45 (alpha 0.982 -> 1.424)", "amplitude"),
        ("amp200", "amplitude x2.00 (alpha 0.982 -> 1.964)", "amplitude")]

snap, SL = board_slopes()
base = read(BASE_221)
ka1 = CO / "ka1-317/base/run/agg_results.csv"

print(f"[317] board scaling map from {snap}")
print("\n=== KA-1: the instrument re-scores a byte-identical candidate ===")
if ka1.exists():
    got = read(ka1)
    print(f"{'member':<6} {'221 recorded':>18} {'317 re-score':>18} {'|diff|':>10}")
    mx = 0.0
    for m in MEMBERS:
        dfr = abs(got[m] - base[m]); mx = max(mx, dfr)
        print(f"{m:<6} {base[m]:18.12f} {got[m]:18.12f} {dfr:10.2e}")
    print(f"KA-1 {'PASS' if mx < 1e-12 else 'FAIL'} — max |diff| {mx:.2e}")
else:
    print("  (not collected yet)")

a15, a15s = read(TWIN[0]), read(TWIN[1])
print("\n=== noise floor: the re-seeded twin a15 vs a15_seed1 (ONE pair, not an sd) ===")
tw = {m: a15s[m] - a15[m] for m in MEMBERS}
tw_ov = sum(tw[m] * SL[m]["slope"] for m in MEMBERS) / 6
print("  " + "  ".join(f"{m} {tw[m]:+.6f}" for m in MEMBERS))
print(f"  priced at transfer factor 1.0: {tw_ov:+.6f} overall")

print("\n=== the sweep: RAW six members on h1-test, one variable per row ===")
hdr = f"{'arm':<9} {'variable':<38} " + " ".join(f"{m:>9}" for m in MEMBERS)
print(hdr); print("-" * len(hdr))
print(f"{'base':<9} {'shipped (= 221 rc_callset)':<38} " +
      " ".join(f"{base[m]:9.6f}" for m in MEMBERS))
rows = {}
for name, desc, cls in ARMS:
    p = CO / f"s317-{name}/{name}/run/agg_results.csv"
    if not p.exists():
        print(f"{name:<9} {desc:<38} " + " ".join(f"{'--':>9}" for _ in MEMBERS)); continue
    v = read(p); rows[name] = (v, cls, desc)
    print(f"{name:<9} {desc:<38} " + " ".join(f"{v[m]:9.6f}" for m in MEMBERS))

print("\n=== deltas vs base, and the board price (transfer factor 1.0 = optimistic) ===")
hdr2 = (f"{'arm':<9} {'class':<10} " + " ".join(f"d{m:>8}" for m in MEMBERS) +
        f" {'d overall':>10} {'x0.4277':>9}")
print(hdr2); print("-" * len(hdr2))
out = {}
for name, desc, cls in ARMS:
    if name not in rows:
        continue
    v = rows[name][0]
    d = {m: v[m] - base[m] for m in MEMBERS}
    ov = sum(d[m] * SL[m]["slope"] for m in MEMBERS) / 6
    out[name] = {"raw": v, "delta_raw": d, "overall_factor1": ov,
                 "overall_factor0.4277": ov * 0.4277, "class": cls, "variable": desc}
    print(f"{name:<9} {cls:<10} " + " ".join(f"{d[m]:+9.6f}" for m in MEMBERS) +
          f" {ov:+10.6f} {ov*0.4277:+9.6f}")

(HERE / "work").mkdir(exist_ok=True)
(HERE / "work/sweep317.json").write_text(json.dumps(
    {"snapshot": snap, "base_raw": base, "twin_delta": tw, "arms": out}, indent=1))
print("\nwrote work/sweep317.json")


# ---------------------------------------------------------------- paired per-target intervals
print("\n=== paired per-TARGET intervals vs base (unit: one perturbation, n<=85, "
      "10,000-resample percentile bootstrap) ===")
from paired317 import paired  # noqa: E402
BASE_RES = L221 / "ladder8-20260904/rc_callset/run/results.csv"
TW = paired(L221 / "ladder6-20260904/a15/run/results.csv",
            L221 / "ladder7-20260904/a15_seed1/run/results.csv")
print("  NOISE FLOOR (re-seeded twin a15 -> a15_seed1, same construction, one seed pair):")
print("    " + "  ".join(f"{m} {TW[m]['mean']:+.5f}[{TW[m]['lo']:+.5f},{TW[m]['hi']:+.5f}]"
                         for m in ("pds", "nmae", "fid", "reach", "jac")))
pt = {}
for name, desc, cls in ARMS:
    p = CO / f"s317-{name}/{name}/run/results.csv"
    if not p.exists():
        continue
    r = paired(BASE_RES, p)
    pt[name] = r
    print(f"\n  {name} — {desc}")
    for m in ("pds", "nmae", "fid", "reach", "jac"):
        if m not in r:
            continue
        v = r[m]
        star = "" if (v["lo"] <= 0 <= v["hi"]) else "  *"
        print(f"    {m:<6} {v['mean']:+.6f} [{v['lo']:+.6f}, {v['hi']:+.6f}]  n={v['n']:2d}  "
              f"up/down {v['n_up']}/{v['n_dn']}{star}")
j = json.loads((HERE / "work/sweep317.json").read_text())
j["paired"] = pt
j["twin_paired"] = TW
(HERE / "work/sweep317.json").write_text(json.dumps(j, indent=1))
print("\n(* = the 95% paired interval excludes zero)")
