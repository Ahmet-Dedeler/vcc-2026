"""317 — the ceiling: what the amplitude axis is worth if `mse` were FREE.

The leader's shape has `score_mse` at (near) the clamp floor, i.e. they pay nothing more for
losing `mse`. So the most generous reading of the hypothesis is: PRETEND `mse` costs nothing
and ask what the other five members are worth over the whole amplitude axis. That is an upper
bound on the trade by construction, and it is the number the verdict should be argued against.
"""
import json
from pathlib import Path
import numpy as np
from read_agg import read
from price317 import board_slopes, MEMBERS

HERE = Path(__file__).resolve().parent
CO = HERE / "cloud-output"
BASE = read(HERE.parents[1] / "experiments/221-cloud-scorer/cloud-output/"
            "ladder8-20260904/rc_callset/run/agg_results.csv")
_, SL = board_slopes()
AMP = {"amp070": 0.70, "amp085": 0.85, "amp120": 1.20, "amp145": 1.45, "amp200": 2.00}
EMI = ["od1", "od8", "ks010", "ts070", "kt16000"]
FREE = [m for m in MEMBERS if m != "mse"]

rows = {n: read(CO / f"s317-{n}/{n}/run/agg_results.csv") for n in list(AMP) + EMI}
d = {n: {m: rows[n][m] - BASE[m] for m in MEMBERS} for n in rows}

print("[317] AMPLITUDE arms, priced with `mse` COUNTED FREE (the generous reading)")
print(f"{'arm':<9} {'s':>5} " + " ".join(f"{m:>9}" for m in FREE) +
      f" {'overall':>10} {'+ mse':>10}")
for n, s in AMP.items():
    free = sum(d[n][m] * SL[m]["slope"] for m in FREE) / 6
    full = free + d[n]["mse"] * SL["mse"]["slope"] / 6
    print(f"{n:<9} {s:5.2f} " + " ".join(f"{d[n][m]:+9.6f}" for m in FREE) +
          f" {free:+10.6f} {full:+10.6f}")

orc = {m: max(d[n][m] * SL[m]["slope"] for n in AMP) for m in FREE}
print(f"\n  FIVE-POINT ORACLE (each member's best over the five amplitude arms, an upper "
      f"bound that no single operating point can reach):")
for m in FREE:
    best = max(AMP, key=lambda n: d[n][m] * SL[m]["slope"])
    print(f"    {m:<6} best {best:<7} d raw {d[best][m]:+.6f}  -> {orc[m]/6:+.6f} overall")
print(f"    SUM {sum(orc.values())/6:+.6f} overall   vs the +0.0089 gate an "
      f"amplitude-changing arm must clear")

print("\n[317] EMISSION arms — `pds` and `mse` are PINNED, so this IS the full price")
print(f"{'arm':<9} " + " ".join(f"{m:>9}" for m in FREE) + f" {'overall':>10}")
for n in EMI:
    ov = sum(d[n][m] * SL[m]["slope"] for m in MEMBERS) / 6
    print(f"{n:<9} " + " ".join(f"{d[n][m]:+9.6f}" for m in FREE) + f" {ov:+10.6f}")
best_emi = max(EMI, key=lambda n: sum(d[n][m] * SL[m]["slope"] for m in MEMBERS))
bo = sum(d[best_emi][m] * SL[m]["slope"] for m in MEMBERS) / 6
orc_e = sum(max(d[n][m] * SL[m]["slope"] for n in EMI) for m in MEMBERS) / 6
print(f"  best single emission point: {best_emi} at {bo:+.6f} overall "
      f"= {bo/0.0019*100:.0f}% of the +0.0019 shape-only gate")
print(f"  five-point emission ORACLE (upper bound, unreachable): {orc_e:+.6f} overall "
      f"= {orc_e/0.0019*100:.0f}% of the shape-only gate")
print(f"  re-seeded-twin noise floor, same pricing: +0.000096 overall")

print("\n[317] the pot on the `mse` side: our live entry's score_mse is 0.18013, i.e. "
      f"{0.18013/6:+.5f} overall is ALL there is to surrender.")
