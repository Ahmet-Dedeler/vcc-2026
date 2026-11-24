"""206 — the high-cosine cohort: who is above us, are they independent of each other, and what does
159's `1/k` law say each of them needs.

The discriminating logic: if ONE party holds a private atlas, exactly one pipeline sits far above
the public-line ceiling. If SEVERAL mutually-independent pipelines sit there, the route is
reachable from public material.

Read-only. Known-answer gates:
  KA-1  197 §2 a15 row: floor(0.87381) = 0.07221 and k_eff = 3.09.
  KA-2  197 §4.2: the best equal consensus of all four panel-covering lines is board-equivalent
        0.0873, which must return k_eff = 4.00 under the same anchor.
"""

import json
import glob
import math
import os
import sys
from collections import defaultdict

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments/206-cluster-method")

NUM, DEN = 1.023611, 2.074495
# 197 §4.1 board<->surface anchors. It reports k_eff as a BAND over these three, not a point.
ANCHORS = {"a15": 2.078, "a23": 2.041, "abc-knee": 1.954}
LAW_C, LAW_B = 66.5, 1920.0
cos_floor = lambda m: (NUM - m) / DEN


def k_from_surface(s):
    d = 1.0 / (s * s) - LAW_C
    return LAW_B / d if d > 0 else float("inf")


def k_band(board_cos):
    ks = [k_from_surface(board_cos / a) for a in ANCHORS.values()]
    return min(ks), max(ks)


def k_eff(board_cos, anchor="a15"):
    return k_from_surface(board_cos / ANCHORS[anchor])


# KA-1: 197 §2's a15 row -> floor 0.07221; §4.2's production surface cosine 0.03814 -> k 3.09
f = cos_floor(0.8738045054832962)
if abs(f - 0.07221) > 5e-6 or abs(k_from_surface(0.03814) - 3.09) > 0.02:
    print(f"GATE FAIL KA-1: floor {f:.5f}, k(0.03814) {k_from_surface(0.03814):.2f}",
          file=sys.stderr)
    sys.exit(1)
# KA-2: 197 §4.2's four-line consensus, surface 0.04278 -> k_eff 4.00
if abs(k_from_surface(0.04278) - 4.00) > 0.02:
    print(f"GATE FAIL KA-2: k(0.04278) = {k_from_surface(0.04278):.2f}, want 4.00", file=sys.stderr)
    sys.exit(1)
print(f"GATES PASS  floor(a15)={f:.5f}  k(production surface 0.03814)="
      f"{k_from_surface(0.03814):.2f}  k(four-line surface 0.04278)="
      f"{k_from_surface(0.04278):.2f}\n")

MEM = ["pds_cosine", "expr_mse_unbiased_capped_norm", "de_wilcoxon_lfc_nmae",
       "de_wilcoxon_direction_fidelity_yield_raw", "de_wilcoxon_direction_reach_raw",
       "de_wilcoxon_sig_jaccard"]
LAB = ["pds", "mse", "nmae", "fid", "reach", "jac"]

live = json.load(open(os.path.join(HERE, "board/leaderboard_live.json")))["entries"]
us = [e for e in live if e["team_name"] == "Ahmet Dedeler"][0]

# ---- the "like us" point estimate: both at claim = K and on the cliff, so
#      Δcos = −Δraw / 2.074495 exactly.
print("=== COSINE: FLOOR vs THE 'LIKE US' POINT ESTIMATE ===")
print("   (we ship claim at 100.0% of the PRED_TRACE_CAP_K cap — 197 §1 rc_hybrid row — so a team")
print("    sitting where we sit inverts to cos_us + (raw_us - raw_them)/2.074495)\n")
print(f"{'team':<22}{'model':<16}{'raw mse':>9}{'floor':>9}{'like-us':>9}{'xours':>7}"
      f"{'k(flr)':>8}{'k(like)':>8}")
COS_US = 0.079249
hi = sorted([e for e in live if cos_floor(e[MEM[1]]) > 0.09], key=lambda e: -cos_floor(e[MEM[1]]))
for e in hi:
    fl = cos_floor(e[MEM[1]])
    lk = COS_US + (us[MEM[1]] - e[MEM[1]]) / DEN
    print(f"{e['team_name']:<22}{str(e['model_name'])[:15]:<16}{e[MEM[1]]:>9.4f}{fl:>9.5f}"
          f"{lk:>9.5f}{lk/COS_US:>7.2f}{k_eff(fl):>8.1f}{k_eff(lk):>8.1f}")
print(f"{'US':<22}{'rc_hybrid':<16}{us[MEM[1]]:>9.4f}{cos_floor(us[MEM[1]]):>9.5f}"
      f"{COS_US:>9.5f}{1.0:>7.2f}{k_eff(cos_floor(us[MEM[1]])):>8.1f}{k_eff(COS_US):>8.1f}")

# ---- independence: pairwise six-tuple distance, normalised by each member's board sd
allv = {m: [e[m] for e in live if e[m] is not None] for m in MEM}
sd = {}
for m in MEM:
    v = sorted(allv[m])
    v = v[len(v)//10: -len(v)//10 or None]          # trimmed, the board has huge mse outliers
    mu = sum(v) / len(v)
    sd[m] = math.sqrt(sum((x - mu) ** 2 for x in v) / len(v))

TOP = ["cqawesome", "Hongjian Zhou", "jw_13288", "Andrew Liu", "MT Zeng", "Garrrbage",
       "Aginglab.com", "Ahmet Dedeler"]
E = {e["team_name"]: e for e in live if e["team_name"] in TOP}
print("\n=== PAIRWISE SIX-TUPLE DISTANCE (root-mean-square of per-member z gaps) ===")
print("    small = the same prediction; large = independent pipelines")
print(f"{'':<16}" + "".join(f"{t[:12]:>13}" for t in TOP))
for a in TOP:
    line = f"{a[:15]:<16}"
    for b in TOP:
        if a == b:
            line += f"{'-':>13}"
            continue
        d = math.sqrt(sum(((E[a][m] - E[b][m]) / sd[m]) ** 2 for m in MEM) / 6)
        line += f"{d:>13.3f}"
    print(line)

# ---- what the 1/k law demands
print("\n=== 159's 1/k LAW: WHAT EACH COSINE COSTS IN EXCHANGEABLE DONOR LINES ===")
print("    k_eff is reported as a BAND over 197 §4.1's three board<->surface anchors")
print(f"{'quantity':<52}{'board cos':>11}{'k_eff band':>16}")
for lab, c in [("production (us, 197's measurement)", COS_US),
               ("best equal consensus of all FOUR public lines (197)", 0.0873),
               ("Garrrbage v14 (floor)", cos_floor(0.7858)),
               ("cluster floor (jw_13288 juniper_46)", cos_floor(0.7515)),
               ("cluster, like-us estimate", COS_US + (us[MEM[1]] - 0.7515) / DEN),
               ("cqawesome m0903c (floor)", cos_floor(0.7495)),
               ("cqawesome, like-us estimate", COS_US + (us[MEM[1]] - 0.7495) / DEN),
               ("context-independent k -> inf ceiling", 0.2548)]:
    lo, hi_ = k_band(c)
    print(f"{lab:<52}{c:>11.5f}{f'{lo:.1f} - {hi_:.1f}':>16}")

# ---- what each hypothesis is worth in overall units
print("\n=== 197's EXCHANGE RATE: +0.00365 OVERALL PER +0.01 COSINE ===")
RATE = 0.00365
for lab, dcos in [("+1 exchangeable public line (3.09 -> 4.00 lines)", 0.0873 - COS_US),
                  ("close HALF the gap to the cluster floor", (0.13117 - COS_US) / 2),
                  ("close the whole gap to the cluster floor", 0.13117 - COS_US),
                  ("close the gap to the like-us cluster estimate",
                   COS_US + (us[MEM[1]] - 0.7515) / DEN - COS_US),
                  ("GISL's published board gain from Marson CD4-T", 0.01325 / RATE / 100)]:
    print(f"  {lab:<52}{dcos:>+9.5f} cos  ->  {dcos * 100 * RATE:>+9.5f} overall")
print(f"\n  gap to 1st (cqawesome 0.22082) from us ({us['score_avg']:.5f}) = "
      f"{0.22082 - us['score_avg']:+.5f} overall")
