"""317 KA-5 — the board scaling map retrodicts OUR OWN scored board rows.

`docs/BOARD-LEDGER.md` records raw AND scaled members for five of our submitted arms. The map
fitted in `price317.py` (from other teams' rows, on each member's unclamped core) is applied to
those raw values and compared to the scaled values the board itself published. Nothing here is
fitted on our rows, so this is an out-of-sample retrodiction of the pricing instrument.
"""
import json
from pathlib import Path
from price317 import board_slopes, MEMBERS

# BOARD-LEDGER.md §"Raw members, where recorded" x §"Scaled members"
ROWS = {
    "a11":      {"raw": {"pds": 0.833307, "mse": 0.909989, "reach": 0.178243},
                 "scaled": {"pds": 0.739889, "mse": 0.082041, "reach": 0.112220}},
    "a14":      {"raw": {"pds": 0.832702, "mse": 0.882775, "nmae": 0.945764,
                         "fid": 0.502881, "reach": 0.180615, "jac": 0.030081},
                 "scaled": {"pds": 0.738614, "mse": 0.110615, "nmae": 0.091972,
                            "fid": -0.032586, "reach": 0.114542, "jac": -0.000612}},
    "a15":      {"raw": {"pds": 0.835987, "mse": 0.873805, "nmae": 0.944251,
                         "fid": 0.503536, "reach": 0.181173, "jac": 0.030080},
                 "scaled": {"pds": 0.745871, "mse": 0.120057, "nmae": 0.094239,
                            "fid": -0.030397, "reach": 0.115368, "jac": -0.000619}},
    "abc-knee": {"raw": {"pds": 0.836786, "mse": 0.874491, "nmae": 0.943171,
                         "fid": 0.503819, "reach": 0.177303, "jac": 0.030145},
                 "scaled": {"pds": 0.747726, "mse": 0.119322, "nmae": 0.096214,
                            "fid": -0.029477, "reach": 0.110910, "jac": -0.000433}},
}
# cd4svd5 and the live entry, from docs/READ-FIRST.md §1 / §1.0 and the snapshot
LIVE = {"cd4svd5": {"raw": {"pds": 0.85064, "mse": 0.82181, "nmae": 0.93397,
                            "fid": 0.50739, "reach": 0.27506, "jac": 0.04090},
                    "scaled": {"pds": 0.77809, "mse": 0.17477, "nmae": 0.11047,
                               "fid": -0.01683, "reach": 0.22139, "jac": 0.02868}}}

name, sl = board_slopes()
print(f"[KA-5] board scaling map from {name}, applied to OUR OWN recorded board rows\n")
print(f"{'arm':<10} {'member':<6} {'raw':>10} {'published':>11} {'retrodicted':>12} {'err':>10}")
worst = {}
for arm, d in list(ROWS.items()) + list(LIVE.items()):
    for m in MEMBERS:
        if m not in d["raw"]:
            continue
        s = sl[m]
        got = (d["raw"][m] - s["b"]) * s["slope"]
        err = got - d["scaled"][m]
        worst[m] = max(worst.get(m, 0.0), abs(err))
        print(f"{arm:<10} {m:<6} {d['raw'][m]:10.6f} {d['scaled'][m]:11.6f} "
              f"{got:12.6f} {err:+10.2e}")
print("\nworst |err| per member: " + "  ".join(f"{m} {worst.get(m, float('nan')):.1e}"
                                               for m in MEMBERS))
print(f"worst overall-point error if all six erred with the same sign: "
      f"{sum(worst.values())/6:.5f}")
