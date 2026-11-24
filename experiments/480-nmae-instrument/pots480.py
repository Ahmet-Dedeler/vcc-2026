"""480 resume-2: the magnitude and sign pots per surface, summarised from decomp480.json / decomp480b.json (no new model).

For each surface and construction: nmae; the magnitude pot (oracle |l| at the construction's signs), the sign pot
(oracle signs at its magnitudes) and the abstain pot (zero exactly where its sign is wrong), each as the nmae it would
recover; the share of the magnitude pot on truth-down genes and in small gates (n_real < 100), and the share of the sign
pot in low/mid-expressed genes. Class contributions add exactly to the totals (decomp480.py asserts it).

    uv run python pots480.py     # -> work/pots480.json
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANCHOR = {"h1": 0.383, "h1_thin": 0.517, "orion_hct116": 0.427, "orion_hek293t": 0.440, "k562_gwps": 0.517}


def tot(rows, key):
    return sum(r[key] for r in rows)


def pick(rows, names, key):
    """Sum over classes whose label starts with one of names (decomp480 labels carry units, e.g. 'low (<0.87/cell)')."""
    return sum(r[key] for r in rows if any(r["cls"].split(" ")[0] == n for n in names))


def share(x, total):
    return x / total if abs(total) > 0.01 else float("nan")          # a pot of ~0 has no meaningful split


out = {}
a = json.load(open(HERE / "work/decomp480.json"))["h1_ours"]
b = json.load(open(HERE / "work/decomp480b.json"))
surfaces = [("h1 (ours, shipped)", a, "h1")] + [(f"{s} (a_analog)", b[s], s) for s in b]
print(f"{'surface':28s} {'anchor':>6s} {'nmae':>6s} {'mag':>7s} {'sign':>7s} {'abst':>7s} {'mag:down':>8s} "
      f"{'mag:small':>9s} {'sign:lowmid':>11s}")
for label, d, s in surfaces:
    tr = d["truth direction"]
    mag, sign, abst = tot(tr, "rec_mag"), tot(tr, "rec_sign"), tot(tr, "rec_abst")
    gs = d["gate size n_real (oracle)"]
    small = pick(gs, {"10-99", "10-29", "30-99"}, "rec_mag")
    ex = d["control expression tertile"]
    lowmid = pick(ex, {"low", "mid"}, "rec_sign")
    row = dict(anchor=ANCHOR[s], nmae=d["_total"], mag=mag, sign=sign, abst=abst,
               mag_down_share=share(pick(tr, {"down"}, "rec_mag"), mag),
               mag_small_gate_share=share(small, mag),
               sign_lowmid_share=share(lowmid, sign))
    out[label] = row
    print(f"{label:28s} {row['anchor']:6.3f} {row['nmae']:6.4f} {mag:+7.4f} {sign:+7.4f} {abst:+7.4f} "
          f"{row['mag_down_share']:8.2f} {row['mag_small_gate_share']:9.2f} {row['sign_lowmid_share']:11.2f}")
json.dump(out, open(HERE / "work/pots480.json", "w"), indent=1)
