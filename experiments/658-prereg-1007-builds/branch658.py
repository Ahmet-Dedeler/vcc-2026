#!/usr/bin/env python3
"""658b: classify 524dr2's resolved recipe against 524dr's (PREREGISTRATION-524dr2.md 1 and 3). Prints ONE line:
  N0 none | N1 mde_sign | N2 tercilec[,h1_dose(<ENV>)] | N3 <bulk...>,mde_sign | STOP <reason>
A dropped switch, changed def_day arguments, or any change other than tercile -> tercilec (same numbers), an added
DEF_DAY_MDE_SIGN=clean4:0.05, or an added env key naming H1 (653 R10's dose, once wired) is a STOP: amend 658b first.
usage: uv run python experiments/658-prereg-1007-builds/branch658.py NEW_RECIPE.json PREV_RECIPE.json"""
import json
import sys


def classify(new: dict, old: dict) -> str:
    if new.get("problems"):
        return "STOP the recipe has problems " + json.dumps(new["problems"])
    ne, oe = new["env"], old["env"]
    if new["args"][3:] != old["args"][3:]:
        return f"STOP def_day arguments changed {old['args'][3:]} -> {new['args'][3:]}"
    dropped = sorted(set(oe) - set(ne))
    if dropped:
        return f"STOP switches dropped since 524dr {dropped}"
    bulk, de, other = [], [], []
    for k in sorted(k for k in oe if k in ne and oe[k] != ne[k]):
        if k == "DEF_DAY_SHARPEN_PT" and ne[k] == oe[k].replace("tercile:", "tercilec:", 1) and oe[k].startswith("tercile:"):
            bulk.append("tercilec")
        else:
            other.append(f"{k}={oe[k]}->{ne[k]}")
    for k in sorted(set(ne) - set(oe)):
        if k == "DEF_DAY_MDE_SIGN" and ne[k] == "clean4:0.05":
            de.append("mde_sign")
        elif "H1" in k:
            bulk.append(f"h1_dose({k})")
        else:
            other.append(f"{k}={ne[k]}")
    if other:
        return "STOP unwired recipe changes " + ", ".join(other) + ": amend 658b first"
    br = "N3" if bulk and de else "N2" if bulk else "N1" if de else "N0"
    return f"{br} {','.join(bulk + de) or 'none'}"


if __name__ == "__main__":
    print(classify(*(json.load(open(p)) for p in sys.argv[1:3])))
