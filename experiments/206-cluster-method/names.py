"""206 — mine every model_name, organization and description ever observed, across all snapshots.

Targets: the `pm-r*` series (Prima Mente?), the tree series (jw_13288), the bird series (MT Zeng),
the `AL-####` series (Andrew Liu), and any org string mentioning Oxford / Prima Mente.

Read-only. Known-answer gate: `pm-r14d` must be attributed to Hongjian Zhou, whose org string must
read 'Oxford / Prima Mente'.
"""

import json
import glob
import os
import re
import sys
from collections import defaultdict

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments/206-cluster-method")

paths = sorted(glob.glob(os.path.join(ROOT, "data/leaderboard/*.json")))
paths.append(os.path.join(ROOT, "experiments/186-intel-0903/board/leaderboard_2026-09-03.json"))
paths.append(os.path.join(HERE, "board/leaderboard_live.json"))

seen = {}
for p in paths:
    with open(p) as fh:
        d = json.load(fh)
    for e in (d["entries"] if isinstance(d, dict) else d):
        key = (e["team_name"], e.get("submission_count"))
        if key not in seen:
            seen[key] = e

rows = sorted(seen.values(), key=lambda e: e["submission_date"])
print(f"distinct observations: {len(rows)}  distinct teams: {len({e['team_name'] for e in rows})}")

by_model = defaultdict(set)
for e in rows:
    by_model[e.get("model_name")].add(e["team_name"])

gate = by_model.get("pm-r14d")
hz = [e for e in rows if e["team_name"] == "Hongjian Zhou"][-1]
if gate != {"Hongjian Zhou"} or hz.get("organization") != "Oxford / Prima Mente":
    print(f"GATE FAIL: pm-r14d -> {gate}, org {hz.get('organization')!r}", file=sys.stderr)
    sys.exit(1)
print("GATE PASS: pm-r14d -> Hongjian Zhou, org 'Oxford / Prima Mente'\n")

# ---- 1. the `pm-` namespace, board-wide
print("=== EVERY MODEL NAME MATCHING /^pm[-_]/i, ANY TEAM, ANY SNAPSHOT ===")
for e in rows:
    m = str(e.get("model_name") or "")
    if re.match(r"^pm[-_]", m, re.I):
        print(f"  {e['submission_date'][:19]}  {e['team_name']:<18}#{e.get('submission_count'):<3}"
              f"{m:<24} ovr {e['score_avg']:>9.5f}  org={e.get('organization')!r}")

# ---- 2. every org string that mentions Oxford / Prima / Mente
print("\n=== ORG STRINGS MENTIONING oxford|prima|mente|clifton ===")
orgs = defaultdict(set)
for e in rows:
    orgs[(e.get("organization") or "").strip()].add(e["team_name"])
for o, ts in sorted(orgs.items()):
    if re.search(r"oxford|prima|mente|clifton", o, re.I):
        print(f"  {o!r:<40} -> {sorted(ts)}")

# ---- 3. the four accounts' full observed name series
print("\n=== THE CLUSTER'S FULL OBSERVED MODEL-NAME SERIES ===")
CL = ["Hongjian Zhou", "jw_13288", "Andrew Liu", "MT Zeng", "cqawesome"]
for t in CL:
    es = [e for e in rows if e["team_name"] == t]
    total = max(e.get("submission_count") or 0 for e in es)
    names = [(e.get("submission_count"), e.get("model_name")) for e in es]
    print(f"  {t:<16} ({len(es)}/{total} seen, org={es[-1].get('organization')!r}, "
          f"members={es[-1].get('team_members')})")
    print(f"      {names}")

# ---- 4. teams sharing a naming convention with any cluster member
print("\n=== OTHER TEAMS USING THE CLUSTER'S NAMING CONVENTIONS ===")
PATS = {
    "trees (jw_13288)": r"^(willow|hazel|cedar|juniper|birch|maple|alder|aspen|rowan|elm|oak|pine|fir|larch|spruce|yew|holly|ash|beech)[-_]?\d*$",
    "birds (MT Zeng)": r"^(gannet|kestrel|petrel|falcon|osprey|merlin|harrier|shrike|swift|tern|skua|curlew|plover|godwit)\d*$",
    "AL-#### (Andrew Liu)": r"^AL-\d{3,4}$",
    "m####/m0903x (cqawesome)": r"^m\d{2,4}[a-z]?$",
}
for label, pat in PATS.items():
    hits = defaultdict(list)
    for e in rows:
        if re.match(pat, str(e.get("model_name") or ""), re.I):
            hits[e["team_name"]].append(e.get("model_name"))
    print(f"  {label}:")
    for t, ms in sorted(hits.items(), key=lambda kv: -len(kv[1])):
        print(f"      {t:<24} {ms}")

# ---- 5. every non-empty description on the LIVE board, ranked by score, for the top 60 teams
live = json.load(open(os.path.join(HERE, "board/leaderboard_live.json")))["entries"]
live.sort(key=lambda e: -(e["score_avg"] or -9))
print("\n=== DESCRIPTIONS ON THE LIVE BOARD, TOP 60 BY SCORE ===")
for e in live[:60]:
    d = e.get("description")
    tag = "" if d else "  (empty)"
    print(f"\n  #{live.index(e)+1:<3} {e['team_name']:<20} {str(e.get('model_name'))[:22]:<24}"
          f"ovr {e['score_avg']:.5f} mse {e['expr_mse_unbiased_capped_norm']:.4f}"
          f"  org={e.get('organization')!r}{tag}")
    if d:
        print("        " + str(d).replace("\n", "\n        ")[:1400])

# ---- 6. any description board-wide naming an external dataset
print("\n\n=== DESCRIPTIONS (ANY SNAPSHOT) NAMING A CANDIDATE EXTERNAL ATLAS ===")
KEY = r"x-atlas|xatlas|pisces|x-cell|xcell|tahoe|scbasecamp|perturb-?base|cellxgene|" \
      r"arc atlas|virtual cell atlas|orion|replogle|gwps|feng|kolf|jurkat|hepg2|rpe1|" \
      r"in-?house|proprietary|internal data|our own (data|screen)|new (atlas|dataset)"
hit = {}
for e in rows:
    d = str(e.get("description") or "")
    if d and re.search(KEY, d, re.I):
        hit.setdefault(e["team_name"], e)
for t, e in sorted(hit.items(), key=lambda kv: -(kv[1]["score_avg"] or -9)):
    print(f"\n  {t:<20} {str(e.get('model_name'))[:22]:<24} ovr {e['score_avg']:.5f} "
          f"mse {e['expr_mse_unbiased_capped_norm']:.4f}")
    print("      " + str(e["description"]).replace("\n", "\n      ")[:900])
