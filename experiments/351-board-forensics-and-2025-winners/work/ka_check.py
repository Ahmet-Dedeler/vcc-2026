"""351 Part A: board forensics. Stdlib only. Reads data/leaderboard/*.json."""
import json, glob, csv, os
from datetime import datetime

REPO = "~/Code/virtual-cell-challenge-2026"
LB = os.path.join(REPO, "data/leaderboard/*.json")
WORK = os.path.join(REPO, "experiments/351-board-forensics-and-2025-winners/work")
os.makedirs(WORK, exist_ok=True)

MEMBERS = ["score_pds", "score_mse", "score_nmae", "score_fid", "score_reach", "score_jac"]
RAW = {"score_pds": "pds_cosine", "score_mse": "expr_mse_unbiased_capped_norm",
       "score_nmae": "de_wilcoxon_lfc_nmae", "score_fid": "de_wilcoxon_direction_fidelity_yield_raw",
       "score_reach": "de_wilcoxon_direction_reach_raw", "score_jac": "de_wilcoxon_sig_jaccard"}

files = sorted(glob.glob(LB))
print(f"snapshots: {len(files)}")
print(f"first: {os.path.basename(files[0])} last: {os.path.basename(files[-1])}")

# Dedupe by id
seen = {}  # id -> record (latest snapshot's version), first_seen, last_seen, n_seen
for f in files:
    d = json.load(open(f))
    snap = os.path.basename(f)
    for e in d["entries"]:
        i = e["id"]
        if i not in seen:
            seen[i] = {"rec": e, "first": snap, "last": snap, "n": 1}
        else:
            seen[i]["rec"] = e
            seen[i]["last"] = snap
            seen[i]["n"] += 1

print(f"distinct ids: {len(seen)}")
latest = json.load(open(files[-1]))
print(f"latest num_entries: {latest['num_entries']} len: {len(latest['entries'])}")
latest_ids = {e["id"] for e in latest["entries"]}
missing = latest_ids - set(seen.keys())
print(f"latest ids missing from dedup: {len(missing)}")

# ---- Known-answer checks ----
ka = json.load(open(os.path.join(REPO, "data/board-entries/xuDelEik7zGWMXdNu0R4.json")))["entry"]
r319 = seen["xuDelEik7zGWMXdNu0R4"]["rec"]
print("\nKNOWN-A (319):")
for k in ["score_avg"] + MEMBERS:
    a, b = ka[k], r319[k]
    print(f"  {k}: board-entry={a!r} snapshot={b!r} diff={abs(a-b):.2e}")
m6 = sum(r319[m] for m in MEMBERS) / 6
print(f"  mean(six)={m6:.16f} score_avg={r319['score_avg']:.16f} diff={abs(m6-r319['score_avg']):.2e}")

m323 = json.load(open(os.path.join(REPO, "experiments/323-current-family-validation/work/metrics-cap15.json")))
print(f"\nKNOWN-B (323 raw.nmae): {m323['raw']['nmae']:.16f} (expect 0.9508364630325662)")
print(f"KNOWN-B (323 raw.reach): {m323['raw']['reach']:.16f} (expect 0.25434726355678233)")

# Leader vs us gap
teams = {}
for i, v in seen.items():
    teams.setdefault(v["rec"]["team_name"], []).append((v["rec"], i))
best = {}
for t, lst in teams.items():
    # best by score_avg using latest record per id; but entries of a team across snapshots:
    # best = max over ids' latest rec
    b = max(lst, key=lambda x: (x[0]["score_avg"] is not None, x[0]["score_avg"] or -9))
    best[t] = b
ranked = sorted(best.items(), key=lambda kv: -(kv[1][0]["score_avg"] or -9))
print(f"\nteams: {len(teams)}")
for rank, (t, (r, i)) in enumerate(ranked[:15], 1):
    print(f"  {rank:2d} {r['score_avg']:.6f} {t[:34]:34s} {r.get('model_name','')[:28]} {r['submission_date'][:10]} {i[:12]}")
ah = [x for x in ranked if x[0] == "Ahmet Dedeler"]
print("us:", [(ranked.index(x) + 1, x[1][0]["score_avg"]) for x in ah])
lv = [x for x in ranked if x[0] == "ViRtUaL cElL builders"]
print("leader:", [(ranked.index(x) + 1, x[1][0]["score_avg"]) for x in lv])

lr, _ = lv[0][1]
ur, _ = ah[0][1]
print("\nKNOWN-C gap leader-minus-us:")
tot = 0
for m in MEMBERS:
    c = (lr[m] - ur[m]) / 6
    tot += c
    print(f"  {m}: {c:+.8f}")
print(f"  sum={tot:.10f} actual diff={lr['score_avg']-ur['score_avg']:.10f}")

json.dump({"n_snapshots": len(files), "distinct": len(seen),
           "latest_n": len(latest_ids), "missing": len(missing)},
          open(os.path.join(WORK, "build_meta.json"), "w"), indent=2)
