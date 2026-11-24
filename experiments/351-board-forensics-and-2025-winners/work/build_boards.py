"""351 Part A full build: entries.csv, top30.csv, improvements, correlations, leader trajectory."""
import json, glob, csv, os, re
from datetime import datetime

REPO = "~/Code/virtual-cell-challenge-2026"
WORK = os.path.join(REPO, "experiments/351-board-forensics-and-2025-winners/work")
MEMBERS = ["score_pds", "score_mse", "score_nmae", "score_fid", "score_reach", "score_jac"]
RAW = {"score_pds": "pds_cosine", "score_mse": "expr_mse_unbiased_capped_norm",
       "score_nmae": "de_wilcoxon_lfc_nmae", "score_fid": "de_wilcoxon_direction_fidelity_yield_raw",
       "score_reach": "de_wilcoxon_direction_reach_raw", "score_jac": "de_wilcoxon_sig_jaccard"}

def snap_time(fn):
    b = os.path.basename(fn)
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})T(\d{2})(\d{2})Z", b)
    if m:
        return datetime(*map(int, m.groups()))
    m = re.match(r"lb_(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})", b)
    if m:
        return datetime(*map(int, m.groups()))
    return datetime.min

files = sorted(glob.glob(os.path.join(REPO, "data/leaderboard/*.json")), key=snap_time)
print("oldest:", os.path.basename(files[0]), "newest:", os.path.basename(files[-1]))
LATEST = files[-1]
print("LATEST =", os.path.basename(LATEST))

seen = {}
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

print(f"distinct={len(seen)}")
latest = json.load(open(LATEST))
print(f"latest entries={len(latest['entries'])} num_entries={latest.get('num_entries')}")

# ---- entries.csv ----
with open(os.path.join(WORK, "entries.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["id", "team_name", "team_id", "model_name", "submission_date", "status",
                "score_avg"] + MEMBERS + [RAW[m] for m in MEMBERS] +
               ["organization", "submission_count", "first_seen_snapshot", "last_seen_snapshot", "n_snapshots_seen"])
    for i, v in sorted(seen.items(), key=lambda kv: (-(kv[1]["rec"].get("score_avg") or -9), kv[0])):
        r = v["rec"]
        w.writerow([i, r.get("team_name"), r.get("team_id"), r.get("model_name"),
                    r.get("submission_date"), r.get("status"), r.get("score_avg")] +
                   [r.get(m) for m in MEMBERS] + [r.get(RAW[m]) for m in MEMBERS] +
                   [r.get("organization"), r.get("submission_count"), v["first"], v["last"], v["n"]])
print("wrote entries.csv")

# ---- top30 from LATEST snapshot ----
le = sorted(latest["entries"], key=lambda e: -(e.get("score_avg") if e.get("score_avg") is not None else -9))[:30]
med = {}
for m in MEMBERS + ["score_avg"]:
    vals = sorted(e[m] for e in le if e.get(m) is not None)
    med[m] = vals[len(vals) // 2]
print("top30 median:", {k: round(v, 6) for k, v in med.items()})
with open(os.path.join(WORK, "top30.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["rank", "team_name", "entry_id", "model_name", "submission_date", "score_avg"] +
               MEMBERS + ["profile"])
    for rank, e in enumerate(le, 1):
        prof = []
        for m in MEMBERS:
            v = e.get(m)
            prof.append(("H" if (v is not None and v >= med[m]) else "L"))
        w.writerow([rank, e["team_name"], e["id"], e.get("model_name"), e.get("submission_date"),
                    e.get("score_avg")] + [e.get(m) for m in MEMBERS] + ["".join(prof)])
print("wrote top30.csv")
# where are we
for rank, e in enumerate(le, 1):
    if e["team_name"] == "Ahmet Dedeler":
        print(f"we are rank {rank} in latest snapshot")
        break
else:
    # find our best rank overall
    all_sorted = sorted(latest["entries"], key=lambda e: -(e.get("score_avg") or -9))
    for rank, e in enumerate(all_sorted, 1):
        if e["team_name"] == "Ahmet Dedeler":
            print(f"we are rank {rank} overall in latest (outside top30)")
            print(json.dumps({k: e[k] for k in ['score_avg'] + MEMBERS}, indent=1))
            break

# ---- consecutive pairs per team (by submission_date, deduped ids) ----
def pdate(s):
    try:
        return datetime.fromisoformat(s.replace("Z", "+00:00"))
    except Exception:
        return None

teams = {}
for i, v in seen.items():
    r = v["rec"]
    if r.get("score_avg") is None:
        continue
    teams.setdefault(r["team_name"], []).append(r)
# dedupe within team by id already (seen); sort by date
pairs = []  # (team, prev, cur)
for t, lst in teams.items():
    lst = sorted(lst, key=lambda e: (pdate(e["submission_date"]) or datetime.min, e["id"]))
    for a, b in zip(lst, lst[1:]):
        if a["id"] == b["id"]:
            continue
        pairs.append((t, a, b))
print(f"teams={len(teams)} consecutive pairs={len(pairs)}")

def delta(a, b, k):
    x, y = a.get(k), b.get(k)
    if x is None or y is None:
        return None
    return y - x

improvements = []
for t, a, b in pairs:
    d = delta(a, b, "score_avg")
    if d is not None and d >= 0.005:
        dm = {m: delta(a, b, m) for m in MEMBERS}
        carry = sorted(dm.items(), key=lambda kv: -((kv[1] or -9) / 6))
        improvements.append((t, a["submission_date"], b["submission_date"], a["id"][:8], b["id"][:8],
                             a.get("model_name"), b.get("model_name"), d, dm))
improvements.sort(key=lambda x: -x[7])
print(f"improvements>=0.005: {len(improvements)}")
with open(os.path.join(WORK, "improvements.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["team", "prev_date", "cur_date", "prev_id8", "cur_id8", "prev_model", "cur_model",
                "d_overall"] + ["d_" + m.replace("score_", "") for m in MEMBERS])
    for t, pd_, cd_, pa, ca, pm, cm, d, dm in improvements:
        w.writerow([t, pd_, cd_, pa, ca, pm, cm, round(d, 8)] + [round(dm[m], 8) if dm[m] is not None else "" for m in MEMBERS])

# correlations over all pairs
import math
def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")

def spearman(xs, ys):
    def ranks(v):
        s = sorted(range(len(v)), key=lambda i: v[i])
        r = [0] * len(v)
        for k, i in enumerate(s):
            r[i] = k
        return r
    return pearson(ranks(xs), ranks(ys))

def corr_table(sub, label):
    def col(k):
        return [(delta(a, b, k[0]), delta(a, b, k[1])) for _, a, b in sub]
    out = {}
    for name, k0, k1 in [("nmae_vs_mse", "score_nmae", "score_mse"),
                         ("reach_vs_pds", "score_reach", "score_pds"),
                         ("nmae_vs_reach", "score_nmae", "score_reach")]:
        pts = [(x, y) for x, y in col((k0, k1)) if x is not None and y is not None]
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        out[name] = {"n": len(pts), "pearson": pearson(xs, ys), "spearman": spearman(xs, ys)}
    print(label, {k: (round(v["pearson"], 4), round(v["spearman"], 4), v["n"]) for k, v in out.items()})
    return out

res = {}
res["all_pairs"] = corr_table(pairs, "ALL")
res["improvers"] = corr_table([(t, a, b) for t, a, b in pairs if (delta(a, b, "score_avg") or -9) > 0], "IMPROVERS")
top30teams = {e["team_name"] for e in le}
res["top30_teams"] = corr_table([(t, a, b) for t, a, b in pairs if t in top30teams], "TOP30")
json.dump(res, open(os.path.join(WORK, "correlations.json"), "w"), indent=2)

# nmae-up pairs: mse behaviour
up = [(t, a, b) for t, a, b in pairs if (delta(a, b, "score_nmae") or -9) > 0]
both_up = sum(1 for t, a, b in up if (delta(a, b, "score_mse") or 0) >= 0)
print(f"pairs with d_nmae>0: {len(up)}, of which d_mse>=0: {both_up} ({both_up/max(1,len(up)):.2f})")
up_top = [(t, a, b) for t, a, b in up if t in top30teams]
both_top = sum(1 for t, a, b in up_top if (delta(a, b, "score_mse") or 0) >= 0)
print(f"top30teams pairs d_nmae>0: {len(up_top)}, d_mse>=0: {both_top} ({both_top/max(1,len(up_top)):.2f})")

# ---- leader trajectory ----
leader = sorted([v["rec"] for v in seen.values() if v["rec"]["team_name"] == "ViRtUaL cElL builders"],
                key=lambda e: (pdate(e["submission_date"]) or datetime.min))
print(f"leader entries: {len(leader)}")
with open(os.path.join(WORK, "leader_trajectory.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["submission_date", "id", "model_name", "score_avg"] + MEMBERS)
    for e in leader:
        w.writerow([e["submission_date"], e["id"], e.get("model_name"), e.get("score_avg")] + [e.get(m) for m in MEMBERS])
for e in leader:
    print(f"  {e['submission_date'][:16]} {e.get('model_name','')[:30]:30s} {e.get('score_avg'):.6f} " +
          " ".join(f"{m[-3:]}={e.get(m):+.4f}" for m in ["score_pds", "score_mse", "score_nmae", "score_fid", "score_reach", "score_jac"]))

# ---- jaccard above ours ----
OUR_JAC = seen["xuDelEik7zGWMXdNu0R4"]["rec"]["score_jac"]
print(f"\nour jac={OUR_JAC}")
above = []
for t, lst in teams.items():
    b = max(lst, key=lambda x: x["score_avg"])
    if b["score_jac"] is not None and b["score_jac"] > OUR_JAC:
        above.append((b["score_jac"], t, b))
above.sort(reverse=True)
print(f"teams with best jac above ours: {len(above)}")
with open(os.path.join(WORK, "jac_above.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["team", "entry_id", "model_name", "submission_date", "score_avg"] + MEMBERS)
    for j, t, b in above[:30]:
        w.writerow([t, b["id"], b.get("model_name"), b["submission_date"], b["score_avg"]] + [b.get(m) for m in MEMBERS])
for j, t, b in above[:15]:
    print(f"  jac={j:.4f} {b['score_avg']:.4f} {t[:30]:30s} " +
          " ".join(f"{m[-3:]}={b.get(m):+.4f}" for m in MEMBERS))
