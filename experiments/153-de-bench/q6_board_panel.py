"""Q6 — the leaderboard as panel data: 11 snapshots, 474 teams, within-team paired deltas.

The board is the only instrument in this project that reads the real objective on the real
contexts.  Each snapshot carries every team's LATEST entry, so consecutive snapshots give
WITHIN-TEAM paired changes -- a natural experiment where the treatment is described in the
team's own `description` field.  094 did this once for amplitude (288 pairs) and it produced
the mse finding that changed the project.  It has never been done for `reach` or `nmae`.
"""
import json, glob, math
from collections import defaultdict

M = ["pds_cosine", "expr_mse_unbiased_capped_norm", "de_wilcoxon_lfc_nmae",
     "de_wilcoxon_direction_fidelity_yield_raw", "de_wilcoxon_direction_reach_raw",
     "de_wilcoxon_sig_jaccard"]
S = ["score_pds", "score_mse", "score_nmae", "score_fid", "score_reach", "score_jac"]
SHORT = ["pds", "mse", "nmae", "fid", "reach", "jac"]

hist = defaultdict(dict)          # team -> entry_id -> row
for f in sorted(glob.glob("../../data/leaderboard/*.json")):
    d = json.load(open(f))
    for r in (d if isinstance(d, list) else d.get("entries") or []):
        if isinstance(r, dict) and r.get("score_avg") is not None:
            hist[r["team_name"]][r["id"]] = r

pairs = []                        # (team, prev, cur) consecutive distinct entries
for team, es in hist.items():
    seq = sorted(es.values(), key=lambda r: r.get("submission_date") or "")
    for a, b in zip(seq, seq[1:]):
        if a["id"] != b["id"]:
            pairs.append((team, a, b))
print(f"{len(hist)} teams, {sum(len(v) for v in hist.values())} distinct entries, "
      f"{len(pairs)} within-team consecutive pairs\n")

def corr(x, y):
    n = len(x); mx, my = sum(x)/n, sum(y)/n
    vx = sum((a-mx)**2 for a in x); vy = sum((b-my)**2 for b in y)
    return sum((a-mx)*(b-my) for a, b in zip(x, y))/math.sqrt(vx*vy) if vx*vy else float("nan")

print("Which member's movement explains a team's own score movement?  (within-team deltas)")
d_avg = [b["score_avg"] - a["score_avg"] for _, a, b in pairs]
print(f"{'member':>8} {'corr with Δscore_avg':>22} {'mean |Δ| contribution':>23}")
for s, nm in zip(S, SHORT):
    d = [b[s] - a[s] for _, a, b in pairs]
    print(f"{nm:>8} {corr(d, d_avg):22.3f} {sum(abs(x) for x in d)/len(d)/6:23.5f}")

print("\nBiggest within-team REACH gains, with what the team said changed")
pairs.sort(key=lambda p: -(p[2]["de_wilcoxon_direction_reach_raw"] - p[1]["de_wilcoxon_direction_reach_raw"]))
for team, a, b in pairs[:10]:
    dr = b["de_wilcoxon_direction_reach_raw"] - a["de_wilcoxon_direction_reach_raw"]
    dp = b["pds_cosine"] - a["pds_cosine"]
    dv = b["score_avg"] - a["score_avg"]
    desc = (b.get("description") or "")[:135].replace("\n", " ")
    print(f"  {team[:17]:17} Δreach{dr:+.4f} Δpds{dp:+.4f} Δavg{dv:+.4f} | {b.get('model_name','')[:20]:20} | {desc}")

print("\nDoes buying reach cost pds?  (within-team, pairs where |Δreach| > 0.01)")
sel = [(b["de_wilcoxon_direction_reach_raw"]-a["de_wilcoxon_direction_reach_raw"],
        b["pds_cosine"]-a["pds_cosine"]) for _, a, b in pairs
       if abs(b["de_wilcoxon_direction_reach_raw"]-a["de_wilcoxon_direction_reach_raw"]) > 0.01]
print(f"  n={len(sel)}  corr(Δreach, Δpds) = {corr([x[0] for x in sel], [x[1] for x in sel]):+.3f}")
up = [p for p in sel if p[0] > 0]
print(f"  of {len(up)} reach-gaining moves, {sum(1 for p in up if p[1] > 0)} also gained pds "
      f"({sum(1 for p in up if p[1] > 0)/max(len(up),1):.1%})")
