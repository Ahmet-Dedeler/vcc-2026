"""Q31 — is the mse/reach trade a frontier, and who is off it?

Q30's "already tuned" split was mislabelled: raw mse ~1.0-1.2 is the BASELINE (0.989), i.e.
scoring ~0 on the member, not tuned.  And in that group every large reach gain came with mse
getting WORSE.  So the board's normal behaviour is: buy reach by over-dispersing past the point
where mse is collectable.  Teams at mse ~1.0 have nothing to lose, so they push hard.

We are at raw mse 0.8828 -- BELOW baseline, actually collecting credit (+0.120 scaled), which
almost nobody does.  The real question is who holds BOTH.
"""
import json
d = json.load(open("/tmp/lb_now.json"))
rows = [r for r in (d if isinstance(d, list) else d.get("entries") or [])
        if isinstance(r, dict) and r.get("score_avg") is not None]
M, R = "expr_mse_unbiased_capped_norm", "de_wilcoxon_direction_reach_raw"
both = [r for r in rows if r[M] < 1.00 and r[R] > 0.20]
mse_only = [r for r in rows if r[M] < 1.00 and r[R] <= 0.20]
print(f"{len(rows)} teams live\n")
print(f"COLLECT mse (raw < 1.00) AND reach > 0.20 — the club we are not in: {len(both)}")
print(f"{'team':24} {'avg':>8} {'raw mse':>8} {'raw reach':>10} {'raw nmae':>9} {'pds':>7} {'subs':>5}")
for r in sorted(both, key=lambda r: -r["score_avg"]):
    print(f"{r['team_name'][:24]:24} {r['score_avg']:8.5f} {r[M]:8.4f} {r[R]:10.4f} "
          f"{r['de_wilcoxon_lfc_nmae']:9.4f} {r['pds_cosine']:7.4f} {r.get('submission_count',0):5}")
print(f"\nCOLLECT mse but reach <= 0.20 — where we sit: {len(mse_only)}")
for r in sorted(mse_only, key=lambda r: -r["score_avg"])[:8]:
    mark = "  <- US" if r["team_name"] == "Ahmet Dedeler" else ""
    print(f"{r['team_name'][:24]:24} {r['score_avg']:8.5f} {r[M]:8.4f} {r[R]:10.4f} "
          f"{r['de_wilcoxon_lfc_nmae']:9.4f} {r['pds_cosine']:7.4f} {r.get('submission_count',0):5}{mark}")
top = sorted(rows, key=lambda r: -r["score_avg"])[:30]
print(f"\nAmong the top 30 by overall: {sum(1 for r in top if r[M] < 1.0)} collect mse, "
      f"{sum(1 for r in top if r[R] > 0.20)} have reach > 0.20, "
      f"{sum(1 for r in top if r[M] < 1.0 and r[R] > 0.20)} have BOTH.")
