"""Per-member world rank of our entry from a public leaderboard fetch (read-only)."""
import json, sys
src = sys.argv[1]; out = sys.argv[2]
d = json.load(open(src))
ents = [e for e in d["entries"] if e.get("status") == "published"]
ours_id = "iUziqE1gewfcPKzYUDaB"
ours = next(e for e in ents if e["id"] == ours_id)
ents.sort(key=lambda e: -e["score_avg"])
overall_rank = 1 + [e["id"] for e in ents].index(ours_id)
raw = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
res = {"source": src, "n_published": len(ents), "n_entries_all": len(d["entries"]),
       "our_entry": ours_id, "our_overall": ours["score_avg"], "overall_rank": overall_rank,
       "our_model_name": ours.get("model_name")}
for m, rk in raw.items():
    sk = f"score_{m}"
    vals = sorted((e[sk] for e in ents if e.get(sk) is not None), reverse=True)
    rank = 1 + sum(1 for v in vals if v > ours[sk])
    best = max(ents, key=lambda e: e[sk] if e.get(sk) is not None else -9e9)
    res[m] = {"raw": ours[rk], "scaled": ours[sk], "rank": rank, "world_best_scaled": best[sk], "holder": best["team_name"]}
res["top5"] = [{"rank": i + 1, "team": e["team_name"], "model_name": e.get("model_name"), "overall": e["score_avg"],
                "submitted": e.get("submission_date"), "raw_mse": e["expr_mse_unbiased_capped_norm"], "raw_pds": e["pds_cosine"],
                "score_mse": e.get("score_mse")} for i, e in enumerate(ents[:5])]
res["leader_gap"] = ents[0]["score_avg"] - ours["score_avg"]
json.dump(res, open(out, "w"), indent=1)
print(json.dumps(res, indent=1))
