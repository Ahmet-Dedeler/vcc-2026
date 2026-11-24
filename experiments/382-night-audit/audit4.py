"""382 item 4: board deltas from score_avg + score_avg == mean(six scaled members)."""
import json

D = "data/board-entries"
FILES = {"319": "xuDelEik7zGWMXdNu0R4.json", "319p": "SvzrTOQJAAaNjawzSyrC.json",
         "319pA": "sEa0LAsWFnXOoXSPJZ59.json", "254": "PVHrJHweHD5gvOf43ftL.json",
         "cd4svd5": "Xv0ayvuoHpsKYWv2iu42.json", "319d3": "5LOjb2VLiIBh0JG6aAgg.json"}
SIX = ["score_pds", "score_mse", "score_fid", "score_reach", "score_jac", "score_nmae"]

out = {"entries": {}, "deltas": {}}
for name, fn in FILES.items():
    d = json.load(open(f"{D}/{fn}"))
    e = d["entry"]
    mean6 = sum(e[k] for k in SIX) / 6
    diff = abs(mean6 - e["score_avg"])
    out["entries"][name] = {"file": fn, "model": e["model_name"],
                            "note": d.get("note"),
                            "score_avg": e["score_avg"],
                            "mean_of_six": mean6,
                            "absdiff": diff, "match_tol1e12": bool(diff <= 1e-12)}

checks = [("319p", "319", "+0.0011593"), ("319pA", "319p", "-0.0012794"),
          ("254", "cd4svd5", "+0.0013325"), ("319d3", "319", "-0.000269")]
for a, b, headline in checks:
    d = out["entries"][a]["score_avg"] - out["entries"][b]["score_avg"]
    out["deltas"][f"{a}_vs_{b}"] = {
        "delta_recomputed": d, "assignment_headline": headline,
        "absdiff_vs_headline": abs(d - float(headline))}

print(json.dumps(out, indent=1))
json.dump(out, open("experiments/382-night-audit/audit4.json", "w"), indent=1)
