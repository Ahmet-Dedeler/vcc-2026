"""487 decision 5: read 446bj2 (446b + the LODO key at matched rho) against 446b from the coordinator's captured board
records (data/board-entries/*.json; read-only), apply 433's PREREGISTRATION-0926 bins, and translate the A/B/C read to
D/E/F with 438's levelled PLOLO ratios (LODO V1 (0.5, 3) vs the 410 key on the stack: abc +0.00274, real +0.00318,
efp +0.00306, efr +0.00169).
usage: uv run python experiments/487-october-decision-review/read487.py [entry_json]  -> work/read487.json
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASE = ROOT / "data/board-entries/UDgjOKL9gI1udAexq6ue.json"          # 446b, 0.2700428
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
       "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard", "nmae": "de_wilcoxon_lfc_nmae"}
SC = {k: f"score_{k}" for k in RAW}
PLOLO = {"abc": 0.00274, "real": 0.00318, "efp": 0.00306, "efr": 0.00169}


def find():
    if len(sys.argv) > 1:
        return Path(sys.argv[1]).resolve()
    cands = []
    for p in (ROOT / "data/board-entries").glob("*.json"):
        d = json.loads(p.read_text())
        txt = (d.get("note") or "") + " " + (d.get("package") or "")
        if "446bj2" in txt and "446bj2z" not in txt:
            cands.append((d.get("captured_at", ""), p))
    if not cands:
        raise SystemExit("no captured 446bj2 entry yet")
    return sorted(cands)[-1][1]


def main():
    b = json.loads(BASE.read_text())["entry"]
    p = find()
    rec = json.loads(p.read_text())
    e = rec["entry"]
    out = {"entry": e["entry_id"], "file": str(p.relative_to(ROOT)), "note": rec.get("note"), "package": rec.get("package"),
           "overall": e["score_avg"], "overall_446b": b["score_avg"], "d_overall": e["score_avg"] - b["score_avg"],
           "d_raw": {k: e[v] - b[v] for k, v in RAW.items()}, "d_scaled": {k: e[v] - b[v] for k, v in SC.items()}}
    d = out["d_overall"]
    bulk = (out["d_scaled"]["pds"] + out["d_scaled"]["mse"]) / 6
    out["d_overall_bulk_members"] = bulk
    key_ok = d >= 0.001 and out["d_raw"]["pds"] >= 0
    out["bin"] = ("key_ok (transfers)" if key_ok else "flat (keep 410's key)" if abs(d) < 0.001 else
                  "close the LODO key" if d <= -0.001 else "pds < 0 with overall >= +0.001: not key_ok")
    out["prediction_hit"] = {"pds +0.001..+0.004": 0.001 <= out["d_raw"]["pds"] <= 0.004,
                             "mse -0.004..-0.012 raw": -0.012 <= out["d_raw"]["mse"] <= -0.004,
                             "overall +0.001..+0.003": 0.001 <= d <= 0.003}
    out["def_projection"] = {k: d * v / PLOLO["abc"] for k, v in PLOLO.items() if k != "abc"}
    print(json.dumps(out, indent=1))
    (HERE / "work/read487.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
