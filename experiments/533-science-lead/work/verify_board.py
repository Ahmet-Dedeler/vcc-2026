"""533 §2: recompute the board arithmetic behind every board number the manuscript will state.

Reads only data/board-entries/<id>.json (the captured records). For each entry it checks that the
published overall equals the unweighted mean of the six scaled members, and for each (arm, base) pair
it recomputes the raw and scaled deltas and the overall delta, then compares them with the numbers
written in the experiments' BOARD-RESULT files (transcribed below with file:line in the comments).

Also prints the implied conversion Δscaled / Δraw per member (the context-weighted 1/(r − b)); it
must be stable across arms, which is an internal-consistency check on the scaled members (the
per-context anchors are not in the captured records, so s = (u − b)/(r − b) cannot be recomputed
from the context-averaged raw values directly).

Run: uv run python experiments/533-science-lead/work/verify_board.py
"""
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[3]
BE = ROOT / "data" / "board-entries"

ENTRIES = {
    "395": "RzcGA26xZEigSuXgUp5t",
    "410": "EETok0iqIbzQFg4CTYMQ",
    "410b": "4WgVi6FIE4iVw9PmRUP3",
    "417e": "wPP8U70kj4mOIvtA83CO",
    "433a": "FCbBcvGxF1HOFPxPISUC",
    "433h": "wMoDAi1mca08FA2Bfqkd",
    "433k": "T76JWD3tOajk5IYvzwUW",
    "446b": "UDgjOKL9gI1udAexq6ue",
    "446bja": "AsHteaDAXCWQo57fLYti",
    "446bj2": "YvVru3DAF89Jiyr34olR",
    "446br2": "wE7usCaR6TAlPp16TR0V",
    "453dr": "q3VuTImh4FRGKXH2waCp",
    "446bj2z": "be90smZUDcIV7j5MgQ03",
    "457dr": "1FfirZsGsAZvaXYTd39J",
    "446bd2": "R5uE5Rve3R9nWXGbeSLM",
}

RAW = {
    "pds": "pds_cosine",
    "mse": "expr_mse_unbiased_capped_norm",
    "nmae": "de_wilcoxon_lfc_nmae",
    "fid": "de_wilcoxon_direction_fidelity_yield_raw",
    "reach": "de_wilcoxon_direction_reach_raw",
    "jac": "de_wilcoxon_sig_jaccard",
}
MEMBERS = list(RAW)

# (arm, base, stated overall delta, {member: stated raw delta}, source)
STATED = [
    ("446bj2", "446b", +0.0022338,
     {"pds": +0.000208, "mse": -0.014473, "nmae": +0.000060, "fid": -0.000022, "reach": -0.001963, "jac": +0.000016},
     "446-q175-level/BOARD-RESULT-446bj2.md"),
    ("446br2", "446b", +0.0003578,
     {"pds": +0.000007, "mse": +0.000070, "nmae": -0.001413, "fid": 0.0, "reach": -0.000039, "jac": 0.0},
     "446-q175-level/BOARD-RESULT-446br2.md"),
    ("453dr", "446b", -0.0036816,
     {"pds": -0.001345, "mse": +0.015033, "nmae": -0.000188, "fid": -0.000036, "reach": -0.003436, "jac": -0.000004},
     "453-abc-dress-rehearsal/BOARD-RESULT-453dr.md"),
    ("446bj2z", "446bj2", +0.0005493,
     {"pds": 0.0, "mse": 0.0, "nmae": -0.002179, "fid": 0.0, "reach": 0.0, "jac": 0.0},
     "456-zs-riders/BOARD-RESULT-446bj2z.md"),
    ("457dr", "453dr", +0.0031484,
     {"pds": +0.001713, "mse": -0.014313, "nmae": +0.000200, "fid": -0.000034, "reach": +0.000711, "jac": +0.000004},
     "457-cd4-half-def/BOARD-RESULT-457dr.md"),
    ("446bd2", "446b", +0.0004200,
     {"pds": +0.002148, "mse": -0.000471, "nmae": +0.000121, "fid": -0.000020, "reach": -0.002314, "jac": +0.000000},
     "445-dld1-bulk-union/BOARD-RESULT-446bd2.md"),
]

# Derived statements to check (arm, base, stated overall delta, source)
DERIVED = [
    ("457dr", "446b", -0.0005332, "457-cd4-half-def/BOARD-RESULT-457dr.md ('What it adds')"),
    ("446bj2z", "446b", None, "chain 446b -> 446bj2 -> 446bj2z"),
    ("446bja", "446b", -0.0010661, "450/BOARD-RESULT-446bja.md (already in the paper, C-63)"),
]


def load(tag):
    d = json.loads((BE / f"{ENTRIES[tag]}.json").read_text())["entry"]
    return d


def main():
    ent = {t: load(t) for t in ENTRIES}
    print("== overall = mean of six scaled members ==")
    worst = 0.0
    for t, d in ent.items():
        m = sum(d[f"score_{k}"] for k in MEMBERS) / 6
        diff = m - d["score_avg"]
        worst = max(worst, abs(diff))
        print(f"{t:8s} {d['entry_id']} score_avg {d['score_avg']:.10f} mean6 {m:.10f} diff {diff:+.1e}")
    print(f"max |diff| {worst:.1e}")

    print("\n== stated deltas vs recomputed ==")
    conv = {k: [] for k in MEMBERS}
    for arm, base, d_over, d_raw, src in STATED:
        a, b = ent[arm], ent[base]
        over = a["score_avg"] - b["score_avg"]
        flag = "OK" if abs(over - d_over) < 6e-8 else "MISMATCH"
        print(f"\n{arm} vs {base} ({src}) overall {over:+.7f} stated {d_over:+.7f} {flag}")
        for k in MEMBERS:
            dr = a[RAW[k]] - b[RAW[k]]
            ds = a[f"score_{k}"] - b[f"score_{k}"]
            st = d_raw[k]
            ok = abs(dr - st) < 6e-7
            bit = "bit" if a[RAW[k]] == b[RAW[k]] else ""
            ratio = ds / dr if dr != 0 else float("nan")
            if abs(dr) > 5e-4:
                conv[k].append((arm, ratio))
            print(f"  {k:5s} raw {dr:+.6f} (stated {st:+.6f}) {'ok' if ok else 'MISMATCH'} {bit:3s}"
                  f" scaled {ds:+.6f} overall {ds / 6:+.7f} ds/dr {ratio:+.3f}")

    print("\n== derived statements ==")
    for arm, base, st, src in DERIVED:
        over = ent[arm]["score_avg"] - ent[base]["score_avg"]
        print(f"{arm} vs {base}: {over:+.7f} stated {st} ({src})")
        bulk = sum(ent[arm][f"score_{k}"] - ent[base][f"score_{k}"] for k in ("pds", "mse")) / 6
        de = sum(ent[arm][f"score_{k}"] - ent[base][f"score_{k}"] for k in ("nmae", "fid", "reach", "jac")) / 6
        print(f"   bulk (pds+mse)/6 {bulk:+.7f}  DE (4 members)/6 {de:+.7f}")

    print("\n== implied Δscaled/Δraw per member (context-weighted 1/(r-b); only |Δraw| > 5e-4) ==")
    for k, v in conv.items():
        print(f"  {k:5s} " + ", ".join(f"{a} {r:+.3f}" for a, r in v))

    print("\n== the October-path decomposition on A/B/C (446b -> 453dr -> 457dr) ==")
    for k in MEMBERS:
        r1 = ent["453dr"][RAW[k]] - ent["446b"][RAW[k]]
        r2 = ent["457dr"][RAW[k]] - ent["453dr"][RAW[k]]
        r3 = ent["457dr"][RAW[k]] - ent["446b"][RAW[k]]
        print(f"  {k:5s} 446b->453dr {r1:+.6f}  453dr->457dr {r2:+.6f}  446b->457dr {r3:+.6f}")

    print("\n== the chain since 446b (overall) ==")
    for t in ("446b", "446bja", "446bj2", "446br2", "453dr", "446bj2z", "457dr", "446bd2"):
        print(f"  {t:8s} {ent[t]['score_avg']:.7f}  {ent[t]['entry_id']}  submitted {ent[t]['submission_date']}")


if __name__ == "__main__":
    main()
