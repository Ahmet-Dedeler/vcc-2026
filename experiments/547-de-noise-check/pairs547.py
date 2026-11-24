"""547: coordinator's independent check of 537's claim that the DE members' board noise between arms that keep the DE
design is far below AGENTS 11.1 rule 13's 0.0045 overall. Pairs = our captured entries (data/board-entries) whose package
metas differ ONLY in bulk keys; `sameline` counts as bulk only when neither side carries a programme (458-family arms put
Arc's programme, a DE block prior, under `sameline`). Prints each pair's raw DE deltas and the rms, raw and in overall units
(raw / (r - b) / 6 with each member's mid-range r - b from AGENTS §4).
usage: uv run python experiments/547-de-noise-check/pairs547.py"""
import glob
import itertools
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"nmae": "de_wilcoxon_lfc_nmae", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
       "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}
SPAN = {"nmae": 1.001 - 0.40, "fid": 0.81 - 0.513, "reach": 0.968 - 0.072, "jac": 0.40 - 0.029}   # |r - b|, mid-range
BULK = {"bulk_ctx_scale", "bulk_load", "bulk_load_sha256", "bulk_sharpen_q", "bulk_sharpen_pt", "sameline"}
rows = []
for f in sorted(glob.glob(str(ROOT / "data/board-entries/*.json"))):
    d = json.load(open(f)); e = d.get("entry", {}); pk = d.get("package") or ""
    mp = ROOT / pk.replace(".vcc", ".h5ad.meta.json")
    if e.get("score_avg") is None or not mp.exists():
        continue
    c = json.load(open(mp)).get("config", {})
    prog = "programme" in json.dumps(c.get("sameline"), sort_keys=True).lower() or "arch1" in os.path.basename(pk)
    rows.append((e["entry_id"], os.path.basename(pk)[:30], c, {m: e[RAW[m]] for m in RAW}, prog))
pairs = []
for a, b in itertools.combinations(rows, 2):
    if a[4] or b[4]:
        continue
    keys = set(a[2]) | set(b[2])
    diff = {k for k in keys if json.dumps(a[2].get(k), sort_keys=True) != json.dumps(b[2].get(k), sort_keys=True)}
    if diff and diff <= BULK:
        pairs.append((a[1], b[1], sorted(diff), {m: b[3][m] - a[3][m] for m in RAW}))
print(f"{len(rows)} entries with metas; {len(pairs)} bulk-only pairs (programme arms excluded)")
for n1, n2, diff, dd in sorted(pairs, key=lambda p: -abs(p[3]["reach"])):
    print(f"  {n1:30s} {n2:30s} reach {dd['reach']:+.5f} fid {dd['fid']:+.6f} jac {dd['jac']:+.6f} nmae {dd['nmae']:+.5f} | {','.join(diff)}")
out = {}
for m in RAW:
    v = [p[3][m] for p in pairs]
    rms = (sum(x * x for x in v) / len(v)) ** 0.5 if v else float("nan")
    out[m] = {"rms_raw": rms, "max_raw": max(abs(x) for x in v) if v else None, "rms_overall": rms / SPAN[m] / 6}
    print(f"{m:5s} rms raw {rms:.6f}  max {out[m]['max_raw']:.6f}  -> overall {out[m]['rms_overall']:.6f}")
tot = sum(x["rms_overall"] ** 2 for x in out.values()) ** 0.5
print(f"DE members together (rms, quadrature): {tot:.6f} overall  vs rule 13's 0.0045")
json.dump({"pairs": [{"a": p[0], "b": p[1], "diff": p[2], "d": p[3]} for p in pairs], "summary": out, "de_total_overall": tot},
          open(ROOT / "experiments/547-de-noise-check/pairs547.json", "w"), indent=1)
