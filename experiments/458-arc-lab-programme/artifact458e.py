"""458e: the dose-curve partner of 458d. Every context's block prior = sign of 0.5*z(Arc-full P_down - 0.5) + 0.5*z(CLEAN4
P_down - 0.5) (Arc alone where CLEAN4 has no value), z-scaled per context on its own CPM block (2,200 top control-CPM
non-panel genes, as the builder's block), artifact458's format (programme_<c>, pdown_<c> = the blend mapped back to a
pseudo-P_down 0.5 - z/4 for provenance only, no rows, w 0). With VCC_SAMELINE=<npz>:<sha>:0 the builder uses it as each
context's block-prior default: 458d's path. usage: uv run python experiments/458-arc-lab-programme/artifact458e.py"""
import csv, hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "scripts"))
import def_sameline as S  # noqa: E402
official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
perts = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))]; panel = set(perts)
d, bp, bg, bsha = S.load_bank("h1_2025_full"); _, arc, nrows, _ = S.programme(d, bp, bg, official, panel)
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
cl = np.asarray(c4["P_down"], float)
assert np.array_equal(np.asarray(c4["genes"]).astype(str), np.asarray(official)), "CLEAN4 gene axis != official"
out = {"genes": np.asarray(official), "perts": np.asarray(perts), "contexts": np.asarray(["A", "B", "C"])}
flips = {}
for c in "ABC":
    b = np.load(ROOT / f"experiments/438-october-stack-repricing/work/abc_frames/bctrl-{c}.npy"); e = np.expm1(b)
    order = np.argsort(-(e / e.sum())); blk = np.asarray([i for i in order if official[i] not in panel][:2200])
    za = (arc - 0.5) / np.nanstd(arc[blk] - 0.5); zc = (cl - 0.5) / np.nanstd(cl[blk] - 0.5)
    z = np.where(np.isfinite(zc) & np.isfinite(za), 0.5 * za + 0.5 * zc, np.where(np.isfinite(za), za, zc))
    psign = np.where(np.isfinite(z) & (z <= 0), 1.0, -1.0)                 # 396's rule: down where P_down > 0.5 / z > 0
    out[f"programme_{c}"] = psign; out[f"programme_ok_{c}"] = np.asarray(True); out[f"w_{c}"] = np.asarray(0.0)
    out[f"pdown_{c}"] = np.where(np.isfinite(z), 0.5 + z / 4, np.nan)
    out[f"targets_{c}"] = np.zeros(0, dtype=np.int64); out[f"rows_{c}"] = np.zeros((0, len(official)), dtype=np.float32)
    out[f"signs_{c}"] = np.zeros((0, len(official)), dtype=np.int8)
    a_sign = np.where(np.isfinite(arc) & (arc <= 0.5), 1.0, -1.0); c_sign = np.asarray(c4["psign"], float)
    flips[c] = dict(vs_458d=float(np.mean(psign[blk] != a_sign[blk])), vs_clean4=float(np.mean(psign[blk] != c_sign[blk])))
prov = {"experiment": "458e", "banks": ["h1_2025_full", "clean4 (391)"], "bank_delta_sha256": bsha, "weight_arc": 0.5,
        "programme_rows_arc": nrows, "block_flips": flips, "note": "0.5*z(Arc-full) + 0.5*z(CLEAN4) block prior; no rows"}
out["provenance"] = np.asarray(json.dumps(prov))
dst = Path(__file__).resolve().parent / "work/arch1e_blend_abc.npz"
np.savez_compressed(dst, **out)
sha = hashlib.sha256(dst.read_bytes()).hexdigest()
Path(str(dst) + ".json").write_text(json.dumps({"sha256": sha, **prov}, indent=1))
print(f"wrote {dst} sha256 {sha}"); print(json.dumps(flips, indent=1))
