"""458: a same-line-format artifact that gives EVERY context (A, B, C) Arc's own H1 programme (P_down over data/donors/h1_2025,
255 perturbations, Flex; the A/B/C panel's rows excluded, def_sameline's 396 rule) and no rows. With
VCC_SAMELINE=<npz>:<sha>:0 (programme on, w 0) the builder uses it as each context's block-prior default, exactly the path
451's KA-2 proved equal to 446bja for context A. No builder change.   usage: uv run python .../artifact458.py"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import def_sameline as S  # noqa: E402

S.BANK_DIRS["h1_2025"] = ROOT / "data/donors/h1_2025"
official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
perts = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))]
panel = set(perts)
delta, bperts, bgenes, bsha = S.load_bank("h1_2025")
psign, pdown, nrows, _ = S.programme(delta, bperts, bgenes, official, panel)
out = {"genes": np.asarray(official), "perts": np.asarray(perts), "contexts": np.asarray(["A", "B", "C"])}
for c in ("A", "B", "C"):
    out[f"programme_{c}"] = psign
    out[f"programme_ok_{c}"] = np.asarray(True)
    out[f"w_{c}"] = np.asarray(0.0)
    out[f"pdown_{c}"] = pdown
    out[f"targets_{c}"] = np.zeros(0, dtype=np.int64)
    out[f"rows_{c}"] = np.zeros((0, len(official)), dtype=np.float32)
    out[f"signs_{c}"] = np.zeros((0, len(official)), dtype=np.int8)
prov = {"experiment": "458", "bank": "h1_2025", "bank_delta_sha256": bsha, "programme_rows": nrows,
        "programme_coverage": float(np.isfinite(pdown).mean()), "programme_frac_down": float((psign < 0).mean()),
        "note": "Arc's own H1 programme for every validation context (lab-signature hypothesis, 450h/450l); no rows"}
out["provenance"] = np.asarray(json.dumps(prov))
dst = Path(__file__).resolve().parent / "work/arch1_programme_abc.npz"
np.savez_compressed(dst, **out)
sha = hashlib.sha256(dst.read_bytes()).hexdigest()
Path(str(dst) + ".json").write_text(json.dumps({"sha256": sha, **prov}, indent=1))
# how different is it from CLEAN4 on each context's block? (CLEAN4 psign: 391's pinned file)
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
print(f"wrote {dst} sha256 {sha}")
print(f"programme rows {nrows}; coverage {prov['programme_coverage']:.3f}; frac down {prov['programme_frac_down']:.3f}")
print("clean4 artifact keys:", c4.files[:8])
if "psign" in c4.files:
    cs = np.asarray(c4["psign"], dtype=float)
    if cs.shape[0] == len(official):
        print(f"sign differs from CLEAN4 on {np.mean(cs != psign):.3f} of all genes")
