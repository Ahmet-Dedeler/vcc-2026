#!/usr/bin/env python3
"""604: the A/B/C board artifacts for 639's identity-augmented CLEAN4 prior (VCC_SAMELINE format; the 458d path:
VCC_SAMELINE=<npz>:<sha256>:0, programme on, w 0, no rows; the hook can read the same file's pdown_<c>).
It runs the SAME function the day-of path runs (def_identity_prior.identity_prior) on the challenge bundle, with the
pinned Flex references (work/flex_refs604.npz): each context vs the other two + H1 (the guard drops its own pinned twin).
KA: lam 0 == CLEAN4 exactly (asserted inside identity_prior).
usage: uv run python experiments/604-identity-prior/build604.py --lam 0 --lam 0.5 --lam 1   (writes work/ident604v3_*)"""
import argparse, csv, json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import def_identity_prior as D
ap = argparse.ArgumentParser(); ap.add_argument("--lam", type=float, action="append", required=True)
ap.add_argument("--flip", default=None, help="625: 8-char 0/1 octile mask (def_identity_prior.flip_vector); all zeros == none")
ap.add_argument("--out-dir", default=None, help="default <604>/work (tests write elsewhere)")
a = ap.parse_args()
OUT = Path(a.out_dir) if a.out_dir else HERE / "work"
FLIP = a.flip if (a.flip and "1" in a.flip) else None
DATA = ROOT / "data/challenge"
with open(DATA / "gene_names.csv") as fh:
    genes = [r["gene_name"] for r in csv.DictReader(fh)]
with open(DATA / "pert_counts.csv") as fh:
    perts = [r["target_gene"] for r in csv.DictReader(fh)]
refs_meta = json.loads((HERE / "work/flex_refs604.npz.json").read_text())
rgenes, rcpm, rcalls = D.load_refs(HERE / "work/flex_refs604.npz", refs_meta["sha256"])
assert rgenes == genes
# the bundle's own contexts: the pinned A/B/C profiles ARE the challenge controls (same computation, same files)
ctx = {c: rcpm[c] for c in "ABC"}
calls = {c: rcalls[c] for c in "ABC"}
P_c4, ps_c4 = D.load_clean4(genes)
for lam in a.lam:
    res = D.identity_prior(genes, ctx, calls, rcpm, rcalls, lam, P_c4, ps_c4, flip_octiles=a.flip)
    art = {"genes": np.asarray(genes), "perts": np.asarray(perts), "contexts": np.asarray(sorted(res))}
    rep = {"lam": lam, "refs_sha256": refs_meta["sha256"], "clean4_sha256": D.CLEAN4_SHA256, "contexts": {}}
    for c, r in res.items():
        art[f"programme_{c}"] = r["psign"]; art[f"programme_ok_{c}"] = np.asarray(True); art[f"w_{c}"] = np.asarray(0.0)
        art[f"pdown_{c}"] = r["pdown"]; art[f"targets_{c}"] = np.zeros(0, dtype=np.int64)
        art[f"rows_{c}"] = np.zeros((0, len(genes)), dtype=np.float32); art[f"signs_{c}"] = np.zeros((0, len(genes)), dtype=np.int8)
        # the context's CPM block proxy: its top 2,200 control-CPM genes (the emitter's block ranks by emit_ref_lin)
        blk = np.argsort(-ctx[c], kind="stable")[:2200]
        tau = np.isfinite(r["pdown"]) & (np.abs(r["pdown"] - 0.5) >= 0.05)
        rep["contexts"][c] = {**r["prov"], "down_frac_block": float((r["psign"][blk] < 0).mean()),
                              "down_frac_block_clean4": float((ps_c4[blk] < 0).mean()),
                              "flips_vs_clean4_block": int((r["psign"][blk] != ps_c4[blk]).sum()),
                              "hook_coverage_tau0.05": int(tau.sum()),
                              "hook_coverage_tau0.05_clean4": int((np.isfinite(P_c4) & (np.abs(P_c4 - 0.5) >= 0.05)).sum())}
    art["provenance"] = np.asarray(json.dumps({k: v for k, v in rep.items() if k != "contexts"} | {"contexts": {c: r["prov"] for c, r in res.items()}}))
    # v3 (639d + d4cd5f 13:30Z): no identity-only extension; side = sign(score), confidence rank-mapped by |score|.
    # v1 (ident604_flex_*) and v2 (ident604v2_flex_*) stay on disk unchanged: 612's KGID prereg and 622k pin v1 by sha.
    dst = OUT / (f"ident604v3_flex_lam{lam:g}_flip{FLIP}_ABC.npz" if FLIP else f"ident604v3_flex_lam{lam:g}_ABC.npz")
    np.savez_compressed(dst, **art)
    rep["path"] = str(dst.relative_to(ROOT)) if dst.is_relative_to(ROOT) else str(dst); rep["sha256"] = D.file_sha256(dst)
    Path(str(dst) + ".json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1), flush=True)
