"""481: the October transfer check. Reads a finished `def_day.sh` work dir (plus its bundle) and prints every quantity this
audit found to change meaning between A/B/C and a new panel, next to its A/B/C value, with a flag. Read-only; < 2 GB.

    uv run python experiments/481-october-audit/transfer_check481.py <def_day work dir> --bundle <bundle> [--json out]

Known-answer rows: on 453's A/B/C-relabelled run every value must equal the A/B/C reference (strength <= 1 share 0.277,
q33 cut 2, NTC constructs 46 x 400, CD4 283/287/284 union 291, key floor share 0.279, ~40 % of blocks at the ballast
cap); on 465's efr run the flags fire where RESULTS.md says they do.

Rows (finding in RESULTS.md):
  F1  CD4 coverage per condition and union vs 047's gates 275 / 290
  F2  NTC construct sizes per context vs the manifest's cells_per_pert (the re-level's claim assumes they are equal)
  gate design strength per context from the stack dump: share <= 1 (q33 cuts at 1 once it exceeds 0.33), q33 threshold
  F4  the programme's delivered dose (median prior_lambda), when a programme is on
  F6  the panel median the carriers' scale is anchored to (whitening before_operator_median_norm; A/B/C 4.1637)
  F7  programme-source rows left after the panel is excluded
  --  key floor/cap/<=0 shares, ballast blocks at KMAX, H1 carrying targets, re-level factors and rho
"""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

import h5py
import numpy as np

REF = {"cd4": {"Rest": 283, "Stim8hr": 287, "Stim48hr": 284, "union": 291}, "strength_le1": 0.277, "q33": 2.0,
       "key_floor": 0.279, "key_cap": 0.143, "key_le0": 0.156, "ballast_at_cap": 0.40, "h1_carrying": 77,
       "median_norm": 4.1637, "lambda_ungated": 0.441, "programme_rows": 275, "factors": "1.086 / 1.138 / 1.166"}


def str_col(h, path):
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            return np.asarray(obj["categories"][:]).astype(str)[np.asarray(obj["codes"][:])]
        return np.asarray(obj["values"][:]).astype(str)
    return np.asarray(obj[:]).astype(str)


def control_gate(path: Path, n_genes: int) -> np.ndarray:
    """The render gate (pooled control CPM > 10 and b_ctrl > 0), streamed."""
    tot = np.zeros(n_genes)
    with h5py.File(path, "r") as h:
        ind, dat = h["X"]["indices"], h["X"]["data"]
        step = 20_000_000
        for s in range(0, ind.shape[0], step):
            tot += np.bincount(np.asarray(ind[s:s + step], np.int64), weights=np.asarray(dat[s:s + step], np.float64),
                               minlength=n_genes)
    return (1e6 * tot / tot.sum() > 10.0) & (np.log1p(5e4 * tot / tot.sum()) > 0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("work", type=Path)
    ap.add_argument("--bundle", type=Path, required=True)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args()
    W, B = a.work, a.bundle
    man = json.loads((B / "manifest.json").read_text())
    ctxs, cpp = list(man["contexts"]), int(man["cells_per_pert"])
    rows, out = [], {}

    def row(fid, name, got, ref, flag):
        rows.append((fid, name, got, ref, "FLAG" if flag else "ok"))
        out[name] = {"got": got, "ref": ref, "flag": bool(flag)}

    # F1 CD4 coverage
    log = W / "logs" / "cd4-extract.log"
    m = re.search(r"panel coverage (\{[^}]*\})", log.read_text()) if log.exists() else None
    cen = json.loads((W / "census.json").read_text()) if (W / "census.json").exists() else {}
    if m:
        cov = json.loads(m.group(1).replace("'", '"'))
        uni = None
        for r in cen.get("coverage_gate", {}).get("rows", []):
            if "CD4 union" in r["row"]:
                uni = r["got"]
        row("F1", "cd4_coverage", f"{cov} union {uni}", f"{REF['cd4']}", min(cov.values()) < 275 or (uni or 0) < 290)
    else:
        row("F1", "cd4_coverage", "no cd4-extract log", REF["cd4"], True)
    # F2 NTC construct sizes
    sizes = {}
    for c in ctxs:
        with h5py.File(B / f"context_{c}.h5ad", "r") as h:
            ids = str_col(h, "obs/ntc_id") if "ntc_id" in h["obs"] else np.array([])
        u, n = np.unique(ids, return_counts=True)
        sizes[c] = (int(u.size), int(n.min()) if n.size else 0, int(n.max()) if n.size else 0)
    row("F2", "ntc_constructs(n,min,max)", sizes, f"46 x {cpp}", any(v[1] != cpp or v[2] != cpp for v in sizes.values()))
    # design strength from the stack dump
    dump = W / "dump_h1x2_cd4x1.npz"
    if dump.exists():
        z = np.load(dump, allow_pickle=False)
        genes = z["genes"].astype(str)
        car = np.asarray(z["carriers"], dtype=np.float64)
        st = {}
        for c in ctxs:
            s = (np.abs(np.where(control_gate(B / f"context_{c}.h5ad", genes.size)[None, :], car, 0.0)) > 0.2).sum(1)
            st[c] = (round(float((s <= 1).mean()), 3), float(np.quantile(s, 0.33)))
        row("gate", "strength_le1_share,q33", st, f"{REF['strength_le1']}, q33 {REF['q33']}",
            any(v[0] > 0.33 or v[1] < 2 for v in st.values()))
    finals = sorted(glob.glob(str(W / "*_final.h5ad.meta.json")))
    meta = json.loads(Path(finals[0]).read_text()) if finals else None
    if meta:
        cfg, eff, rc = meta["config"], meta["effects"], meta["output"]["rc_callset"]
        prog = any(v.get("programme_default") for v in (cfg.get("sameline") or {}).get("contexts", {}).values())
        lam = rc.get("prior_lambda", {}).get("median")
        row("F4", "programme_dose_median_lambda" + ("" if prog else " (no programme)"), lam,
            f"{REF['lambda_ungated']} ungated A/B/C", prog and lam is not None and lam < 0.35)
        mn = eff.get("whitening", {}).get("before_operator_median_norm")
        row("F6", "carrier_anchor_median_norm", round(mn, 4) if mn else None, REF["median_norm"],
            mn is not None and abs(mn / REF["median_norm"] - 1) > 0.05)
        row("--", "h1_carrying_targets", eff.get("h1_union", {}).get("carrying_targets"), REF["h1_carrying"], False)
        row("--", "bulk_ctx_scale", cfg.get("bulk_ctx_scale"), REF["factors"],
            any(not 0.8 <= v <= 1.3 for v in (cfg.get("bulk_ctx_scale") or {}).values()))
        if "de_prior_gate" in cfg:
            row("gate", "builder_prior_gate", json.dumps(cfg["de_prior_gate"])[:160], "le1 (472) or q33 at A/B/C: cut 2", False)
    # F7 programme rows
    slj = W / "sameline.npz.json"
    if slj.exists():
        prov = json.loads(slj.read_text())
        pr = {c: v.get("programme_source", {}).get("programme_rows") for c, v in prov.get("contexts", {}).items()}
        pr = {c: v for c, v in pr.items() if v is not None}
        if pr:
            row("F7", "programme_source_rows", pr, REF["programme_rows"], any(v < 150 for v in pr.values()))
    # keys
    for kj in sorted(glob.glob(str(W / "pert_scale_*def*.json"))):
        k = json.loads(Path(kj).read_text())
        row("--", f"key {Path(kj).stem}: floor/cap/<=0", (round(k["frac_at_floor"], 3), round(k["frac_at_cap"], 3),
            round(k["frac_key_le0"], 3)), (REF["key_floor"], REF["key_cap"], REF["key_le0"]), k["frac_at_floor"] > 0.40)
    # ballast blocks at KMAX (ballast cells: nnz > 8,000; 481 ballastK)
    fin = sorted(glob.glob(str(W / "*_final.h5ad")))
    if fin:
        with h5py.File(fin[0], "r") as h:
            nnz = np.diff(h["X"]["indptr"][:])
        K = (nnz > 8000).reshape(-1, cpp).sum(1)
        row("--", "ballast_blocks_at_KMAX16", round(float((K >= 16).mean()), 3), REF["ballast_at_cap"], (K >= 16).mean() > 0.6)
    print(f"{'':5s}{'quantity':38s} {'this run':60s} {'A/B/C':34s} flag")
    for fid, name, got, ref, flag in rows:
        print(f"{fid:5s}{name[:38]:38s} {str(got)[:60]:60s} {str(ref)[:34]:34s} {flag}")
    if a.json:
        a.json.write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
