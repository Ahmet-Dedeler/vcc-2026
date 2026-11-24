"""223 Phase 2a(ii) — build panel-covering CMap donor tables from the two LINCS-2020
level-5 files already on disk.

176 measured CMap only on **h1-test's 85-target DE pool**, as a SIGN predictor, and stopped
at its own gate 3 ("only if 3 is positive: pds/cos in delivered space").  Gate 3 came out
marginal (+0.005 on 867 landmark genes, nothing on 7,101 inferred), so the delivered-space
measurement was never made -- and CMap is the largest partial-coverage source we hold
(union 127/300 = 42.3%, i.e. squarely inside the median-rule's void band).

Three tables, all on the official 18,533-gene axis:
  cmap_sh_cgs        9 core trt_sh.cgs lines, mean z where a line covers a target
  cmap_sh_cgs_lm     the same, with every INFERRED coordinate zeroed (176 gate 2 measured
                     the inferred axis at chance, 0.49-0.53, against landmark 0.50-0.55)
  cmap_xpr           trt_xpr CRISPR-KO, mean z per target across every line that has it

Written only into this experiment's work dir.  Sign convention is 176's, asserted below by
the own-target check: a knockdown's own gene must read NEGATIVE z on median.

    uv run python build_cmap.py
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import h5py
import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "work/tables"
OUT.mkdir(parents=True, exist_ok=True)

OFFICIAL = [line.strip() for line in
            (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:]
            if line.strip()]
assert len(OFFICIAL) == 18533
GPOS = {g: i for i, g in enumerate(OFFICIAL)}
PANEL = [line.split(",")[0]
         for line in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
         if line.strip()]
PANEL_SET = set(PANEL)

# 176's INVENTORY.md: the nine trt_sh lines with >=100 panel targets.
SH_LINES = ("VCAP", "PC3", "A375", "HT29", "HA1E", "A549", "MCF7", "HCC515", "HEPG2")
T0 = time.time()


def gene_axis(rid: np.ndarray):
    gi = pl.read_csv(ROOT / "data/external/cmap/geneinfo_beta.txt", separator="\t",
                     infer_schema_length=0)
    id2sym = dict(zip(gi["gene_id"].to_list(), gi["gene_symbol"].to_list()))
    id2fs = dict(zip(gi["gene_id"].to_list(), gi["feature_space"].to_list()))
    sym = np.asarray([id2sym.get(r, "") for r in rid])
    fs = np.asarray([id2fs.get(r, "") for r in rid])
    on = np.asarray([s in GPOS for s in sym])
    col = np.asarray([GPOS[s] if s in GPOS else -1 for s in sym])
    return on, col, (fs == "landmark")


def emit(name: str, tot: np.ndarray, cnt: np.ndarray, perts: list[str],
         extra: dict | None = None) -> dict:
    delta = np.zeros_like(tot, dtype=np.float32)
    live = cnt > 0
    delta[live] = (tot[live] / cnt[live]).astype(np.float32)
    # never let a panel-target coordinate leak into a residual
    pcols = np.asarray([GPOS[g] for g in PANEL if g in GPOS])
    delta[:, pcols] = 0.0
    keep = (delta != 0).any(axis=1)
    delta = np.ascontiguousarray(delta[keep])
    perts = [p for p, k in zip(perts, keep) if k]
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "delta.npy", delta)
    (d / "perts.txt").write_text("\n".join(perts) + "\n")
    (d / "genes.txt").write_text("\n".join(OFFICIAL) + "\n")
    np.save(d / "b_ctrl.npy", np.zeros(len(OFFICIAL), dtype=np.float32))
    rn = (delta != 0).sum(1)
    meta = dict(name=name, shape=list(delta.shape),
                panel_targets=len(set(perts) & PANEL_SET),
                median_nonzero_genes_per_row=int(np.median(rn)),
                rows_all_zero=int((rn == 0).sum()), **(extra or {}))
    (d / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[{name}] {delta.shape} panel {meta['panel_targets']}/300 "
          f"median nz/row {meta['median_nonzero_genes_per_row']} [{time.time()-T0:.0f}s]",
          flush=True)
    return meta


def own_target_check(tot, cnt, perts, name) -> dict:
    """176's sign-convention gate, re-run here: a knockdown's own gene must read z < 0."""
    vals = []
    for i, p in enumerate(perts):
        j = GPOS.get(p)
        if j is None or cnt[i, j] == 0:
            continue
        vals.append(tot[i, j] / cnt[i, j])
    v = np.asarray(vals)
    r = dict(source=name, n=int(v.size), median_own_z=float(np.median(v)) if v.size else None,
             frac_negative=float((v < 0).mean()) if v.size else None)
    print(f"    own-target check {name}: n={r['n']} median z {r['median_own_z']} "
          f"frac neg {r['frac_negative']}", flush=True)
    assert v.size == 0 or r["median_own_z"] < 0, \
        f"{name}: own-target z is not negative -- sign convention is wrong"
    return r


def build_sh() -> tuple[dict, dict, dict]:
    p = ROOT / "data/external/cmap/level5_beta_trt_sh.gctx"
    with h5py.File(p, "r") as h:
        cid = np.asarray([x.decode() for x in h["0/META/COL/id"][:]])
        rid = np.asarray([x.decode() for x in h["0/META/ROW/id"][:]])
        on, col, is_lm = gene_axis(rid)
        is_cgs = np.char.startswith(cid, "CGS")
        idx_all = np.flatnonzero(is_cgs)
        parts = [c.split(":") for c in cid[is_cgs]]
        line = np.asarray([q[0].split("_")[1] for q in parts])
        tgt = np.asarray([q[1] for q in parts])
        keep_line = np.isin(line, SH_LINES)
        perts = sorted(set(tgt[keep_line]))
        tpos = {t: i for i, t in enumerate(perts)}
        print(f"trt_sh CGS: {is_cgs.sum()} sigs, {keep_line.sum()} on the 9 core lines, "
              f"{len(perts)} distinct targets [{time.time()-T0:.0f}s]", flush=True)
        tot = np.zeros((len(perts), len(OFFICIAL)), dtype=np.float64)
        cnt = np.zeros((len(perts), len(OFFICIAL)), dtype=np.uint8)
        oncols, ocols = np.flatnonzero(on), col[on]
        for L in SH_LINES:
            sel = np.flatnonzero(keep_line & (line == L))
            rows = idx_all[sel]
            order = np.argsort(rows)
            M = np.empty((rows.size, rid.size), dtype=np.float32)
            for k in range(0, rows.size, 2000):
                rr = rows[order[k:k + 2000]]
                M[order[k:k + 2000]] = h["0/DATA/0/matrix"][rr, :]
            for i, t in enumerate(tgt[sel]):
                r = tpos[t]
                tot[r, ocols] += M[i, oncols]
                cnt[r, ocols] += 1
            del M
            print(f"  {L:<8} {rows.size:5d} CGS sigs [{time.time()-T0:.0f}s]", flush=True)
    chk = own_target_check(tot, cnt, perts, "cmap_sh_cgs")
    lm_official = np.zeros(len(OFFICIAL), dtype=bool)
    lm_official[col[on & is_lm]] = True
    a = emit("cmap_sh_cgs", tot, cnt, perts,
             dict(lines=list(SH_LINES), own_target=chk,
                  landmark_genes_on_axis=int(lm_official.sum())))
    cnt_lm = cnt.copy()
    cnt_lm[:, ~lm_official] = 0
    tot_lm = tot.copy()
    tot_lm[:, ~lm_official] = 0.0
    b = emit("cmap_sh_cgs_lm", tot_lm, cnt_lm, perts,
             dict(lines=list(SH_LINES), own_target=chk,
                  note="inferred coordinates zeroed; 176 gate 2 measured them at chance"))
    return a, b, chk


def build_xpr() -> dict:
    sig = pl.read_parquet(ROOT / "data/external/cmap/genetic_sigs.parquet").filter(
        pl.col("pert_type") == "trt_xpr")
    want = dict(zip(sig["sig_id"].to_list(),
                    zip(sig["cell_iname"].to_list(), sig["cmap_name"].to_list())))
    p = ROOT / "data/external/cmap/level5_beta_trt_xpr.gctx"
    with h5py.File(p, "r") as h:
        cid = np.asarray([x.decode() for x in h["0/META/COL/id"][:]])
        rid = np.asarray([x.decode() for x in h["0/META/ROW/id"][:]])
        on, col, _ = gene_axis(rid)
        rows = np.asarray([i for i, c in enumerate(cid) if c in want])
        tg = np.asarray([want[cid[r]][1] for r in rows])
        perts = sorted(set(tg))
        tpos = {t: i for i, t in enumerate(perts)}
        print(f"trt_xpr: {rows.size} sigs, {len(perts)} distinct targets "
              f"[{time.time()-T0:.0f}s]", flush=True)
        tot = np.zeros((len(perts), len(OFFICIAL)), dtype=np.float64)
        cnt = np.zeros((len(perts), len(OFFICIAL)), dtype=np.uint16)
        oncols, ocols = np.flatnonzero(on), col[on]
        for k in range(0, rows.size, 1000):
            blk = rows[k:k + 1000]
            M = h["0/DATA/0/matrix"][blk, :]
            for i, t in enumerate(tg[k:k + 1000]):
                r = tpos[t]
                tot[r, ocols] += M[i, oncols]
                cnt[r, ocols] += 1
            del M
            if (k // 1000) % 20 == 0:
                print(f"  xpr {k}/{rows.size} [{time.time()-T0:.0f}s]", flush=True)
    chk = own_target_check(tot, cnt, perts, "cmap_xpr")
    return emit("cmap_xpr", tot, cnt.astype(np.float64), perts, dict(own_target=chk))


def main() -> None:
    out = {}
    a, b, _ = build_sh()
    out["cmap_sh_cgs"], out["cmap_sh_cgs_lm"] = a, b
    out["cmap_xpr"] = build_xpr()
    (HERE / "work/build_cmap.json").write_text(json.dumps(out, indent=2))
    print(f"-> {HERE / 'work/build_cmap.json'}")


if __name__ == "__main__":
    main()
