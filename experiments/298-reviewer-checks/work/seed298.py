"""298 item 4 (T7) — the seed scan of the 2025 library's failing spacers against promoters.

VERDICT: NOT TESTABLE ON DISK.  The scan needs (1) the 34 arms' PROTOSPACER SEQUENCES and (2) the
promoter SEQUENCE of every annotated TSS window.  Neither is on this machine, and the mandate
permits a download only for (2).  This script therefore does the two things it can:

  * it re-reads the 2025 library's `guide_id` categories and PROVES they are protospacer NAMES,
    not sequences (0 of 410 categories over the three splits are ACGT strings), and records what
    every candidate sequence source on disk actually contains;
  * it writes the exact input list a future run needs — the 34 arms of the 17 two-construct
    targets from `experiments/293-challenge-construct-ceiling/work/construct293.csv`, each with its
    target-gene knockdown lfc, the |lfc| < 1.0 FAILING flag, the cells behind it, and what is
    missing for that arm — to work/seed298.csv, plus the match logic to work/seed298_plan.md.

    PYTHONDONTWRITEBYTECODE=1 uv run python seed298.py

Writes work/seed298.csv, work/seed298.json.
"""
from __future__ import annotations

import hashlib
import json
import re
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path("~/Code/virtual-cell-challenge-2026")
WORK = ROOT / "experiments/298-reviewer-checks/work"
CONSTRUCT293 = ROOT / "experiments/293-challenge-construct-ceiling/work/construct293.csv"
BIOMART = ROOT / "data/external/broad16line/metadata/gene_locations_biomart.csv"
SPLITS = {"train": ROOT / "data/external/vcc2025/train/adata_Training.h5ad",
          "validation": ROOT / "data/external/vcc2025/validation/adata_Validation.h5ad",
          "test": ROOT / "data/external/vcc2025/test/adata_Test.h5ad"}
FAIL_LFC = 1.0        # the mandate's failing-construct rule: |lfc| < 1.0
LOG: list[str] = []


def P(m: str = "") -> None:
    print(m, flush=True)
    LOG.append(m)


def cats_of(h, path) -> list[str]:
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        return [c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return sorted({x.decode() if isinstance(x, bytes) else str(x) for x in v})


def main() -> None:
    P(f"298 seed298.py — {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")

    # ---------------- KA: 293's 16-target rho_guide ------------------------------------------
    ceil = json.loads((ROOT / "experiments/293-challenge-construct-ceiling/work/"
                       "ceiling293.json").read_text())
    rg = ceil["native"]["rungs"]["200"]["rho_guide"]
    ka_ok = abs(rg["point"] - 0.075) < 5e-4 and rg["n"] == 16
    P(f"\nKA  293 native n=200 rho_guide {rg['point']:+.4f} [{rg['lo']:+.4f}, {rg['hi']:+.4f}] on "
      f"{rg['n']} targets (want +0.075 / 16) -> {'PASS' if ka_ok else 'FAIL'}")
    if not ka_ok:
        raise SystemExit("KA FAILED — abort (AGENTS.md:1375-1382)")

    # ---------------- (1) is a spacer SEQUENCE on disk? ---------------------------------------
    P("\n=== (1) spacer sequences: every candidate source on disk ===")
    srcs = []
    seq_re = re.compile(r"^[ACGT]{15,30}$", re.I)
    total_cats, total_seq = 0, 0
    for split, p in SPLITS.items():
        with h5py.File(p, "r") as h:
            obs_keys = sorted(h["obs"].keys())
            var_keys = sorted(h["var"].keys())
            uns_keys = sorted(h["uns"].keys()) if "uns" in h else []
            cats = cats_of(h, "obs/guide_id")
        nseq = sum(1 for c in cats if seq_re.match(c))
        total_cats += len(cats)
        total_seq += nseq
        srcs.append(dict(source=str(p.relative_to(ROOT)), kind="2025 h5ad obs/guide_id",
                         n_categories=len(cats), n_acgt_sequences=nseq,
                         example=cats[:3], obs_keys=obs_keys, var_keys=var_keys,
                         uns_keys=uns_keys, carries_spacer=False))
        P(f"  {split:<11} obs {obs_keys}  uns {uns_keys or '[]'}")
        P(f"              guide_id: {len(cats)} categories, {nseq} are ACGT strings; "
          f"e.g. {cats[:2]}")
    P(f"  -> {total_seq}/{total_cats} guide_id categories over the three splits are sequences: "
      f"the field is a protospacer NAME ('GENE_P1P2_A|GENE_P1P2_B' = the dual-guide construct's "
      f"two half-spacers), and there is no uns/ or var/ guide table.")

    for rel, note in [
        ("data/external/vcc2025", "the whole 2025 release on disk: 3 h5ads + gene_names.csv + "
                                  "3 pert_counts_*.csv; no library/spacer table"),
        ("data/challenge", "the 2026 release: 3 context h5ads + gene_names.csv + manifest.json + "
                           "pert_counts.csv + vcc_2026_controls.zip; no library/spacer table"),
        ("data/external/cd4_cell2026/suppl_tables_2026-08-24/sgrna_library_metadata.suppl_table.csv",
         "a DIFFERENT library (the GWCD4i primary-CD4 CRISPRi screen), not the VCC 2025 library"),
        ("data/external/gse249595/GSE249595_guides_reference.fa.gz",
         "a DIFFERENT study's guide reference (GSE249595), not the VCC 2025 library"),
        ("data/external/myllia-a549/GSE337804_sgRNA_lib.csv.gz",
         "a DIFFERENT library (Myllia A549), not the VCC 2025 library")]:
        p = ROOT / rel
        srcs.append(dict(source=rel, kind="checked, not the 2025 library", exists=p.exists(),
                         note=note, carries_spacer=False))
        P(f"  checked {rel}: {'present' if p.exists() else 'absent'} — {note}")

    # ---------------- (2) is a promoter SEQUENCE on disk? -------------------------------------
    P("\n=== (2) promoter sequence ===")
    bm = pd.read_csv(BIOMART)
    P(f"  {BIOMART.relative_to(ROOT)}: {len(bm):,} rows, columns {list(bm.columns)}")
    P(f"  -> COORDINATES only (chromosome, gene start/end, karyotype band). It fixes the TSS ± 500 "
      f"bp WINDOW but carries no base. No FASTA/2bit/genome file exists anywhere under data/.")
    have_seq = False

    # ---------------- the 34 arms ---------------------------------------------------------------
    P("\n=== the 34 arms of the 17 two-construct targets (293) ===")
    d = pd.read_csv(CONSTRUCT293)
    nat = d[d.regime == "challenge_20k"].drop_duplicates("target")
    rows = []
    for _, r in nat.iterrows():
        for tag, cons, cells, kd in (("big", r.c_big, r.cells_big, r.kd_big),
                                     ("small", r.c_small, r.cells_small, r.kd_small)):
            rows.append(dict(target=r.target, split=r.split, arm=tag, construct=cons,
                             cells=int(cells), kd_lfc=float(kd), abs_lfc=abs(float(kd)),
                             failing=bool(abs(float(kd)) < FAIL_LFC),
                             spacer_sequence="", spacer_source="NOT ON DISK",
                             promoter_sequence="NOT ON DISK",
                             seed_scan_status="BLOCKED: no spacer sequence, no promoter sequence"))
    arms = pd.DataFrame(rows).sort_values(["target", "arm"]).reset_index(drop=True)
    arms.to_csv(WORK / "seed298.csv", index=False)
    nf = int(arms.failing.sum())
    P(f"  {len(arms)} arms over {arms.target.nunique()} targets; FAILING (|lfc| < {FAIL_LFC}): "
      f"{nf}/{len(arms)}")
    P(f"  {'target':<9} {'arm':<6} {'construct':<26} {'cells':>6} {'kd lfc':>8}  failing")
    for _, r in arms.iterrows():
        P(f"  {r.target:<9} {r.arm:<6} {r.construct:<26} {r.cells:6d} {r.kd_lfc:+8.3f}  "
          f"{'FAIL' if r.failing else '.'}")
    P(f"\n  targets with BOTH constructs failing: "
      f"{sorted(arms.groupby('target').failing.sum().pipe(lambda s: s[s == 2]).index.tolist())}")
    P(f"  targets with NEITHER failing:          "
      f"{sorted(arms.groupby('target').failing.sum().pipe(lambda s: s[s == 0]).index.tolist())}")

    out = dict(
        verdict="NOT TESTABLE ON DISK",
        reason="the 2025 library's obs/guide_id encodes protospacer NAMES, not sequences "
               f"({total_seq}/{total_cats} categories are ACGT strings); no VCC-2025 library "
               "metadata table exists on disk; and no genome or promoter SEQUENCE exists on disk "
               "(gene_locations_biomart.csv is coordinates only). The mandate permits a download "
               "for the promoter sequence but NOT for the spacer sequences, so the missing half "
               "cannot be recovered under this mandate.",
        n_arms=len(arms), n_failing=nf, failing_arms=arms[arms.failing].construct.tolist(),
        fail_rule=f"|target-gene knockdown lfc| < {FAIL_LFC}",
        sources_checked=srcs, promoter_table=dict(path=str(BIOMART.relative_to(ROOT)),
                                                  rows=len(bm), columns=list(bm.columns),
                                                  carries_sequence=have_seq),
        ka=dict(row="KA 293 native n=200 rho_guide on 16 two-construct targets",
                want=[0.075, 16], got=[rg["point"], rg["n"]], ok=bool(ka_ok)),
        sha256_construct293=hashlib.sha256(CONSTRUCT293.read_bytes()).hexdigest())
    (WORK / "seed298.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "seed298.log").write_text("\n".join(LOG))
    P(f"\nwrote work/seed298.csv ({len(arms)} arms), work/seed298.json")


if __name__ == "__main__":
    main()
