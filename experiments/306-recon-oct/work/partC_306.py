"""306 Part C — known-answer check on the one numeric claim a rival published.

Claim (public, leaderboard `description` field, entry `v6loco_plus_H1hESC_nb_a2`,
team "Testh", score_avg 0.14954, read live 2026-09-07T15:48Z):
    "... with the VCC2025 H1 hESC context added to the evidence pool (21/300 panel targets)."

Repo's own established number (READ-FIRST §6, 285's 16:15Z check):
    "the 2026 panel shares only 25/300 targets with the 2025 H1 library
     (`data/challenge/pert_counts.csv` vs the 2025 train/val/test targets)"

KNOWN-ANSWER ROW: recompute 25/300 from the two files named in that sentence.
Then enumerate the other defensible readings of "the 2025 H1 library" to see whether
21 is one of them.

Read-only. No writes outside experiments/306-recon-oct/. Tiny CSVs only; RSS < 200 MB.
"""

import csv
import hashlib
import pathlib

ROOT = pathlib.Path("~/Code/virtual-cell-challenge-2026")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def col(path, name):
    """Read one named column out of a CSV as a list of stripped strings."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    assert name in rows[0], (path, list(rows[0]))
    return [r[name].strip() for r in rows if r[name].strip()]


def lines(path):
    return [l.strip() for l in open(path) if l.strip()]


# ---------------------------------------------------------------- the 2026 panel
panel_path = ROOT / "data/challenge/pert_counts.csv"
panel = col(panel_path, "target_gene")
panel_set = set(panel)
print(f"2026 panel      {panel_path}  sha256[:16]={sha(panel_path)}")
print(f"                rows={len(panel)}  unique={len(panel_set)}")

# --------------------------------------------- reading 1: the 2025 train/val/test
# The sentence under test names exactly these files.
tvt = {}
for split, fn in [
    ("train", "train/pert_counts_Training.csv"),
    ("validation", "validation/pert_counts_Validation.csv"),
    ("test", "test/pert_counts_Test.csv"),
]:
    p = ROOT / "data/external/vcc2025" / fn
    with open(p, newline="") as f:
        hdr = next(csv.reader(f))
    tgt_col = "target_gene" if "target_gene" in hdr else hdr[0]
    g = [x for x in col(p, tgt_col) if x.lower() not in {"non-targeting", "ntc", "control"}]
    tvt[split] = set(g)
    print(f"2025 {split:<11} {p.name}  sha256[:16]={sha(p)}  header={hdr}  n={len(g)}  unique={len(set(g))}")

union_2025 = tvt["train"] | tvt["validation"] | tvt["test"]
print(f"\n2025 train/val/test union: {len(union_2025)} targets")

ov_tvt = panel_set & union_2025
print(f"READING 1  panel n 2025 train/val/test  =  {len(ov_tvt)} / {len(panel_set)}")
print(f"           {sorted(ov_tvt)}")
for split in ("train", "validation", "test"):
    print(f"           of which in {split}: {len(panel_set & tvt[split])}")

# ------------------------------- reading 2: the h1_2025 donor bank actually built
p = ROOT / "data/donors/h1_2025/perts.txt"
donor = set(lines(p))
print(f"\nh1_2025 donor bank {p} sha256[:16]={sha(p)} n={len(donor)}")
ov_donor = panel_set & donor
print(f"READING 2  panel n h1_2025 donor rows   =  {len(ov_donor)} / {len(panel_set)}")

# ------------------------- reading 3: the 2025 library as resolved spacers (301)
p = ROOT / "experiments/301-metadata-acquisitions/work/vcc2025_spacers.csv"
with open(p, newline="") as f:
    rows = list(csv.DictReader(f))
spacer_all = {r["guide_id"].split("_")[0] for r in rows if r["guide_id"].strip()}
spacer_res = {r["gene"].strip() for r in rows if r["gene"].strip()}
print(f"\n301 spacer table {p.name} sha256[:16]={sha(p)} rows={len(rows)}"
      f"  guide_id-derived genes={len(spacer_all)}  resolved gene column={len(spacer_res)}")
print(f"READING 3a panel n 2025 library (all constructs)      = {len(panel_set & spacer_all)} / {len(panel_set)}")
print(f"READING 3b panel n 2025 library (spacer-resolved only) = {len(panel_set & spacer_res)} / {len(panel_set)}")

# ---------------------------------------------------------------------- verdict
print("\n" + "=" * 72)
print("KNOWN-ANSWER ROW")
print(f"  expected (READ-FIRST §6 / 285, 16:15Z): 25 / 300")
print(f"  measured (reading 1, the files that sentence names): {len(ov_tvt)} / {len(panel_set)}")
print(f"  verdict: {'PASS' if (len(ov_tvt), len(panel_set)) == (25, 300) else 'FAIL'}")
print(f"\nRIVAL CLAIM under test: 21 / 300 (Testh, v6loco_plus_H1hESC_nb_a2)")
cands = {
    "1  train+val+test union": len(ov_tvt),
    "2  h1_2025 donor bank": len(ov_donor),
    "3a 2025 library, all constructs": len(panel_set & spacer_all),
    "3b 2025 library, spacer-resolved": len(panel_set & spacer_res),
    "1t train only": len(panel_set & tvt["train"]),
    "1v validation only": len(panel_set & tvt["validation"]),
    "1e test only": len(panel_set & tvt["test"]),
}
for k, v in cands.items():
    print(f"  {k:<36} {v:>3} / 300   {'<-- matches the rival number' if v == 21 else ''}")
