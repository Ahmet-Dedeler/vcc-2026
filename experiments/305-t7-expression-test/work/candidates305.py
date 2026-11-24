"""305 step 1 -- rebuild 301's FROZEN candidate table from its saved hits, with 301's exact filters.

Reads (read-only) 301's `seed_hits301_cage.csv`, `promoters_cage_tss_hg19.fa.gz`,
`vcc2025_spacers.csv` and 298's `seed298.csv`.  Reproduces the primary reading of
`seed_contrast.py` -- 12-nt seed, CRISPRi active window (TSS -50 .. +300), DISTAL filter
(candidate promoter window >= 10 kb from any window of the arm's own gene) -- and asserts the
published 4/13, 4/21, 4/31 counts.  Emits EVERY evaluated construct-gene pair (failing arms,
working arms AND the 31 non-targeting constructs), duplicate spacers/hits collapsed to one
construct-gene pair.

    PYTHONDONTWRITEBYTECODE=1 uv run python candidates305.py
"""
from __future__ import annotations

import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path("~/Code/virtual-cell-challenge-2026")
E301 = ROOT / "experiments/301-metadata-acquisitions/work"
E298 = ROOT / "experiments/298-reviewer-checks/work"
HERE = Path(__file__).resolve().parent

L_PRIMARY = 12
WIN = (-50, 300)          # 301's CRISPRi active window
DISTAL_KB = 10            # 301's self-locus filter


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> None:
    # ---- 301's promoter windows, for the distal filter (verbatim logic from seed_contrast.py)
    wins: dict[str, list] = collections.defaultdict(list)
    with gzip.open(E301 / "promoters_cage_tss_hg19.fa.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                g, t, c, s, tss, lo, hi, n = line[1:].strip().split("|")
                wins[g].append((c, int(lo), int(hi)))

    def distal(tgt: str, cg: str, kb: int = DISTAL_KB) -> bool:
        if not tgt:
            return True
        for c1, l1, h1 in wins.get(tgt, []):
            for c2, l2, h2 in wins.get(cg, []):
                if c1 == c2 and max(l2 - h1, l1 - h2, 0) < kb * 1000:
                    return False
        return True

    hits = list(csv.DictReader(open(E301 / "seed_hits301_cage.csv")))
    for h in hits:
        h["L"] = int(h["L"])
        h["off"] = int(h["offset_from_tss"])

    arms = list(csv.DictReader(open(E298 / "seed298.csv")))
    fail = [a["construct"] for a in arms if a["failing"] == "True"]
    work = [a["construct"] for a in arms if a["failing"] != "True"]
    split_of = {a["construct"]: a["split"] for a in arms}
    target_of = {a["construct"]: a["target"] for a in arms}
    kd_of = {a["construct"]: float(a["kd_lfc"]) for a in arms}
    cells_of = {a["construct"]: int(a["cells"]) for a in arms}
    ntc = [r["guide_id"] for r in csv.DictReader(open(E301 / "vcc2025_spacers.csv"))
           if r["guide_id"].startswith("non-targeting") and r["route"] != "UNRESOLVED"]
    assert (len(fail), len(work), len(ntc)) == (13, 21, 31), (len(fail), len(work), len(ntc))

    dropped = [h for h in hits if not distal(h["target"], h["gene"])]
    kept = [h for h in hits if distal(h["target"], h["gene"])]
    sub = [h for h in kept if h["L"] == L_PRIMARY and WIN[0] <= h["off"] <= WIN[1]]

    cand: dict[str, dict[str, list]] = collections.defaultdict(lambda: collections.defaultdict(list))
    for h in sub:
        cand[h["construct"]][h["gene"]].append(h)

    counts = {
        "failing": [sum(1 for c in fail if cand.get(c)), len(fail)],
        "working": [sum(1 for c in work if cand.get(c)), len(work)],
        "ntc": [sum(1 for c in ntc if cand.get(c)), len(ntc)],
    }
    print("301 primary reading reproduced (L=12, CRISPRi window, distal):")
    for k, (a, b) in counts.items():
        print(f"  {k:9s} {a:2d} / {b}")
    assert counts["failing"] == [4, 13] and counts["working"] == [4, 21] and counts["ntc"] == [4, 31], counts
    print(f"  self-locus hits dropped: {len(dropped)} over "
          f"{len({(h['construct'], h['gene']) for h in dropped})} construct-gene pairs")

    rows = []
    for grp, members in (("failing", fail), ("working", work), ("ntc", ntc)):
        for c in members:
            for g, hs in sorted(cand.get(c, {}).items()):
                rows.append(dict(
                    pair_id=f"{c}|{g}",
                    group=grp,
                    construct=c,
                    intended_target=target_of.get(c, ""),
                    split=split_of.get(c, "ALL"),
                    candidate_gene=g,
                    n_hits=len(hs),
                    spacers=",".join(sorted({h["spacer"] for h in hs})),
                    offsets=",".join(str(h["off"]) for h in sorted(hs, key=lambda x: x["off"])),
                    orients=",".join(sorted({h["orient"] for h in hs})),
                    pams=",".join(sorted({h["pam"] for h in hs})),
                    kd_lfc=kd_of.get(c, ""),
                    cells_full=cells_of.get(c, ""),
                ))
    with open(HERE / "candidates305.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    meta = dict(
        seed_length=L_PRIMARY, crispri_window=list(WIN), distal_kb=DISTAL_KB,
        counts=counts, n_pairs=len(rows),
        n_pairs_by_group={g: sum(1 for r in rows if r["group"] == g) for g in ("failing", "working", "ntc")},
        primary_family=[r["pair_id"] for r in rows if r["group"] in ("failing", "working")],
        comparator_family=[r["pair_id"] for r in rows if r["group"] == "ntc"],
        inputs={str(p.relative_to(ROOT)): dict(bytes=p.stat().st_size, sha256=sha256_file(p)) for p in (
            E301 / "seed_hits301_cage.csv", E301 / "promoters_cage_tss_hg19.fa.gz",
            E301 / "vcc2025_spacers.csv", E298 / "seed298.csv")},
    )
    (HERE / "candidates305_meta.json").write_text(json.dumps(meta, indent=1))
    print(f"\nwrote candidates305.csv: {len(rows)} construct-gene pairs "
          f"({meta['n_pairs_by_group']})")
    for r in rows:
        print(f"  {r['group']:8s} {r['construct']:32s} -> {r['candidate_gene']:10s} "
              f"({r['split']}, {r['n_hits']} hit(s), offsets {r['offsets']})")


if __name__ == "__main__":
    main()
