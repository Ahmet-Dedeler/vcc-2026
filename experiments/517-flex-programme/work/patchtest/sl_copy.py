"""The same-line artifact for `VCC_SAMELINE` (451): per context, what a public bank of the SAME cell line says.

450 (H1, 2026-09-25) measured what a same-line bank from another lab and chemistry carries (Pan 2025, H1 hESC, 10x 5′,
against the 2025 VCC H1 Flex truth): block-sign programme 0.625 vs CLEAN4 0.581, per-target gate signs 0.764 vs our
0.551, per-target centred cosine 0.260 vs 0.069. So when the day's identity call (`def_identity.py`) names a line we hold
a bank for, that context should read the bank.

For each context whose consensus identity call is in REGISTRY, this writes (official gene axis, panel order):
  programme_ok_<c>  whether a bank in PROGRAMME_OK supplied the programme (else the builder keeps CLEAN4 there)
  w_<c>           the context's blend weight when the knob's w is `auto` (451b): --w-strong (2) if every bank that supplied
                  a covered row is in STRONG_BANKS, else --w-weak (0.5)
  programme_<c>   the bank's per-gene block sign: -1 where P_down > 0.5, +1 where P_down <= 0.5, -1 where the bank has no
                  value (396's `psign_A` rule: P_down = fraction of the bank's perturbations with delta < 0 among its
                  non-zero entries; the panel's own targets are excluded from the programme)
  pdown_<c>       P_down itself (NaN where uncovered)
  targets_<c>     indices (into the panel order) of the targets the bank covers
  rows_<c>        their RAW bank rows, unit-normalised over the non-panel genes, panel columns 0 (float32). `--centre`
                  subtracts the bank's own mean over its non-panel perturbations first; 450's emitter-level H1 price
                  prefers raw rows (w 1 + per-target prior: +0.0806 raw vs +0.0742 centred, board units), because the
                  centred rows lose mse and nmae for a small pds gain
  signs_<c>       sign of the RAW rows (int8; 0 where the bank does not measure the gene): per-target DE signs keep the
                  programme, which is informative for signs (450 §1-§2)
plus `contexts`, `genes`, `perts` and a JSON `provenance` (identity sample, unanimity, banks, bank file sha256s).
A context whose call is not in REGISTRY is simply absent: the builder leaves it exactly as without the knob.

`--programme-source <bank>` (460, default off; for the case where 458 reads positive, i.e. the lab-signature hypothesis
holds): EVERY context in the call gets that bank's programme (programme_<c>, pdown_<c>, the panel's rows excluded, 396's
rule) with programme_ok true. Matched contexts keep their same-line rows, targets, signs and w. The rest get 0 rows and w 0,
the construction `experiments/458-arc-lab-programme/artifact458.py` scored on A/B/C. The builder reads the programme only
without `noprog`, so the day runs sameline_flags `pt`. `h1_2025` = Arc's own H1 screen (Flex, 255 perturbations).

Rules boundary (AGENTS.md §7): public Perturb-seq conditioning a fitted transfer model, like every donor we ship. The
identity selects which public bank conditions the context; nothing here reads a literature effect for a specific target.

usage:
  uv run python scripts/def_sameline.py --bundle data/challenge --identity <identity.json> --out <work>/sameline.npz
  uv run python scripts/def_sameline.py --bundle data/challenge --force A=JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE \
      --out <work>/sameline_ka.npz            # known-answer: A's programme must equal 396's psign_A
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")  # 517 test copy only

# CCLE sample name -> banks, primary first. The primary bank alone defines the programme (the construction 450 and 396
# validated); a covered target reads the first bank in this order that has a row for it.
REGISTRY: dict[str, tuple[str, ...]] = {
    "HCT116_LARGE_INTESTINE": ("orion_hct116",),
    "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE": ("k562_gwps", "viperturb_k562", "k562_essential"),
    "DLD1_LARGE_INTESTINE": ("dld1_lowmoi",),
    "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE": ("nadig_jurkat",),
    "HEPG2_LIVER": ("nadig_hepg2",),
    "A549_LUNG": ("jiang_A549", "myllia_a549"),
    "HT29_LARGE_INTESTINE": ("jiang_HT29",),
    "MCF7_BREAST": ("jiang_MCF7",),
    "BXPC3_PANCREAS": ("jiang_BXPC3",),
    "HAP1_ENGINEERED": ("jiang_HAP1",),
    "CALU3_LUNG": ("sunshine_calu3",),
    "THP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE": ("wessels_thp1",),
}
# Banks whose per-gene programme may replace CLEAN4: unstimulated, broad CRISPRi screens (the kind 450 validated: Pan's
# genome-scale H1 screen, and its essential half). Targeted, cytokine-stimulated or infection screens (Jiang, Sunshine,
# Wessels, Myllia) contribute per-target rows only: their mean response is the stimulation's, not the line's.
PROGRAMME_OK = {"orion_hct116", "k562_gwps", "viperturb_k562", "k562_essential", "dld1_lowmoi", "nadig_jurkat", "nadig_hepg2"}
# 451b: bank quality classes for the per-context blend weight (`w_<c>`, read by the builder when VCC_SAMELINE's w is
# `auto`). STRONG = own-target knockdown visible in >= 0.9 of rows (451's bank audit); everything else is weak. 450f (H1):
# the best w is 2 for a Pan-quality bank and 0.5 for a 0.61-sign bank, where w 2 is neutral overall and costs mse.
STRONG_BANKS = {"orion_hct116", "k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2", "sunshine_calu3"}
BANK_DIRS: dict[str, Path] = {
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "viperturb_k562": ROOT / "data/donors/viperturb_k562",
    "k562_essential": ROOT / "data/donors/k562_essential",
    "dld1_lowmoi": ROOT / "experiments/223-silent-zero-audit/work/tables/dld1_lowmoi",
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
    "nadig_hepg2": ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2",
    "jiang_A549": ROOT / "experiments/225-source-sweep/work/tables/jiang_A549",
    "myllia_a549": ROOT / "experiments/120-myllia-a549/work/compact-donor",
    "jiang_HT29": ROOT / "experiments/225-source-sweep/work/tables/jiang_HT29",
    "jiang_MCF7": ROOT / "experiments/225-source-sweep/work/tables/jiang_MCF7",
    "jiang_BXPC3": ROOT / "experiments/225-source-sweep/work/tables/jiang_BXPC3",
    "jiang_HAP1": ROOT / "experiments/225-source-sweep/work/tables/jiang_HAP1",
    "sunshine_calu3": ROOT / "experiments/225-source-sweep/work/tables/sunshine_calu3",
    "wessels_thp1": ROOT / "experiments/225-source-sweep/work/tables/wessels_thp1",
    "h1_2025": ROOT / "data/donors/h1_2025",          # 460: --programme-source only (Arc's lab, 458); not in REGISTRY
    # 458d: Arc's H1 from ALL public cells (300 perts x ~1,071 cells vs 38,176 controls; denoise458d.py): the source to use
    "h1_2025_full": ROOT / "data/donors/h1_2025_full",
}


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def read_column(path: Path, name: str) -> list[str]:
    with open(path, newline="") as f:
        return [row[name] for row in csv.DictReader(f)]


def load_bank(name: str):
    d = BANK_DIRS[name]
    delta = np.load(d / "delta.npy", mmap_mode="r")
    perts = [line.strip() for line in open(d / "perts.txt")]
    genes = [line.strip() for line in open(d / "genes.txt")]
    if delta.shape != (len(perts), len(genes)):
        raise SystemExit(f"{name}: delta {delta.shape} != perts x genes {(len(perts), len(genes))}")
    return delta, perts, genes, file_sha256(d / "delta.npy")


ESSENTIAL_LIST = ROOT / "data/donors/k562_essential/perts.txt"   # 517: meanz_ne drops these rows


def programme_meanz_ne(delta, perts, genes, official, panel) -> tuple[np.ndarray, np.ndarray, int]:
    """517 (artifact517.py): mean over the non-essential, non-panel rows of delta / sd (sd over all non-panel rows),
    as a pseudo-P whose |P - 0.5| >= 0.05 coverage equals the P_down rule's (so the hook's tau keeps its meaning)."""
    ess = {l.strip() for l in open(ESSENTIAL_LIST) if l.strip()}
    keep = np.array([p not in panel for p in perts])
    X = np.asarray(delta[np.flatnonzero(keep)], dtype=np.float64)
    ne = np.array([p not in ess for p in np.asarray(perts)[keep]])
    nz = (X != 0).sum(0); dn = (X < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pd_all = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    pos = {g: j for j, g in enumerate(genes)}
    col = np.array([pos.get(g, -1) for g in official])
    P_all = np.where(col >= 0, pd_all[np.maximum(col, 0)], np.nan)
    sd = X.std(0); sd[sd == 0] = 1.0
    zb = (X[ne] / sd).mean(0)
    z = np.where(np.isfinite(P_all) & (col >= 0), zb[np.maximum(col, 0)], np.nan)
    n_cov = int((np.isfinite(P_all) & (np.abs(P_all - 0.5) >= 0.05)).sum())
    z_q = float(np.sort(np.abs(np.nan_to_num(z)))[::-1][n_cov - 1])
    pdown = np.where(np.isfinite(z), np.clip(0.5 - 0.05 * z / z_q, 0.001, 0.999), np.nan)
    psign = np.where(np.isfinite(pdown) & (pdown <= 0.5), 1.0, -1.0)
    return psign, pdown, int(ne.sum())


def programme(delta, perts, genes, official, panel) -> tuple[np.ndarray, np.ndarray, int, np.ndarray]:
    """396's rule on the official axis: P_down over non-zero entries, panel targets' rows excluded. Also returns the bank's
    mean delta over the same rows (its own gene axis), which the per-target rows are centred by."""
    keep = np.array([p not in panel for p in perts])
    dn = np.zeros(len(genes))
    nz = np.zeros(len(genes))
    tot = np.zeros(len(genes))
    rows = np.where(keep)[0]
    for s in range(0, len(rows), 512):
        blk = np.asarray(delta[rows[s:s + 512]], dtype=np.float64)
        dn += (blk < 0).sum(0)
        nz += (blk != 0).sum(0)
        tot += blk.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pd_bank = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    pos = {g: j for j, g in enumerate(genes)}
    col = np.array([pos.get(g, -1) for g in official])
    pdown = np.where(col >= 0, pd_bank[np.maximum(col, 0)], np.nan)
    psign = np.where(np.isfinite(pdown) & (pdown <= 0.5), 1.0, -1.0)
    return psign, pdown, int(keep.sum()), tot / max(int(keep.sum()), 1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--identity", help="def_identity.py's JSON (uses its consensus `call` per context)")
    ap.add_argument("--force", default="", help="C=CCLE_SAMPLE[,..]: bypass the identity file (known-answer rows only)")
    ap.add_argument("--contexts", default="", help="restrict to these contexts (default: every context in the call)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--centre", action="store_true", help="centre each bank row by the bank's own mean (not the default)")
    ap.add_argument("--w-strong", type=float, default=2.0, help="w_<c> when every row-supplying bank is strong (450d)")
    ap.add_argument("--w-weak", type=float, default=0.5, help="w_<c> when any row-supplying bank is weak (450f)")
    ap.add_argument("--min-programme-rows", type=int, default=150,
                    help="481 F7: refuse a --programme-source with fewer rows left after the panel is excluded")
    ap.add_argument("--programme-estimator", choices=["pdown", "meanz_ne"], default="pdown",
                    help="517: meanz_ne = the hook-source encoding for --programme-source (not for the block prior)")
    ap.add_argument("--programme-source", default="",
                    help="460: a bank whose programme becomes EVERY context's block prior (h1_2025 = Arc's H1, 458); "
                         "matched contexts keep their same-line rows; default off")
    a = ap.parse_args()
    t0 = time.time()
    bundle = Path(a.bundle)
    official = read_column(bundle / "gene_names.csv", "gene_name")
    perts = read_column(bundle / "pert_counts.csv", "target_gene")
    panel = set(perts)
    panel_cols = np.array([g in panel for g in official])
    if a.force:
        call = dict(p.split("=", 1) for p in a.force.split(","))
        ident = {"call": call, "call_unanimous_across_methods": {c: None for c in call}, "source": "--force"}
    else:
        ident = json.loads(Path(a.identity).read_text())
        call = dict(ident["call"])
    if a.contexts:
        call = {c: v for c, v in call.items() if c in set(a.contexts.split(","))}
    out: dict[str, np.ndarray] = {"genes": np.asarray(official), "perts": np.asarray(perts)}
    prov: dict[str, object] = {"identity_source": str(a.identity or "--force"), "contexts": {}}
    matched = []
    for c, sample in sorted(call.items()):
        banks = REGISTRY.get(str(sample))
        rec: dict[str, object] = {"call": sample, "unanimous": ident.get("call_unanimous_across_methods", {}).get(c)}
        if banks and not a.force and rec["unanimous"] is not True:
            # 041 calibrated the matcher on known controls; a split call can be a near neighbour (a sister line of the
            # same tissue), where a bank is no longer same-line information. Only a unanimous call switches the mode on.
            rec["mode"] = "none (registered line, but the call is not unanimous across the 8 methods)"
            prov["contexts"][c] = rec
            print(f"[sameline] {c}: {sample} is registered but the call is not unanimous; context unchanged", flush=True)
            continue
        if not banks:
            rec["mode"] = "none (no registered bank)"
            prov["contexts"][c] = rec
            print(f"[sameline] {c}: {sample} -> no registered bank; context unchanged", flush=True)
            continue
        loaded = {b: load_bank(b) for b in banks}
        pbank = next((b for b in banks if b in PROGRAMME_OK), None)
        if pbank is not None:
            psign, pdown, n_prog, _ = programme(*loaded[pbank][:3], official, panel)
        else:
            psign, pdown, n_prog = -np.ones(len(official)), np.full(len(official), np.nan), 0
        centre = {b: programme(*loaded[b][:3], official, panel)[3] for b in banks} if a.centre else {}
        tgt_idx, rows, signs, src = [], [], [], []
        for i, t in enumerate(perts):
            for b in banks:
                delta, bperts, bgenes, _ = loaded[b]
                if t in bperts:
                    r = np.asarray(delta[bperts.index(t)], dtype=np.float64)
                    pos = {g: j for j, g in enumerate(bgenes)}
                    col = np.array([pos.get(g, -1) for g in official])
                    raw = np.where(col >= 0, r[np.maximum(col, 0)], 0.0)
                    v = np.where(col >= 0, ((r - centre[b]) if a.centre else r)[np.maximum(col, 0)], 0.0)
                    raw[panel_cols] = 0.0
                    v[panel_cols] = 0.0
                    n = np.linalg.norm(v)
                    if n > 0:
                        tgt_idx.append(i); rows.append((v / n).astype(np.float32)); signs.append(np.sign(raw).astype(np.int8))
                        src.append(b)
                    break
        matched.append(c)
        out[f"programme_{c}"] = psign
        out[f"programme_ok_{c}"] = np.asarray(pbank is not None)
        weak = sorted({b for b in src if b not in STRONG_BANKS})
        out[f"w_{c}"] = np.asarray(a.w_weak if weak else a.w_strong)
        out[f"pdown_{c}"] = pdown
        out[f"targets_{c}"] = np.asarray(tgt_idx, dtype=np.int64)
        out[f"rows_{c}"] = np.asarray(rows, dtype=np.float32).reshape(len(rows), len(official))
        out[f"signs_{c}"] = np.asarray(signs, dtype=np.int8).reshape(len(signs), len(official))
        rec.update({
            "mode": "same-line", "banks": list(banks),
            "bank_delta_sha256": {b: loaded[b][3] for b in banks},
            "programme_bank": pbank, "programme_rows": n_prog, "rows_centred": bool(a.centre),
            "w": float(out[f"w_{c}"]), "weak_row_banks": weak,
            "programme_coverage": float(np.isfinite(pdown).mean()),
            "programme_frac_down": float((psign < 0).mean()),
            "covered_targets": len(tgt_idx),
            "covered_by_bank": {b: src.count(b) for b in banks},
        })
        prov["contexts"][c] = rec
        print(f"[sameline] {c}: {sample} -> {list(banks)}; programme from {pbank or 'NONE (stimulated/targeted banks only)'} ({n_prog} rows, coverage "
              f"{rec['programme_coverage']:.3f}); {len(tgt_idx)}/{len(perts)} panel targets covered {rec['covered_by_bank']}",
              flush=True)
    if a.programme_source:
        if a.programme_source not in BANK_DIRS:
            raise SystemExit(f"--programme-source {a.programme_source!r} is not in BANK_DIRS")
        sdelta, sperts, sgenes, ssha = load_bank(a.programme_source)
        spsign, spdown, sn, _ = programme(sdelta, sperts, sgenes, official, panel)
        if a.programme_estimator == "meanz_ne":
            spsign, spdown, sn_ne = programme_meanz_ne(sdelta, sperts, sgenes, official, panel)
            print(f"[sameline] --programme-estimator meanz_ne: {sn_ne} of {sn} rows (non-essential)", flush=True)
        scov = float(np.isfinite(spdown).mean())
        if sn < a.min_programme_rows or scov < 0.5:
            # 481 F7: the programme excludes the panel's own rows; a panel that overlaps the source leaves too few, and at
            # 0 rows every gene silently defaults to down while programme_ok stays true
            raise SystemExit(f"--programme-source {a.programme_source}: {sn} rows left after excluding the panel "
                             f"(gene coverage {scov:.3f}); need >= {a.min_programme_rows} rows and coverage >= 0.5")
        for c in sorted(call):
            out[f"programme_{c}"] = spsign
            out[f"pdown_{c}"] = spdown
            out[f"programme_ok_{c}"] = np.asarray(True)
            if c not in matched:
                out[f"w_{c}"] = np.asarray(0.0)
                out[f"targets_{c}"] = np.zeros(0, dtype=np.int64)
                out[f"rows_{c}"] = np.zeros((0, len(official)), dtype=np.float32)
                out[f"signs_{c}"] = np.zeros((0, len(official)), dtype=np.int8)
                matched.append(c)
            prov["contexts"].setdefault(c, {"call": call[c]})["programme_source"] = {
                "bank": a.programme_source, "bank_delta_sha256": ssha, "programme_rows": sn,
                **({} if a.programme_estimator == "pdown" else {"programme_estimator": a.programme_estimator}),
                "programme_coverage": float(np.isfinite(spdown).mean()), "programme_frac_down": float((spsign < 0).mean())}
        matched = sorted(matched)
        print(f"[sameline] --programme-source {a.programme_source}: every context's programme from it ({sn} rows, coverage "
              f"{float(np.isfinite(spdown).mean()):.3f}); contexts {matched}", flush=True)
    out["contexts"] = np.asarray(matched)
    out["provenance"] = np.asarray(json.dumps(prov))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, **out)
    sha = file_sha256(Path(a.out))
    Path(str(a.out) + ".json").write_text(json.dumps({"sha256": sha, **prov}, indent=1))
    print(f"[sameline] wrote {a.out} sha256={sha} contexts={matched} [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
