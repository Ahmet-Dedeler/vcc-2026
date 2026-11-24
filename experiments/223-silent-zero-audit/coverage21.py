"""223 Phase 1b — COMPUTED panel coverage for every one of the 21 experiments' sources.

The mandate requires coverage computed, not eyeballed.  For each experiment the source's
own target list is read off disk and intersected with `data/challenge/pert_counts.csv`.
Where the source has no local artifact (values sealed, or never acquired), that is stated
explicitly and the experiment is classified UNKNOWN-COVERAGE rather than guessed.

Three coverage numbers matter and they are different:
  * `panel_targets`   -- rows of the 300-gene panel this source carries.  This is what the
                         MEDIAN sizing rule reads: <= 50% and the median row norm is 0.
  * `gene_axis`       -- response coordinates the source measures.  116's NO-GO was a
                         GENE-axis coverage failure (13.91% of the non-panel axis), a
                         different thing entirely.
  * `fold_coverage`   -- the same intersection restricted to each PLOLO fold's targets.

    uv run python coverage21.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

PANEL = [line.split(",")[0]
         for line in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
         if line.strip()]
assert len(PANEL) == 300
PANEL_SET = set(PANEL)


def lines(p: Path) -> list[str]:
    return [x.strip() for x in p.read_text().splitlines() if x.strip()]


def from_perts(rel: str) -> dict:
    p = ROOT / rel
    perts = set(lines(p / "perts.txt"))
    genes = lines(p / "genes.txt")
    return dict(kind="donor-table", path=rel, n_targets=len(perts),
                n_genes=len(genes), panel_targets=len(perts & PANEL_SET))


def from_symbol_file(rel: str, kind: str, n_genes=None) -> dict:
    p = ROOT / rel
    syms = set(lines(p))
    return dict(kind=kind, path=rel, n_targets=len(syms), n_genes=n_genes,
                panel_targets=len(syms & PANEL_SET))


def dld1() -> dict:
    p = ROOT / "experiments/116-dld1-lowmoi/work/compact-source-v1"
    syms = set(lines(p / "target_symbols.txt"))
    genes = lines(p / "canonical_genes.txt")
    mask = np.load(p / "panel_response_mask.npy")
    return dict(kind="compact-source (npy, no perts.txt)", path=str(p.relative_to(ROOT)),
                n_targets=len(syms), n_genes=len(genes),
                panel_targets=len(syms & PANEL_SET),
                panel_response_coords_zeroed=int(np.count_nonzero(mask)))


def cmap() -> dict:
    """176's per-line CGS matrices are already SUBSET to h1-test's 85 targets, so the
    panel coverage has to come from its own inventory table, which counts distinct target
    genes per line over the whole L1000 CGS block.  Parsed, not eyeballed."""
    inv = (ROOT / "experiments/176-cmap-donors/INVENTORY.md").read_text()
    per = {}
    for ln in inv.splitlines():
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if len(cells) < 5 or not cells[3].replace("*", "").strip().rstrip("/300").isdigit():
            continue
        try:
            per[f"{cells[0]}:{cells[1]}"] = dict(
                n_targets=int(cells[2]),
                panel_targets=int(cells[3].replace("*", "").strip()))
        except ValueError:
            continue
    union = 127  # INVENTORY.md's own union line, quoted below and asserted present
    assert "Union panel coverage across all lines with" in inv and "127/300" in inv
    d = ROOT / "experiments/176-cmap-donors/work"
    scored = {}
    for f in sorted(d.glob("cgs_*.npz")):
        with np.load(f, allow_pickle=True) as z:
            scored[f.stem] = dict(scored_targets=int(len(z["perts"])),
                                  scored_genes=int(len(z["genes"])),
                                  have=int(np.count_nonzero(z["have"])))
    return dict(kind="L1000 CGS; scored on h1-test's 85-target DE pool, NOT the panel",
                path="experiments/176-cmap-donors/work", per_line=per,
                scored_matrices=scored, panel_targets=union, n_targets=None,
                n_genes=None,
                note="union over the 9 core lines; no single line exceeds 118/300")


def hcasmc() -> dict:
    p = ROOT / "experiments/126-hcasmc-panel-confirmation/VCC_TARGET_OVERLAP.v1.txt"
    syms = set(lines(p))
    return dict(kind="overlap list only (effect values SEALED, not opened)",
                path=str(p.relative_to(ROOT)), n_targets=None, n_genes=None,
                panel_targets=len(syms & PANEL_SET))


def hues8() -> dict:
    p = ROOT / "experiments/131-hues8-hepatocyte-metadata/vcc_overlap.csv"
    rows = [r.split(",") for r in p.read_text().splitlines()[1:] if r.strip()]
    hdr = p.read_text().splitlines()[0].split(",")
    col = hdr.index("target_gene")
    per = {}
    for r in rows:
        per.setdefault(r[0].strip(), set()).add(r[col].strip())
    syms = set().union(*per.values()) & PANEL_SET
    return dict(kind="metadata only (NEVER ACQUIRED -- HTTP 403/404 on every object)",
                path=str(p.relative_to(ROOT)), n_targets=1950, n_genes=None,
                panel_targets=len(syms),
                per_dataset={k: len(v & PANEL_SET) for k, v in per.items()})


def jurkat107() -> dict:
    j = json.loads((ROOT / "experiments/107-jurkat-high-moi/work/source-model-input.json")
                   .read_text())
    z = np.load(ROOT / "experiments/107-jurkat-high-moi/work/source-model-input.npz",
                allow_pickle=True)
    ts = np.asarray(z["target_symbols"])
    elig = np.asarray(z["precomputed_eligible"], dtype=bool)
    all_t = set(map(str, ts))
    eligible = set(map(str, ts[elig]))
    return dict(kind="joint high-MOI deconvolution; estimable target set",
                path="experiments/107-jurkat-high-moi/work/source-model-input.npz",
                n_targets=len(eligible), n_genes=int(j["modeled_response_genes"]),
                panel_targets=len(eligible & PANEL_SET),
                guide_library_targets=len(all_t),
                guide_library_panel_overlap=len(all_t & PANEL_SET),
                note=("TARGET coverage 299/300 -> the median row-norm rule is LIVE; the "
                      "GENE axis is only 349/18,533 = 1.9%, a different limitation"))


SOURCES = {
    "055-kolf-direct": lambda: from_perts("data/external/kolf2"),
    "057-kolf-union-whiten": lambda: from_perts("data/external/kolf2"),
    "059-feng-ipsc": lambda: from_perts("experiments/099-donors/work/feng_ipsc"),
    "079-kolf-local-centering": lambda: from_perts("data/external/kolf2"),
    "085-kolf-pm-gate": lambda: from_perts("data/external/kolf2"),
    "105-telohaec-residual": lambda: from_perts("data/external/igvf-telohaec"),
    "107-jurkat-high-moi": jurkat107,
    "115-gse295214-h1": lambda: from_perts("experiments/115-gse295214-h1/work/compact-donor"),
    "116-dld1-lowmoi": dld1,
    "119-altair-a549": lambda: from_perts("experiments/119-altair-a549/work/compact-source"),
    "120-myllia-a549": lambda: from_perts("experiments/120-myllia-a549/work/compact-donor"),
    "121-a549-cross-source-concordance": lambda: dict(
        kind="two compact sources, 20-target common gallery",
        path="119 compact-source AND 120 compact-donor",
        n_targets=20, n_genes=7685,
        panel_targets=len(set(lines(ROOT / "experiments/119-altair-a549/work/compact-source"
                                    / "perts.txt")) & PANEL_SET
                          & set(lines(ROOT / "experiments/120-myllia-a549/work/compact-donor"
                                      / "perts.txt")))),
    "126-hcasmc-panel-confirmation": hcasmc,
    "129-gse295214-transfer": lambda: from_perts(
        "experiments/115-gse295214-h1/work/compact-donor"),
    "130-telohaec-panel-confirmation": lambda: from_perts("data/external/igvf-telohaec"),
    "131-hues8-hepatocyte-metadata": hues8,
    "135-telohaec-atlas-confirmation": lambda: from_perts("data/external/igvf-telohaec"),
    "137-jurkat-tangent-cone": jurkat107,
    "138-dld1-replicate-sign": dld1,
    "139-jurkat-donor-simplex": jurkat107,
    "176-cmap-donors": cmap,
}


def main() -> None:
    out = {}
    for exp, fn in SOURCES.items():
        try:
            r = fn()
        except Exception as e:  # noqa: BLE001
            r = dict(kind="ERROR", error=f"{type(e).__name__}: {e}")
        pt = r.get("panel_targets")
        r["panel_pct"] = (round(100 * pt / 300, 1) if isinstance(pt, int) else None)
        r["median_rule_vacuous"] = (None if not isinstance(pt, int) else bool(pt <= 150))
        out[exp] = r
    (WORK / "coverage21.json").write_text(json.dumps(out, indent=2))

    hdr = f"{'experiment':38s}{'panel':>7s}{'pct':>7s}{'targets':>9s}{'genes':>7s}  median-rule"
    print(hdr)
    print("-" * len(hdr))
    for exp, r in out.items():
        pt = r.get("panel_targets")
        vac = r.get("median_rule_vacuous")
        tag = ("VACUOUS (<=50%)" if vac is True else
               "live (>50%)" if vac is False else "UNKNOWN")
        print(f"{exp:38s}{(pt if pt is not None else '?'):>7}"
              f"{(r['panel_pct'] if r['panel_pct'] is not None else '?'):>7}"
              f"{(r.get('n_targets') or '?'):>9}{(r.get('n_genes') or '?'):>7}  {tag}")
    print(f"\n-> {WORK / 'coverage21.json'}")


if __name__ == "__main__":
    main()
