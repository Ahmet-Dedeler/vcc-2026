"""481: donor/source coverage of every panel we hold, from the banks' perts lists alone (no matrices, < 0.3 GB).

Panels: A/B/C (Arc 2026 validation), the 2025 VCC H1 panel (Arc's other 300-target panel: the best stand-in we have for
how Arc builds a panel), and the synthetic October draws (238 efr/efp, 199 realistic/pessimistic).
Sources: the three production donors (carrier histogram 3/2/1/0 = the builder's `ncar`), the shipped H1 union pair, the
programme source (h1_2025_full: its programme excludes the panel's own rows), PerturbMatch, and the same-line registry
banks. CD4 per-condition coverage comes from extract_cd4.py's dry runs (RESULTS.md F1), not from here.

usage: uv run python census481.py  -> work/census481.json and a table on stdout
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PANELS = ["abc", "h1_2025", "efr", "efp", "real199", "pess199"]
BANKS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps/perts.txt",
    "orion_hct116": ROOT / "data/donors/orion_hct116/perts.txt",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t/perts.txt",
    "gse295214_h1": ROOT / "experiments/115-gse295214-h1/work/compact-donor/perts.txt",
    "h1_2025": ROOT / "data/donors/h1_2025/perts.txt",
    "h1_2025_full": ROOT / "data/donors/h1_2025_full/perts.txt",
    "perturbmatch": ROOT / "data/external/perturbmatch/perts.txt",
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat/perts.txt",
    "nadig_hepg2": ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2/perts.txt",
    "k562_essential": ROOT / "data/donors/k562_essential/perts.txt",
    "dld1_lowmoi": ROOT / "experiments/223-silent-zero-audit/work/tables/dld1_lowmoi/perts.txt",
}


def read_list(p: Path) -> set[str]:
    return {ln.strip() for ln in p.read_text().split("\n") if ln.strip()}


def main() -> None:
    banks = {k: read_list(v) for k, v in BANKS.items() if v.exists()}
    out = {}
    for name in PANELS:
        panel = [str(t) for t in pd.read_csv(HERE / "work/panels" / f"{name}.csv")["target_gene"]]
        pset = set(panel)
        rec = {k: len(pset & v) for k, v in banks.items()}
        ncar = [sum(t in banks[d] for d in ("k562_gwps", "orion_hct116", "orion_hek293t")) for t in panel]
        rec["carriers_3_2_1_0"] = [ncar.count(k) for k in (3, 2, 1, 0)]
        rec["h1_union"] = len(pset & (banks["gse295214_h1"] | banks["h1_2025"]))
        rec["pm_support_k562_and_pm"] = len(pset & banks["k562_gwps"] & banks.get("perturbmatch", set()))
        rec["programme_rows_h1_2025_full"] = len(banks["h1_2025_full"] - pset)
        out[name] = rec
    cols = ["k562_gwps", "orion_hct116", "orion_hek293t", "carriers_3_2_1_0", "gse295214_h1", "h1_2025", "h1_union",
            "h1_2025_full", "programme_rows_h1_2025_full", "pm_support_k562_and_pm", "nadig_jurkat", "nadig_hepg2",
            "k562_essential", "dld1_lowmoi"]
    print("panel".ljust(9) + " ".join(c[:14].rjust(14) for c in cols))
    for name, rec in out.items():
        print(name.ljust(9) + " ".join(str(rec.get(c, "-")).rjust(14) for c in cols))
    (HERE / "work/census481.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
