"""223 Phase 1a — mechanical scan: which experiments contain a MEDIAN-ROW-NORM sizing path?

The silent zero has exactly one signature: a scalar computed as
`median(norm(residual_rows))` over EVERY row, used as a divisor or as the numerator of a
dose ratio.  Three known spellings:

  P1  production  `scripts/make_submission_pm_cd4.py::solve_cd4_fraction`
                  -> f_per_unit = median_all / base_median ; RAISES when 0 and TARGET_F set
  P2  harness     `162 plolo.solve_fraction` (and copies)
                  -> returns fraction 0.0 SILENTLY when per_unit == 0
  P3  099-lineage `nV = float(np.median(np.linalg.norm(V, axis=1))) or 1.0 ; V *= nE / nV`
                  -> does NOT zero; it MIS-SIZES by falling back to a divisor of 1.0

Anything else (an explicit blend `fraction`, a rotation, a router, a source-only audit)
cannot express the bug and is reported as such.  The scan is textual and deliberately
over-inclusive; every hit is read by hand afterwards.

    uv run python scan_paths.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

EXPS = [
    "055-kolf-direct", "057-kolf-union-whiten", "059-feng-ipsc",
    "079-kolf-local-centering", "085-kolf-pm-gate", "105-telohaec-residual",
    "107-jurkat-high-moi", "115-gse295214-h1", "116-dld1-lowmoi",
    "119-altair-a549", "120-myllia-a549", "121-a549-cross-source-concordance",
    "126-hcasmc-panel-confirmation", "129-gse295214-transfer",
    "130-telohaec-panel-confirmation", "131-hues8-hepatocyte-metadata",
    "135-telohaec-atlas-confirmation", "137-jurkat-tangent-cone",
    "138-dld1-replicate-sign", "139-jurkat-donor-simplex", "176-cmap-donors",
]
# not in the mandate's 21 but named by 218 section 7.5 as the residual-path closures
BONUS = ["117-pan-h1-direct-residual", "118-pds-null-panel-residual",
         "103-residual-only-production", "099-donors", "062-cd4-post-whiten"]

PATTERNS = {
    "median_row_norm": re.compile(
        r"median\s*\(\s*np\.linalg\.norm|np\.median\s*\(\s*np\.linalg\.norm"),
    "median_or_one_fallback": re.compile(r"np\.median\([^)]*\)\s*\)?\s*or\s+1\.0"),
    "solve_cd4_fraction": re.compile(r"solve_cd4_fraction"),
    "solve_fraction": re.compile(r"\bsolve_fraction\b"),
    "plolo_import": re.compile(r"\bimport plolo\b|from plolo import"),
    "postwhiten_residual": re.compile(r"add_post_whiten_residual|post_whiten|postwhiten"),
    "explicit_fraction_knob": re.compile(
        r"\bfraction\b|\bblend\b|\bFRACTIONS?\b|\bfracs?\b"),
    "makes_submission_builder": re.compile(r"make_submission_pm_cd4"),
}


def scan(exp: str) -> dict:
    d = ROOT / "experiments" / exp
    hits: dict[str, list[str]] = {k: [] for k in PATTERNS}
    n_py = 0
    for p in sorted(d.rglob("*.py")):
        if "__pycache__" in p.parts:
            continue
        n_py += 1
        try:
            txt = p.read_text(errors="replace")
        except Exception:
            continue
        for k, rx in PATTERNS.items():
            if rx.search(txt):
                hits[k].append(str(p.relative_to(d)))
    return dict(exp=exp, n_py=n_py, hits={k: v for k, v in hits.items() if v})


def main() -> None:
    rows = [scan(e) for e in EXPS]
    bonus = [scan(e) for e in BONUS]
    (WORK / "scan_paths.json").write_text(
        json.dumps(dict(mandate=rows, bonus=bonus), indent=2))

    keys = ["median_row_norm", "median_or_one_fallback", "solve_cd4_fraction",
            "solve_fraction", "plolo_import", "postwhiten_residual",
            "explicit_fraction_knob", "makes_submission_builder"]
    hdr = f"{'experiment':38s}{'py':>4s}" + "".join(f"{k[:9]:>10s}" for k in keys)
    for title, group in (("THE 21", rows), ("BONUS (residual-path closures)", bonus)):
        print(f"\n### {title}")
        print(hdr)
        print("-" * len(hdr))
        for r in group:
            marks = "".join(f"{('YES' if k in r['hits'] else '.'):>10s}" for k in keys)
            print(f"{r['exp']:38s}{r['n_py']:>4d}{marks}")
    print(f"\n-> {WORK / 'scan_paths.json'}")


if __name__ == "__main__":
    main()
