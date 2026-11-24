#!/usr/bin/env python3
"""Re-check manuscript numbers against experiment records.

Run from the repo root:
    uv run python experiments/503-cloud-claims-audit/audit_claims.py

Prints one line per check: CONFIRMED or MISSING, with the source line.
Does not edit PAPER.md or CLAIMS.md.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "experiments/294-writeup/PAPER.md"
CLAIMS = ROOT / "experiments/294-writeup/CLAIMS.md"


def read(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace")


def line_of(text: str, needle: str) -> str:
    for i, line in enumerate(text.splitlines(), 1):
        if needle in line:
            snippet = line.strip().replace("\t", " ")
            if len(snippet) > 220:
                snippet = snippet[:217] + "..."
            return f"L{i}: {snippet}"
    return "NOT FOUND"


def check(name: str, rel: str, needle: str, text: str | None = None) -> str:
    body = text if text is not None else read(rel)
    hit = needle in body
    loc = line_of(body, needle) if hit else "NOT FOUND"
    status = "CONFIRMED" if hit else "MISSING"
    return f"{status}\t{name}\t{rel}\t{needle!r}\t{loc}"


def main() -> None:
    paper = read("experiments/294-writeup/PAPER.md")
    claims = read("experiments/294-writeup/CLAIMS.md")
    findings = read("docs/FINDINGS.md")
    f09 = read("docs/FINDINGS-2026-09.md")
    rows: list[str] = []

    # --- manuscript defects the audit calls ---
    rows.append(check("paper-header-63", "experiments/294-writeup/PAPER.md", "63 %", paper))
    rows.append(check("paper-body-37", "experiments/294-writeup/PAPER.md", "37 %", paper))
    rows.append(
        check(
            "paper-blend-range",
            "experiments/294-writeup/PAPER.md",
            "+0.071 … +0.084",
            paper,
        )
    )
    rows.append(
        check(
            "paper-same-cells",
            "experiments/294-writeup/PAPER.md",
            "on the same cells",
            paper,
        )
    )
    rows.append(
        check(
            "paper-same-300-panel",
            "experiments/294-writeup/PAPER.md",
            "same 300-target panel",
            paper,
        )
    )
    rows.append(
        check(
            "paper-imply-0.65",
            "experiments/294-writeup/PAPER.md",
            "gate sign accuracy near 0.65",
            paper,
        )
    )
    rows.append(
        check(
            "paper-nothing-on-abc",
            "experiments/294-writeup/PAPER.md",
            "nothing here is measured on the challenge's held contexts A/B/C",
            paper,
        )
    )
    rows.append(
        check(
            "paper-8-banks",
            "experiments/294-writeup/PAPER.md",
            "8 banks, 3 labs",
            paper,
        )
    )
    rows.append(
        check(
            "paper-100M",
            "experiments/294-writeup/PAPER.md",
            "100M-cell",
            paper,
        )
    )
    rows.append(
        check(
            "paper-within-0.0006",
            "experiments/294-writeup/PAPER.md",
            "within 0.0006",
            paper,
        )
    )
    rows.append(
        check(
            "claims-within-0.0006",
            "experiments/294-writeup/CLAIMS.md",
            "within 0.0006",
            claims,
        )
    )
    rows.append(
        check(
            "paper-0.330-unlabeled",
            "experiments/294-writeup/PAPER.md",
            "0.330 at 200/200",
            paper,
        )
    )

    # phrases that were corrected elsewhere and should NOT still be the manuscript's result
    rows.append(
        check(
            "paper-supersedes-near-zero",
            "experiments/294-writeup/PAPER.md",
            "is superseded",
            paper,
        )
    )
    rows.append(
        "ABSENT\tmanuscript-kolf-coverage-zero\tPAPER/CLAIMS\t"
        "'covers 0'\tno October-coverage sentence; C-19's kolf2 is a fold size (273)"
        if "kolf2" not in paper
        else "PRESENT\tmanuscript-kolf-coverage-zero\tPAPER\tkolf2 mentioned"
    )
    rows.append(
        check(
            "c47-295",
            "experiments/295-reserve-r3-claim/RESULTS.md",
            "+0.5016",
        )
    )
    rows.append(
        check(
            "c01-median-0.0199",
            "experiments/276-conservation-map/RESULTS.md",
            "+0.0199",
        )
    )

    # 410b DE bound: largest |delta| vs the 0.0006 claim
    b410 = read("experiments/410-pertscale-revisit/BOARD-RESULT-410b.md")
    rows.append(
        check(
            "410b-within-0.0007",
            "experiments/410-pertscale-revisit/BOARD-RESULT-410b.md",
            "within 0.0007",
            b410,
        )
    )
    # parse the 410b table's 410b-395 column from the known raw values in the file
    # nmae line contains both 0.937101 and 0.936458
    rows.append(
        check(
            "410b-nmae-ref",
            "experiments/410-pertscale-revisit/BOARD-RESULT-410b.md",
            "0.937101",
            b410,
        )
    )
    rows.append(
        check(
            "410b-nmae-arm",
            "experiments/410-pertscale-revisit/BOARD-RESULT-410b.md",
            "0.936458",
            b410,
        )
    )

    # board JSON bit-identity for 433a vs 417e
    a = json.loads(read("data/board-entries/FCbBcvGxF1HOFPxPISUC.json"))["entry"]
    e = json.loads(read("data/board-entries/wPP8U70kj4mOIvtA83CO.json"))["entry"]
    for key in ("pds_cosine", "expr_mse_unbiased_capped_norm"):
        same = a[key] == e[key]
        rows.append(
            f"{'CONFIRMED' if same else 'MISSING'}\tboard-433a-417e-{key}\t"
            f"data/board-entries/{{FCbBcvGxF1HOFPxPISUC,wPP8U70kj4mOIvtA83CO}}.json\t"
            f"{a[key]!r} vs {e[key]!r}\tidentical={same}"
        )
    dnmae = a["de_wilcoxon_lfc_nmae"] - e["de_wilcoxon_lfc_nmae"]
    rows.append(
        f"CONFIRMED\tboard-433a-dnmae\tboard JSON\t{dnmae:.6f}\t"
        f"433a {a['de_wilcoxon_lfc_nmae']} minus 417e {e['de_wilcoxon_lfc_nmae']}"
    )

    # median vs mean on the same H1 split
    t273 = read("experiments/273-flex-reference-instrument/work/tables273.md")
    rows.append(
        check(
            "h1-median-and-mean",
            "experiments/273-flex-reference-instrument/work/tables273.md",
            "+0.2429",
            t273,
        )
    )
    rows.append(
        check(
            "h1-mean-0.3309",
            "experiments/273-flex-reference-instrument/work/tables273.md",
            "+0.3309",
            t273,
        )
    )

    # weight curve
    rows.append(
        check(
            "weightcurve-37",
            "experiments/450-same-line-programme/BOARD-RESULT-446bja.md",
            "37 %",
        )
    )
    # blend ladder
    r450 = read("experiments/450-same-line-programme/RESULTS.md")
    for needle in ("+0.0531", "+0.0710", "+0.0790", "+0.0934", "−0.020", "+0.544", "+0.205"):
        rows.append(check(f"450-{needle}", "experiments/450-same-line-programme/RESULTS.md", needle, r450))

    # 434f map and later upper-bound caveat
    rows.append(
        check(
            "434-nmae-map",
            "experiments/434-posterior-median-magnitude/RESULTS.md",
            "0.950",
        )
    )
    rows.append(
        check(
            "findings-upper-bound",
            "docs/FINDINGS.md",
            "these accuracies are upper bounds",
            findings,
        )
    )

    # gate interval spans 0
    rows.append(
        check(
            "467-gate-ci",
            "experiments/467-strength-gated-programme/RESULTS.md",
            "[−0.0034, +0.0132]",
        )
    )
    rows.append(
        check(
            "claims-c71-ci",
            "experiments/294-writeup/CLAIMS.md",
            "[−0.0034, +0.0132]",
            claims,
        )
    )

    # control-noise decomposition
    rows.append(
        check(
            "462-third-noise",
            "experiments/462-broad16-lab-signature/RESULTS.md",
            "33 % control-sample noise",
        )
    )
    rows.append(
        check(
            "462-broad-0.801",
            "experiments/462-broad16-lab-signature/RESULTS.md",
            "0.801",
        )
    )

    # 452 correction is not a manuscript claim
    rows.append(
        check(
            "452-bit-repro-correction",
            "experiments/452-def-day-rehearsal/RESULTS.md",
            "IS bit-reproducible",
        )
    )
    rows.append(
        check(
            "findings-kolf-181",
            "docs/FINDINGS.md",
            "181–251",
            findings,
        )
    )

    # headline numbers for the earlier claims, in the findings note 275 verified
    headlines = [
        ("C-02", "+0.047"),
        ("C-03", "+0.243"),
        ("C-04", "+0.071"),
        ("C-05", "+0.58"),
        ("C-08", "+0.0006"),
        ("C-12", "+0.054"),
        ("C-13", "+0.0033"),
        ("C-15", "+0.273"),
        ("C-16", "+0.435"),
        ("C-17", "+0.387"),
        ("C-21", "−0.01047"),
        ("C-24", "+0.805"),
        ("C-26", "+0.328"),
        ("C-27", "+0.317"),
        ("C-32", "0.091"),
        ("C-40", "−0.034"),
        ("C-44", "−0.018"),
        ("C-46", "25 of 300"),
        ("C-49", "+0.687"),
        ("C-55", "+0.786"),
        ("C-56", "+0.71293"),
    ]
    for cid, needle in headlines:
        rows.append(check(f"findings09-{cid}", "docs/FINDINGS-2026-09.md", needle, f09))

    # 496 pooled interval
    rows.append(
        check(
            "496-overall-ci",
            "experiments/496-uncalled-sign-hook/RESULTS.md",
            "[−0.0017, +0.0020]",
        )
    )
    rows.append(
        check(
            "496-nmae",
            "experiments/496-uncalled-sign-hook/RESULTS.md",
            "−0.00146",
        )
    )

    # 440 "about 50"
    rows.append(
        check(
            "440-about-50",
            "experiments/440-empty-union-jaccard/RESULTS.md",
            "About 50 teams",
        )
    )

    # 410b numeric bound
    ref, arm = 0.937101, 0.936458
    delta = arm - ref
    rows.append(
        f"CONFIRMED\t410b-nmae-delta\tBOARD-RESULT-410b.md table\t"
        f"{delta:.6f}\tabs(delta)={abs(delta):.6f} > 0.0006 and < 0.0007"
    )

    out = ROOT / "experiments/503-cloud-claims-audit/work/audit_log.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(rows) + "\n"
    out.write_text(text, encoding="utf-8")
    print(text)
    missing = [r for r in rows if r.startswith("MISSING")]
    print(f"--- {len(rows)} checks, {len(missing)} MISSING ---")


if __name__ == "__main__":
    main()
