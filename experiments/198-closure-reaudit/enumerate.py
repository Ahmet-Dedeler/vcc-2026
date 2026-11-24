#!/usr/bin/env python3
"""198 — enumerate every experiment's verdict, the endpoint it rests on, and the surface.

Read-only sweep of experiments/*/RESULTS.md (+ *BOARD-RESULT.md).
Writes work/enumeration.json and a markdown table.
"""
from __future__ import annotations
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXP = ROOT / "experiments"
OUT = Path(__file__).resolve().parent / "work"
OUT.mkdir(exist_ok=True)

# ---- lexicons -------------------------------------------------------------

VERDICT_PATTERNS = [
    ("NO-GO", re.compile(r"\bno[\s\-_]?go\b", re.I)),
    ("GO", re.compile(r"(?<![a-z\-])\bGO\b(?!way|od|es|ing)")),
    ("NEGATIVE", re.compile(r"\bnegative result\b|\brefuted\b|\bfalsified\b|\bdead end\b", re.I)),
    ("CLOSED", re.compile(r"\bclosed\b|\bdo not reopen\b|\bterminal\b", re.I)),
    ("SHIP", re.compile(r"\bship(ped)?\b|\bsubmitted\b", re.I)),
    ("VETO", re.compile(r"\bveto(ed)?\b", re.I)),
]

# Endpoint = which scored member (or proxy) the verdict was priced on
ENDPOINTS = {
    "pds": re.compile(r"\bpds\b|pds_cosine|macro PDS", re.I),
    "mse": re.compile(r"\bmse\b|expr_mse", re.I),
    "nmae": re.compile(r"\bnmae\b|lfc_nmae", re.I),
    "fid": re.compile(r"\bfid\b|direction_fidelity", re.I),
    "reach": re.compile(r"\breach\b|direction_reach", re.I),
    "jac": re.compile(r"\bjac\b|jaccard", re.I),
    "cosine": re.compile(r"\bcosine\b|\bcos ?θ\b|\bcos theta\b", re.I),
    "sign_acc": re.compile(r"sign accuracy|sign purity|sign precision|signed credit", re.I),
    "n_real": re.compile(r"\bn_real\b", re.I),
    "overall": re.compile(r"\boverall\b", re.I),
}

SURFACES = {
    "h1-test": re.compile(r"h1[\s\-_]?test", re.I),
    "orion-hct116": re.compile(r"orion[\s\-_]?hct116|hct116[\s\-_]?thinned|hct-116", re.I),
    "orion-hek293t": re.compile(r"orion[\s\-_]?hek293t|hek293t[\s\-_]?thinned|hek-293", re.I),
    "board": re.compile(r"\bboard\b|leaderboard|live score", re.I),
    "lodo": re.compile(r"\blodo\b|leave-one-donor", re.I),
    "plolo": re.compile(r"\bplolo\b", re.I),
    "debench": re.compile(r"debench|de-bench|DEBench", re.I),
    "k562-gwps": re.compile(r"k562[\s\-_]?gwps", re.I),
    "h1-2025": re.compile(r"h1[\s\-_]?2025|2025 H1|adata_(Training|Validation|Test)", re.I),
    "somatic": re.compile(r"\bsomatic\b", re.I),
}

# Invalidating-finding classes from the mandate
CLASSES = {
    "C1_magnitude_pool": dict(
        label="nmae / DE-pool built by top-N by magnitude (188: pool is significance-selected)",
        pat=re.compile(
            r"top[\s\-_]?\d+.{0,40}(by |on )?(\|?lfc\|?|magnitude|effect)|"
            r"(magnitude|\|lfc\||abs.{0,4}lfc)[\s\-_]?selected|"
            r"top[\s\-_]?n .{0,20}magnitude", re.I),
        needs=("nmae",),
    ),
    "C2_pds_h1test_detable": dict(
        label="pds verdict sized from an h1-test DE table rather than PLOLO (188 Phase 3)",
        pat=re.compile(r"Surface\s*\(|debench.*pds|pds.*debench|\.pds\(\)", re.I),
        needs=("pds",),
    ),
    "C3_orion_absolute": dict(
        label="absolute pds/cosine/fid/reach/jac/n_real read off an Orion harness",
        pat=re.compile(r"orion|hct116|hek293t", re.I),
        needs=("pds", "fid", "reach", "jac", "cosine", "n_real"),
    ),
    "C4_uniform_nreal": dict(
        label="assumed uniform / single-valued n_real, or used Surface(n_real_target=N)",
        pat=re.compile(r"n_real_target|uniform n_real|fixed n_real|n_real\s*=\s*\d+|"
                       r"power[\s\-]match", re.I),
        needs=(),
    ),
    "C5_nreal2670": dict(
        label="downstream of 093's refuted n_real ≈ 2,670",
        pat=re.compile(r"2,?670|2670", re.I),
        needs=(),
    ),
    "C6_reach_emission_h1test": dict(
        label="reach verdict on an EMISSION-LAW arm taken from h1-test (not validated)",
        pat=re.compile(r"emission|emit_counts|resample|dispersion|overdispersion|"
                       r"emitter|library size|depth", re.I),
        needs=("reach",),
    ),
    "C7_rule13_toll": dict(
        label="profile-class arm vetoed by the −0.0146 emission-law reach toll",
        pat=re.compile(r"0\.0146|reach toll|composed slot gate|slot gate", re.I),
        needs=(),
    ),
    "C8_crossline_adaptation": dict(
        label="cross-line adaptation family (closed strong-form by 184)",
        pat=re.compile(r"cross[\s\-_]?line|context adapt|basal[\s\-_]?condition|adaptation", re.I),
        needs=(),
    ),
    "C9_amplitude": dict(
        label="verdict rested on a global-amplitude argument (187 closed calibration)",
        pat=re.compile(r"\balpha\b|\bα\b|amplitude|shrink", re.I),
        needs=("mse", "nmae"),
    ),
    "C10_key_vs_delivery": dict(
        label="reach verdict conflating ordering-KEY quality with DELIVERY (189/190/185)",
        pat=re.compile(r"ordering|rank(ing)? key|confidence key|\|H\||prefix", re.I),
        needs=("reach",),
    ),
}


def sections(text: str) -> list[str]:
    return text.split("\n")


def scan(path: Path) -> dict:
    text = path.read_text(errors="replace")
    lines = sections(text)
    name = path.parent.name

    title = ""
    for ln in lines[:5]:
        if ln.startswith("#"):
            title = ln.lstrip("# ").strip()
            break

    # verdict lines: any line carrying a verdict token
    vlines = []
    for i, ln in enumerate(lines):
        for vname, pat in VERDICT_PATTERNS:
            if pat.search(ln):
                vlines.append({"kind": vname, "line": i + 1, "text": ln.strip()[:400]})
                break

    kinds = {v["kind"] for v in vlines}
    if "NO-GO" in kinds:
        verdict = "NO-GO"
    elif "VETO" in kinds:
        verdict = "VETO"
    elif "NEGATIVE" in kinds:
        verdict = "NEGATIVE"
    elif "CLOSED" in kinds:
        verdict = "CLOSED"
    elif "GO" in kinds or "SHIP" in kinds:
        verdict = "GO/SHIP"
    else:
        verdict = "—"

    # endpoints and surfaces: count mentions across the whole file, weight the
    # verdict-bearing lines 5x since that is what the verdict rests on
    ep = {}
    for k, pat in ENDPOINTS.items():
        n = len(pat.findall(text))
        nv = sum(len(pat.findall(v["text"])) for v in vlines)
        if n:
            ep[k] = n + 4 * nv
    sf = {}
    for k, pat in SURFACES.items():
        n = len(pat.findall(text))
        nv = sum(len(pat.findall(v["text"])) for v in vlines)
        if n:
            sf[k] = n + 4 * nv

    scored = [k for k in ("pds", "mse", "nmae", "fid", "reach", "jac") if k in ep]
    primary_ep = max(ep, key=ep.get) if ep else "—"
    primary_sf = max(sf, key=sf.get) if sf else "—"

    flags = []
    for cid, spec in CLASSES.items():
        if not spec["pat"].search(text):
            continue
        if spec["needs"] and not any(n in ep for n in spec["needs"]):
            continue
        hits = spec["pat"].findall(text)
        flags.append({"class": cid, "n_hits": len(hits)})

    return dict(
        exp=name, title=title, path=str(path.relative_to(ROOT)),
        n_lines=len(lines), verdict=verdict, verdict_kinds=sorted(kinds),
        verdict_lines=vlines[:12],
        endpoints=ep, primary_endpoint=primary_ep, scored_members=scored,
        surfaces=sf, primary_surface=primary_sf,
        flags=flags,
    )


def main():
    rows = []
    for p in sorted(EXP.glob("*/RESULTS.md")):
        rows.append(scan(p))
    # also board results
    boards = sorted(EXP.glob("*/*BOARD-RESULT.md"))

    (OUT / "enumeration.json").write_text(json.dumps(
        {"rows": rows, "board_results": [str(b.relative_to(ROOT)) for b in boards]},
        indent=1))

    from collections import Counter
    print(f"{len(rows)} RESULTS.md scanned; {len(boards)} BOARD-RESULT.md")
    print("\nverdict counts:", dict(Counter(r["verdict"] for r in rows)))
    print("\nprimary endpoint:", dict(Counter(r["primary_endpoint"] for r in rows).most_common()))
    print("\nprimary surface:", dict(Counter(r["primary_surface"] for r in rows).most_common()))
    fc = Counter()
    for r in rows:
        for f in r["flags"]:
            fc[f["class"]] += 1
    print("\nflag counts (all verdicts):")
    for k, v in fc.most_common():
        print(f"  {k:28s} {v:3d}   {CLASSES[k]['label']}")

    term = [r for r in rows if r["verdict"] in ("NO-GO", "VETO", "NEGATIVE", "CLOSED")]
    print(f"\nterminal-verdict experiments: {len(term)}")
    fc2 = Counter()
    for r in term:
        for f in r["flags"]:
            fc2[f["class"]] += 1
    print("flag counts (terminal only):")
    for k, v in fc2.most_common():
        print(f"  {k:28s} {v:3d}")


if __name__ == "__main__":
    main()
