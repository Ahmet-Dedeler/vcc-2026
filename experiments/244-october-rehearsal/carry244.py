"""244 -- 235's coverage/sizing table (carry235.py, imported and run verbatim) on the rehearsal panel.

Known-answer rows, abort on failure: every `sized` row carry235.json holds for 199's pessimistic
bundle and the shipped-pair frame on A/B/C (7 fields, bit-exact).  Then the same measurement on
238's ef-bundle: every registry source, and the shipped H1 pair as a union, both sizing rules,
stub base = the shipped h1union base median (235's construction: base-independent comparison).

    uv run python experiments/244-october-rehearsal/carry244.py
"""
from __future__ import annotations

import json
import resource
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/235-october-hardening"))
import carry235 as C  # noqa: E402  (235's functions; its main() is never called)

M = C.M
EF = ROOT / "experiments/238-essential-free-panel/work/ef-bundle"
DEF199 = ROOT / "experiments/199-def-runbook/work/def-bundle"
CARRY235 = json.loads((ROOT / "experiments/235-october-hardening/work/carry235.json").read_text())


def strip(row: dict) -> dict:
    return {k: v for k, v in row.items() if k in ("covered_targets", "all_rows_f_per_unit",
                                                    "all_rows_would_solve", "solved_fraction", "carrying_raised")}


def main() -> None:
    checks = []

    def chk(name, got, want):
        checks.append({"check": name, "got": got, "want": want, "verdict": "PASS" if got == want else "FAIL"})

    # KA 1: the shipped pair on A/B/C reproduces the shipped frame (235's own KA, re-run)
    panel, genes = C.panel_and_genes(ROOT / "data/challenge")
    ka = C.size_sources(M.UNION_SOURCES_DEFAULT, panel, genes)
    for key, want in CARRY235["known_answer_h1_pair_abc"].items():
        chk(f"KA 235 h1union abc frame/{key}", ka["frame"][key], want)
    # KA 2: every sized row 235 recorded on 199's pessimistic bundle, re-measured
    panel, genes = C.panel_and_genes(DEF199)
    for key, want in CARRY235["sized"].items():
        tag, name = key.split("/")
        if tag != "def199_pessimistic":
            continue
        row = C.size_sources((name,), panel, genes)
        for f in ("solved_fraction", "all_rows_would_solve", "all_rows_f_per_unit"):
            chk(f"KA 235 sized {key}/{f}", row[f], want[f])
        chk(f"KA 235 sized {key}/carrying_rows", row.get("frame", {}).get("carrying_rows"),
            want.get("frame", {}).get("carrying_rows"))
    # the rehearsal panel
    panel, genes = C.panel_and_genes(EF)
    cov = {n: C.coverage(n, panel) for n in sorted(M.UNION_SOURCE_REGISTRY)}
    sized = {}
    for n in sorted(M.UNION_SOURCE_REGISTRY):
        if cov[n] == 0:
            sized[n] = {"covered": 0, "note": "carries nothing; 235's raise would fire if selected"}
            continue
        row = C.size_sources((n,), panel, genes)
        sized[n] = strip(row) | {"frame": row.get("frame")}
        print(f"[ef238/{n}] covered={cov[n]:3d} carrying_solved={row['solved_fraction']:.6f} "
              f"all_rows_would_solve={row['all_rows_would_solve']:.6f}", flush=True)
        checks.append({"check": f"ef238/{n}: carrying rule > 0", "got": row["solved_fraction"],
                       "verdict": "PASS" if row["solved_fraction"] > 0 else "FAIL"})
        checks.append({"check": f"ef238/{n}: all-rows rule == 0.0", "got": row["all_rows_would_solve"],
                       "verdict": "PASS" if row["all_rows_would_solve"] == 0.0 else "FAIL"})
    pair = C.size_sources(M.UNION_SOURCES_DEFAULT, panel, genes)
    print(f"[ef238/H1 pair {list(M.UNION_SOURCES_DEFAULT)}] covered={pair['covered_targets']} "
          f"carrying_rows={pair['frame']['carrying_rows']} coverage={pair['frame']['coverage_frac']:.4f} "
          f"g_per_unit={pair['frame']['g_per_unit_fraction']!r} solved={pair['solved_fraction']!r} "
          f"all_rows_would_solve={pair['all_rows_would_solve']!r}")
    out = {"ef238_coverage": cov, "ef238_sized": sized,
           "ef238_h1_pair_stub_base": strip(pair) | {"frame": pair["frame"]},
           "stub_base_median_row_norm": C.SHIPPED_BASE_NORM, "target_g": C.G,
           "checks": checks, "peak_rss_gib": round(C.peak_gib(), 3)}
    out["all_pass"] = all(c["verdict"] == "PASS" for c in checks)
    (HERE / "work/carry244.json").write_text(json.dumps(out, indent=1, sort_keys=True, default=str))
    for c in checks:
        if c["verdict"] != "PASS" or c["check"].startswith("KA 235 h1union"):
            print(f"[{c['verdict']}] {c['check']}: got={c['got']!r}" + (f" want={c['want']!r}" if "want" in c else ""))
    n_ka = sum(1 for c in checks if c["check"].startswith("KA"))
    print(f"\nCARRY {'PASS' if out['all_pass'] else 'FAIL'}  {len(checks)} checks ({n_ka} known-answer), "
          f"peak {out['peak_rss_gib']} GiB")
    raise SystemExit(0 if out["all_pass"] else 1)


if __name__ == "__main__":
    main()
