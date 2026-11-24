#!/usr/bin/env python3
"""708 AMENDMENT 5: the Part A known answer, comparing like-labelled stages (science [d80a10], 10-16).

The rehearsal runs dayof438 on 180's relabelled A/B/C bundle, i.e. the D/E/F code path, where every context is unfitted:
claims f = 1 from the bundle's own NTC cells and the EV margin sd 0.01. The reference (work/expected_abcdef.txt) is 438's
stored abcdef run on A/B/C labels, which uses 187's board-fitted claims (ntclevel417: f = 1.00 / 1.00 / 0.9825, 187 identity
ratios 1.0000 / 0.9927 / 1.0000) and sd 0.002. So every LEVELLED quantity differs by design, and the known answer is:
  (i)   every claim-independent cell is identical: the design table's two unlevelled columns (per fold; mean, n+, CI) on every
        row, and the retention line (ret_f per fold);
  (ii)  every claim-dependent rho is the reference's rho times claim_ref / claim_got per context (rho = B / claim, B unchanged),
        to rounding: the design table's 3-decimal "rho unlevelled" column within 0.002 relative, the header's 4-decimal
        "rho unlevelled" lines within 0.0004 relative;
  (iii) the RECOMMENDATION line is identical (tag aside).
Levelled cells (levelled cos / pds_bd / mse, levelled prices, movers, argmax, and each q's RULE CHOICE, which the rule reads
from LEVELLED prices) are reported, never compared; (iii) is the one levelled output required to agree, because it is the
decision (the board's q 1.25).
usage: uv run python experiments/708-def-panel-repricing/kaA708.py <dayof438 printout> [<expected>] [--out JSON]
Prints "PART A KA PASS (Amendment 5)" or "PART A KA FAIL (Amendment 5): <why>"; exit 0 / 1.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def parse(p: Path) -> dict:
    lines = p.read_text().splitlines()
    out = {"claims": None, "ret": None, "rho_hdr": {}, "table": {}, "choices": [], "rec": None}
    for ln in lines:
        m = re.search(r"retention (.*); claims ([0-9.]+)/([0-9.]+)/([0-9.]+)", ln)
        if m:
            out["ret"], out["claims"] = m.group(1), [float(m.group(k)) for k in (2, 3, 4)]
        m = re.match(r"(\S+): rho unlevelled ([0-9.]+)/([0-9.]+)/([0-9.]+); levelled", ln)
        if m:
            out["rho_hdr"][m.group(1)] = [float(m.group(k)) for k in (2, 3, 4)]
        if ln.startswith("| ") and not ln.startswith("|---") and "design (vs" not in ln:
            cells = [c.strip() for c in ln.strip().strip("|").split("|")]
            if len(cells) == 9:                       # the design table (the surfaces have 7 cells)
                out["table"].setdefault(cells[0], []).append(cells)
        m = re.match(r"RULE CHOICE at (q [0-9.]+)[^:]*: (H1 x[0-9.]+, CD4 x[0-9.]+)", ln)
        if m:
            out["choices"].append((m.group(1), m.group(2)))
        if ln.startswith("RECOMMENDATION"):
            out["rec"] = re.sub(r"\([^)]*\)", "(TAG)", ln, count=1)
    return out


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("--out")]
    outp = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--out=")), None)
    got_p = Path(args[0])
    exp_p = Path(args[1]) if len(args) > 1 else HERE / "work" / "expected_abcdef.txt"
    g, e = parse(got_p), parse(exp_p)
    bad, info = [], {}
    if g["claims"] is None or e["claims"] is None:
        bad.append("claims line missing")
    ratio = [ce / cg for ce, cg in zip(e["claims"], g["claims"])] if not bad else [1, 1, 1]
    info["claims_got"], info["claims_ref"], info["claim_ratio_ref_over_got"] = g["claims"], e["claims"], ratio
    if g["ret"] != e["ret"]:
        bad.append(f"retention line differs: {g['ret']!r} vs {e['ret']!r}")
    n_rows, unlev_same, rho_dev = 0, 0, 0.0
    for name, rows in e["table"].items():
        grows = g["table"].get(name, [])
        if len(grows) != len(rows):
            bad.append(f"row {name!r}: {len(grows)} got vs {len(rows)} expected")
            continue
        for gr, er in zip(grows, rows):
            n_rows += 1
            if gr[4] != er[4] or gr[5] != er[5]:
                bad.append(f"row {name!r}: unlevelled cells differ: {gr[4:6]} vs {er[4:6]}")
            else:
                unlev_same += 1
            rg = [float(v) for v in gr[8].split("/")]
            re_ = [float(v) for v in er[8].split("/")]
            for k in range(3):
                dev = abs(rg[k] / (re_[k] * ratio[k]) - 1.0)
                rho_dev = max(rho_dev, dev)
                if dev > 0.002:
                    bad.append(f"row {name!r} context {k + 1}: rho {rg[k]} vs ref {re_[k]} x ratio {ratio[k]:.5f} (dev {dev:.4f})")
    hdr_dev = 0.0
    for key, rv in e["rho_hdr"].items():
        gv = g["rho_hdr"].get(key)
        if gv is None:
            bad.append(f"header rho line {key!r} missing")
            continue
        for k in range(3):
            dev = abs(gv[k] / (rv[k] * ratio[k]) - 1.0)
            hdr_dev = max(hdr_dev, dev)
            if dev > 0.0004:
                bad.append(f"header {key!r} context {k + 1}: rho {gv[k]} vs ref {rv[k]} x {ratio[k]:.5f} (dev {dev:.5f})")
    info["rule_choice_ref"] = e["choices"]                 # levelled: reported, not compared (Amendment 5)
    if not g["rec"] or g["rec"] != e["rec"]:
        bad.append(f"RECOMMENDATION differs: {g['rec']!r} vs {e['rec']!r}")
    info.update(rows_compared=n_rows, rows_unlevelled_identical=unlev_same, max_rho_table_dev=rho_dev,
                max_rho_header_dev=hdr_dev, rule_choices=g["choices"], recommendation=g["rec"],
                got=str(got_p), ref=str(exp_p), verdict="PASS" if not bad else "FAIL", problems=bad[:20])
    if outp:
        Path(outp).write_text(json.dumps(info, indent=1))
    print(f"claims got {g['claims']} ref {e['claims']}: ratio ref/got {[round(r, 5) for r in ratio]}")
    print(f"{unlev_same}/{n_rows} design rows: unlevelled cells identical; max rho deviation from the claim ratio: table "
          f"{rho_dev:.5f}, header {hdr_dev:.6f}; RECOMMENDATION identical: {g['rec'] == e['rec']}")
    print(f"reported (levelled, not compared): RULE CHOICE got {g['choices']} / ref {e['choices']}")
    if bad:
        print("PART A KA FAIL (Amendment 5): " + "; ".join(bad[:5]))
        return 1
    print("PART A KA PASS (Amendment 5)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
