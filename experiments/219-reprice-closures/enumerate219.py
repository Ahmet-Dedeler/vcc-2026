"""219 Phase 1 — enumerate every experiment with a terminal verdict and record, for each,
the machine-readable facts the re-price needs:

  * the verdict token and the headline line it came from
  * every "overall"-unit number the file quotes near a gate word (the MARGIN candidates)
  * which of the six members the file actually prices (a member name within 60 chars of a
    number), and which it says are zero / not priced / not readable
  * which instrument(s) it quotes: PLOLO, h1-test, board, orion, lodo, somatic, k562
  * whether it uses the 096 `claim` route, rule 13's floors (+0.0019 / +0.0089), and
    whether it carries a permutation / shuffle control

This is deliberately mechanical.  Nothing here is a verdict; the classification in
RESULTS.md is made by reading the flagged files.  Method follows 198's `enumerate.py`
(verdict-bearing lines weighted) but the extracted fields are different: 198 recorded the
primary ENDPOINT and SURFACE, this records the PRICED MEMBERS and the MARGIN.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

MEMBERS = ("pds", "mse", "nmae", "fid", "reach", "jac")

VERDICTS = [
    ("VETO", re.compile(r"\bVETO(?:ED|ES)?\b")),
    ("NO-GO", re.compile(r"\bNO[- ]?GO\b", re.I)),
    ("CLOSED", re.compile(r"\bCLOSED?\b|\bSTAYS? CLOSED\b|\bRETIRED\b", re.I)),
    ("NEGATIVE", re.compile(r"\bNEGATIVE\b", re.I)),
    ("GO/SHIP", re.compile(r"\bSHIP IT\b|\bGO\b(?!/)|\bSHIP\b", re.I)),
]

INSTRUMENTS = {
    "plolo": re.compile(r"\bPLOLO\b", re.I),
    "h1-test": re.compile(r"h1[-_ ]?test", re.I),
    "board": re.compile(r"\bboard\b", re.I),
    "orion": re.compile(r"orion[-_]?(hct116|hek293t)", re.I),
    "lodo": re.compile(r"\bLODO\b"),
    "somatic": re.compile(r"\bsomatic\b", re.I),
    "k562-gwps": re.compile(r"k562[-_]gwps", re.I),
    "kolf2": re.compile(r"\bkolf2\b", re.I),
}

ZERO_LANG = re.compile(
    r"carried at (?:~?)zero|not (?:locally )?priceab|unmeasured|not readable|unreadable"
    r"|cannot (?:be )?(?:priced|read)|held at zero|assumed zero|not priced", re.I)
CLAIM_ROUTE = re.compile(r"096|claim route|`claim`|claim\b", re.I)
PERM_CTRL = re.compile(r"permut|shuffl|scrambl", re.I)
RATIO_ARG = re.compile(r"\bratio\b|exchange rate|per unit of", re.I)
FLOOR_19 = re.compile(r"\+?0\.0019\b")
FLOOR_89 = re.compile(r"\+?0\.0089\b")
TOLL = re.compile(r"0\.0146\b")
DIRECTION = re.compile(r"cosine|transfer cos|direction (?:arm|channel|gain)|\bcos\b", re.I)

NUM = re.compile(r"[-+−]?\d+\.\d{3,8}")


def _f(tok: str) -> float:
    return float(tok.replace("−", "-"))


def scan(path: Path) -> dict:
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    head = "\n".join(lines[:1])

    # --- verdict: first match, searching the headline and any "verdict" line first ------
    verdict, vline = None, ""
    pri = [l for l in lines if re.search(r"verdict|GO / NO-GO|GO/NO-GO|^#", l, re.I)]
    for pool in (pri, lines):
        for name, rx in VERDICTS:
            for l in pool:
                if rx.search(l):
                    verdict, vline = name, l.strip()[:220]
                    break
            if verdict:
                break
        if verdict:
            break

    # --- overall-unit numbers appearing on a gate/price line ----------------------------
    margins = []
    for l in lines:
        if not re.search(r"overall|predicted|gate|floor|price|Δovr|TOTAL|net\b", l, re.I):
            continue
        for tok in NUM.findall(l):
            v = _f(tok)
            if abs(v) < 0.2:                      # overall units live here
                margins.append(dict(v=v, line=l.strip()[:200]))

    priced, zeroed = set(), set()
    for l in lines:
        low = l.lower()
        for m in MEMBERS:
            if re.search(rf"[`\b]{m}[`\b]", low) and NUM.search(l):
                priced.add(m)
        if ZERO_LANG.search(l):
            for m in MEMBERS:
                if re.search(rf"[`\b]{m}[`\b]", low):
                    zeroed.add(m)

    instruments = [k for k, rx in INSTRUMENTS.items() if rx.search(text)]

    return dict(
        exp=path.parent.name,
        title=head.lstrip("# ").strip()[:200],
        verdict=verdict,
        verdict_line=vline,
        n_lines=len(lines),
        priced_members=sorted(priced),
        zeroed_members=sorted(zeroed),
        instruments=instruments,
        uses_096_claim=bool(CLAIM_ROUTE.search(text)),
        cites_floor_0019=bool(FLOOR_19.search(text)),
        cites_floor_0089=bool(FLOOR_89.search(text)),
        cites_toll_0146=bool(TOLL.search(text)),
        has_perm_control=bool(PERM_CTRL.search(text)),
        has_ratio_argument=bool(RATIO_ARG.search(text)),
        is_direction_flavoured=bool(DIRECTION.search(text)),
        n_direction_mentions=len(DIRECTION.findall(text)),
        margins=margins[:80],
    )


def main() -> None:
    rows = [scan(p) for p in sorted((ROOT / "experiments").glob("*/RESULTS.md"))]
    rows.sort(key=lambda r: r["exp"])
    (WORK / "enumeration219.json").write_text(json.dumps(rows, indent=1))

    term = [r for r in rows if r["verdict"] in ("NO-GO", "VETO", "CLOSED", "NEGATIVE")]
    print(f"RESULTS.md files            : {len(rows)}")
    from collections import Counter
    print("verdicts                    :", dict(Counter(r["verdict"] for r in rows)))
    print(f"terminal verdicts           : {len(term)}")

    # the mandate's Class-A screen, mechanically: a direction-flavoured terminal verdict
    # that prices mse and/or uses PLOLO, and carries DE members at zero or never prices them
    cand = []
    for r in term:
        de = {"nmae", "fid", "reach", "jac"}
        de_missing = de - set(r["priced_members"])
        direction = r["n_direction_mentions"] >= 3
        mse_or_plolo = ("mse" in r["priced_members"]) or ("plolo" in r["instruments"])
        if direction and mse_or_plolo and (de_missing or set(r["zeroed_members"]) & de):
            cand.append((r["exp"], len(de_missing), sorted(de_missing),
                         r["n_direction_mentions"]))
    cand.sort(key=lambda t: (-t[3], t[0]))
    print(f"\nmechanical Class-A candidates: {len(cand)}")
    for e, n, miss, dm in cand:
        print(f"  {e:38s} de_unpriced={n} {miss}  dir_mentions={dm}")
    (WORK / "candidates219.json").write_text(json.dumps(cand, indent=1))


if __name__ == "__main__":
    main()
