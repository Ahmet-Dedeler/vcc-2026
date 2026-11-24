"""Read a D/E/F entry pair against docs/FINAL-PHASE-FEEDBACK-PLAN.md section 2.

  uv run python experiments/531-cloud-final-phase-tooling/readdef.py ROW ENTRY BASE [ALT]

ROW is a row id from that table: 1, 2, 3, 4, 5, 6, 7, 8-9, 10-12, 13-14.
8, 9, 10, 11, 12, 13 and 14 are accepted as aliases of the grouped rows.
ENTRY, BASE and the optional ALT are captured-entry JSON files (the wrapper
written by scripts/capture_board_entry.py, or a bare entry object) or entry
ids resolved under data/board-entries/<id>.json.

Bulk rows (2, 3, 4, and the 10-12 combination) read exactly: any strictly
positive difference adopts. There is no noise floor. Direction members
reproduce to 16 significant figures (plan section 1 rule 4).

DE rows that are not prior-only (8-9) need an effect > 0.009 overall on top
of the row's own bin, unless a second seed is being read. Prior-only rows
(6, and 7 which is read468's rule) do not apply that floor: they share
emission draws with the base (plan section 1 rule 4, the 468 case).

Row 5 is bulk and DE. Overall is read exactly, and a w that raises the bulk
and costs the DE channel more than 0.009 overall is not adopted.

"Overall contribution" of a set of members is the sum of their scaled deltas
divided by 6, the same convention as experiments/456-zs-riders/read0928.py.
score_avg is the board's overall and is what "overall better" reads.

Exit 0 when a verdict was printed, including a decision not to adopt.
Exit 2 when a row that reads members was given an entry with members hidden
or absent. Exit 1 on usage errors.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ENTRIES = Path(__file__).resolve().parents[2] / "data" / "board-entries"

RAW = {
    "pds": "pds_cosine",
    "mse": "expr_mse_unbiased_capped_norm",
    "nmae": "de_wilcoxon_lfc_nmae",
    "fid": "de_wilcoxon_direction_fidelity_yield_raw",
    "reach": "de_wilcoxon_direction_reach_raw",
    "jac": "de_wilcoxon_sig_jaccard",
}
SCALED = {m: f"score_{m}" for m in RAW}
NEEDED = ["score_avg", *RAW.values(), *SCALED.values()]

# Plan section 1 rule 4. Strict: an effect of exactly 0.009 does not clear it,
# and a DE cost of exactly 0.009 is not "beyond" 0.009.
DE_NOISE = 0.009
# read0928 / read468: pds and mse "identical" at 1e-6 raw. jac at 0.0010.
# The slack absorbs the ulp from subtracting two nearby floats (0.5 + 1e-6 is not
# exactly 1e-6 above 0.5). It is far below the next bin edge the tests use (1e-7).
PIN = 1e-6
JAC_PIN = 0.0010
PIN_SLACK = 1e-15

ALIASES = {
    "1": "1",
    "2": "2",
    "3": "3",
    "4": "4",
    "5": "5",
    "6": "6",
    "7": "7",
    "8": "8-9",
    "9": "8-9",
    "8-9": "8-9",
    "10": "10-12",
    "11": "10-12",
    "12": "10-12",
    "10-12": "10-12",
    "13": "13-14",
    "14": "13-14",
    "13-14": "13-14",
}


def _find(x):
    if isinstance(x, dict):
        if "score_avg" in x or "entry_id" in x:
            return x
        for v in x.values():
            y = _find(v)
            if y:
                return y
    if isinstance(x, list):
        for v in x:
            y = _find(v)
            if y:
                return y
    return None


def load_entry(spec: str) -> dict:
    """A path to a capture, or an entry id under data/board-entries/."""
    p = Path(spec)
    if not p.is_file():
        cand = ENTRIES / f"{spec}.json"
        if cand.is_file():
            p = cand
        else:
            raise SystemExit(f"no capture at {spec!r} and no data/board-entries/{spec}.json")
    found = _find(json.loads(p.read_text()))
    if not isinstance(found, dict):
        raise SystemExit(f"{p} has no entry object")
    out = dict(found)
    out["_source"] = str(p)
    out.setdefault("entry_id", p.stem)
    return out


def _num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def members_hidden(e: dict) -> list[str]:
    return [k for k in NEEDED if not _num(e.get(k))]


def delta(arm: dict, base: dict) -> dict:
    raw = {m: arm[RAW[m]] - base[RAW[m]] for m in RAW}
    scaled = {m: arm[SCALED[m]] - base[SCALED[m]] for m in RAW}
    bulk = (scaled["pds"] + scaled["mse"]) / 6.0
    de = (scaled["fid"] + scaled["reach"] + scaled["jac"] + scaled["nmae"]) / 6.0
    fnr = (scaled["fid"] + scaled["nmae"] + scaled["reach"]) / 6.0
    return {
        "raw": raw,
        "scaled": scaled,
        "overall": arm["score_avg"] - base["score_avg"],
        "bulk": bulk,
        "de": de,
        "fnr": fnr,
        "arm_overall": arm["score_avg"],
        "base_overall": base["score_avg"],
    }


def _fmt_block(tag: str, arm: dict, base: dict, d: dict) -> list[str]:
    lines = [
        f"== {tag} {arm.get('entry_id')} vs {base.get('entry_id')}: "
        f"overall {d['arm_overall']:.7f} vs {d['base_overall']:.7f}  d {d['overall']:+.7f}"
    ]
    for m in RAW:
        lines.append(
            f"  {m:5s} raw d {d['raw'][m]:+.6f}  scaled d {d['scaled'][m]:+.6f}  "
            f"-> overall {d['scaled'][m] / 6:+.6f}"
        )
    lines.append(
        f"  bulk (pds+mse)/6 {d['bulk']:+.12g}; DE (fid+reach+jac+nmae)/6 {d['de']:+.12g}; "
        f"fid+nmae+reach overall {d['fnr']:+.12g}"
    )
    return lines


def _hidden(row: str, which: str) -> dict:
    return {
        "row": row,
        "code": "hidden",
        "adopt": None,
        "which": None,
        "bin": f"MEMBERS HIDDEN on {which}",
        "adoption": (
            "none. Members are hidden or absent, so this row is not read. "
            "The dark plan (DEF-RUNBOOK section 6) stands for this pair."
        ),
        "lines": [],
    }


def _pack(row, code, adopt, which, bin_s, adoption, lines) -> dict:
    return {
        "row": row,
        "code": code,
        "adopt": adopt,
        "which": which,
        "bin": bin_s,
        "adoption": adoption,
        "lines": lines,
    }


def _better(value: float) -> bool:
    """Bulk rows read exactly. Strictly positive adopts; a tie does not."""
    return value > 0.0


def _row2(arm_d, alt_d) -> tuple[str, bool | None, str | None, str, str]:
    """Returns code, adopt, which, bin, adoption."""
    a_win = _better(arm_d["bulk"])
    if alt_d is None:
        if a_win:
            return (
                "adopt_arm", True, "arm",
                "this side beats the base on mse+pds overall",
                "ADOPT this side. One side only: the other scale was not supplied, so this is not a both-lose read.",
            )
        return (
            "base_stands", False, None,
            "this side does not beat the base on mse+pds overall",
            "do not adopt this side. R0's ρ* stands for this side.",
        )
    b_win = _better(alt_d["bulk"])
    if a_win and b_win:
        if arm_d["bulk"] > alt_d["bulk"]:
            return (
                "adopt_arm", True, "arm",
                "the better side beats the base on mse+pds overall",
                "ADOPT the first arm (the better mse+pds side).",
            )
        if alt_d["bulk"] > arm_d["bulk"]:
            return (
                "adopt_alt", True, "alt",
                "the better side beats the base on mse+pds overall",
                "ADOPT the second arm (the better mse+pds side).",
            )
        return (
            "tie", False, None,
            "both sides beat the base by the same mse+pds overall",
            "do not adopt either side. A tie is not a choice of ρ*.",
        )
    if a_win:
        return (
            "adopt_arm", True, "arm",
            "the better side beats the base on mse+pds overall; the other side loses",
            "ADOPT the first arm.",
        )
    if b_win:
        return (
            "adopt_alt", True, "alt",
            "the better side beats the base on mse+pds overall; the other side loses",
            "ADOPT the second arm.",
        )
    return (
        "base_stands", False, None,
        "both sides lose on mse+pds overall",
        "do not adopt either side. R0's ρ* stands.",
    )


def _row3(arm_d, alt_d) -> tuple[str, bool | None, str | None, str, str]:
    """ENTRY is the lower q (1.5), ALT the higher q (2.0), each vs the day-2 best."""
    a_win = _better(arm_d["overall"])
    if alt_d is None:
        if a_win:
            return (
                "adopt_arm", True, "arm",
                "overall better than the day-2 best",
                "ADOPT this q. One q only: the curve-turn read needs both.",
            )
        return (
            "base_stands", False, None,
            "not overall better than the day-2 best",
            "do not adopt this q. The day-2 best stands. Stop the q curve here.",
        )
    b_win = _better(alt_d["overall"])
    if not a_win and not b_win:
        return (
            "base_stands", False, None,
            "neither q is overall better than the day-2 best",
            "do not adopt either q. The day-2 best stands. Stop the q curve.",
        )
    if a_win and b_win:
        if alt_d["overall"] > arm_d["overall"]:
            return (
                "adopt_alt", True, "alt",
                "both q beat the day-2 best and the higher q is better: the curve has not turned",
                "ADOPT the higher q. No further q is in the plan.",
            )
        if arm_d["overall"] > alt_d["overall"]:
            return (
                "adopt_arm", True, "arm",
                "the higher q is worse than the lower q: the q curve has turned",
                "ADOPT the lower q. Stop the q curve. Do not adopt the higher q.",
            )
        return (
            "tie", False, None,
            "both q beat the day-2 best by the same overall",
            "do not adopt either q. A tie is not a step up the curve.",
        )
    if a_win:
        return (
            "adopt_arm", True, "arm",
            "the lower q is overall better than the day-2 best and the higher q is not: the q curve has turned",
            "ADOPT the lower q. Stop the q curve. Do not adopt the higher q.",
        )
    return (
        "adopt_alt", True, "alt",
        "the lower q loses and the higher q is overall better than the day-2 best",
        "ADOPT the higher q. It beats the day-2 best. The curve is not monotone.",
    )


def _row4_one(d) -> tuple[bool, str]:
    if _better(d["overall"]):
        return True, "overall better"
    return False, "not overall better"


def _row5(d) -> tuple[str, bool, str, str]:
    raises_bulk = d["bulk"] > 0.0
    de_cost = -d["de"]
    veto = raises_bulk and d["de"] < -DE_NOISE
    better = _better(d["overall"])
    if veto:
        return (
            "veto_de", False,
            f"bulk up ({d['bulk']:+.12g}) and DE cost {de_cost:.12g} is beyond 0.009",
            "NOT adopted. A w that raises the bulk and costs DE beyond 0.009 is not adopted.",
        )
    if better:
        return (
            "adopt", True,
            "overall better, and the DE cost is not beyond 0.009 with the bulk up",
            "ADOPT this w.",
        )
    return (
        "not_better", False,
        "not overall better",
        "do not adopt this w. The base stands.",
    )


def _row6(d) -> tuple[str, bool, str, str]:
    # Prior-only: the 0.009 floor does not apply, even when overall itself is under it.
    if d["fnr"] >= 0.0015:
        return (
            "adopt", True,
            f"fid+nmae+reach overall {d['fnr']:+.12g} >= +0.0015",
            "ADOPT the programme source. Prior-only: the 0.009 noise rule does not apply.",
        )
    return (
        "below_bin", False,
        f"fid+nmae+reach overall {d['fnr']:+.12g} < +0.0015",
        "do not adopt the programme source. Prior-only: the 0.009 noise rule does not apply and does not rescue a miss.",
    )


def _row7(d) -> tuple[str, bool, str, str]:
    """read468's rule, on this pair. Prior-only, so 0.009 does not apply.

    Falsifiers match read468.py: |raw pds| or |raw mse| > 1e-6, or |raw jac| > 0.0010.
    Bins match read468.py's PRIMARY and the 2026-09-27 addendum: only HARM closes
    the gate; NULL keeps le1.
    """
    bad = [m for m in ("pds", "mse") if _exceeds(d["raw"][m], PIN)]
    jac_bad = _exceeds(d["raw"]["jac"], JAC_PIN)
    if d["raw"]["reach"] >= 0.003 and d["raw"]["fid"] <= 0.0:
        mechanism = "the H1 mechanism transfers (reach up, fid down)"
    elif d["raw"]["reach"] <= 0.0 and d["raw"]["fid"] < 0.0:
        mechanism = "the gate only removes signal (reach and fid both down)"
    else:
        mechanism = f"mixed (d reach {d['raw']['reach']:+.6f}, d fid {d['raw']['fid']:+.6f})"
    if bad or jac_bad:
        fired = list(bad) + (["jac"] if jac_bad else [])
        return (
            "bug", False,
            f"FALSIFIER FIRED {fired}. {mechanism}",
            "do not adopt. A falsifier fired: diagnose before reading the gate.",
        )
    ov = d["overall"]
    if ov >= 0.0015:
        primary = "TRANSFERS (>= +0.0015)"
        code, adopt = "adopt_transfers", True
        adoption = "ADOPT le1."
    elif ov >= 0.0005:
        primary = "WEAK POSITIVE (+0.0005 .. +0.0015): keep the gate, flagged unconfirmed"
        code, adopt = "adopt_weak", True
        adoption = "ADOPT le1, flagged unconfirmed."
    elif ov <= -0.0005:
        primary = "HARM (<= -0.0005): close the gate"
        code, adopt = "harm", False
        adoption = "do not adopt. HARM closes the gate."
    else:
        primary = "NULL (|d| < 0.0005): keep le1 (only HARM closes the gate)"
        code, adopt = "adopt_null", True
        adoption = "ADOPT le1. A NULL read keeps the gate."
    return code, adopt, f"{primary}. {mechanism}", adoption + " Prior-only: the 0.009 noise rule does not apply."


def _exceeds(value: float, limit: float) -> bool:
    return abs(value) > limit + PIN_SLACK


def _pinned_bulk(d) -> list[str]:
    return [m for m in ("pds", "mse") if _exceeds(d["raw"][m], PIN)]


def _row89(d) -> tuple[str, bool, str, str]:
    """Called-gene scale. Not prior-only: nmae bin AND overall > 0.009."""
    moved = _pinned_bulk(d)
    nmae_ok = d["raw"]["nmae"] <= -0.002
    noise_ok = d["overall"] > DE_NOISE
    if moved:
        return (
            "not_pinned", False,
            f"pds/mse not pinned {moved} (raw |d| > 1e-6); nmae raw {d['raw']['nmae']:+.6f}",
            "do not adopt. The plan requires pinned pds/mse.",
        )
    if not nmae_ok:
        return (
            "nmae_miss", False,
            f"nmae raw {d['raw']['nmae']:+.6f} is not <= -0.002",
            "do not adopt.",
        )
    if not noise_ok:
        return (
            "below_noise", False,
            f"nmae raw {d['raw']['nmae']:+.6f} <= -0.002 at pinned pds/mse, "
            f"but overall {d['overall']:+.12g} is not > 0.009",
            "do not adopt from this single entry. A DE arm needs an effect > ~0.009 overall, or two seeds.",
        )
    return (
        "adopt", True,
        f"nmae raw {d['raw']['nmae']:+.6f} <= -0.002 at pinned pds/mse, overall {d['overall']:+.12g} > 0.009",
        "ADOPT this scale.",
    )


def _row89_pair(arm_d, alt_d) -> tuple[str, bool | None, str | None, str, str]:
    ca, aa, ba, da = _row89(arm_d)
    if alt_d is None:
        return ca, aa, ("arm" if aa else None), ba, da
    cb, ab, bb, db = _row89(alt_d)
    if aa and ab:
        # One scale ships. Prefer the more negative raw nmae.
        if arm_d["raw"]["nmae"] < alt_d["raw"]["nmae"]:
            which, extra = "arm", "ADOPT the first scale (the more negative nmae). Do not ship both."
        elif alt_d["raw"]["nmae"] < arm_d["raw"]["nmae"]:
            which, extra = "alt", "ADOPT the second scale (the more negative nmae). Do not ship both."
        else:
            return (
                "tie", False, None,
                f"both scales clear the bin at the same nmae. first: {ba}. second: {bb}",
                "do not adopt either. A tie of scales is not a choice.",
            )
        return "adopt_" + which, True, which, f"both clear. first: {ba}. second: {bb}", extra
    if aa:
        return "adopt_arm", True, "arm", f"first clears ({ba}). second does not ({bb}).", "ADOPT the first scale."
    if ab:
        return "adopt_alt", True, "alt", f"second clears ({bb}). first does not ({ba}).", "ADOPT the second scale."
    return "none", False, None, f"neither clears. first {ca}: {ba}. second {cb}: {bb}.", "do not adopt either scale."


def _row10(d) -> tuple[str, bool, str, str]:
    if _better(d["overall"]):
        return (
            "adopt", True,
            "the combination is scored and overall better than the base",
            "ADOPT the combination as the new best. It was scored, so rule 2 allows it to be the last entry.",
        )
    return (
        "scored_not_better", False,
        "the combination is scored and is not overall better than the base",
        "KEEP the base. The combination was scored (rule 2 is met) and it is not the best-supported package.",
    )


def verdict(row: str, arm: dict, base: dict, alt: dict | None = None) -> dict:
    key = ALIASES.get(str(row).strip())
    if key is None:
        raise SystemExit(f"unknown row {row!r}; known: {', '.join(sorted(set(ALIASES)))}")
    lines: list[str] = []
    if key in ("1", "13-14"):
        # These rows do not read a bin. Still show deltas when members exist.
        if not members_hidden(arm) and not members_hidden(base):
            lines += _fmt_block("entry", arm, base, delta(arm, base))
        if key == "1":
            return _pack(
                key, "baseline", None, None,
                "day 1 baseline on D/E/F",
                "no adoption decision. R0 is the recipe, not an arm.",
                lines,
            )
        return _pack(
            key, "reserve", None, None,
            "reserve: re-upload the best scored package as the LAST entry",
            "no new arm. The last upload is the best-supported package already scored, unchanged.",
            lines,
        )

    missing = []
    ha, hb = members_hidden(arm), members_hidden(base)
    if ha:
        missing.append(f"entry {arm.get('entry_id')} ({', '.join(ha)})")
    if hb:
        missing.append(f"base {base.get('entry_id')} ({', '.join(hb)})")
    if alt is not None:
        hc = members_hidden(alt)
        if hc:
            missing.append(f"alt {alt.get('entry_id')} ({', '.join(hc)})")
    if missing:
        out = _hidden(key, "; ".join(missing))
        out["lines"] = lines
        return out

    arm_d = delta(arm, base)
    alt_d = delta(alt, base) if alt is not None else None
    lines += _fmt_block("arm", arm, base, arm_d)
    if alt_d is not None:
        lines += _fmt_block("alt", alt, base, alt_d)

    if key == "2":
        code, adopt, which, bin_s, adoption = _row2(arm_d, alt_d)
    elif key == "3":
        code, adopt, which, bin_s, adoption = _row3(arm_d, alt_d)
    elif key == "4":
        a_ok, a_bin = _row4_one(arm_d)
        if alt_d is None:
            code = "adopt_arm" if a_ok else "base_stands"
            which = "arm" if a_ok else None
            bin_s = a_bin
            adoption = "ADOPT this arm." if a_ok else "do not adopt. The base stands."
            adopt = a_ok
        else:
            b_ok, b_bin = _row4_one(alt_d)
            bin_s = f"first: {a_bin}. second: {b_bin}. The two arms are different switches; each is read on its own."
            if a_ok and b_ok:
                code, adopt, which = "adopt_both", True, "both"
                adoption = "ADOPT both switches. They are not alternatives."
            elif a_ok:
                code, adopt, which = "adopt_arm", True, "arm"
                adoption = "ADOPT the first arm only."
            elif b_ok:
                code, adopt, which = "adopt_alt", True, "alt"
                adoption = "ADOPT the second arm only."
            else:
                code, adopt, which = "base_stands", False, None
                adoption = "do not adopt either arm. The base stands."
    elif key == "5":
        if alt_d is not None:
            # w = 0 and w = 1 are alternatives. Each must clear the row-5 rule; ship the better overall among those that do.
            c1, a1, b1, d1 = _row5(arm_d)
            c2, a2, b2, d2 = _row5(alt_d)
            bin_s = f"first ({c1}): {b1}. second ({c2}): {b2}."
            if a1 and a2:
                if arm_d["overall"] > alt_d["overall"]:
                    code, adopt, which = "adopt_arm", True, "arm"
                    adoption = "ADOPT the first w (the better overall among the two that clear)."
                elif alt_d["overall"] > arm_d["overall"]:
                    code, adopt, which = "adopt_alt", True, "alt"
                    adoption = "ADOPT the second w (the better overall among the two that clear)."
                else:
                    code, adopt, which = "tie", False, None
                    adoption = "do not adopt either w. Both clear and the overalls tie."
            elif a1:
                code, adopt, which = "adopt_arm", True, "arm"
                adoption = "ADOPT the first w. " + d2
            elif a2:
                code, adopt, which = "adopt_alt", True, "alt"
                adoption = "ADOPT the second w. " + d1
            else:
                code, adopt, which = "none", False, None
                adoption = "do not adopt either w. " + d1 + " " + d2
        else:
            code, adopt, bin_s, adoption = _row5(arm_d)
            which = "arm" if adopt else None
    elif key == "6":
        code, adopt, bin_s, adoption = _row6(arm_d)
        which = "arm" if adopt else None
    elif key == "7":
        code, adopt, bin_s, adoption = _row7(arm_d)
        which = "arm" if adopt else None
    elif key == "8-9":
        code, adopt, which, bin_s, adoption = _row89_pair(arm_d, alt_d)
    elif key == "10-12":
        code, adopt, bin_s, adoption = _row10(arm_d)
        which = "arm" if adopt else None
        adoption += " The ρ* re-check of this combination is a separate row-2 read against it."
    else:
        raise SystemExit(f"unhandled row {key}")
    out = _pack(key, code, adopt, which, bin_s, adoption, lines)
    return out


def render(v: dict) -> str:
    head = [f"row {v['row']}  code {v['code']}"]
    body = list(v["lines"])
    body.append(f"BIN: {v['bin']}")
    body.append(f"ADOPTION: {v['adoption']}")
    return "\n".join(head + body)


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("-h", "--help"):
        print(__doc__.strip())
        return 0 if args else 1
    if args[0] == "--row":
        args = args[1:]
    if len(args) not in (3, 4):
        print(__doc__.strip(), file=sys.stderr)
        return 1
    row, entry_s, base_s = args[0], args[1], args[2]
    alt_s = args[3] if len(args) == 4 else None
    arm = load_entry(entry_s)
    base = load_entry(base_s)
    alt = load_entry(alt_s) if alt_s else None
    v = verdict(row, arm, base, alt)
    print(render(v))
    return 2 if v["code"] == "hidden" else 0


if __name__ == "__main__":
    sys.exit(main())
