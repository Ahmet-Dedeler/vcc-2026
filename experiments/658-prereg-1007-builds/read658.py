#!/usr/bin/env python3
"""658 — the reader for the 10-07 builds' board entries (PREREGISTRATION-ALLSW.md, PREREGISTRATION-524dr2.md), and for 524dr
itself (524's PREREGISTRATION, unchanged; 653 R6's "524's reader", which did not exist).

  uv run python experiments/658-prereg-1007-builds/read658.py allsw ALLSW [--hook ENTRY_496A] [--out J]
  uv run python experiments/658-prereg-1007-builds/read658.py dr ENTRY --recipe RECIPE.json
        [--prev ENTRY_524DR --prev-recipe RECIPE524.json] [--hook ENTRY_496A] [--h1 ENTRY_459DR]
        [--alias ENVKEY=h1_dose] [--extra ENVKEY=ENTRY:REF:bulk|de ...] [--out J]

ENTRY = an entry id under data/board-entries/ or a capture path (531's readdef.load_entry). Exit 0 = read (whatever the bin),
2 = members hidden or a refusal (an unknown recipe switch, a missing read), 3 = a FALSIFIER fired (diagnose before reading).

allsw: falsifiers F1 (bulk pinned to 498p2g, 1e-6 raw) and F2 (the DE switches' own leak rules), then R1 = 637's --validate as
registered (FACT = ALLSW, X = 498p2g, Y = 446bj2z, or a synthetic zs+hook entry when the hook is ON; BASE = 446bj2) and R2 = the
exact interaction I = dDE(ALLSW - 498p2g) - sum dDE(single reads vs 446bj2), binned at 1 and 2 same-law re-roll sd (0.00025).

dr: the secondary read vs 453dr with 524's Sigma rule and bins (|d - Sigma| <= 0.0020 PASS / LAG / SUPER, member falsifiers),
for 524dr and 524dr2 alike; with --prev, 524dr2's primary read vs 524dr by branch (N0 none new, N1 DE-only new, N2 bulk-only new,
N3 both). Every line that misses names the members and the switches that own them (the culprit map)."""
from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/531-cloud-final-phase-tooling"))
sys.path.insert(0, str(ROOT / "experiments/637-final-phase-factorial"))
import readdef as R  # noqa: E402  (read-only)
import readdef_factorial637 as F637  # noqa: E402  (read-only)

E446B, E446BJ2, E446BJ2Z = "UDgjOKL9gI1udAexq6ue", "YvVru3DAF89Jiyr34olR", "be90smZUDcIV7j5MgQ03"
E479P2G, E495, E498P2G, E453DR = "ioNMsGwrG2EvYuoB5Zu7", "9Hxef2F6Clf4G2hIkNPq", "JWSghc0o7MN5D6o8Iadc", "q3VuTImh4FRGKXH2waCp"
BULK, DE = ("pds", "mse"), ("nmae", "fid", "reach", "jac")
REROLL = 0.00025                       # same-law DE re-roll sd, overall units (537 section 2.5, FINDINGS 10-03)
PIN = 1e-6                             # raw; bulk members of a DE-only change, and 637's validation pin
DR_TOL = 0.0020                        # 524's overall tolerance
DR_MEMBER_TOL = {"nmae": 0.002, "fid": 0.003, "reach": 0.012, "jac": 0.001}   # 524's DE falsifier tolerances (raw)
HOOK_LEAK = {"fid": 0.001, "jac": 0.001, "reach": 0.01}                       # 496a's leak check (raw)
TERCILEC_BULK_MAX = 0.0003             # tercilec moves 4 A/B/C targets' tiers (598's KA); more than this is a defect
# which switch owns which members on A/B/C (the culprit map; 524dr2 prereg table)
OWNER = {
    "lodo": "mse (pds small): the LODO key's stack dump. Check meta bulk_load / bulk_load_sha256 and def_day's 'LODO key ON' line",
    "lodo_hi_pct": "mse, pds: the cap-6 percentile of the LODO load. Check bulk_load (hi6 dump) and the hi-pct export",
    "sharpen_pt": "pds up, mse down: per-target key-tercile sharpening. Check meta bulk_sharpen_pt (mode, edges, counts) and the def_day 'per-target sharpening' line",
    "tercilec": "pds (4 A/B/C targets change tier): carried-only cut points. Check bulk_sharpen_pt.mode == 'tercilec' and the tier counts",
    "de_mag": "nmae only (zs; the other five bit-identical on 446bj2z). Check meta de_mag and the DEF_DAY_DE_MAG pass-through",
    "mde_sign": "nmae (fid / jac pinned; reach unpriced): the uncalled-gene hook. Check meta de_mde_sign and the DEF_DAY_MDE_SIGN pass-through",
    "h1_dose": "pds, mse: the stack's H1 dose (459 / dayof438). Check bulk_load and the dump's VCC_H1_UNION_G",
}
ENV2SWITCH = {"DEF_DAY_LODO": "lodo", "DEF_DAY_LODO_HI_PCT": "lodo_hi_pct", "DEF_DAY_SHARPEN_PT": "sharpen_pt",
              "DEF_DAY_DE_MAG": "de_mag", "DEF_DAY_MDE_SIGN": "mde_sign"}
CHANNEL = {"lodo": "bulk", "lodo_hi_pct": "bulk", "sharpen_pt": "bulk", "tercilec": "bulk", "h1_dose": "bulk",
           "de_mag": "de", "mde_sign": "de"}


def load(spec: str) -> dict:
    e = R.load_entry(spec)
    h = R.members_hidden(e)
    if h:
        raise Hidden(f"{spec}: members hidden {h}")
    return e


class Hidden(Exception):
    pass


class Refuse(Exception):
    pass


def synth(base: dict, *pairs: tuple[dict, dict], tag: str = "synthetic") -> dict:
    """base + sum of (arm - ref) deltas, on score_avg, raw and scaled members (all linear inside the clamps)."""
    out = copy.deepcopy(base)
    for arm, ref in pairs:
        out["score_avg"] += arm["score_avg"] - ref["score_avg"]
        for m in R.RAW:
            out[R.RAW[m]] += arm[R.RAW[m]] - ref[R.RAW[m]]
            out[R.SCALED[m]] += arm[R.SCALED[m]] - ref[R.SCALED[m]]
    out["entry_id"] = tag
    return out


def dsum(ds: list[dict]) -> dict:
    """the sum of several R.delta dicts (raw, scaled, overall, bulk, de)."""
    out = {"raw": {m: 0.0 for m in R.RAW}, "scaled": {m: 0.0 for m in R.RAW}, "overall": 0.0, "bulk": 0.0, "de": 0.0}
    for d in ds:
        for m in R.RAW:
            out["raw"][m] += d["raw"][m]
            out["scaled"][m] += d["scaled"][m]
        for k in ("overall", "bulk", "de"):
            out[k] += d[k]
    return out


def fmt(d: dict, tag: str) -> str:
    return (f"{tag}: overall {d['overall']:+.7f} (bulk {d['bulk']:+.7f}, DE {d['de']:+.7f}) | raw "
            + " ".join(f"{m} {d['raw'][m]:+.6f}" for m in R.RAW))


# ---------------------------------------------------------------------------------------------------------------- allsw
ALLSW_LEAK_BAND = "v2"                 # Amendment 1 (10-04, before ALLSW is built): 637's v2 band (0b484be2), centre -0.0002


def read_allsw(allsw: str, hook: str | None) -> dict:
    lines, fals = [], []
    F637.set_leak_band(ALLSW_LEAK_BAND)
    A, X, Y, B = load(allsw), load(E498P2G), load(E446BJ2Z), load(E446BJ2)
    H = load(hook) if hook else None
    dX, dY = R.delta(X, B), R.delta(Y, B)
    singles = [dY] + ([R.delta(H, B)] if H else [])
    dYsum = dsum(singles)
    dA_X = R.delta(A, X)                                   # exact: ALLSW shares 498p2g's bulk and emission draws
    lines.append(f"ALLSW = 498p2g's bulk + DE switches {{zs{', hook' if H else ''}}}; X = 498p2g, Y = zs{' + hook' if H else ''}, BASE = 446bj2")
    lines.append(fmt(dX, "X  498p2g - 446bj2 (measured leak = its DE part)"))
    lines.append(fmt(dYsum, "Y  sum of the DE switches' single reads vs 446bj2"))
    lines.append(fmt(dA_X, "ALLSW - 498p2g (exact read of the DE switches on the new bulk)"))
    pred = dsum([dX, dYsum])
    lines.append(f"prediction (additive): ALLSW = 446bj2 + {pred['overall']:+.7f} = {B['score_avg'] + pred['overall']:.7f}; "
                 f"vs 498p2g {dYsum['overall']:+.7f}")
    # F1: bulk pinned to 498p2g
    for m in BULK:
        if abs(dA_X["raw"][m]) > PIN + R.PIN_SLACK:
            fals.append(f"F1 {m} raw moved {dA_X['raw'][m]:+.3e} vs 498p2g (> 1e-6): ALLSW's bulk is not 498p2g's. Diagnose the build "
                        "(bulk_load sha, bulk_ctx_scale, bulk_sharpen_pt, builder changes since 498p2g) before reading anything")
    # F2: the DE switches' own leak rules
    if not H:
        for m in ("fid", "reach", "jac"):
            if abs(dA_X["raw"][m]) > PIN + R.PIN_SLACK:
                fals.append(f"F2 {m} raw moved {dA_X['raw'][m]:+.3e} vs 498p2g: zs alone is call-neutral (bit-identical on 446bj2z)")
    else:
        dH = singles[1]
        for m, tol in HOOK_LEAK.items():
            if abs(dA_X["raw"][m] - dH["raw"][m] - dY["raw"][m]) > tol:
                fals.append(f"F2 {m} raw {dA_X['raw'][m]:+.6f} vs the single reads' {dH['raw'][m] + dY['raw'][m]:+.6f} (> {tol}): 496a's leak rule")
    out = {"mode": "allsw", "hook": bool(H), "falsifiers": fals, "prediction_vs_446bj2": pred["overall"],
           "prediction_score": B["score_avg"] + pred["overall"], "d_allsw_minus_498p2g": dA_X, "lines": lines}
    if fals:
        out["verdict"] = "FALSIFIED"
        return out
    # R1: 637's validation as registered
    Ysyn = synth(B, (Y, B), *([(H, B)] if H else []), tag="zs+hook (synthetic)") if H else Y
    v = F637.validate(A, X, Ysyn, B)
    out["R1_637_validate"] = {"pass": v["pass"], "lines": v["lines"], "leak_band": F637.LEAK_VERSION}
    lines += [f"R1 (637 --validate, leak band {F637.LEAK_VERSION}: centre {F637.LEAK_CENTRE}, tol {F637.VALIDATE_DE_TOL}): "
              + ("PASS" if v["pass"] else "FAIL")] + ["   " + x for x in v["lines"]]
    # R2: exact interaction
    I = dA_X["de"] - dYsum["de"]
    Im = {m: dA_X["raw"][m] - dYsum["raw"][m] for m in DE}
    b = "ADDITIVE" if abs(I) <= REROLL else ("WEAK INTERACTION" if abs(I) <= 2 * REROLL else "INTERACTION")
    out["R2"] = {"I_overall": I, "I_raw": Im, "bin": b}
    lines.append(f"R2 exact interaction I = {I:+.7f} overall (1 sd = {REROLL}): {b}; raw " + " ".join(f"{m} {Im[m]:+.6f}" for m in DE))
    lo = F637.LEAK_CENTRE - F637.VALIDATE_DE_TOL - dX["de"]      # |leak + I - centre| <= tol
    hi = F637.LEAK_CENTRE + F637.VALIDATE_DE_TOL - dX["de"]
    out["R1_window_I"] = [lo, hi]
    lines.append(f"   R1 passes iff I in [{lo:+.7f}, {hi:+.7f}] (X's measured leak {dX['de']:+.7f} vs 637's centre {F637.LEAK_CENTRE})")
    if dA_X["de"] <= 0:
        lines.append("   the DE switches give nothing on the new bulk (dDE(ALLSW - 498p2g) <= 0): their October values are re-opened "
                     "by the coordinator with a written reason (prereg section 6)")
    out["verdict"] = b
    return out


# ------------------------------------------------------------------------------------------------------------------- dr
def recipe_switches(path: str, extra: dict, alias: dict) -> dict[str, str]:
    """ON switches of a resolved recipe (october_recipe.py command --json): {switch: value}. tercilec is reported as BOTH
    sharpen_pt (its tercile read) and tercilec (no board read; Sigma 0). Unknown exports refuse."""
    rec = json.loads(Path(path).read_text())
    if rec.get("problems"):
        raise Refuse(f"{path}: the recipe has problems {rec['problems']}")
    sw = {}
    for k, v in sorted(rec["env"].items()):
        name = ENV2SWITCH.get(k) or alias.get(k) or (k if k in extra else None)
        if name is None:
            raise Refuse(f"{path}: env {k}={v} has no wired read in 658 (pass --alias / --extra, or amend the prereg first)")
        sw[name] = v
        if name == "sharpen_pt" and str(v).startswith("tercilec"):
            sw["tercilec"] = v
    args = rec.get("args", [])
    if len(args) >= 6 and (args[3], args[4], args[5]) != ("1.75", "0.5", "noprog"):
        raise Refuse(f"{path}: q / w / flags {args[3:6]} differ from 524's (1.75, 0.5, noprog); Sigma has no read for that change")
    return sw


def switch_delta(name: str, entries: dict, extra: dict) -> dict | None:
    """a switch's own A/B/C board read (the entry that decided it minus that entry's reference); None = no read (tercilec)."""
    reads = {"lodo": (E446BJ2, E446B), "de_mag": (E446BJ2Z, E446BJ2), "sharpen_pt": (E479P2G, E446BJ2),
             "lodo_hi_pct": (E495, E446BJ2)}
    if name in reads:
        a, b = reads[name]
        return R.delta(load(a), load(b))
    if name == "tercilec":
        return None
    if name == "mde_sign":
        if not entries.get("hook"):
            raise Refuse("mde_sign is ON: pass --hook <496a entry>")
        return R.delta(load(entries["hook"]), load(E446BJ2))
    if name == "h1_dose":
        if not entries.get("h1"):
            raise Refuse("h1_dose is ON: pass --h1 <459dr entry>")
        return R.delta(load(entries["h1"]), load(E453DR))
    if name in extra:
        a, b, _ch = extra[name]
        return R.delta(load(a), load(b))
    raise Refuse(f"no read wired for switch {name}")


def chan(s: str, extra: dict) -> str | None:
    return CHANNEL.get(s) or (extra[s][2] if s in extra else None)


def culprits(dev: dict, on: list[str], extra: dict) -> list[str]:
    out = []
    for m in R.RAW:
        owners = [s for s in on if (m in BULK and chan(s, extra) == "bulk") or (m in DE and chan(s, extra) == "de")]
        out.append(f"   {m:5s} d - Sigma {dev[m]:+.6f} raw; owned by {owners or ['no ON switch']}")
    return out


def read_dr(entry: str, recipe: str, prev: str | None, prev_recipe: str | None, hook: str | None, h1: str | None,
            extra: dict, alias: dict) -> dict:
    lines, fals = [], []
    entries = {"hook": hook, "h1": h1}
    sw = recipe_switches(recipe, extra, alias)
    on = list(sw)
    E, D453 = load(entry), load(E453DR)
    d = R.delta(E, D453)
    reads = {s: switch_delta(s, entries, extra) for s in on}
    S = dsum([x for x in reads.values() if x is not None])
    lines.append(f"ON switches: {', '.join(f'{s}={sw[s]}' for s in on)}")
    for s, x in reads.items():
        lines.append(fmt(x, f"   read {s}") if x is not None else f"   read {s}: none on the board (Sigma 0; A/B/C moves 4 targets)")
    lines.append(fmt(S, "Sigma (524's rule)"))
    lines.append(fmt(d, f"{entry} - 453dr"))
    dev = {m: d["raw"][m] - S["raw"][m] for m in R.RAW}
    # 524's falsifiers
    for m, tol in DR_MEMBER_TOL.items():
        if abs(dev[m]) > tol:
            fals.append(f"DE member {m} moved {d['raw'][m]:+.6f} raw vs Sigma {S['raw'][m]:+.6f} (> {tol}): integration defect")
    if not any(chan(s, extra) == "bulk" for s in on) and any(abs(d["raw"][m]) > PIN for m in BULK):
        fals.append("pds or mse moved while every ON switch is DE-only")
    gap = d["overall"] - S["overall"]
    sec = "PASS" if abs(gap) <= DR_TOL else ("LAG" if gap < -DR_TOL else "SUPER-ADDITIVE")
    out = {"mode": "dr", "on": on, "sigma": S, "delta_vs_453dr": d, "gap": gap, "falsifiers_secondary": fals}
    lines.append(f"SECONDARY (vs 453dr, 524's bins): d - Sigma = {gap:+.7f} -> {sec}" + (" [falsifiers fired]" if fals else ""))
    if sec != "PASS" or fals:
        lines += culprits(dev, on, extra)
        lines += [f"   {s}: {OWNER[s]}" for s in on if s in OWNER]
    out["secondary"] = sec
    if prev:
        out["primary"] = read_primary(E, sw, prev, prev_recipe, entries, extra, alias, lines)
    out["lines"] = lines
    out["verdict"] = "FALSIFIED" if fals or (out.get("primary", {}).get("falsifiers")) else (out.get("primary", {}).get("bin") or sec)
    return out


def read_primary(E: dict, sw: dict, prev: str, prev_recipe: str | None, entries: dict, extra: dict, alias: dict,
                 lines: list[str]) -> dict:
    if not prev_recipe:
        raise Refuse("--prev needs --prev-recipe (524dr's work/recipe524.json)")
    psw = recipe_switches(prev_recipe, extra, alias)
    dropped = [s for s in psw if s not in sw]
    changed = [s for s in psw if s in sw and psw[s] != sw[s]
               and not (s == "sharpen_pt" and sw[s] == str(psw[s]).replace("tercile:", "tercilec:", 1))]
    if dropped or changed:
        raise Refuse(f"switches dropped {dropped} or changed {changed} since 524dr: the primary read has no rule for that")
    new = [s for s in sw if s not in psw]
    P = load(prev)
    dP = R.delta(E, P)
    newreads = {s: switch_delta(s, entries, extra) for s in new}
    Sn = dsum([x for x in newreads.values() if x is not None])
    bulk_new = [s for s in new if chan(s, extra) == "bulk"]
    de_new = [s for s in new if chan(s, extra) == "de"]
    branch = "N3" if bulk_new and de_new else "N2" if bulk_new else "N1" if de_new else "N0"
    res = {"branch": branch, "new": new, "delta_vs_524dr": dP, "falsifiers": []}
    lines.append(f"PRIMARY (vs 524dr {prev}): new switches {new or 'none'} -> branch {branch}")
    lines.append(fmt(dP, "   524dr2 - 524dr"))
    if branch == "N0":
        moved = [m for m in R.RAW if abs(dP["raw"][m]) > PIN]
        res["falsifiers"] = [f"N0 member {m} moved {dP['raw'][m]:+.3e}: the applied tree changed the October path" for m in moved]
        res["bin"] = "TREE-IDENTICAL" if not moved else "FALSIFIED"
        lines.append("   N0: no new switch. This should have been a LOCAL content check, never a slot (identical content). "
                     + ("Members identical." if not moved else "; ".join(res["falsifiers"])))
        return res
    # bulk side
    if not bulk_new:
        for m in BULK:
            if abs(dP["raw"][m]) > PIN + R.PIN_SLACK:
                res["falsifiers"].append(f"bulk {m} moved {dP['raw'][m]:+.3e} with only DE switches new: tree or determinism defect")
    else:
        h1r = newreads.get("h1_dose")
        if h1r is None:                                 # tercilec alone on the bulk side
            if abs(dP["bulk"]) > TERCILEC_BULK_MAX:
                res["falsifiers"].append(f"bulk moved {dP['bulk']:+.7f} overall (> {TERCILEC_BULK_MAX}) with only tercilec new: "
                                         "it moves 4 A/B/C targets; diagnose the tiers (meta bulk_sharpen_pt) and the re-level")
            res["bulk_bin"] = "CARRIED (tercilec)" if abs(dP["bulk"]) <= TERCILEC_BULK_MAX else "FALSIFIED"
        else:
            half = 0.5 * h1r["bulk"]
            res["bulk_bin"] = ("CONFIRMED (H1 dose)" if dP["bulk"] >= half else
                               "NOT CONFIRMED (H1 dose): October keeps x2" if dP["bulk"] < 0 else "WEAK (H1 dose): keep x2")
        lines.append(f"   bulk: d {dP['bulk']:+.7f} overall vs predicted {Sn['bulk']:+.7f} -> {res.get('bulk_bin', 'pinned')}")
    # DE side
    tol = REROLL if not bulk_new else 2 * REROLL
    if de_new:
        miss = dP["de"] - Sn["de"]
        res["de_bin"] = ("CARRIED" if abs(miss) <= tol else
                         "LAG (the day does not deliver the DE switch)" if dP["de"] <= 0.25 * Sn["de"] else "PARTIAL")
        if "mde_sign" in de_new:
            dH = newreads["mde_sign"]
            for m, t in HOOK_LEAK.items():
                if abs(dP["raw"][m] - dH["raw"][m]) > t:
                    res["falsifiers"].append(f"hook leak: {m} {dP['raw'][m]:+.6f} vs 496a's {dH['raw'][m]:+.6f} (> {t})")
        lines.append(f"   DE: d {dP['de']:+.7f} vs predicted {Sn['de']:+.7f} (tol {tol}) -> {res['de_bin']}")
    else:
        if abs(dP["de"]) > tol:
            res["falsifiers"].append(f"DE moved {dP['de']:+.7f} overall with only bulk switches new (> {tol}, 2 re-roll sd)")
        lines.append(f"   DE: d {dP['de']:+.7f} (a bulk step only re-rolls it; tol {tol})")
    res["bin"] = "FALSIFIED" if res["falsifiers"] else " / ".join(x for x in (res.get("bulk_bin"), res.get("de_bin")) if x)
    for f in res["falsifiers"]:
        lines.append("   FALSIFIER: " + f)
    if res["falsifiers"] or "LAG" in res["bin"] or "PARTIAL" in res["bin"]:
        lines += [f"   {s}: {OWNER[s]}" for s in new if s in OWNER]
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a1 = sub.add_parser("allsw")
    a1.add_argument("allsw")
    a1.add_argument("--hook")
    a1.add_argument("--out")
    a2 = sub.add_parser("dr")
    a2.add_argument("entry")
    a2.add_argument("--recipe", required=True)
    a2.add_argument("--prev")
    a2.add_argument("--prev-recipe")
    a2.add_argument("--hook")
    a2.add_argument("--h1")
    a2.add_argument("--alias", action="append", default=[])
    a2.add_argument("--extra", action="append", default=[])
    a2.add_argument("--out")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "allsw":
            r = read_allsw(a.allsw, a.hook)
        else:
            extra, alias = {}, {}
            for x in a.extra:
                k, v = x.split("=", 1)
                ea, eb, ch = v.split(":", 2)
                if ch not in ("bulk", "de"):
                    raise Refuse(f"--extra {x}: channel must be bulk or de")
                extra[k] = (ea, eb, ch)
            for x in a.alias:
                k, v = x.split("=", 1)
                if v not in CHANNEL:
                    raise Refuse(f"--alias {x}: {v} is not a known switch")
                alias[k] = v
            r = read_dr(a.entry, a.recipe, a.prev, a.prev_recipe, a.hook, a.h1, extra, alias)
    except Hidden as e:
        print(f"MEMBERS HIDDEN: {e}")
        return 2
    except Refuse as e:
        print(f"REFUSED: {e}")
        return 2
    print("\n".join(r["lines"]))
    for f in r.get("falsifiers", []) + r.get("falsifiers_secondary", []):
        print("FALSIFIER: " + f)
    print(f"VERDICT: {r['verdict']}")
    if a.out:
        Path(a.out).write_text(json.dumps({k: v for k, v in r.items() if k != "lines"}, indent=1, default=str))
    return 3 if r["verdict"] == "FALSIFIED" else 0


if __name__ == "__main__":
    sys.exit(main())
