#!/usr/bin/env python3
"""625 board reader: captured entries -> the module decode -> the octile flip mask (column space -> octile space, asserted).

For the A/B/C rehearsal (transfer rule, PREREGISTRATION §6) and the approved D/E/F decode (in-sample rule, final phase).

Manifest (JSON):
  {"design": "v8" | "single",
   "base": "<entry id | path>",
   "rows":    ["<id>", ...],              # v8: the 8 Hadamard rows, in modules625_v8.npz `flips` order (COLUMN space)
   "checks":  ["<id>", "<id>"],           # v8: the 2 additivity arms, `check_flips` order (octile 0, then octile 7)
   "octiles": {"<k>": "<id>", ...},       # single: one entry per flipped octile k (OCTILE space), any subset of 0..7
   "combo":   {"mask": "<8 chars, octile order>", "entry": "<id>"}}   # optional: the scored combination
Entries: captured-entry JSON (scripts/capture_board_entry.py's wrapper, or a bare entry object), or an entry id resolved
under data/board-entries/<id>.json.

Before any decode:
  * every entry is terminal, members not hidden, and on the base's panel and partition;
  * PRIOR-ONLY: pds_cosine and expr_mse_unbiased_capped_norm equal the base's to 1e-12. That makes the reads exact
    (468, shared emission draws). A changed bulk is refused unless --allow-bulk-change, in which case the bulk delta is
    printed and the DE reads carry the re-roll.
Decode (raw members, board-pooled; overall = score_avg):
  * v8: decode625.decode on nmae / fid / reach gives x per COLUMN, plus the additivity check (§5). score_avg and raw jac
    are decoded with the same algebra (reported). Translation to octile order uses the npz's colmap, asserted equal to
    [3, 0, 1, 2, 4, 5, 6, 7].
  * single: x_k = m(entry_k) - m(base), directly in octile space, with no additivity assumption.
Rules:
  * --rule transfer (A/B/C; PREREGISTRATION §6): flip k iff x^nmae < 0 AND x^fid > 0 AND |x^nmae| >= 0.25 × median|x^nmae|
    (the median is over the decoded / probed octiles). v8 with additivity failing: decode625's fallback.
  * --rule transfer3 (Amendment 3, the CURRENT A/B/C rule): flip k iff x^nmae < 0 AND |x^nmae| >= max(0.25 × median|x^nmae|,
    2 σ_k) AND x^fid >= -0.002. σ_k is each octile's realisation sd of Δnmae on the board's scale (--sigma <json {octile: σ}>;
    H1's 625d seed sd × sqrt(85/900)). fid is a no-harm check only, because DOWN -> UP flips expand the call set.
  * --rule insample (D/E/F, the final's own truth): flip k iff x^overall > 0 (the board's score_avg delta). The final
    package is the combination, and it must itself be scored before Nov 5 (FINAL-PHASE-FEEDBACK-PLAN §1 rule 2).
  * combo (optional): its measured delta vs the base, and vs the sum of the chosen x's (additivity of the chosen set).
Output: a per-octile table, the mask (8 chars, octile order), the reading, and with --out a JSON record for
docs/october-decisions.json (entry ids, modules625_v8.npz sha, this file's sha).

  uv run python experiments/625-module-readout/readdecode625.py --manifest <json> --rule transfer|insample [--out <json>]
  uv run python experiments/625-module-readout/readdecode625.py --selftest
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ENTRIES = ROOT / "data/board-entries"
sys.path.insert(0, str(HERE))
import decode625 as D  # noqa: E402

V8_PATH = HERE / "work/modules625_v8.npz"
V8_SHA256 = "c24e3304ffd63d77d76b4a74341e09a2ee8396d7e372955e9c6feee167133dc1"
COLMAP = [3, 0, 1, 2, 4, 5, 6, 7]
RAW = {"nmae": "de_wilcoxon_lfc_nmae", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
       "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard", "overall": "score_avg"}
BULK = ("pds_cosine", "expr_mse_unbiased_capped_norm")


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_entry(ref):
    """A captured-entry wrapper, a bare entry object, a path, or an id under data/board-entries/."""
    if isinstance(ref, dict):
        obj = ref
    else:
        p = Path(ref)
        if not p.exists():
            p = ENTRIES / f"{ref}.json"
        obj = json.loads(p.read_text())
    e = obj.get("entry", obj)
    hidden = bool(obj.get("members_hidden", False))
    return e, hidden


def members(e) -> dict:
    out = {}
    for k, col in RAW.items():
        v = e.get(col)
        if v is None:
            raise SystemExit(f"entry {e.get('entry_id')}: member {col} absent (members hidden or not scored)")
        out[k] = float(v)
    return out


def gate(base_e, e, allow_bulk: bool, label: str) -> dict:
    if not e.get("is_terminal", True):
        raise SystemExit(f"{label} {e.get('entry_id')}: not terminal")
    for k in ("panel_id", "partition"):
        if e.get(k) != base_e.get(k):
            raise SystemExit(f"{label} {e.get('entry_id')}: {k} {e.get(k)!r} != base {base_e.get(k)!r}")
    db = {k: float(e[k]) - float(base_e[k]) for k in BULK}
    if any(abs(v) > 1e-12 for v in db.values()) and not allow_bulk:
        raise SystemExit(f"{label} {e.get('entry_id')}: NOT prior-only (bulk moved {db}); refusing (use --allow-bulk-change "
                         f"for a fused entry; its DE reads then carry the re-roll)")
    return db


def decide(x: dict, probed: list[int], rule: str, fallback_mask=None, sigma=None) -> list[int]:
    xn, xf, xo = (np.array([x[m][k] for k in probed]) for m in ("nmae", "fid", "overall"))
    if rule == "insample":
        sel = xo > 0
    elif rule == "transfer3":
        sg = np.array([float((sigma or {}).get(str(k), 0.0)) for k in probed])
        bar = np.maximum(0.25 * float(np.median(np.abs(xn))), 2.0 * sg)
        sel = (xn < 0) & (np.abs(xn) >= bar) & (xf >= -0.002)
    elif fallback_mask is not None:
        sel = np.array([fallback_mask[k] for k in probed], bool)
    else:
        bar = 0.25 * float(np.median(np.abs(xn)))
        sel = (xn < 0) & (xf > 0) & (np.abs(xn) >= bar)
    mask = [0] * 8
    for k, s in zip(probed, sel):
        mask[k] = int(bool(s))
    return mask


def read(man: dict, rule: str, allow_bulk: bool = False, quiet: bool = False, sigma=None) -> dict:
    assert sha256(V8_PATH) == V8_SHA256, "modules625_v8.npz changed"
    z = np.load(V8_PATH, allow_pickle=True)
    assert list(map(int, z["colmap"])) == COLMAP, "colmap changed"
    base_e, hid = load_entry(man["base"])
    if hid:
        raise SystemExit("base: members hidden")
    mb = members(base_e)
    rec = {"design": man["design"], "rule": rule, "base": base_e.get("entry_id"), "bulk_deltas": {}}
    x = {m: [np.nan] * 8 for m in RAW}
    fallback = None
    if man["design"] == "v8":
        D.use_design("v8")
        rows = [load_entry(r) for r in man["rows"]] + [load_entry(r) for r in man["checks"]]
        if len(rows) != 10:
            raise SystemExit("v8 needs 8 rows + 2 checks")
        arms = []
        for i, (e, h) in enumerate(rows):
            if h:
                raise SystemExit(f"arm {i}: members hidden")
            rec["bulk_deltas"][str(e.get("entry_id"))] = gate(base_e, e, allow_bulk, f"arm {i}")
            arms.append(members(e))
        dec = D.decode({"base": mb, "arms": arms})                 # nmae / fid / reach + additivity + §6 (column space)
        xc = {m: dec[m]["x"] for m in ("nmae", "fid", "reach")}
        for m in ("jac", "overall"):
            xc[m] = D.decode_member(mb[m], [a[m] for a in arms]).tolist()
        for m in RAW:
            for j in range(8):
                x[m][COLMAP[j]] = float(xc[m][j])                 # column j -> octile COLMAP[j]
        rec["additive"] = dec["additive"]
        rec["additivity"] = {m: dec[m]["additivity"] for m in ("nmae", "fid")}
        if not dec["additive"]:
            fallback = [0] * 8
            for j in range(8):
                fallback[COLMAP[j]] = int(dec["transfer_mask"][j])
        probed = list(range(8))
        rec["decode625_mask_column_order"] = dec["transfer_mask"]
    elif man["design"] == "single":
        probed = sorted(int(k) for k in man["octiles"])
        for k in probed:
            e, h = load_entry(man["octiles"][str(k)])
            if h:
                raise SystemExit(f"octile {k}: members hidden")
            rec["bulk_deltas"][str(e.get("entry_id"))] = gate(base_e, e, allow_bulk, f"octile {k}")
            mk = members(e)
            for m in RAW:
                x[m][k] = mk[m] - mb[m]
        rec["additive"] = None
    else:
        raise SystemExit(f"unknown design {man['design']!r}")
    mask = decide(x, probed, rule, fallback if rule == "transfer" else None, sigma=sigma)
    rec["sigma"] = sigma
    rec.update(x={m: x[m] for m in RAW}, probed=probed, mask="".join(map(str, mask)), n_flip=int(sum(mask)),
               predicted_overall_gain=float(sum(x["overall"][k] for k in range(8) if mask[k])),
               reading=("ABORT-additivity (fallback mask)" if (rec["additive"] is False and rule in ("transfer", "transfer3"))
                        else ("GO" if any(mask) else "NULL")))
    if "combo" in man:
        ce, h = load_entry(man["combo"]["entry"])
        if h:
            raise SystemExit("combo: members hidden")
        rec["bulk_deltas"][str(ce.get("entry_id"))] = gate(base_e, ce, allow_bulk, "combo")
        cm = members(ce)
        cmask = [int(c) for c in man["combo"]["mask"]]
        pred = {m: float(sum(x[m][k] for k in range(8) if cmask[k])) for m in RAW}
        rec["combo"] = dict(mask=man["combo"]["mask"], measured={m: cm[m] - mb[m] for m in RAW}, predicted_additive=pred,
                            matches_decoded_mask=(man["combo"]["mask"] == rec["mask"]))
    rec["sha256"] = {"modules625_v8": V8_SHA256, "readdecode625": sha256(Path(__file__))}
    if not quiet:
        print(f"design {man['design']}, rule {rule}, base {rec['base']}; additive {rec['additive']}")
        print("octile | x nmae (raw) | x fid (raw) | x reach | x overall (score_avg) | flip")
        for k in range(8):
            f = "-" if k not in probed else ("FLIP" if mask[k] else "keep")
            print(f"  {k}    | {x['nmae'][k]:+.6f} | {x['fid'][k]:+.6f} | {x['reach'][k]:+.6f} | {x['overall'][k]:+.6f} | {f}")
        print(f"mask (octile order 0..7): {rec['mask']}  reading {rec['reading']}  predicted overall gain "
              f"{rec['predicted_overall_gain']:+.6f}")
        if "combo" in rec:
            c = rec["combo"]
            print(f"combo {c['mask']}: measured overall {c['measured']['overall']:+.6f} vs additive prediction "
                  f"{c['predicted_additive']['overall']:+.6f}; matches the decoded mask: {c['matches_decoded_mask']}")
    return rec


def _fake(base, d, eid, bulk_shift=0.0):
    e = {"entry_id": eid, "is_terminal": True, "panel_id": "p", "partition": "val",
         "pds_cosine": 0.86 + bulk_shift, "expr_mse_unbiased_capped_norm": 0.62}
    for m, col in RAW.items():
        e[col] = base[m] + d.get(m, 0.0)
    return {"entry": e, "members_hidden": False}


def selftest() -> None:
    rng = np.random.default_rng(6255)
    z = np.load(V8_PATH, allow_pickle=True)
    flips, checks = z["flips"].astype(bool), z["check_flips"].astype(bool)
    base = {"nmae": 0.90, "fid": 0.51, "reach": 0.28, "jac": 0.04, "overall": 0.27}
    n_ok = 0
    for t in range(300):
        xo = {m: rng.normal(0, 1e-3, 8) for m in RAW}           # true effects in OCTILE space
        if t % 3 == 0:                                         # a designed case: octiles 4 and 5 disagree with CLEAN4
            xo["nmae"][:] = 2e-4; xo["fid"][:] = -1e-4; xo["overall"][:] = -1e-4
            for k in (4, 5):
                xo["nmae"][k] = -6e-4; xo["fid"][k] = 3e-4; xo["overall"][k] = 2e-4
        xcol = {m: np.array([xo[m][COLMAP[j]] for j in range(8)]) for m in RAW}   # column j holds octile COLMAP[j]
        bE = _fake(base, {}, "B")
        rows = [_fake(base, {m: float(r @ xcol[m]) for m in RAW}, f"r{i}") for i, r in enumerate(flips)]
        chks = [_fake(base, {m: float(r @ xcol[m]) for m in RAW}, f"c{i}") for i, r in enumerate(checks)]
        man = {"design": "v8", "base": bE, "rows": rows, "checks": chks}
        r = read(man, "transfer", quiet=True)
        for m in RAW:
            assert np.allclose(r["x"][m], xo[m], atol=1e-12), (t, m)
        assert r["additive"] is True
        exp = (xo["nmae"] < 0) & (xo["fid"] > 0) & (np.abs(xo["nmae"]) >= 0.25 * np.median(np.abs(xo["nmae"])))
        assert r["mask"] == "".join(map(str, exp.astype(int))), (t, r["mask"])
        if t % 3 == 0:
            assert r["mask"] == "00001100", r["mask"]
            ri = read(man, "insample", quiet=True)
            assert ri["mask"] == "00001100" and abs(ri["predicted_overall_gain"] - 4e-4) < 1e-12
        # single design on the same truth, a subset of octiles
        sub = sorted(rng.choice(8, 5, replace=False).tolist())
        mans = {"design": "single", "base": bE, "octiles": {str(k): _fake(base, {m: float(xo[m][k]) for m in RAW}, f"s{k}")
                                                             for k in sub}}
        rs = read(mans, "insample", quiet=True)
        assert rs["mask"] == "".join("1" if (k in sub and xo["overall"][k] > 0) else "0" for k in range(8)), rs["mask"]
        n_ok += 1
    # refusals: a fused entry (bulk moved) and hidden members
    bE = _fake(base, {}, "B")
    bad = {"design": "single", "base": bE, "octiles": {"3": _fake(base, {"nmae": -1e-4}, "x", bulk_shift=1e-6)}}
    try:
        read(bad, "insample", quiet=True)
        raise AssertionError("a fused entry was not refused")
    except SystemExit as e:
        assert "NOT prior-only" in str(e)
    r = read(bad, "insample", allow_bulk=True, quiet=True)
    assert abs(list(r["bulk_deltas"].values())[0]["pds_cosine"] - 1e-6) < 1e-12
    hid = _fake(base, {}, "h"); hid["members_hidden"] = True
    try:
        read({"design": "single", "base": bE, "octiles": {"1": hid}}, "insample", quiet=True)
        raise AssertionError("hidden members were not refused")
    except SystemExit as e:
        assert "hidden" in str(e)
    # combo bookkeeping
    xo = {m: np.zeros(8) for m in RAW}; xo["overall"][[2, 6]] = [3e-4, 1e-4]
    mans = {"design": "single", "base": bE, "octiles": {str(k): _fake(base, {m: float(xo[m][k]) for m in RAW}, f"s{k}") for k in range(8)},
            "combo": {"mask": "00100010", "entry": _fake(base, {"overall": 3.5e-4}, "combo")}}
    rc = read(mans, "insample", quiet=True)
    assert rc["mask"] == "00100010" and rc["combo"]["matches_decoded_mask"]
    assert abs(rc["combo"]["predicted_additive"]["overall"] - 4e-4) < 1e-12 and abs(rc["combo"]["measured"]["overall"] - 3.5e-4) < 1e-12
    # Amendment 3 rule: nmae-only with a sigma bar and a fid no-harm check
    xo = {m: np.zeros(8) for m in RAW}
    xo["nmae"][:] = [2e-4, -8e-4, -3e-4, 1e-4, -9e-4, -6e-4, 2e-4, -1e-3]
    xo["fid"][:] = [0, -1e-3, 1e-4, 0, -3e-3, 5e-4, 0, -1e-4]
    mans = {"design": "single", "base": bE, "octiles": {str(k): _fake(base, {m: float(xo[m][k]) for m in RAW}, f"s{k}") for k in range(8)}}
    sig = {"1": 1e-4, "2": 2e-4, "4": 1e-4, "5": 1e-4, "7": 6e-4}
    r3 = read(mans, "transfer3", quiet=True, sigma=sig)
    # bar = max(0.25*median|x| = 0.25*4.5e-4 = 1.125e-4, 2 sigma): oct1 1.6e-4 ok; oct2 4e-4 > 3e-4 -> no; oct4 fid -3e-3 -> no;
    # oct5 2e-4 ok; oct7 1.2e-3 > 1e-3 -> no
    assert r3["mask"] == "01000100", r3["mask"]
    print(f"SELFTEST PASS: {n_ok} random v8 decodes recover the octile-space truth to 1e-12 through the column map; the §6 "
          f"rule and the in-sample rule exact; the single design exact on random subsets; fused / hidden entries refused; "
          f"combo bookkeeping exact; Amendment 3 transfer3 rule exact (sigma bar, fid no-harm)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path)
    ap.add_argument("--rule", choices=("transfer", "transfer3", "insample"))
    ap.add_argument("--sigma", type=Path, help="transfer3: JSON {octile: sigma of the board-scale Δnmae}")
    ap.add_argument("--allow-bulk-change", action="store_true")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    if a.manifest is None or a.rule is None:
        ap.error("--manifest and --rule are required (or --selftest)")
    sigma = json.loads(a.sigma.read_text()) if a.sigma else None
    if a.rule == "transfer3" and sigma is None:
        ap.error("--rule transfer3 needs --sigma (Amendment 3's realisation noise per octile)")
    rec = read(json.loads(a.manifest.read_text()), a.rule, a.allow_bulk_change, sigma=sigma)
    if a.out:
        a.out.write_text(json.dumps(rec, indent=1, default=float))
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
