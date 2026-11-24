"""804 — the checker's calibration (experiments/800-checker/ka800.json) from the graded known answers. Light.

Per class and surface set: rho = checker point / board delta over the graded rows with |board| >= 3e-4 and matching sign
(stats800.calibrate); sign failures counted. Per member (bulk rows): rho_pds (literal), rho_pds_lam05 (board-emulated, a
diagnostic chosen AFTER the first read, so not a validated unit), rho_mse (518's identity). Prospective rows never enter.
    uv run python experiments/804-known-answers/calib804.py
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/800-checker"))
import stats800 as S  # noqa: E402

G = json.loads((HERE / "work/grade804.json").read_text())
CLS = {"495": "bulk-key", "495pct": "bulk-key", "479p2g": "bulk-shape", "479sp2": "bulk-shape", "498p2g": "bulk-shape",
       "702": "bulk-external", "707": "bulk-direction"}   # 707 added 10-17 after its board read (prospective row graded first)
out = {"calibration": {}, "members": {}, "rows": {}, "note": "rho = checker / board; projection = checker / rho (DESIGN 7)"}
for surf in ("W", "ABC"):
    rows = []
    for k, r in G.items():
        row, s = k.split("|")
        if s != surf or row not in CLS or not r.get("board") or row == "495pct":
            continue
        rows.append(dict(**{"class": CLS[row]}, row=row, checker=r["point"], board=r["board"]["overall"]))
        out["rows"][k] = dict(row=row, surface=surf, cls=CLS[row], checker=r["point"], board=r["board"]["overall"],
                              rho=r.get("rho_overall"), rho_pds=r.get("rho_pds"), rho_mse=r.get("rho_mse"),
                              decision=r["decision"], grade=r["grade"])
    for cls, v in S.calibrate(rows).items():
        out["calibration"][f"{cls}|{surf}"] = v
    allr = [dict(**{"class": "bulk-any"}, row=x["row"], checker=x["checker"], board=x["board"]) for x in rows]
    out["calibration"][f"bulk-any|{surf}"] = S.calibrate(allr)["bulk-any"]
# ABCp (801's A/B/C folds with day-recipe held-out operators; grade804p.py): its own calibration, keyed "|ABCp"
GP = json.loads((HERE / "work/grade804p.json").read_text()) if (HERE / "work/grade804p.json").exists() else {}
rows = []
for row, r in GP.items():
    if row not in CLS or not r.get("board") or row == "495pct":
        continue
    rows.append(dict(**{"class": CLS[row]}, row=row, checker=r["point"], board=r["board"]["overall"]))
    out["rows"][f"{row}|ABCp"] = dict(row=row, surface="ABCp", cls=CLS[row], checker=r["point"], board=r["board"]["overall"],
                                      rho=r.get("rho_overall"), decision=r["decision"], grade=r["grade"])
# board-pointed composite arms read through 801's ABCp folds by 811's composite811 (6ec25788): the checker numbers are 811's
# information reads (811/RESULTS.md §000000), the board numbers the def_day reads of 10-20 (BOARD-LEDGER; PREREG-C1/C0-ABC)
EXTRA_ABCP = [dict(row="C1-811 (C1 - B1)", cls="bulk-direction", checker=0.00153, board=0.0028713),
              dict(row="C0-811 (C0 - B0)", cls="bulk-direction", checker=0.00269, board=0.0040775)]
for x in EXTRA_ABCP:
    rows.append(dict(**{"class": x["cls"]}, row=x["row"], checker=x["checker"], board=x["board"]))
    out["rows"][f"{x['row']}|ABCp"] = dict(row=x["row"], surface="ABCp", cls=x["cls"], checker=x["checker"], board=x["board"],
                                           rho=x["checker"] / x["board"], source="811 composite811 info read; board def_day 10-20")
if rows:
    for cls, v in S.calibrate(rows).items():
        out["calibration"][f"{cls}|ABCp"] = v
    allr = [dict(**{"class": "bulk-any"}, row=x["row"], checker=x["checker"], board=x["board"]) for x in rows]
    out["calibration"]["bulk-any|ABCp"] = S.calibrate(allr)["bulk-any"]
# DE rows (803's engine on h1test + h1tv): the frozen regime f 0.35, and f 0.61 for call-set size (662 / 619's amendment;
# chosen after the f 0.35 read, so recorded as a post-hoc class rule, not a validated unit)
GD = json.loads((HERE / "work/grade_de804.json").read_text()) if (HERE / "work/grade_de804.json").exists() else {}
DECLS = {"446bj2z zs": "de-mag", "DN": "de-mag", "630arm UPRN": "de-mag", "496a hook (618 form)": "de-sign",
         "idall3 vs hook": "de-sign", "K20": "de-callset", "K16": "de-callset"}
for regime in ("f0.35", "f0.61"):
    rows = []
    for row, cls in DECLS.items():
        r = GD.get(f"{row}|{regime}")
        if not r:
            continue
        rows.append(dict(**{"class": cls}, row=row, checker=r["point"], board=r["board"]))
        out["rows"][f"{row}|DE-H1|{regime}"] = dict(row=row, surface="DE-H1", regime=regime, cls=cls, checker=r["point"],
                                                     board=r["board"], decision=r["decision"], grade=r.get("grade"))
    for cls, v in S.calibrate(rows).items():
        out["calibration"][f"{cls}|DE-H1|{regime}"] = v
_prev = ROOT / "experiments/800-checker/ka800.json"   # carry the hand-written post-hoc notes across regenerations (10-18)
if _prev.exists() and json.loads(_prev.read_text()).get("post_hoc"):
    out["post_hoc"] = json.loads(_prev.read_text())["post_hoc"]
(ROOT / "experiments/800-checker/ka800.json").write_text(json.dumps(out, indent=1))
for k, v in out["calibration"].items():
    print(f"{k:22s} {v}")
