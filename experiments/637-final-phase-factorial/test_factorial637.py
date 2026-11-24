"""Tests for 637's factorial reader. Synthetic captures built with 531's own helpers + one smoke read on real captures,
then 531's full suite to prove its verdicts are untouched.

  uv run python experiments/637-final-phase-factorial/test_factorial637.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/531-cloud-final-phase-tooling"))
import readdef as R  # noqa: E402
import readdef_factorial637 as F  # noqa: E402
from test_readdef import entry, shifted  # noqa: E402  (531's synthetic-capture helpers)

N = 0


# The original synthetic cases (2026-10-02) were built around the v1 band; they run under v1 unchanged.
F.set_leak_band("v1")


def check(name, cond, detail=""):
    global N
    if not cond:
        raise SystemExit(f"FAIL {name} {detail}")
    N += 1
    print(f"PASS {name}")


def combo(base, eid, bulk_scaled=None, de_scaled=None, leak=-0.0003, raw=None):
    """base + a bulk switch (scaled pds/mse) + a DE switch (scaled DE members) + the bulk switch's leak, put on nmae."""
    sd = {"pds": 0.0, "mse": 0.0, "fid": 0.0, "reach": 0.0, "jac": 0.0, "nmae": 0.0}
    for k, v in (bulk_scaled or {}).items():
        sd[k] += v
    for k, v in (de_scaled or {}).items():
        sd[k] += v
    sd["nmae"] += 6 * leak            # leak in overall units -> scaled nmae
    return shifted(base, eid, draw=raw or {}, sdraw=sd)


def main() -> int:
    base = entry("base", overall=0.30, raw={"pds": 0.86, "mse": 0.62, "fid": 0.514, "reach": 0.286, "jac": 0.041,
                                            "nmae": 0.897},
                 scaled={"pds": 0.81, "mse": 0.38, "fid": 0.006, "reach": 0.234, "jac": 0.028, "nmae": 0.168})

    # 1. additive: bulk +0.001 overall (scaled pds+mse = 0.006), DE fnr +0.003 overall (scaled fid 0.018)
    e = combo(base, "f1", bulk_scaled={"pds": 0.002, "mse": 0.004}, de_scaled={"fid": 0.018},
              raw={"pds": 0.001, "mse": -0.004})
    v = F.factorial_verdict("4", "6", e, base)
    check("additive: bulk row 4 adopts", v["bulk"]["adopt"] is True and v["bulk"]["code"] == "adopt_arm", v["bulk"])
    check("additive: DE row 6 adopts", v["de"]["adopt"] is True and v["de"]["code"] == "adopt", v["de"])
    check("additive: no inconclusive", "inconclusive" not in v["bulk"]["code"] + v["de"]["code"])

    # 2. DE near the row-6 bin edge: true fnr +0.0013 -> fnr - leak spans the +0.0015 bin inside the band
    e = combo(base, "f2", bulk_scaled={"pds": 0.002, "mse": 0.004}, de_scaled={"fid": 6 * 0.0013 + 6 * 0.0003 - 6 * 0.0003},
              leak=0.0, raw={"pds": 0.001, "mse": -0.004})
    v = F.factorial_verdict("4", "6", e, base)
    check("edge: DE inconclusive inside the leak band", v["de"]["code"] == "inconclusive_leak_band", v["de"])
    check("edge: inconclusive never adopts", v["de"]["adopt"] is False)

    # 3. both lose: bulk -0.0005 overall, DE harm (row 7: overall <= -0.0005 after removing the leak)
    e = combo(base, "f3", bulk_scaled={"pds": -0.001, "mse": -0.002}, de_scaled={"reach": -0.012},
              raw={"pds": -0.0005, "mse": 0.002})
    v = F.factorial_verdict("4", "7", e, base)
    check("lose: bulk base stands", v["bulk"]["adopt"] is False, v["bulk"])
    check("lose: DE row 7 harm", v["de"]["code"] == "harm", v["de"])

    # 4. row 2 reads bulk only (no leak charged); row 3 reads overall = bulk + leak
    e = combo(base, "f4", bulk_scaled={"pds": 0.0003, "mse": 0.0003}, de_scaled={}, raw={"pds": 0.0002})
    v2 = F.factorial_verdict("2", "6", e, base)
    v3 = F.factorial_verdict("3", "6", e, base)
    check("row 2 adopts a +0.0001 bulk step (exact read)", v2["bulk"]["adopt"] is True, v2["bulk"])
    check("row 3 charges the leak: +0.0001 bulk - 0.0003..-0.0005 leak does not adopt", v3["bulk"]["adopt"] is False,
          v3["bulk"])

    # 5. row 7 falsifier on jac still fires in factorial mode
    e = combo(base, "f5", bulk_scaled={"pds": 0.002}, de_scaled={"fid": 0.03}, raw={"pds": 0.001, "jac": 0.002})
    v = F.factorial_verdict("4", "7", e, base)
    check("row 7 jac falsifier still fires", v["de"]["code"] == "bug", v["de"])

    # 6. refusals
    for row in ("5", "8", "9", "8-9", "10", "1", "13"):
        try:
            F.factorial_verdict("4", row, base, base)
            check(f"refuse DE row {row}", False)
        except SystemExit as ex:
            check(f"refuse DE row {row}", "REFUSED" in str(ex) or "must be" in str(ex), str(ex))
    for row in ("6", "7", "5"):
        try:
            F.factorial_verdict(row, "6", base, base)
            check(f"refuse bulk row {row}", False)
        except SystemExit as ex:
            check(f"refuse bulk row {row}", True, str(ex))

    # 7. hidden members
    hid = entry("hid", overall=0.3, hide=("score_fid",))
    v = F.factorial_verdict("4", "6", hid, base)
    check("hidden members -> no read", v["code"] == "hidden")

    # 8. validation PASS / FAIL
    X = combo(base, "X", bulk_scaled={"pds": 0.002, "mse": 0.004}, de_scaled={}, raw={"pds": 0.001, "mse": -0.004})
    Y = combo(base, "Y", bulk_scaled={}, de_scaled={"fid": 0.018}, leak=0.0)
    FCT = combo(base, "FCT", bulk_scaled={"pds": 0.002, "mse": 0.004}, de_scaled={"fid": 0.018},
                raw={"pds": 0.001, "mse": -0.004})
    check("validate PASS on an additive triple", F.validate(FCT, X, Y, base)["pass"] is True)
    bad_bulk = combo(base, "Fb", bulk_scaled={"pds": 0.002, "mse": 0.004}, de_scaled={"fid": 0.018},
                     raw={"pds": 0.001 + 2e-6, "mse": -0.004})
    check("validate FAIL when pds is not additive (2e-6)", F.validate(bad_bulk, X, Y, base)["pass"] is False)
    bad_de = combo(base, "Fd", bulk_scaled={"pds": 0.002, "mse": 0.004}, de_scaled={"fid": 0.018 + 6 * 0.001},
                   raw={"pds": 0.001, "mse": -0.004})
    check("validate FAIL when DE misses by 0.001 overall", F.validate(bad_de, X, Y, base)["pass"] is False)

    # 9. smoke read on REAL captures: 479p2g (bulk-only step) vs 446bj2 read as factorial rows 4 x 6
    import json
    pk = {}  # locate captures by package name (entry ids are opaque)
    for p in (ROOT / "data/board-entries").glob("*.json"):
        d = json.loads(p.read_text())
        pk[Path(d.get("package") or "").name[:7]] = p
    if "479p2g-" in pk and "446bj2-" in pk:
        a, b = R.load_entry(str(pk["479p2g-"])), R.load_entry(str(pk["446bj2-"]))
        v = F.factorial_verdict("4", "6", a, b)
        check("real 479p2g vs 446bj2: bulk switch adopts (it was +0.0012 bulk)", v["bulk"]["adopt"] is True, v["bulk"])
        check("real 479p2g vs 446bj2: no DE switch in it reads below the row-6 bin", v["de"]["adopt"] is False, v["de"])
    else:
        print("SKIP real-capture smoke (captures not found)")

    # 11. the v2 leak band (2026-10-03; 658's 9 bulk-only board pairs): constants, recording, and the validate window
    F.set_leak_band("v2")
    check("v2 band: centre -0.0002, band (+0.0002, -0.0002, -0.0006)",
          F.LEAK_CENTRE == -0.0002 and F.LEAK_BAND == (0.0002, -0.0002, -0.0006))
    Fm = combo(base, "Fm", bulk_scaled={"pds": 0.002, "mse": 0.004}, de_scaled={"fid": 0.018 + 6 * 0.00035},
               raw={"pds": 0.001, "mse": -0.004}, leak=0.0)
    rv2 = F.validate(Fm, X, Y, base)
    check("v2: a +0.00035 overall measured leak validates (|0.00035 + 0.0002| = 0.00055 <= 0.0006)",
          rv2["pass"] is True and rv2["leak_band"] == "v2", rv2["lines"][-1])
    F.set_leak_band("v1")
    rv1 = F.validate(Fm, X, Y, base)
    check("v1: the same triple fails (|0.00035 + 0.0003| = 0.00065 > 0.0006) and records band v1",
          rv1["pass"] is False and rv1["leak_band"] == "v1", rv1["lines"][-1])
    F.set_leak_band("v2")
    vv = F.factorial_verdict("4", "6", combo(base, "e2", bulk_scaled={"pds": 0.003, "mse": 0.005}, de_scaled={"fid": 0.03},
                                             raw={"pds": 0.0015, "mse": -0.005}, leak=-0.0002), base)
    check("v2: a factorial read records its leak band", vv.get("leak_band") == "v2", vv.get("leak_band"))
    try:
        F.set_leak_band("v9")
        check("unknown leak band refused", False)
    except SystemExit:
        check("unknown leak band refused", True)

    # 10. 531's suite is untouched
    r = subprocess.run([sys.executable, str(ROOT / "experiments/531-cloud-final-phase-tooling/test_readdef.py")],
                       capture_output=True, text=True)
    check("531's suite still passes 74 row-verdicts", r.returncode == 0 and "74 row-verdicts passed" in r.stdout,
          r.stdout[-300:] + r.stderr[-300:])
    print(f"\n{N} factorial checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
