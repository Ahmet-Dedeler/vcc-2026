"""429 known answers: compare every run429 output tag that has a reference row (stats arrays bit-for-bit, members to
<= 1e-12). Exit 1 if a MUST row fails; SHOULD rows (zero repair == no-op on H1) are reported, never abort.
usage: uv run python ka429.py [TAG ...]   (default: every metrics-*.json in work/h1)"""
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
W = HERE / "work/h1"
MDE = "mde0.1cap0.1min0.3"


def refs(tag):
    """(kind, reference metrics path) pairs for one of our tags."""
    out = []
    m = re.fullmatch(r"base_s(\d)_bps", tag)
    if m:
        out.append(("MUST 416 base", EXP / f"416-de-carrier-menu/work/h1/metrics-base_s{m[1]}_bps.json"))
    m = re.fullmatch(r"base_s(\d)_bps_dnk250zs", tag)
    if m:
        out.append(("MUST 423 dn_k250_zs", EXP / f"423-rank-preserving-magnitude/work/h1/metrics-dn_k250_zs_s{m[1]}_bps.json"))
    m = re.fullmatch(rf"base_s(\d)_bps_{MDE}", tag)
    if m:
        out.append(("MUST 426 MDE", EXP / f"426-mde-magnitude/work/h1/metrics-base_s{m[1]}_bps_{MDE}.json"))
    m = re.fullmatch(rf"base_s(\d)_bps_zr-5_{MDE}(_dnk250zs)?", tag)
    if m:
        out.append(("SHOULD zr no-op", W / f"metrics-base_s{m[1]}_bps_{MDE}{m[2] or ''}.json"))
    m = re.fullmatch(r"(.*)_dnk250zsrp", tag)
    if m:
        out.append(("MUST vendored == rp423", W / f"metrics-{m[1]}_dnk250zs.json"))
    return out


def compare(tag, ref):
    a = json.loads((W / f"metrics-{tag}.json").read_text())["raw"]
    b = json.loads(ref.read_text())["raw"]
    dmax = max(abs(a[k] - b[k]) for k in a)
    sa = np.load(W / f"stats-{tag}.npz", allow_pickle=True)
    sb = np.load(ref.with_name(ref.name.replace("metrics-", "stats-").replace(".json", ".npz")), allow_pickle=True)
    bad = [k for k in ("targets", "sums", "jk", "p", "q", "lfc")
           if not np.array_equal(sa[k], sb[k], **({} if k == "targets" else {"equal_nan": True}))]
    return dmax, bad


def main():
    tags = sys.argv[1:] or sorted(p.name[len("metrics-"):-len(".json")] for p in W.glob("metrics-*.json"))
    fail = 0
    for t in tags:
        for kind, ref in refs(t):
            if not ref.exists():
                print(f"[SKIP] {t} vs {kind}: {ref.name} not there yet")
                continue
            dmax, bad = compare(t, ref)
            ok = dmax <= 1e-12 and not bad
            print(f"[{'PASS' if ok else 'FAIL'}] {t} vs {kind} ({ref.relative_to(EXP)}): members max|d| {dmax:.3g}; "
                  f"stats arrays differing: {bad or 'none'}")
            fail += (not ok) and kind.startswith("MUST")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
