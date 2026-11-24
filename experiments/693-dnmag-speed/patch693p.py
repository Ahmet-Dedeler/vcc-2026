#!/usr/bin/env python3
"""693p — pin bump: the builder uses dn_mag v2 (src/vcc2026/dn_mag2.py, byte-identical outputs, 3-6x faster). UNAPPLIED.

Refuses unless 693's six KA records exist (work/ka693_s{0,1,2}_{dn,updn}.json), each says PASS (blocks 85, v2 == v1), and
each names the same v2 sha256 as this directory's dn_mag2.py. v1 (src/vcc2026/dn_mag.py, 0a8cc560) stays on disk untouched,
so every earlier build stays reproducible.
Changes:
  src/vcc2026/dn_mag2.py                      new file (a byte copy of experiments/693-dnmag-speed/dn_mag2.py)
  scripts/make_submission_pm_cd4.py           DE_DNMAG_MODULE_PATH -> dn_mag2.py, DE_DNMAG_MODULE_SHA256 -> v2, both imports
                                              `from vcc2026 import dn_mag as _dn_mag` -> `dn_mag2 as _dn_mag`, the comment
  scripts/def_day.sh                          690a's final-meta module check (sha + message) and its doc line -> v2
After apply (the coordinator's KA): the all-unset KA (== ec3a5ae6) and a DN rebuild whose prediction content equals dn631's
(3c6a51b8's content digest; the meta's de_dnmag.module_sha256 moves to v2 by design).
    uv run python experiments/693-dnmag-speed/patch693p.py --root <repo> [--check]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
V1 = "0a8cc560231a16f0fff7fbbf316e3a5e190cbdf20490239ef354b49b3ff60dac"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--ka-dir", type=Path, default=HERE / "work", help="where ka693_s*_*.json live (tests only)")
    a = ap.parse_args()
    r = a.root.resolve()
    v2 = sha(HERE / "dn_mag2.py")
    errs = []
    for s in (0, 1, 2):
        for m in ("dn", "updn"):
            f = a.ka_dir / f"ka693_s{s}_{m}.json"
            if not f.exists():
                errs.append(f"missing KA record {f.name}")
                continue
            k = json.loads(f.read_text())
            if k.get("blocks") != 85 or k.get("v2_sha") != v2 or k.get("v1_sha") != V1:
                errs.append(f"KA record {f.name} does not cover 85 blocks of this v2 ({str(k.get('v2_sha'))[:16]}) vs v1")
    bu, dd = r / "scripts/make_submission_pm_cd4.py", r / "scripts/def_day.sh"
    b, d = bu.read_text(), dd.read_text()
    if (r / "src/vcc2026/dn_mag2.py").exists() or "dn_mag2" in b:
        errs.append("693p is already applied")
    edits_b = [
        ('DE_DNMAG_MODULE_PATH = ROOT / "src/vcc2026/dn_mag.py"', 'DE_DNMAG_MODULE_PATH = ROOT / "src/vcc2026/dn_mag2.py"', 1),
        (f'DE_DNMAG_MODULE_SHA256 = "{V1}"', f'DE_DNMAG_MODULE_SHA256 = "{v2}"', 1),
        ("from vcc2026 import dn_mag as _dn_mag  # noqa: PLC0415", "from vcc2026 import dn_mag2 as _dn_mag  # noqa: PLC0415", 2),
        ("Module src/vcc2026/dn_mag.py, sha pinned.", "Module src/vcc2026/dn_mag2.py (693: v1's outputs byte for byte, faster), sha pinned.", 1),
    ]
    edits_d = [
        (f'"$VCC_DE_DNMAG" "{V1}"', f'"$VCC_DE_DNMAG" "{v2}"', 1),
        ("with dn_mag module 0a8cc560", f"with dn_mag module {v2[:8]} (693 v2)", 1),
        ("Builder >= 631 (dn_mag module 0a8cc560).", f"Builder >= 631 + 693p (dn_mag2 module {v2[:8]}).", 1),
    ]
    for old, new, n in edits_b:
        if b.count(old) != n:
            errs.append(f"builder anchor {old[:50]!r} found {b.count(old)} times (want {n})")
        b = b.replace(old, new)
    for old, new, n in edits_d:
        if d.count(old) != n:
            errs.append(f"def_day anchor {old[:50]!r} found {d.count(old)} times (want {n})")
        d = d.replace(old, new)
    if errs:
        print("REFUSED:\n  " + "\n  ".join(errs), file=sys.stderr)
        return 1
    if a.check:
        print(f"CHECK OK: v2 {v2[:16]} covered by 6 KA records; 5 builder + 3 def_day anchors; nothing written")
        return 0
    shutil.copy2(HERE / "dn_mag2.py", r / "src/vcc2026/dn_mag2.py")
    bu.write_text(b)
    dd.write_text(d)
    print(f"APPLIED 693p to {r}: builder -> dn_mag2 {v2[:16]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
