"""319 -- the built artifact's config against the object it is ONE variable away from.

  uv run python experiments/319-cap-variant/configdiff319.py <built.h5ad.meta.json>

Exits non-zero unless the built config differs from `submissions/312-normrel15.h5ad.meta.json`
(entry `Q519I2pLWRCW7PAEx5ul`, the shipped `VCC_NORMREL_EXP=-1.5` arm) in EXACTLY the one key
`normrel_mode`.  AGENTS rule 10: diff the finished artifact's config, never the command you typed.
The diff against `249-cd4svd5-f034` (the incumbent the postbuild check uses) is printed for the
record but does not gate.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REFS = {
    "312-normrel15": ROOT / "submissions/312-normrel15.h5ad.meta.json",
    "249-cd4svd5-f034": ROOT / "submissions/249-cd4svd5-f034.h5ad.meta.json",
}
EXPECTED = {"normrel_mode"}


def diff(a: dict, b: dict) -> dict[str, tuple[object, object]]:
    out: dict[str, tuple[object, object]] = {}
    for k in sorted(set(a) | set(b)):
        va, vb = a.get(k, "<absent>"), b.get(k, "<absent>")
        if va != vb:
            out[k] = (va, vb)
    return out


def main() -> int:
    built = json.loads(Path(sys.argv[1]).read_text())["config"]
    rc = 0
    for label, path in REFS.items():
        if not path.exists():
            print(f"[319] config vs {label}: MISSING {path}")
            if label == "312-normrel15":
                rc = 2
            continue
        ref = json.loads(path.read_text())["config"]
        d = diff(ref, built)
        print(f"[319] config vs {label}: {len(d)} key(s) differ")
        for k, (va, vb) in d.items():
            print(f"[319]     {k}: {va!r} -> {vb!r}")
        if label == "312-normrel15":
            if set(d) != EXPECTED:
                print(f"[319] FAIL: expected exactly {sorted(EXPECTED)}, got {sorted(d)}")
                rc = 1
            elif d["normrel_mode"] != ("<absent>", "cap"):
                print(f"[319] FAIL: normrel_mode moved {d['normrel_mode']}, expected "
                      f"('<absent>', 'cap')")
                rc = 1
            else:
                print("[319] PASS: exactly one key, `normrel_mode`: <absent> -> 'cap'")
    return rc


if __name__ == "__main__":
    sys.exit(main())
