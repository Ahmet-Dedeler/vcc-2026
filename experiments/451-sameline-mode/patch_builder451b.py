"""451b: a PER-CONTEXT same-line blend weight. `VCC_SAMELINE=<npz>:<sha>:auto[:pt][:noprog]` reads each matched context's
w from the artifact (`w_<c>`, written by `scripts/def_sameline.py` from the bank's quality class). A numeric w keeps today's
global behaviour, so every existing call and both known answers are unchanged.

Why (450f, H1): the best w falls with the bank's sign quality: 2 for a Pan-quality bank (+0.093), 0.5 for a 0.61 bank,
where w 2 is neutral overall and costs mse. Two October contexts can match banks of different quality.

usage: uv run python experiments/451-sameline-mode/patch_builder451b.py <builder.py> [--check]   (applies on top of 451)
"""
import sys
from pathlib import Path

path = Path(sys.argv[1])
src = path.read_text()
edits = [
    ('''    SAMELINE = {"path": Path(_sl[0]) if Path(_sl[0]).is_absolute() else ROOT / _sl[0], "sha256": _sl[1],
                "w": float(_sl[2]), "pt": "pt" in _flags, "prog": "noprog" not in _flags}
    if not 0.0 <= float(SAMELINE["w"]) <= 10.0:
        raise SystemExit(f"VCC_SAMELINE blend weight must be in [0, 10]: {SAMELINE['w']!r}")''',
     '''    SAMELINE = {"path": Path(_sl[0]) if Path(_sl[0]).is_absolute() else ROOT / _sl[0], "sha256": _sl[1],
                # 451b: "auto" = each context's w from the artifact (`w_<c>`, by bank quality); a number = one global w
                "w": None if _sl[2] == "auto" else float(_sl[2]), "pt": "pt" in _flags, "prog": "noprog" not in _flags}
    if SAMELINE["w"] is not None and not 0.0 <= float(SAMELINE["w"]) <= 10.0:
        raise SystemExit(f"VCC_SAMELINE blend weight must be in [0, 10] or auto: {SAMELINE['w']!r}")'''),
    ('''            rec["programme_ok"] = np.asarray(z[f"programme_ok_{c}"]) if f"programme_ok_{c}" in z.files else np.asarray(True)''',
     '''            rec["programme_ok"] = np.asarray(z[f"programme_ok_{c}"]) if f"programme_ok_{c}" in z.files else np.asarray(True)
            if f"w_{c}" in z.files:
                rec["w"] = np.asarray(z[f"w_{c}"])
                if not 0.0 <= float(rec["w"]) <= 10.0:
                    raise SystemExit(f"VCC_SAMELINE artifact context {c}: w_{c} {float(rec['w'])!r} outside [0, 10]")
            elif SAMELINE["w"] is None:
                raise SystemExit(f"VCC_SAMELINE w is auto but the artifact carries no w_{c} (rebuild it with def_sameline.py)")'''),
    ('''                    _w = float(SAMELINE["w"])''',
     '''                    _w = float(SAMELINE["w"]) if SAMELINE["w"] is not None else float(_slr["w"])'''),
    ('''"w": SAMELINE["w"], "pt": SAMELINE["pt"], "prog": SAMELINE["prog"],''',
     '''"w": "auto" if SAMELINE["w"] is None else SAMELINE["w"], "pt": SAMELINE["pt"], "prog": SAMELINE["prog"],'''),
    ('''# `VCC_SAMELINE="<npz>:<sha256>:<w>[:pt][:noprog]"` -- 451''',
     '''# `VCC_SAMELINE="<npz>:<sha256>:<w|auto>[:pt][:noprog]"` -- 451 (451b: `auto` = per-context w from the artifact)'''),
]
for old, new in edits:
    if src.count(old) != 1:
        raise SystemExit(f"anchor count {src.count(old)} != 1: {old[:80]!r}")
for old, new in edits:
    src = src.replace(old, new)
if "--check" in sys.argv:
    print("all anchors found; not writing")
else:
    path.write_text(src)
    print(f"patched {path}")
