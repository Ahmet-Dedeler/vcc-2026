"""590 — UNAPPLIED patch: a similarity floor on the day-of identity call before same-line rows switch on.

Why (590/RESULTS.md): on 402 labelled public profiles the matcher is right whenever the line is in CCLE (every wrong
unanimous call is a Cellosaurus-documented duplicate: PC-14 = PC-9, KCI-MOH1 = HPAC). A line ABSENT from CCLE is still
called unanimously: both HEK293 profiles (Orion HEK293T, Xu HEK293) -> HAP1_ENGINEERED, 8/8 methods, a registry line,
so today's rule would switch Jiang's HAP1 rows on. What separates them is the absolute match quality: the lowest top-1
correlation over the 8 methods is 0.51-0.61 for those false calls (plain and Flex-simulated), >= 0.68 for 99 % of correct
calls, and 0.719 / 0.733 / 0.772 for A/B/C (Flex, 18,400 cells, exactly how D/E/F will be measured).

Two edits, both inert on every known-answer row:
1. `scripts/def_identity.py`: `match()` also returns `min_top1_score` per context (min over the 8 methods of the top-1
   score it already computes) and prints it. No existing key or printed table changes; `--verify` is unchanged.
2. `scripts/def_sameline.py`: `--min-top1` (default 0.65; 0 disables). A registered, unanimous call whose
   `min_top1_score` is below it leaves the context unchanged, with the value in the provenance record. `--force` rows
   are untouched, an identity file without the key (pre-590) skips the gate with a WARNING, and a passing call writes
   nothing new (so a matched context's npz content is unchanged).

Order on the day it is applied: AFTER 538's applier (530 edits def_sameline.py; the anchors here are disjoint from
530's and match both the live and the post-530 file: `--check` proves it) and after 554.

    uv run python experiments/590-identity-validation/patch_identity590.py --check      # live + post-530, no writes
    uv run python experiments/590-identity-validation/patch_identity590.py --apply      # writes the two scripts
    uv run python experiments/590-identity-validation/patch_identity590.py --root <dir> --apply   # a copy
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT_DEFAULT = HERE.parents[1]
PATCH530 = ROOT_DEFAULT / "experiments/530-cloud-october-bundle/patch_bundle530.py"


def _apply(src: str, edits: list[tuple[str, str]], name: str) -> str:
    for old, new in edits:
        n = src.count(old)
        if n != 1:
            raise SystemExit(f"{name}: anchor found {n} times (need exactly 1):\n{old[:200]}")
        src = src.replace(old, new)
    return src


def patch_identity(src: str) -> str:
    if "min_top1_score" in src:
        raise SystemExit("def_identity.py: already carries min_top1_score (590 applied?)")
    edits = [
        (
            '        "call_unanimous_across_methods": {\n'
            '            ctx: len({results[k][ctx][0]["sample"] for k in results}) == 1 for ctx in contexts\n'
            '        },\n'
            '    }\n',
            '        "call_unanimous_across_methods": {\n'
            '            ctx: len({results[k][ctx][0]["sample"] for k in results}) == 1 for ctx in contexts\n'
            '        },\n'
            '        # 590: the lowest top-1 correlation over the 8 methods. A line absent from CCLE is still called\n'
            '        # unanimously (HEK293 -> HAP1, 8/8) at 0.51-0.61; correct calls sit >= 0.68 (q01 over 390 labelled\n'
            '        # profiles), A/B/C at 0.719-0.772. def_sameline.py --min-top1 gates the same-line mode on it.\n'
            '        "min_top1_score": {ctx: min(results[k][ctx][0]["score"] for k in results) for ctx in contexts},\n'
            '    }\n',
        ),
        (
            "        for row in consensus[ctx][:10]:\n"
            "            print(f\"  {row['rank']:2d} {row['sample']:<45} {str(row['lineage']):<20} {row['score']:.3f}\")\n"
            "    return out\n",
            "        for row in consensus[ctx][:10]:\n"
            "            print(f\"  {row['rank']:2d} {row['sample']:<45} {str(row['lineage']):<20} {row['score']:.3f}\")\n"
            "        print(f\"  lowest top-1 correlation over the 8 methods: {out['min_top1_score'][ctx]:.3f} \"\n"
            "              f\"(590: same-line rows need >= 0.65; A/B/C 0.719-0.772)\")\n"
            "    return out\n",
        ),
    ]
    return _apply(src, edits, "def_identity.py")


def patch_sameline(src: str) -> str:
    if "--min-top1" in src:
        raise SystemExit("def_sameline.py: already carries --min-top1 (590 applied?)")
    edits = [
        (
            '    ap.add_argument("--w-weak", type=float, default=0.5, help="w_<c> when any row-supplying bank is weak (450f)")\n',
            '    ap.add_argument("--w-weak", type=float, default=0.5, help="w_<c> when any row-supplying bank is weak (450f)")\n'
            '    ap.add_argument("--min-top1", type=float, default=0.65,\n'
            '                    help="590: a registered, unanimous call also needs the identity file\'s min_top1_score (the lowest "\n'
            '                         "top-1 correlation over the 8 methods) >= this. Lines absent from CCLE are called unanimously "\n'
            '                         "at 0.51-0.61 (HEK293 -> HAP1); A/B/C sit at 0.719-0.772. 0 disables")\n',
        ),
        (
            '            print(f"[sameline] {c}: {sample} is registered but the call is not unanimous; context unchanged", flush=True)\n'
            '            continue\n',
            '            print(f"[sameline] {c}: {sample} is registered but the call is not unanimous; context unchanged", flush=True)\n'
            '            continue\n'
            '        if banks and not a.force and a.min_top1 > 0:\n'
            '            s1 = ident.get("min_top1_score", {}).get(c)\n'
            '            if s1 is None:\n'
            '                print(f"[sameline] {c}: WARNING the identity file has no min_top1_score (pre-590 def_identity.py); "\n'
            '                      "the similarity gate is skipped", flush=True)\n'
            '            elif s1 < a.min_top1:\n'
            '                # 590: a line absent from CCLE is matched to its nearest CCLE line, unanimously; a registered nearest\n'
            '                # line would then inject another line\'s rows. Correct calls sit well above this floor.\n'
            '                rec["mode"] = (f"none (registered and unanimous, but the lowest top-1 correlation {s1:.3f} < "\n'
            '                               f"{a.min_top1}: the line may be absent from CCLE (590))")\n'
            '                prov["contexts"][c] = rec\n'
            '                print(f"[sameline] {c}: {sample} is registered and unanimous, but min top-1 correlation {s1:.3f} < "\n'
            '                      f"{a.min_top1}; context unchanged (590)", flush=True)\n'
            '                continue\n'
            '            else:\n'
            '                print(f"[sameline] {c}: min top-1 correlation {s1:.3f} >= {a.min_top1} (590 gate passed)", flush=True)\n',
        ),
    ]
    return _apply(src, edits, "def_sameline.py")


def _load_530():
    spec = importlib.util.spec_from_file_location("patch530", PATCH530)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=ROOT_DEFAULT)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--check", action="store_true")
    g.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    ident = a.root / "scripts/def_identity.py"
    same = a.root / "scripts/def_sameline.py"
    new_ident = patch_identity(ident.read_text())
    new_same = patch_sameline(same.read_text())
    for name, text in (("def_identity.py", new_ident), ("def_sameline.py", new_same)):
        compile(text, name, "exec")
    if a.check:
        print(f"[590] def_identity.py {hashlib.sha256(ident.read_bytes()).hexdigest()[:12]} -> "
              f"{hashlib.sha256(new_ident.encode()).hexdigest()[:12]} (anchors OK, compiles)")
        print(f"[590] def_sameline.py {hashlib.sha256(same.read_bytes()).hexdigest()[:12]} -> "
              f"{hashlib.sha256(new_same.encode()).hexdigest()[:12]} (anchors OK, compiles)")
        if "--programme-estimator" not in same.read_text():
            post530 = _load_530().patch_sameline(same.read_text())
            both = patch_sameline(post530)
            compile(both, "def_sameline.py (530 then 590)", "exec")
            print("[590] def_sameline.py after 530's patch_sameline: 590's anchors OK, compiles")
        print("[590] --check: nothing written")
        return
    ident.write_text(new_ident)
    same.write_text(new_same)
    print(f"[590] applied: def_identity.py {hashlib.sha256(ident.read_bytes()).hexdigest()}")
    print(f"[590] applied: def_sameline.py {hashlib.sha256(same.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
