"""533: build Draft v15's PAPER.md from the untouched v14 copy and the text blocks in v15_blocks_{a,b,c,d}.py.

Idempotent: always starts from work/PAPER.v14.orig.md (a byte copy of experiments/294-writeup/PAPER.md at f31ad98a+) and
asserts that every anchor occurs exactly once before it is used, so a changed source fails loudly instead of misplacing
text. Run: uv run python experiments/533-science-lead/work/apply_v15.py
"""
import hashlib
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from v15_blocks_a import ANCHORS_A  # noqa: E402
from v15_blocks_b import ANCHORS_B  # noqa: E402
from v15_blocks_c import ANCHORS_C  # noqa: E402
from v15_blocks_d import ANCHORS_D  # noqa: E402

SRC = HERE / "PAPER.v14.orig.md"
DST = HERE.parent / "paper" / "PAPER.md"


def apply(text: str, anchors) -> str:
    for anchor, mode, block in anchors:
        n = text.count(anchor)
        if n != 1:
            raise SystemExit(f"anchor occurs {n}x (need 1): {anchor[:90]!r}")
        if mode == "after":
            text = text.replace(anchor, anchor + block)
        elif mode == "before":
            text = text.replace(anchor, block + anchor)
        elif mode == "replace":
            text = text.replace(anchor, block)
        else:
            raise SystemExit(f"bad mode {mode}")
    return text


def main() -> None:
    src = SRC.read_text()
    out = src
    for group in (ANCHORS_A, ANCHORS_B, ANCHORS_C, ANCHORS_D):
        out = apply(out, group)
    DST.write_text(out)
    print(f"v14 {len(src):,} chars sha {hashlib.sha256(src.encode()).hexdigest()[:12]} -> v15 {len(out):,} chars "
          f"sha {hashlib.sha256(out.encode()).hexdigest()[:12]}; pending markers {out.count('[[PENDING 09-30')}")


if __name__ == "__main__":
    main()
