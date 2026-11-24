"""533: insert the Round 14 addendum into CLAIMS.md (from the untouched v14 copy), before the historical summary table."""
import pathlib, sys
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from v15_claims import ADDENDUM  # noqa: E402
src = (HERE / "CLAIMS.v14.orig.md").read_text()
anchor = "## Summary of confirmation status (historical through Round 4)"
assert src.count(anchor) == 1
out = src.replace(anchor, ADDENDUM + anchor)
(HERE.parent / "paper" / "CLAIMS.md").write_text(out)
print(f"CLAIMS v14 {len(src):,} -> v15 {len(out):,} chars; rows C-74..C-93: {sum(out.count(f'**C-{i}**') for i in range(74, 94))}")
