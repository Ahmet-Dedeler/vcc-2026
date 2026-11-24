"""420: assemble RESULTS.md from results_head.md + results_tail.md, filling {{T-xxx}} placeholders with the sections of
work/table420.md (rendered from the analysis JSON by table420.py), so no table number is transcribed by hand.

    uv run python assemble420.py
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
tab = (HERE / "work/table420.md").read_text()
sections = {}
for block in re.split(r"(?m)^(?=### T-)", tab):
    m = re.match(r"### (T-[a-z0-9]+)\.", block)
    if m:
        sections[m.group(1)] = block.strip()
text = (HERE / "results_head.md").read_text() + "\n" + (HERE / "results_tail.md").read_text()
missing = [k for k in re.findall(r"\{\{(T-[a-z0-9]+)\}\}", text) if k not in sections]
assert not missing, f"missing table sections: {missing}"
text = re.sub(r"\{\{(T-[a-z0-9]+)\}\}", lambda m: sections[m.group(1)], text)
(HERE / "RESULTS.md").write_text(text)
print(f"RESULTS.md: {len(text.splitlines())} lines; sections used: {sorted(set(re.findall(r'### (T-[a-z0-9]+)', text)))}")
