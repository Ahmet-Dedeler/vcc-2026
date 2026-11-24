#!/usr/bin/env python3
"""Move one SHIPPED_ENV default in scripts/def_pipeline.py after a board read adopts a rider (455's pattern, generalised).

    uv run python scripts/shipped_env_update.py VCC_DE_ZERO_REPAIR=-2 "zero-repair floor -2: +0.00xx (446br2, <entry>)" [--file F] [--check]

- An existing key: its value literal is replaced; the key must occur exactly once in SHIPPED_ENV's literal.
- A new key: it is appended to SHIPPED_ENV's dict literal.
- The ARM_BOARD_DELTA text for the key is set (added if absent), so an arm-drift message can always name it.
- The patched module must import and expose the new value; otherwise the file is left unchanged.
After it, re-prove row 1: the driver's defaults plus the adopted arm's two day-of exports must rebuild that arm byte for byte
(experiments/455-driver-q175/chain455.sh is the template)."""
import argparse
import ast
import importlib.util
import re
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("assign"); ap.add_argument("delta_text")
ap.add_argument("--file", default=str(ROOT / "scripts/def_pipeline.py")); ap.add_argument("--check", action="store_true")
a = ap.parse_args()
key, value = a.assign.split("=", 1)
p = Path(a.file); s = p.read_text()
m = re.search(r"SHIPPED_ENV = \{(.*?)\}\n", s, re.S)
if not m:
    raise SystemExit("SHIPPED_ENV literal not found")
body = m.group(1)
pat = re.compile(r'"%s": "([^"]*)"' % re.escape(key))
hits = pat.findall(body)
if len(hits) > 1:
    raise SystemExit(f"{key} occurs {len(hits)} times in SHIPPED_ENV")
if hits:
    new_body = pat.sub(f'"{key}": "{value}"', body)
else:
    new_body = body.rstrip() + f',\n               "{key}": "{value}"'
s2 = s[:m.start(1)] + new_body + s[m.end(1):]
d = re.search(r"ARM_BOARD_DELTA = \{", s2)
if not d:
    raise SystemExit("ARM_BOARD_DELTA not found")
# the key's value: one or more adjacent string literals (implicit concatenation across lines), then "," or "}"
dpat = re.compile(r'"%s":\s*((?:"(?:[^"\\]|\\.)*"\s*)+)(?=[,}])' % re.escape(key))
lit = '"' + a.delta_text.replace("\\", "\\\\").replace('"', '\\"') + '"'
if dpat.search(s2[d.end():]):
    head, tail = s2[:d.end()], s2[d.end():]
    s2 = head + dpat.sub(lambda mm: f'"{key}": {lit}', tail, count=1)
else:
    s2 = s2[:d.end()] + f'"{key}": {lit},\n' + "                   " + s2[d.end():]
with tempfile.TemporaryDirectory() as td:
    tp = Path(td) / "def_pipeline_patched.py"; tp.write_text(s2)
    spec = importlib.util.spec_from_file_location("dpp", tp); mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "scripts"))
    spec.loader.exec_module(mod)
    assert mod.SHIPPED_ENV[key] == value, "patched module does not expose the new value"
    assert mod.ARM_BOARD_DELTA.get(key) == a.delta_text, "ARM_BOARD_DELTA does not carry the new text"
print(f"{key}: {hits[0] if hits else '<absent>'} -> {value}; module imports and exposes it")
if not a.check:
    p.write_text(s2); print("patched", p)
