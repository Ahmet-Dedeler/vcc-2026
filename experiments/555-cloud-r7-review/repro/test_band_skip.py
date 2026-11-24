"""band553's missing-effect path: a scored effect absent from a half is dropped, not raised.

Replays the loop at band553.py:50-52 against in-memory dicts. Does not read a half access.
"""
from __future__ import annotations


def kept(base_effects, halves, splits):
    out = []
    for k_ in base_effects:
        if not all(k_ in halves[(k, h)] for k in splits for h in "AB"):
            continue
        out.append(k_)
    return out


base = {"lcap-kept vs cpm": {"mean": 0.001}, "cpm-n1200 vs cpm": {"mean": -0.002}}
halves = {(k, h): {"lcap-kept vs cpm": {"mean": 0.0}} for k in (0, 1) for h in "AB"}
got = kept(base, halves, (0, 1))
assert got == ["lcap-kept vs cpm"], got
print("base effects", list(base))
print("present on every half", got)
print("silently dropped", [k for k in base if k not in got])
print('{"silent_skip": true}')
