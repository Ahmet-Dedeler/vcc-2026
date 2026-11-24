"""252 — diff the BLIND (pre-252) and FIRE (guarded) probe runs into one verdict table.

Rules, per probe key:
  * a key tagged [PASS] (a correct input) must not raise in either run and every recorded
    number / digest must be IDENTICAL across the two runs (bit-identity of the default path);
  * a key tagged [FIRE] / [SWAP] / [TAUTOLOGY] / [TRAP] / "must miss" is a broken input: the
    guarded run must raise a ValueError subclass (HarnessGuardError) where the blind run
    returned a value -- or, where the blind harness could not even express the check
    (`name=` did not exist: TypeError), the guarded run must raise the guard;
  * everything else is reported as-is.

    uv run python compare252.py        -> work/compare252.json, work/compare252.md ; exit 1 on FAIL
"""
from __future__ import annotations

import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
BROKEN_TAGS = ("[FIRE", "[SWAP", "[TAUTOLOGY", "[TRAP", "must miss", "[empty expectation")


def flat(v, prefix=""):
    if isinstance(v, dict):
        out = {}
        for k, x in v.items():
            out.update(flat(x, f"{prefix}{k}."))
        return out
    if isinstance(v, str) and "/" in v and prefix.rstrip(".").endswith("_op"):
        v = v.rsplit("/", 1)[-1]          # a path: the blind run resolves it under its faux root
    return {prefix.rstrip("."): v}


def same(a, b) -> bool:
    if isinstance(a, float) and isinstance(b, float):
        return (math.isnan(a) and math.isnan(b)) or a == b
    return a == b


def main() -> None:
    blind = json.loads((WORK / "probes-blind.json").read_text())
    fire = json.loads((WORK / "probes-fire.json").read_text())
    rows, fails = [], 0
    for section in ("identity", "e064", "axis_trap", "preserve_endpoints_zero", "cache"):
        b_sec, f_sec = blind[section], fire[section]
        for key in f_sec:
            if key.startswith("_"):
                continue
            b, f = b_sec.get(key), f_sec[key]
            if not (isinstance(f, dict) and "raised" in f):
                rows.append(dict(section=section, probe=key, verdict="info",
                                 blind=str(b)[:160], fire=str(f)[:160]))
                continue
            broken = any(t in key for t in BROKEN_TAGS)
            is_pass = "[PASS" in key
            b_raised = b.get("raised") if isinstance(b, dict) else "<absent in blind>"
            f_raised = f.get("raised")
            if "endpoints(all-zero truth)" in key:
                # designed to return NaN (not raise): blind crashed (ZeroDivisionError) or read 0.0
                v = f.get("value") or {}
                ok = f_raised is None and all(isinstance(v.get(k), float) and math.isnan(v[k])
                                              for k in ("sgn", "sgn_w"))
                verdict = "FIRE" if ok else "FAIL"
                detail = f"blind: {b_raised or 'silent 0.0'}; fire: sgn/sgn_w = NaN (mean_endpoints skips them)"
            elif is_pass:
                if b is None:
                    ok = f_raised is None
                    detail = "fire-only PASS row"
                elif b_raised == "TypeError" and f_raised is None:
                    ok = True
                    detail = "blind: `name=` not expressible (TypeError); fire: accepted"
                else:
                    fb, ff = flat(b.get("value", {})), flat(f.get("value", {}))
                    diff = {k: (fb.get(k), ff.get(k)) for k in set(fb) | set(ff)
                            if k not in ("name", "sha256", "axis_sha256")
                            and not k.startswith("audit_cache")      # 252's new provenance keys
                            and not same(fb.get(k), ff.get(k))}
                    ok = b_raised is None and f_raised is None and not diff
                    detail = "identical" if not diff else f"DIFFER {diff}"
                verdict = "PASS" if ok else "FAIL"
            elif broken:
                fired = f_raised is not None and (f.get("is_value_error") or f_raised == "RebuildAttempted")
                silent_before = (b is None) or (b_raised is None) or (b_raised == "TypeError")
                already = isinstance(b, dict) and b.get("is_value_error")   # the OLD check fired too
                ok = fired and (silent_before or already)
                verdict = ("FIRE" if silent_before else "FIRE (pre-252 too)") if ok else "FAIL"
                detail = (f"blind: {'silent' if b_raised is None else b_raised}; "
                          f"fire: {f_raised}: {str(f.get('message', ''))[:110]}")
            else:
                verdict, ok, detail = "info", True, ""
            fails += int(not ok)
            rows.append(dict(section=section, probe=key, verdict=verdict, detail=detail,
                             blind_raised=b_raised, fire_raised=f_raised,
                             blind_seconds=(b or {}).get("seconds"), fire_seconds=f.get("seconds")))
    md = ["| section | probe | blind (pre-252) | fire (252) | verdict |", "|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['section']} | `{r['probe']}` | {r.get('blind_raised') or 'silent'} | "
                  f"{r.get('fire_raised') or 'accepted'} | **{r['verdict']}** {r.get('detail', '')[:120]} |")
    (WORK / "compare252.md").write_text("\n".join(md) + "\n")
    (WORK / "compare252.json").write_text(json.dumps(dict(rows=rows, fails=fails,
                                                          blind_plolo=blind["plolo_sha256"],
                                                          fire_plolo=fire["plolo_sha256"],
                                                          blind_e064=blind["e064_sha256"],
                                                          fire_e064=fire["e064_sha256"]), indent=1))
    print("\n".join(md))
    print(f"\n{len(rows)} rows, {fails} failures -> {'PASS' if fails == 0 else 'FAIL'}")
    if fails:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
