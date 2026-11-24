#!/usr/bin/env python3
"""690h — the day-of DNMAG delivery report (coordinator 10-15, after DN read PARTIAL). UNAPPLIED; the coordinator applies it.

When VCC_DE_DNMAG is set, def_day.sh reads the final build's rc_audit dump (554's per-block record) and reports per context:
targets, mean dn_frac (delivered / needed DOWN magnitude), mean dn_moved (called DOWN genes moved per target), targets
under 200 moved genes, and the fail-safe counts (dn_invariant_failed, dn_guard_failed). A context is LOW if its mean dn_frac
< 0.10 or its mean dn_moved < 200 (the coordinator's PARTIAL condition; on A/B/C the DN build reads 0.19-0.22 and ~1,500).
REPORT only, never fatal: a LOW context is logged, written to <work>/dnmag_low.txt and to <work>/dnmag_audit.json, and is
repeated as a WARNING right before def_day's DONE line; it is NOT rebuilt (691's watcher copies the lines into NEXT.txt).
ONE greppable line always ends the report: `690h: dnmag LOW_CONTEXTS <ctx,...|none>`, or `-1` if the dump could not be
read or held no dn_ rows.

One file, scripts/def_day.sh, two anchored insertions (after 690a and 690g): the report before 690g's final-meta block,
the warning after 630dd's warning line.
    uv run python experiments/693-dnmag-speed/patch690h.py --root <repo> [--check]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

MARK = "690h: DNMAG delivery report"
R_OLD = '''if [ -n "${BLC:-}" ]; then   # 690g: dump hook -- the final meta must carry the per-context bulk it was built with
'''
R_NEW = '''if [ -n "${VCC_DE_DNMAG:-}" ]; then   # 690h: DNMAG delivery report -- per context from the final rc_audit dump; REPORT, never fatal
  python3 - "$FINAL.rc_audit.jsonl" "$W/dnmag_audit.json" "$W/dnmag_low.txt" <<'PY' 2>&1 | while IFS= read -r _l; do log "690h: $_l"; done || true
import collections, json, sys
src, out_json, out_low = sys.argv[1:4]
try:
    rows = [json.loads(l) for l in open(src) if l.strip()]
    rows = [r for r in rows if "dn_frac" in r and "dn_moved" in r]
except Exception as exc:
    print(f"dnmag: could not read {src}: {exc}")
    rows = []
if not rows:
    print("dnmag LOW_CONTEXTS -1")
    sys.exit(0)
by = collections.defaultdict(list)
for r in rows:
    by[str(r.get("context"))].append(r)
rep, low = {}, []
for c in sorted(by):
    rs = by[c]
    fr = sum(float(r["dn_frac"]) for r in rs) / len(rs)
    mv = sum(float(r["dn_moved"]) for r in rs) / len(rs)
    rec = {"targets": len(rs), "mean_dn_frac": round(fr, 4), "mean_dn_moved": round(mv, 1),
           "targets_moved_lt_200": sum(1 for r in rs if float(r["dn_moved"]) < 200),
           "invariant_failed": int(sum(float(r.get("dn_invariant_failed", 0)) for r in rs)),
           "guard_failed": int(sum(float(r.get("dn_guard_failed", 0)) for r in rs)),
           "low": bool(fr < 0.10 or mv < 200)}
    rep[c] = rec
    if rec["low"]:
        low.append(c)
    print(f"dnmag {c}: " + json.dumps(rec))
json.dump(rep, open(out_json, "w"), indent=1)
if low:
    open(out_low, "w").write("\\n".join(low) + "\\n")
print("dnmag LOW_CONTEXTS " + (",".join(low) if low else "none"))
PY
fi
''' + R_OLD
W_ANCHOR = '''[ -s "$W/uprn_effectively_off.txt" ] && log "630dd: WARNING 630 EFFECTIVELY OFF'''
W_LINE = '''[ -s "$W/dnmag_low.txt" ] && log "690h: WARNING DNMAG under-delivered (mean dn_frac < 0.10 or < 200 genes/target) in: $(tr '\\n' ' ' < "$W/dnmag_low.txt")"   # 690h: DNMAG delivery report
'''


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    p = a.root.resolve() / "scripts/def_day.sh"
    t = p.read_text()
    if MARK in t:
        print("REFUSED: 690h is already applied", file=sys.stderr)
        return 1
    if "690a: DOWN magnitude pass-through" not in t or "690g: dump hook" not in t:
        print("REFUSED: def_day.sh lacks 690a / 690g", file=sys.stderr)
        return 1
    if t.count(R_OLD) != 1:
        print(f"REFUSED: report anchor found {t.count(R_OLD)} times", file=sys.stderr)
        return 1
    lines = t.splitlines(keepends=True)
    hits = [i for i, ln in enumerate(lines) if ln.startswith(W_ANCHOR)]
    if len(hits) != 1:
        print(f"REFUSED: warning anchor found {len(hits)} times", file=sys.stderr)
        return 1
    t2 = "".join(lines[:hits[0] + 1]) + W_LINE + "".join(lines[hits[0] + 1:])
    t2 = t2.replace(R_OLD, R_NEW, 1)
    if a.check:
        print("CHECK OK: 2 anchors; nothing written")
        return 0
    p.write_text(t2)
    print(f"APPLIED 690h to {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
