"""472: an absolute form of the prior gate, `VCC_DE_PRIOR_GATE=le<N>` (weak = design strength <= N), beside 468's share
form `q<NN>`.

Why (peer review 2026-09-26): strength is an integer count, so q33's cut jumps between classes. On H1 and A/B/C it lands
at "strength <= 1" (thresholds 1.72 and 2). On the efr synthetic October panel, 54-58 % of targets have strength <= 1, so
q33's threshold drops to 1 and the gate protects only strength 0. That leaves the strength-1 class, the source of all of
H1's reach loss, on the programme. `le1` names the rule the evidence supports: on A/B/C it selects exactly q33's weak
sets, so it must rebuild 458g byte for byte (chain472's KA).
q-mode is untouched: same parse, same threshold, same printed line, same config keys in the same order.
usage: uv run python experiments/472-gate-absolute/patch_builder472.py <builder.py> [--check]
"""
import sys
from pathlib import Path

path = Path(sys.argv[1])
src = path.read_text()
edits = [
    # 1. parse: q<NN> as before, or le<N>
    ('''    _pg = _re.fullmatch(r"q(\\d{1,2})", DE_PRIOR_GATE_ENV)
    if _pg is None or not 1 <= int(_pg.group(1)) <= 99:
        raise SystemExit(f"VCC_DE_PRIOR_GATE must be q<1..99>, e.g. q33: {DE_PRIOR_GATE_ENV!r}")
    if SAMELINE is None or not SAMELINE["prog"]:
        raise SystemExit("VCC_DE_PRIOR_GATE gates a programme default: it needs VCC_SAMELINE without noprog")
    DE_PRIOR_GATE = {"q": int(_pg.group(1)) / 100.0, "tau": 0.2}
''',
     '''    # 472: `le<N>` is the absolute form (weak = strength <= N); `q<NN>` the share form (468). A share cut on an integer
    # strength jumps between classes on a panel whose strength distribution differs from H1's / A/B/C's.
    _pg = _re.fullmatch(r"q(\\d{1,2})", DE_PRIOR_GATE_ENV)
    _pl = _re.fullmatch(r"le(\\d{1,2})", DE_PRIOR_GATE_ENV)
    if not ((_pg is not None and 1 <= int(_pg.group(1)) <= 99) or _pl is not None):
        raise SystemExit(f"VCC_DE_PRIOR_GATE must be q<1..99> (a share, e.g. q33) or le<0..99> (an absolute strength, "
                         f"e.g. le1): {DE_PRIOR_GATE_ENV!r}")
    if SAMELINE is None or not SAMELINE["prog"]:
        raise SystemExit("VCC_DE_PRIOR_GATE gates a programme default: it needs VCC_SAMELINE without noprog")
    DE_PRIOR_GATE = ({"q": int(_pg.group(1)) / 100.0, "tau": 0.2} if _pg is not None
                     else {"le": int(_pl.group(1)), "tau": 0.2})
'''),
    # 2. the threshold
    ('''                        _pg_thr = float(np.quantile(_pg_s, DE_PRIOR_GATE["q"]))
''',
     '''                        _pg_thr = (float(DE_PRIOR_GATE["le"]) + 1.0 if "le" in DE_PRIOR_GATE   # integer strength <= N
                                   else float(np.quantile(_pg_s, DE_PRIOR_GATE["q"])))
'''),
    # 3. the printed line (q-mode text unchanged)
    ('''                        print(f"[prior_gate] {context}: q {DE_PRIOR_GATE['q']:g} threshold {_pg_thr:g}; programme on "''',
     '''                        _pg_mode = (f"q {DE_PRIOR_GATE['q']:g}" if "q" in DE_PRIOR_GATE else f"le {DE_PRIOR_GATE['le']}")
                        print(f"[prior_gate] {context}: {_pg_mode} threshold {_pg_thr:g}; programme on "'''),
    # 4. config (q-mode: same keys, same order)
    ('''                        "q": DE_PRIOR_GATE["q"], "tau": DE_PRIOR_GATE["tau"], "weak": "clean4",''',
     '''                        **({"q": DE_PRIOR_GATE["q"]} if "q" in DE_PRIOR_GATE else {"le": DE_PRIOR_GATE["le"]}),
                        "tau": DE_PRIOR_GATE["tau"], "weak": "clean4",'''),
]
for old, new in edits:
    if src.count(old) != 1:
        raise SystemExit(f"anchor count {src.count(old)} != 1: {old[:90]!r}")
for old, new in edits:
    src = src.replace(old, new)
compile(src, str(path), "exec")
if "--check" in sys.argv:
    print("all anchors found; compiles; not writing")
else:
    path.write_text(src)
    print(f"patched {path}")
