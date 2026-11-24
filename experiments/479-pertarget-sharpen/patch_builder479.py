"""479: per-target bulk sharpening, `VCC_BULK_SHARPEN_Q=tercile:<low>,<mid>,<high>` or `stercile:<low>,<mid>,<high>`, beside
the scalar `VCC_BULK_SHARPEN_Q=<q>` (unchanged). Third session's 492/493 (492/RESULTS.md "Knob definition"):
  tercile  : the key tercile. r_t = ||bulk_t|| / ||carriers_t|| over NON-PANEL columns of the LOADED stack dump
             (= c * w_t for the key the dump was built with; exact to 1.6e-14 on 446bj2's LODO dump). Edges
             np.quantile(r, [1/3, 2/3]), np.digitize. Floor guard: targets at min(r) (within 1e-9 relative) -> LOW.
             Tie block: >= 2 targets with identical r (1e-9 relative) that contain a cut -> MID (non-carried targets
             sit at w = 1).
  stercile : the strength tercile, s_t = ||carriers_t|| over non-panel columns (the unkeyed stack row). Continuous.
LOW (weak / unreliable rows) gets the first q, MID the second, HIGH the third. Rows keep their own gated non-panel norm
(bulk_sharpen), panel columns are untouched, and the carriers (DE channel) are untouched, exactly as the scalar knob.
Needs VCC_BULK_LOAD (the tercile is read from the dump) and VCC_PERT_SCALE_BALLAST. Unset or scalar: no code path changes.
usage: uv run python experiments/479-pertarget-sharpen/patch_builder479.py <builder.py> [--check]
"""
import sys
from pathlib import Path

path = Path(sys.argv[1]); src = path.read_text()
edits = [
    # 1. parse: the per-target forms, before the scalar float parse
    ('''if BULK_SHARPEN_Q is not None:
    BULK_SHARPEN_Q = float(BULK_SHARPEN_Q)
    if not 1.0 < BULK_SHARPEN_Q <= 2.0:
        raise SystemExit(f"VCC_BULK_SHARPEN_Q must be in (1, 2] when set: {BULK_SHARPEN_Q!r}")
''',
     '''# 479 (2026-09-27): `tercile:<low>,<mid>,<high>` / `stercile:<low>,<mid>,<high>` sharpen each target's bulk row with
# the q of its tercile (492/493, third session): weak or unreliable rows gain most from a sharper q. The tercile is read
# from the LOADED stack dump (sharpen_q_by_tercile). The scalar form below is unchanged.
BULK_SHARPEN_PT: dict[str, object] | None = None
BULK_SHARPEN_PT_STATE: dict[str, object] = {}
if BULK_SHARPEN_Q is not None and BULK_SHARPEN_Q.split(":", 1)[0].strip().lower() in ("tercile", "stercile"):
    _pt_mode, _pt_q = BULK_SHARPEN_Q.split(":", 1)
    try:
        _pt_qv = [float(x) for x in _pt_q.split(",")]
    except ValueError:
        raise SystemExit(f"VCC_BULK_SHARPEN_Q per-target form must be tercile:<q>,<q>,<q> or stercile:...: {BULK_SHARPEN_Q!r}")
    if len(_pt_qv) != 3 or not all(1.0 < x <= 3.5 for x in _pt_qv):
        raise SystemExit(f"VCC_BULK_SHARPEN_Q per-target qs must be three values in (1, 3.5]: {BULK_SHARPEN_Q!r}")
    if PERT_SCALE_BALLAST is None:
        raise SystemExit("VCC_BULK_SHARPEN_Q is a bulk-only arm and needs VCC_PERT_SCALE_BALLAST (welded it reaches the DE channel)")
    BULK_SHARPEN_PT = {"mode": _pt_mode.strip().lower(), "q": _pt_qv}
    BULK_SHARPEN_Q = None
if BULK_SHARPEN_Q is not None:
    BULK_SHARPEN_Q = float(BULK_SHARPEN_Q)
    if not 1.0 < BULK_SHARPEN_Q <= 2.0:
        raise SystemExit(f"VCC_BULK_SHARPEN_Q must be in (1, 2] when set: {BULK_SHARPEN_Q!r}")
'''),
    # 2. the tercile function, after bulk_sharpen
    ('''    if not np.isfinite(n1) or n1 <= 0.0:
        return h
    return np.where(on, s * (n0 / n1), h)
''',
     '''    if not np.isfinite(n1) or n1 <= 0.0:
        return h
    return np.where(on, s * (n0 / n1), h)


def sharpen_q_by_tercile(bulk, carriers, genes, perts, pt) -> tuple[np.ndarray, dict[str, object]]:
    """479: a per-target sharpening q from the loaded stack dump's arrays (see VCC_BULK_SHARPEN_Q's per-target forms)."""
    panel = set(map(str, perts))
    nonpanel = np.fromiter((str(g) not in panel for g in genes), dtype=bool, count=len(genes))
    cn = np.linalg.norm(np.asarray(carriers, dtype=np.float64)[:, nonpanel], axis=1)
    floor = np.zeros(len(cn), dtype=bool)
    if pt["mode"] == "stercile":
        v = cn
    else:
        bn = np.linalg.norm(np.asarray(bulk, dtype=np.float64)[:, nonpanel], axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.where(cn > 0, bn / np.where(cn > 0, cn, 1.0), np.nan)
        if not np.isfinite(v).any():
            raise SystemExit("VCC_BULK_SHARPEN_Q tercile: no finite key ratio in the loaded dump")
        vmin = float(np.nanmin(v))
        floor = np.isfinite(v) & (v <= vmin * (1.0 + 1e-9))
    fill = float(np.nanmedian(v))
    vv = np.where(np.isfinite(v), v, fill)
    edges = np.quantile(vv, [1.0 / 3.0, 2.0 / 3.0])
    t = np.digitize(vv, edges)
    tie_mid = np.zeros(len(vv), dtype=bool)
    if pt["mode"] == "tercile":
        for e in edges:   # a block of identical ratios containing a cut goes to MID (non-carried targets sit at w = 1)
            blk = np.abs(vv - e) <= abs(e) * 1e-9
            if int(blk.sum()) >= 2:
                tie_mid |= blk
        t = np.where(tie_mid, 1, t)
        t = np.where(floor, 0, t)
    q = np.asarray(pt["q"], dtype=np.float64)[t]
    audit = {"mode": pt["mode"], "q": list(pt["q"]), "edges": [float(e) for e in edges],
             "counts": [int((t == i).sum()) for i in range(3)], "floor_to_low": int(floor.sum()),
             "tie_to_mid": int(tie_mid.sum()), "nonfinite": int((~np.isfinite(v)).sum())}
    return q, audit
'''),
    # 3. read the carriers with the loaded bulk
    ('''            loaded = np.asarray(z["bulk"], dtype=np.float64)
            dump_env = json.loads(str(z["env"]))
''',
     '''            loaded = np.asarray(z["bulk"], dtype=np.float64)
            dump_env = json.loads(str(z["env"]))
            _pt_carriers = np.asarray(z["carriers"], dtype=np.float64) if BULK_SHARPEN_PT is not None else None
'''),
    # 4. compute the per-target q once the bulk is swapped in
    ('''        effects = loaded
''',
     '''        effects = loaded
        if BULK_SHARPEN_PT is not None:
            _qt, _pt_audit = sharpen_q_by_tercile(loaded, _pt_carriers, spec.genes, spec.perts, BULK_SHARPEN_PT)
            BULK_SHARPEN_PT_STATE.update({"qt": _qt, "audit": _pt_audit})
            print(f"[bulk-sharpen-pt] {_pt_audit}", flush=True)
'''),
    # 5. the render call
    ('''                                sharpen_q=BULK_SHARPEN_Q,
''',
     '''                                sharpen_q=(BULK_SHARPEN_Q if BULK_SHARPEN_PT is None
                                           else float(BULK_SHARPEN_PT_STATE["qt"][target_index])),
'''),
    # 6. the guard before emission: per-target needs the loaded dump and the ballast path
    ('''    if BULK_SHARPEN_Q is not None:
        if bulk_effects is None:
            raise ValueError("VCC_BULK_SHARPEN_Q needs the ballast path (bulk_effects)")
''',
     '''    if BULK_SHARPEN_Q is not None:
        if bulk_effects is None:
            raise ValueError("VCC_BULK_SHARPEN_Q needs the ballast path (bulk_effects)")
    if BULK_SHARPEN_PT is not None:
        if bulk_effects is None or "qt" not in BULK_SHARPEN_PT_STATE:
            raise ValueError("VCC_BULK_SHARPEN_Q per-target needs VCC_BULK_LOAD and the ballast path")
'''),
    # 7. config
    ('''            **({} if BULK_SHARPEN_Q is None else {"bulk_sharpen_q": BULK_SHARPEN_Q}),
''',
     '''            **({} if BULK_SHARPEN_Q is None else {"bulk_sharpen_q": BULK_SHARPEN_Q}),
            **({} if BULK_SHARPEN_PT is None else {"bulk_sharpen_pt": dict(BULK_SHARPEN_PT_STATE.get("audit", BULK_SHARPEN_PT))}),
'''),
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
    path.write_text(src); print(f"patched {path}")
