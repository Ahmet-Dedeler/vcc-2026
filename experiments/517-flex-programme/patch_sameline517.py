"""517: `scripts/def_sameline.py --programme-estimator {pdown,meanz_ne}` (default pdown = byte-identical behaviour). NOT
APPLIED. With meanz_ne, the --programme-source programme's pdown_<c> is 517's hook-source encoding (artifact517.py): the mean
over the source's NON-ESSENTIAL rows (K562-essential list out; the panel's own rows out as before) of delta / per-gene sd,
mapped to a pseudo-P whose tau-0.05 coverage equals the P_down rule's; programme_<c> is its sign. It is meant for the hook's
own artifact (def_day's `hookbank:` form, DEF_DAY_HOOK_ESTIMATOR=meanz_ne), not the block prior. KA: artifact517.py's
meanz_ne artifact is reproduced (arrays equal) from the A/B/C bundle with --programme-source h1_2025_full.
usage: uv run python experiments/517-flex-programme/patch_sameline517.py <def_sameline.py> [--check]
"""
import os
import sys
from pathlib import Path

p = Path(sys.argv[1]); s = p.read_text()
if "--programme-estimator" in s:
    raise SystemExit("already patched")
E = [
    ('def programme(delta, perts, genes, official, panel) -> tuple[np.ndarray, np.ndarray, int, np.ndarray]:\n',
     'ESSENTIAL_LIST = ROOT / "data/donors/k562_essential/perts.txt"   # 517: meanz_ne drops these rows\n\n\n'
     'def programme_meanz_ne(delta, perts, genes, official, panel) -> tuple[np.ndarray, np.ndarray, int]:\n'
     '    """517 (artifact517.py): mean over the non-essential, non-panel rows of delta / sd (sd over all non-panel rows),\n'
     '    as a pseudo-P whose |P - 0.5| >= 0.05 coverage equals the P_down rule\'s (so the hook\'s tau keeps its meaning)."""\n'
     '    ess = {l.strip() for l in open(ESSENTIAL_LIST) if l.strip()}\n'
     '    keep = np.array([p not in panel for p in perts])\n'
     '    X = np.asarray(delta[np.flatnonzero(keep)], dtype=np.float64)\n'
     '    ne = np.array([p not in ess for p in np.asarray(perts)[keep]])\n'
     '    nz = (X != 0).sum(0); dn = (X < 0).sum(0)\n'
     '    with np.errstate(invalid="ignore", divide="ignore"):\n'
     '        pd_all = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)\n'
     '    pos = {g: j for j, g in enumerate(genes)}\n'
     '    col = np.array([pos.get(g, -1) for g in official])\n'
     '    P_all = np.where(col >= 0, pd_all[np.maximum(col, 0)], np.nan)\n'
     '    sd = X.std(0); sd[sd == 0] = 1.0\n'
     '    zb = (X[ne] / sd).mean(0)\n'
     '    z = np.where(np.isfinite(P_all) & (col >= 0), zb[np.maximum(col, 0)], np.nan)\n'
     '    n_cov = int((np.isfinite(P_all) & (np.abs(P_all - 0.5) >= 0.05)).sum())\n'
     '    z_q = float(np.sort(np.abs(np.nan_to_num(z)))[::-1][n_cov - 1])\n'
     '    pdown = np.where(np.isfinite(z), np.clip(0.5 - 0.05 * z / z_q, 0.001, 0.999), np.nan)\n'
     '    psign = np.where(np.isfinite(pdown) & (pdown <= 0.5), 1.0, -1.0)\n'
     '    return psign, pdown, int(ne.sum())\n\n\n'
     'def programme(delta, perts, genes, official, panel) -> tuple[np.ndarray, np.ndarray, int, np.ndarray]:\n'),
    ('    ap.add_argument("--programme-source", default="",\n',
     '    ap.add_argument("--programme-estimator", choices=["pdown", "meanz_ne"], default="pdown",\n'
     '                    help="517: meanz_ne = the hook-source encoding for --programme-source (not for the block prior)")\n'
     '    ap.add_argument("--programme-source", default="",\n'),
    ('        spsign, spdown, sn, _ = programme(sdelta, sperts, sgenes, official, panel)\n',
     '        spsign, spdown, sn, _ = programme(sdelta, sperts, sgenes, official, panel)\n'
     '        if a.programme_estimator == "meanz_ne":\n'
     '            spsign, spdown, sn_ne = programme_meanz_ne(sdelta, sperts, sgenes, official, panel)\n'
     '            print(f"[sameline] --programme-estimator meanz_ne: {sn_ne} of {sn} rows (non-essential)", flush=True)\n'),
    ('                "bank": a.programme_source, "bank_delta_sha256": ssha, "programme_rows": sn,\n',
     '                "bank": a.programme_source, "bank_delta_sha256": ssha, "programme_rows": sn,\n'
     '                **({} if a.programme_estimator == "pdown" else {"programme_estimator": a.programme_estimator}),\n'),
]
for old, new in E:
    if s.count(old) != 1:
        raise SystemExit(f"anchor count {s.count(old)} != 1: {old[:80]!r}")
for old, new in E:
    s = s.replace(old, new)
compile(s, str(p), "exec")
if "--check" in sys.argv:
    print("patch applies cleanly (not written)")
else:
    tmp = p.with_name(p.name + ".517tmp"); tmp.write_text(s); os.replace(tmp, p)
    print(f"patched {p}")
