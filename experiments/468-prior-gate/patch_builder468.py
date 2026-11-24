"""468: `VCC_DE_PRIOR_GATE=q<NN>` in the production builder -- a same-line/lab PROGRAMME as a context's block-prior default
goes only to the context's design-STRONG targets; the weakest NN % keep CLEAN4.

Why (467, H1, 3 seeds, production emitter): Arc's clean programme on every target buys fid +0.0198 / nmae -0.0148 raw but
costs reach -0.0388 (four small-gate targets whose few DE genes follow CLEAN4's cross-line consensus, not the line's own
programme). On the top 67 % by strength it keeps fid +0.0161 / nmae -0.0102 with reach +0.0015: net +0.0110 vs +0.0067.
Giving weak targets NO prior instead loses reach (467b, c4off33: reach -0.0114), so weak targets fall back to CLEAN4.
Strength = the number of gated genes whose |design row| > 0.2 (385's n_real predictor; Spearman 0.48 on H1), with the
context's own targets setting the NN-th percentile. A weak target that already has a same-line per-target prior (`pt`)
keeps it. Unset is a no-op: no code path runs, no config key is written.

usage: uv run python experiments/468-prior-gate/patch_builder468.py <builder.py> [--check]
"""
import sys
from pathlib import Path

path = Path(sys.argv[1])
src = path.read_text()
edits = [
    # 1. parse, after VCC_SAMELINE's parse (it needs SAMELINE)
    ('''    if DE_SIGN_PRIOR not in ("clean4", "clean4q80"):
        raise SystemExit("VCC_SAMELINE replaces the CLEAN4 default per context: it needs VCC_DE_SIGN_PRIOR=clean4 or clean4q80 "
                         f"(not the per-context / per-target variants); got {DE_SIGN_PRIOR!r}")
''',
     '''    if DE_SIGN_PRIOR not in ("clean4", "clean4q80"):
        raise SystemExit("VCC_SAMELINE replaces the CLEAN4 default per context: it needs VCC_DE_SIGN_PRIOR=clean4 or clean4q80 "
                         f"(not the per-context / per-target variants); got {DE_SIGN_PRIOR!r}")
# `VCC_DE_PRIOR_GATE="q<NN>"` -- 468 (2026-09-26): with a same-line/lab PROGRAMME as a context's block-prior default (VCC_SAMELINE
# without noprog), the programme goes only to the context's design-STRONG targets; the weakest NN % keep CLEAN4. Strength =
# the number of gated genes whose |design row| > 0.2 (385's n_real predictor), the threshold its NN-th percentile over the
# context's targets. H1 (467, Arc's programme, production emitter, 3 seeds): ungated fid +0.0198 / nmae -0.0148 / reach
# -0.0388 raw; q33 fid +0.0161 / nmae -0.0102 / reach +0.0015. No prior on weak targets loses reach instead (467b), hence
# CLEAN4. A weak target with a same-line per-target prior (`pt`) keeps it. Unset is a no-op.
DE_PRIOR_GATE_ENV = os.environ.get("VCC_DE_PRIOR_GATE", "").strip().lower()
DE_PRIOR_GATE: dict[str, float] | None = None
DE_PRIOR_GATE_AUDIT: dict[str, dict[str, object]] = {}
if DE_PRIOR_GATE_ENV:
    import re as _re  # noqa: PLC0415
    _pg = _re.fullmatch(r"q(\\d{1,2})", DE_PRIOR_GATE_ENV)
    if _pg is None or not 1 <= int(_pg.group(1)) <= 99:
        raise SystemExit(f"VCC_DE_PRIOR_GATE must be q<1..99>, e.g. q33: {DE_PRIOR_GATE_ENV!r}")
    if SAMELINE is None or not SAMELINE["prog"]:
        raise SystemExit("VCC_DE_PRIOR_GATE gates a programme default: it needs VCC_SAMELINE without noprog")
    DE_PRIOR_GATE = {"q": int(_pg.group(1)) / 100.0, "tau": 0.2}
'''),
    # 2. the per-context pre-pass, once the programme has replaced the default (inside the matched-context branch)
    ('''                    SAMELINE_AUDIT[str(context)] = {
                        "covered_targets": int(len(_slr["targets"])), "w": _w, "pt": bool(SAMELINE["pt"]),''',
     '''                    if DE_PRIOR_GATE is not None and _prog:
                        # 468: design strength per target (the emitter's own key: gated design row), the context's NN-th
                        # percentile as the threshold; weak targets without a same-line per-target prior keep CLEAN4
                        _pg_c4 = de_sign_prior_full(spec.genes, context)
                        _pg_s = np.empty(len(spec.perts), dtype=np.float64)
                        for _pg_t in range(len(spec.perts)):
                            _pg_g = emission_rank_axis(context_rows[_pg_t], b_ctrl, min_cpm=min_cpm)[0]
                            _pg_k = np.where(_pg_g, np.asarray(context_rows[_pg_t], dtype=np.float64), 0.0)
                            _pg_s[_pg_t] = float((np.abs(_pg_k) > DE_PRIOR_GATE["tau"]).sum())
                        _pg_thr = float(np.quantile(_pg_s, DE_PRIOR_GATE["q"]))
                        _pg_weak = [str(spec.perts[_i]) for _i in np.flatnonzero(_pg_s < _pg_thr)]
                        _pg_kept = [t for t in _pg_weak if t in sl_prior]
                        for _t in _pg_weak:
                            if _t not in sl_prior:
                                sl_prior[_t] = _pg_c4
                        DE_PRIOR_GATE_AUDIT[str(context)] = {
                            "threshold": _pg_thr, "weak_clean4": len(_pg_weak) - len(_pg_kept),
                            "weak_kept_pt": len(_pg_kept), "strong_programme": int(len(spec.perts) - len(_pg_weak)),
                            "strength_q": [float(x) for x in np.quantile(_pg_s, [0, .25, .5, .75, 1])]}
                        print(f"[prior_gate] {context}: q {DE_PRIOR_GATE['q']:g} threshold {_pg_thr:g}; programme on "
                              f"{len(spec.perts) - len(_pg_weak)} targets, CLEAN4 on {len(_pg_weak) - len(_pg_kept)} "
                              f"(+{len(_pg_kept)} weak with a pt prior kept)", flush=True)
                    SAMELINE_AUDIT[str(context)] = {
                        "covered_targets": int(len(_slr["targets"])), "w": _w, "pt": bool(SAMELINE["pt"]),'''),
    # 3. config: only when set
    ('''                                                                 "contexts": dict(SAMELINE_AUDIT)}}),''',
     '''                                                                 "contexts": dict(SAMELINE_AUDIT)}}),
                    **({} if DE_PRIOR_GATE is None else {"de_prior_gate": {
                        "q": DE_PRIOR_GATE["q"], "tau": DE_PRIOR_GATE["tau"], "weak": "clean4",
                        "strength": "count(|gated design row| > tau)", "contexts": dict(DE_PRIOR_GATE_AUDIT)}}),'''),
]
for old, new in edits:
    if src.count(old) != 1:
        raise SystemExit(f"anchor count {src.count(old)} != 1: {old[:80]!r}")
for old, new in edits:
    src = src.replace(old, new)
compile(src, str(path), "exec")
if "--check" in sys.argv:
    print("all anchors found; compiles; not writing")
else:
    path.write_text(src)
    print(f"patched {path}")
