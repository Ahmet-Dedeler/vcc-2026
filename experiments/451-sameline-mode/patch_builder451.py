"""451: add `VCC_SAMELINE` (the October same-line mode) to the production builder, by exact anchored replacements.

usage: uv run python experiments/451-sameline-mode/patch_builder451.py <builder.py> [--check]
Every anchor must occur exactly the expected number of times; otherwise nothing is written. Unset, the knob changes no
operation: every context reads `context_rows is context_effect`, `bulk_rows is bulk_effects` and an empty per-target
override, and the config gains no key.
"""
import sys
from pathlib import Path

path = Path(sys.argv[1])
src = path.read_text()

KNOB = '''DE_SIGN_SHRINK_FLOOR = 0.002
# `VCC_SAMELINE="<npz>:<sha256>:<w>[:pt][:noprog]"` -- 451 (2026-09-25): the SAME-LINE mode for a context whose line has a public bank.
# `scripts/def_sameline.py` writes the artifact from the day's identity call: per matched context, the bank's per-gene
# programme (396's P_down rule) and unit rows (non-panel genes) for the panel targets the bank covers. In a matched context
# only: (1) the DE block prior's default is the bank's programme instead of CLEAN4 (A = Nadig Jurkat reproduces
# `clean4q80ja`); (2) with w > 0, each covered target's design row and bulk row become unit(ours) + w * unit(bank row) over the
# non-panel genes, rescaled to our row's own non-panel norm (magnitude kept, direction moved; panel columns untouched);
# (3) with `pt`, a covered target's block prior is the bank row's own sign wherever the bank measures the gene; `noprog`
# keeps CLEAN4 as the context's default instead of (1). 450 (H1, Pan
# 2025 H1 -> h1-test, other lab and chemistry): programme 0.625 vs CLEAN4 0.581, per-target gate signs 0.764 vs 0.551,
# per-target centred cosine 0.260 vs 0.069. Unmatched contexts and the unset knob change nothing.
SAMELINE_ENV = os.environ.get("VCC_SAMELINE", "").strip()
SAMELINE: dict[str, object] | None = None
SAMELINE_AUDIT: dict[str, dict[str, object]] = {}
if SAMELINE_ENV:
    _sl = SAMELINE_ENV.split(":")
    _flags = _sl[3:]
    if len(_sl) < 3 or len(_sl[1]) != 64 or any(f not in ("pt", "noprog") for f in _flags) or len(set(_flags)) != len(_flags):
        raise SystemExit(f"VCC_SAMELINE must be <npz>:<sha256>:<w>[:pt][:noprog]; got {SAMELINE_ENV!r}")
    SAMELINE = {"path": Path(_sl[0]) if Path(_sl[0]).is_absolute() else ROOT / _sl[0], "sha256": _sl[1],
                "w": float(_sl[2]), "pt": "pt" in _flags, "prog": "noprog" not in _flags}
    if not 0.0 <= float(SAMELINE["w"]) <= 10.0:
        raise SystemExit(f"VCC_SAMELINE blend weight must be in [0, 10]: {SAMELINE['w']!r}")
    if DE_SIGN_PRIOR not in ("clean4", "clean4q80"):
        raise SystemExit("VCC_SAMELINE replaces the CLEAN4 default per context: it needs VCC_DE_SIGN_PRIOR=clean4 or clean4q80 "
                         f"(not the per-context / per-target variants); got {DE_SIGN_PRIOR!r}")'''

FUNCS = '''def sameline_load(genes, perts) -> dict[str, dict[str, np.ndarray]]:
    """`VCC_SAMELINE`'s artifact, verified: sha256, the official gene order and the panel order (it is panel-specific)."""
    import hashlib as _hashlib  # noqa: PLC0415
    got = _hashlib.sha256(Path(SAMELINE["path"]).read_bytes()).hexdigest()
    if got != SAMELINE["sha256"]:
        raise SystemExit(f"VCC_SAMELINE sha256 {got} != {SAMELINE['sha256']}")
    out: dict[str, dict[str, np.ndarray]] = {}
    with np.load(SAMELINE["path"], allow_pickle=False) as z:
        if list(map(str, z["genes"])) != list(map(str, genes)):
            raise SystemExit("VCC_SAMELINE artifact gene axis differs from the official gene order")
        if list(map(str, z["perts"])) != list(map(str, perts)):
            raise SystemExit("VCC_SAMELINE artifact panel differs from this bundle's (it is panel-specific)")
        for c in map(str, z["contexts"]):
            rec = {k: np.asarray(z[f"{k}_{c}"]) for k in ("programme", "targets", "rows", "signs")}
            rec["programme_ok"] = np.asarray(z[f"programme_ok_{c}"]) if f"programme_ok_{c}" in z.files else np.asarray(True)
            n = len(genes)
            if rec["programme"].shape != (n,) or rec["rows"].shape != (len(rec["targets"]), n) \\
                    or rec["signs"].shape != rec["rows"].shape or not np.isin(rec["programme"], (-1.0, 1.0)).all():
                raise SystemExit(f"VCC_SAMELINE artifact context {c}: malformed arrays")
            out[c] = rec
    return out


def sameline_blend(matrix, rec, panel_mask, w) -> np.ndarray:
    """Each covered target's row -> unit(ours) + w * unit(bank row) over the non-panel genes, at our row's own non-panel
    norm; panel columns and uncovered targets untouched. w == 0 returns `matrix` itself (no arithmetic at all)."""
    if w == 0.0:
        return matrix
    out = np.array(matrix, dtype=np.float64, copy=True)
    npm = ~np.asarray(panel_mask, dtype=bool)
    for k, ti in enumerate(np.asarray(rec["targets"], dtype=np.int64)):
        o = out[ti, npm]
        no = float(np.linalg.norm(o))
        b = np.asarray(rec["rows"][k], dtype=np.float64)[npm]
        nb = float(np.linalg.norm(b))
        if no == 0.0 or nb == 0.0:
            continue
        mix = o / no + w * b / nb
        nm = float(np.linalg.norm(mix))
        if nm > 0.0:
            out[ti, npm] = mix * (no / nm)
    return out


if NORMREL_PM is not None and NORMREL_EXP is None:'''

CTX = '''            # 451 (`VCC_SAMELINE`): a context the day's identity call matched to a public same-line bank reads the bank's
            # programme as its block-prior default and, for the targets the bank covers, blended rows. Unset or unmatched:
            # the aliases below ARE the shared matrices and the override is empty.
            context_rows = context_effect
            bulk_rows = bulk_effects
            sl_prior: dict[str, np.ndarray] = {}
            if SAMELINE is not None:
                if context_effects is not None or de_sign_full is None:
                    raise SystemExit("VCC_SAMELINE needs one shared effect matrix and VCC_DE_SIGN_PRIOR")
                if not _SAMELINE_CACHE:
                    _SAMELINE_CACHE.update(sameline_load(spec.genes, spec.perts))
                _slr = _SAMELINE_CACHE.get(str(context))
                if _slr is not None:
                    _w = float(SAMELINE["w"])
                    _prog = bool(SAMELINE["prog"]) and bool(_slr["programme_ok"])
                    if _prog:
                        de_sign_full = np.asarray(_slr["programme"], dtype=np.float64)
                    context_rows = sameline_blend(context_effect, _slr, panel_mask, _w)
                    bulk_rows = None if bulk_effects is None else sameline_blend(bulk_effects, _slr, panel_mask, _w)
                    if SAMELINE["pt"]:
                        for _k, _ti in enumerate(np.asarray(_slr["targets"], dtype=np.int64)):
                            _s = np.asarray(_slr["signs"][_k], dtype=np.float64)
                            sl_prior[str(spec.perts[int(_ti)])] = np.where(_s != 0.0, _s, de_sign_full)
                    SAMELINE_AUDIT[str(context)] = {
                        "covered_targets": int(len(_slr["targets"])), "w": _w, "pt": bool(SAMELINE["pt"]),
                        "programme_default": _prog,
                        "programme_down_frac": float((de_sign_full < 0).mean()),
                        "rows_changed": int(0 if _w == 0.0 else len(_slr["targets"])),
                    }
                    print(f"[sameline] {context}: {'programme' if _prog else 'CLEAN4'} default (down_frac {SAMELINE_AUDIT[str(context)]['programme_down_frac']:.4f}); "
                          f"{len(_slr['targets'])} covered targets blended at w {_w:g}; per-target prior {bool(SAMELINE['pt'])}",
                          flush=True)
                else:
                    print(f"[sameline] {context}: not in the artifact -- unchanged", flush=True)
            control_rows = int(control_matrix.shape[0])'''

CONF = '''"de_sign_emitter_sha256": DE_SIGN_EMITTER_SHA256,
                    **({} if SAMELINE is None else {"sameline": {"path": str(SAMELINE["path"]), "sha256": SAMELINE["sha256"],
                                                                 "w": SAMELINE["w"], "pt": SAMELINE["pt"], "prog": SAMELINE["prog"],
                                                                 "contexts": dict(SAMELINE_AUDIT)}}),'''

edits = [
    ("DE_SIGN_SHRINK_FLOOR = 0.002", KNOB, 1),
    ("if NORMREL_PM is not None and NORMREL_EXP is None:", FUNCS, 1),
    ("            control_rows = int(control_matrix.shape[0])", CTX, 1),
    ("context_effect[target_index]", "context_rows[target_index]", 6),
    ("bulk_effects[target_index] * (1.0 if BULK_CTX_SCALE is None else BULK_CTX_SCALE[str(context)])",
     "bulk_rows[target_index] * (1.0 if BULK_CTX_SCALE is None else BULK_CTX_SCALE[str(context)])", 1),
    ("DE_SIGN_PRIOR_BY_TARGET.get(str(target_value), de_sign_full)[rc_control.keep]",
     "sl_prior.get(str(target_value), DE_SIGN_PRIOR_BY_TARGET.get(str(target_value), de_sign_full))[rc_control.keep]", 1),
    ('"de_sign_emitter_sha256": DE_SIGN_EMITTER_SHA256,', CONF, 1),
    ("_ZR_GENE_POS: dict[str, int] = {}", "_ZR_GENE_POS: dict[str, int] = {}\n_SAMELINE_CACHE: dict[str, dict[str, np.ndarray]] = {}   # 451", 1),
]
for old, new, n in edits:
    got = src.count(old)
    if got != n:
        raise SystemExit(f"anchor count {got} != {n}: {old[:70]!r}")
for old, new, n in edits:
    src = src.replace(old, new)
if "--check" in sys.argv:
    print("all anchors found; not writing")
else:
    path.write_text(src)
    print(f"patched {path}")
