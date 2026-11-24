"""373 — `blend_rows` mode for the `VCC_CTX_A_PRIOR` artifact, as exact-anchor
replacements chained AFTER 366's patch.

Applies to a COPY of `scripts/make_submission_pm_cd4.py` placed under
`work/patched/`; never writes to `scripts/`. Chain: 366's four replacements
(E1 env+loader, E2 applicators, E3 run() block, E4 config keys), then 373's
E5 loader branch, E6 applicator + dispatcher, E7 run() call-site, E8
mode-conditional config. Every anchor must occur exactly once.

    uv run python experiments/373-blend-rows-patch/patch373.py            # (re)build work/patched copy
    uv run python experiments/373-blend-rows-patch/patch373.py --check    # chain on HEAD bytes + pre/mid/post sha256
"""
from __future__ import annotations

import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WT = ROOT / "scripts/make_submission_pm_cd4.py"
DEST = HERE / "work/patched/make_submission_pm_cd4.py"
P366 = ROOT / "experiments/366-ctx-a-prior-patch/patch366.py"


def replace_once(text: str, old: str, new: str) -> str:
    n = text.count(old)
    assert n == 1, f"anchor found {n} times:\n{old[:300]}"
    return text.replace(old, new)


def load_366():
    spec = importlib.util.spec_from_file_location("patch366_chained", P366)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def patch(src: str) -> str:
    mod366 = load_366()
    src = mod366.patch(src)  # E1-E4; asserts CTX_A_PRIOR absent, anchors x1
    assert "blend_rows" not in src, "already 373-patched"

    # ------------------------------------------------------- E5a. mode tuple
    OLD5A = '''CTX_A_PRIOR_MODES = ("sign_block", "add_vector")'''
    NEW5A = '''CTX_A_PRIOR_MODES = ("sign_block", "add_vector", "blend_rows")'''
    src = replace_once(src, OLD5A, NEW5A)

    # ------------------------------------------------------- E5b. mode-aware key sets
    OLD5B = '''    allowed = {"genes", "mode", "vector", "p_down", "tau", "alpha", "scale_rule", "context"}
    if keys - allowed:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown keys {sorted(keys - allowed)}")
    lacking = {"genes", "mode", "vector"} - keys
    if lacking:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: missing keys {sorted(lacking)}")'''
    NEW5B = '''    _mode0 = str(np.asarray(z["mode"]).reshape(-1)[0]) if "mode" in keys else None
    if _mode0 is not None and _mode0 not in CTX_A_PRIOR_MODES:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown mode {_mode0!r}")
    if _mode0 == "blend_rows":
        allowed = {"targets", "genes", "matrix", "beta", "mode", "context", "norm_rule"}
    else:
        allowed = {"genes", "mode", "vector", "p_down", "tau", "alpha", "scale_rule", "context"}
    if keys - allowed:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown keys {sorted(keys - allowed)}")
    if _mode0 == "blend_rows":
        lacking = {"targets", "genes", "matrix", "beta", "mode"} - keys
    else:
        lacking = {"genes", "mode", "vector"} - keys
    if lacking:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: missing keys {sorted(lacking)}")'''
    src = replace_once(src, OLD5B, NEW5B)

    # ------------------------------------------------------- E5c. inline blend_rows parse
    # (inline, not a later-defined helper: this loader runs AT IMPORT, before any
    # function defined 2000 lines below exists. Panel-membership of `targets` is
    # unknowable here — the official/panel axes are apply-time inputs — so it is
    # checked in `apply_ctx_a_blend_rows`, which raises on any target outside the panel.)
    OLD5C = '''    genes = np.asarray(z["genes"]).astype(str).ravel()
    vector = np.asarray(z["vector"], dtype=np.float64).ravel()
    mode = str(np.asarray(z["mode"]).reshape(-1)[0])
    if mode not in CTX_A_PRIOR_MODES:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown mode {mode!r}")'''
    NEW5C = '''    mode = str(np.asarray(z["mode"]).reshape(-1)[0])
    if mode not in CTX_A_PRIOR_MODES:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown mode {mode!r}")
    if mode == "blend_rows":
        # --- 373 `blend_rows`: per-row fixed-norm blend toward an artifact matrix.
        t_raw, g_raw = np.asarray(z["targets"]), np.asarray(z["genes"])
        if t_raw.dtype.kind not in ("U", "S") or t_raw.ndim != 1:
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: targets must be a 1-D unicode array")
        if g_raw.ndim != 1:
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: genes must be a 1-D array")
        targets = t_raw.astype(str).ravel()
        genes = g_raw.astype(str).ravel()
        try:
            matrix = np.asarray(z["matrix"], dtype=np.float64)
        except (ValueError, TypeError):
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: matrix must be numeric")
        if targets.size == 0:
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: blend_rows needs >=1 target")
        if genes.size == 0:
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: blend_rows needs >=1 gene")
        if matrix.ndim != 2 or matrix.shape != (targets.size, genes.size):
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: matrix shape {matrix.shape} vs {(targets.size, genes.size)}")
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: matrix must be finite (NaN refused)")
        if len(set(targets.tolist())) != targets.size:
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: duplicate targets refused")
        if len(set(genes.tolist())) != genes.size:
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: duplicate genes refused")
        try:
            _b = np.asarray(z["beta"], dtype=np.float64)
        except (ValueError, TypeError):
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: beta must be numeric")
        if _b.size != 1:
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: beta must be scalar")
        beta = float(_b.reshape(-1)[0])
        if not np.isfinite(beta) or not 0.0 < beta <= 1.0:
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: beta must be in (0,1]; got {beta!r}")
        norm_rule = str(np.asarray(z["norm_rule"]).reshape(-1)[0]) if "norm_rule" in keys else "fixed_norm_on_covered"
        if norm_rule != "fixed_norm_on_covered":
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown norm_rule {norm_rule!r}")
        _ctx = str(np.asarray(z["context"]).reshape(-1)[0]) if "context" in keys else "A"
        if _ctx not in CTX_A_PRIOR_CONTEXTS:
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown context {_ctx!r}")
        return {"path": str(src), "sha256": file_sha256(src), "targets": targets, "genes": genes,
                "matrix": matrix, "beta": beta, "mode": mode, "context": _ctx, "norm_rule": norm_rule}
    genes = np.asarray(z["genes"]).astype(str).ravel()
    vector = np.asarray(z["vector"], dtype=np.float64).ravel()'''
    src = replace_once(src, OLD5C, NEW5C)

    # ------------------------------------------------------- E6. applicator + dispatcher
    OLD6 = '''def build_ctx_a_context_effects(effects, prior, b_ctrl, official_genes, panel_positions, contexts, **kw):
    """None when the flag is off (the writer path stays bit-identical); else one corrected context."""
    if prior is None:
        return None, None
    name = str(prior["context"])
    if name not in [str(c) for c in contexts]:
        raise ValueError(f"ctx-a-prior context {name!r} not in bundle contexts {[str(c) for c in contexts]}")
    corrected, audit = apply_ctx_a_prior(effects, prior, b_ctrl, official_genes, panel_positions, **kw)
    return {str(c): (corrected if str(c) == name else effects) for c in contexts}, audit'''
    NEW6 = '''def apply_ctx_a_blend_rows(design, prior, official_genes, panel_positions, panel_targets):
    """373 `blend_rows`: per-row fixed-norm blend toward the artifact matrix (346's rule).

    For each artifact target resolving to a design row, on the covered genes C
    (`prior["genes"]` mapped BY NAME): `r = design[row, C]`, `m = matrix[k]`;
    `m` is scaled to `||r||` on C (`||m|| == 0` raises); `blended = (1-b)*r +
    b*m_scaled`; then `blended` is rescaled to `||r||` on C and written back on
    C only. Genes outside C and targets outside the artifact are untouched."""
    genes = [str(x) for x in official_genes]
    pos_of = {x: i for i, x in enumerate(genes)}
    absent = [str(x) for x in prior["genes"] if str(x) not in pos_of]
    if absent:
        raise ValueError(f"ctx-a-prior blend_rows: {len(absent)} artifact genes not on the official axis (e.g. {absent[:3]})")
    at = np.array([pos_of[str(x)] for x in prior["genes"]], dtype=np.int64)
    if set(int(i) for i in np.asarray(panel_positions).ravel().tolist()).intersection(at.tolist()):
        raise ValueError("ctx-a-prior blend_rows: artifact carries panel-coordinate support; refused")
    if panel_targets is None:
        raise ValueError("ctx-a-prior blend_rows: panel_targets (row labels) required")
    panel = [str(t) for t in panel_targets]
    row_of = {t: i for i, t in enumerate(panel)}
    art_targets = [str(t) for t in np.asarray(prior["targets"]).ravel().tolist()]
    outside = [t for t in art_targets if t not in row_of]
    if outside:
        raise ValueError(f"ctx-a-prior blend_rows: {len(outside)} artifact targets not in the panel (e.g. {outside[:3]})")
    D = np.asarray(design, dtype=np.float64)
    if D.ndim != 2 or D.shape[1] != len(genes) or D.shape[0] != len(panel):
        raise ValueError(f"ctx-a-prior blend_rows: design {D.shape} vs ({len(panel)}, {len(genes)})")
    M = np.asarray(prior["matrix"], dtype=np.float64)
    if M.shape != (len(art_targets), len(at)):
        raise ValueError(f"ctx-a-prior blend_rows: matrix {M.shape} vs ({len(art_targets)}, {len(at)})")
    beta = float(prior["beta"])
    out = D.copy()
    cosines: list[float] = []
    written = 0
    for k, t in enumerate(art_targets):
        i = row_of[t]
        r = D[i, at]
        m = M[k]
        rn = float(np.linalg.norm(r))
        mn = float(np.linalg.norm(m))
        if mn == 0.0:
            raise ValueError(f"ctx-a-prior blend_rows: zero artifact row for target {t!r}")
        if rn == 0.0:
            continue  # row already zero on C; the fixed-norm blend is the identity there
        m_scaled = m * (rn / mn)
        blended = (1.0 - beta) * r + beta * m_scaled
        bn = float(np.linalg.norm(blended))
        if bn == 0.0:
            raise ValueError(f"ctx-a-prior blend_rows: degenerate blend for target {t!r}")
        out[i, at] = blended * (rn / bn)
        cosines.append(float(np.dot(r, m) / (rn * mn)))
        written += 1
    audit: dict[str, object] = {"mode": "blend_rows", "context": prior["context"], "rows": int(D.shape[0]),
                                "rows_touched": int(written), "artifact_targets": int(len(art_targets)),
                                "covered_genes": int(len(at)), "beta": beta,
                                "mean_cosine": float(np.mean(cosines)) if cosines else float("nan"),
                                "artifact_sha256": prior["sha256"]}
    return out, audit


def build_ctx_a_context_effects(effects, prior, b_ctrl, official_genes, panel_positions, contexts, **kw):
    """None when the flag is off (the writer path stays bit-identical); else one corrected context."""
    if prior is None:
        return None, None
    name = str(prior["context"])
    if name not in [str(c) for c in contexts]:
        raise ValueError(f"ctx-a-prior context {name!r} not in bundle contexts {[str(c) for c in contexts]}")
    panel_targets = kw.pop("panel_targets", None)
    if str(prior.get("mode")) == "blend_rows":
        corrected, audit = apply_ctx_a_blend_rows(effects, prior, official_genes, panel_positions, panel_targets)
    else:
        corrected, audit = apply_ctx_a_prior(effects, prior, b_ctrl, official_genes, panel_positions, **kw)
    return {str(c): (corrected if str(c) == name else effects) for c in contexts}, audit'''
    src = replace_once(src, OLD6, NEW6)

    # ------------------------------------------------------- E7a. run() call-site
    OLD7A = '''            tuple(map(str, spec.contexts)), head_k=RC_CALLSET_HEAD, block_k=RC_CALLSET_CPM, min_cpm=min_cpm)'''
    NEW7A = '''            tuple(map(str, spec.contexts)), head_k=RC_CALLSET_HEAD, block_k=RC_CALLSET_CPM, min_cpm=min_cpm,
            panel_targets=tuple(map(str, spec.perts)))'''
    src = replace_once(src, OLD7A, NEW7A)

    # ------------------------------------------------------- E7b. mode-aware print
    OLD7B = '''touched={ctx_a_audit['touched_entries']}'''
    NEW7B = '''touched={ctx_a_audit.get('touched_entries', ctx_a_audit.get('rows_touched'))}'''
    src = replace_once(src, OLD7B, NEW7B)

    # ------------------------------------------------------- E8. mode-conditional config
    OLD8 = '''            **({} if CTX_A_PRIOR is None else {"ctx_a_prior_mode": CTX_A_PRIOR["mode"], "ctx_a_prior_sha256": CTX_A_PRIOR["sha256"], "ctx_a_prior_alpha": CTX_A_PRIOR["alpha"], "ctx_a_prior_tau": CTX_A_PRIOR["tau"], "ctx_a_prior_context": CTX_A_PRIOR["context"]}),'''
    NEW8 = '''            **({} if CTX_A_PRIOR is None else {"ctx_a_prior_mode": CTX_A_PRIOR["mode"], "ctx_a_prior_sha256": CTX_A_PRIOR["sha256"], "ctx_a_prior_context": CTX_A_PRIOR["context"], **({"ctx_a_prior_alpha": CTX_A_PRIOR["alpha"], "ctx_a_prior_tau": CTX_A_PRIOR["tau"]} if CTX_A_PRIOR["mode"] in ("sign_block", "add_vector") else {"ctx_a_prior_beta": CTX_A_PRIOR["beta"], "ctx_a_prior_targets": int(np.asarray(CTX_A_PRIOR["targets"]).size)})}),'''
    src = replace_once(src, OLD8, NEW8)
    return src


def head_bytes() -> bytes:
    return subprocess.check_output(["git", "show", f"HEAD:scripts/{WT.name}"], cwd=ROOT)


def main() -> None:
    if len(sys.argv) >= 2 and sys.argv[1] == "--check":
        pre = head_bytes()
        mod366 = load_366()
        mid = mod366.patch(pre.decode()).encode()
        post = patch(pre.decode()).encode()
        print(f"--check: HEAD bytes {WT.name} {hashlib.sha256(pre).hexdigest()}")
        print(f"--check: 366-only  {WT.name} {hashlib.sha256(mid).hexdigest()}")
        print(f"--check: 366+373   {WT.name} {hashlib.sha256(post).hexdigest()}")
        n_pre = len(pre.decode().splitlines())
        print(f"--check: +{len(mid.decode().splitlines()) - n_pre} lines (366); "
              f"+{len(post.decode().splitlines()) - len(mid.decode().splitlines())} lines (373); "
              "all 4+7 anchors found exactly once")
        return
    DEST.parent.mkdir(parents=True, exist_ok=True)
    pre = WT.read_bytes()
    post = patch(pre.decode()).encode()
    DEST.write_bytes(post)
    print(f"patched copy {DEST}")
    print(f"pre  sha256 {hashlib.sha256(pre).hexdigest()}")
    print(f"post sha256 {hashlib.sha256(post).hexdigest()}")
    print(f"+{len(post.decode().splitlines()) - len(pre.decode().splitlines())} lines")


if __name__ == "__main__":
    main()
