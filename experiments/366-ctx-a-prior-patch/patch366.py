"""366 — the context-A prior arm (`VCC_CTX_A_PRIOR`), as exact-anchor replacements.

Applies to a COPY of `scripts/make_submission_pm_cd4.py` placed under `work/patched/`;
never writes to `scripts/`. Four replacements (E1 env+loader, E2 applicators,
E3 run() block, E4 config keys); every anchor must occur exactly once.

    uv run python experiments/366-ctx-a-prior-patch/patch366.py            # (re)build work/patched copy
    uv run python experiments/366-ctx-a-prior-patch/patch366.py --check    # anchors in HEAD bytes + pre/post sha256
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WT = ROOT / "scripts/make_submission_pm_cd4.py"
DEST = HERE / "work/patched/make_submission_pm_cd4.py"


def replace_once(text: str, old: str, new: str) -> str:
    n = text.count(old)
    assert n == 1, f"anchor found {n} times:\n{old[:300]}"
    return text.replace(old, new)


def patch(src: str) -> str:
    assert "CTX_A_PRIOR" not in src, "already patched"

    # ------------------------------------------------------- E1. env flag + loader
    OLD1 = '''    if not 0.0 < CD4_PARTIAL_G < 5.0:
        raise ValueError(f"VCC_CD4_PARTIAL_G must be in (0, 5); got {CD4_PARTIAL_G}")
else:
    CD4_PARTIAL_G = None
'''
    NEW1 = '''    if not 0.0 < CD4_PARTIAL_G < 5.0:
        raise ValueError(f"VCC_CD4_PARTIAL_G must be in (0, 5); got {CD4_PARTIAL_G}")
else:
    CD4_PARTIAL_G = None
# --- `VCC_CTX_A_PRIOR` -- one-variable context prior arm (366; DESIGN 365).
# Unset/blank is BYTE-INERT. Set to a path, the `.npz` is loaded + validated AT IMPORT;
# the named context's rows are corrected after the shared preserve, before rendering.
CTX_A_PRIOR_MODES = ("sign_block", "add_vector")
CTX_A_PRIOR_CONTEXTS = ("A", "B", "C", "D", "E", "F")
CTX_A_PRIOR_PATH = os.environ.get("VCC_CTX_A_PRIOR", "").strip() or None


def _load_ctx_a_prior(path):
    """Load + structurally validate the 366 prior artifact (import time)."""
    src = Path(path).expanduser()
    if not src.is_file():
        raise ValueError(f"VCC_CTX_A_PRIOR must be an existing .npz path; got {path!r}")
    z = np.load(src, allow_pickle=False)
    keys = set(z.files)
    allowed = {"genes", "mode", "vector", "p_down", "tau", "alpha", "scale_rule", "context"}
    if keys - allowed:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown keys {sorted(keys - allowed)}")
    lacking = {"genes", "mode", "vector"} - keys
    if lacking:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: missing keys {sorted(lacking)}")
    genes = np.asarray(z["genes"]).astype(str).ravel()
    vector = np.asarray(z["vector"], dtype=np.float64).ravel()
    mode = str(np.asarray(z["mode"]).reshape(-1)[0])
    if mode not in CTX_A_PRIOR_MODES:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown mode {mode!r}")
    if genes.shape != vector.shape or genes.size == 0 or not np.all(np.isfinite(vector)):
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: genes/vector shape {genes.shape}/{vector.shape} or non-finite")
    out = {"path": str(src), "sha256": file_sha256(src), "genes": genes, "vector": vector, "mode": mode}
    if "p_down" in keys:
        down = np.asarray(z["p_down"], dtype=np.float64).ravel()
        if down.shape != genes.shape or not np.all(np.isfinite(down)) or down.min() < 0.0 or down.max() > 1.0:
            raise ValueError(f"VCC_CTX_A_PRIOR {path}: p_down must be finite in [0,1] on the genes axis")
        out["p_down"] = down
    out["tau"] = float(np.asarray(z["tau"]).reshape(-1)[0]) if "tau" in keys else 0.0
    out["alpha"] = float(np.asarray(z["alpha"]).reshape(-1)[0]) if "alpha" in keys else 1.0
    if not np.isfinite(out["tau"]) or not 0.0 <= out["tau"] <= 0.5:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: tau must be in [0,0.5]; got {out['tau']!r}")
    if not np.isfinite(out["alpha"]) or out["alpha"] == 0.0:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: alpha must be finite nonzero; got {out['alpha']!r}")
    out["scale_rule"] = str(np.asarray(z["scale_rule"]).reshape(-1)[0]) if "scale_rule" in keys else "gate_median"
    if out["scale_rule"] != "gate_median":
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown scale_rule {out['scale_rule']!r}")
    out["context"] = str(np.asarray(z["context"]).reshape(-1)[0]) if "context" in keys else "A"
    if out["context"] not in CTX_A_PRIOR_CONTEXTS:
        raise ValueError(f"VCC_CTX_A_PRIOR {path}: unknown context {out['context']!r}")
    return out


CTX_A_PRIOR = _load_ctx_a_prior(CTX_A_PRIOR_PATH) if CTX_A_PRIOR_PATH is not None else None
if CTX_A_PRIOR is not None and ALPHA_MODE == "knee":
    raise ValueError("VCC_CTX_A_PRIOR and VCC_ALPHA=knee are two per-context mechanisms; set exactly one")
'''
    src = replace_once(src, OLD1, NEW1)

    # ------------------------------------------------------- E2. the applicators
    OLD2 = '''def build_production_effects(
'''
    NEW2 = '''def _ctx_a_head_block(row, b_ctrl, head_k, block_k):
    """361's partition: head = top-head_k by |design|; block = top-block_k by control among non-head."""
    g = int(np.asarray(row).size)
    by_abs = np.argsort(-np.abs(np.asarray(row, dtype=np.float64)), kind="stable")
    head = by_abs[: max(0, min(int(head_k), g))]
    in_head = np.zeros(g, dtype=bool)
    in_head[head] = True
    by_ctrl = np.argsort(-np.asarray(b_ctrl, dtype=np.float64), kind="stable")
    rest = by_ctrl[~in_head[by_ctrl]]
    return head, rest[: max(0, min(int(block_k), rest.size))]


def apply_ctx_a_prior(design, prior, b_ctrl, official_genes, panel_positions, *, head_k=400, block_k=2200,
                      min_cpm=10.0):
    """Correct ONE context's post-preserve design matrix -> (corrected, audit). Genes map
    BY NAME (zero misses); panel support raises. `sign_block` = sign*|design| on the block
    (p_down tau-rule else down; sign(vector), 0 = untouched; artifact-absent untouched).
    `add_vector` = row + alpha*s*vector, per-row s over b_ctrl>0 ('gate_median')."""
    genes = [str(x) for x in official_genes]
    pos_of = {x: i for i, x in enumerate(genes)}
    absent = [str(x) for x in prior["genes"] if str(x) not in pos_of]
    if absent:
        raise ValueError(f"ctx-a-prior: {len(absent)} artifact genes not on the official axis (e.g. {absent[:3]})")
    at = np.array([pos_of[str(x)] for x in prior["genes"]], dtype=np.int64)
    if set(int(i) for i in np.asarray(panel_positions).ravel().tolist()).intersection(at.tolist()):
        raise ValueError("ctx-a-prior: artifact carries panel-coordinate support; refused")
    D = np.asarray(design, dtype=np.float64)
    ctrl = np.asarray(b_ctrl, dtype=np.float64)
    if D.ndim != 2 or D.shape[1] != len(genes) or ctrl.shape != (len(genes),):
        raise ValueError(f"ctx-a-prior: design {D.shape} / b_ctrl {ctrl.shape} vs {len(genes)} genes")
    full = np.zeros(len(genes), dtype=np.float64)
    full[at] = np.asarray(prior["vector"], dtype=np.float64)
    if "p_down" in prior:
        down, have = np.zeros(len(genes)), np.zeros(len(genes), dtype=bool)
        down[at] = np.asarray(prior["p_down"], dtype=np.float64)
        have[at] = True
        sign = np.where(have & (np.abs(down - 0.5) > float(prior["tau"])), np.sign(down - 0.5), -1.0)
        sign[~have] = 0.0
    else:
        sign = np.where(full != 0.0, np.sign(full), 0.0)
    out = D.copy()
    expr = np.expm1(ctrl)
    gate = (1e6 * expr / max(float(expr.sum()), 1e-12) > float(min_cpm)) & (ctrl > 0)
    univ = ctrl > 0
    touched, scales = 0, []
    for t in range(D.shape[0]):
        head, block = _ctx_a_head_block(D[t], ctrl, head_k, block_k)
        if prior["mode"] == "sign_block":
            use = block[sign[block] != 0.0]
            out[t, use] = sign[use] * np.abs(D[t, use])
            touched += int(use.size)
        else:
            if not bool(univ.any()):
                raise ValueError("ctx-a-prior add_vector: empty control-expressed universe")
            denom = float(np.median(np.abs(full[univ])))
            if denom == 0.0:
                raise ValueError("ctx-a-prior add_vector: zero vector median on the universe (vacuous)")
            s = float(np.median(np.abs(D[t, univ]))) / denom
            scales.append(s)
            out[t] += float(prior["alpha"]) * s * full
            touched += int(np.count_nonzero(full[block]))
    audit: dict[str, object] = {"mode": prior["mode"], "context": prior["context"], "head_k": int(head_k), "block_k": int(block_k),
                                "rows": int(D.shape[0]), "touched_entries": int(touched), "support_entries": int(np.count_nonzero(full)),
                                "artifact_sha256": prior["sha256"], "gated_out_frac": float(1.0 - gate[full != 0.0].mean()) if bool((full != 0.0).any()) else 0.0}
    if scales:
        audit.update(scale_median=float(np.median(scales)), scale_min=float(np.min(scales)), scale_max=float(np.max(scales)))
    return out, audit


def build_ctx_a_context_effects(effects, prior, b_ctrl, official_genes, panel_positions, contexts, **kw):
    """None when the flag is off (the writer path stays bit-identical); else one corrected context."""
    if prior is None:
        return None, None
    name = str(prior["context"])
    if name not in [str(c) for c in contexts]:
        raise ValueError(f"ctx-a-prior context {name!r} not in bundle contexts {[str(c) for c in contexts]}")
    corrected, audit = apply_ctx_a_prior(effects, prior, b_ctrl, official_genes, panel_positions, **kw)
    return {str(c): (corrected if str(c) == name else effects) for c in contexts}, audit


def build_production_effects(
'''
    src = replace_once(src, OLD2, NEW2)

    # ------------------------------------------------------- E3. run() block
    OLD3 = '''        del frames
        gc.collect()

    del donors, pm, cd4_banks, cd4, h1_union, union2, operator, model
'''
    NEW3 = '''        del frames
        gc.collect()

    # --- `VCC_CTX_A_PRIOR` (366): correct the named context's rows AFTER the shared
    # preserve, BEFORE rendering (`context_effects` is knee's entry point; refused with it).
    ctx_a_audit: dict[str, object] | None = None
    if CTX_A_PRIOR is not None:
        if tuple(map(str, spec.contexts)) != ("A", "B", "C"):
            raise ValueError(
                "VCC_CTX_A_PRIOR is identity-bound to A/B/C (041: A = Jurkat); "
                f"got {list(map(str, spec.contexts))}"
            )
        if str(CTX_A_PRIOR["context"]) not in set(map(str, spec.contexts)):
            raise ValueError(f"VCC_CTX_A_PRIOR context {CTX_A_PRIOR['context']!r} not in this bundle")
        ctx_frame = load_control_frame(bundle, str(CTX_A_PRIOR["context"]), spec.genes, cells_per_pert=spec.cells_per_pert)
        memory_checkpoint("ctx-a-prior control frame", args.max_rss_gib)
        context_effects, ctx_a_audit = build_ctx_a_context_effects(
            effects, CTX_A_PRIOR, ctx_frame.b_ctrl, tuple(map(str, spec.genes)), np.asarray(panel_positions),
            tuple(map(str, spec.contexts)), head_k=RC_CALLSET_HEAD, block_k=RC_CALLSET_CPM, min_cpm=min_cpm)
        effect_audit["ctx_a"] = ctx_a_audit
        print(f"[ctx-a-prior] mode={CTX_A_PRIOR['mode']} context={CTX_A_PRIOR['context']} sha={CTX_A_PRIOR['sha256'][:16]} touched={ctx_a_audit['touched_entries']}", flush=True)
        del ctx_frame
        gc.collect()

    del donors, pm, cd4_banks, cd4, h1_union, union2, operator, model
'''
    src = replace_once(src, OLD3, NEW3)

    # ------------------------------------------------------- E4. the config keys
    OLD4 = '''            **({} if CD4_PARTIAL_G is None else {"cd4_partial_g": CD4_PARTIAL_G}),
'''
    NEW4 = '''            **({} if CD4_PARTIAL_G is None else {"cd4_partial_g": CD4_PARTIAL_G}),
            # Written only when the context-A prior is ON (366), the h1_union_g convention.
            **({} if CTX_A_PRIOR is None else {"ctx_a_prior_mode": CTX_A_PRIOR["mode"], "ctx_a_prior_sha256": CTX_A_PRIOR["sha256"], "ctx_a_prior_alpha": CTX_A_PRIOR["alpha"], "ctx_a_prior_tau": CTX_A_PRIOR["tau"], "ctx_a_prior_context": CTX_A_PRIOR["context"]}),
'''
    src = replace_once(src, OLD4, NEW4)
    return src


def head_bytes() -> bytes:
    return subprocess.check_output(["git", "show", f"HEAD:scripts/{WT.name}"], cwd=ROOT)


def main() -> None:
    if len(sys.argv) >= 2 and sys.argv[1] == "--check":
        pre = head_bytes()
        post = patch(pre.decode()).encode()
        print(f"--check: HEAD bytes {WT.name} {hashlib.sha256(pre).hexdigest()}")
        print(f"--check: patched   {WT.name} {hashlib.sha256(post).hexdigest()}")
        print(f"--check: +{len(post.decode().splitlines()) - len(pre.decode().splitlines())} lines; "
              "all 4 anchors found exactly once")
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
