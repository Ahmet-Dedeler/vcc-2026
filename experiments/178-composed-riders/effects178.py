"""178 — the production EFFECTS stage, reproduced arm by arm on the real challenge panel.

WHY THIS EXISTS (AGENTS.md 11.1 rule 10, and its corollary)
-----------------------------------------------------------
Three new switches were added to `scripts/make_submission_pm_cd4.py` — `VCC_WHITEN_GAMMA`
(already there), `VCC_DONOR_WEIGHTS` and `VCC_SIGN_REPAIR`. Each is claimed to default to the
shipped behaviour bit-for-bit. That claim is only worth something if a build with the patched
builder at a15's config reproduces a15's own audit, so the KNOWN-ANSWER row here is the `a15`
arm and it is checked against `submissions/a15-cd4f010.vcc.meta.json`:

    effects.baseline_pre_sha256   7678bab9...      effects.candidate_pre_sha256 0601336f...
    effects.post_whiten? (not in a15's meta)       effects.final_sha256         2980d706...
    effects.whitening.after_all_median_norm        4.2292956044543
    effects.whitening.after_operator_median_norm   25.73792992728807
    operator sha256                                e2f65ecf...

It runs the effects stage ONLY. It emits no cells, writes nothing into `submissions/`, and
never calls `vcc`. The count build is what makes a submission; this is the part of the build
where all three riders live, and it is the part whose audit rule 10 diffs.

MEMORY
------
`main()` holds the three `Donor.delta` arrays (5.84 GB) alive next to the model's `s` copies
(another 5.84 GB) and then loads a 2.58 GB operator on top — a23's build log records
**15.9 GiB** at the three-donor model. The cap here is 10 GB with two sibling agents running,
so the stage is split in two processes:

  stage `pre`   donors -> `MultiDonorModel` -> pre-whiten effects -> PerturbMatch K562-slot
                replacement, for BOTH donor weightings, cached to `work/pre.<tag>.npy`.
                `LeanDonor` drops each donor's `delta` the moment `decompose()` has consumed
                it, so the peak is one donor's delta plus the accumulated `s`.
  stage `arms`  operator + CD4 only: whitening at the arm's gamma, the optional sign repair,
                the CD4 residual, and the audit. Cheap and re-runnable.

Nothing about the arithmetic changes: `LeanDonor` only frees memory, and the `a15` known-answer
row is what proves it (a changed `m`, `s` or ordering would move `baseline_pre_sha256`).

    uv run python effects178.py pre
    uv run python effects178.py arms
"""
from __future__ import annotations

import gc
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "src"))

# The builder is the object under test: it is imported, never copied (156's lesson, a19's
# cost).  The env is left at the shipped defaults here; the arms below drive the switches
# through the builder's own functions rather than through a second import.
os.environ.setdefault("VCC_ALPHA", "0.982")          # a15/a23
os.environ.setdefault("VCC_CD4_TARGET_F", "0.1")     # a15/a23
import importlib.util  # noqa: E402


def load_builder(**env: str):
    saved = {k: os.environ.get(k) for k in env}
    os.environ.update(env)
    try:
        name = "prod_builder_178_" + "_".join(f"{k}{v}" for k, v in sorted(env.items()))
        spec = importlib.util.spec_from_file_location(
            name or "prod_builder_178", ROOT / "scripts/make_submission_pm_cd4.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    return mod


BUNDLE = ROOT / "data/challenge"
DONOR_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
PM_PATH = ROOT / "data/external/perturbmatch"
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
OP_PATH = ROOT / "data/donors/whitener_union_pruned.npz"
OP_META = ROOT / "data/donors/whitener_union_pruned.json"

# tag -> VCC_DONOR_WEIGHTS value (None = the shipped default, i.e. the env unset)
WEIGHTINGS = {"w111": None, "w211": "2,1,1"}

# arm -> (weighting tag, gamma, sign repair)
ARMS: dict[str, tuple[str, float, bool]] = {
    "a15":            ("w111", 0.65, False),   # KNOWN ANSWER: must reproduce the shipped meta
    "g045":           ("w111", 0.45, False),
    "g030":           ("w111", 0.30, False),
    "flip":           ("w111", 0.65, True),
    "w211":           ("w211", 0.65, False),
    "flip_w211":      ("w211", 0.65, True),   # 170 killed gamma; this is the gamma-free composite
    "g045_flip":      ("w111", 0.45, True),
    "g045_w211":      ("w211", 0.45, False),
    "all3_g045":      ("w211", 0.45, True),
    "all3_g030":      ("w211", 0.30, True),
}

A15_META = json.loads((ROOT / "submissions/a15-cd4f010.vcc.meta.json").read_text())


def rss_gib() -> float:
    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (1024 ** 3) if sys.platform == "darwin" else peak / (1024 ** 2)


# --------------------------------------------------------------------------- stage `pre`
def stage_pre() -> None:
    t0 = time.time()
    M = load_builder()
    from vcc2026.donor_transfer import Donor
    from vcc2026.multi_donor import MultiDonorModel

    class LeanDonor(Donor):
        """`Donor` that frees `delta` as soon as `decompose()` has read it.

        `MultiDonorModel.__init__` touches `delta` exactly once, through `decompose()`, and
        then only uses `b_ctrl`, `genes` and `perts`. Dropping the 2.7 GB array there halves
        the peak. The `a15` known-answer sha is what proves the arithmetic is untouched.
        """

        @classmethod
        def lazy(cls, d: Path) -> "LeanDonor":
            return cls(delta=np.load(d / "delta.npy", mmap_mode="r"),
                       perts=np.asarray((d / "perts.txt").read_text().split("\n")),
                       genes=np.asarray((d / "genes.txt").read_text().split("\n")),
                       b_ctrl=np.load(d / "b_ctrl.npy"))

        def decompose(self):
            block = np.asarray(self.delta, dtype=np.float32)
            m = block.mean(axis=0)
            s = block - m
            del block
            self.delta = None
            gc.collect()
            return m, s

    spec = M.ChallengeSpec.load(BUNDLE)  # re-exported from scripts/build_submission.py
    assert len(spec.perts) == 300 and len(spec.genes) == 18533
    donors = {n: LeanDonor.lazy(p) for n, p in DONOR_PATHS.items()}
    model = MultiDonorModel(donors, spec.genes, gene_norm=False, eb=False)
    del donors
    gc.collect()
    print(f"[pre] three-donor model built, peak RSS {rss_gib():.2f} GiB "
          f"[{time.time()-t0:.0f}s]", flush=True)

    # The operator's columns are a function of the MODEL, so they are built here (where the
    # model lives) and cached with their sha for stage `arms` to re-check.
    columns = M.production_operator_columns(model, tuple(map(str, spec.perts)))
    np.save(WORK / "opcolumns.npy", np.asarray(columns, dtype=np.int64))
    col_sha = hashlib.sha256(
        np.ascontiguousarray(columns, dtype="<i8").tobytes()).hexdigest()
    print(f"[pre] operator columns {len(columns)} sha {col_sha[:16]} vs a15 "
          f"{A15_META['operator']['columns_sha256'][:16]}  "
          f"{'PASS' if col_sha == A15_META['operator']['columns_sha256'] else 'FAIL'}")
    assert col_sha == A15_META["operator"]["columns_sha256"]

    pm = M.EXP061.prepare_perturbmatch(PM_PATH, tuple(map(str, spec.genes)),
                                       chunk_rows=128, max_rss_gib=32.0)
    print(f"[pre] PerturbMatch prepared, peak RSS {rss_gib():.2f} GiB "
          f"[{time.time()-t0:.0f}s]", flush=True)

    audit: dict[str, object] = {}
    for tag, env in WEIGHTINGS.items():
        V = M if env is None else load_builder(VCC_DONOR_WEIGHTS=env)
        weights = V.donor_weights(model)
        baseline_pre = V.build_v9_pre_whiten_effects(
            model, tuple(map(str, spec.perts)), alpha=V.ALPHA, beta=V.BETA, weights=weights)
        candidate_pre, pm_audit = V.apply_perturbmatch_k562_replacement(
            baseline_pre, model, pm, tuple(map(str, spec.perts)),
            alpha=V.ALPHA, weights=weights)
        np.save(WORK / f"pre.{tag}.npy", candidate_pre)
        audit[tag] = dict(
            weights=weights, alpha=V.ALPHA,
            baseline_pre_sha256=V.array_sha256(baseline_pre),
            candidate_pre_sha256=V.array_sha256(candidate_pre),
            pm={"enabled": True, **{k: v for k, v in pm_audit.items()
                                    if not isinstance(v, list)}},
        )
        print(f"[pre] {tag:<5} weights {weights} "
              f"baseline {audit[tag]['baseline_pre_sha256'][:8]} "
              f"candidate {audit[tag]['candidate_pre_sha256'][:8]}", flush=True)
        del baseline_pre, candidate_pre
        gc.collect()

    (WORK / "pre-audit.json").write_text(json.dumps(audit, indent=1, default=str))
    # ---- rule 10, half 1: the pre-whitening shas must be a15's -----------------------
    for key in ("baseline_pre_sha256", "candidate_pre_sha256"):
        got, want = audit["w111"][key], A15_META["effects"][key]
        print(f"[KA] {key:<24} {got[:16]} vs a15 {want[:16]}  "
              f"{'PASS' if got == want else 'FAIL'}")
        assert got == want, f"{key} does not reproduce a15"
    assert audit["w111"]["candidate_pre_sha256"] != audit["w211"]["candidate_pre_sha256"]
    print(f"[pre] done, peak RSS {rss_gib():.2f} GiB [{time.time()-t0:.0f}s]", flush=True)


# -------------------------------------------------------------------------- stage `arms`
def stage_arms() -> None:
    t0 = time.time()
    M = load_builder()
    from vcc2026.whiten import Whitener

    spec = M.ChallengeSpec.load(BUNDLE)  # re-exported from scripts/build_submission.py
    genes = tuple(map(str, spec.genes))
    perts = tuple(map(str, spec.perts))
    gpos = {g: i for i, g in enumerate(genes)}
    panel_positions = np.asarray([gpos[t] for t in perts], dtype=np.int64)

    # The operator, by PATH and SHA (rule 11) — never by description.
    op_sha = M.file_sha256(OP_PATH)
    op_meta = json.loads(OP_META.read_text())
    print(f"[arms] operator {OP_PATH.relative_to(ROOT)} sha {op_sha[:16]} "
          f"n_cols {op_meta['n_cols']} prune {op_meta.get('basis_prune')}", flush=True)
    assert op_sha == A15_META["operator"]["sha256"], "operator is not the one a15 shipped"

    # `columns` as `production_operator_columns` built them in stage `pre`, re-checked here
    # against the SHIPPED meta's `columns_sha256` — a19's near-miss was a silently different
    # operator, and the column set is the other half of the operator's identity.
    columns = np.load(WORK / "opcolumns.npy")
    col_sha = hashlib.sha256(
        np.ascontiguousarray(columns, dtype="<i8").tobytes()).hexdigest()
    assert col_sha == A15_META["operator"]["columns_sha256"], "operator columns are not a15's"
    assert len(columns) == int(op_meta["n_cols"]) == A15_META["operator"]["columns"]

    whitener = Whitener.load(OP_PATH)
    print(f"[arms] whitener loaded {whitener.evecs.shape} peak RSS {rss_gib():.2f} GiB "
          f"[{time.time()-t0:.0f}s]", flush=True)

    cd4_banks, cd4 = M.prepare_frozen_cd4(CD4_ROOT, perts, genes, perts,
                                          chunk_rows=128, max_rss_gib=32.0)
    del cd4_banks
    gc.collect()
    print(f"[arms] CD4 prepared, peak RSS {rss_gib():.2f} GiB [{time.time()-t0:.0f}s]",
          flush=True)

    pre_cache = {tag: np.load(WORK / f"pre.{tag}.npy") for tag in WEIGHTINGS}
    out: dict[str, object] = {}
    for arm, (tag, gamma, repair) in ARMS.items():
        env = {}
        if WEIGHTINGS[tag] is not None:
            env["VCC_DONOR_WEIGHTS"] = WEIGHTINGS[tag]
        if gamma != 0.65:
            env["VCC_WHITEN_GAMMA"] = f"{gamma}"
        if repair:
            env["VCC_SIGN_REPAIR"] = "cons"
        V = M if not env else load_builder(**env)
        assert V.WHITEN_GAMMA == gamma and (V.SIGN_REPAIR == "cons") == repair

        candidate_pre = pre_cache[tag]
        whitened, whitening_audit = V.apply_v9_whitening(candidate_pre, whitener, columns)
        if V.SIGN_REPAIR is None:
            sr_audit = {"enabled": False}
        else:
            whitened, sr_audit = V.sign_repair_to_consensus(whitened, candidate_pre, columns)
            sr_audit = {"enabled": True, **sr_audit}
        final, cd4_audit = V.apply_frozen_cd4_residual(whitened, cd4, columns,
                                                       panel_positions)
        np.save(WORK / f"eff.{arm}.npy", final.astype(np.float64))
        out[arm] = dict(
            weights=WEIGHTINGS[tag] or "1,1,1",
            gamma=gamma, sign_repair=bool(repair),
            candidate_pre_sha256=V.array_sha256(candidate_pre),
            final_sha256=V.array_sha256(final),
            whitening=whitening_audit, sign_repair_audit=sr_audit,
            cd4_fraction=float(cd4_audit["fraction"]),
            cd4_effective_f=float(
                cd4_audit["residual_frame"]["effective_f_at_cd4_fraction"]),
            cd4_whitened_base_median_row_norm=float(
                cd4_audit["residual_frame"]["whitened_base_median_row_norm"]),
        )
        print(f"[arms] {arm:<11} gamma {gamma:<5} w {WEIGHTINGS[tag] or '1,1,1':<6} "
              f"repair {int(repair)}  after_all {whitening_audit['after_all_median_norm']:.10f}  "
              f"final {out[arm]['final_sha256'][:8]}  "
              f"flips {sr_audit.get('flipped_entries', 0)}", flush=True)
        del whitened, final
        gc.collect()

    (WORK / "arms-audit.json").write_text(json.dumps(out, indent=1, default=str))

    # ---- rule 10, half 2: the a15 arm's audit against the SHIPPED metas ---------------
    # Both are diffed. a15 was built on AWS (`/work/exp089-.../input/...`) and a23 on this
    # machine, and the two shipped metas ALREADY disagree with each other in the last ULP of
    # the whitening audit and in `post_whiten_sha256` while agreeing exactly on
    # `candidate_pre_sha256` — i.e. the difference is the operator matmul's summation order
    # on different hardware, upstream of anything this patch touches. a23 is the incumbent
    # the riders would be added to, so a23 is the binding comparison and a15 is reported.
    A23_META = json.loads(
        (ROOT / "submissions/a23-massfix.h5ad.meta.json").read_text())
    a = out["a15"]
    for label, meta in (("a23-massfix (incumbent, built HERE)", A23_META),
                        ("a15-cd4f010 (built on AWS)", A15_META)):
        e = meta["effects"]
        print(f"\n=== rule 10: the a15-config arm vs {label} ===")
        checks = [
            ("effects.candidate_pre_sha256", a["candidate_pre_sha256"],
             e["candidate_pre_sha256"]),
            ("effects.final_sha256", a["final_sha256"], e["final_sha256"]),
            ("whitening.after_all_median_norm", a["whitening"]["after_all_median_norm"],
             e["whitening"]["after_all_median_norm"]),
            ("whitening.after_operator_median_norm",
             a["whitening"]["after_operator_median_norm"],
             e["whitening"]["after_operator_median_norm"]),
            ("whitening.before_operator_median_norm",
             a["whitening"]["before_operator_median_norm"],
             e["whitening"]["before_operator_median_norm"]),
            ("whitening.operator_preserve_scale",
             a["whitening"]["operator_preserve_scale"],
             e["whitening"]["operator_preserve_scale"]),
            ("cd4.fraction", a["cd4_fraction"], e["cd4"]["fraction"]),
            ("cd4.whitened_base_median_row_norm",
             a["cd4_whitened_base_median_row_norm"],
             e["cd4"]["residual_frame"]["whitened_base_median_row_norm"]),
            ("operator.sha256", op_sha, meta["operator"]["sha256"]),
        ]
        ndiff = 0
        for name, got, want in checks:
            same = got == want
            ndiff += 0 if same else 1
            shown_got = got[:16] if isinstance(got, str) else f"{got!r}"
            shown_want = want[:16] if isinstance(want, str) else f"{want!r}"
            print(f"  {name:<40} {shown_got:<24} vs {shown_want:<24} "
                  f"{'EXACT' if same else 'differs'}")
        print(f"  EXACT-match differences: {ndiff}  "
              f"{'(rule 10 clean)' if ndiff == 0 else ''}")

    # every arm's difference count against the a15 arm, in the audit the build writes
    print("\n=== per-arm difference count against the a15 arm's audit ===")
    base = out["a15"]
    for arm, row in out.items():
        diffs = [k for k in ("weights", "gamma", "sign_repair")
                 if row[k] != base[k]]
        moved = row["final_sha256"] != base["final_sha256"]
        print(f"  {arm:<11} config differences {len(diffs)} {diffs}  "
              f"effect matrix moved: {moved}")
    print(f"[arms] done, peak RSS {rss_gib():.2f} GiB [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    {"pre": stage_pre, "arms": stage_arms}[sys.argv[1]]()
