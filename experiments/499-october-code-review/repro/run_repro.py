"""499 reproductions. Synthetic inputs only. Does not modify tracked files.

Run from the repo root:
  uv run python experiments/499-october-code-review/repro/run_repro.py
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "experiments/499-october-code-review"
PY = sys.executable


def section(title: str) -> None:
    print(f"\n===== {title} =====")


def run(cmd: list[str], *, cwd: Path | None = None, check: bool = True, input_bytes: bytes | None = None) -> subprocess.CompletedProcess:
    r = subprocess.run(cmd, cwd=cwd or ROOT, input=input_bytes, capture_output=True)
    out = (r.stdout or b"") + (r.stderr or b"")
    text = out.decode("utf-8", "replace")
    print(f"$ {' '.join(cmd)}")
    if text.strip():
        print(text.rstrip())
    print(f"exit {r.returncode}")
    if check and r.returncode != 0:
        raise SystemExit(f"command failed: {cmd}")
    return r


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def blank_resolved(sharpen: str) -> dict[str, str]:
    keys = [
        "VCC_NORMREL_EXP", "VCC_NORMREL_MODE", "VCC_CD4_PARTIAL_G", "VCC_DE_SIGN_PRIOR",
        "VCC_PERT_SCALE", "VCC_BULK_SHARPEN_Q", "VCC_PERT_SCALE_BALLAST", "VCC_DE_MDE",
        "VCC_DE_ZERO_REPAIR", "VCC_DE_MAG",
    ]
    resolved = {k: "" for k in keys}
    resolved["VCC_BULK_SHARPEN_Q"] = sharpen
    return resolved


def call_rows(mod, cfg: dict, sharpen: str):
    saved = {k: os.environ.pop(k) for k in ("VCC_BULK_LOAD", "VCC_BULK_CTX_SCALE") if k in os.environ}
    try:
        return mod.september_rows(cfg, blank_resolved(sharpen))
    finally:
        os.environ.update(saved)


def sharpen_row(rows):
    return next(r for r in rows if r[0] in ("bulk_sharpen_q", "bulk_sharpen_pt"))


def apply_patch(patch: Path, src: Path) -> Path:
    td = Path(tempfile.mkdtemp(prefix="499-"))
    header = next(line[4:].strip() for line in patch.read_text().splitlines() if line.startswith("+++ "))
    rel = header.split("\t", 1)[0]
    if rel.startswith("b/"):
        rel = rel[2:]
    target = td / rel
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, target)
    r = subprocess.run(["patch", "-p1", "--forward", "--batch"], cwd=td, input=patch.read_bytes(), capture_output=True)
    text = (r.stdout + r.stderr).decode()
    print(f"patch {patch.name} onto {src.relative_to(ROOT) if str(src).startswith(str(ROOT)) else src}")
    print(text.rstrip())
    print(f"exit {r.returncode}")
    if r.returncode != 0:
        raise SystemExit("patch failed")
    return target


def main() -> None:
    os.chdir(ROOT)
    section("env")
    run(["bash", "--version"])
    run(["bash", "-c", "echo BASH_VERSINFO=${BASH_VERSINFO[0]}.${BASH_VERSINFO[1]}"])
    docker = shutil.which("docker") or shutil.which("podman")
    print(f"container runtime: {docker or 'none'}")

    section("F1 unpatched september_rows raises ValueError on the per-target form")
    pipe = load_module(ROOT / "scripts/def_pipeline.py", "def_pipeline_head")
    try:
        call_rows(pipe, {}, "stercile:2.5,2.0,1.75")
        print("NO RAISE")
    except ValueError as exc:
        print(f"ValueError: {exc}")
    row = sharpen_row(call_rows(pipe, {"bulk_sharpen_q": 1.75}, "1.75"))
    print(f"scalar 1.75 row: {row}")

    section("F1 patch: PT meta passes, scalar 1.75 still passes, wrong mode fails")
    patched = apply_patch(EXP / "patches/def_pipeline.py.patch", ROOT / "scripts/def_pipeline.py")
    pmod = load_module(patched, "def_pipeline_499")
    good = {"bulk_sharpen_pt": {"mode": "stercile", "q": [2.5, 2.0, 1.75]}}
    print("PT match", sharpen_row(call_rows(pmod, good, "stercile:2.5,2.0,1.75")))
    print("scalar 1.75", sharpen_row(call_rows(pmod, {"bulk_sharpen_q": 1.75}, "1.75")))
    bad = {"bulk_sharpen_pt": {"mode": "tercile", "q": [2.5, 2.0, 1.75]}}
    print("wrong mode", sharpen_row(call_rows(pmod, bad, "stercile:2.5,2.0,1.75")))
    print("unset empty", sharpen_row(call_rows(pmod, {}, "")))
    # 485's copy does not touch line 206; the same hunk must apply there too.
    apply_patch(EXP / "patches/def_pipeline.py.patch",
                ROOT / "experiments/476-def-day-merge/repro/def_pipeline.py.new")

    section("F1 shell: nonzero driver rc deletes the h5ad the driver just wrote")
    run(["bash", "-c", textwrap.dedent(r"""
        set -uo pipefail
        W=$(mktemp -d)
        h5="$W/day-lvl0.h5ad"
        driver() { printf 'counts\n' > "$h5"; printf '{"config":{"bulk_sharpen_pt":{"mode":"stercile"}}}\n' > "$h5.meta.json"; return 1; }
        rm -f "$h5" "$h5.meta.json"
        driver
        rc=$?
        echo "after driver rc=$rc h5=$([ -s "$h5" ] && echo nonempty || echo missing)"
        if [ "$rc" -ne 0 ] || [ ! -s "$h5" ]; then
          rm -f "$h5" "$h5.meta.json"
          echo "deleted h5=$([ -e "$h5" ] && echo still-there || echo gone) meta=$([ -e "$h5.meta.json" ] && echo still-there || echo gone)"
          exit 1
        fi
    """)], check=False)

    section("F1 meta_matches would accept the PT meta the driver deletes")
    meta = Path(tempfile.mkdtemp()) / "m.json"
    meta.write_text(json.dumps({"config": {
        "bulk_sharpen_pt": {"mode": "stercile", "q": [2.5, 2.0, 1.75]},
        "bulk_ctx_scale": {"D": 1.0, "E": 1.0, "F": 1.0},
    }}))
    run([PY, "experiments/497-october-switches/work/def_day_helpers.497.py",
         "meta_matches", str(meta), "D=1,E=1,F=1", "stercile:2.5,2.0,1.75"], check=False)
    run(["bash", "-c", "grep -n 'def sharpen_q_by_tercile' scripts/make_submission_pm_cd4.py || echo 'HEAD builder has no sharpen_q_by_tercile'"])
    run(["bash", "-c", "grep -n 'assert_build_carries_arms' experiments/476-def-day-merge/repro/def_pipeline.py.new | head"])

    section("F2 496 pin moves even when the knob is unset")
    td = Path(tempfile.mkdtemp(prefix="499-496-"))
    emitter = td / "rc_emit_k.py"
    builder = td / "make_submission_pm_cd4.py"
    shutil.copy(ROOT / "src/vcc2026/rc_emit_k.py", emitter)
    shutil.copy(ROOT / "scripts/make_submission_pm_cd4.py", builder)
    before = hashlib.sha256(emitter.read_bytes()).hexdigest()
    print(f"unpatched emitter sha {before}")
    run([PY, "experiments/496-uncalled-sign-hook/patch_emitter496.py", str(emitter)])
    run([PY, "experiments/496-uncalled-sign-hook/patch_builder496.py", str(builder), str(emitter)])
    after = hashlib.sha256(emitter.read_bytes()).hexdigest()
    pin = re.search(r'^DE_SIGN_EMITTER_SHA256 = "([0-9a-f]{64})"', builder.read_text(), re.M)
    print(f"patched emitter sha {after}")
    print(f"builder pin          {pin.group(1) if pin else 'MISSING'}")
    print(f"pin equals patched emitter: {pin.group(1) == after}")
    print(f"pin differs from HEAD pin: {pin.group(1) != before}")
    print(f"VCC_DE_MDE_SIGN in patched builder: {'VCC_DE_MDE_SIGN' in builder.read_text()}")
    # knob-unset: the config spread writes de_mde_sign only inside `if DE_MDE_SIGN is None` false branch
    print("de_mde_sign spread present:", 'if DE_MDE_SIGN is None else {"de_mde_sign"' in builder.read_text())

    section("F2 rule-10: a moved emitter sha is unexpected under today's ALLOW")
    inc_path = ROOT / "experiments/455-driver-q175/work/ka446b/455-ka-446b.h5ad.meta.json"
    inc = json.loads(inc_path.read_text())["config"]
    meta_cfg = dict(inc)
    meta_cfg["de_sign_emitter_sha256"] = after
    moved = sorted(k for k in set(meta_cfg) | set(inc)
                   if meta_cfg.get(k, "<absent>") != inc.get(k, "<absent>"))
    allow = set("whitener_id,pert_scale,pert_scale_sha256,bulk_load,bulk_load_sha256,bulk_ctx_scale".split(","))
    unexpected = [k for k in moved if k not in allow]
    print(f"incumbent sha {inc.get('de_sign_emitter_sha256')}")
    print(f"moved {moved}")
    print(f"unexpected {unexpected}")
    # postbuild exits 1 when unexpected is non-empty (scripts/def_postbuild.py:524)
    print(f"would SystemExit(1): {bool(unexpected)}")

    section("F2 fix: allow the key only when the source pin differs; quote-safe regex")
    day = apply_patch(EXP / "patches/def_day.497.sh.patch",
                      ROOT / "experiments/497-october-switches/work/def_day.497.sh")
    run(["bash", "-n", str(day)])
    # The broken form (backslash-quote inside bash single quotes) does not match the source line.
    run(["bash", "-c", r"""
        python3 -c 'import re,sys; t=open(sys.argv[1]).read(); m=re.search(r"^DE_SIGN_EMITTER_SHA256 = \"([0-9a-f]{64})\"", t, re.M); print("quoted", repr(m.group(1) if m else ""))' scripts/make_submission_pm_cd4.py
        python3 -c 'import re,sys; t=open(sys.argv[1]).read(); m=re.search(r"^DE_SIGN_EMITTER_SHA256 = .([0-9a-f]{64}).", t, re.M); print("dot   ", repr(m.group(1) if m else ""))' scripts/make_submission_pm_cd4.py
    """])
    allow_script = day.read_text().split("# 496's emitter pin")[1].split("uv run python scripts/def_postbuild.py")[0]
    probe = td / "probe.sh"
    probe.write_text(textwrap.dedent(f"""\
        set -uo pipefail
        cd {ROOT}
        ALLOW="whitener_id,pert_scale,pert_scale_sha256,bulk_load,bulk_load_sha256,bulk_ctx_scale"
        SQ=1.75; Q=1.75; SL=""; unset VCC_DE_PRIOR_GATE
        [ "$SQ" != "$Q" ] && ALLOW="$ALLOW,bulk_sharpen_q,bulk_sharpen_pt"
        INC={inc_path}
        {allow_script}
        echo "HEAD_ALLOW=$ALLOW"
        # point the extractor at the patched builder by shadowing the relative path
        mkdir -p "$PWD/scripts" 2>/dev/null || true
    """))
    # Run the extracted lines twice: once against HEAD builder, once with a temp copy
    # substituted by running python ourselves (the shell line hardcodes the relative path).
    run(["bash", "-c", textwrap.dedent(f"""
        set -uo pipefail
        cd {ROOT}
        extract() {{
          INC="$1"; SRC="$2"
          ALLOW="whitener_id,pert_scale,pert_scale_sha256,bulk_load,bulk_load_sha256,bulk_ctx_scale"
          SQ=1.75; Q=1.75; SL=""; unset VCC_DE_PRIOR_GATE || true
          [ "$SQ" != "$Q" ] && ALLOW="$ALLOW,bulk_sharpen_q,bulk_sharpen_pt"
          _inc_em="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["config"].get("de_sign_emitter_sha256",""))' "$INC" 2>/dev/null || true)"
          _src_em="$(python3 -c 'import re,sys; t=open(sys.argv[1]).read(); m=re.search(r"^DE_SIGN_EMITTER_SHA256 = .([0-9a-f]{{64}}).", t, re.M); print(m.group(1) if m else "")' "$SRC" 2>/dev/null || true)"
          [ -n "$_src_em" ] && [ -n "$_inc_em" ] && [ "$_src_em" != "$_inc_em" ] && ALLOW="$ALLOW,de_sign_emitter_sha256"
          echo "src=$(printf %s "$_src_em" | cut -c1-16) inc=$(printf %s "$_inc_em" | cut -c1-16) ALLOW=$ALLOW"
        }}
        extract {inc_path} scripts/make_submission_pm_cd4.py
        extract {inc_path} {builder}
    """)])

    section("F3 closing log prints q None and omits pt")
    run(["bash", "-c", textwrap.dedent(f"""
        python3 - {meta} <<'PY'
        import json, sys
        c = json.load(open(sys.argv[1]))["config"]
        print("OLD", "meta: bulk_ctx_scale", c.get("bulk_ctx_scale"), "| q", c.get("bulk_sharpen_q"), "| sameline", json.dumps(c.get("sameline"))[:300])
        print("NEW", "meta: bulk_ctx_scale", c.get("bulk_ctx_scale"), "| q", c.get("bulk_sharpen_q"), "| pt", c.get("bulk_sharpen_pt"), "| sameline", json.dumps(c.get("sameline"))[:300])
        PY
    """)])

    section("496 emitter: mde_add present, mde_add_sign absent, move_m and mde_rec unchanged")
    import numpy as np
    old_src = (ROOT / "src/vcc2026/rc_emit_k.py").read_text()
    new_src = emitter.read_text()
    start = old_src.index('            if "mde_add" in lfc_target:')
    end = old_src.index("            w_ = (nk * Lk).astype(np.float64)", start)
    old_block = textwrap.dedent(old_src[start:end])
    new_start = new_src.index('            if "mde_add" in lfc_target:')
    new_end = new_src.index("            w_ = (nk * Lk).astype(np.float64)", new_start)
    new_block = textwrap.dedent(new_src[new_start:new_end])
    n_keep = 12
    def ns(with_sign: bool):
        rng = np.random.default_rng(496)
        E0v = rng.normal(size=n_keep)
        blk = np.array([True, True, False, False, False, False, True, False, False, False, True, False])
        head = np.array([False, True, False, False, True, False, False, False, False, True, False, False])
        madd = np.abs(rng.normal(size=n_keep)) + 0.05
        lfc = {"mde_add": madd, "mde_scope": "all"}
        if with_sign:
            lfc["mde_add_sign"] = np.array([0, 0, -1, 1, 0, -1, 0, 1, 0, 0, 0, -1], dtype=np.float64)
        return dict(
            np=np, n_keep=n_keep, blk=blk, head=head, E0v=E0v, Ms=np.full(n_keep, 3.0),
            A0=np.full(n_keep, 10.0), signs=np.sign(E0v), kf=0.4,
            _t2=lambda k, e=E0v: e + 0.2, lfc_target=lfc, move_m=None, mde_rec={},
        )
    def exec_block(block, with_sign):
        g = ns(with_sign)
        exec(block, g, g)
        return g["move_m"], g["mde_rec"]
    om, od = exec_block(old_block, False)
    nm, nd = exec_block(new_block, False)
    print(f"move_m equal {np.array_equal(om, nm)}")
    print(f"mde_rec equal {od == nd} keys {sorted(od)}")
    om2, _ = exec_block(old_block, True)  # old ignores the key; same as absent
    nm2, nd2 = exec_block(new_block, True)
    print(f"with sign present, move_m differs from old {not np.array_equal(om2, nm2)}")
    print(f"with sign present, diagnostic { {k: nd2[k] for k in nd2 if k.startswith('mde_n_sign')} }")

    section("496 gene axis: P_down[keep] matches integer index and ref_lin")
    rng = np.random.default_rng(7)
    n = 18533
    P = rng.random(n)
    P[rng.choice(n, 12000, replace=False)] = np.nan
    keep = rng.random(n) > 0.7
    a = P[keep]
    b = P[np.flatnonzero(keep)]
    cpm = rng.random(n)
    ref = np.ascontiguousarray(cpm[keep])
    print(f"equal_nan {np.array_equal(a, b, equal_nan=True)} naive_equal {np.array_equal(a, b)}")
    print(f"shapes ref {ref.shape} P_down[keep] {a.shape} n_keep {int(keep.sum())}")
    prior = ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz"
    if prior.exists():
        z = np.load(prior, allow_pickle=False)
        pdn = np.asarray(z["P_down"], dtype=np.float64)
        psign = np.asarray(z["psign"])
        finite = np.isfinite(pdn)
        agree = np.where(pdn > 0.5, -1, 1)
        print(f"prior genes {z['genes'].shape} P_down {pdn.shape} finite {int(finite.sum())} "
              f"agree {int(np.sum(psign[finite] == agree[finite]))}")
    else:
        print("prior npz absent")

    section("497 knobs inert when unset; hi-pct argv; pin id")
    run(["bash", "-c", r"""
        set -uo pipefail
        Q=1.75
        unset DEF_DAY_SHARPEN_PT DEF_DAY_LODO_HI_PCT
        SQ="$Q"
        if [ -n "${DEF_DAY_SHARPEN_PT:-}" ]; then SQ="$DEF_DAY_SHARPEN_PT"; fi
        set -- ${DEF_DAY_LODO_HI_PCT:+--hi-pct "$DEF_DAY_LODO_HI_PCT"}
        echo "unset SQ=$SQ hi_argc=$#"
        DEF_DAY_LODO_HI_PCT=
        set -- ${DEF_DAY_LODO_HI_PCT:+--hi-pct "$DEF_DAY_LODO_HI_PCT"}
        echo "empty hi_argc=$#"
        DEF_DAY_LODO_HI_PCT=96.64697522769474
        set -- ${DEF_DAY_LODO_HI_PCT:+--hi-pct "$DEF_DAY_LODO_HI_PCT"}
        echo "set hi_argc=$# a1=$1 a2=$2"
        DEF_DAY_SHARPEN_PT='stercile:2.5,2.0,1.75'
        SQ="$Q"
        if [ -n "${DEF_DAY_SHARPEN_PT:-}" ]; then SQ="$DEF_DAY_SHARPEN_PT"; fi
        echo "set SQ=$SQ"
        printf 'upper [%s]\n' "$(echo 'VCC_BULK_CTX_SCALE=D=1.01,E=0.99,F=1.1' | grep -o 'VCC_BULK_CTX_SCALE=[A-Z0-9=.,]*' | cut -d= -f2-)"
        printf 'lower [%s]\n' "$(echo 'VCC_BULK_CTX_SCALE=d=1.01' | grep -o 'VCC_BULK_CTX_SCALE=[A-Z0-9=.,]*' | cut -d= -f2- || true)"
    """])
    run([PY, "-c", "import re; print(bool(re.fullmatch(r'[a-z0-9_.\\-]+', 'lodo_lo0.5_pct_def')))"])

    section("476 merge markers and rehearsed forms still present")
    resolved = (ROOT / "experiments/476-def-day-merge/work/def_day.resolved.sh").read_text()
    print("conflict markers", any(line.startswith(("<<<<<<<", "=======", ">>>>>>>")) for line in resolved.splitlines()))
    for needle in ("--cd4-thresholds \"$CD4T\"", "rm -f \"$h5\"", "past_cliff", "meta_matches", "best_safe", "UNIT="):
        print(f"resolved contains {needle!r}: {needle.replace(chr(92), '') in resolved or needle in resolved}")
    run(["bash", "-n", "experiments/476-def-day-merge/work/def_day.resolved.sh"])
    run(["bash", "-n", "experiments/497-october-switches/work/def_day.497.sh"])
    print("done")


if __name__ == "__main__":
    main()
