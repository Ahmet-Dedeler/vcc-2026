"""The Oct-22 path as ONE command, stage by stage, timed, fail-closed. Experiment 199.

    preflight -> identity -> cd4-extract -> census -> whitener -> repin-frozen -> pin-artifacts
             -> repin-panel-support -> build -> postbuild -> prep

Every stage is a subprocess wrapped in `/usr/bin/time -l`, so its wall clock and peak RSS are
captured into `<work>/timings.json` as it goes; the run stops at the first non-zero exit.
Nothing here submits: `vcc submit` is never invoked and the driver refuses to run it.

The two heavy stages (`repin-panel-support`, `build`) need ~17-20 GiB (AGENTS.md: the three-donor
model load is the floor; 180 §11). The driver refuses to start them while another builder is
running or while free+inactive memory is under `--min-free-gib`.

Rehearsal (synthetic bundle, throwaway frozen table, throwaway whitener id, builder restored):

    uv run python scripts/def_pipeline.py \
        --bundle experiments/199-def-runbook/work/def-bundle \
        --work experiments/199-def-runbook/work/run1 \
        --frozen-inputs experiments/199-def-runbook/work/run1/frozen_inputs.json \
        --whitener-id union_pruned_def199 --name def199-rehearsal \
        --cd4-thresholds relaxed --restore-builder-after

The round (real bundle, the production table, the builder edit COMMITTED afterwards):

    uv run python scripts/def_pipeline.py --bundle data/challenge \
        --work experiments/<round-dir>/work --name <opaque-name> --whitener-id union_pruned_def

Stage subsets: `--stages preflight,identity,census` or `--from build`. The exact command of
every stage is printed and logged (`<work>/logs/<stage>.log`), so the runbook and the run can
never disagree.

The configuration it builds is `SHIPPED_ENV` + `SHIPPED_FLAGS` below: the currently shipped
`h1union` (board 0.2118310), carrying BOTH direction arms. Those environment variables are
DEFAULTS -- export one and yours wins (experiment 240; before it they were overrides, and the
driver silently rebuilt the two-arms-ago incumbent). After the build the driver re-reads the
arms off the finished artifact and refuses the package if they are not there.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGES = ["preflight", "identity", "cd4-extract", "census", "whitener", "key", "repin-frozen",
          "pin-artifacts", "repin-panel-support", "build", "postbuild", "prep"]
# `vcc prep` is heavy too: measured 21.7 GiB max RSS on a 13.2 GB artifact (the CLI reads the
# whole X into memory), so it gets the same "alone, and enough free memory" guard as the builds.
HEAVY = {"repin-panel-support", "build", "prep"}
CD4_SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
CD4_SOURCE_SHA256 = "c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62"
CD4_SOURCE_BYTES = 16_786_240_107
ABC_OPERATOR_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"

# --- THE CONFIGURATION THE DAY MUST CARRY (240; the defect is experiment 239 section 2.1) -----
#
# These are DEFAULTS, not overrides. `resolve_env` puts `os.environ` ON TOP of them, so an
# operator who exports one of these WINS. Before 240 this dict went last in the subprocess env
# merge and silently reset the operator's value; run as the runbook documented it, that dropped
# BOTH board-validated direction arms from the only prize-eligible entry -- `-0.0199671` overall,
# larger than our whole gap to first place -- and `def_postbuild.py` passed, because it diffed
# against `rc_callset`, the incumbent from before both arms.
#
# Pinned to the SHIPPED entry `h1union` (board 0.2118310, rank 3/656, entry
# `iUziqE1gewfcPKzYUDaB`, 2026-09-04) -- NOT to 191's `rc_hybrid` or 202's `rc_callset`. The two
# arms this carries, and what the board paid for them:
#
#   VCC_CD4_TARGET_F  0.10 -> 0.34   `cd4f034`  +0.0130557 overall  (entry QSe3JlNUYZ2BzvlcsL14)
#   VCC_H1_UNION_G         = 0.28    `h1union`  +0.0069114 overall  (entry iUziqE1gewfcPKzYUDaB)
#
# `VCC_UNION_SOURCES` is deliberately NOT in this dict: unset selects the H1 pair, which is the
# source selection 229 landed (0204168) and 235 hardened (b540148), and is 229's known-answer
# row 2. Setting it here would make the default path depend on a registry ordering; leaving it
# unset is what the shipped arm did. `SHIPPED_UNION_SOURCES` is what "unset" resolves to, and
# the post-build assertion checks the built config against it.
SHIPPED_ENV = {"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.34", "VCC_H1_UNION_G": "0.28",
               # 432 (2026-09-23): the September knobs that need no panel artifact -- 433h's config
               # (417e stack LOAD + coarse MDE + zero repair + re-level) minus the two DAY-OF values
               # (VCC_BULK_LOAD, VCC_BULK_CTX_SCALE), which are exported on the day, never defaulted.
               "VCC_NORMREL_EXP": "-1.5", "VCC_NORMREL_MODE": "cap", "VCC_CD4_PARTIAL_G": "0.5",
               "VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_PERT_SCALE_BALLAST": "1", "VCC_BULK_SHARPEN_Q": "1.75",
               "VCC_DE_MDE": "0.1,0.1,0.3,coarse", "VCC_DE_ZERO_REPAIR": "-5"}
SHIPPED_UNION_SOURCES = ["gse295214_h1", "h1_2025"]
# --- 247: CD4's dose on a NEW panel, and ONLY there --------------------------------------------
#
# `solve_cd4_fraction` sizes CD4 by the median row norm over ALL 300 rows, which over-doses the
# rows CD4 actually carries by `median_carrying / median_all` whenever coverage is below ~92%
# (242 #2): 1.196x on DEF-RUNBOOK section 5's pessimistic draw, measured on the real banks (247).
# On a non-A/B/C bundle the driver therefore sizes CD4 on its CARRYING rows, at the dose that
# reproduces the A/B/C DELIVERED carrying-row dose exactly:
#     g_cd4 = 0.34 * (median_carrying / median_all on A/B/C) = 0.34 * 1.0128397878381297
# At this `g` the A/B/C `effects.cd4.fraction` reproduces 0.32207201067578145 to 5.6e-17
# (numeric identity, one ULP). A/B/C itself stays on the legacy path -- these two variables are
# NOT set there -- so 240's row 1 (the driver's defaults rebuild `8df4dadd.../f13ce40b...`
# byte-for-byte) is untouched. Both are defaults an exported value overrides, like SHIPPED_ENV.
DEF_ENV = {"VCC_CD4_SIZING": "carrying", "VCC_CD4_TARGET_G": "0.3443655278649641",
           # 256 (BR249, 2026-09-05): the CD4 deflation SHIPS on D/E/F as the PINNED A/B/C basis
           # (251), never the day-of re-fit (250: the re-fit reads <= 0 on both low-coverage draws).
           "VCC_CD4_DEFLATE_BASIS": "V5_abc",
           # 432: 410's key RE-FIT for this panel on the day's operator (key432.py), registered under
           # this id by the proposed `def_pin_panel_artifacts.py --pert-scale` (DEF-RUNBOOK 410 addition).
           "VCC_PERT_SCALE": "wshare_lo0.3_hi3_def"}
# --- 256: the CD4 deflation on A/B/C ONLY -----------------------------------------------------
#
# `cd4svd5` (249; board 0.2160944, +0.0042634, entry Xv0ayvuoHpsKYWv2iu42, 2026-09-05) was built
# with the day-of RE-FIT `VCC_CD4_DEFLATE_K=5`, and that is what reproduces it byte-for-byte
# (240's row 1, re-run by 256). On A/B/C the pinned basis equals the re-fit only to 2.4e-15
# (251 KA 2), so A/B/C keeps the re-fit and every other panel carries the pin. The builder
# refuses both variables at once (rule 10), which is why they live in two tables.
ABC_ENV = {"VCC_CD4_DEFLATE_K": "5", "VCC_PERT_SCALE": "wshare_lo0.3_hi3"}
# What each arm is worth on the board, for the message the guard raises with.
ARM_BOARD_DELTA = {"VCC_CD4_TARGET_F": "+0.0130557 (`cd4f034`)",
                   "VCC_H1_UNION_G": "+0.0069114 (`h1union`)",
                   "VCC_ALPHA": "the a11->a15 amplitude series; 211 closed it at 0.982",
                   "VCC_CD4_SIZING": "carrying-row CD4 sizing off A/B/C (247): without it the "
                                     "cd4f034 dose inflates by median_carrying/median_all",
                   "VCC_CD4_TARGET_G": "A/B/C's delivered carrying-row dose 0.3443655278649641 (247)",
                   "VCC_CD4_DEFLATE_K": "+0.0042634 (`cd4svd5`, the A/B/C re-fit, entry Xv0ayvuoHpsKYWv2iu42)",
                   "VCC_NORMREL_EXP": "per-target donor weighting, 097 route (a) a6 / 308 `normrel1.5`; "
                                       "local PLOLO +0.00259 [+0.00000, +0.00527] 3/3 folds, "
                                       "priced +0.00104 pessimistic -- UNTESTED on the board",
                   "VCC_NORMREL_PM": "the PARITY-CORRECT rendering of VCC_NORMREL_EXP (316, to 313 "
                                     "section 8): PerturbMatch folded in at the arm's own per-target K562 "
                                     "weight and INSIDE the norm match, so the row is exactly shape-only. "
                                     "The shipped weight-first form scored +0.0016245 (Q519I2pLWRCW7PAEx5ul) "
                                     "with score_nmae -0.007705 -- UNTESTED on the board",
                   "VCC_NORMREL_MODE": "`cap` = the SUPPRESSION-ONLY half of VCC_NORMREL_EXP (319): the "
                                       "ratio is clamped at 1 before the exponent, so a donor row at or "
                                       "below its own panel median keeps weight 1 and only OVERSIZED rows "
                                       "are down-weighted. `plain` (the default) is the shipped two-sided "
                                       "rule, +0.0016245 at exp -1.5 (Q519I2pLWRCW7PAEx5ul) -- `cap` is "
                                       "UNTESTED on the board",
                   "VCC_CD4_DEFLATE_BASIS": "+0.0042634 (`cd4svd5`) carried to a new panel as the PINNED "
                                            "A/B/C basis (251); never the re-fit off A/B/C (250)",
                   "VCC_CD4_PARTIAL_G": "partial CD4 rows, board steps 319p/319p2 (+0.0011593 / +0.0002969, FINDINGS s0)",
                   "VCC_DE_SIGN_PRIOR": "CLEAN4 sign prior q80 (392 +0.0025062, 395 +0.0013142)",
                   "VCC_PERT_SCALE": "410's per-target bulk key (410b +0.0170570 with ballast); the artifact is "
                                     "panel- and operator-specific: re-fit on the day (432)",
                   "VCC_PERT_SCALE_BALLAST": "410b ballast cells: +0.0170570 over the welded 410",
                   "VCC_BULK_SHARPEN_Q": "bulk sharpening, levelled: q 1.25 +0.0069612 (417a), q 1.5 +0.0046475 (433k), q 1.75 +0.0022834 (446b, UDgjOKL9gI1udAexq6ue)",
                   "VCC_DE_MDE": "426/433 coarse MDE (H1 nmae -0.0095; board read pending)",
                   "VCC_DE_ZERO_REPAIR": "428 zero-gene repair (cannot hurt; A/B/C 24 genes at lfc ~ -35)"}
PRODUCTION_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")


def day_env(is_abc: bool) -> dict[str, str]:
    """The day's DEFAULTS: the shipped arms, plus the CD4 deflation as the re-fit on A/B/C and as
    the pinned basis + CD4's panel-invariant sizing on any other panel (256)."""
    return {**SHIPPED_ENV, **ABC_ENV} if is_abc else {**SHIPPED_ENV, **DEF_ENV}
# The build flags, byte-for-byte from experiments/222-h1-union-build/build222.sh, which is
# 210's, which is 202's and 191's -- the flags never moved, only the env did.
SHIPPED_FLAGS = ["--min-cpm", "10", "--tail-spread", "1.0", "--keep-top", "8000",
                 "--keep-spread", "0.25", "--overdispersion", "4.0"]


def resolve_env(defaults: dict[str, str] = SHIPPED_ENV) -> tuple[dict[str, str], list[str]]:
    """The day's configuration, as DEFAULTS an operator can override.

    A variable PRESENT in `os.environ` is the operator's and wins -- including when it is
    empty, which the builder reads as "this arm is off" (`VCC_H1_UNION_G=""` is how 222's own
    known-answer row turns the arm off). Returns the resolved env and the keys the operator
    supplied, so every caller can say where each value came from.
    """
    resolved, from_operator = {}, []
    for key, default in defaults.items():
        if key in os.environ:
            resolved[key] = os.environ[key]
            from_operator.append(key)
        else:
            resolved[key] = default
    return resolved, from_operator


def _as_float(value: str) -> float | None:
    """The builder's own reading of these variables: blank means unset/off."""
    value = (value or "").strip()
    return float(value) if value else None


def september_rows(cfg: dict, resolved: dict[str, str]) -> list[tuple[str, object, object, bool]]:
    """432: the September knobs, re-read off the written config by value. A knob the operator exported
    EMPTY (the builder's "off") must be absent from the config. VCC_BULK_LOAD / VCC_BULK_CTX_SCALE are day-of
    values: checked only when exported, never defaulted."""
    def val(key: str) -> str:
        return (resolved.get(key, os.environ.get(key, "")) or "").strip()
    rows: list[tuple[str, object, object, bool]] = []
    def scalar(env: str, key: str, conv) -> None:
        v = val(env)
        want = conv(v) if v else "<absent>"
        got = cfg.get(key, "<absent>")
        rows.append((key, want, got, got == want))
    scalar("VCC_NORMREL_EXP", "normrel_exp", float)
    if val("VCC_NORMREL_EXP"):   # the builder always writes normrel_mode alongside the exponent
        scalar("VCC_NORMREL_MODE", "normrel_mode", lambda s: s.lower())
    scalar("VCC_CD4_PARTIAL_G", "cd4_partial_g", float)
    scalar("VCC_DE_SIGN_PRIOR", "de_sign_prior", lambda s: s.lower())
    scalar("VCC_PERT_SCALE", "pert_scale", lambda s: s.lower())
    # 499: `tercile:` / `stercile:` is 479's per-target form (DEF_DAY_SHARPEN_PT). It is not a
    # float, and the scalar config key is absent. float(sq) raised ValueError after a good
    # h5ad existed; def_day.sh then deleted that h5ad and stopped.
    sq = val("VCC_BULK_SHARPEN_Q")
    pt_mode = sq.split(":", 1)[0].strip().lower()
    if pt_mode in ("tercile", "stercile", "tercilec"):   # 598: tercilec, the carried-only key tercile
        _, _, qs = sq.partition(":")
        try:
            want_q = [float(x) for x in qs.split(",")]
        except ValueError:
            want_q = []
        pt = cfg.get("bulk_sharpen_pt")
        qv = pt.get("q") if isinstance(pt, dict) else None
        ok = ( "bulk_sharpen_q" not in cfg and isinstance(pt, dict) and pt.get("mode") == pt_mode
              and len(want_q) == 3 and isinstance(qv, (list, tuple)) and len(qv) == 3
              and all(abs(float(a) - b) < 5e-5 for a, b in zip(qv, want_q)))
        rows.append(("bulk_sharpen_pt", {"mode": pt_mode, "q": want_q},
                     {"mode": pt.get("mode"), "q": list(qv)} if isinstance(qv, (list, tuple)) else pt, ok))
    else:
        scalar("VCC_BULK_SHARPEN_Q", "bulk_sharpen_q", float)
    b = val("VCC_PERT_SCALE_BALLAST")
    rows.append(("pert_scale_ballast", "<present>" if b else "<absent>", cfg.get("pert_scale_ballast", "<absent>"),
                 ("pert_scale_ballast" in cfg) == bool(b)))
    m = val("VCC_DE_MDE")
    got = cfg.get("de_mde", "<absent>")
    if m:
        p = m.split(",")
        want = {"c": float(p[0]), "cap": float(p[1]), "min_count": float(p[2]), "coarse": len(p) == 4}
        # the builder writes "coarse": True only for the coarse form (426a's plain meta has no key)
        ok = isinstance(got, dict) and all(got.get(k, False if k == "coarse" else None) == v for k, v in want.items())
    else:
        want, ok = "<absent>", got == "<absent>"
    rows.append(("de_mde", want, got, ok))
    z = val("VCC_DE_ZERO_REPAIR")
    got = cfg.get("de_zero_repair", "<absent>")
    rows.append(("de_zero_repair.floor_log2", float(z) if z else "<absent>",
                 got.get("floor_log2") if isinstance(got, dict) else got,
                 (isinstance(got, dict) and got.get("floor_log2") == float(z)) if z else got == "<absent>"))
    g = val("VCC_DE_MAG")
    got = cfg.get("de_mag", "<absent>")
    rows.append(("de_mag.variant", g or "<absent>", got.get("variant") if isinstance(got, dict) else got,
                 (isinstance(got, dict) and got.get("variant") == g) if g else got == "<absent>"))
    if os.environ.get("VCC_BULK_LOAD", "").strip():
        want = os.environ["VCC_BULK_LOAD"].strip().rsplit(":", 1)[1].lower()
        rows.append(("bulk_load_sha256", want, cfg.get("bulk_load_sha256", "<absent>"),
                     cfg.get("bulk_load_sha256") == want))
    if os.environ.get("VCC_BULK_CTX_SCALE", "").strip():
        want = {k.strip(): float(v) for k, v in (p.split("=") for p in os.environ["VCC_BULK_CTX_SCALE"].split(","))}
        rows.append(("bulk_ctx_scale", want, cfg.get("bulk_ctx_scale", "<absent>"), cfg.get("bulk_ctx_scale") == want))
    return rows


def assert_build_carries_arms(meta_path: Path, resolved: dict[str, str]) -> dict[str, str]:
    """RAISE unless the config the builder WROTE carries the arms the driver ASKED for.

    235 applied this discipline to the operator axis: do not trust that the thing you passed is
    the thing that got built -- re-read it off the finished artifact. This is the same rule on
    the env axis, and it is the check that would have caught 239 section 2.1 on the day: the
    pre-240 driver asked for `f = 0.34`, built `f = 0.10`, and nothing looked.

    Compared against the RESOLVED intent, not against `SHIPPED_ENV`, so an operator override is
    honoured and still verified end to end. `def_postbuild.py` is the layer that additionally
    requires the arms the runbook says we ship.
    """
    cfg = json.loads(Path(meta_path).read_text())["config"]
    want_g = _as_float(resolved.get("VCC_H1_UNION_G", ""))
    want_sources = [s for s in (os.environ.get("VCC_UNION_SOURCES", "").strip().split(","))
                    if s] or SHIPPED_UNION_SOURCES
    rows: list[tuple[str, object, object, bool]] = []

    alpha_env = (resolved.get("VCC_ALPHA", "") or "").strip()
    if alpha_env.lower() == "knee":
        rows.append(("alpha_mode", "knee", cfg.get("alpha_mode", "<absent>"),
                     cfg.get("alpha_mode") == "knee"))
    else:
        want_alpha = _as_float(alpha_env)
        rows.append(("alpha", want_alpha if want_alpha is not None else "<builder default>",
                     cfg.get("alpha", "<absent>"),
                     want_alpha is None or cfg.get("alpha") == want_alpha))

    want_f = _as_float(resolved.get("VCC_CD4_TARGET_F", ""))
    rows.append(("cd4_target_f", want_f, cfg.get("cd4_target_f", "<absent>"),
                 cfg.get("cd4_target_f") == want_f))

    rows.append(("h1_union_g", want_g, cfg.get("h1_union_g", "<absent>"),
                 cfg.get("h1_union_g", None) == want_g))
    if want_g is not None:
        rows.append(("h1_union_sources", want_sources, cfg.get("h1_union_sources", "<absent>"),
                     list(cfg.get("h1_union_sources") or []) == want_sources))

    # 247: CD4's sizing rule. On A/B/C the driver does not resolve these two, so an operator's
    # export is read straight from the environment; on D/E/F `resolved` carries DEF_ENV.
    want_sizing = (resolved.get("VCC_CD4_SIZING", os.environ.get("VCC_CD4_SIZING", ""))
                   or "").strip().lower() or None
    rows.append(("cd4_sizing", want_sizing if want_sizing else "<absent: the all-rows median>",
                 cfg.get("cd4_sizing", "<absent>"),
                 (cfg.get("cd4_sizing") == want_sizing) if want_sizing else "cd4_sizing" not in cfg))
    if want_sizing:
        want_gc = _as_float(resolved.get("VCC_CD4_TARGET_G", os.environ.get("VCC_CD4_TARGET_G", "")))
        rows.append(("cd4_target_g", want_gc, cfg.get("cd4_target_g", "<absent>"),
                     cfg.get("cd4_target_g") == want_gc))

    # 256: exactly one CD4 deflation variable, by panel (the builder refuses both at once). On
    # A/B/C `resolved` carries ABC_ENV's re-fit k; elsewhere DEF_ENV's pinned basis; an operator's
    # export of the other one is read straight from the environment, as for the sizing rows.
    want_k = (resolved.get("VCC_CD4_DEFLATE_K", os.environ.get("VCC_CD4_DEFLATE_K", "")) or "").strip()
    want_k_int = int(want_k) if want_k.isdigit() else None
    rows.append(("cd4_deflate_k", want_k_int if want_k else "<absent>", cfg.get("cd4_deflate_k", "<absent>"),
                 (cfg.get("cd4_deflate_k") == want_k_int) if want_k else "cd4_deflate_k" not in cfg))
    want_basis = (resolved.get("VCC_CD4_DEFLATE_BASIS", os.environ.get("VCC_CD4_DEFLATE_BASIS", "")) or "").strip()
    rows.append(("cd4_deflate_basis", want_basis or "<absent>", cfg.get("cd4_deflate_basis", "<absent>"),
                 (cfg.get("cd4_deflate_basis") == want_basis) if want_basis else "cd4_deflate_basis" not in cfg))

    rows.extend(september_rows(cfg, resolved))
    bad = [r for r in rows if not r[3]]
    lines = [f"  {'OK ' if ok else 'NO '} config.{key}: asked {want!r}, the artifact carries {got!r}"
             for key, want, got, ok in rows]
    print("[arms] the built config against the configuration asked for:", flush=True)
    print("\n".join(lines), flush=True)
    if bad:
        env_arms = {k: v for k, v in resolved.items() if k in ARM_BOARD_DELTA}
        raise SystemExit(
            "REFUSING THE BUILD: the artifact does not carry the configuration it was built for.\n"
            + "\n".join(lines)
            + f"\n  meta        {meta_path}\n"
            f"  asked for   {env_arms}\n"
            "  A dropped arm is silent everywhere else -- the artifact builds, `vcc prep` passes,\n"
            "  and the board is dark on the day. What each arm was measured to be worth:\n"
            + "".join(f"    {k}: {ARM_BOARD_DELTA[k]}\n" for k in sorted(env_arms))
            + "  Do not submit this package. DEF-RUNBOOK section 2 step 9.")
    return {key: str(got) for key, _, got, _ in rows}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def free_gib() -> float:
    out = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
    pages = {}
    for line in out.splitlines():
        m = re.match(r"Pages (free|inactive|speculative|purgeable):\s+(\d+)\.", line)
        if m:
            pages[m.group(1)] = int(m.group(2))
    return sum(pages.values()) * 16384 / 2**30


# The builder and the operator fit, as they appear on a command line: the bare file name,
# preceded by a path separator, a quote, or the start of an argument -- NOT by a word character,
# so `test_make_submission_pm_cd4.py` (the test suite) does not read as a build.
BUILDER_MARKER = re.compile(r"(?<![\w-])(make_submission_pm_cd4|build_union_whitener)\.py\b")
# Executables that can be RUNNING the builder: an interpreter (any flags), the `uv` shim before
# and after it execs the interpreter, and the wrappers the runbook and the rehearsals put in
# front of them. Shells are deliberately absent: a shell whose command line mentions the builder
# is a launcher, not a build, and the first rehearsal tripped on its own (`bash build191.sh`).
BUILDER_LAUNCHERS = re.compile(r"(^|/)(python[0-9.]*|uv|time|caffeinate|nice|nohup|timeout)$")


def other_builders() -> list[str]:
    """Processes actually executing the builder or the operator fit (242 #8; fixed by 247).

    The pre-247 predicate required `argv[1]` to be an interpreter and the script to sit in
    `argv[2:4]`, so it missed `python -X dev -u <script>` (two interpreter flags),
    `python -c "import runpy; runpy.run_path('<script>')"` (how the test suite starts it), the
    `uv run python <script>` shim before it execs the interpreter, and `/usr/bin/time -l ...`
    -- and a miss prints nothing, which is indistinguishable from "nothing running". Now:
    the executable must be a launcher (above) and the builder's file name must appear
    ANYWHERE in the remaining arguments. Demonstrated on all six launch forms with a sleeping
    dummy in experiments/247-oct22-guards (forms247.py).
    """
    out = subprocess.run(["ps", "-axo", "pid=,args="], capture_output=True, text=True).stdout
    me = os.getpid()
    hits = []
    for ln in out.splitlines():
        parts = ln.split(None, 1)
        if len(parts) < 2 or not parts[0].isdigit() or int(parts[0]) == me:
            continue
        argv = parts[1].split()
        if not argv or not BUILDER_LAUNCHERS.search(argv[0].lstrip("-")):
            continue
        if any(BUILDER_MARKER.search(a) for a in argv[1:]):
            hits.append(ln.strip())
    return hits


def builder_literal(name: str, source: Path = ROOT / "scripts/make_submission_pm_cd4.py") -> dict:
    """A module-level dict literal out of the builder's SOURCE, never an import."""
    import ast
    for node in ast.parse(source.read_text()).body:
        targets: list[str] = []
        if isinstance(node, ast.Assign):
            targets = [t.id for t in node.targets if isinstance(t, ast.Name)]
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            targets = [node.target.id]
        if name in targets and node.value is not None:
            return ast.literal_eval(node.value)
    raise SystemExit(f"{source}: no module-level literal named {name}")


def bundle_axis(bundle: Path) -> tuple[str, int]:
    """The operator axis THIS bundle implies -- `sorted_by_official_index(union(donors) - panel)`
    as ordered symbols -- and its 235 digest, from the donors' `genes.txt` alone (axis235.py)."""
    import pandas as pd
    official = pd.read_csv(bundle / "gene_names.csv")["gene_name"].astype(str).tolist()
    panel = set(pd.read_csv(bundle / "pert_counts.csv")["target_gene"].astype(str))
    tindex = {g: i for i, g in enumerate(official)}
    union: set[str] = set()
    for donor in PRODUCTION_DONORS:
        union |= {g for g in (ROOT / "data/donors" / donor / "genes.txt").read_text().split("\n")
                  if g in tindex}
    columns = sorted(tindex[g] for g in union - panel)
    digest = hashlib.sha256("\n".join(official[i] for i in columns).encode("utf-8")).hexdigest()
    return digest, len(columns)


def verify_existing_operator(wout: Path, wmeta: Path, bundle: Path) -> str:
    """RAISE unless an operator handed to `--whitener-out` was fitted for THIS bundle (242 #18).

    The skip branch used to reuse whatever sat at `<work>/whitener_<id>.npz` with no check at
    all -- exactly the case `--whitener-out`'s help text recommends -- and the driver's own
    `sha256(wout) == ABC_OPERATOR_SHA` abort never ran. The builder catches it 20 minutes later
    (235); this catches it before the stage is skipped. Expectation, in the builder's order:
    the metadata json's own `axis_sha256` (every fit since 235 writes it), else the builder's
    `LEGACY_OPERATOR_AXIS_SHA256[operator sha256]`, else refuse.
    """
    meta = json.loads(wmeta.read_text())
    expect, n_cols = bundle_axis(bundle)
    recorded = meta.get("axis_sha256")
    source = f"metadata:{wmeta.name}"
    if not isinstance(recorded, str) or not recorded:
        recorded = builder_literal("LEGACY_OPERATOR_AXIS_SHA256").get(sha256(wout))
        source = "LEGACY_OPERATOR_AXIS_SHA256[operator sha256]"
    if not isinstance(recorded, str) or not recorded:
        raise SystemExit(
            "REFUSING TO SKIP the whitener stage: the existing operator carries no axis digest, "
            "so the panel it was fitted for cannot be checked.\n"
            f"  operator  {wout}\n  metadata  {wmeta}\n"
            f"  this bundle's axis_sha256 {expect} over {n_cols} columns\n"
            "  Re-fit it (delete the two files and re-run the stage), or point --whitener-out "
            "at an operator whose json carries 'axis_sha256'.")
    if recorded != expect:
        raise SystemExit(
            "REFUSING TO SKIP the whitener stage: the existing operator was fitted for a "
            "DIFFERENT PANEL.\n"
            f"  operator                  {wout}\n"
            f"  axis_sha256 recorded ({source}) {recorded}\n"
            f"  axis_sha256 of THIS bundle                {expect}\n"
            f"  columns                   {n_cols} (not the check: n_cols is 17,942 on every panel measured)\n"
            "  A stale operator in the work directory pairs ~98.8% of its coordinates with the "
            "wrong gene, silently (DEF-RUNBOOK section 0.1). Delete it and re-run the stage, or "
            "pass the operator fitted for this bundle.")
    if meta.get("n_cols") != n_cols:
        raise SystemExit(f"REFUSING TO SKIP the whitener stage: metadata n_cols {meta.get('n_cols')} "
                         f"!= {n_cols} recomputed from this bundle ({wmeta})")
    return source


class Runner:
    def __init__(self, work: Path, dry: bool):
        self.work, self.dry = work, dry
        self.logs = work / "logs"
        self.logs.mkdir(parents=True, exist_ok=True)
        self.timings_path = work / "timings.json"
        self.timings = json.loads(self.timings_path.read_text()) if self.timings_path.exists() else {}

    def run(self, stage: str, cmd: list[str], env_defaults: dict[str, str] | None = None,
            bare_env: bool = False) -> None:
        """`env_defaults` are DEFAULTS, never overrides (240).

        `bare_env` (437, 2026-09-23): run the stage with EVERY `VCC_*` variable removed, exports included. The pin
        stages (`repin-frozen`, `pin-artifacts`) need no arm knob, and an exported day-of value whose partner is a
        driver DEFAULT (e.g. `VCC_BULK_CTX_SCALE` needs `VCC_PERT_SCALE_BALLAST`) otherwise makes the builder refuse
        at import there: defaults never reach these stages, exports do (432 T3). 432's rehearsal ran them with nothing
        exported; this makes that the driver's behaviour whatever the operator exported.

        `resolve_env` is the single implementation of that rule: for a key the operator
        exported, `resolved[k]` IS `os.environ[k]`, so the merge below is `defaults first,
        operator on top` however you read it. The printed and logged command carries the
        RESOLVED values -- the ones the child will actually see -- and names which of them the
        operator supplied, so a run and its log can never disagree about the configuration.
        """
        assert "submit" not in cmd, "this driver never submits"
        resolved, from_operator = resolve_env(env_defaults or {})
        if bare_env:
            resolved, from_operator = {}, []
        full = ["/usr/bin/time", "-l", *cmd]
        shown = " ".join(f"{k}={v}" for k, v in resolved.items()) + (" " if resolved else "") + " ".join(map(str, cmd))
        print(f"\n=== [{stage}] {shown}" + ("   [bare env: every VCC_* removed]" if bare_env else ""), flush=True)
        if from_operator:
            print(f"    [env] operator-supplied, and it WINS: "
                  + ", ".join(f"{k}={resolved[k]!r}" for k in from_operator), flush=True)
        if self.dry:
            return
        log = self.logs / f"{stage}.log"
        t0 = time.time()
        with open(log, "w") as fh:
            fh.write(f"$ {shown}\n")
            if from_operator:
                fh.write(f"# operator-supplied (wins over the driver's default): {from_operator}\n")
            fh.write("\n")
            fh.flush()
            base_env = ({k: v for k, v in os.environ.items() if not k.startswith("VCC_")} if bare_env
                        else dict(os.environ))
            proc = subprocess.run(full, cwd=ROOT, env={**base_env, **resolved},
                                  stdout=fh, stderr=subprocess.STDOUT)
        wall = time.time() - t0
        text = log.read_text()
        m = re.search(r"(\d+)\s+maximum resident set size", text)
        rss = int(m.group(1)) / 2**30 if m else None
        self.timings[stage] = {"seconds": round(wall, 1), "peak_rss_gib": round(rss, 2) if rss else None,
                               "returncode": proc.returncode, "log": str(log.relative_to(ROOT)),
                               "command": shown}
        self.timings_path.write_text(json.dumps(self.timings, indent=2) + "\n")
        tail = "\n".join(ln for ln in text.splitlines()[-25:]
                         if not re.match(r"^\s*\d+\s+(average|page|swaps|block|messages|signals|voluntary|involuntary|instructions|cycles|peak memory)", ln))
        print(tail, flush=True)
        print(f"--- [{stage}] {wall:.1f}s, peak RSS {rss:.2f} GiB, rc={proc.returncode}" if rss
              else f"--- [{stage}] {wall:.1f}s rc={proc.returncode}", flush=True)
        if proc.returncode != 0:
            raise SystemExit(f"stage {stage} failed (rc={proc.returncode}); see {log}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", type=Path, default=ROOT / "data/challenge")
    ap.add_argument("--work", type=Path, required=True, help="everything this run writes goes here")
    ap.add_argument("--name", required=True, help="artifact stem; the .h5ad is <work>/<name>.h5ad")
    ap.add_argument("--frozen-inputs", type=Path, default=ROOT / "scripts/frozen_inputs.json",
                    help="the frozen table to re-pin and build against (rehearsal: a throwaway copy)")
    ap.add_argument("--whitener-id", default="union_pruned_def")
    ap.add_argument("--whitener-out", type=Path, default=None,
                    help="default <work>/whitener_<id>.npz; pass an existing fit to skip the stage")
    ap.add_argument("--cd4-source", type=Path, default=CD4_SOURCE)
    ap.add_argument("--cd4-out", type=Path, default=None, help="default <work>/cd4_compact")
    ap.add_argument("--cd4-thresholds", choices=["default", "relaxed"], default="default",
                    help="047's gates (275 per condition, 290 union) or none. Relax ONLY for a synthetic panel")
    ap.add_argument("--census-max-zero", default=None,
                    help="538: forward --max-zero-carrier N to def_census.py (integer 6..30). "
                         "Omitted, the census keeps its floor of 5 and this argument adds no argv word")
    ap.add_argument("--emitter", default="rc_callset",
                    help="the incumbent's emitter: rc_callset (202, and unchanged through "
                         "`cd4f034` and `h1union`, board 0.2118310, entry iUziqE1gewfcPKzYUDaB); "
                         "rc_hybrid is 191's (0.1834190)")
    ap.add_argument("--incumbent-meta", type=Path,
                    default=ROOT / "submissions/249-cd4svd5-f034.h5ad.meta.json",
                    help="the CURRENT incumbent's meta -- the config the D/E/F build must differ "
                         "from in exactly contexts + whitener_id (+ the re-fit -> pinned-basis swap). "
                         "Default is the shipped `cd4svd5` (board 0.2160944, 2026-09-05), NOT "
                         "`h1union` (0.2118310) or `rc_callset` (0.1918638), which predate arms "
                         "and would let a build that dropped them pass")
    ap.add_argument("--allow-arm-drift", default="",
                    help="comma-separated arms this run DELIBERATELY does not ship "
                         "(e.g. VCC_H1_UNION_G); forwarded to def_postbuild.py. Requires a "
                         "matching exported value and a note in the round's RESULTS.md")
    ap.add_argument("--max-rss-gib", type=float, default=24.0, help="builder guard, as build191.sh")
    ap.add_argument("--min-free-gib", type=float, default=26.0)
    ap.add_argument("--stages", default=None, help="comma-separated subset, in pipeline order")
    ap.add_argument("--from", dest="from_stage", default=None)
    ap.add_argument("--skip-source-hash", action="store_true", help="skip the 30 s sha256 of the CD4 source")
    ap.add_argument("--restore-builder-after", action="store_true",
                    help="REHEARSAL ONLY: byte-restore the builder from git HEAD after `prep`")
    ap.add_argument("--dry-run", action="store_true", help="print every command, run nothing")
    args = ap.parse_args()

    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    bundle = args.bundle.resolve()
    manifest = json.loads((bundle / "manifest.json").read_text())
    contexts = list(manifest["contexts"])
    ctx = ",".join(contexts)
    is_abc = contexts == ["A", "B", "C"]
    frozen = args.frozen_inputs.resolve()
    wout = (args.whitener_out or work / f"whitener_{args.whitener_id}.npz").resolve()
    wmeta = wout.with_suffix(".json")
    cd4_out = (args.cd4_out or work / "cd4_compact").resolve()
    key_out = work / "pert_scale_def_lo0.3_hi3.npz"   # 432: 410's key re-fit for this panel (non-A/B/C only)
    out_h5ad = work / f"{args.name}.h5ad"
    out_meta = Path(str(out_h5ad) + ".meta.json")
    py = ["uv", "run", "python"]

    stages = STAGES if not args.stages else [s.strip() for s in args.stages.split(",")]
    for s in stages:
        if s not in STAGES:
            ap.error(f"unknown stage {s}; order is {STAGES}")
    if args.from_stage:
        stages = [s for s in STAGES if STAGES.index(s) >= STAGES.index(args.from_stage)]
    r = Runner(work, args.dry_run)
    defaults = day_env(is_abc)
    resolved, from_operator = resolve_env(defaults)
    union_sources = os.environ.get("VCC_UNION_SOURCES", "").strip() or "<unset: the H1 pair " \
        f"{SHIPPED_UNION_SOURCES}>"
    print(f"bundle {bundle}  panel {manifest.get('panel_id')}  contexts {contexts}\n"
          f"work {work}\nfrozen table {frozen}\nwhitener {wout} (id {args.whitener_id})\n"
          f"cd4 banks {cd4_out}\nartifact {out_h5ad}\nstages {stages}", flush=True)
    print("\nthe configuration this run will build (driver defaults; an exported value wins):")
    for key, value in resolved.items():
        origin = "OPERATOR" if key in from_operator else "driver default"
        drift = "" if value == defaults[key] else f"   <-- NOT the shipped {defaults[key]!r}"
        print(f"  {key:<20} {value!r:<10} ({origin}){drift}")
    if is_abc:
        print(f"  {'VCC_CD4_SIZING':<20} <unset on A/B/C: the legacy all-rows median, so 240 row 1 "
              "stays byte-for-byte>")
        print(f"  {'VCC_CD4_DEFLATE_BASIS':<20} <unset on A/B/C: the re-fit above IS the incumbent's; "
              "the pinned basis is for other panels>")
    print(f"  {'VCC_UNION_SOURCES':<20} {union_sources}")
    print(f"  emitter              {args.emitter!r}\n  flags                {' '.join(SHIPPED_FLAGS)}")
    drifted = [k for k, v in resolved.items() if v != defaults[k]]
    if drifted:
        print("\n  *** ARM DRIFT: the configuration above is NOT the shipped one. The arms and what\n"
              "      the board paid for them:\n"
              + "".join(f"        {k}: shipped {defaults[k]!r}, this run {resolved[k]!r} -- "
                        f"{ARM_BOARD_DELTA[k]}\n" for k in drifted)
              + "      Deliberate? then `def_postbuild.py --allow-arm-drift "
              + ",".join(sorted(drifted)) + "` must be passed too, and the round's RESULTS.md\n"
              "      must say why. Otherwise stop and unset it. ***", flush=True)
    t_all = time.time()

    # 437 (2026-09-23): the restore used to live inside the `prep` stage, so a run whose stage list stops before
    # prep (every arm build) or that fails mid-way left the production builder PINNED. It now runs in `finally`
    # whenever --restore-builder-after is set and a pin record exists.
    try:
        for stage in stages:
            if stage == "preflight":
                r.run(stage, [*py, "scripts/def_preflight.py", "--bundle", str(bundle),
                              *(["--verify"] if is_abc else []), "--out", str(work / "preflight.json")])
            elif stage == "identity":
                r.run(stage, [*py, "scripts/def_identity.py", "--bundle", str(bundle),
                              *(["--verify"] if is_abc else ["--contexts", ctx]),
                              "--out", str(work / "identity.json")])
            elif stage == "census":
                cd4_root = cd4_out if (cd4_out / "cd4_rest" / "perts.txt").exists() else None
                if cd4_root is None and not is_abc:
                    # 242 #10: on a new panel the census must count against the banks extracted
                    # FOR that panel; counting against the A/B/C banks (the script's default)
                    # measures the wrong object. The census's own coverage gate is fail-closed
                    # either way; this makes the reason for a failure the right one.
                    raise SystemExit(
                        f"census on a non-A/B/C bundle needs the CD4 banks re-extracted for this panel "
                        f"at {cd4_out}; run the cd4-extract stage first (DEF-RUNBOOK section 2 step 3)")
                zero_arg: list[str] = []
                if args.census_max_zero is not None:
                    tok = str(args.census_max_zero)
                    if not (tok.isdigit() and not (len(tok) > 1 and tok[0] == "0") and 6 <= int(tok) <= 30):
                        raise SystemExit(f"--census-max-zero must be an integer in [6, 30] (got {args.census_max_zero})")
                    zero_arg = ["--max-zero-carrier", str(int(tok))]
                r.run(stage, [*py, "scripts/def_census.py", "--bundle", str(bundle),
                              *(["--verify"] if is_abc else []),
                              *(["--cd4-root", str(cd4_root)] if cd4_root else []),
                              *zero_arg,
                              "--out", str(work / "census.json")])
            elif stage == "cd4-extract":
                src = args.cd4_source.resolve()
                if not src.is_file() or src.stat().st_size != CD4_SOURCE_BYTES:
                    raise SystemExit(f"CD4 source missing or wrong size: {src} (expect {CD4_SOURCE_BYTES} bytes)")
                if not args.skip_source_hash and not args.dry_run:
                    t0 = time.time()
                    got = sha256(src)
                    if got != CD4_SOURCE_SHA256:
                        raise SystemExit(f"CD4 source sha256 {got} != pinned {CD4_SOURCE_SHA256}")
                    print(f"CD4 source sha256 verified in {time.time() - t0:.0f}s", flush=True)
                thr = [] if args.cd4_thresholds == "default" else ["--min-targets", "0", "--min-union-targets", "0"]
                r.run(stage, [*py, "experiments/047-cd4-direct/extract_cd4.py", "--input", str(src),
                              "--output", str(cd4_out), "--panel-csv", str(bundle / "pert_counts.csv"),
                              "--genes-csv", str(bundle / "gene_names.csv"), "--overwrite", *thr])
            elif stage == "whitener":
                if wout.exists() and wmeta.exists():
                    # 242 #18: never reuse an operator without proving it was fitted for THIS bundle.
                    source = verify_existing_operator(wout, wmeta, bundle)
                    print(f"\n=== [whitener] exists and its axis matches this bundle ({source}); "
                          f"skipping the fit: {wout}", flush=True)
                    continue
                r.run(stage, [*py, "scripts/build_union_whitener.py", "--bundle", str(bundle),
                              "--basis-prune", "essential", "--out", str(wout),
                              "--scratch", str(work / ".whitener_scratch"), "--max-rss-gib", "14"])
                if not args.dry_run and not is_abc and sha256(wout) == ABC_OPERATOR_SHA:
                    raise SystemExit("the re-fit operator is byte-identical to the A/B/C one: the panel did not change")
                if not args.dry_run:
                    # The fresh fit must carry this bundle's axis too (its json is what step 9 reads).
                    verify_existing_operator(wout, wmeta, bundle)
            elif stage == "key":
                # 432: 410's per-target key is panel- AND operator-specific (DEF-RUNBOOK 410 addition). Same arithmetic
                # as build_weights410.py; its --ka mode reproduces the pinned A/B/C artifact d7f76a22 byte for byte.
                if is_abc:
                    print("\n=== [key] A/B/C: the pinned wshare_lo0.3_hi3 artifact IS the key; nothing to fit", flush=True)
                    continue
                r.run(stage, [*py, "scripts/def_pert_scale_key.py", "--bundle", str(bundle),
                              "--op", str(wout), "--lo", "0.3", "--hi", "3", "--out", str(key_out)])
            elif stage == "repin-frozen":
                if not frozen.exists():
                    shutil.copyfile(ROOT / "scripts/frozen_inputs.json", frozen)
                    print(f"seeded throwaway table {frozen} from scripts/frozen_inputs.json", flush=True)
                r.run(stage, [*py, "scripts/make_submission_pm_cd4.py", "--repin-frozen-inputs",
                              "--bundle", str(bundle), "--frozen-inputs", str(frozen)], bare_env=True)
            elif stage == "pin-artifacts":
                r.run(stage, [*py, "scripts/def_pin_panel_artifacts.py", "--whitener", str(wout),
                              "--whitener-metadata", str(wmeta), "--whitener-id", args.whitener_id,
                              "--cd4-root", str(cd4_out), "--apply", "--out", str(work / "pin-artifacts.json"),
                              *([] if is_abc else ["--pert-scale", str(key_out)])], bare_env=True)
            elif stage in ("repin-panel-support", "build"):
                if not args.dry_run:
                    others = other_builders()
                    if others:
                        raise SystemExit("another builder is running; refusing to start a ~19 GiB stage:\n  " + "\n  ".join(others))
                    fg = free_gib()
                    if fg < args.min_free_gib:
                        raise SystemExit(f"free+inactive memory {fg:.1f} GiB < {args.min_free_gib} GiB; refusing")
                    print(f"memory check: {fg:.1f} GiB free+inactive, no other builder", flush=True)
                common = [*py, "scripts/make_submission_pm_cd4.py", "--bundle", str(bundle),
                          "--frozen-inputs", str(frozen), "--expect-contexts", ctx,
                          "--whitener", str(wout), "--whitener-metadata", str(wmeta),
                          "--whitener-id", args.whitener_id, "--cd4-root", str(cd4_out),
                          "--max-rss-gib", str(args.max_rss_gib)]
                if stage == "repin-panel-support":
                    r.run(stage, [*common, "--repin-panel-support", "--min-cpm", "10"], defaults)
                else:
                    if out_h5ad.exists():
                        raise SystemExit(f"{out_h5ad} exists; the builder refuses to overwrite. Pick a new --name")
                    r.run(stage, [*common, *SHIPPED_FLAGS, "--emitter", args.emitter, "-o", str(out_h5ad)],
                          defaults)
                    # 240's guard: re-read the arms off the finished artifact. An assertion that has
                    # never been observed to fire is not a guard -- this one is demonstrated firing
                    # in experiments/240-driver-fix (row 4).
                    # 485: VCC_BULK_DUMP writes the npz and SystemExit(0)s BEFORE write_submission, so
                    # there is no meta for the arm assert. Requiring the assert here makes every
                    # successful dump exit 1, and def_day.sh then keeps whatever npz is on disk.
                    # Skip the assert only when a dump was requested and no h5ad was published.
                    if not args.dry_run:
                        dump = os.environ.get("VCC_BULK_DUMP", "").strip()
                        if dump and not out_h5ad.exists():
                            dump_path = Path(dump)
                            if not dump_path.is_file() or dump_path.stat().st_size == 0:
                                raise SystemExit(
                                    f"dump stage finished without a non-empty dump file: {dump_path}")
                            print(f"[dump] {dump_path} is non-empty and no h5ad was written; "
                                  "skipping the arm assert (the dump exits before meta is written)",
                                  flush=True)
                        else:
                            assert_build_carries_arms(out_meta, resolved)
            elif stage == "postbuild":
                r.run(stage, [*py, "scripts/def_postbuild.py", "--meta", str(out_meta), "--bundle", str(bundle),
                              "--incumbent", str(args.incumbent_meta.resolve()), "--cd4-root", str(cd4_out),
                              "--allow-config-diff", "whitener_id" if is_abc else
                              "whitener_id,pert_scale,pert_scale_sha256,bulk_load_sha256,bulk_ctx_scale",
                              *(["--allow-arm-drift", args.allow_arm_drift] if args.allow_arm_drift else []),
                              "--out", str(work / "postbuild.json")])
            elif stage == "prep":
                if not args.dry_run:
                    others = other_builders()
                    if others:
                        raise SystemExit("another builder is running; refusing to start vcc prep (~22 GiB):\n  " + "\n  ".join(others))
                    fg = free_gib()
                    if fg < args.min_free_gib:
                        raise SystemExit(f"free+inactive memory {fg:.1f} GiB < {args.min_free_gib} GiB; refusing")
                r.run(stage, ["vcc", "prep", str(out_h5ad), "-g", str(bundle / "gene_names.csv"),
                              "--perts", str(bundle / "pert_counts.csv"), "--contexts", ctx,
                              "-o", str(work / f"{args.name}.vcc"), "--dry-run"])
    finally:
        rec_path = work / "pin-artifacts.json"
        if args.restore_builder_after and not args.dry_run and rec_path.exists():
            # Restore only if the builder on disk is still THIS pin's output (`builder_sha256_after`): a stale record
            # from an earlier run in the same work dir must never revert a later, committed builder.
            now = sha256(ROOT / "scripts/make_submission_pm_cd4.py")
            rec = json.loads(rec_path.read_text())
            if now == rec.get("builder_sha256_after"):
                r.run("restore-builder", [*py, "scripts/def_pin_panel_artifacts.py", "--restore",
                                          "--record", str(rec_path)])
            else:
                print(f"[restore-builder] builder on disk {now[:16]}... is not this pin's output "
                      f"{str(rec.get('builder_sha256_after'))[:16]}...; not restoring", flush=True)

    total = time.time() - t_all
    if not args.dry_run:
        r.timings["_total"] = {"seconds": round(total, 1), "stages": stages}
        r.timings_path.write_text(json.dumps(r.timings, indent=2) + "\n")
        print("\nstage                   seconds   peak RSS GiB")
        for s, t in r.timings.items():
            if s.startswith("_"):
                continue
            print(f"{s:<22} {t['seconds']:>9.1f}   {t.get('peak_rss_gib') or 0:>6.2f}")
        print(f"{'TOTAL (this run)':<22} {total:>9.1f}")
        print(f"timings -> {r.timings_path}")


if __name__ == "__main__":
    main()
