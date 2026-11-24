"""432: write the PROPOSED `scripts/def_pipeline.py` (never applied here) and its unified diff.

The production driver's defaults still describe early September (`h1union` + cd4svd5). Run on Oct 22 as DEF-RUNBOOK
section 2 prints it, it would silently drop every knob added since 319 (normrel cap, partial CD4 rows, the CLEAN4 prior,
410's key + ballast, sharpening, the MDE, the zero repair) -- and its arms guard and default postbuild incumbent
(`249-cd4svd5`) would both pass the package, the 239 failure mode again.

This script makes four edits, each anchored on text that must occur exactly once:
  1. SHIPPED_ENV gains the panel-agnostic September knobs (433h's config minus the day-of LOAD/CTX);
  2. ABC_ENV / DEF_ENV gain the panel-specific key id (`wshare_lo0.3_hi3` on A/B/C, the day's `wshare_lo0.3_hi3_def`
     registered by the proposed `def_pin_panel_artifacts.py --pert-scale`);
  3. ARM_BOARD_DELTA gains a line per new key (the driver indexes it for every drifted key: KeyError otherwise);
  4. assert_build_carries_arms re-reads the new knobs off the written meta, plus VCC_BULK_LOAD / VCC_BULK_CTX_SCALE when
     exported (day-of values, never defaults), and the postbuild stage allows the panel-specific artifact keys.
usage: uv run python make_def_pipeline_proposal.py     (writes def_pipeline.proposed.py + def_pipeline.proposed.diff)
"""
from __future__ import annotations

import difflib
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "scripts/def_pipeline.py"

EDITS = [
    (  # 1
        'SHIPPED_ENV = {"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.34", "VCC_H1_UNION_G": "0.28"}\n',
        'SHIPPED_ENV = {"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.34", "VCC_H1_UNION_G": "0.28",\n'
        '               # 432 (2026-09-23): the September knobs that need no panel artifact -- 433h\'s config\n'
        '               # (417e stack LOAD + coarse MDE + zero repair + re-level) minus the two DAY-OF values\n'
        '               # (VCC_BULK_LOAD, VCC_BULK_CTX_SCALE), which are exported on the day, never defaulted.\n'
        '               "VCC_NORMREL_EXP": "-1.5", "VCC_NORMREL_MODE": "cap", "VCC_CD4_PARTIAL_G": "0.5",\n'
        '               "VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_PERT_SCALE_BALLAST": "1", "VCC_BULK_SHARPEN_Q": "1.25",\n'
        '               "VCC_DE_MDE": "0.1,0.1,0.3,coarse", "VCC_DE_ZERO_REPAIR": "-5"}\n',
    ),
    (  # 2a
        '           "VCC_CD4_DEFLATE_BASIS": "V5_abc"}\n',
        '           "VCC_CD4_DEFLATE_BASIS": "V5_abc",\n'
        '           # 432: 410\'s key RE-FIT for this panel on the day\'s operator (key432.py), registered under\n'
        '           # this id by the proposed `def_pin_panel_artifacts.py --pert-scale` (DEF-RUNBOOK 410 addition).\n'
        '           "VCC_PERT_SCALE": "wshare_lo0.3_hi3_def"}\n',
    ),
    (  # 2b
        'ABC_ENV = {"VCC_CD4_DEFLATE_K": "5"}\n',
        'ABC_ENV = {"VCC_CD4_DEFLATE_K": "5", "VCC_PERT_SCALE": "wshare_lo0.3_hi3"}\n',
    ),
    (  # 3
        '                   "VCC_CD4_DEFLATE_BASIS": "+0.0042634 (`cd4svd5`) carried to a new panel as the PINNED "\n'
        '                                            "A/B/C basis (251); never the re-fit off A/B/C (250)"}\n',
        '                   "VCC_CD4_DEFLATE_BASIS": "+0.0042634 (`cd4svd5`) carried to a new panel as the PINNED "\n'
        '                                            "A/B/C basis (251); never the re-fit off A/B/C (250)",\n'
        '                   "VCC_CD4_PARTIAL_G": "partial CD4 rows, board steps 319p/319p2 (+0.0011593 / +0.0002969, FINDINGS s0)",\n'
        '                   "VCC_DE_SIGN_PRIOR": "CLEAN4 sign prior q80 (392 +0.0025062, 395 +0.0013142)",\n'
        '                   "VCC_PERT_SCALE": "410\'s per-target bulk key (410b +0.0170570 with ballast); the artifact is "\n'
        '                                     "panel- and operator-specific: re-fit on the day (432)",\n'
        '                   "VCC_PERT_SCALE_BALLAST": "410b ballast cells: +0.0170570 over the welded 410",\n'
        '                   "VCC_BULK_SHARPEN_Q": "417a bulk sharpening q 1.25: +0.0069612",\n'
        '                   "VCC_DE_MDE": "426/433 coarse MDE (H1 nmae -0.0095; board read pending)",\n'
        '                   "VCC_DE_ZERO_REPAIR": "428 zero-gene repair (cannot hurt; A/B/C 24 genes at lfc ~ -35)"}\n',
    ),
    (  # 4a: the generic September-knob rows, just before the verdict
        '    bad = [r for r in rows if not r[3]]\n'
        '    lines = [f"  {\'OK \' if ok else \'NO \'} config.{key}: asked {want!r}, the artifact carries {got!r}"\n',
        '    rows.extend(september_rows(cfg, resolved))\n'
        '    bad = [r for r in rows if not r[3]]\n'
        '    lines = [f"  {\'OK \' if ok else \'NO \'} config.{key}: asked {want!r}, the artifact carries {got!r}"\n',
    ),
    (  # 4b: the helper itself, before assert_build_carries_arms
        'def assert_build_carries_arms(meta_path: Path, resolved: dict[str, str]) -> dict[str, str]:\n',
        'def september_rows(cfg: dict, resolved: dict[str, str]) -> list[tuple[str, object, object, bool]]:\n'
        '    """432: the September knobs, re-read off the written config by value. A knob the operator exported\n'
        '    EMPTY (the builder\'s "off") must be absent from the config. VCC_BULK_LOAD / VCC_BULK_CTX_SCALE are day-of\n'
        '    values: checked only when exported, never defaulted."""\n'
        '    def val(key: str) -> str:\n'
        '        return (resolved.get(key, os.environ.get(key, "")) or "").strip()\n'
        '    rows: list[tuple[str, object, object, bool]] = []\n'
        '    def scalar(env: str, key: str, conv) -> None:\n'
        '        v = val(env)\n'
        '        want = conv(v) if v else "<absent>"\n'
        '        got = cfg.get(key, "<absent>")\n'
        '        rows.append((key, want, got, got == want))\n'
        '    scalar("VCC_NORMREL_EXP", "normrel_exp", float)\n'
        '    if val("VCC_NORMREL_EXP"):   # the builder always writes normrel_mode alongside the exponent\n'
        '        scalar("VCC_NORMREL_MODE", "normrel_mode", lambda s: s.lower())\n'
        '    scalar("VCC_CD4_PARTIAL_G", "cd4_partial_g", float)\n'
        '    scalar("VCC_DE_SIGN_PRIOR", "de_sign_prior", lambda s: s.lower())\n'
        '    scalar("VCC_PERT_SCALE", "pert_scale", lambda s: s.lower())\n'
        '    scalar("VCC_BULK_SHARPEN_Q", "bulk_sharpen_q", float)\n'
        '    b = val("VCC_PERT_SCALE_BALLAST")\n'
        '    rows.append(("pert_scale_ballast", "<present>" if b else "<absent>", cfg.get("pert_scale_ballast", "<absent>"),\n'
        '                 ("pert_scale_ballast" in cfg) == bool(b)))\n'
        '    m = val("VCC_DE_MDE")\n'
        '    got = cfg.get("de_mde", "<absent>")\n'
        '    if m:\n'
        '        p = m.split(",")\n'
        '        want = {"c": float(p[0]), "cap": float(p[1]), "min_count": float(p[2]), "coarse": len(p) == 4}\n'
        '        # the builder writes "coarse": True only for the coarse form (426a\'s plain meta has no key)\n'
        '        ok = isinstance(got, dict) and all(got.get(k, False if k == "coarse" else None) == v for k, v in want.items())\n'
        '    else:\n'
        '        want, ok = "<absent>", got == "<absent>"\n'
        '    rows.append(("de_mde", want, got, ok))\n'
        '    z = val("VCC_DE_ZERO_REPAIR")\n'
        '    got = cfg.get("de_zero_repair", "<absent>")\n'
        '    rows.append(("de_zero_repair.floor_log2", float(z) if z else "<absent>",\n'
        '                 got.get("floor_log2") if isinstance(got, dict) else got,\n'
        '                 (isinstance(got, dict) and got.get("floor_log2") == float(z)) if z else got == "<absent>"))\n'
        '    g = val("VCC_DE_MAG")\n'
        '    got = cfg.get("de_mag", "<absent>")\n'
        '    rows.append(("de_mag.variant", g or "<absent>", got.get("variant") if isinstance(got, dict) else got,\n'
        '                 (isinstance(got, dict) and got.get("variant") == g) if g else got == "<absent>"))\n'
        '    if os.environ.get("VCC_BULK_LOAD", "").strip():\n'
        '        want = os.environ["VCC_BULK_LOAD"].strip().rsplit(":", 1)[1].lower()\n'
        '        rows.append(("bulk_load_sha256", want, cfg.get("bulk_load_sha256", "<absent>"),\n'
        '                     cfg.get("bulk_load_sha256") == want))\n'
        '    if os.environ.get("VCC_BULK_CTX_SCALE", "").strip():\n'
        '        want = {k.strip(): float(v) for k, v in (p.split("=") for p in os.environ["VCC_BULK_CTX_SCALE"].split(","))}\n'
        '        rows.append(("bulk_ctx_scale", want, cfg.get("bulk_ctx_scale", "<absent>"), cfg.get("bulk_ctx_scale") == want))\n'
        '    return rows\n'
        '\n'
        '\n'
        'def assert_build_carries_arms(meta_path: Path, resolved: dict[str, str]) -> dict[str, str]:\n',
    ),
    (  # 4c: postbuild on a new panel: the panel-specific artifact keys differ from any A/B/C incumbent by construction
        '                          "--allow-config-diff", "whitener_id",\n',
        '                          "--allow-config-diff", "whitener_id" if is_abc else\n'
        '                          "whitener_id,pert_scale,pert_scale_sha256,bulk_load_sha256,bulk_ctx_scale",\n',
    ),
]


def main() -> None:
    src = SRC.read_text()
    out = src
    for i, (old, new) in enumerate(EDITS):
        n = out.count(old)
        if n != 1:
            raise SystemExit(f"edit {i}: anchor found {n} times")
        out = out.replace(old, new, 1)
    (HERE / "def_pipeline.proposed.py").write_text(out)
    diff = "".join(difflib.unified_diff(src.splitlines(keepends=True), out.splitlines(keepends=True),
                                        fromfile="a/scripts/def_pipeline.py", tofile="b/scripts/def_pipeline.py"))
    (HERE / "def_pipeline.proposed.diff").write_text(diff)
    print(f"wrote def_pipeline.proposed.py ({len(out)} bytes) and def_pipeline.proposed.diff ({diff.count(chr(10))} lines)")


if __name__ == "__main__":
    main()
