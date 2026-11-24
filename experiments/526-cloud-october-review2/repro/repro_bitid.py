#!/usr/bin/env python3
"""Bit-identity of the knobs added since 2d1fc214, read off the source.

Unset means the branch is not entered. The one config key that moves while its knob
is unset is de_sign_emitter_sha256: the constant changed in 8960c102 and is written
whenever VCC_DE_SIGN_PRIOR is set (it always is on this path). The file hash matches
the new constant. Exits 0.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BUILDER = (ROOT / "scripts/make_submission_pm_cd4.py").read_text()
EMITTER = ROOT / "src/vcc2026/rc_emit_k.py"
DEF_DAY = (ROOT / "scripts/def_day.sh").read_text()
INC = ROOT / "submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"


def main() -> int:
    pin = re.search(r'^DE_SIGN_EMITTER_SHA256 = "([0-9a-f]{64})"', BUILDER, re.M).group(1)
    file_sha = hashlib.sha256(EMITTER.read_bytes()).hexdigest()
    inc = json.loads(INC.read_text())["config"]["de_sign_emitter_sha256"]
    print(f"emitter file {file_sha}")
    print(f"builder pin  {pin}")
    print(f"446b meta    {inc}")
    assert file_sha == pin
    assert inc.startswith("edb5b8dc") and inc != pin

    # The allow line fires exactly when those two differ. Today they do.
    allow = "de_sign_emitter_sha256" if (pin and inc and pin != inc) else ""
    print(f"def_day would add to --allow-config-diff: {allow or '<nothing>'}")
    assert allow == "de_sign_emitter_sha256"
    assert "de_sign_emitter_sha256" in DEF_DAY
    # postbuild then allows the key by name. A third value is not compared to the pin.
    post = (ROOT / "scripts/def_postbuild.py").read_text()
    assert "not unexpected" in post
    assert "de_sign_emitter_sha256" not in post  # no pin equality there

    # Unset MDE sign: the config splat contributes nothing. No RNG in that branch.
    assert '**({} if DE_MDE_SIGN is None else {"de_mde_sign": dict(DE_MDE_SIGN)})' in BUILDER
    emit = EMITTER.read_text()
    assert 'if "mde_add_sign" in lfc_target:' in emit
    # The sign override is inside that if. Unset never enters it.
    start = emit.index('if "mde_add_sign" in lfc_target:')
    window = emit[start:start + 600]
    assert "asg = np.asarray" in window

    # Unset DE_MAG: same splat shape, and the scorer-control block is inside `if DE_MAG is not None`.
    assert '**({} if DE_MAG is None else {"de_mag":' in BUILDER

    # Call-set CPM: assignment stays 2200 unless the env is non-empty. The key is unconditional.
    assert re.search(r"^RC_CALLSET_CPM = 2200$", BUILDER, re.M)
    assert 'if os.environ.get("VCC_RC_CALLSET_CPM", "").strip():' in BUILDER
    assert '"rc_callset_cpm": RC_CALLSET_CPM' in BUILDER
    # def_day wipes every VCC_* and never exports this one, so a leaked shell value cannot ride along.
    assert "VCC_RC_CALLSET_CPM" not in DEF_DAY
    assert "unset" in DEF_DAY and "VCC_" in DEF_DAY

    # Scalar sharpen does not call the tercile function. The call is guarded.
    assert "if BULK_SHARPEN_PT is not None else None" in BUILDER
    print("BIT-IDENTITY: unset MDE sign, unset DE mag and unset call-set CPM do not enter their branches")
    print("CONFIG KEY: de_sign_emitter_sha256 moves on every build (pin != 446b) even with MDE sign unset")
    return 0


if __name__ == "__main__":
    sys.exit(main())
