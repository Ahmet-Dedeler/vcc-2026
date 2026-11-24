"""595 — run the control-dependent day-of stages on bundles of FOREIGN control cells (build595.py), on a scratch worktree
of HEAD with the 10-05 patches applied (538 → 554 → 590; shas asserted equal to 594's final table). No builds, no
def_day.sh execution, live files untouched (hashed before/after). Stage commands are exactly the driver's / day's:

  preflight  scripts/def_preflight.py --bundle B --out …           (def_pipeline.py's preflight stage on a non-A/B/C bundle)
             → if it fails on the 554 library-spread row: the same with DEF_DAY_ALLOW_WIDE_LIBS=1 (554's override)
  identity   scripts/def_identity.py --bundle B --contexts D,E,F --out …
  census     scripts/def_census.py --bundle B --cd4-root <the A/B/C-panel CD4 extraction> --out …
             (the mock bundles carry the CURRENT panel, so its 199 KA extraction is the panel's own); the override path
             is exercised with --max-zero-carrier below the observed count (must fail) and at the default (must pass)
  sameline   scripts/def_sameline.py --bundle B --identity … --out …  (def_day.sh's call; 590's floor at its default)

Usage: uv run python experiments/595-foreign-control-rehearsal/rehearse595.py <bundles dir>
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = HERE.parents[1]
os.environ.setdefault("SCRATCH594", str(Path(os.environ.get("TMP595", "/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-"
                                                              "challenge-2026/e132c515-fa77-477b-acf1-db04d2a4f831/scratchpad/r595"))))
sys.path.insert(0, str(R / "experiments/594-def-bundle-composition"))
import compose594 as C  # noqa: E402  (worktree, overlay, appliers, sha table)
sys.path.insert(0, str(R / "experiments/594-def-bundle-composition"))
from post594 import FINAL  # noqa: E402

LOG = HERE / "work/rehearse595.log"
CD4_ABC = "experiments/199-def-runbook/work/cd4_compact_abc_ka"


def log(m: str) -> None:
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {m}"
    print(line, flush=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")


def tail(t: str, n: int = 8) -> str:
    return "\n".join(t.strip().splitlines()[-n:])


def main() -> None:
    bundles = Path(sys.argv[1])
    W = HERE / "work"
    W.mkdir(parents=True, exist_ok=True)
    LOG.write_text("")
    S = C.S
    Path(os.environ["SCRATCH594"]).mkdir(parents=True, exist_ok=True)
    live_before = C.shas(C.R)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=R, capture_output=True, text=True).stdout.strip()
    if S.exists():
        subprocess.run(["git", "worktree", "remove", "--force", str(S)], cwd=R, capture_output=True)
    subprocess.run(["git", "worktree", "add", "--detach", str(S), head], cwd=R, check=True, capture_output=True)
    links = sum(C.fill(r) for r in C.OVERLAY + [CD4_ABC])
    out: dict = {"head": head, "bundles": str(bundles), "build": json.loads((bundles / "build595.json").read_text()), "runs": {}}
    try:
        for p in ("538", "554", "590"):
            rc, t = C.run(C.APPLY[p](S, False), S)
            if rc != 0:
                raise SystemExit(f"apply {p} failed: {tail(t)}")
        got = C.shas(S)
        out["patched_equals_594_final"] = got == FINAL
        log(f"worktree {head[:12]}, {links} overlays; 538+554+590 applied; shas == 594 final: {out['patched_equals_594_final']}")
        if not out["patched_equals_594_final"]:
            raise SystemExit("patched tree differs from 594's final table")
        for b in ("X", "Y"):
            B = bundles / b
            res: dict = {}
            # preflight
            rc, t = C.run([C.PY, "scripts/def_preflight.py", "--bundle", str(B), "--out", str(W / f"preflight_{b}.json")], S)
            res["preflight"] = {"rc": rc, "tail": tail(t, 12)}
            spread_lines = [ln.strip() for ln in t.splitlines() if "q75/q25" in ln or "spread" in ln.lower()]
            res["preflight"]["spread_lines"] = spread_lines
            log(f"[{b}] preflight rc {rc}; " + " | ".join(spread_lines[:6]))
            if rc != 0:
                rc2, t2 = C.run([C.PY, "scripts/def_preflight.py", "--bundle", str(B), "--out", str(W / f"preflight_{b}_allow.json")],
                                S, env={"DEF_DAY_ALLOW_WIDE_LIBS": "1"})
                res["preflight_allow_wide_libs"] = {"rc": rc2, "tail": tail(t2, 12),
                                                    "spread_lines": [ln.strip() for ln in t2.splitlines() if "q75/q25" in ln or "spread" in ln.lower()]}
                log(f"[{b}] preflight with DEF_DAY_ALLOW_WIDE_LIBS=1: rc {rc2}")
            # identity
            idj = W / f"identity_{b}.json"
            rc, t = C.run([C.PY, "scripts/def_identity.py", "--bundle", str(B), "--contexts", "D,E,F", "--out", str(idj)], S)
            res["identity"] = {"rc": rc}
            if rc == 0:
                ident = json.loads(idj.read_text())
                per = {}
                for c in "DEF":
                    per[c] = {"call": ident["call"][c], "unanimous": ident["call_unanimous_across_methods"][c],
                              "min_top1_score": round(float(ident["min_top1_score"][c]), 4),
                              "top5": [(r["sample"], round(r["score"], 3), r.get("lineage")) for r in ident["consensus"][c][:5]],
                              "per_method_top1": sorted(set(ident["top1_per_method"][c].values()))}
                    log(f"[{b}] identity {c}: {per[c]['call']} unanimous={per[c]['unanimous']} min_top1={per[c]['min_top1_score']} "
                        f"top5={[s for s, _, _ in per[c]['top5']]}")
                res["identity"]["contexts"] = per
                res["identity"]["anchor_quality"] = ident.get("method_quality")
            else:
                log(f"[{b}] identity FAILED rc {rc}: {tail(t, 4)}")
            # census
            cj = W / f"census_{b}.json"
            rc, t = C.run([C.PY, "scripts/def_census.py", "--bundle", str(B), "--cd4-root", str(S / CD4_ABC), "--out", str(cj)], S)
            gate = [ln.strip() for ln in t.splitlines() if ln.strip().startswith(("OK", "FAIL", "PASS")) or "gate" in ln.lower()]
            res["census"] = {"rc": rc, "gate_lines": gate[-8:]}
            zc = None
            if cj.exists():
                cen = json.loads(cj.read_text())
                zc = (cen.get("consensus") or {}).get("targets_with_0_donors")
                res["census"]["targets_with_0_donors"] = zc
            log(f"[{b}] census rc {rc}; targets with 0 production donors: {zc}")
            if zc is not None and zc > 0:
                rcs, ts = C.run([C.PY, "scripts/def_census.py", "--bundle", str(B), "--cd4-root", str(S / CD4_ABC),
                                 "--max-zero-carrier", str(zc - 1), "--out", str(W / f"census_{b}_strict.json")], S)
                res["census"]["override_below_count_rc"] = rcs
                log(f"[{b}] census with --max-zero-carrier {zc - 1} (below the count): rc {rcs} (must be 1)")
            # sameline, exactly as def_day.sh calls it (run even if a gate above failed, to see its behaviour)
            sl = W / f"sameline_{b}.npz"
            rcl, tl = C.run([C.PY, "scripts/def_sameline.py", "--bundle", str(B), "--identity", str(idj), "--out", str(sl)], S)
            res["sameline"] = {"rc": rcl, "lines": [ln.strip() for ln in tl.splitlines() if ln.startswith("[sameline]")]}
            if (Path(str(sl) + ".json")).exists():
                prov = json.loads(Path(str(sl) + ".json").read_text())
                res["sameline"]["contexts"] = {c: {k: v for k, v in r.items() if k in ("call", "mode", "banks", "w", "covered_targets",
                                                                                       "programme_bank", "weak_row_banks")}
                                               for c, r in prov.get("contexts", {}).items()}
            for ln in res["sameline"]["lines"]:
                log(f"[{b}] {ln}")
            out["runs"][b] = res
    finally:
        out["live_untouched"] = C.shas(C.R) == live_before
        log(f"live files untouched: {out['live_untouched']}")
        subprocess.run(["git", "worktree", "remove", "--force", str(S)], cwd=R, capture_output=True)
        (W / "rehearse595.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
        log("wrote work/rehearse595.json; worktree removed")


if __name__ == "__main__":
    main()
