"""551 -- the operator-deletion step of drive551.panel and r7a551.panel resumes after an interrupted run (555 finding 5).

Stubbed: every heavy step and fit is already 'done' in a scratch state file; comp jsons, operator sidecar jsons and ONE of the
four operator npz files exist (three were 'deleted' by an interrupted run). The panel must record all four (bytes None for the
three) and delete the remaining one. Runs in work/test-resume555/ (removed after); touches no R6 or R7 panel directory.

    uv run python experiments/551-r6-c10b/test_resume555.py
"""
import sys

sys.dont_write_bytecode = True
import guard551  # noqa: E402,F401
import json  # noqa: E402
import shutil  # noqa: E402
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
T = HERE / "work" / "test-resume555"


def fake_panel_tree(work: Path, p: str, tags, cancer, config, repro):
    pw = work / p
    ops = pw / "ops"
    ops.mkdir(parents=True, exist_ok=True)
    for i, tag in enumerate(tags):
        (ops / f"{p}__{tag}.json").write_text(json.dumps({"sha256": f"{i:064x}"}))
        if i == 0:
            (ops / f"{p}__{tag}.npz").write_bytes(b"x" * 1000)
    for f in cancer:
        c = {"peak_rss_gib": 1.0, "op_axis_checked": True, "config": config}
        if repro:
            c["ka"] = {"repro": {"worst": 0.0}}
        (pw / f"comp551-{f}.json").write_text(json.dumps(c))
    return pw


def rm_tree(d: Path):
    """bottom-up with absolute paths: shutil.rmtree's dir_fd calls trip guard551 (r7a551.ka_cd4 does the same)"""
    if d.exists():
        for f in sorted(d.rglob("*"), key=lambda x: len(x.parts), reverse=True):
            f.unlink() if f.is_file() else f.rmdir()
        d.rmdir()


def run(mod, p, state_steps, config, repro, sub):
    w = T / sub
    w.mkdir(parents=True)
    mod.WORK = w
    mod.STATE = w / "state.json"
    if hasattr(mod, "R7A"):
        mod.R7A = w
    fake_panel_tree(w, p, list(mod.FITS), mod.CANCER, config, repro)
    mod.STATE.write_text(json.dumps({s: {"done": "x", "rc": 0} for s in state_steps(p, mod)}))
    mod.heavy = lambda step, cmd, log: {"rc": 0}
    return mod


def check(mod, p):
    s = json.loads(mod.STATE.read_text())[f"{p}:ops_deleted"]
    assert set(s) - {"done"} == set(mod.FITS) and sum(r["bytes"] is None for t, r in s.items() if t != "done") == 3, s
    assert not list((mod.WORK / p / "ops").glob("*.npz")), "an operator was left"


def main():
    import drive551 as D1
    steps = lambda p, m: [f"{p}:fit:{t}" for t in m.FITS] + [f"{p}:full551", f"{p}:comp551"]
    rm_tree(T)
    run(D1, "S1", steps, "M6+P2g", True, "d551")
    D1.panel("S1", "M6+P2g")
    check(D1, "S1")
    print("drive551: interrupted deletion resumes (3 already gone -> bytes None, the last one deleted): PASS")
    import r7a551 as R
    steps_r = lambda p, m: [f"{p}:cd4"] + [f"{p}:fit:{t}" for t in m.FITS] + [f"{p}:axes", f"{p}:full551", f"{p}:comp551"]
    run(R, "T1", steps_r, "M6+P2g", False, "r7a")
    R.DRAWS = R.DRAWS                                   # the real draws.json (read-only): T1's bundle path
    R.panel("T1", "M6+P2g")
    check(R, "T1")
    print("r7a551: interrupted deletion resumes: PASS")
    rm_tree(T)
    print("ALL PASS")


if __name__ == "__main__":
    main()
