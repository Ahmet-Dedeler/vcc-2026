"""Small helpers for scripts/def_day.sh, kept out of inline `python3 -c` strings (bash brace-expands `{a, b}` in them:
452's first rehearsal logged an empty rho line three times because of it).

usage:
  def_day_helpers.py rho <build.ntclevel.json>            -> "D 0.6861 E 0.6206 F 0.5911"
  def_day_helpers.py past_cliff <build.ntclevel.json>     -> 0 if any rho > 1, 1 if all safe, 2 if unreadable
  def_day_helpers.py same <D=..,E=..> <D=..,E=..>         -> exit 0 if the factor strings agree within 5e-5
  def_day_helpers.py inbound <D=..,E=..,F=..>             -> exit 0 if every factor is inside the builder bound [0.7, 1.4]
  def_day_helpers.py nctx <sameline.npz>                  -> the number of contexts the same-line artifact switches on
  def_day_helpers.py inband <build.ntclevel.json>         -> 0 if every rho is finite, <= 1, and within 1% of rho*
  def_day_helpers.py best_safe <json:factors> [...]       -> stem of the nearest safe build (rho* from the same formula as relevel_iter.py)
  def_day_helpers.py meta_matches <meta.json> <factors> <q | tercile:q,q,q | stercile:q,q,q>
        -> 0 if the artifact's bulk_sharpen_q matches and bulk_ctx_scale matches (absent only when every factor is 1)

rho* is exp(2 * sd * Phi^-1(q)) with q = 0.31 / (0.31 + 1.78) and sd = 0.01. That formula is copied from
experiments/448-def-rehearsal-q175/relevel_iter.py (the q and sd literals). Keep the two in step.
"""
import json
import math
import sys


def factors(s: str) -> dict[str, float]:
    return {k: float(v) for k, v in (p.split("=") for p in s.split(","))}


def rho_star(sd: float = 0.01) -> float:
    """Match relevel_iter.py: q = 0.31 / (0.31 + 1.78), rho* = exp(2 * sd * Phi^-1(q))."""
    from scipy.stats import norm
    q = 0.31 / (0.31 + 1.78)
    return math.exp(2.0 * sd * float(norm.ppf(q)))


def _rhos(path: str) -> dict[str, float]:
    per = json.load(open(path))["per_context"]
    return {c: float(v["rho"]) for c, v in per.items()}


def main() -> int:
    cmd, *args = sys.argv[1:]
    if cmd == "rho":
        per = json.load(open(args[0]))["per_context"]
        print(" ".join(f"{c} {per[c]['rho']:.4f}" for c in sorted(per)))
        return 0
    if cmd == "past_cliff":
        try:
            rho = _rhos(args[0])
            if not rho or any(not math.isfinite(r) for r in rho.values()):
                print("past_cliff: rho missing or non-finite", file=sys.stderr)
                return 2
            return 0 if max(rho.values()) > 1.0 else 1
        except Exception as exc:
            print(f"past_cliff failed: {exc}", file=sys.stderr)
            return 2
    if cmd == "inband":
        try:
            rho = _rhos(args[0])
            star = rho_star()
            ok = bool(rho) and all(
                math.isfinite(r) and r > 0.0 and r <= 1.0 and abs(r / star - 1.0) <= 0.01
                for r in rho.values()
            )
            return 0 if ok else 1
        except Exception as exc:
            print(f"inband failed: {exc}", file=sys.stderr)
            return 2
    if cmd == "best_safe":
        import os
        star = rho_star()
        best = None
        for spec in args:
            path = spec.rsplit(":", 1)[0]
            try:
                rho = _rhos(path)
            except Exception as exc:
                print(f"best_safe skip {path}: {exc}", file=sys.stderr)
                continue
            h5 = path[: -len(".ntclevel.json")] + ".h5ad"
            if (not rho or not os.path.exists(h5)
                    or any((not math.isfinite(r)) or r <= 0.0 or r > 1.0 for r in rho.values())):
                continue
            gap = max(abs(math.log(r / star)) for r in rho.values())
            stem = path[: -len(".ntclevel.json")]
            if best is None or gap < best[0]:
                best = (gap, stem)
        print(best[1] if best else "")
        return 0
    if cmd == "meta_matches":
        meta_path, fact, qarg = args[0], args[1], args[2]
        try:
            cfg = json.load(open(meta_path))["config"]
            want = factors(fact)
            got_q = cfg.get("bulk_sharpen_q")
            if qarg.split(":", 1)[0] in ("tercile", "stercile"):
                # 497: the per-target form (479): the scalar q is absent and bulk_sharpen_pt carries mode and qs
                mode, qs = qarg.split(":", 1)
                want_q = [float(x) for x in qs.split(",")]
                pt = cfg.get("bulk_sharpen_pt")
                if (got_q is not None or not isinstance(pt, dict) or pt.get("mode") != mode or len(pt.get("q", [])) != len(want_q)
                        or any(abs(float(a) - b) >= 5e-5 for a, b in zip(pt["q"], want_q))):
                    print(f"per-target sharpening mismatch: meta q {got_q}, pt {pt} vs {qarg}", file=sys.stderr)
                    return 1
            else:
                q = float(qarg)
                if got_q is None or abs(float(got_q) - q) >= 5e-5:
                    print(f"q mismatch: meta {got_q} vs {q}", file=sys.stderr)
                    return 1
            got = cfg.get("bulk_ctx_scale")
            if got is None:
                if all(abs(v - 1.0) < 5e-5 for v in want.values()):
                    return 0
                print("meta has no bulk_ctx_scale but factors are not all 1", file=sys.stderr)
                return 1
            if set(map(str, got)) != set(want):
                print(f"scale keys {sorted(map(str, got))} != {sorted(want)}", file=sys.stderr)
                return 1
            for k, v in want.items():
                if abs(float(got[k]) - v) >= 5e-5:
                    print(f"scale {k}: meta {got[k]} vs {v}", file=sys.stderr)
                    return 1
            return 0
        except Exception as exc:
            print(f"meta_matches failed: {exc}", file=sys.stderr)
            return 2
    if cmd == "same":
        a, b = factors(args[0]), factors(args[1])
        return 0 if a.keys() == b.keys() and all(abs(a[k] - b[k]) < 5e-5 for k in a) else 1
    if cmd == "inbound":
        return 0 if all(0.7 <= v <= 1.4 for v in factors(args[0]).values()) else 1
    if cmd == "nctx":
        import numpy as np  # noqa: PLC0415
        print(len(np.load(args[0])["contexts"]))
        return 0
    raise SystemExit(f"unknown command {cmd!r}")


if __name__ == "__main__":
    sys.exit(main())
