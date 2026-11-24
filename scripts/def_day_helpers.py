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
  def_day_helpers.py meta_matches <meta.json> <factors> <q | tercile:q,q,q | stercile:q,q,q | tercilec:q,q,q>
  def_day_helpers.py riders_match <meta.json> <de_mag variant | -> <clean4:tau | programme:tau[:tau_c4] |
      file:npz:sha:tau[:tau_c4[:leN]] | ->   [rc=<-|N>,stack=<x0.75..x2>,rho=<-|s>,strength=<-|theta>[,lodo_hi=<h>][,callset=<prop[:swap[:steer]]>][,h1=<x1..x3>]]
      (502; 517 adds programme / file; 530's optional 4th arg checks the October meta rows; 538 adds
      optional lodo_hi and callset. "-" = absent. A 4-key spec with no lodo_hi/callset is 530's spec.)
        -> 0 if the artifact's bulk_sharpen_q matches and bulk_ctx_scale matches (absent only when every factor is 1)

rho* is exp(2 * sd * Phi^-1(q)) with q = 0.31 / (0.31 + 1.78) and sd = 0.01. That formula is copied from
experiments/448-def-rehearsal-q175/relevel_iter.py (the q and sd literals). Keep the two in step.
"""
import json
import math
import sys


def factors(s: str) -> dict[str, float]:
    return {k: float(v) for k, v in (p.split("=") for p in s.split(","))}


def rho_star(sd: float = 0.01, scale: float = 1.0) -> float:
    """Match relevel_iter.py: q = 0.31 / (0.31 + 1.78), rho* = exp(2 * sd * Phi^-1(q)).
    scale 1 (the default, and DEF_DAY_RHO_SCALE unset) returns that value with no multiply.
    """
    from scipy.stats import norm
    q = 0.31 / (0.31 + 1.78)
    base = math.exp(2.0 * sd * float(norm.ppf(q)))
    return base if scale == 1.0 else base * scale


def _rhos(path: str) -> dict[str, float]:
    per = json.load(open(path))["per_context"]
    return {c: float(v["rho"]) for c, v in per.items()}


X1_G = 0.3443655278649641  # DEF_ENV's VCC_CD4_TARGET_G; the level build records this even when the dump is x1.5


def _mde_sign_request(sign: str):
    """517: the builder's VCC_DE_MDE_SIGN grammar -> the fields its config dict must carry (None = unparseable)."""
    p = [x.strip() for x in sign.split(":")]
    try:
        if p[0] == "clean4" and len(p) == 2:
            return {"source": "clean4", "tau": float(p[1])}
        if p[0] == "programme" and len(p) in (2, 3):
            return {"source": "programme", "tau": float(p[1]), "tau_c4": float(p[2]) if len(p) == 3 else 0.05}
        if p[0] == "file" and len(p) in (4, 5, 6):
            out = {"source": "file", "sha256": p[2], "tau": float(p[3]), "tau_c4": float(p[4]) if len(p) >= 5 else 0.05}
            if len(p) == 6:
                out["le"] = float(p[5][2:])
            return out
    except (ValueError, IndexError):
        return None
    return None


# 690i: 438's dose grid; the dump-name token drops the dot (def_day.sh's case tables)
H1_TOK690I = {"x1": "h1x1", "x1.5": "h1x15", "x2": "h1x2", "x2.5": "h1x25", "x3": "h1x3"}
CD4_TOK690I = {"x0.75": "cd4x075", "x1": "cd4x1", "x1.25": "cd4x125", "x1.5": "cd4x15", "x1.75": "cd4x175", "x2": "cd4x2"}


def _dump_name530(stack: str, strength: str, h1: str = "x2") -> str:
    """The stack-dump filename def_day.sh writes. Unset stack, strength and h1 -> today's name (690i: any 438 grid dose)."""
    base = f"dump_{H1_TOK690I[h1]}_{CD4_TOK690I[stack]}"
    if strength not in ("-", ""):
        base += "_s" + strength.replace("-", "m")
    return base + ".npz"


def _callset538(raw: str):
    """538: DEF_DAY_CALLSET_BLOCK grammar -> the level-meta fields riders_match requires."""
    if raw == "prop":
        return {"source": "prop", "swap": 1.0, "steer": "kept"}
    p = raw.split(":")
    if not p or p[0] != "prop" or len(p) not in (2, 3):
        return None
    try:
        swap = float(p[1])
    except ValueError:
        return None
    steer = p[2].lower() if len(p) == 3 else "kept"
    if not (0.0 < swap <= 1.0) or steer not in {"kept", "key", "cpm"}:
        return None
    return {"source": "prop", "swap": swap, "steer": steer}


def _october530_spec(spec: str):
    required = {"rc", "stack", "rho", "strength"}
    optional = {"lodo_hi", "callset", "h1"}   # 690i: h1
    out = {}
    for part in spec.split(","):
        k, sep, v = part.partition("=")
        if sep != "=" or k not in required | optional or v == "":
            return None
        out[k] = v
    if not required <= set(out):
        return None
    return out


def _october530_mismatch(cfg: dict, spec: dict) -> str | None:
    rc, stack, rho, strength = spec["rc"], spec["stack"], spec["rho"], spec["strength"]
    if stack not in CD4_TOK690I:
        return f"stack {stack!r} is not one of {sorted(CD4_TOK690I)}"
    h1 = spec.get("h1", "x2")   # 690i
    if h1 not in H1_TOK690I:
        return f"h1 {h1!r} is not one of {sorted(H1_TOK690I)}"
    want_name = _dump_name530(stack, strength, h1)
    if cfg.get("bulk_load") != want_name:
        return f"bulk_load {cfg.get('bulk_load')!r} vs {want_name}"
    g = cfg.get("cd4_target_g")
    if g is None or abs(float(g) - X1_G) >= 1e-12:
        return f"cd4_target_g {g!r} vs {X1_G} (the level build keeps DEF_ENV even when the dump is x1.5)"
    got_cpm = cfg.get("rc_callset_cpm")
    want_cpm = 2200 if rc == "-" else int(rc) if rc.isdigit() else None
    if want_cpm is None or got_cpm is None or int(got_cpm) != want_cpm:
        return f"rc_callset_cpm {got_cpm!r} vs {rc}"
    got_rho = cfg.get("rho_scale")
    if rho == "-":
        if got_rho is not None:
            return f"rho_scale present but not requested: {got_rho}"
    else:
        s = float(rho)
        if not 0.9 <= s <= 1.1:
            return f"rho scale {rho} is outside [0.9, 1.1]"
        # The builder does not record the loop target (530 does not touch it). A value that
        # IS recorded must be the one that was requested.
        if got_rho is not None and abs(float(got_rho) - s) >= 1e-12:
            return f"rho_scale meta {got_rho} vs {rho}"
    # 538: the level meta does not record the LODO cap (the pinned key's json does). A value that
    # IS recorded must be the one requested. Absence passes. callset IS a level-meta row.
    if "lodo_hi" in spec:
        got = cfg.get("lodo_hi")
        if got is not None and str(got) != spec["lodo_hi"]:
            try:
                same = abs(float(got) - float(spec["lodo_hi"])) < 1e-12
            except (TypeError, ValueError):
                same = False
            if not same:
                return f"lodo_hi meta {got!r} vs {spec['lodo_hi']}"
    if "callset" in spec:
        want = _callset538(spec["callset"])
        got = cfg.get("rc_callset_block")
        if want is None:
            return f"callset {spec['callset']!r} is not prop[:swap[:steer]]"
        if not isinstance(got, dict) or got.get("source") != "prop" or got.get("steer") != want["steer"]:
            return f"rc_callset_block {got!r} vs {spec['callset']}"
        sha = got.get("sha256")
        if not isinstance(sha, str) or len(sha) != 64:
            return f"rc_callset_block sha256 {sha!r}"
        try:
            swap_ok = abs(float(got.get("swap", -1.0)) - want["swap"]) < 1e-12
        except (TypeError, ValueError):
            swap_ok = False
        if not swap_ok:
            return f"rc_callset_block swap {got.get('swap')!r} vs {want['swap']}"
    return None


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
            scale = float(args[1]) if len(args) > 1 else 1.0
            star = rho_star(scale=scale)
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
        scale = 1.0
        specs = list(args)
        if specs and specs[0] == "--scale":
            scale = float(specs[1])
            specs = specs[2:]
        star = rho_star(scale=scale)
        best = None
        for spec in specs:
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
    if cmd == "riders_match":
        if len(args) not in (3, 4):
            print(f"riders_match: expected 3 or 4 arguments, got {len(args)}", file=sys.stderr)
            return 2
        meta_path, mag, sign = args[0], args[1], args[2]
        spec = args[3] if len(args) == 4 else None
        try:
            cfg = json.load(open(meta_path))["config"]
        except Exception as exc:
            print(f"riders_match failed: {exc}", file=sys.stderr)
            return 2
        got_mag = (cfg.get("de_mag") or {}).get("variant") if isinstance(cfg.get("de_mag"), dict) else None
        if (mag == "-") != (got_mag is None) or (mag != "-" and got_mag != mag):
            print(f"de_mag mismatch: meta {got_mag!r} vs {mag!r}", file=sys.stderr)
            return 1
        ms = cfg.get("de_mde_sign")
        if sign == "-":
            if ms is not None:
                print(f"de_mde_sign present but not requested: {ms}", file=sys.stderr)
                return 1
        else:
            want = _mde_sign_request(sign)
            bad = (not isinstance(ms, dict) or want is None
                   or any(ms.get(k) != v if isinstance(v, str) else (ms.get(k) is None or abs(float(ms[k]) - v) >= 1e-12)
                          for k, v in want.items()))
            if bad:
                print(f"de_mde_sign mismatch: meta {ms} vs {sign}", file=sys.stderr)
                return 1
        if spec is None:
            return 0
        parsed = _october530_spec(spec)
        if parsed is None:
            print(f"riders_match: bad October spec {spec!r}", file=sys.stderr)
            return 2
        msg = _october530_mismatch(cfg, parsed)
        if msg:
            print(msg, file=sys.stderr)
            return 1
        return 0
    if cmd == "meta_matches":
        meta_path, fact, qarg = args[0], args[1], args[2]
        try:
            cfg = json.load(open(meta_path))["config"]
            want = factors(fact)
            got_q = cfg.get("bulk_sharpen_q")
            if qarg.split(":", 1)[0] in ("tercile", "stercile", "tercilec"):   # 598: tercilec, the carried-only key tercile
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
