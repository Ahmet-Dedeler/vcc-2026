"""Static scan: repo modules imported inside functions (not at module level) by the R7a instrument.

The instrument pin (r7a551.instrument_files) imports lib551/designs551/price551/gate551 and exec's only the
source of full551/comp551/aggregate551 above the first top-level ``def``. Imports that live inside functions
are invisible to that probe. This script walks the AST of those modules and of every repo module they import
at module level (the closure), and reports Import / ImportFrom nodes whose parent is not the module body.
"""
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PINNED = {
    "experiments/062-cd4-post-whiten/run.py",
    "experiments/064-pm-cd4-composition/run.py",
    "experiments/162-retrodicting-validator/plolo.py",
    "experiments/162-retrodicting-validator/run_arms.py",
    "experiments/169-sharpening-six-member/sharp.py",
    "experiments/187-joint-amplitude/cliff187.py",
    "experiments/187-joint-amplitude/joint187.py",
    "experiments/187-joint-amplitude/mse187.py",
    "experiments/218-independent-donors/compose218.py",
    "experiments/218-independent-donors/resid218.py",
    "experiments/218-independent-donors/sources218.py",
    "experiments/228-theta-gain-plolo/arms228.py",
    "experiments/410-pertscale-revisit/key410.py",
    "experiments/410-pertscale-revisit/size410.py",
    "experiments/414-bulk-reprice/arms414.py",
    "experiments/414-bulk-reprice/ballast414.py",
    "experiments/415-lodo-key/lodo415.py",
    "experiments/415-lodo-key/sens415.py",
    "experiments/418-bulk-composition/analyse418.py",
    "experiments/418-bulk-composition/lib418.py",
    "experiments/420-bulk-stack/lib420.py",
    "experiments/551-r6-c10b/comp551.py",
    "experiments/551-r6-c10b/designs551.py",
    "experiments/551-r6-c10b/full551.py",
    "experiments/551-r6-c10b/gate551.py",
    "experiments/551-r6-c10b/guard551.py",
    "experiments/551-r6-c10b/lib551.py",
    "experiments/551-r6-c10b/price551.py",
}
SEEDS = [
    "experiments/551-r6-c10b/lib551.py",
    "experiments/551-r6-c10b/full551.py",
    "experiments/551-r6-c10b/comp551.py",
    "experiments/551-r6-c10b/price551.py",
    "experiments/551-r6-c10b/aggregate551.py",
    "experiments/551-r6-c10b/designs551.py",
    "experiments/551-r6-c10b/gate551.py",
]
# directories the instrument puts on sys.path (lib551, price551, analyse418)
PATH_DIRS = [
    ROOT / "experiments/420-bulk-stack",
    ROOT / "experiments/418-bulk-composition",
    ROOT / "experiments/417-bulk-sharpen",
    ROOT / "experiments/415-lodo-key",
    ROOT / "experiments/414-bulk-reprice",
    ROOT / "experiments/410-pertscale-revisit",
    ROOT / "experiments/228-theta-gain-plolo",
    ROOT / "experiments/218-independent-donors",
    ROOT / "experiments/187-joint-amplitude",
    ROOT / "experiments/169-sharpening-six-member",
    ROOT / "experiments/162-retrodicting-validator",
    ROOT / "experiments/096-mse-floor",
    ROOT / "experiments/064-pm-cd4-composition",
    ROOT / "experiments/062-cd4-post-whiten",
    ROOT / "experiments/551-r6-c10b",
    ROOT / "src",
    ROOT,
]


def resolve(name: str, level: int, src: Path) -> str | None:
    """Best-effort map of an imported name to a repo-relative .py path. None if not a repo file."""
    if level and level > 0:
        return None  # no relative imports in this stack
    parts = name.split(".")
    candidates = []
    for d in PATH_DIRS:
        p = d.joinpath(*parts).with_suffix(".py")
        candidates.append(p)
        p2 = d.joinpath(*parts, "__init__.py")
        candidates.append(p2)
    # also experiments/*/<name>.py by glob of the first component only when it is a bare module
    if len(parts) == 1:
        for hit in (ROOT / "experiments").glob(f"*/{parts[0]}.py"):
            candidates.append(hit)
        for hit in (ROOT / "src").glob(f"**/{parts[0]}.py"):
            candidates.append(hit)
    for c in candidates:
        if c.is_file():
            try:
                return str(c.relative_to(ROOT))
            except ValueError:
                return None
    return None


def module_level_imports(tree: ast.AST) -> list[tuple[str, int]]:
    out = []
    for node in tree.body if isinstance(tree, ast.Module) else []:
        if isinstance(node, ast.Import):
            for a in node.names:
                out.append((a.name, node.lineno))
        elif isinstance(node, ast.ImportFrom) and node.module:
            out.append((node.module, node.lineno))
    return out


def lazy_imports(tree: ast.AST) -> list[dict]:
    """Import nodes that are not direct children of the module body."""
    top = set(id(n) for n in (tree.body if isinstance(tree, ast.Module) else []))
    found = []

    class V(ast.NodeVisitor):
        def __init__(self):
            self.stack: list[str] = []

        def visit_FunctionDef(self, node):
            self.stack.append(node.name)
            self.generic_visit(node)
            self.stack.pop()

        visit_AsyncFunctionDef = visit_FunctionDef

        def visit_ClassDef(self, node):
            self.stack.append(node.name)
            self.generic_visit(node)
            self.stack.pop()

        def visit_Import(self, node):
            if id(node) not in top and self.stack:
                for a in node.names:
                    found.append(dict(name=a.name, line=node.lineno, where=".".join(self.stack), kind="import"))
            self.generic_visit(node)

        def visit_ImportFrom(self, node):
            if id(node) not in top and self.stack and node.module:
                found.append(dict(name=node.module, line=node.lineno, where=".".join(self.stack),
                                  names=[a.name for a in node.names], kind="from"))
            self.generic_visit(node)

    V().visit(tree)
    return found


def main() -> None:
    seen: set[str] = set()
    queue = list(SEEDS)
    closure: list[str] = []
    while queue:
        rel = queue.pop(0)
        if rel in seen:
            continue
        path = ROOT / rel
        if not path.is_file():
            continue
        seen.add(rel)
        closure.append(rel)
        tree = ast.parse(path.read_text(), filename=str(path))
        for name, _ln in module_level_imports(tree):
            hit = resolve(name, 0, path)
            if hit and hit not in seen and hit.startswith("experiments/"):
                queue.append(hit)
            elif hit and hit not in seen and (hit.startswith("src/") or hit.startswith("scripts/")):
                queue.append(hit)

    rows = []
    for rel in closure:
        path = ROOT / rel
        tree = ast.parse(path.read_text(), filename=str(path))
        for item in lazy_imports(tree):
            hit = resolve(item["name"], 0, path)
            item.update(file=rel, resolved=hit, pinned=bool(hit and hit in PINNED),
                        repo=bool(hit and (hit.startswith("experiments/") or hit.startswith("src/") or hit.startswith("scripts/"))))
            rows.append(item)

    uncovered = [r for r in rows if r["repo"] and not r["pinned"]]
    covered = [r for r in rows if r["repo"] and r["pinned"]]
    external = [r for r in rows if not r["repo"]]
    report = dict(closure_n=len(closure), closure=closure, n_lazy=len(rows),
                  uncovered_repo=uncovered, covered_repo=covered, non_repo=external)
    out = Path(__file__).resolve().parent / "lazy_imports.json"
    out.write_text(json.dumps(report, indent=1) + "\n")
    print(f"closure {len(closure)} files; lazy import sites {len(rows)}; "
          f"repo and NOT in the 28 pins: {len(uncovered)}")
    for r in uncovered:
        print(f"  UNCOVERED {r['file']}:{r['line']} in {r['where']} imports {r['name']} -> {r['resolved']}")
    print("repo lazy imports that ARE pinned:")
    for r in covered:
        print(f"  pinned    {r['file']}:{r['line']} in {r['where']} imports {r['name']}")
    print("non-repo lazy imports:")
    for r in external:
        print(f"  external  {r['file']}:{r['line']} in {r['where']} imports {r['name']}")


if __name__ == "__main__":
    sys.exit(main())
