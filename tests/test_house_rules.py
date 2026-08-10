"""The house rules, enforced against the source rather than against a reviewer.

Four of Mjolnir's rules are structural, and structural rules decay silently: the
first bare ``except`` is added in a hurry, the first ``from __future__ import
annotations`` is forgotten in a new module, and the failure shows up months later
as a run that produced a clean-looking result out of a swallowed exception.

So they are checked here by walking the package's own syntax trees. Nothing is
imported for the AST checks — a module that fails to import for want of an
optional dependency still has its source read — and nothing runs.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

PACKAGE = Path(__file__).resolve().parent.parent / "src" / "mjolnir"
MODULES = sorted(PACKAGE.rglob("*.py"))


def _tree(path):
    return ast.parse(path.read_text(), filename=str(path))


def test_the_package_was_found():
    assert MODULES, "no modules under src/mjolnir"


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_no_bare_except(path):
    """A bare ``except`` turns a failure into a clean-looking result."""
    for node in ast.walk(_tree(path)):
        if isinstance(node, ast.ExceptHandler) and node.type is None:
            pytest.fail("{0}:{1} has a bare except".format(path.name, node.lineno))


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_no_match_statement(path):
    """``match`` is 3.10+; the package supports 3.9."""
    match_node = getattr(ast, "Match", None)
    if match_node is None:  # pragma: no cover - only on very old interpreters
        pytest.skip("this interpreter cannot parse match statements at all")
    for node in ast.walk(_tree(path)):
        if isinstance(node, match_node):
            pytest.fail("{0}:{1} uses a match statement".format(path.name, node.lineno))


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_future_annotations_or_no_annotations_at_all(path):
    """PEP-604 unions and builtin generics are only safe behind the future import."""
    tree = _tree(path)
    has_future = any(
        isinstance(node, ast.ImportFrom) and node.module == "__future__"
        and any(alias.name == "annotations" for alias in node.names)
        for node in tree.body)
    annotated = any(
        isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.AnnAssign))
        for node in ast.walk(tree))
    if annotated and len(tree.body) > 1:
        assert has_future, \
            "{0} carries annotations without 'from __future__ import annotations'".format(
                path.name)


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_every_module_compiles_under_this_interpreter(path):
    compile(path.read_text(), str(path), "exec")


def test_no_module_raises_a_bare_exception():
    """Errors are explicit: MjolnirError, or a standard exception with a message."""
    offenders = []
    for path in MODULES:
        for node in ast.walk(_tree(path)):
            if isinstance(node, ast.Raise) and isinstance(node.exc, ast.Name):
                offenders.append("{0}:{1}".format(path.name, node.lineno))
    assert not offenders, "raise without an instance: {0}".format(", ".join(offenders))


def test_every_module_imports_without_any_optional_dependency():
    """``mjolnir doctor`` has to report a missing tool, not die on it.

    Imported through :func:`importlib.import_module` rather than as bare import
    statements. ``# noqa`` is a flake8 convention and pyflakes does not honour
    it, so a block of deliberately-unused imports fails the lint step that CI
    runs over this directory — the test for the house rules would have been the
    thing that broke the build.

    Walking the package also means a module added later is covered without
    anyone remembering to add it here.
    """
    import importlib
    import pkgutil

    import mjolnir

    failures = []
    for module in pkgutil.walk_packages(mjolnir.__path__, "mjolnir."):
        try:
            importlib.import_module(module.name)
        except Exception as exc:  # noqa: BLE001 - the failure is the finding
            failures.append("{0}: {1}: {2}".format(
                module.name, type(exc).__name__, exc))
    assert not failures, "modules that do not import: {0}".format("; ".join(failures))
