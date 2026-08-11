"""The conda recipe must keep up with the code.

The recipe's own comment explains why it lists every module rather than every
package: *a missing submodule imports fine through its parent and fails at the
first run that needs it*. That is a real failure mode and the recipe was written
against it — and then two new modules were added and neither was listed, so the
package test would have passed on a build that shipped a broken tool.

A comment cannot enforce itself. This does.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
RECIPE = ROOT / "conda-recipe" / "meta.yaml"
PACKAGE = ROOT / "src" / "mjolnir"


def shipped_modules():
    """Every importable module under src/mjolnir, in dotted form."""
    found = set()
    for path in PACKAGE.rglob("*.py"):
        parts = list(path.relative_to(ROOT / "src").with_suffix("").parts)
        if parts[-1] == "__init__":
            parts = parts[:-1]
        found.add(".".join(parts))
    return found


def listed_modules():
    """Every module named in the recipe's import test."""
    return set(re.findall(r"^\s+- (mjolnir[\w.]*)$", RECIPE.read_text(), re.M))


@pytest.mark.skipif(not RECIPE.exists(), reason="no conda recipe in this checkout")
def test_every_shipped_module_is_in_the_recipe_import_test():
    missing = sorted(shipped_modules() - listed_modules())
    assert not missing, (
        "these modules ship but are not imported by the conda package test, so a "
        "build that broke them would pass: {0}".format(", ".join(missing)))


@pytest.mark.skipif(not RECIPE.exists(), reason="no conda recipe in this checkout")
def test_the_recipe_does_not_name_modules_that_no_longer_exist():
    """A stale entry fails the build for a module nobody removed on purpose."""
    extra = sorted(listed_modules() - shipped_modules())
    assert not extra, (
        "the recipe imports modules that are not in the package: {0}".format(
            ", ".join(extra)))


@pytest.mark.skipif(not RECIPE.exists(), reason="no conda recipe in this checkout")
def test_the_package_data_is_declared_in_pyproject_as_well_as_the_manifest():
    """The playbooks reach the wheel through include-package-data defaulting true.

    Naming them in pyproject makes the guarantee explicit rather than incidental:
    a wheel that ships no playbooks installs cleanly and then runs with no
    organism knowledge, which is indistinguishable from an organism that has no
    playbook.
    """
    pyproject = (ROOT / "pyproject.toml").read_text()
    assert "agent/playbooks/*.yaml" in pyproject
