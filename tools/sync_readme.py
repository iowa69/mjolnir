#!/usr/bin/env python3
"""Write the README validation block from metrics.json, or check it still matches.

Every number this project publishes about its own accuracy is generated from a
measured file. Nothing is typed into the README by hand, because a hand-typed
figure survives the run being redone and then disagrees with the file it came
from — and the reader has no way to tell which one is current.

  sync_readme.py METRICS_JSON README.md          rewrite the block
  sync_readme.py METRICS_JSON README.md --check  exit 1 if it would change

--check is what CI and the test suite run, so a hand-edited number fails the
build rather than reaching a clinician.

The renderer deliberately refuses to print a bare accuracy figure. Every row
carries what the measurement establishes — agreement with another tool is not
correctness, and a row measured against a tool rather than against a phenotype
says so in the table, not in a footnote nobody reads.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

BEGIN = "<!-- METRICS:BEGIN -->"
END = "<!-- METRICS:END -->"

#: What a row was measured against. The distinction is the point: only
#: ``phenotype`` supports the word "accuracy".
ESTABLISHES = {
    "phenotype": "vs measured phenotype",
    "truth": "vs recorded truth",
    "agreement": "agreement only, not correctness",
    "structure": "vs known collection structure",
    "self": "internal consistency",
}


def _row(entry: dict) -> str:
    kind = entry.get("establishes", "agreement")
    if kind not in ESTABLISHES:
        raise SystemExit(
            "metrics.json: unknown 'establishes' value {0!r}; expected one of {1}".format(
                kind, ", ".join(sorted(ESTABLISHES))))
    return "| {0} | **{1}** | {2} | {3} |".format(
        entry["metric"], entry["mjolnir"], entry.get("comparator", "—"),
        ESTABLISHES[kind])


def render(metrics: dict) -> str:
    lines = [BEGIN, "", "## Validation", ""]
    lines.append(metrics["blurb"])
    lines.append("")
    lines.append("| | Mjolnir | Comparator | What this establishes |")
    lines.append("|---|---|---|---|")
    for entry in metrics["rows"]:
        lines.append(_row(entry))
    lines.append("")

    untested = metrics.get("untested") or []
    if untested:
        lines.append("**Not tested, and therefore not claimed:**")
        lines.append("")
        for item in untested:
            lines.append("- {0}".format(item))
        lines.append("")

    for figure in metrics.get("figures") or []:
        lines.append('<p align="center"><img src="{0}" alt="{1}" width="820"></p>'.format(
            figure["src"], figure["alt"]))
        lines.append("")

    lines.append("Generated from `{0}` by `tools/sync_readme.py`. ".format(
        metrics.get("source", "analysis/metrics.json"))
        + "The test suite re-runs it with `--check`, so a hand-edited figure fails the build.")
    lines.append("")
    lines.append(END)
    return "\n".join(lines)


def splice(readme: str, block: str) -> str:
    start = readme.index(BEGIN)
    stop = readme.index(END) + len(END)
    return readme[:start] + block + readme[stop:]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("metrics", type=Path)
    parser.add_argument("readme", type=Path)
    parser.add_argument("--check", action="store_true",
                        help="report whether the block is current, change nothing")
    args = parser.parse_args()

    if not args.metrics.exists():
        # Before the first validation run there are no measurements, and that is
        # a legitimate state: --check must not fail a build for an absence the
        # README already describes honestly.
        print("{0} does not exist yet; nothing to sync".format(args.metrics))
        return 0

    metrics = json.loads(args.metrics.read_text())
    readme = args.readme.read_text()
    if BEGIN not in readme or END not in readme:
        print("error: {0} has no {1} / {2} markers".format(args.readme, BEGIN, END),
              file=sys.stderr)
        return 2

    updated = splice(readme, render(metrics))
    if args.check:
        if updated == readme:
            print("README validation block matches metrics.json")
            return 0
        print("README validation block is stale; run sync_readme.py without --check",
              file=sys.stderr)
        return 1
    if updated != readme:
        args.readme.write_text(updated)
        print("updated {0}".format(args.readme))
    else:
        print("{0} already current".format(args.readme))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
