"""``python -m mjolnir`` entry point.

The console script and this module resolve to the same :func:`mjolnir.cli.main`,
so a run started either way produces byte-identical output — which matters
because the command line is recorded in the report.
"""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
