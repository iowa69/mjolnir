"""Database management: what each source is, how it is obtained, what landed.

:mod:`~mjolnir.db.registry` declares; :mod:`~mjolnir.db.fetch` acts. Consumers
elsewhere in the package want three things from here and should take them
through this namespace rather than reaching into either module: the path to an
installed file (:func:`database_file`, which fails with the command that would
fetch it), the versions and checksums to print in the report
(:func:`installed`), and the licence facts (:func:`attributions`).
"""

from .registry import (
    DATABASES,
    DB_GROUPS,
    Database,
    DatabaseFile,
    Licence,
    attributions,
    catalogue_database,
    check_redistribution,
    fetch_hint,
    latest_in_family,
    must_fetch,
    resolve_names,
    spec_for,
)
from .fetch import (
    database_file,
    database_version,
    fetch_database,
    fetch_databases,
    format_listing,
    installed,
    is_installed,
    list_databases,
    missing,
    verify_installed,
)

__all__ = [
    "DATABASES",
    "DB_GROUPS",
    "Database",
    "DatabaseFile",
    "Licence",
    "attributions",
    "catalogue_database",
    "check_redistribution",
    "database_file",
    "database_version",
    "fetch_database",
    "fetch_databases",
    "fetch_hint",
    "format_listing",
    "installed",
    "is_installed",
    "latest_in_family",
    "list_databases",
    "missing",
    "must_fetch",
    "resolve_names",
    "spec_for",
    "verify_installed",
]
