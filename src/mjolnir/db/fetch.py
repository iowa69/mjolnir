"""Obtain the databases, prove they are what the registry says, write down what landed.

``mjolnir db list`` is the part most people will use, and it touches the network
not at all: it prints every database with its licence, its citation and its
source, because someone deciding whether they may use these sources at all
should not have to download tens of megabytes to find out. The size of the
default set is computed from the registry rather than written down, so it cannot
go stale.

The fetching is deliberately unclever. One file at a time, no mirrors, no
partial-content resume, no fallback to a different source when one fails. A
resistance report is only comparable between two installations if both know
exactly which bytes produced it (design §5.5), and every "clever" recovery path
is a way of ending up with a catalogue nobody can name.

Verification is honest about what it can prove. Where the registry pins a
sha256, a mismatch is fatal. Where it pins only a git blob id — which is what
GitHub's API can give without transferring the file — the blob id is recomputed
from the received bytes and a mismatch means upstream changed the file since the
registry snapshot: loud, recorded in the manifest, and carried into the report's
database table, but fatal only under ``strict``, because upstream moving on is a
normal event and refusing to run would be the wrong answer to it. Where nothing
byte-stable exists at all — NCBI's efetch rewrites its own FASTA headers — the
sequence length and the accession are checked instead, and those *are* fatal,
since the wrong genome is never something to warn about and continue.

Whatever is measured goes into ``mjolnir_databases.json`` at the database root,
and every run reads it back as :class:`~mjolnir.records.DatabaseVersion` records
so the report can print the version and checksum of everything it consulted.
"""

from __future__ import annotations

import hashlib
import json
import os
import tarfile
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

from ..config import DB_ENV_VAR, default_db_dir
from ..records import DatabaseVersion
from ..utils import (
    LOG,
    MjolnirError,
    PathLike,
    ensure_dir,
    free_bytes,
    human_bytes,
    sha256sum,
    smart_open,
)
from .registry import (
    DATABASES,
    MANIFEST_NAME,
    NO_DATABASES_TEXT,
    REDISTRIBUTION_POLICY,
    SNAPSHOT_DATE,
    VERIFY_CHECKSUM,
    VERIFY_FASTA_LENGTH,
    VERIFY_NON_EMPTY,
    Database,
    DatabaseFile,
    attributions,
    fetch_hint,
    resolve_names,
    spec_for,
)

#: Upstream is told who is asking, so an operator reading their own logs — or
#: GitHub's rate limiter — can see what the traffic is.
USER_AGENT = "mjolnir (+https://github.com/iowa69/mjolnir)"

#: Seconds. The largest file here is 31 MB; anything slower than this timeout
#: allows is a broken connection rather than a slow one.
TIMEOUT = 300

#: Transport failures are retried, HTTP 4xx is not. A 404 means the registry is
#: wrong about where the file lives, and asking again cannot fix that.
RETRY_ATTEMPTS = 3
RETRY_BACKOFF = 2.0

#: Headroom demanded on top of the declared download size before starting.
#: Downloads are written whole, and a full filesystem halfway through a 31 MB
#: xlsx leaves a truncated file that verification would then have to catch.
FREE_SPACE_MARGIN = 1.5

_ProgressFn = Optional[Callable[[str], None]]


# ---------------------------------------------------------------------------
# Transport
# ---------------------------------------------------------------------------

def download(url: str, dest: PathLike, progress: _ProgressFn = None,
             attempts: int = RETRY_ATTEMPTS) -> Path:
    """Fetch *url* to *dest*, retrying transport failures only.

    The file is written to a ``.part`` sibling and moved into place at the end,
    so an interrupted download can never be mistaken for a complete one by the
    presence check that decides whether a database is installed.
    """
    target = Path(dest)
    last: Optional[Exception] = None
    for attempt in range(1, attempts + 1):
        try:
            return _download_once(url, target, progress)
        except MjolnirError as exc:
            last = exc
            if "HTTP 4" in str(exc) and "HTTP 429" not in str(exc):
                raise
            if attempt < attempts:
                delay = RETRY_BACKOFF * attempt
                LOG.warning("download failed (%s); retrying in %.0fs", exc, delay)
                time.sleep(delay)
    raise last if last is not None else MjolnirError("could not download " + url)


def _download_once(url: str, dest: Path, progress: _ProgressFn = None) -> Path:
    ensure_dir(dest.parent)
    partial = dest.with_name(dest.name + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            total = int(response.headers.get("Content-Length") or 0)
            done = 0
            step = 0
            with open(str(partial), "wb") as handle:
                while True:
                    block = response.read(1 << 16)
                    if not block:
                        break
                    handle.write(block)
                    done += len(block)
                    if progress and total and done // (1 << 22) > step:
                        step = done // (1 << 22)
                        progress("{0}: {1} of {2}".format(
                            dest.name, human_bytes(done), human_bytes(total)))
    except urllib.error.HTTPError as exc:
        _discard(partial)
        raise MjolnirError("could not download {0}: HTTP {1} {2}".format(
            url, exc.code, exc.reason)) from exc
    except (urllib.error.URLError, OSError) as exc:
        _discard(partial)
        raise MjolnirError("could not download {0}: {1}".format(url, exc)) from exc
    if not partial.exists() or partial.stat().st_size == 0:
        _discard(partial)
        raise MjolnirError("{0} returned an empty file".format(url))
    os.replace(str(partial), str(dest))
    return dest


def _discard(partial: Path) -> None:
    """Remove a half-written download.

    A ``.part`` file left behind is harmless — nothing looks for one — but it is
    also a lie about disk usage, and the next attempt writes over it anyway.
    """
    try:
        if partial.exists():
            partial.unlink()
    except OSError as exc:
        LOG.debug("could not remove %s: %s", partial, exc)


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def git_blob_sha1(path: PathLike) -> str:
    """The git object id of a file's contents.

    This is the one integrity anchor GitHub's API hands over without handing
    over the file, so it is the only pin the registry could be written with
    offline. Reproducing it locally is just sha1 over ``blob <size>\\0`` and the
    bytes — the same value ``git hash-object`` prints.
    """
    target = Path(path)
    size = target.stat().st_size
    digest = hashlib.sha1()
    digest.update("blob {0}\0".format(size).encode("ascii"))
    with open(str(target), "rb") as handle:
        while True:
            block = handle.read(1 << 20)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def fasta_sequence_length(path: PathLike) -> Tuple[int, str]:
    """Total sequence length and first header of a FASTA file.

    Used where the bytes are not stable but the sequence is: NCBI's efetch
    rewrites description text and line wrapping between fetches, and a checksum
    over that would fail for reasons that have nothing to do with the genome.
    """
    length = 0
    header = ""
    with smart_open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if not header:
                    header = line[1:].strip()
                continue
            length += len(line.strip())
    return length, header


def verify_file(spec: DatabaseFile, path: PathLike, strict: bool = False) -> Dict[str, Any]:
    """Check a downloaded file against whatever the registry could pin.

    Returns the measurements — sha256 always, blob id when one was pinned — for
    the manifest, plus a ``note`` that is empty when everything matched. The
    note is not decoration: it is what the report prints beside the database
    version when the file on disk is not the file the registry describes.
    """
    target = Path(path)
    if not target.exists():
        raise MjolnirError("{0} was not written to {1}".format(spec.name, target))
    size = target.stat().st_size
    if size == 0:
        raise MjolnirError("{0} downloaded as an empty file".format(spec.name))

    measured: Dict[str, Any] = {
        "size": size,
        "sha256": sha256sum(target),
        "verify": spec.verify,
        "note": "",
    }

    if spec.verify == VERIFY_NON_EMPTY:
        measured["note"] = (
            "{0} has no published integrity anchor; only its presence and "
            "non-emptiness were checked".format(spec.name))
        return measured

    if spec.verify == VERIFY_FASTA_LENGTH:
        length, header = fasta_sequence_length(target)
        measured["sequence_length"] = length
        measured["header"] = header
        if spec.expect_sequence_length and length != spec.expect_sequence_length:
            raise MjolnirError(
                "{0} holds {1} bases but {2} were expected. This is the wrong "
                "genome or a truncated download, and every catalogue coordinate "
                "would be wrong against it.\n  source: {3}".format(
                    spec.name, length, spec.expect_sequence_length, spec.url))
        if spec.expect_header and spec.expect_header not in header:
            raise MjolnirError(
                "{0} does not name {1!r} in its header (got {2!r}); coordinates "
                "against another assembly are silently wrong".format(
                    spec.name, spec.expect_header, header[:100]))
        return measured

    if spec.verify != VERIFY_CHECKSUM:  # pragma: no cover - guarded in the registry
        raise MjolnirError(
            "{0} declares verification strategy {1!r}, which fetch.py does not "
            "implement".format(spec.name, spec.verify))

    # sha256 when the registry pins one, otherwise the git blob id.
    if spec.sha256:
        if measured["sha256"] != spec.sha256:
            raise MjolnirError(
                "{0} failed its checksum.\n  expected sha256 {1}\n  measured    "
                " {2}\n  source: {3}\nDelete it and fetch again; if it fails "
                "twice the pinned digest and upstream disagree and the registry "
                "needs updating.".format(
                    spec.name, spec.sha256, measured["sha256"], spec.url))
        return measured

    if spec.git_blob_sha1:
        blob = git_blob_sha1(target)
        measured["git_blob_sha1"] = blob
        if blob != spec.git_blob_sha1:
            message = (
                "{0} differs from the registry snapshot taken {1}: upstream git "
                "blob {2}, received {3} ({4} on disk, {5} expected). Upstream has "
                "changed the file. The data may be a legitimate update; the "
                "checksum recorded for this run will not match another "
                "installation's.".format(
                    spec.name, SNAPSHOT_DATE, spec.git_blob_sha1, blob,
                    human_bytes(size), human_bytes(spec.size_bytes or size)))
            if strict:
                raise MjolnirError(message + "\n  --strict was requested, so this is fatal.")
            LOG.warning("%s", message)
            measured["note"] = message
        return measured

    measured["note"] = (
        "no integrity pin exists for {0}; only its presence and non-emptiness "
        "were checked".format(spec.name))
    LOG.warning("%s", measured["note"])
    return measured


# ---------------------------------------------------------------------------
# Unpacking
# ---------------------------------------------------------------------------

def unpack_archive(archive: PathLike, into: PathLike) -> Path:
    """Expand a verified tar or zip archive, refusing members that escape *into*.

    Verification happens before this, never after: an archive is checked as the
    bytes that arrived, because a member rewritten during extraction is no
    longer the thing the checksum described.
    """
    source = Path(archive)
    destination = ensure_dir(into)
    root = destination.resolve()

    def _safe(name: str) -> None:
        resolved = (destination / name).resolve()
        if root != resolved and root not in resolved.parents:
            raise MjolnirError(
                "{0} contains a member that would be written outside the "
                "database directory: {1}".format(source.name, name))

    if tarfile.is_tarfile(str(source)):
        with tarfile.open(str(source)) as tar:
            for member in tar.getmembers():
                if member.issym() or member.islnk():
                    raise MjolnirError(
                        "{0} contains a link member ({1}); Mjolnir does not "
                        "extract links from downloaded archives".format(
                            source.name, member.name))
                _safe(member.name)
            try:
                # 3.12 wants the filter named; 3.14 makes it the default. Ask for
                # it where it exists so the behaviour does not change under us.
                tar.extractall(str(destination), filter="data")
            except TypeError:
                tar.extractall(str(destination))
        return destination

    if zipfile.is_zipfile(str(source)):
        with zipfile.ZipFile(str(source)) as archive_file:
            for name in archive_file.namelist():
                _safe(name)
            archive_file.extractall(str(destination))
        return destination

    raise MjolnirError(
        "{0} is marked for unpacking but is neither a tar nor a zip "
        "archive".format(source))


# ---------------------------------------------------------------------------
# The manifest
# ---------------------------------------------------------------------------

def manifest_path(db_root: PathLike) -> Path:
    return Path(db_root).expanduser() / MANIFEST_NAME


def read_manifest(db_root: PathLike) -> Dict[str, Any]:
    """The record of what was fetched into *db_root*, or an empty record.

    A missing manifest is an empty result, not an error — a database root that
    has never been populated is an ordinary state. A manifest that exists and
    cannot be parsed *is* an error: it means something wrote over it, and
    guessing at what used to be installed there would put unverifiable versions
    into a clinical report.
    """
    path = manifest_path(db_root)
    if not path.exists():
        return {"databases": {}}
    try:
        with open(str(path), "rt") as handle:
            data = json.load(handle)
    except (OSError, ValueError) as exc:
        raise MjolnirError(
            "the database manifest {0} could not be read: {1}\nMove it aside and "
            "run 'mjolnir db fetch' again to rebuild it.".format(path, exc))
    if not isinstance(data, dict) or not isinstance(data.get("databases"), dict):
        raise MjolnirError(
            "the database manifest {0} is not in the expected form; move it "
            "aside and run 'mjolnir db fetch' again".format(path))
    return data


def write_manifest(db_root: PathLike, records: Dict[str, Any]) -> Path:
    """Merge *records* into the manifest at *db_root* and write it back."""
    root = ensure_dir(db_root)
    data = read_manifest(root)
    data["databases"].update(records)
    data["written"] = _now()
    data["db_root"] = str(Path(root).resolve())
    data["registry_snapshot"] = SNAPSHOT_DATE
    path = manifest_path(root)
    temporary = path.with_name(path.name + ".part")
    with open(str(temporary), "wt") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(str(temporary), str(path))
    return path


def _now() -> str:
    """UTC, because a database root is shared between machines in two timezones."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def database_checksum(files: Dict[str, Any]) -> str:
    """One digest standing for a whole database.

    A sha256 over ``name sha256`` lines, sorted, so that two installations can
    be compared with a single string in the report and a difference points at
    the database rather than at a file nobody printed.
    """
    digest = hashlib.sha256()
    for name in sorted(files):
        digest.update("{0} {1}\n".format(name, files[name].get("sha256", "")).encode("utf-8"))
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------

def plan(names: Optional[Sequence[str]] = None, db_root: Optional[PathLike] = None,
         force: bool = False) -> List[Dict[str, Any]]:
    """What :func:`fetch_databases` would do, without doing any of it.

    One entry per database, with the files that would be downloaded and the
    bytes that implies. The CLI prints this before starting and the free-space
    check is computed from it.
    """
    root = Path(db_root) if db_root is not None else default_db_dir()
    out: List[Dict[str, Any]] = []
    for name in resolve_names(names):
        spec = spec_for(name)
        if not spec.fetchable:
            continue
        wanted = [item for item in spec.files
                  if force or not _file_present(spec, item, root)]
        out.append({
            "name": spec.name,
            "version": spec.version,
            "directory": str(spec.directory(root)),
            "files": [item.name for item in wanted],
            "bytes": sum(item.size_bytes for item in wanted),
            "licence": spec.licence.spdx,
            "already_present": not wanted,
        })
    return out


def _file_present(spec: Database, item: DatabaseFile, db_root: PathLike) -> bool:
    path = spec.directory(db_root) / item.name
    return path.exists() and path.stat().st_size > 0


def _check_free_space(db_root: PathLike, needed: int) -> None:
    """Refuse before downloading rather than halfway through."""
    if needed <= 0:
        return
    probe = Path(db_root).expanduser()
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    available = free_bytes(probe)
    if available and available < needed * FREE_SPACE_MARGIN:
        raise MjolnirError(
            "not enough free space under {0}: {1} available, {2} needed for the "
            "download plus working room. Free some space or point ${3} at "
            "another filesystem.".format(
                probe, human_bytes(available), human_bytes(int(needed * FREE_SPACE_MARGIN)),
                DB_ENV_VAR))


def fetch_database(name: str, db_root: Optional[PathLike] = None, force: bool = False,
                   strict: bool = False, progress: _ProgressFn = None) -> DatabaseVersion:
    """Download, verify and record one database. Returns what the report prints."""
    spec = spec_for(name)
    root = Path(db_root) if db_root is not None else default_db_dir()
    if not spec.fetchable:
        raise MjolnirError(
            "{0} is not something Mjolnir downloads.\n  {1}\n  see: {2}".format(
                spec.name, spec.note or spec.what, spec.homepage))

    directory = ensure_dir(spec.directory(root))
    todo = set(item.name for item in spec.files
               if force or not _file_present(spec, item, root))
    _check_free_space(directory, sum(item.size_bytes for item in spec.files
                                     if item.name in todo))

    measured: Dict[str, Any] = {}
    notes: List[str] = []
    existing = read_manifest(root)["databases"].get(spec.name, {})
    previous_files = existing.get("files", {}) if isinstance(existing, dict) else {}

    for item in spec.files:
        target = directory / item.name
        if item.name in todo:
            LOG.info("fetching %s -> %s (%s)", item.name, target,
                     human_bytes(item.size_bytes) if item.size_bytes else "size unknown")
            try:
                download(item.url, target, progress)
            except MjolnirError:
                if item.required:
                    raise
                # An optional file that upstream no longer publishes is a gap in
                # the record, not a failed install; it is named in the manifest
                # so the report does not imply the database is complete.
                LOG.warning("optional file %s could not be fetched from %s",
                            item.name, item.url)
                notes.append("optional file {0} was not obtained".format(item.name))
                continue
            check = verify_file(item, target, strict=strict)
            if item.unpack:
                # Verified first, expanded second: a member rewritten during
                # extraction is no longer the bytes the checksum described.
                unpack_archive(target, directory)
                check["unpacked_into"] = str(directory)
        elif target.exists():
            check = dict(previous_files.get(item.name) or {})
            if not check:
                check = verify_file(item, target, strict=strict)
        elif item.required:
            raise MjolnirError(
                "{0} is missing from {1} and was not re-fetched; run with "
                "--force".format(item.name, directory))
        else:
            continue
        measured[item.name] = check
        if check.get("note"):
            notes.append(check["note"])

    version = _resolve_version(spec, directory)
    licence_note = _record_licence(spec, directory)
    if licence_note:
        notes.append(licence_note)
    if spec.must_fetch:
        notes.append(
            "fetched at install time rather than distributed with Mjolnir: {0}".format(
                spec.licence.describe()))

    record = {
        "name": spec.name,
        "version": version,
        "checksum": database_checksum(measured),
        "path": str(directory),
        "licence": spec.licence.spdx,
        "licence_name": spec.licence.name,
        "citation": spec.citation,
        "url": spec.homepage,
        "fetched": _now(),
        "attribution": spec.licence.attribution,
        "redistributable": spec.redistributable,
        "registry_snapshot": SNAPSHOT_DATE,
        "note": "; ".join(notes),
        "files": measured,
    }
    write_manifest(root, {spec.name: record})
    return _as_version(record)


def fetch_databases(names: Optional[Sequence[str]] = None,
                    db_root: Optional[PathLike] = None, force: bool = False,
                    strict: bool = False,
                    progress: _ProgressFn = None) -> List[DatabaseVersion]:
    """Fetch a set of databases, or the default set when none is named.

    Named ``fetch_databases`` rather than ``fetch`` so that ``from mjolnir.db
    import fetch`` keeps meaning the module: a package that exports a function
    with the same name as one of its submodules shadows it on import, and the
    resulting AttributeError arrives a long way from its cause.

    Failures are not swallowed: the first database that cannot be installed
    stops the run. A partially-populated database root that reports success is
    the failure mode this whole project is written against.
    """
    wanted = [name for name in resolve_names(names) if spec_for(name).fetchable]
    if not wanted:
        raise MjolnirError(
            "nothing to fetch: every name given is a database Mjolnir does not "
            "download. Run 'mjolnir db list' to see what each one needs.")
    total = sum(entry["bytes"] for entry in plan(wanted, db_root, force))
    LOG.info("fetching %d database(s), about %s to download",
             len(wanted), human_bytes(total))
    return [fetch_database(name, db_root, force=force, strict=strict, progress=progress)
            for name in wanted]


def _resolve_version(spec: Database, directory: Path) -> str:
    """The version to record, read from the data where the data carries one.

    tbdb publishes ``db-schema-version`` in ``variables.json``; that is a fact
    about the files on disk, where the registry's version string is only a fact
    about the day the registry was written. Where a database states its own
    version, the stated one wins and the registry's is kept beside it — the
    report needs both, since neither alone identifies the files.
    """
    if not spec.version_file:
        return spec.version
    stated = directory / spec.version_file
    if not stated.exists():
        return spec.version
    try:
        with open(str(stated), "rt") as handle:
            declared = json.load(handle).get(spec.version_key)
    except (OSError, ValueError) as exc:
        raise MjolnirError(
            "{0} states its version in {1}, which could not be read: {2}".format(
                spec.name, stated, exc))
    if not declared:
        raise MjolnirError(
            "{0} states its version in {1} under {2!r}, which is absent or "
            "empty; the file upstream publishes has changed shape".format(
                spec.name, stated, spec.version_key))
    return "{0} {1} ({2})".format(spec.version_key, declared, spec.version)


#: Written into every database directory, so a directory copied to another
#: machine carries its own terms rather than depending on this source tree.
LICENCE_STAMP = "mjolnir_licence.txt"


def _record_licence(spec: Database, directory: Path) -> str:
    """Write the licence terms beside the data and check upstream's own copy.

    ODC-By is an attribution licence: the WHO data may be redistributed, and may
    not be redistributed anonymously. A directory holding the catalogue without
    the terms that came with it is the thing that turns a permitted mirror into
    a breach, so the terms are written next to it every time.
    """
    lines = [
        "{0} {1}".format(spec.name, spec.version),
        spec.title,
        "",
        "Licence:    {0}".format(spec.licence.describe()),
        "Terms:      {0}".format(spec.licence.url or "see the provider"),
        "Provider:   {0}".format(spec.provider),
        "Source:     {0}".format(spec.homepage),
        "Citation:   {0}".format(spec.citation),
    ]
    if spec.licence.attribution:
        lines += ["", "Attribution required with any redistribution or publication:",
                  "  " + spec.licence.attribution]
    if spec.licence.note:
        lines += ["", "Note: " + spec.licence.note]
    lines += ["", "Fetched by Mjolnir on {0}. {1}".format(_now(), REDISTRIBUTION_POLICY)]
    (directory / LICENCE_STAMP).write_text("\n".join(lines) + "\n")

    if not spec.licence_file:
        return ""
    upstream = directory / spec.licence_file
    if not upstream.exists():
        raise MjolnirError(
            "{0} declares its licence must be verified on this machine but "
            "{1} is not present at {2}".format(spec.name, spec.licence_file, upstream))
    first = ""
    with smart_open(upstream, "rt") as handle:
        for line in handle:
            if line.strip():
                first = line.strip()
                break
    LOG.info("%s licence file %s: %s", spec.name, spec.licence_file, first[:120])
    return "upstream licence file {0} reads: {1}".format(spec.licence_file, first[:200])


# ---------------------------------------------------------------------------
# Reading back what is installed
# ---------------------------------------------------------------------------

def _as_version(record: Dict[str, Any]) -> DatabaseVersion:
    return DatabaseVersion(
        name=str(record.get("name", "")),
        version=str(record.get("version") or "unknown"),
        checksum=str(record.get("checksum", "")),
        path=str(record.get("path", "")),
        licence=str(record.get("licence", "")),
        citation=str(record.get("citation", "")),
        url=str(record.get("url", "")),
        fetched=str(record.get("fetched", "")),
        note=str(record.get("note", "")),
    )


def installed(db_root: Optional[PathLike] = None) -> List[DatabaseVersion]:
    """Every database recorded as installed under *db_root*.

    This is what the report's database table is built from, so it reports what
    the manifest says was fetched — not what the registry says exists.
    """
    root = Path(db_root) if db_root is not None else default_db_dir()
    entries = read_manifest(root)["databases"]
    return [_as_version(entries[name]) for name in sorted(entries)]


def database_version(name: str, db_root: Optional[PathLike] = None) -> DatabaseVersion:
    """The installed record for one database, or a MjolnirError saying how to get it."""
    spec = spec_for(name)
    root = Path(db_root) if db_root is not None else default_db_dir()
    record = read_manifest(root)["databases"].get(spec.name)
    if record is None:
        raise MjolnirError(
            "{0} is not installed under {1}.\n  fetch it with: {2}".format(
                spec.name, root, fetch_hint(spec.name)))
    return _as_version(record)


def is_installed(name: str, db_root: Optional[PathLike] = None) -> bool:
    """Whether every required file of *name* is present and non-empty."""
    spec = spec_for(name)
    root = Path(db_root) if db_root is not None else default_db_dir()
    required = spec.required_files()
    if not required:
        return False
    return all(_file_present(spec, item, root) for item in required)


def missing(names: Optional[Sequence[str]] = None,
            db_root: Optional[PathLike] = None) -> List[str]:
    """The fetchable databases among *names* that are not installed."""
    root = Path(db_root) if db_root is not None else default_db_dir()
    return [name for name in resolve_names(names)
            if spec_for(name).fetchable and not is_installed(name, root)]


def database_file(name: str, filename: str,
                  db_root: Optional[PathLike] = None) -> Path:
    """The path to one file of an installed database.

    Every consumer — the catalogue loaders, the barcode reader, the masker —
    asks for its input through here, so that "the file is not there" is answered
    once, with the command that would fetch it, instead of once per module with
    a different message each time.
    """
    spec = spec_for(name)
    root = Path(db_root) if db_root is not None else default_db_dir()
    item = spec.file(filename)
    path = spec.directory(root) / item.name
    if not path.exists():
        raise MjolnirError(
            "{0} ({1}) is not present at {2}.\n  fetch it with: {3}\n  it is "
            "needed for: {4}".format(
                item.name, spec.title, path, fetch_hint(spec.name),
                ", ".join(spec.required_for) or spec.what))
    if path.stat().st_size == 0:
        raise MjolnirError(
            "{0} is present but empty at {1}; re-fetch with: {2} --force".format(
                item.name, path, fetch_hint(spec.name)))
    return path


def verify_installed(name: str, db_root: Optional[PathLike] = None,
                     deep: bool = False) -> List[str]:
    """Problems with an installed database: missing files, or changed contents.

    ``deep`` re-reads every file and compares its sha256 against the manifest,
    which is how a database edited in place after installation is caught. It is
    off by default because it re-hashes every file in the database, and
    ``mjolnir doctor`` is expected to be cheap.
    """
    spec = spec_for(name)
    root = Path(db_root) if db_root is not None else default_db_dir()
    record = read_manifest(root)["databases"].get(spec.name)
    problems: List[str] = []
    if record is None:
        return ["{0} is not recorded as installed; {1}".format(
            spec.name, fetch_hint(spec.name))]
    if record.get("note"):
        problems.append("{0}: {1}".format(spec.name, record["note"]))
    recorded = record.get("files", {})
    for item in spec.files:
        path = spec.directory(root) / item.name
        if not path.exists():
            if item.required:
                problems.append("{0}: {1} is missing from {2}".format(
                    spec.name, item.name, path.parent))
            continue
        if not deep:
            continue
        expected = (recorded.get(item.name) or {}).get("sha256")
        if not expected:
            problems.append("{0}: {1} has no recorded checksum to compare "
                            "against".format(spec.name, item.name))
            continue
        if sha256sum(path) != expected:
            problems.append(
                "{0}: {1} has changed since it was fetched (sha256 differs from "
                "the manifest); re-fetch with {2} --force".format(
                    spec.name, item.name, fetch_hint(spec.name)))
    return problems


# ---------------------------------------------------------------------------
# Listing - the mode that downloads nothing
# ---------------------------------------------------------------------------

def list_databases(db_root: Optional[PathLike] = None) -> List[Dict[str, Any]]:
    """Every registered database with its licence, citation and install state.

    No network access, by construction: this is the mode an operator runs before
    deciding whether their institution may use these sources at all.
    """
    root = Path(db_root) if db_root is not None else default_db_dir()
    entries = read_manifest(root)["databases"]
    out: List[Dict[str, Any]] = []
    for name in sorted(DATABASES):
        spec = DATABASES[name]
        record = entries.get(name, {})
        out.append({
            "name": spec.name,
            "title": spec.title,
            "what": spec.what,
            "family": spec.family,
            "version": spec.version,
            "installed_version": record.get("version", ""),
            "version_date": spec.version_date,
            "provider": spec.provider,
            "homepage": spec.homepage,
            "citation": spec.citation,
            "licence": spec.licence.spdx,
            "licence_name": spec.licence.name,
            "licence_url": spec.licence.url,
            "redistributable": spec.redistributable,
            "distribution": ("may be redistributed with attribution"
                             if spec.redistributable else
                             "fetched at install time; not redistributable"),
            "attribution": spec.licence.attribution,
            "fetchable": spec.fetchable,
            "auto": spec.auto,
            "bytes": spec.total_bytes,
            "files": [item.name for item in spec.files],
            "required_for": list(spec.required_for),
            "superseded_by": spec.superseded_by,
            "successor_watch": spec.successor_watch,
            "installed": is_installed(name, root) if spec.fetchable else False,
            "path": str(spec.directory(root)),
            "checksum": record.get("checksum", ""),
            "fetched": record.get("fetched", ""),
            "note": spec.note,
        })
    return out


def format_listing(db_root: Optional[PathLike] = None, verbose: bool = False) -> str:
    """The human-readable form of :func:`list_databases`, for ``mjolnir db list``."""
    root = Path(db_root) if db_root is not None else default_db_dir()
    rows = list_databases(root)
    lines: List[str] = [
        "Mjolnir databases",
        "  database root:     {0}".format(root),
        "  registry snapshot: {0}".format(SNAPSHOT_DATE),
        "  " + REDISTRIBUTION_POLICY,
        "",
    ]
    if not any(row["installed"] for row in rows):
        lines += ["  " + NO_DATABASES_TEXT, ""]

    for row in rows:
        if row["fetchable"]:
            state = "installed" if row["installed"] else "not installed"
        else:
            state = "not fetched by Mjolnir"
        lines.append("{0}  [{1}]".format(row["name"], state))
        lines.append("  {0}".format(row["title"]))
        lines.append("  version:      {0}{1}".format(
            row["installed_version"] or row["version"],
            " (registry: {0})".format(row["version"]) if row["installed_version"]
            and row["installed_version"] != row["version"] else ""))
        lines.append("  licence:      {0} ({1})".format(row["licence"], row["licence_name"]))
        lines.append("  distribution: {0}".format(row["distribution"]))
        if row["licence_url"]:
            lines.append("  terms:        {0}".format(row["licence_url"]))
        lines.append("  source:       {0}".format(row["homepage"]))
        lines.append("  citation:     {0}".format(row["citation"]))
        if row["bytes"]:
            lines.append("  download:     {0} in {1} file(s)".format(
                human_bytes(row["bytes"]), len(row["files"])))
        if row["required_for"]:
            lines.append("  needed for:   {0}".format(", ".join(row["required_for"])))
        if row["checksum"]:
            lines.append("  checksum:     {0}".format(row["checksum"]))
            lines.append("  fetched:      {0}".format(row["fetched"]))
        if verbose:
            lines.append("  what:         {0}".format(row["what"]))
            if row["files"]:
                lines.append("  files:        {0}".format(", ".join(row["files"])))
            if row["note"]:
                lines.append("  note:         {0}".format(row["note"]))
            if row["successor_watch"]:
                lines.append("  versioning:   {0}".format(row["successor_watch"]))
        if row["superseded_by"]:
            lines.append("  superseded by: {0}".format(row["superseded_by"]))
        lines.append("")

    owed = attributions()
    if owed:
        lines.append("Attribution owed by any run that uses these databases:")
        for text in owed:
            lines.append("  - {0}".format(text))
        lines.append("")
    return "\n".join(lines)


def format_plan(names: Optional[Sequence[str]] = None,
                db_root: Optional[PathLike] = None, force: bool = False) -> str:
    """What a fetch would download, for printing before it starts."""
    entries = plan(names, db_root, force)
    total = sum(entry["bytes"] for entry in entries)
    lines = ["would fetch {0} database(s), {1}:".format(len(entries), human_bytes(total))]
    for entry in entries:
        if entry["already_present"]:
            lines.append("  {0}: already present at {1}".format(
                entry["name"], entry["directory"]))
            continue
        lines.append("  {0} {1}: {2} file(s), {3} -> {4}".format(
            entry["name"], entry["version"], len(entry["files"]),
            human_bytes(entry["bytes"]), entry["directory"]))
    return "\n".join(lines)


def iter_files(names: Optional[Sequence[str]] = None) -> Iterable[Tuple[str, DatabaseFile]]:
    """(database name, file) for everything the named databases would fetch."""
    for name in resolve_names(names):
        spec = spec_for(name)
        for item in spec.files:
            yield spec.name, item
