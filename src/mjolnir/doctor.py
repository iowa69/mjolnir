"""The whole environment, reported before anything runs.

The failure this module exists to prevent is specific and has happened to every
one of the tools Mjolnir replaces: a run starts, spends forty minutes mapping,
and dies at the variant-calling step because Clair3 was never installed. Or
worse, it does not die — a missing database turns into an empty catalogue turns
into a report with no resistance determinants, which reads exactly like a
susceptible isolate.

So ``mjolnir doctor`` probes everything up front: every external tool, every
Python dependency, every database in the registry, and the derived question that
actually matters — which of the pipeline's capabilities can run at all. A
missing optional tool is reported as a lost capability with its consequence
named, never as a footnote.

Nothing in here raises. A probe that cannot be completed becomes a row saying so
with the reason attached; a doctor that aborted at the first missing tool would
report one problem when there are five, and the operator would install them one
run at a time.

The one thing it will not do is call anything "fine" that it did not measure.
Absent tools are absent, an unresolvable database is unresolved, and a Kraken2
index that is not a mycobacterial pangenome is reported as an uninformative
screen rather than as a working contamination check.
"""

from __future__ import annotations

import importlib
import importlib.util
import os
import platform as platform_module
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from . import config, shell
from .records import DatabaseVersion
from .utils import (
    MjolnirError,
    PathLike,
    conda_package,
    cpu_count,
    free_bytes,
    human_bytes,
    sha256sum,
)

# ---------------------------------------------------------------------------
# Requirement levels
# ---------------------------------------------------------------------------

#: Without it, nothing that reads sequence works.
LEVEL_REQUIRED = "required"
#: One member of a group must be present; which one is a real choice, not a
#: degradation (bwa-mem2 or bwa; skani or mash).
LEVEL_ALTERNATIVE = "alternative"
#: Absence costs a named capability, and the report says which.
LEVEL_OPTIONAL = "optional"

LEVELS: Tuple[str, ...] = (LEVEL_REQUIRED, LEVEL_ALTERNATIVE, LEVEL_OPTIONAL)


def _register(name: str, value: Any, source: str, unit: str = "",
              note: str = "") -> Any:
    """File a doctor threshold in config.py's registry and return its value."""
    return config._define(name, value, source, unit=unit, note=note)


#: SOURCE: Mjolnir policy. Files larger than this are listed without a checksum.
#: The registry checksums a catalogue so two installations can be compared; a
#: Kraken2 index is tens of gigabytes and hashing it would turn ``doctor`` into
#: the slowest command in the tool.
CHECKSUM_MAX_BYTES = _register(
    "doctor_checksum_max_bytes", 64 << 20, config.SRC_POLICY, unit="bytes",
    note="largest database file doctor will checksum inline")

#: Seconds allowed for the model host to answer. Plumbing, not science: the LLM
#: is optional by design and a slow host must not delay an environment report.
LLM_PROBE_TIMEOUT_SECONDS = 2.0


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ToolSpec:
    """One external program, why Mjolnir wants it, and how badly."""

    label: str
    #: Executable names to look for, in order of preference. More than one
    #: because upstream renames happen: IQ-TREE ships as iqtree3, iqtree2 or
    #: iqtree depending on the channel, and Clair3's entry point is a shell
    #: script, not a binary called "clair3".
    binaries: Tuple[str, ...]
    purpose: str
    level: str
    group: str = ""
    note: str = ""


#: Every tool Mjolnir can drive. Ordered by the stage that uses it, because the
#: report reads top to bottom as the pipeline runs.
TOOLS: Tuple[ToolSpec, ...] = (
    ToolSpec("samtools", ("samtools",),
             "BAM handling, indexing and the direct pileup used for barcode "
             "sites and catalogue positions",
             LEVEL_REQUIRED),
    ToolSpec("bcftools", ("bcftools",),
             "variant calling on Illumina, and VCF normalisation on every "
             "platform",
             LEVEL_REQUIRED),
    ToolSpec("minimap2", ("minimap2",),
             "ONT read mapping (-x map-ont) and assembly-to-reference "
             "alignment for FASTA input",
             LEVEL_REQUIRED),
    ToolSpec("bwa-mem2", ("bwa-mem2",),
             "Illumina read mapping",
             LEVEL_ALTERNATIVE, group="illumina-mapper"),
    ToolSpec("bwa", ("bwa",),
             "Illumina read mapping, the fallback when bwa-mem2 is absent",
             LEVEL_ALTERNATIVE, group="illumina-mapper",
             note="bwa-mem2 and bwa produce equivalent alignments; the "
                  "fallback is a speed difference, not a sensitivity one"),
    ToolSpec("skani", ("skani",),
             "whole-genome ANI for species identification",
             LEVEL_ALTERNATIVE, group="ani"),
    ToolSpec("mash", ("mash",),
             "sketch-based ANI for species identification, and the cheap "
             "first pass before skani",
             LEVEL_ALTERNATIVE, group="ani"),
    ToolSpec("clair3", ("run_clair3.sh",),
             "the preferred ONT variant caller",
             LEVEL_OPTIONAL,
             note="without it ONT calling falls back to bcftools, which is "
                  "specifically weak on ONT indels; the report states which "
                  "caller produced the numbers"),
    ToolSpec("freebayes", ("freebayes",),
             "alternative Illumina variant caller, for cross-checking a call",
             LEVEL_OPTIONAL),
    ToolSpec("kraken2", ("kraken2",),
             "read composition screen",
             LEVEL_OPTIONAL,
             note="only informative with a mycobacterial pangenome index; "
                  "measured sensitivity for M. tuberculosis reads against a "
                  "standard index is 0.0731"),
    ToolSpec("seqkit", ("seqkit",),
             "FASTA/FASTQ manipulation in the database build and subsetting "
             "steps",
             LEVEL_OPTIONAL),
    ToolSpec("iqtree", ("iqtree3", "iqtree2", "iqtree"),
             "maximum-likelihood tree over the cohort's joint variant "
             "alignment",
             LEVEL_OPTIONAL,
             note="cohort distances and clusters do not need it; the tree "
                  "figure does"),
)


@dataclass
class ToolReport:
    """One tool's presence, path and version — or the reason none of those."""

    label: str
    level: str
    purpose: str
    found: bool = False
    binary: str = ""
    path: str = ""
    version: str = ""
    group: str = ""
    note: str = ""
    install: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "label": self.label, "level": self.level, "purpose": self.purpose,
            "found": self.found, "binary": self.binary, "path": self.path,
            "version": self.version, "group": self.group, "note": self.note,
            "install": self.install,
        }


def check_tool(spec: ToolSpec) -> ToolReport:
    """Probe one tool. Never raises; an unprobeable tool is a report, not a stop."""
    report = ToolReport(label=spec.label, level=spec.level, purpose=spec.purpose,
                        group=spec.group, note=spec.note,
                        install="conda install -c conda-forge -c bioconda {0}".format(
                            conda_package(spec.binaries[0])))
    for binary in spec.binaries:
        try:
            path = shell.tool_path(binary, required=False)
        except (OSError, MjolnirError):
            # tool_path with required=False does not raise for absence; this
            # catches an unreadable PATH entry, which is a broken environment
            # rather than a missing tool, and is worth saying out loud.
            report.note = "could not search PATH for {0}".format(binary)
            continue
        if not path:
            continue
        report.found = True
        report.binary = binary
        report.path = path
        try:
            report.version = shell.tool_version(binary) or shell.VERSION_UNKNOWN
        except (OSError, MjolnirError):
            report.version = shell.VERSION_UNKNOWN
        return report
    return report


def check_tools(specs: Sequence[ToolSpec] = TOOLS) -> List[ToolReport]:
    """Probe every tool. Always returns one row per spec."""
    return [check_tool(spec) for spec in specs]


def _alternative_groups(reports: Sequence[ToolReport]) -> Dict[str, List[ToolReport]]:
    groups: Dict[str, List[ToolReport]] = {}
    for report in reports:
        if report.level == LEVEL_ALTERNATIVE and report.group:
            groups.setdefault(report.group, []).append(report)
    return groups


def missing_required_tools(reports: Sequence[ToolReport]) -> List[ToolReport]:
    """Required tools that are absent, plus every member of an empty group."""
    missing = [r for r in reports if r.level == LEVEL_REQUIRED and not r.found]
    for members in _alternative_groups(reports).values():
        if not any(m.found for m in members):
            missing.extend(members)
    return missing


def tool_problems(reports: Sequence[ToolReport]) -> List[str]:
    """One problem line per unmet requirement.

    An alternative group produces a single line naming the choice rather than
    one line per member: telling an operator that both bwa-mem2 and bwa are
    "required" when either would do is how a doctor's output stops being read.
    """
    problems: List[str] = []
    for report in reports:
        if report.level == LEVEL_REQUIRED and not report.found:
            problems.append("required tool missing: {0} ({1})\n    {2}".format(
                report.label, report.purpose, report.install))
    for group, members in sorted(_alternative_groups(reports).items()):
        if any(m.found for m in members):
            continue
        problems.append(
            "no {0} installed: one of {1} is required ({2})\n    {3}".format(
                group, " or ".join(m.label for m in members),
                members[0].purpose,
                "\n    ".join(m.install for m in members)))
    return problems


# ---------------------------------------------------------------------------
# Python dependencies
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PythonDep:
    """One import Mjolnir needs, and what stops working without it."""

    module: str
    distribution: str
    purpose: str
    level: str
    #: The pyproject extra that installs it, for the error message.
    extra: str = ""


PYTHON_DEPS: Tuple[PythonDep, ...] = (
    PythonDep("pandas", "pandas",
              "catalogue tables and the joint variant table", LEVEL_REQUIRED),
    PythonDep("numpy", "numpy",
              "depth, breadth and evenness arithmetic over the coverage array",
              LEVEL_REQUIRED),
    PythonDep("openpyxl", "openpyxl",
              "reading the WHO catalogue, which is distributed only as .xlsx - "
              "the repository's .txt master file is missing a drug and four "
              "genes, so it is not a substitute",
              LEVEL_REQUIRED, extra="excel"),
    PythonDep("pysam", "pysam",
              "direct pileup at barcode and catalogue positions without "
              "shelling out per site",
              LEVEL_OPTIONAL, extra="reads"),
    PythonDep("scipy", "scipy",
              "hierarchical clustering and the cohort dendrogram",
              LEVEL_OPTIONAL, extra="cluster"),
    PythonDep("reportlab", "reportlab",
              "the PDF deliverable; without it the report is HTML, TSV and "
              "JSON only",
              LEVEL_OPTIONAL, extra="report"),
    PythonDep("matplotlib", "matplotlib",
              "the drug grid, coverage strip, allele-fraction plot and cluster "
              "dendrogram",
              LEVEL_OPTIONAL, extra="report"),
)


@dataclass
class PythonDepReport:
    module: str
    level: str
    purpose: str
    found: bool = False
    version: str = ""
    extra: str = ""
    note: str = ""
    install: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "module": self.module, "level": self.level, "purpose": self.purpose,
            "found": self.found, "version": self.version, "extra": self.extra,
            "note": self.note, "install": self.install,
        }


def _distribution_version(distribution: str) -> str:
    """Version from installed metadata, without importing the package.

    Importing pandas to read ``pandas.__version__`` costs the better part of a
    second and can fail for reasons that have nothing to do with the version, so
    the metadata is read instead and the import is only checked for existence.
    """
    from importlib import metadata

    try:
        return metadata.version(distribution)
    except metadata.PackageNotFoundError:
        # Importable but with no distribution metadata: a vendored copy, or a
        # source checkout on PYTHONPATH. The module is there; its version is
        # genuinely unknown, and saying so is the honest answer.
        return ""
    except (OSError, ValueError) as exc:
        return "version metadata unreadable: {0}".format(exc)


def check_python_deps(deps: Sequence[PythonDep] = PYTHON_DEPS) -> List[PythonDepReport]:
    """Probe every Python dependency by spec lookup, not by importing it."""
    out: List[PythonDepReport] = []
    for dep in deps:
        install = "pip install {0}".format(dep.distribution)
        if dep.extra:
            install = "pip install 'mjolnir-myco[{0}]'".format(dep.extra)
        report = PythonDepReport(module=dep.module, level=dep.level,
                                 purpose=dep.purpose, extra=dep.extra,
                                 install=install)
        try:
            spec = importlib.util.find_spec(dep.module)
        except (ImportError, ValueError) as exc:
            report.note = "import system could not resolve {0}: {1}".format(
                dep.module, exc)
            out.append(report)
            continue
        if spec is None:
            out.append(report)
            continue
        report.found = True
        report.version = _distribution_version(dep.distribution) or "installed"
        out.append(report)
    return out


# ---------------------------------------------------------------------------
# Databases
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DatabaseSpec:
    """One database Mjolnir consults, and how to get it.

    ``patterns`` are globs relative to the database root, tried in order and to
    a bounded depth. The layout is ``db/fetch.py``'s to decide, so the doctor
    looks for the file by name rather than assuming a directory structure.
    """

    name: str
    patterns: Tuple[str, ...]
    purpose: str
    level: str
    licence: str = ""
    citation: str = ""
    fetch_hint: str = "mjolnir db fetch"


#: The design's §12 table. This is the fallback: if ``db/registry.py`` is
#: importable and exposes a database list, that is authoritative and is used
#: instead, since it is the module that records what was actually fetched.
DATABASES: Tuple[DatabaseSpec, ...] = (
    DatabaseSpec(
        "WHO catalogue v2 (xlsx)",
        ("WHO-UCN-TB-2023.7-eng.xlsx", "*WHO*2023.7*.xlsx", "who/*.xlsx"),
        "the graded resistance catalogue and the anchor of the consensus rule",
        LEVEL_REQUIRED,
        licence="ODC-By v1.0 (redistributable with attribution)",
        citation=config.SRC_WHO_V2,
        fetch_hint="mjolnir db fetch who"),
    DatabaseSpec(
        "H37Rv reference (NC_000962.3)",
        ("NC_000962.3.fasta", "*NC_000962*.fa*", "reference/*NC_000962*"),
        "the coordinate system for every MTBC variant, barcode site and mask "
        "interval in this tool",
        LEVEL_REQUIRED,
        licence="public",
        citation=config.SRC_NC_000962,
        fetch_hint="mjolnir db fetch reference"),
    DatabaseSpec(
        "tbdb mutations.csv",
        ("mutations.csv", "tbdb/mutations.csv"),
        "the tbdb catalogue arm of the three-catalogue consensus",
        LEVEL_OPTIONAL,
        licence="verify at fetch time",
        citation=config.SRC_TBDB,
        fetch_hint="mjolnir db fetch tbdb"),
    DatabaseSpec(
        "tbdb barcode.bed",
        ("barcode.bed", "tbdb/barcode.bed"),
        "the lineage-defining SNP barcode; without it no MTBC lineage, "
        "sublineage, BCG or animal-lineage call can be made at all",
        LEVEL_OPTIONAL,
        licence="verify at fetch time",
        citation=config.SRC_TBDB,
        fetch_hint="mjolnir db fetch tbdb"),
    DatabaseSpec(
        "tbdb mask.bed",
        ("mask.bed", "tbdb/mask.bed"),
        "the repetitive and error-prone regions excluded from cohort SNP "
        "distances; masking is mandatory and the mask used is printed in the "
        "report",
        LEVEL_OPTIONAL,
        licence="verify at fetch time",
        citation=config.SRC_TBDB,
        fetch_hint="mjolnir db fetch tbdb"),
    DatabaseSpec(
        "MTBseq resistance lists",
        ("MTB_Resistance_Mediating.txt", "mtbseq/MTB_Resistance_Mediating.txt"),
        "the MTBseq catalogue arm of the consensus, and the comparison "
        "baseline for a lab migrating from MTBseq",
        LEVEL_OPTIONAL,
        licence="GPL-3.0 (MTBseq)",
        citation="ngs-fzb/MTBseq_source var/res/",
        fetch_hint="mjolnir db fetch mtbseq"),
    DatabaseSpec(
        "mycobacterial ANI reference set",
        ("ani/*.msh", "ani/*.sketch", "ani/references.txt", "*.msh"),
        "species identification; without it Mjolnir cannot name a species and "
        "will not guess one from a read classifier",
        LEVEL_OPTIONAL,
        licence="per source assembly",
        citation=config.SRC_ANI_SPECIES,
        fetch_hint="mjolnir db fetch ani"),
)


@dataclass
class DatabaseReport:
    """One database as found on disk, or the command that would fetch it."""

    name: str
    level: str
    purpose: str
    present: bool = False
    #: Present is not the same as fit for purpose. A standard Kraken2 index is
    #: on disk and cannot answer the question it would be used for, so that
    #: distinction is carried rather than collapsed into a tick.
    usable: Optional[bool] = None
    version: Optional[DatabaseVersion] = None
    fetch_hint: str = ""
    note: str = ""

    @property
    def ready(self) -> bool:
        return self.present and (self.usable is not False)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name, "level": self.level, "purpose": self.purpose,
            "present": self.present, "usable": self.usable,
            "ready": self.ready,
            "version": self.version.to_dict() if self.version else None,
            "fetch_hint": self.fetch_hint, "note": self.note,
        }


def _find_in_db(db_dir: Path, patterns: Sequence[str]) -> Optional[Path]:
    """First match for any pattern, searched to three levels.

    Bounded rather than recursive on purpose: a database root may contain a
    Kraken2 index of several million files, and walking it to find a 5 MB
    spreadsheet is not a reasonable thing for an environment check to do.
    """
    for pattern in patterns:
        for prefix in ("", "*/", "*/*/"):
            try:
                matches = sorted(db_dir.glob(prefix + pattern))
            except OSError:
                continue
            for match in matches:
                if match.is_file():
                    return match
    return None


def _describe_file(path: Path, spec: DatabaseSpec) -> DatabaseVersion:
    """Build the registry record for a database file that is present."""
    checksum = ""
    note = ""
    try:
        size = path.stat().st_size
        if size <= CHECKSUM_MAX_BYTES:
            checksum = sha256sum(path)
        else:
            note = "not checksummed: {0} exceeds the {1} inline limit".format(
                human_bytes(size), human_bytes(CHECKSUM_MAX_BYTES))
    except OSError as exc:
        note = "could not read {0}: {1}".format(path, exc)
    return DatabaseVersion(
        name=spec.name,
        version="unknown",
        checksum=checksum,
        path=str(path),
        licence=spec.licence,
        citation=spec.citation,
        note=note,
    )


def _registry_databases() -> Optional[List[DatabaseReport]]:
    """Ask ``db/registry.py`` for the authoritative list, if it exists yet.

    The registry module owns what was actually fetched, including versions and
    checksums recorded at fetch time, which is strictly better information than
    a filename glob. Its API is not fixed at the time this was written, so three
    plausible shapes are accepted and anything else falls through to the glob —
    reported as a fallback rather than passed off as the registry's answer.
    """
    try:
        registry = importlib.import_module("mjolnir.db.registry")
    except ImportError:
        return None
    for attribute in ("doctor_reports", "all_databases", "entries"):
        function = getattr(registry, attribute, None)
        if callable(function):
            try:
                raw = function()
            # Deliberately broad, and only here. This calls into a module the
            # doctor does not own, at exactly the moment the operator is asking
            # what is broken; whatever it raises becomes a visible finding
            # rather than a traceback that hides the rest of the report.
            except Exception as exc:  # noqa: BLE001 - see comment above
                return [DatabaseReport(
                    name="database registry",
                    level=LEVEL_REQUIRED,
                    purpose="the authoritative record of what was fetched",
                    present=False,
                    note="mjolnir.db.registry.{0}() raised {1}: {2}".format(
                        attribute, exc.__class__.__name__, exc),
                    fetch_hint="mjolnir db fetch")]
            return _normalise_registry(raw)
    return None


def _normalise_registry(raw: Any) -> List[DatabaseReport]:
    """Turn whatever the registry returned into DatabaseReport rows."""
    out: List[DatabaseReport] = []
    items = raw.values() if isinstance(raw, dict) else raw
    for item in items or []:
        if isinstance(item, DatabaseReport):
            out.append(item)
            continue
        get = item.get if isinstance(item, dict) else (
            lambda key, default=None: getattr(item, key, default))
        path = get("path", "") or ""
        version = DatabaseVersion(
            name=str(get("name", "unnamed")),
            version=str(get("version", "unknown")),
            checksum=str(get("checksum", "") or ""),
            path=str(path),
            licence=str(get("licence", "") or ""),
            citation=str(get("citation", "") or ""),
            url=str(get("url", "") or ""),
            fetched=str(get("fetched", "") or ""),
            note=str(get("note", "") or ""),
        )
        present = bool(path) and Path(str(path)).exists()
        out.append(DatabaseReport(
            name=version.name,
            level=LEVEL_REQUIRED if get("required", False) else LEVEL_OPTIONAL,
            purpose=str(get("purpose", "") or ""),
            present=present,
            version=version,
            fetch_hint=str(get("fetch_hint", "mjolnir db fetch")),
        ))
    return out


def check_databases(db_dir: PathLike,
                    specs: Sequence[DatabaseSpec] = DATABASES
                    ) -> List[DatabaseReport]:
    """Report every database: present with its checksum, or absent with its fetch line."""
    from_registry = _registry_databases()
    if from_registry:
        return from_registry

    root = Path(db_dir).expanduser()
    reports: List[DatabaseReport] = []
    if not root.exists():
        for spec in specs:
            reports.append(DatabaseReport(
                name=spec.name, level=spec.level, purpose=spec.purpose,
                present=False, fetch_hint=spec.fetch_hint,
                note="database root {0} does not exist; set ${1} or run "
                     "'mjolnir db fetch'".format(root, config.DB_ENV_VAR)))
        return reports

    for spec in specs:
        found = _find_in_db(root, spec.patterns)
        if found is None:
            reports.append(DatabaseReport(
                name=spec.name, level=spec.level, purpose=spec.purpose,
                present=False, fetch_hint=spec.fetch_hint,
                note="not found under {0} (looked for {1})".format(
                    root, ", ".join(spec.patterns))))
        else:
            reports.append(DatabaseReport(
                name=spec.name, level=spec.level, purpose=spec.purpose,
                present=True, version=_describe_file(found, spec),
                fetch_hint=spec.fetch_hint,
                note="version and checksum recorded here are of the file on "
                     "disk; db/registry.py records what was fetched"))
    return reports


def check_kraken2_index(db_path: Optional[PathLike]) -> DatabaseReport:
    """The Kraken2 index, and whether it is a mycobacterial screen at all.

    Separated from the other databases because the answer is not "present or
    absent". A standard or capped index is present and useless for this
    organism: measured sensitivity for *M. tuberculosis* reads is 0.0731, so
    reporting a clean screen from one would be an invented result.
    """
    report = DatabaseReport(
        name="Kraken2 index",
        level=LEVEL_OPTIONAL,
        purpose="read composition screen",
        fetch_hint="set ${0} to a mycobacterial pangenome index".format(
            config.KRAKEN2_DB_ENV_VAR),
    )
    if not db_path:
        report.note = (
            "no index configured; the read composition screen will not run, "
            "and the report will say the screen was not performed rather than "
            "that the sample was clean")
        return report
    path = Path(db_path).expanduser()
    if not path.exists():
        report.note = "configured index {0} does not exist".format(path)
        return report
    report.present = True
    try:
        informative, why = config.kraken2_index_informative(path)
    except MjolnirError as exc:
        report.note = str(exc)
        report.usable = False
        return report
    report.usable = informative
    report.version = DatabaseVersion(name="Kraken2 index", path=str(path),
                                     note=why)
    # `why` already carries the refusal text when the index is uninformative;
    # config.kraken2_index_informative composes it. Repeating it here would
    # print the same paragraph twice.
    report.note = why
    return report


# ---------------------------------------------------------------------------
# Capabilities
# ---------------------------------------------------------------------------

@dataclass
class Capability:
    """Something the pipeline can or cannot do, and what is stopping it."""

    name: str
    available: bool
    missing: List[str] = field(default_factory=list)
    consequence: str = ""
    note: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {"name": self.name, "available": self.available,
                "missing": list(self.missing), "consequence": self.consequence,
                "note": self.note}


def _present(reports: Sequence[ToolReport], label: str) -> bool:
    return any(r.label == label and r.found for r in reports)


def _dep_present(deps: Sequence[PythonDepReport], module: str) -> bool:
    return any(d.module == module and d.found for d in deps)


def _db_present(dbs: Sequence[DatabaseReport], name: str) -> bool:
    return any(name.lower() in d.name.lower() and d.ready for d in dbs)


def capabilities(tools: Sequence[ToolReport],
                 deps: Sequence[PythonDepReport],
                 databases: Sequence[DatabaseReport],
                 kraken2: Optional[DatabaseReport] = None) -> List[Capability]:
    """Roll the inventory up into what can actually be run.

    This is the part an operator reads. "freebayes: missing" means nothing on
    its own; "Illumina reads: available" and "MTBC lineage: unavailable, needs
    barcode.bed" are decisions they can act on.
    """
    out: List[Capability] = []

    def need(name: str, requirements: Sequence[Tuple[bool, str]],
             consequence: str, note: str = "") -> None:
        missing = [what for ok, what in requirements if not ok]
        out.append(Capability(name=name, available=not missing, missing=missing,
                              consequence=consequence, note=note))

    need("Illumina read input",
         [(_present(tools, "bwa-mem2") or _present(tools, "bwa"),
           "bwa-mem2 or bwa"),
          (_present(tools, "samtools"), "samtools"),
          (_present(tools, "bcftools"), "bcftools")],
         "paired Illumina FASTQ cannot be mapped or called")

    ont_caller = "clair3" if _present(tools, "clair3") else (
        "bcftools" if _present(tools, "bcftools") else "")
    need("ONT read input",
         [(_present(tools, "minimap2"), "minimap2"),
          (_present(tools, "samtools"), "samtools"),
          (bool(ont_caller), "clair3 or bcftools")],
         "ONT FASTQ cannot be mapped or called",
         note="" if _present(tools, "clair3") else
              "Clair3 is absent, so ONT calling would fall back to bcftools, "
              "which is specifically weak on ONT indels; the report states "
              "which caller was used")

    need("assembly (FASTA) input",
         [(_present(tools, "minimap2"), "minimap2")],
         "assemblies cannot be aligned to the reference",
         note=config.FASTA_CAPABILITY_LOSS)

    need("species identification",
         [(_present(tools, "skani") or _present(tools, "mash"), "skani or mash"),
          (_db_present(databases, "ANI reference"),
           "the mycobacterial ANI reference set")],
         "no species can be named. Mjolnir will not substitute a taxonomic "
         "read classifier for an ANI call: in current NCBI taxonomy the MTBC "
         "members are not at species rank, so a classifier row naming one is "
         "not a species identification")

    need("MTBC lineage and sublineage",
         [(_db_present(databases, "barcode"), "tbdb barcode.bed"),
          (_present(tools, "samtools"), "samtools")],
         "no lineage, sublineage, BCG or animal-lineage call, and the "
         "intrinsic pyrazinamide resistance of BCG would go unflagged")

    need("resistance calling (WHO anchor)",
         [(_db_present(databases, "WHO catalogue"), "the WHO v2 xlsx"),
          (_dep_present(deps, "openpyxl"), "openpyxl")],
         "no graded resistance calls. Absence of a catalogue is reported as "
         "absence of evidence, never as susceptibility")

    need("read composition screen",
         [(_present(tools, "kraken2"), "kraken2"),
          (bool(kraken2 and kraken2.ready), "a mycobacterial pangenome index")],
         "no read-level composition screen runs. The heterozygosity and "
         "coverage-based contamination measures are unaffected; what is lost "
         "is the non-target read fraction, and its absence is reported as "
         "absence rather than as a clean sample",
         note=(kraken2.note if kraken2 and kraken2.usable is False else ""))

    need("cohort distances and clusters",
         [(_db_present(databases, "mask"), "tbdb mask.bed"),
          (_present(tools, "bcftools"), "bcftools")],
         "pairwise SNP distances cannot be masked, and an unmasked distance "
         "over repetitive and error-prone regions is not comparable to a "
         "published threshold")

    need("cohort tree figure",
         [(_present(tools, "iqtree"), "iqtree")],
         "the cohort report loses its maximum-likelihood tree; distances and "
         "clusters are unaffected")

    need("PDF report",
         [(_dep_present(deps, "reportlab"), "reportlab"),
          (_dep_present(deps, "matplotlib"), "matplotlib")],
         "output is HTML, TSV and JSON only")

    return out


# ---------------------------------------------------------------------------
# The model host
# ---------------------------------------------------------------------------

def llm_reachable(host: str, timeout: float = LLM_PROBE_TIMEOUT_SECONDS
                  ) -> Tuple[bool, str]:
    """Whether the configured model host answers, and what it said.

    A False here is never fatal. The design's fourth house rule is that Mjolnir
    runs without the model: every gate takes its declared default and the report
    states that the interpretation is rule-only. This probe exists so the
    operator learns that before the run, not from the finished PDF.
    """
    if not host:
        return False, (
            "no model host configured (${0}); the report will be rule-only, "
            "which is a supported mode and is stated in the output".format(
                config.LLM_HOST_ENV_VAR))
    import urllib.error
    import urllib.request

    base = host.rstrip("/")
    for endpoint, flavour in (("/api/tags", "ollama"),
                              ("/v1/models", "OpenAI-compatible")):
        url = base + endpoint
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                if 200 <= response.status < 300:
                    return True, "{0} answers the {1} API".format(host, flavour)
        except (urllib.error.URLError, OSError, ValueError):
            continue
    return False, (
        "{0} did not answer on /api/tags or /v1/models within {1:g}s; the "
        "report will be rule-only".format(host, timeout))


# ---------------------------------------------------------------------------
# The whole picture
# ---------------------------------------------------------------------------

@dataclass
class Diagnosis:
    """Everything ``mjolnir doctor`` found, in one object."""

    tools: List[ToolReport] = field(default_factory=list)
    python_deps: List[PythonDepReport] = field(default_factory=list)
    databases: List[DatabaseReport] = field(default_factory=list)
    kraken2: Optional[DatabaseReport] = None
    capabilities: List[Capability] = field(default_factory=list)
    environment: Dict[str, Any] = field(default_factory=dict)
    unverified_thresholds: List[str] = field(default_factory=list)
    llm_ok: bool = False
    llm_note: str = ""
    problems: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """True when nothing required is missing.

        Optional tools and optional databases do not make a diagnosis fail —
        they cost capabilities, which are listed separately — so that an
        Illumina-only lab is not told its environment is broken because it has
        no Clair3.
        """
        return not self.problems

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ok": self.ok,
            "environment": dict(self.environment),
            "tools": [t.to_dict() for t in self.tools],
            "python_deps": [d.to_dict() for d in self.python_deps],
            "databases": [d.to_dict() for d in self.databases],
            "kraken2": self.kraken2.to_dict() if self.kraken2 else None,
            "capabilities": [c.to_dict() for c in self.capabilities],
            "unverified_thresholds": list(self.unverified_thresholds),
            "llm": {"reachable": self.llm_ok, "note": self.llm_note},
            "problems": list(self.problems),
        }

    def tool_versions(self) -> Dict[str, str]:
        """Installed tools and versions, for the methods annex.

        Distinct from ``shell.captured_versions()``, which lists what actually
        ran. A report must never claim to have used a caller it merely had
        installed.
        """
        return dict((t.label, t.version) for t in self.tools if t.found)


def _environment(cfg: config.Config) -> Dict[str, Any]:
    env: Dict[str, Any] = {
        "mjolnir_version": _mjolnir_version(),
        "python": sys.version.split()[0],
        "python_executable": sys.executable,
        "platform": platform_module.platform(),
        "cpus_available": cpu_count(),
        "threads_configured": cfg.threads,
        "db_dir": str(cfg.db_dir),
        "out_dir": str(cfg.out_dir),
        "profile": cfg.profile,
    }
    for label, path in (("db_dir_free", cfg.db_dir), ("out_dir_free", cfg.out_dir)):
        # Neither directory need exist yet, so walk up to the first ancestor
        # that does: the free space that matters is the filesystem they will be
        # created on, and "unknown" here would be a needless blank in the
        # report of a machine that is simply not set up yet.
        target = Path(path).expanduser().absolute()
        while not target.exists() and target != target.parent:
            target = target.parent
        free = free_bytes(target)
        env[label] = human_bytes(free) if free else "unknown"
    for name in (config.DB_ENV_VAR, config.KRAKEN2_DB_ENV_VAR,
                 config.LLM_HOST_ENV_VAR, config.LLM_MODEL_ENV_VAR,
                 config.THREADS_ENV_VAR):
        env[name] = os.environ.get(name, "")
    return env


def _mjolnir_version() -> str:
    try:
        package = importlib.import_module("mjolnir")
        return str(getattr(package, "__version__", "unknown"))
    except ImportError:
        return "unknown"


def diagnose(cfg: Optional[config.Config] = None) -> Diagnosis:
    """Probe everything and return the result. Never raises.

    Each section is guarded independently: a database root that cannot be read
    must not cost the operator the tool inventory, because the two problems are
    usually fixed in the same sitting and they should be shown together.
    """
    if cfg is None:
        try:
            cfg = config.Config()
        except MjolnirError as exc:
            # A configuration that will not even construct — a non-numeric
            # $MJOLNIR_THREADS, say — is itself the finding. Report it as the
            # single problem rather than probing an environment the run could
            # not use anyway.
            diagnosis = Diagnosis()
            diagnosis.problems.append("configuration is invalid: {0}".format(exc))
            return diagnosis

    diagnosis = Diagnosis()

    try:
        diagnosis.environment = _environment(cfg)
    except (OSError, MjolnirError) as exc:
        diagnosis.problems.append("could not describe the environment: {0}".format(exc))

    diagnosis.tools = check_tools()
    diagnosis.python_deps = check_python_deps()

    try:
        diagnosis.databases = check_databases(cfg.db_dir)
    except (OSError, MjolnirError) as exc:
        diagnosis.problems.append(
            "could not inspect the database root {0}: {1}".format(cfg.db_dir, exc))

    try:
        diagnosis.kraken2 = check_kraken2_index(cfg.kraken2_db)
    except (OSError, MjolnirError) as exc:
        diagnosis.problems.append("could not inspect the Kraken2 index: {0}".format(exc))

    diagnosis.capabilities = capabilities(
        diagnosis.tools, diagnosis.python_deps, diagnosis.databases,
        diagnosis.kraken2)

    try:
        diagnosis.unverified_thresholds = [t.name for t in config.unverified()]
    except MjolnirError as exc:  # pragma: no cover - registry is built at import
        diagnosis.problems.append("threshold registry unreadable: {0}".format(exc))

    if cfg.use_llm:
        diagnosis.llm_ok, diagnosis.llm_note = llm_reachable(cfg.llm_host)
    else:
        diagnosis.llm_note = "model interpretation disabled for this run"

    diagnosis.problems.extend(tool_problems(diagnosis.tools))
    for dep in diagnosis.python_deps:
        if dep.level == LEVEL_REQUIRED and not dep.found:
            diagnosis.problems.append(
                "required Python package missing: {0} ({1})\n    {2}".format(
                    dep.module, dep.purpose, dep.install))
    for database in diagnosis.databases:
        if database.level == LEVEL_REQUIRED and not database.ready:
            diagnosis.problems.append(
                "required database missing: {0} ({1})\n    {2}".format(
                    database.name, database.purpose, database.fetch_hint))
    return diagnosis


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

_MARK_PRESENT = "ok"
_MARK_ABSENT = "--"


def _status_mark(present: bool, level: str) -> str:
    if present:
        return _MARK_PRESENT
    return "MISSING" if level == LEVEL_REQUIRED else _MARK_ABSENT


def render(diagnosis: Diagnosis, verbose: bool = False) -> str:
    """The environment report as text, for the CLI.

    Returns a string rather than printing, so the same content can go to a log,
    to the JSON artefact and into the report's methods annex without being
    rebuilt in three places.
    """
    lines: List[str] = []
    env = diagnosis.environment

    lines.append("Mjolnir {0} - environment".format(env.get("mjolnir_version", "?")))
    lines.append("  python      {0}  ({1})".format(
        env.get("python", "?"), env.get("python_executable", "?")))
    lines.append("  platform    {0}".format(env.get("platform", "?")))
    lines.append("  cpus        {0} available, {1} configured".format(
        env.get("cpus_available", "?"), env.get("threads_configured", "?")))
    lines.append("  database    {0}  ({1} free)".format(
        env.get("db_dir", "?"), env.get("db_dir_free", "?")))
    lines.append("  output      {0}  ({1} free)".format(
        env.get("out_dir", "?"), env.get("out_dir_free", "?")))

    lines.append("")
    lines.append("External tools")
    for tool in diagnosis.tools:
        lines.append("  {0:<8} {1:<12} {2:<11} {3}".format(
            _status_mark(tool.found, tool.level), tool.label, tool.level,
            tool.version if tool.found else tool.install))
        if verbose or not tool.found:
            lines.append("           {0}".format(tool.purpose))
        if tool.note and (verbose or not tool.found):
            lines.append("           note: {0}".format(tool.note))
        if verbose and tool.found:
            lines.append("           {0}".format(tool.path))

    lines.append("")
    lines.append("Python packages")
    for dep in diagnosis.python_deps:
        lines.append("  {0:<8} {1:<12} {2:<11} {3}".format(
            _status_mark(dep.found, dep.level), dep.module, dep.level,
            dep.version if dep.found else dep.install))
        if verbose or not dep.found:
            lines.append("           {0}".format(dep.purpose))
        if dep.note:
            lines.append("           note: {0}".format(dep.note))

    lines.append("")
    lines.append("Databases")
    for database in diagnosis.databases:
        detail = ""
        if database.ready and database.version:
            detail = database.version.path
            if database.version.checksum:
                detail += "  sha256:{0}".format(database.version.checksum[:12])
        else:
            detail = database.fetch_hint
        lines.append("  {0:<8} {1:<32} {2:<11} {3}".format(
            _status_mark(database.ready, database.level), database.name,
            database.level, detail))
        if verbose or not database.ready:
            lines.append("           {0}".format(database.purpose))
        if database.note and (verbose or not database.ready):
            lines.append("           note: {0}".format(database.note))
        if verbose and database.ready and database.version:
            if database.version.licence:
                lines.append("           licence: {0}".format(database.version.licence))

    if diagnosis.kraken2 is not None:
        kraken2 = diagnosis.kraken2
        lines.append("")
        lines.append("Contamination screen")
        if kraken2.present and kraken2.usable is False:
            mark = "REFUSED"
        else:
            mark = _status_mark(kraken2.present, LEVEL_OPTIONAL)
        lines.append("  {0:<8} {1}".format(mark, kraken2.name))
        lines.append("           {0}".format(kraken2.note))

    lines.append("")
    lines.append("Capabilities")
    for capability in diagnosis.capabilities:
        lines.append("  {0:<8} {1}".format(
            "ok" if capability.available else _MARK_ABSENT, capability.name))
        if not capability.available:
            lines.append("           needs: {0}".format(", ".join(capability.missing)))
            lines.append("           without it: {0}".format(capability.consequence))
        if capability.note and (verbose or not capability.available):
            lines.append("           note: {0}".format(capability.note))

    lines.append("")
    lines.append("Interpretation model")
    lines.append("  {0:<8} {1}".format("ok" if diagnosis.llm_ok else _MARK_ABSENT,
                                       diagnosis.llm_note))

    if diagnosis.unverified_thresholds:
        lines.append("")
        lines.append("Thresholds whose citation has not been checked against the "
                     "primary document ({0} of {1}):".format(
                         len(diagnosis.unverified_thresholds),
                         len(config.THRESHOLDS)))
        for name in diagnosis.unverified_thresholds:
            lines.append("  {0}".format(config.threshold(name).describe()))

    lines.append("")
    if diagnosis.problems:
        lines.append("NOT READY - {0} problem(s):".format(len(diagnosis.problems)))
        for problem in diagnosis.problems:
            lines.append("  - {0}".format(problem))
    else:
        unavailable = [c.name for c in diagnosis.capabilities if not c.available]
        if unavailable:
            lines.append("READY, with {0} capability/capabilities unavailable: "
                         "{1}".format(len(unavailable), ", ".join(unavailable)))
        else:
            lines.append("READY - every capability is available.")
    return "\n".join(lines)
