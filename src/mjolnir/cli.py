"""Mjolnir's command line.

Five commands, and the boundary between them is what the operator is asking
for rather than what the code does: ``run`` for one sample, ``cohort`` for a set
compared against one another, ``report`` to render a finished run again, ``db``
for the reference data with its licence and citation printed, and ``doctor``
for the environment before any of it.

Three things here are deliberate.

**A threshold that changed is announced.** Every flag that overrides a published
number goes through :meth:`~mjolnir.config.Config.set_explicit`, so the report
can say "this run used a depth floor of 8, not the published 25". A silently
changed threshold and a threshold with no source are the same problem.

**The model is a flag, never a gate.** ``--llm-host`` points at any local
server; ``--no-llm`` produces a rule-only report. Neither changes a single
verdict — the verdicts are computed before the model is asked anything — and the
report states which of the two happened.

**Nothing here decides anything clinical.** The CLI collects inputs, builds a
:class:`~mjolnir.config.Config`, and hands both to :mod:`mjolnir.pipeline`. When
it prints a number it prints one a record already carried.
"""

from __future__ import annotations

import argparse
import os
import shlex
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from . import __version__, config, doctor
from .config import Config, PROFILES, default_db_dir
from .db import fetch as db_fetch
from .db import registry as db_registry
from .pipeline import (Pipeline, RunOptions, load_cohort_result, load_sample_result,
                       write_outputs)
from .records import PLATFORMS
from .seqio import describe_inputs, detect_inputs
from .utils import LOG, MjolnirError, cpu_count, ensure_dir, setup_logging, tempdir

EPILOG = """
examples:
  # one Illumina isolate, everything on, PDF + HTML + tables
  mjolnir run -1 30-20_S1_R1_001.fastq.gz -2 30-20_S1_R2_001.fastq.gz -o results/

  # ONT reads, forcing the platform rather than detecting it
  mjolnir run --reads barcode07.fastq.gz --platform ont -o results/

  # an assembly: no allele fractions, and the report says so
  mjolnir run -a isolate.fasta -o results/

  # a directory of paired FASTQ as one cohort, clustered at 6 SNPs
  mjolnir cohort reads/ -o outbreak/ --distance 6

  # render a finished run again, research profile
  mjolnir report results/30-20.json -o results/research/ --profile research

  # what is installed, what it costs, and under which licence
  mjolnir db list
  mjolnir db info who-catalogue-v2
  mjolnir db fetch

  # the environment, required versus optional
  mjolnir doctor

Mjolnir v{version}
""".format(version=__version__)


# ---------------------------------------------------------------------------
# Argument groups
# ---------------------------------------------------------------------------

def _add_inputs(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("inputs")
    group.add_argument("inputs", nargs="*", type=Path, metavar="INPUT",
                       help="FASTQ files, assemblies, or directories to scan")
    group.add_argument("-1", "--r1", type=Path, default=None, metavar="FASTQ",
                       help="forward reads (with --r2)")
    group.add_argument("-2", "--r2", type=Path, default=None, metavar="FASTQ",
                       help="reverse reads")
    group.add_argument("--reads", type=Path, default=None, metavar="FASTQ",
                       help="single-end or ONT reads")
    group.add_argument("-a", "--assembly", type=Path, default=None, metavar="FASTA",
                       help="an assembled genome")
    group.add_argument("--sample", default=None, metavar="NAME",
                       help="sample name (default: derived from the filename)")
    group.add_argument("-r", "--recursive", action="store_true",
                       help="descend into subdirectories when scanning a directory")
    group.add_argument("--platform", default=None, metavar="NAME",
                       choices=sorted(set(PLATFORMS)),
                       help="force the platform instead of detecting it: {0}".format(
                           ", ".join(PLATFORMS)))
    group.add_argument("--ref", "--reference", dest="reference", type=Path, default=None,
                       metavar="FASTA",
                       help="reference genome (default: H37Rv for MTBC, the ANI "
                            "set's own genome for anything else)")


def _add_thresholds(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group(
        "thresholds (every default is a published number; an override is printed "
        "in the report)")
    group.add_argument("--min-depth", type=int, default=None, metavar="X",
                       help="target mean depth (default: {0}; {1})".format(
                           config.MIN_DEPTH, config.source_for("min_depth")))
    group.add_argument("--degraded-depth-floor", type=int, default=None, metavar="X",
                       help="depth below which a sample is not callable "
                            "(default: {0})".format(config.DEGRADED_DEPTH_FLOOR))
    group.add_argument("--min-breadth", type=float, default=None, metavar="FRAC",
                       help="fraction of the reference that must reach the floor "
                            "(default: {0})".format(config.MIN_BREADTH))
    group.add_argument("--min-reads-illumina", type=int, default=None, metavar="N",
                       help="reads supporting a variant on Illumina (default: {0})".format(
                           config.MIN_READS_ILLUMINA))
    group.add_argument("--min-reads-ont", type=int, default=None, metavar="N",
                       help="reads supporting a variant on ONT (default: {0})".format(
                           config.MIN_READS_ONT))
    group.add_argument("--major-variant-fraction", type=float, default=None,
                       metavar="FRAC",
                       help="read support at or above which a variant is major "
                            "(default: {0})".format(config.MAJOR_VARIANT_FRACTION))
    group.add_argument("--min-minor-variant-fraction", type=float, default=None,
                       metavar="FRAC",
                       help="lowest minor-allele fraction reported on Illumina "
                            "(default: {0})".format(config.MIN_MINOR_VARIANT_FRACTION))
    group.add_argument("--distance", type=int, default=None, metavar="SNPS",
                       help="cohort clustering distance (default: {0}; the TB "
                            "conventions are {1} and {2}, and the local M. chimaera "
                            "baseline was {3})".format(
                                config.DEFAULT_CLUSTER_DISTANCE,
                                config.CLUSTER_SNP_STRICT, config.CLUSTER_SNP_RELAXED,
                                config.CHIMAERA_LOCAL_DISTANCE))
    group.add_argument("--mask", type=Path, default=None, metavar="BED",
                       help="repeat/low-complexity mask for cohort distances "
                            "(default: tbdb mask.bed from the database directory)")
    group.add_argument("--mtbseq-compat", action="store_true",
                       help="reproduce MTBseq's denominator, tie-break and filter "
                            "stack so a run can be reconciled with a legacy one")


def _add_analysis(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("analysis")
    group.add_argument("--mapper", default=None, metavar="TOOL",
                       help="force the aligner (default: bwa-mem2 or bwa on "
                            "Illumina, minimap2 on ONT)")
    group.add_argument("--caller", default=None, metavar="TOOL",
                       help="force the variant caller (default: bcftools on "
                            "Illumina, Clair3 on ONT)")
    group.add_argument("--clair3-model", default=None, metavar="NAME",
                       help="Clair3 model directory name under the database root")
    group.add_argument("--allow-degraded-ont-calling", action="store_true",
                       help="allow the bcftools fallback on ONT when Clair3 is "
                            "absent. It disables indel calling entirely and the "
                            "report says the numbers came from it")
    group.add_argument("--gff", metavar="GFF3",
                       help="gene models used to name variants (default: beside "
                            "the reference, or tbdb's H37Rv models for MTBC)")
    group.add_argument("--build-index", action="store_true",
                       help="build a missing reference index instead of refusing")
    group.add_argument("--no-markdup", dest="mark_duplicates", action="store_false",
                       default=None, help="do not mark duplicates on Illumina")
    group.add_argument("--kraken2-db", type=Path, default=None, metavar="DIR",
                       help="Kraken2 index for the read composition screen. Only a "
                            "mycobacterial pangenome index is treated as "
                            "informative (measured sensitivity with a standard "
                            "index is {0} for M. tuberculosis reads)".format(
                                config.KRAKEN2_MTB_SENSITIVITY_STANDARD))
    group.add_argument("--kraken2-confidence", type=float, default=None, metavar="FRAC",
                       help="Kraken2 --confidence (default: {0}; Kraken2's own 0.0 "
                            "default is refused)".format(config.KRAKEN2_MIN_CONFIDENCE))
    group.add_argument("--kraken2-report", type=Path, default=None, metavar="FILE",
                       help="use a Kraken2 report that already exists instead of "
                            "running kraken2")
    group.add_argument("--no-typing", dest="typing", action="store_false", default=True,
                       help="skip species and lineage typing")
    group.add_argument("--no-resistance", dest="resistance", action="store_false",
                       default=True, help="skip resistance calling")
    group.add_argument("--no-contamination", dest="contamination", action="store_false",
                       default=True, help="skip the contamination panel")


def _add_llm(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("interpretation")
    group.add_argument("--llm-host", default=None, metavar="URL",
                       help="local model host: ollama, vLLM, SGLang or llama.cpp "
                            "(default: ${0} or {1})".format(
                                config.LLM_HOST_ENV_VAR, "http://127.0.0.1:11434"))
    group.add_argument("--llm-model", default=None, metavar="NAME",
                       help="model tag (default: the first the host advertises)")
    group.add_argument("--no-llm", dest="use_llm", action="store_false", default=True,
                       help="do not contact a model; the report is rule-derived "
                            "and says so")


def _add_outputs(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("output")
    group.add_argument("-o", "--outdir", type=Path, default=None, metavar="DIR",
                       help="write results here (created if missing)")
    group.add_argument("-f", "--format", dest="formats", action="append", default=None,
                       metavar="FMT",
                       help="output formats, repeatable or comma-separated: "
                            "tsv, json, html, pdf (default: all; tsv and json are "
                            "written together)")
    group.add_argument("--profile", default=None, choices=list(PROFILES),
                       help="report profile: {0} (default: clinical)".format(
                           ", ".join(PROFILES)))


def _add_common(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("runtime")
    group.add_argument("-t", "--threads", type=int, default=None,
                       help="CPU threads (default: all {0})".format(cpu_count()))
    group.add_argument("--db-dir", type=Path, default=None, metavar="DIR",
                       help="database directory (default: ${0} or {1})".format(
                           config.DB_ENV_VAR, default_db_dir()))
    group.add_argument("--tmpdir", type=Path, default=None, metavar="DIR",
                       help="directory for intermediate files (default: system temp)")
    group.add_argument("--keep-temp", action="store_true",
                       help="keep BAMs, VCFs and logs for debugging")
    group.add_argument("-v", "--verbose", action="count", default=0,
                       help="verbose logging (repeat for more)")
    group.add_argument("-q", "--quiet", action="store_true",
                       help="only warnings and errors")


# ---------------------------------------------------------------------------
# Turning arguments into a Config
# ---------------------------------------------------------------------------

#: Command-line flag to ``Config`` field. Every one of these overrides a
#: published number, so they all go through ``set_explicit`` and all appear in
#: the report's list of what this run changed.
_THRESHOLD_FLAGS: Dict[str, str] = {
    "min_depth": "min_depth",
    "degraded_depth_floor": "degraded_depth_floor",
    "min_breadth": "min_breadth",
    "min_reads_illumina": "min_reads_illumina",
    "min_reads_ont": "min_reads_ont",
    "major_variant_fraction": "major_variant_fraction",
    "min_minor_variant_fraction": "min_minor_variant_fraction",
    "distance": "cluster_distance",
    "kraken2_confidence": "kraken2_confidence",
}


def build_config(args: Any) -> Config:
    """A validated :class:`~mjolnir.config.Config` from the parsed arguments."""
    cfg = Config(
        db_dir=args.db_dir or default_db_dir(),
        out_dir=args.outdir or Path("mjolnir_out"),
        threads=args.threads if args.threads else cpu_count(),
        platform=getattr(args, "platform", None),
        reference=getattr(args, "reference", None),
        profile=args.profile or "clinical",
        mask_bed=getattr(args, "mask", None),
        kraken2_db=getattr(args, "kraken2_db", None),
        llm_host=getattr(args, "llm_host", None) or "",
        llm_model=getattr(args, "llm_model", None) or "",
        use_llm=bool(getattr(args, "use_llm", True)),
        keep_temp=bool(getattr(args, "keep_temp", False)),
        tmp_dir=getattr(args, "tmpdir", None),
        mtbseq_compat=bool(getattr(args, "mtbseq_compat", False)),
    )
    if args.threads is not None:
        cfg.set_explicit("threads", cfg.threads)
    for flag, field_name in _THRESHOLD_FLAGS.items():
        value = getattr(args, flag, None)
        if value is not None:
            cfg.set_explicit(field_name, value)
    cfg.validate()
    if cfg.tmp_dir is not None:
        ensure_dir(cfg.tmp_dir)
    return cfg


def build_options(args: Any, *, callable_regions: bool = False) -> RunOptions:
    return RunOptions(
        mapper=getattr(args, "mapper", None) or "",
        caller=getattr(args, "caller", None) or "",
        clair3_model=getattr(args, "clair3_model", None) or "",
        build_index=bool(getattr(args, "build_index", False)),
        gff=getattr(args, "gff", "") or "",
        mark_duplicates=getattr(args, "mark_duplicates", None),
        allow_degraded_ont_calling=bool(
            getattr(args, "allow_degraded_ont_calling", False)),
        kraken2_report=getattr(args, "kraken2_report", None),
        callable_regions=callable_regions,
        typing=bool(getattr(args, "typing", True)),
        resistance=bool(getattr(args, "resistance", True)),
        contamination=bool(getattr(args, "contamination", True)),
        interpret=bool(getattr(args, "use_llm", True)),
    )


def collect_samples(args: Any) -> List[Any]:
    """Detect every sample the flags describe, or say why nothing was found.

    Detection is :mod:`mjolnir.seqio`'s job — pairing, platform evidence and the
    checks that come with it — so this function only assembles the path list and
    hands the naming decisions over.
    """
    paths: List[Path] = list(args.inputs or [])
    explicit: List[Path] = []
    if args.r1 is not None:
        explicit.append(args.r1)
    if args.r2 is not None:
        if args.r1 is None:
            raise MjolnirError("--r2 was given without --r1")
        explicit.append(args.r2)
    if args.reads is not None:
        explicit.append(args.reads)
    if args.assembly is not None:
        explicit.append(args.assembly)
    paths.extend(explicit)
    if not paths:
        raise MjolnirError(
            "no inputs given. Pass FASTQ files, an assembly or a directory; see "
            "'mjolnir run --help'")

    detected = detect_inputs(
        paths, platform=args.platform, sample_id=args.sample,
        reference=getattr(args, "reference", None),
        recursive=bool(getattr(args, "recursive", False)))
    if not detected:
        raise MjolnirError(
            "no readable FASTQ or FASTA input was found in: {0}".format(
                ", ".join(str(p) for p in paths)))
    return detected


def _command_line() -> str:
    return "mjolnir " + " ".join(shlex.quote(a) for a in sys.argv[1:])


def _require_outdir(args: Any) -> Path:
    if args.outdir is None:
        raise MjolnirError(
            "no output directory set; pass -o/--outdir. Mjolnir writes a PDF, an "
            "HTML report and the TSV/JSON artefacts, so it needs somewhere to "
            "put them")
    return ensure_dir(args.outdir)


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_run(args: Any) -> int:
    """One sample (or several, independently) through the whole pipeline."""
    out_dir = _require_outdir(args)
    cfg = build_config(args)
    detected = collect_samples(args)
    if len(detected) > 1:
        LOG.info("%d samples detected; each is analysed independently. Use "
                 "'mjolnir cohort' to compare them to one another.", len(detected))
    LOG.info("Mjolnir v%s | %d sample(s) | %d threads | db %s",
             __version__, len(detected), cfg.threads, cfg.db_dir)
    LOG.debug("%s", _command_line())
    print(describe_inputs(detected), file=sys.stderr)

    pipeline = Pipeline(cfg, build_options(args))
    evidence = dict((d.sample_id, d.checks) for d in detected)
    with tempdir(prefix="mjolnir.", keep=cfg.keep_temp, parent=cfg.tmp_dir) as workdir:
        results = pipeline.run([d.sample for d in detected], workdir, evidence)
        if cfg.keep_temp:
            LOG.info("intermediate files kept in %s", workdir)
        written = write_outputs(results, out_dir, formats=args.formats,
                                profile=cfg.profile)
    _recap(results, None, written)
    return 0 if all(r.status != "fail" for r in results) else 1


def cmd_cohort(args: Any) -> int:
    """Several samples, then the joint table, distances and clusters."""
    out_dir = _require_outdir(args)
    cfg = build_config(args)
    detected = collect_samples(args)
    LOG.info("Mjolnir v%s | cohort of %d | %d threads | clustering at %d SNPs",
             __version__, len(detected), cfg.threads, cfg.cluster_distance)
    print(describe_inputs(detected), file=sys.stderr)

    pipeline = Pipeline(cfg, build_options(args, callable_regions=True))
    evidence = dict((d.sample_id, d.checks) for d in detected)
    cohort = None
    with tempdir(prefix="mjolnir.", keep=cfg.keep_temp, parent=cfg.tmp_dir) as workdir:
        results = pipeline.run([d.sample for d in detected], workdir, evidence)
        try:
            cohort = pipeline.build_cohort(results)
        except MjolnirError as exc:
            # The per-sample results are finished and are worth writing even when
            # the comparison between them is not possible.
            LOG.error("no cohort comparison: %s", exc)
            for result in results:
                result.warnings.append("cohort comparison not made: {0}".format(exc))
        written = write_outputs(results, out_dir, cohort=cohort,
                                formats=args.formats, profile=cfg.profile)
    _recap(results, cohort, written)
    return 0 if cohort is not None else 1


def cmd_report(args: Any) -> int:
    """Render a finished run again, without re-deriving anything.

    Regenerating a report must not be able to change a call, so this reads the
    JSON artefacts and renders them. A different profile, a newly installed
    reportlab or a fixed template are the reasons to run it; a different answer
    is not one of them.
    """
    out_dir = _require_outdir(args)
    results = []
    cohort = None
    for path in args.results:
        results.append(load_sample_result(Path(path)))
    if args.cohort is not None:
        cohort = load_cohort_result(Path(args.cohort))
    if not results and cohort is None:
        raise MjolnirError("nothing to render: give at least one <sample>.json")
    profile = args.profile or (results[0].profile if results else "clinical")
    written = write_outputs(results, out_dir, cohort=cohort, formats=args.formats,
                            profile=profile)
    for path in written:
        print(path)
    return 0


def cmd_db(args: Any) -> int:
    """List, fetch and describe the reference databases, licences included."""
    db_root = args.db_dir or default_db_dir()
    action = args.db_action

    if action == "list":
        print(db_fetch.format_listing(db_root, verbose=bool(args.verbose)))
        print("\nlicences and attribution:")
        for line in db_registry.attributions():
            print("  " + line)
        return 0

    if action == "info":
        names = db_registry.resolve_names(args.names) if args.names else sorted(
            db_registry.DATABASES)
        for name in names:
            spec = db_registry.spec_for(name)
            print(spec.describe())
            try:
                installed = db_fetch.database_version(name, db_root=db_root)
            except MjolnirError:
                installed = None
            if installed is not None:
                print("  installed   {0}".format(installed.path))
                print("  version     {0}".format(installed.version))
                print("  checksum    {0}".format(installed.checksum or "-"))
                print("  fetched     {0}".format(installed.fetched or "-"))
            else:
                print("  installed   no ({0})".format(db_registry.fetch_hint(name)))
            print()
        return 0

    if action == "fetch":
        names = db_registry.resolve_names(args.names) if args.names else None
        if args.dry_run:
            print(db_fetch.format_plan(names, db_root=db_root, force=args.force))
            return 0
        try:
            versions = db_fetch.fetch_databases(names, db_root=db_root,
                                                force=args.force, strict=args.strict)
        except MjolnirError:
            raise
        for version in versions:
            print("ok {0:<28} {1}".format(version.name, version.version))
            if version.licence:
                print("   licence   {0}".format(version.licence))
            if version.citation:
                print("   cite      {0}".format(version.citation))
        if not versions:
            print("nothing fetched; everything requested is already installed "
                  "(--force to replace it)")
        return 0

    if action == "verify":
        names = db_registry.resolve_names(args.names) if args.names else [
            v.name for v in db_fetch.installed(db_root)]
        if not names:
            print(db_registry.NO_DATABASES_TEXT)
            return 1
        problems = 0
        for name in names:
            issues = db_fetch.verify_installed(name, db_root=db_root, deep=args.deep)
            if issues:
                problems += len(issues)
                for issue in issues:
                    print("PROBLEM {0}: {1}".format(name, issue))
            else:
                print("ok      {0}".format(name))
        return 1 if problems else 0

    raise MjolnirError("unknown db action {0!r}".format(action))


def cmd_doctor(args: Any) -> int:
    """The whole environment: tools, Python deps, databases, capabilities."""
    cfg = Config(db_dir=args.db_dir or default_db_dir(),
                 kraken2_db=getattr(args, "kraken2_db", None),
                 llm_host=getattr(args, "llm_host", None) or "",
                 use_llm=bool(getattr(args, "use_llm", True)))
    diagnosis = doctor.diagnose(cfg)
    if args.json:
        import json

        print(json.dumps(diagnosis.to_dict(), indent=2, sort_keys=True))
    else:
        print(doctor.render(diagnosis, verbose=bool(args.verbose)))
    return 0 if diagnosis.ok else 1


# ---------------------------------------------------------------------------
# Recap
# ---------------------------------------------------------------------------

def _recap(results: Sequence[Any], cohort: Optional[Any],
           written: Sequence[Path]) -> None:
    """One screen of what happened, on stderr so a piped table stays clean."""
    out = sys.stderr
    print("", file=out)
    for result in results:
        resistant = [d.drug for d in result.resistant_drugs()]
        print("{0:<20} {1:<34} {2}".format(
            result.sample_id[:20], result.species.display[:34],
            result.lineage.display), file=out)
        print("{0:<20} validity: {1:<14} {2}".format(
            "", result.contamination.verdict,
            "resistance determinants: " + (", ".join(resistant) if resistant
                                           else "none detected")), file=out)
        for warning in result.warnings[:5]:
            print("{0:<20} ! {1}".format("", warning), file=out)
        unmeasured = result.unmeasured()
        if unmeasured:
            print("{0:<20} not measured: {1}".format(
                "", ", ".join(unmeasured[:8])
                + (" (+{0} more)".format(len(unmeasured) - 8)
                   if len(unmeasured) > 8 else "")), file=out)
    if cohort is not None:
        print("\ncohort: {0} samples, {1} clusters at {2} SNPs ({3})".format(
            len(cohort.samples), len(cohort.clusters), cohort.threshold,
            cohort.mask_name or "no mask"), file=out)
        for cluster in cohort.clusters:
            print("  {0}: {1}".format(cluster.cluster_id,
                                      ", ".join(cluster.members)), file=out)
    print("\n{0} file(s) written".format(len(written)), file=out)


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mjolnir",
        description="Mjolnir - resistance, lineage, species and contamination "
                    "calling for the M. tuberculosis complex and non-tuberculous "
                    "mycobacteria, from Illumina reads, ONT reads or assemblies.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version",
                        version="mjolnir {0}".format(__version__))
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    run_parser = subparsers.add_parser(
        "run", help="one sample: full outputs and a PDF",
        description="Analyse one sample end to end: species, lineage, resistance, "
                    "contamination and a clinician-first report.",
        epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    _add_inputs(run_parser)
    _add_analysis(run_parser)
    _add_thresholds(run_parser)
    _add_llm(run_parser)
    _add_outputs(run_parser)
    _add_common(run_parser)
    run_parser.set_defaults(func=cmd_run)

    cohort_parser = subparsers.add_parser(
        "cohort", help="many samples: joint analysis, distances, clusters",
        description="Analyse several samples, then compare them: a joint variant "
                    "table, masked pairwise SNP distances with their shared "
                    "callable-site denominators, and clusters at a stated "
                    "threshold.",
        epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    _add_inputs(cohort_parser)
    _add_analysis(cohort_parser)
    _add_thresholds(cohort_parser)
    _add_llm(cohort_parser)
    _add_outputs(cohort_parser)
    _add_common(cohort_parser)
    cohort_parser.set_defaults(func=cmd_cohort)

    report_parser = subparsers.add_parser(
        "report", help="regenerate the report from a finished run",
        description="Render finished results again. Reads the JSON artefacts and "
                    "renders them; it re-derives nothing, so a regenerated report "
                    "cannot disagree with the run that produced it.")
    report_parser.add_argument("results", nargs="*", type=Path, metavar="JSON",
                               help="<sample>.json files written by 'mjolnir run'")
    report_parser.add_argument("--cohort", type=Path, default=None, metavar="JSON",
                               help="cohort.json written by 'mjolnir cohort'")
    _add_outputs(report_parser)
    report_parser.add_argument("-v", "--verbose", action="count", default=0)
    report_parser.add_argument("-q", "--quiet", action="store_true")
    report_parser.set_defaults(func=cmd_report)

    db_parser = subparsers.add_parser(
        "db", help="fetch and inspect the reference databases",
        description="Mjolnir keeps its databases in ${0} (default {1}). Every one "
                    "of them records its version, checksum, licence and citation, "
                    "and the report prints them.".format(
                        config.DB_ENV_VAR, default_db_dir()))
    db_sub = db_parser.add_subparsers(dest="db_action", metavar="ACTION")
    for action, help_text in (
        ("list", "show what is installed, with licences"),
        ("fetch", "download and verify databases"),
        ("info", "describe one or more databases in full"),
        ("verify", "re-check installed files against their recorded checksums"),
    ):
        sub = db_sub.add_parser(action, help=help_text)
        sub.add_argument("names", nargs="*", metavar="NAME",
                         help="database or group names (default: {0})".format(
                             ", ".join(sorted(db_registry.DB_GROUPS))))
        if action == "fetch":
            sub.add_argument("--force", action="store_true",
                             help="re-fetch databases that are already installed")
            sub.add_argument("--dry-run", action="store_true",
                             help="print what would be fetched, and its size")
            sub.add_argument("--strict", action="store_true",
                             help="treat a checksum mismatch as fatal rather than "
                                  "as a recorded warning")
        if action == "verify":
            sub.add_argument("--deep", action="store_true",
                             help="re-checksum every file rather than checking "
                                  "sizes and presence")
        sub.add_argument("--db-dir", type=Path, default=None, metavar="DIR",
                         help="database directory (default: {0})".format(default_db_dir()))
        sub.add_argument("-v", "--verbose", action="count", default=0)
        sub.add_argument("-q", "--quiet", action="store_true")
        sub.set_defaults(func=cmd_db, db_action=action)
    db_parser.set_defaults(func=cmd_db, db_action=None)

    doctor_parser = subparsers.add_parser(
        "doctor", help="the environment: tools, databases, capabilities",
        description="Report the whole environment before running anything: which "
                    "tools and databases are present, which are missing, and "
                    "exactly which capability each absence removes.")
    doctor_parser.add_argument("--db-dir", type=Path, default=None, metavar="DIR")
    doctor_parser.add_argument("--kraken2-db", type=Path, default=None, metavar="DIR",
                               help="Kraken2 index to interrogate")
    doctor_parser.add_argument("--llm-host", default=None, metavar="URL",
                               help="model host to probe")
    doctor_parser.add_argument("--no-llm", dest="use_llm", action="store_false",
                               default=True, help="do not probe a model host")
    doctor_parser.add_argument("--json", action="store_true",
                               help="machine-readable output")
    doctor_parser.add_argument("-v", "--verbose", action="count", default=0)
    doctor_parser.add_argument("-q", "--quiet", action="store_true")
    doctor_parser.set_defaults(func=cmd_doctor)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(getattr(args, "verbose", 0), getattr(args, "quiet", False))

    if not getattr(args, "command", None):
        parser.print_help()
        return 1
    if args.command == "db" and not getattr(args, "db_action", None):
        print("usage: mjolnir db {list,fetch,info,verify}\n\n"
              "Start with:  mjolnir db fetch", file=sys.stderr)
        return 1

    try:
        return int(args.func(args))
    except MjolnirError as exc:
        # These are the errors that name what to install or fetch. A traceback
        # would bury that line in twenty frames of no interest to the reader.
        LOG.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        LOG.error("interrupted")
        return 130
    except BrokenPipeError:
        os._exit(0)


if __name__ == "__main__":
    sys.exit(main())
