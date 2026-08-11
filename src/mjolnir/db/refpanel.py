"""The mycobacterial reference panel: genomes, gene models, and the manifest.

Species identification is ANI against a set of reference genomes, and design §14
left the composition of that set open. This module answers it, and does the part
that turned out to matter more than the choice of species: **it fetches the gene
models too**.

Without a GFF a variant is a coordinate and nothing else. The NTM resistance
rules — ``erm(41)`` sequevar typing, ``rrl`` 2058/2059, ``rrs`` 1408 — are all
keyed on gene names, so on a reference with no annotation they are written,
tested, and completely dead. That was the state of every NTM reference Mjolnir
had until this existed.

The panel is deliberately wider than the organisms Mjolnir reports on. A species
that is *not* in the set cannot be excluded by ANI, so the laboratory
contaminants and the near neighbours are there to be matched and named rather
than to be silently absorbed into the closest thing that happens to be present.

Fetched from the NCBI datasets API rather than declared in ``registry.py``
because that API serves a zip per accession, not a URL per file, and the
alternative is hard-coding an assembly-name path segment that changes whenever
NCBI re-annotates.
"""

from __future__ import annotations

import io
import json
import time
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ..utils import LOG, MjolnirError, PathLike, ensure_dir

DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

#: Where the panel lives, relative to the database root.
PANEL_DIRNAME = "ani"
MANIFEST_NAME = "references.tsv"


@dataclass(frozen=True)
class PanelSpecies:
    """One species Mjolnir wants a reference genome for."""

    taxid: int
    name: str
    complex: str = ""
    why: str = ""


#: The panel. MTBC members are here so that ANI can place a sample in the
#: complex; it cannot separate them, and :mod:`mjolnir.typing.species` refuses to
#: try. Everything else is a species a mycobacterial isolate might actually be,
#: plus the laboratory contaminants that must be nameable rather than absorbed.
PANEL: Tuple[PanelSpecies, ...] = (
    PanelSpecies(1773, "Mycobacterium tuberculosis", "MTBC",
                 "the MTBC anchor; members are separated by barcode, not by ANI"),
    PanelSpecies(1765, "Mycobacterium tuberculosis variant bovis", "MTBC",
                 "animal lineage La1, and BCG's parent"),
    PanelSpecies(33894, "Mycobacterium tuberculosis variant africanum", "MTBC",
                 "lineages 5 and 6"),
    PanelSpecies(1806, "Mycobacterium tuberculosis variant microti", "MTBC"),
    PanelSpecies(78331, "Mycobacterium canettii", "MTBC",
                 "the most divergent MTBC member; a useful outgroup"),
    PanelSpecies(222805, "Mycobacterium intracellulare subsp. chimaera", "MAC",
                 "the heater-cooler outbreak organism"),
    PanelSpecies(1767, "Mycobacterium intracellulare", "MAC",
                 "chimaera's nearest neighbour; the hard MAC call"),
    PanelSpecies(1764, "Mycobacterium avium", "MAC"),
    PanelSpecies(1782, "Mycobacterium scrofulaceum", "MAC"),
    PanelSpecies(1305738, "Mycobacterium paraintracellulare", "MAC"),
    PanelSpecies(1157943, "Mycobacterium yongonense", "MAC"),
    PanelSpecies(36809, "Mycobacteroides abscessus", "abscessus",
                 "erm(41) sequevar typing applies to this group"),
    PanelSpecies(1698, "Mycobacteroides chelonae", "abscessus"),
    PanelSpecies(1768, "Mycobacterium kansasii"),
    PanelSpecies(1781, "Mycobacterium marinum"),
    PanelSpecies(1769, "Mycolicibacterium fortuitum"),
    PanelSpecies(1778, "Mycobacterium gordonae",
                 "the commonest contaminant of mycobacterial culture"),
    PanelSpecies(1789, "Mycobacterium xenopi"),
    PanelSpecies(1780, "Mycobacterium malmoense"),
    PanelSpecies(1784, "Mycobacterium simiae"),
    PanelSpecies(1787, "Mycobacterium szulgai"),
    PanelSpecies(1809, "Mycobacterium ulcerans"),
    PanelSpecies(29311, "Mycobacterium haemophilum"),
    PanelSpecies(56689, "Mycolicibacterium mucogenicum"),
    PanelSpecies(1772, "Mycolicibacterium smegmatis",
                 why="a laboratory strain, and a contaminant worth naming"),
    PanelSpecies(1096, "Mycobacterium leprae"),
    PanelSpecies(1771, "Mycolicibacterium phlei"),
)


def _get(url: str, timeout: int = 120) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "mjolnir"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as handle:
            return handle.read()
    except urllib.error.URLError as exc:
        raise MjolnirError("could not reach the NCBI datasets API ({0}): {1}\n"
                           "  the reference panel needs network access once; "
                           "after that it is a local directory".format(url, exc))


def best_accession(taxid: int) -> Tuple[Optional[str], str]:
    """The assembly to use for a taxon: RefSeq and complete, where one exists.

    Preferring a complete genome is not fussiness. ANI over a draft assembly is
    computed across contig boundaries that are artefacts of assembly rather than
    of biology, and the margin that separates *M. chimaera* from
    *M. intracellulare* is two percent.
    """
    for level in ("&filters.assembly_level=complete_genome", ""):
        url = "{0}/genome/taxon/{1}/dataset_report?page_size=20{2}".format(
            DATASETS_API, taxid, level)
        try:
            reports = (json.loads(_get(url)) or {}).get("reports") or []
        except (ValueError, MjolnirError):
            continue
        ranked = sorted(reports, key=lambda r: (
            not str(r.get("accession", "")).startswith("GCF"),
            (r.get("assembly_info") or {}).get("refseq_category") != "reference genome"))
        for report in ranked:
            accession = report.get("accession")
            if accession:
                return accession, (report.get("organism") or {}).get("organism_name", "")
    return None, ""


def fetch_assembly(accession: str) -> Tuple[Optional[bytes], Optional[bytes]]:
    """``(genome FASTA, GFF)`` for an accession; either may be None."""
    url = ("{0}/genome/accession/{1}/download?include_annotation_type=GENOME_FASTA"
           "&include_annotation_type=GENOME_GFF").format(DATASETS_API, accession)
    archive = zipfile.ZipFile(io.BytesIO(_get(url, timeout=300)))
    fasta = gff = None
    for name in archive.namelist():
        if name.endswith(".fna") and fasta is None:
            fasta = archive.read(name)
        elif name.endswith(".gff") and gff is None:
            gff = archive.read(name)
    return fasta, gff


def panel_dir(db_root: PathLike) -> Path:
    return Path(str(db_root)).expanduser() / PANEL_DIRNAME


def _same_organism(requested: str, reported: str) -> bool:
    """Whether NCBI's organism name is the species that was asked for.

    Compared on the first two words, so a subspecies or a strain suffix is
    accepted - "Mycobacterium avium subsp. hominissuis" answers a request for
    "Mycobacterium avium" - while a different species is not.
    """
    def head(text):
        parts = str(text or "").replace("[", "").replace("]", "").split()
        return " ".join(parts[:2]).lower()
    return bool(reported) and head(requested) == head(reported)


def build_panel(db_root: PathLike, *, species: Sequence[PanelSpecies] = PANEL,
                overwrite: bool = False, pause: float = 0.4) -> Dict[str, int]:
    """Fetch the panel and write its manifest. Returns a small tally.

    A species that cannot be fetched is skipped and counted rather than raising:
    a panel of 26 is worth having while one species is unavailable, and the
    manifest records exactly which genomes are behind the ANI call.
    """
    destination = ensure_dir(panel_dir(db_root))
    rows: List[Tuple[str, str, str, str, bool]] = []
    tally = {"fetched": 0, "with_gene_models": 0, "skipped": 0, "reused": 0,
             "mislabelled": 0}

    for entry in species:
        accession, organism = best_accession(entry.taxid)
        if not accession:
            LOG.warning("no assembly found for %s (taxid %d)", entry.name, entry.taxid)
            tally["skipped"] += 1
            continue
        # The taxid is a number a human typed, and a wrong one returns a real
        # genome for the wrong organism. Taxid 1770 is M. avium subsp.
        # paratuberculosis, not M. haemophilum, and the panel carried that
        # genome under the wrong name until a sample matched it at 98.5% and the
        # number made no biological sense. What NCBI says the genome is wins.
        if not _same_organism(entry.name, organism):
            LOG.warning(
                "taxid %d returned %r, which is not %r; skipping rather than "
                "adding a genome under a name that is not its own",
                entry.taxid, organism, entry.name)
            tally["mislabelled"] = tally.get("mislabelled", 0) + 1
            tally["skipped"] += 1
            continue
        fasta_path = destination / "{0}.fna".format(accession)
        gff_path = destination / "{0}.gff".format(accession)
        if fasta_path.exists() and not overwrite:
            rows.append((fasta_path.name, organism or entry.name, entry.complex,
                         accession, gff_path.exists()))
            tally["reused"] += 1
            if gff_path.exists():
                tally["with_gene_models"] += 1
            continue
        try:
            fasta, gff = fetch_assembly(accession)
        except (MjolnirError, zipfile.BadZipFile) as exc:
            LOG.warning("could not fetch %s for %s: %s", accession, entry.name, exc)
            tally["skipped"] += 1
            continue
        if not fasta:
            LOG.warning("%s returned no genome sequence for %s", accession, entry.name)
            tally["skipped"] += 1
            continue
        fasta_path.write_bytes(fasta)
        if gff:
            gff_path.write_bytes(gff)
            tally["with_gene_models"] += 1
        # NCBI's name, not the one that was asked for: they agree by the check
        # above, and where they differ in detail the authority is the record.
        rows.append((fasta_path.name, organism or entry.name, entry.complex,
                     accession, bool(gff)))
        tally["fetched"] += 1
        LOG.info("panel: %s -> %s%s", entry.name, accession,
                 "" if gff else " (NO gene models)")
        time.sleep(pause)

    if not rows:
        raise MjolnirError(
            "no reference genome could be fetched, so the panel would name "
            "nothing. Check network access to {0}".format(DATASETS_API))
    write_manifest(destination / MANIFEST_NAME, rows)
    return tally


def write_manifest(path: PathLike, rows: Sequence[Tuple[str, str, str, str, bool]]) -> Path:
    """The manifest ``typing/species.py`` reads, plus what each genome can do.

    The ``note`` column is load-bearing: a reference with no gene models can
    still carry an ANI species call and cannot carry a resistance call, and a
    reader comparing two runs needs to see which they had.
    """
    target = Path(path)
    ensure_dir(target.parent)
    with open(str(target), "w") as handle:
        handle.write("file\tname\tcomplex\taccession\tsubspecies\tnote\tsource\n")
        for filename, name, complex_name, accession, has_gff in rows:
            note = ("gene models present" if has_gff else
                    "no gene models: variants against this reference cannot be "
                    "named, so no resistance rule can fire")
            handle.write("{0}\t{1}\t{2}\t{3}\t\t{4}\tNCBI RefSeq\n".format(
                filename, name, complex_name, accession, note))
    LOG.info("wrote %s: %d references, %d with gene models",
             target, len(rows), sum(1 for r in rows if r[4]))
    return target


def panel_status(db_root: PathLike) -> Dict[str, int]:
    """What is installed: genomes, gene models, and species without models."""
    directory = panel_dir(db_root)
    manifest = directory / MANIFEST_NAME
    if not manifest.exists():
        return {"references": 0, "with_gene_models": 0}
    references = models = 0
    with open(str(manifest), "rt", errors="replace") as handle:
        next(handle, None)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0]:
                continue
            references += 1
            if (directory / fields[0]).with_suffix(".gff").exists():
                models += 1
    return {"references": references, "with_gene_models": models}


__all__ = [
    "MANIFEST_NAME", "PANEL", "PANEL_DIRNAME", "PanelSpecies", "best_accession",
    "build_panel", "fetch_assembly", "panel_dir", "panel_status", "write_manifest",
]
