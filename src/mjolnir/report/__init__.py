"""The deliverable: a clinician-first PDF, the same report as HTML, and the tables.

Three modules, one data layer. :mod:`mjolnir.report.tables` turns a
:class:`~mjolnir.records.SampleResult` into rows and panels;
:mod:`mjolnir.report.pdf` and :mod:`mjolnir.report.html` render those same rows
and the same vector scenes into two media. Nothing in the PDF is computed in the
PDF, which is why the HTML cannot disagree with it and why the golden-file test
in the design's §13 only has to cover ``tables``.

reportlab and matplotlib are optional. They are imported inside the functions
that need them, and their absence raises :class:`~mjolnir.utils.MjolnirError`
with the install line — a report step that quietly wrote nothing is the failure
this project was written against. The HTML report and every TSV/JSON artefact
need neither.
"""

from .html import (
    render_cohort_html,
    render_html,
    scene_to_svg,
    write_cohort_html,
    write_html,
    write_scene_svgs,
)
from .pdf import (
    Scene,
    build_cohort_scenes,
    build_scenes,
    call_style,
    generated_stamp,
    scene_to_matplotlib,
    write_cohort_pdf,
    write_figures,
    write_pdf,
)
from .tables import (
    CALL_GLYPH,
    CALL_LEGEND,
    Grid,
    GridCell,
    all_checks,
    catalogue_call_rows,
    check_rows,
    cohort_drug_grid,
    cohort_headline,
    cohort_json,
    cohort_pairs,
    contamination_panel,
    contamination_rows,
    database_rows,
    disagreement_rows,
    distance_matrix_rows,
    distance_rows,
    drug_grid,
    drug_rows,
    headline_sentence,
    identity_pairs,
    lineage_pairs,
    lineage_support_rows,
    methods_pairs,
    observation_rows,
    qc_panel,
    qc_rows,
    rule_headline,
    sample_json,
    species_pairs,
    threshold_rows,
    tool_version_rows,
    unverified_rows,
    validity_pairs,
    variant_catalogue_rows,
    variant_rows,
    write_cohort_tables,
    write_json,
    write_sample_tables,
    write_tables,
    write_tsv,
)

__all__ = [
    # tables: the data layer
    "CALL_GLYPH", "CALL_LEGEND", "Grid", "GridCell",
    "all_checks", "qc_panel", "contamination_panel",
    "identity_pairs", "species_pairs", "lineage_pairs", "validity_pairs",
    "methods_pairs", "cohort_pairs",
    "headline_sentence", "rule_headline", "cohort_headline",
    "drug_grid", "cohort_drug_grid",
    "drug_rows", "catalogue_call_rows", "disagreement_rows", "variant_rows",
    "variant_catalogue_rows", "check_rows", "qc_rows", "contamination_rows",
    "observation_rows", "lineage_support_rows", "database_rows", "tool_version_rows",
    "threshold_rows", "unverified_rows", "distance_rows", "distance_matrix_rows",
    "sample_json", "cohort_json",
    "write_tsv", "write_json", "write_sample_tables", "write_cohort_tables", "write_tables",
    # pdf: the deliverable
    "Scene", "build_scenes", "build_cohort_scenes", "call_style", "generated_stamp",
    "write_pdf", "write_cohort_pdf", "write_figures", "scene_to_matplotlib",
    # html: the same content, self-contained
    "render_html", "render_cohort_html", "write_html", "write_cohort_html",
    "scene_to_svg", "write_scene_svgs",
]
