"""The deliverable: a clinician-first PDF, with the charts drawn as vectors.

Design §11. Page 1 is what a clinician reads: who the sample is, what the
organism is and how confidently, the drug-by-drug table with the WHO grade and
the catalogue agreement beside it, the sample-validity verdict, and the single
most important sentence the agent produced. Page 2 is the typing and purity
evidence. Everything a bioinformatician needs is behind that, in annexes, in
full. ``--profile research`` moves the annexes forward without changing a word
of their content.

**Why the charts are not matplotlib figures pasted in.** They are drawn straight
onto the PDF canvas through :class:`Scene`, a tiny vector model of rectangles,
lines, circles and text that this module renders with reportlab's own drawing
primitives and that :mod:`mjolnir.report.html` renders as inline SVG. The reason
is mechanical: reportlab has no vector import path for matplotlib output —
svglib and pypdf would both be new dependencies — so the only way to embed a
matplotlib figure directly is to rasterise it, and a clinical report gets
printed, so a 300 dpi bitmap of a drug grid is exactly the wrong trade. One
scene, two renderers, real vectors in both media, and the PDF and the HTML
cannot drift apart because they draw the same objects.

matplotlib still has a job here: :func:`write_figures` replays the same scenes
through it to produce standalone SVG, PDF or PNG figures for a slide deck or a
manuscript. It is imported lazily, and its absence raises with the install line
rather than silently producing nothing.

**Colour carries meaning and survives a monochrome printer.** Every call cell
carries an ASCII glyph as well as a fill, the fills are ordered by luminance so
that severity still reads as darkness in greyscale, and the one call that must
never be mistaken for another — ``R (outside WHO catalogue)`` — is additionally
hatched. "No resistance determinant detected" is drawn as a near-empty cell with
a dotted border and the glyph ``ND``, which is visually nothing like the solid,
bordered cell that means a catalogue actively graded a variant as susceptible.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .. import __version__
from ..config import (
    CATALOGUES,
    MAJOR_VARIANT_FRACTION,
    MIN_DEPTH,
    MIN_MINOR_VARIANT_FRACTION,
    MIN_SHARED_CALLABLE_SITES,
    PROFILES,
    source_for,
)
from ..records import (
    CALL_NO_CALL,
    CALL_R,
    CALL_R_INTERIM,
    CALL_R_OUTSIDE_WHO,
    CALL_S,
    CALL_S_INTERIM,
    CALL_UNCERTAIN,
    CohortResult,
    SampleResult,
    STATUS_FAIL,
    STATUS_PASS,
    STATUS_WARN,
    VALIDITY_INVALID,
    VALIDITY_NOT_ASSESSED,
    VALIDITY_SUSPECT,
    VALIDITY_VALID,
)
from ..utils import LOG, MjolnirError, ensure_dir, natural_key
from . import tables as T

# ---------------------------------------------------------------------------
# Palette
# ---------------------------------------------------------------------------

PAPER = "#ffffff"
INK = "#16181d"
MUTED = "#5f6672"
LINE = "#c9ced8"
PANEL = "#f3f5f8"
ACCENT = "#1f4e9c"
HATCH_INK = "#7c828c"

#: Fill, ink and texture per resistance call.
#:
#: The luminance ordering is deliberate: R is by far the darkest cell on the
#: page and stays the darkest in greyscale, ND is the lightest, and everything
#: else falls between them in severity order. The hatch on ``R-outside-WHO``
#: exists because a photocopy must not be able to turn it into a WHO Group 1
#: call — design §5.5 rule 3 forbids the two being presented as equivalent.
CALL_STYLE: Dict[str, Dict[str, Any]] = {
    CALL_R: {"fill": "#9b1c1c", "ink": "#ffffff", "border": "#6d1212",
             "hatch": "", "dash": None, "short": "resistance (WHO Group 1)"},
    CALL_R_INTERIM: {"fill": "#d4744a", "ink": INK, "border": "#a3512e",
                     "hatch": "", "dash": None, "short": "resistance (interim grade)"},
    CALL_R_OUTSIDE_WHO: {"fill": "#dda04a", "ink": INK, "border": "#a3701f",
                         "hatch": "diag", "dash": None,
                         "short": "resistance, outside the WHO catalogue"},
    CALL_UNCERTAIN: {"fill": "#e8dba0", "ink": INK, "border": "#b3a25e",
                     "hatch": "", "dash": None, "short": "uncertain significance"},
    CALL_S_INTERIM: {"fill": "#dfeae1", "ink": "#2f5d40", "border": "#7ba98c",
                     "hatch": "", "dash": (1.4, 1.4),
                     "short": "graded not associated with resistance (interim)"},
    CALL_S: {"fill": "#a9cdb6", "ink": "#1f4a30", "border": "#3f7d55",
             "hatch": "", "dash": None,
             "short": "graded not associated with resistance"},
    CALL_NO_CALL: {"fill": "#fbfbfc", "ink": "#6b7280", "border": "#b6bcc6",
                   "hatch": "", "dash": (0.9, 1.6),
                   "short": "no resistance determinant detected"},
}

STATUS_STYLE: Dict[str, Dict[str, str]] = {
    STATUS_PASS: {"fill": "#dff0e4", "ink": "#1f5c33", "border": "#3f7d55"},
    STATUS_WARN: {"fill": "#fbe9d0", "ink": "#8a4a08", "border": "#b45309"},
    STATUS_FAIL: {"fill": "#f6d9d9", "ink": "#8a1414", "border": "#b91c1c"},
}
NOT_MEASURED_STYLE = {"fill": "#eceef1", "ink": MUTED, "border": HATCH_INK}

#: Sequential ramp for the cohort distance heatmap: dark where samples are close,
#: pale where they are far. Steps rather than a continuous gradient so that a
#: greyscale print still shows discrete bands.
DISTANCE_RAMP: Tuple[str, ...] = ("#1b3a63", "#27528a", "#3a6fae", "#6f9bcb",
                                  "#a8c3de", "#d5e1ee")
FAR_RAMP: Tuple[str, ...] = ("#d9dce1", "#e6e8ec", "#f1f2f5")

#: Legend order: most alarming first. The vocabulary's own order is declaration
#: order, which reads as arbitrary in a swatch block.
LEGEND_ORDER: Tuple[str, ...] = (CALL_R, CALL_R_OUTSIDE_WHO, CALL_R_INTERIM,
                                 CALL_UNCERTAIN, CALL_S, CALL_S_INTERIM, CALL_NO_CALL)

#: The longest annex table Mjolnir will typeset. Beyond this the PDF says how
#: many rows it dropped and where the complete table is; it never truncates
#: silently.
MAX_ANNEX_ROWS = 400

CONTENT_WIDTH = 523.0


def call_style(call: str) -> Dict[str, Any]:
    """Fill, ink, border and texture for a call, with a loud default.

    An unknown call gets the uncertain style and keeps its own text, rather than
    being coerced into one of the seven — a call Mjolnir does not recognise is a
    bug, and it should look odd on the page.
    """
    return CALL_STYLE.get(call, CALL_STYLE[CALL_UNCERTAIN])


# ---------------------------------------------------------------------------
# The vector scene
# ---------------------------------------------------------------------------
#
# Coordinates are points, origin top-left, y increasing downward — the same
# convention as SVG, so the HTML renderer is a direct transcription and the PDF
# renderer flips y once, in one place.

@dataclass
class Rect:
    x: float
    y: float
    w: float
    h: float
    fill: Optional[str] = None
    stroke: Optional[str] = None
    stroke_width: float = 0.6
    dash: Optional[Tuple[float, float]] = None
    radius: float = 0.0
    title: str = ""


@dataclass
class Line:
    x1: float
    y1: float
    x2: float
    y2: float
    stroke: str = INK
    width: float = 0.6
    dash: Optional[Tuple[float, float]] = None


@dataclass
class Text:
    x: float
    y: float
    text: str
    size: float = 7.0
    fill: str = INK
    #: "start", "middle" or "end" — the same vocabulary SVG uses.
    anchor: str = "start"
    bold: bool = False
    #: Degrees counter-clockwise, for vertical column headers.
    rotate: float = 0.0


@dataclass
class Circle:
    cx: float
    cy: float
    r: float
    fill: Optional[str] = None
    stroke: Optional[str] = None
    stroke_width: float = 0.5
    title: str = ""


@dataclass
class Poly:
    points: List[Tuple[float, float]] = field(default_factory=list)
    fill: Optional[str] = None
    stroke: Optional[str] = None
    stroke_width: float = 0.5
    close: bool = True


@dataclass
class Scene:
    """A finished drawing, in points, ready for either renderer."""

    width: float
    height: float
    items: List[Any] = field(default_factory=list)
    title: str = ""
    caption: str = ""

    def add(self, *items: Any) -> None:
        self.items.extend(items)


def _hatch(scene: Scene, x: float, y: float, w: float, h: float,
           colour: str = HATCH_INK, step: float = 2.6, width: float = 0.4) -> None:
    """Diagonal hatching, emitted as clipped line segments.

    Hatching is built here rather than expressed as a fill attribute so that
    neither renderer needs a pattern concept: the PDF and the SVG both receive
    the same handful of line segments.
    """
    x1, y1 = x + w, y + h
    offset = x - y1
    end = x1 - y
    while offset <= end:
        ax = max(x, y + offset)
        ay = ax - offset
        bx = min(x1, y1 + offset)
        by = bx - offset
        if bx > ax:
            scene.add(Line(ax, ay, bx, by, stroke=colour, width=width))
        offset += step


def _cell(scene: Scene, x: float, y: float, w: float, h: float, call: str,
          glyph: str, title: str = "", not_evaluable: bool = False) -> None:
    """One call cell: fill, border, optional hatch, and the glyph on top."""
    style = call_style(call)
    if not_evaluable:
        scene.add(Rect(x, y, w, h, fill=NOT_MEASURED_STYLE["fill"],
                       stroke=NOT_MEASURED_STYLE["border"], stroke_width=0.5,
                       title=title or "target regions not callable"))
        _hatch(scene, x, y, w, h, colour=HATCH_INK, step=3.0)
        scene.add(Text(x + w / 2.0, y + h / 2.0 + 2.4, T.FLAG_NOT_EVALUABLE,
                       size=6.4, fill=MUTED, anchor="middle", bold=True))
        return
    scene.add(Rect(x, y, w, h, fill=style["fill"], stroke=style["border"],
                   stroke_width=0.5, dash=style["dash"], title=title))
    if style["hatch"] == "diag":
        _hatch(scene, x, y, w, h, colour=style["ink"], step=2.6, width=0.35)
    scene.add(Text(x + w / 2.0, y + h / 2.0 + 2.4, glyph, size=6.8,
                   fill=style["ink"], anchor="middle", bold=True))


def _legend(scene: Scene, x: float, y: float, width: float,
            calls: Sequence[str] = LEGEND_ORDER, columns: int = 2) -> float:
    """Swatch legend. Returns the y below it."""
    entry_h = 12.0
    col_w = width / float(columns)
    for index, call in enumerate(calls):
        style = call_style(call)
        col = index % columns
        row = index // columns
        cx = x + col * col_w
        cy = y + row * entry_h
        _cell(scene, cx, cy, 16.0, 9.0, call, T.CALL_GLYPH.get(call, "?"))
        scene.add(Text(cx + 20.0, cy + 7.0, style["short"], size=6.4, fill=MUTED))
    rows = int(math.ceil(len(calls) / float(columns)))
    return y + rows * entry_h


# ---------------------------------------------------------------------------
# Charts
# ---------------------------------------------------------------------------

def drug_grid_scene(grid: T.Grid, cell_w: float = 44.0, cell_h: float = 14.0,
                    label_w: float = 108.0, flag_w: float = 150.0) -> Scene:
    """The drug grid: every catalogue's answer beside the consensus.

    The consensus column is separated by a heavier rule and its cells are drawn
    one point taller, because it is the call that goes in the patient's record
    and the three columns to its left are the evidence for it.
    """
    if grid.is_empty:
        return _empty_scene("No drug calls were made for this sample.")
    header_h = 26.0
    n_rows = len(grid.row_labels)
    width = label_w + cell_w * len(grid.column_labels) + flag_w
    body_h = header_h + n_rows * cell_h
    legend_h = 12.0 * int(math.ceil(len(LEGEND_ORDER) / 2.0)) + 10.0
    scene = Scene(width=width, height=body_h + legend_h + 6.0,
                  title="Drug by catalogue", caption=grid.caption)

    for index, column in enumerate(grid.column_labels):
        cx = label_w + index * cell_w
        bold = column == T.CONSENSUS_COLUMN
        scene.add(Text(cx + cell_w / 2.0, header_h - 6.0, column, size=6.6,
                       fill=INK if bold else MUTED, anchor="middle", bold=bold))
    scene.add(Text(label_w + cell_w * len(grid.column_labels) + 4.0, header_h - 6.0,
                   "flags", size=6.6, fill=MUTED))
    scene.add(Line(0, header_h - 2.0, width, header_h - 2.0, stroke=LINE, width=0.7))

    has_consensus = grid.column_labels and grid.column_labels[-1] == T.CONSENSUS_COLUMN
    if has_consensus:
        consensus_x = label_w + cell_w * (len(grid.column_labels) - 1)
        scene.add(Line(consensus_x - 2.0, 2.0, consensus_x - 2.0, body_h,
                       stroke=LINE, width=1.1))

    for row_index in range(n_rows):
        y = header_h + row_index * cell_h
        if row_index % 2 == 1:
            scene.add(Rect(0, y, width, cell_h, fill=PANEL, stroke=None))
        flags = grid.row_flags[row_index] if row_index < len(grid.row_flags) else []
        disagreeing = any(f.startswith(T.FLAG_DISAGREEMENT) for f in flags)
        if disagreeing:
            # The front page must show a catalogue disagreement, not only the
            # annex: a bar in the margin plus the mark in the flag column.
            scene.add(Rect(0, y, 2.4, cell_h, fill=ACCENT, stroke=None))
        label = grid.row_labels[row_index]
        sub = grid.row_sublabels[row_index] if row_index < len(grid.row_sublabels) else ""
        scene.add(Text(6.0, y + cell_h - 4.5, label, size=7.0, fill=INK, bold=True))
        if sub:
            scene.add(Text(34.0, y + cell_h - 4.5, _ellipsis(sub, 18), size=6.2, fill=MUTED))
        for col_index, cell in enumerate(grid.cells[row_index]):
            cx = label_w + col_index * cell_w
            inset = 0.8 if (has_consensus and col_index == len(grid.column_labels) - 1) else 1.5
            _cell(scene, cx + inset, y + inset, cell_w - 2 * inset, cell_h - 2 * inset,
                  cell.call, cell.glyph, title=cell.detail,
                  not_evaluable=cell.not_evaluable)
        if flags:
            scene.add(Text(label_w + cell_w * len(grid.column_labels) + 4.0,
                           y + cell_h - 4.5, _ellipsis("; ".join(flags), 46),
                           size=6.0, fill=ACCENT if disagreeing else MUTED))

    _legend(scene, 0.0, body_h + 6.0, width)
    return scene


def coverage_strip_scene(result: SampleResult, width: float = CONTENT_WIDTH) -> Scene:
    """Depth and coverage against the thresholds that judge them.

    Two axes, because depth and fractions are not the same quantity and a single
    normalised bar would hide that 12x is a different kind of problem from 88%
    breadth. Every threshold gets a tick with its value, and a metric that was
    not measured gets a hatched track and the reason — never a bar at zero.
    """
    panel = T.qc_panel(result)
    depth_specs = [s for s in T.QC_METRIC_SPECS if s.unit == "x"]
    frac_specs = [s for s in T.QC_METRIC_SPECS if s.unit == "fraction"]
    by_label = dict((s.label, T.find_check(panel, s.key, s.label))
                    for s in T.QC_METRIC_SPECS)

    label_w, value_w = 148.0, 96.0
    track_w = width - label_w - value_w
    row_h, track_h = 19.0, 9.0
    header_h = 12.0

    observed = [float(by_label[spec.label].value) * 1.25 for spec in depth_specs
                if by_label.get(spec.label) is not None
                and isinstance(by_label[spec.label].value, (int, float))]
    # Round the axis up to a multiple of ten so the tick labels are readable
    # numbers; the bars and the thresholds are drawn from the real values.
    depth_max = math.ceil(max([MIN_DEPTH * 1.6] + observed) / 10.0) * 10.0

    height = (header_h + len(depth_specs) * row_h + 12.0
              + header_h + len(frac_specs) * row_h + 8.0)
    scene = Scene(width=width, height=height, title="Coverage and depth",
                  caption="Bars are the measured value; the caret under each track is the "
                          "threshold and the source is in the QC annex. A hatched track "
                          "means the metric was not measured.")

    y = 0.0
    y = _metric_block(scene, "read depth (x)", depth_specs, by_label, y, label_w,
                      track_w, row_h, track_h, header_h,
                      scale_max=depth_max, is_fraction=False)
    y += 12.0
    _metric_block(scene, "coverage and composition (%)", frac_specs, by_label, y,
                  label_w, track_w, row_h, track_h, header_h,
                  scale_max=1.0, is_fraction=True)
    return scene


def _metric_block(scene: Scene, title: str, specs: Sequence[T.MetricSpec],
                  by_label: Dict[str, Any], y: float, label_w: float, track_w: float,
                  row_h: float, track_h: float, header_h: float,
                  scale_max: float, is_fraction: bool) -> float:
    """One axis and its rows. Returns the y below the block."""
    scene.add(Text(0.0, y + 8.0, title, size=6.6, fill=MUTED, bold=True))
    axis_y = y + header_h - 3.0
    scene.add(Line(label_w, axis_y, label_w + track_w, axis_y, stroke=LINE, width=0.5))
    for step in range(0, 5):
        fraction = step / 4.0
        tx = label_w + fraction * track_w
        scene.add(Line(tx, axis_y - 2.0, tx, axis_y, stroke=LINE, width=0.5))
        label = ("{0:.0f}%".format(fraction * 100.0) if is_fraction
                 else "{0:.0f}".format(fraction * scale_max))
        scene.add(Text(tx, axis_y - 3.5, label, size=5.6, fill=MUTED, anchor="middle"))

    row_y = y + header_h
    for spec in specs:
        check = by_label.get(spec.label)
        scene.add(Text(0.0, row_y + track_h - 1.0, _ellipsis(spec.label, 34),
                       size=6.4, fill=INK))
        track_x = label_w
        scene.add(Rect(track_x, row_y, track_w, track_h, fill="#f7f8fa",
                       stroke=LINE, stroke_width=0.4))
        value = check.value if (check is not None and check.measured) else None
        if value is None or not isinstance(value, (int, float)):
            _hatch(scene, track_x, row_y, track_w, track_h, colour=HATCH_INK, step=4.5)
            why = (check.reading if check is not None else "") or "not measured"
            # The reason goes inside the track, where there is room for it. A
            # hatched bar on its own says "no data"; the sentence says which
            # capability was absent, which is the part a reader needs.
            scene.add(Text(track_x + 5.0, row_y + track_h - 2.0,
                           _ellipsis(why, 88), size=5.4, fill=INK))
            scene.add(Text(track_x + track_w + 4.0, row_y + track_h - 1.0,
                           "not measured", size=6.0, fill=MUTED, bold=True))
        else:
            style = STATUS_STYLE.get(check.status, STATUS_STYLE[STATUS_WARN])
            fraction = min(1.0, float(value) / scale_max) if scale_max else 0.0
            scene.add(Rect(track_x, row_y, max(0.8, fraction * track_w), track_h,
                           fill=style["fill"], stroke=style["border"], stroke_width=0.5,
                           title="{0} = {1}".format(spec.label, value)))
            printed = (T.fmt_fraction(float(value)) if is_fraction
                       else "{0:.1f}x".format(float(value)))
            scene.add(Text(track_x + track_w + 4.0, row_y + track_h - 1.0,
                           "{0}  {1}".format(printed, check.status), size=6.0,
                           fill=style["ink"], bold=True))
        for bound, _kind in spec.bounds():
            if not scale_max:
                continue
            bx = track_x + min(1.0, float(bound) / scale_max) * track_w
            scene.add(Poly([(bx - 2.2, row_y + track_h + 3.4), (bx + 2.2, row_y + track_h + 3.4),
                            (bx, row_y + track_h + 0.6)], fill=INK, stroke=None))
            printed = ("{0:.0f}%".format(bound * 100.0) if is_fraction
                       else "{0:g}x".format(bound))
            scene.add(Text(bx, row_y + track_h + 8.6, printed, size=5.2, fill=MUTED,
                           anchor="middle"))
        row_y += row_h
    return row_y


def allele_fraction_scene(result: SampleResult, width: float = CONTENT_WIDTH,
                          height: float = 168.0) -> Scene:
    """Allele fraction against genome position, coloured by catalogue call.

    The two horizontal rules are the thresholds that decide what a variant *is*:
    the major-variant line at 0.90 and the floor below which Mjolnir will not
    report a minor allele at all. Variants with no allele fraction — every
    variant on assembly input — are drawn in their own lane below the axis
    rather than dropped, so a FASTA sample does not look like a sample with no
    variants.
    """
    points = T.catalogue_variant_points(result)
    if not points:
        return _empty_scene("No variants were called for this sample.")
    missing = [p for p in points if p["allele_fraction"] is None]
    left, right, top = 34.0, 12.0, 14.0
    # The lane for variants with no allele fraction only takes space when there
    # are any; on Illumina it usually collapses to nothing.
    lane_h = 24.0 if missing else 0.0
    axis_h = 26.0
    plot_h = height - top - axis_h - lane_h
    plot_w = width - left - right

    positions = [p["pos"] for p in points]
    span_max = result.qc.reference_length or (max(positions) if positions else 1)
    span_min = 0 if result.qc.reference_length else min(positions)
    if span_max <= span_min:
        span_max = span_min + 1
    axis_note = ("position on {0}".format(result.qc.reference or result.reference or "the reference")
                 if result.qc.reference_length else
                 "position (axis spans the observed variants; no reference length recorded)")

    scene = Scene(width=width, height=height, title="Allele fraction at variant positions",
                  caption="Filled circles are variants a catalogue graded; open grey circles "
                          "are variants no catalogue grades. Circle area grows with read depth.")

    def px(pos: float) -> float:
        return left + (float(pos) - span_min) / float(span_max - span_min) * plot_w

    def py(fraction: float) -> float:
        return top + (1.0 - float(fraction)) * plot_h

    scene.add(Rect(left, top, plot_w, plot_h, fill="#fcfcfd", stroke=LINE, stroke_width=0.5))
    for value in (0.0, 0.25, 0.5, 0.75, 1.0):
        gy = py(value)
        scene.add(Line(left, gy, left + plot_w, gy, stroke="#e6e9ee", width=0.4))
        scene.add(Text(left - 3.0, gy + 2.0, "{0:.2f}".format(value), size=5.6,
                       fill=MUTED, anchor="end"))
    major_y = py(MAJOR_VARIANT_FRACTION)
    scene.add(Line(left, major_y, left + plot_w, major_y, stroke=INK, width=0.7))
    scene.add(Text(left + plot_w, major_y - 2.5,
                   "major variant >= {0:.2f}".format(MAJOR_VARIANT_FRACTION),
                   size=5.6, fill=INK, anchor="end"))
    minor_y = py(MIN_MINOR_VARIANT_FRACTION)
    scene.add(Line(left, minor_y, left + plot_w, minor_y, stroke=MUTED, width=0.5,
                   dash=(2.0, 2.0)))
    scene.add(Text(left + plot_w, minor_y + 6.0,
                   "minor-variant reporting floor {0:.2f}".format(MIN_MINOR_VARIANT_FRACTION),
                   size=5.6, fill=MUTED, anchor="end"))

    # Position axis. Without ticks the horizontal placement of a point carries
    # no information at all, which is worse than not drawing it.
    axis_y = top + plot_h
    scene.add(Line(left, axis_y, left + plot_w, axis_y, stroke=LINE, width=0.5))
    for step in range(0, 5):
        position = span_min + (span_max - span_min) * step / 4.0
        tx = px(position)
        scene.add(Line(tx, axis_y, tx, axis_y + 2.5, stroke=LINE, width=0.5))
        scene.add(Text(tx, axis_y + 9.0, _position_label(position), size=5.6,
                       fill=MUTED, anchor="middle"))
    lane_y = axis_y + axis_h - 6.0
    labelled = 0
    for point in points:
        depth = point["depth"] or 0
        radius = 1.9 + min(3.2, math.sqrt(float(depth)) / 5.0 if depth else 0.0)
        style = call_style(point["call"])
        graded = point["catalogued"]
        title = "{0} | {1} | depth {2} | AF {3}".format(
            point["label"], T.call_label(point["call"]),
            point["depth"] if point["depth"] is not None else "NA",
            T.fmt_number(point["allele_fraction"], 3, na="NA"))
        if point["allele_fraction"] is None:
            x = px(point["pos"])
            scene.add(Line(x, lane_y, x, lane_y + 7.0,
                           stroke=style["border"] if graded else "#b6bcc6", width=0.9))
            continue
        scene.add(Circle(px(point["pos"]), py(point["allele_fraction"]), radius,
                         fill=style["fill"] if graded else "#e9ebef",
                         stroke=style["border"] if graded else "#b6bcc6",
                         stroke_width=0.5, title=title))
        if graded and labelled < 6 and point["gene"]:
            # Stagger alternate labels so two variants a few kilobases apart do
            # not overprint each other.
            offset = -(radius + 2.0) if labelled % 2 == 0 else (radius + 6.0)
            labelled += 1
            scene.add(Text(px(point["pos"]) + radius + 1.5,
                           py(point["allele_fraction"]) + offset,
                           _ellipsis(point["label"], 22), size=5.4, fill=INK))

    if missing:
        scene.add(Text(left, lane_y + 17.0,
                       "{0} variant(s) with no allele fraction, marked as ticks: "
                       "{1}".format(
                           len(missing),
                           "assembly input carries none, which is a capability loss "
                           "and not a clean result" if result.platform == "fasta"
                           else "the caller reported none"),
                       size=5.8, fill=MUTED))
    scene.add(Text(left + plot_w, top - 4.0, axis_note, size=5.8,
                   fill=MUTED, anchor="end"))
    return scene


def _position_label(position: float) -> str:
    """A genome coordinate at a readable magnitude."""
    if position >= 1e6:
        return "{0:.1f} Mb".format(position / 1e6)
    if position >= 1e3:
        return "{0:.0f} kb".format(position / 1e3)
    return "{0:.0f}".format(position)


def distance_matrix_scene(cohort: CohortResult, width: float = CONTENT_WIDTH) -> Scene:
    """The pairwise SNP matrix, with uncompared pairs hatched rather than zeroed."""
    samples = list(cohort.samples)
    if len(samples) < 2:
        return _empty_scene("A distance matrix needs at least two samples.")
    label_w = 96.0
    cell = min(26.0, (width - label_w) / float(len(samples)))
    header_h = 62.0
    size = cell * len(samples)
    legend_w = 470.0
    scene = Scene(width=max(label_w + size + 60.0, legend_w),
                  height=header_h + size + 46.0,
                  title="Pairwise masked SNP distance",
                  caption="Hatched cells were never compared; they are absent distances, "
                          "not zeros. The shared callable denominator for every pair is in "
                          "the cohort annex.")
    matrix = cohort.distance_matrix()
    values = [v for row in matrix.values() for v in row.values() if v is not None]
    vmax = max(values) if values else 1
    threshold = cohort.threshold if cohort.threshold is not None else 0

    for index, sample in enumerate(samples):
        scene.add(Text(label_w - 4.0, header_h + index * cell + cell / 2.0 + 2.2,
                       _ellipsis(sample, 18), size=6.0, fill=INK, anchor="end"))
        scene.add(Text(label_w + index * cell + cell / 2.0, header_h - 4.0,
                       _ellipsis(sample, 16), size=6.0, fill=INK, anchor="start",
                       rotate=90.0))

    for row_index, a in enumerate(samples):
        for col_index, b in enumerate(samples):
            x = label_w + col_index * cell
            y = header_h + row_index * cell
            value = matrix[a][b]
            if value is None:
                scene.add(Rect(x, y, cell, cell, fill="#f6f7f9", stroke=PAPER,
                               stroke_width=0.6,
                               title="{0} vs {1}: never compared".format(a, b)))
                _hatch(scene, x, y, cell, cell, colour=HATCH_INK, step=3.0)
                continue
            fill, ink = _distance_colour(value, threshold, vmax)
            shared = cohort.shared_callable_sites(a, b)
            title = "{0} vs {1}: {2} SNPs over {3} shared callable bases".format(
                a, b, value, T.fmt_number(shared, na="an unrecorded number of"))
            scene.add(Rect(x, y, cell, cell, fill=fill, stroke=PAPER, stroke_width=0.6,
                           title=title))
            if cell >= 17.0:
                scene.add(Text(x + cell / 2.0, y + cell / 2.0 + 2.2, str(value),
                               size=5.8, fill=ink, anchor="middle"))

    legend_y = header_h + size + 12.0
    cursor = 0.0
    scene.add(Text(cursor, legend_y + 7.0,
                   "<= threshold ({0} SNPs)".format(threshold), size=6.0, fill=MUTED))
    cursor += 104.0
    for colour in DISTANCE_RAMP:
        scene.add(Rect(cursor, legend_y, 13.0, 8.0, fill=colour, stroke=LINE,
                       stroke_width=0.3))
        cursor += 14.0
    cursor += 8.0
    scene.add(Text(cursor, legend_y + 7.0, "above threshold", size=6.0, fill=MUTED))
    cursor += 76.0
    for colour in FAR_RAMP:
        scene.add(Rect(cursor, legend_y, 13.0, 8.0, fill=colour, stroke=LINE,
                       stroke_width=0.3))
        cursor += 14.0
    cursor += 8.0
    scene.add(Rect(cursor, legend_y, 13.0, 8.0, fill="#f6f7f9", stroke=LINE,
                   stroke_width=0.3))
    _hatch(scene, cursor, legend_y, 13.0, 8.0, colour=HATCH_INK, step=3.0)
    cursor += 18.0
    scene.add(Text(cursor, legend_y + 7.0, "never compared", size=6.0, fill=MUTED))
    return scene


def _distance_colour(value: int, threshold: int, vmax: int) -> Tuple[str, str]:
    if threshold and value <= threshold:
        step = int(round((value / float(max(1, threshold))) * (len(DISTANCE_RAMP) - 1)))
        colour = DISTANCE_RAMP[max(0, min(len(DISTANCE_RAMP) - 1, step))]
        return colour, "#ffffff" if step <= 2 else INK
    span = max(1, vmax - threshold)
    step = int(round(((value - threshold) / float(span)) * (len(FAR_RAMP) - 1)))
    return FAR_RAMP[max(0, min(len(FAR_RAMP) - 1, step))], MUTED


def dendrogram_scene(cohort: CohortResult, width: float = CONTENT_WIDTH) -> Scene:
    """Average-linkage tree over the masked SNP distances.

    Samples whose distances are not all present are left out of the tree and
    named beneath it. Imputing a missing comparison would draw a branch that no
    measurement supports, and a dendrogram is read as if every branch were
    measured.
    """
    distances: Dict[Tuple[str, str], Optional[int]] = {}
    for pair in cohort.pairs:
        distances[(pair.sample_a, pair.sample_b)] = pair.snps
    merges, leaves, excluded = T.upgma_tree(cohort.samples, distances)
    if not merges:
        note = "Not enough complete pairwise distances to build a tree."
        if excluded:
            note += " Missing comparisons for: {0}.".format(
                ", ".join(sorted(excluded, key=natural_key)))
        return _empty_scene(note)

    order: List[int] = []

    def walk(node_id: int) -> None:
        node = next((m for m in merges if m["id"] == node_id), None)
        if node is None:
            order.append(node_id)
            return
        walk(node["left"])
        walk(node["right"])

    walk(merges[-1]["id"])

    label_w = 118.0
    row_h = 14.0
    plot_w = width - label_w - 46.0
    height = 26.0 + len(order) * row_h + 34.0
    max_distance = max([m["distance"] for m in merges] + [1.0])
    threshold = float(cohort.threshold) if cohort.threshold is not None else None
    axis_max = max(max_distance * 1.1, (threshold or 0) * 1.2, 1.0)

    scene = Scene(width=width, height=height, title="Cluster dendrogram",
                  caption="Average linkage over masked pairwise SNP distances. The dashed "
                          "line is the clustering threshold; its basis is printed with the "
                          "cluster table.")

    def px(distance: float) -> float:
        return label_w + plot_w * min(1.0, distance / axis_max)

    y_of: Dict[int, float] = {}
    for index, leaf_id in enumerate(order):
        y = 26.0 + index * row_h + row_h / 2.0
        y_of[leaf_id] = y
        name = leaves[leaf_id] if leaf_id < len(leaves) else str(leaf_id)
        cluster = cohort.cluster_of(name)
        scene.add(Text(label_w - 6.0, y + 2.2, _ellipsis(name, 22), size=6.4,
                       fill=INK, anchor="end"))
        if cluster is not None:
            scene.add(Text(label_w + plot_w + 6.0, y + 2.2, cluster.cluster_id,
                           size=5.8, fill=ACCENT))

    for node in merges:
        x = px(node["distance"])
        y1 = y_of[node["left"]]
        y2 = y_of[node["right"]]
        inside = threshold is not None and node["distance"] <= threshold
        stroke = ACCENT if inside else MUTED
        scene.add(Line(x, y1, x, y2, stroke=stroke, width=0.9))
        scene.add(Line(x, y1, px(_distance_of(merges, node["left"])), y1,
                       stroke=stroke, width=0.9))
        scene.add(Line(x, y2, px(_distance_of(merges, node["right"])), y2,
                       stroke=stroke, width=0.9))
        y_of[node["id"]] = (y1 + y2) / 2.0
        if node["distance"] > 0:
            scene.add(Text(x, min(y1, y2) - 2.0, "{0:g}".format(node["distance"]),
                           size=5.2, fill=MUTED, anchor="middle"))

    axis_y = 26.0 + len(order) * row_h + 6.0
    scene.add(Line(label_w, axis_y, label_w + plot_w, axis_y, stroke=LINE, width=0.5))
    for step in range(0, 5):
        value = axis_max * step / 4.0
        tx = px(value)
        scene.add(Line(tx, axis_y, tx, axis_y + 2.5, stroke=LINE, width=0.5))
        scene.add(Text(tx, axis_y + 9.0, "{0:g}".format(round(value)), size=5.6,
                       fill=MUTED, anchor="middle"))
    scene.add(Text(label_w + plot_w / 2.0, axis_y + 18.0, "masked SNP distance",
                   size=5.8, fill=MUTED, anchor="middle"))
    if threshold is not None:
        tx = px(threshold)
        scene.add(Line(tx, 20.0, tx, axis_y, stroke=ACCENT, width=0.8, dash=(2.5, 2.0)))
        scene.add(Text(tx, 16.0, "threshold {0:g}".format(threshold), size=5.8,
                       fill=ACCENT, anchor="middle"))
    if excluded:
        scene.caption += (" Excluded for missing comparisons: {0}.".format(
            ", ".join(sorted(excluded, key=natural_key))))
    return scene


def _distance_of(merges: Sequence[Dict[str, Any]], node_id: int) -> float:
    for node in merges:
        if node["id"] == node_id:
            return float(node["distance"])
    return 0.0


def _empty_scene(message: str) -> Scene:
    scene = Scene(width=CONTENT_WIDTH, height=26.0)
    scene.add(Rect(0, 0, CONTENT_WIDTH, 24.0, fill=PANEL, stroke=LINE, stroke_width=0.4))
    scene.add(Text(8.0, 15.0, message, size=7.0, fill=MUTED))
    return scene


def _ellipsis(text: str, limit: int) -> str:
    """Shorten a label for a fixed-width slot, with the cut made visible."""
    text = str(text or "")
    if len(text) <= limit:
        return text
    return text[: max(1, limit - 3)] + "..."


def build_scenes(result: SampleResult,
                 cohort: Optional[CohortResult] = None) -> List[Tuple[str, Scene]]:
    """Every figure for one sample, named, in report order.

    Shared with the HTML renderer so that both media show the same figures built
    from the same numbers.
    """
    scenes: List[Tuple[str, Scene]] = [
        ("drug-grid", drug_grid_scene(T.drug_grid(result))),
        ("coverage", coverage_strip_scene(result)),
        ("allele-fraction", allele_fraction_scene(result)),
    ]
    if cohort is not None and len(cohort.samples) > 1:
        scenes.append(("distance-matrix", distance_matrix_scene(cohort)))
        scenes.append(("dendrogram", dendrogram_scene(cohort)))
    return scenes


def build_cohort_scenes(cohort: CohortResult,
                        results: Sequence[SampleResult] = ()) -> List[Tuple[str, Scene]]:
    scenes: List[Tuple[str, Scene]] = []
    if results:
        scenes.append(("cohort-drug-grid", drug_grid_scene(T.cohort_drug_grid(results),
                                                           cell_w=26.0, label_w=118.0,
                                                           flag_w=120.0)))
    scenes.append(("distance-matrix", distance_matrix_scene(cohort)))
    scenes.append(("dendrogram", dendrogram_scene(cohort)))
    return scenes


# ---------------------------------------------------------------------------
# reportlab
# ---------------------------------------------------------------------------

def _require_reportlab():
    """Import reportlab or say exactly how to get it.

    A report step that quietly wrote nothing because a library was missing is
    the failure this project exists to prevent, so the absence is fatal and the
    message carries the command.
    """
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.platypus import (KeepTogether, PageBreak, Paragraph,
                                        SimpleDocTemplate, Spacer, Table, TableStyle)
        from reportlab.platypus.flowables import Flowable
    except ImportError as exc:
        raise MjolnirError(
            "reportlab is required to write the PDF report but is not installed "
            "({0}).\n  conda install -c conda-forge reportlab\n"
            "  (or: pip install 'mjolnir-myco[report]')".format(exc)
        )
    return {
        "colors": colors, "A4": A4, "ParagraphStyle": ParagraphStyle,
        "KeepTogether": KeepTogether, "PageBreak": PageBreak, "Paragraph": Paragraph,
        "SimpleDocTemplate": SimpleDocTemplate, "Spacer": Spacer, "Table": Table,
        "TableStyle": TableStyle, "Flowable": Flowable,
    }


_SCENE_FLOWABLE_CACHE: Dict[str, Any] = {}


def _scene_flowable_class(rl: Dict[str, Any]):
    """A platypus flowable that draws a :class:`Scene` at whatever width it gets.

    Built on first use rather than at import time, because subclassing
    ``Flowable`` requires reportlab to be importable and this module must import
    cleanly on a machine that only ever wanted the HTML.
    """
    cached = _SCENE_FLOWABLE_CACHE.get("class")
    if cached is not None:
        return cached

    class _SceneFlowable(rl["Flowable"]):
        def __init__(self, scene: Scene) -> None:
            rl["Flowable"].__init__(self)
            self.scene = scene
            self.scale = 1.0

        def wrap(self, availWidth, availHeight):  # noqa: N803 - reportlab's signature
            self.scale = min(1.0, float(availWidth) / max(1.0, self.scene.width))
            return (self.scene.width * self.scale, self.scene.height * self.scale)

        def draw(self) -> None:
            canvas = self.canv
            canvas.saveState()
            canvas.scale(self.scale, self.scale)
            _render_scene(canvas, self.scene)
            canvas.restoreState()

    _SCENE_FLOWABLE_CACHE["class"] = _SceneFlowable
    return _SceneFlowable


def _render_scene(canvas: Any, scene: Scene) -> None:
    """Draw a scene onto a reportlab canvas, flipping y exactly once."""
    from reportlab.lib.colors import HexColor

    def flip(y: float) -> float:
        return scene.height - y

    for item in scene.items:
        if isinstance(item, Rect):
            canvas.saveState()
            if item.dash:
                canvas.setDash(list(item.dash))
            if item.stroke:
                canvas.setStrokeColor(HexColor(item.stroke))
                canvas.setLineWidth(item.stroke_width)
            if item.fill:
                canvas.setFillColor(HexColor(item.fill))
            mode = (1 if item.fill else 0, 1 if item.stroke else 0)
            if mode != (0, 0):
                if item.radius:
                    canvas.roundRect(item.x, flip(item.y + item.h), item.w, item.h,
                                     item.radius, stroke=mode[1], fill=mode[0])
                else:
                    canvas.rect(item.x, flip(item.y + item.h), item.w, item.h,
                                stroke=mode[1], fill=mode[0])
            canvas.restoreState()
        elif isinstance(item, Line):
            canvas.saveState()
            canvas.setStrokeColor(HexColor(item.stroke))
            canvas.setLineWidth(item.width)
            if item.dash:
                canvas.setDash(list(item.dash))
            canvas.line(item.x1, flip(item.y1), item.x2, flip(item.y2))
            canvas.restoreState()
        elif isinstance(item, Circle):
            canvas.saveState()
            if item.fill:
                canvas.setFillColor(HexColor(item.fill))
            if item.stroke:
                canvas.setStrokeColor(HexColor(item.stroke))
                canvas.setLineWidth(item.stroke_width)
            canvas.circle(item.cx, flip(item.cy), item.r,
                          stroke=1 if item.stroke else 0, fill=1 if item.fill else 0)
            canvas.restoreState()
        elif isinstance(item, Poly):
            if not item.points:
                continue
            canvas.saveState()
            path = canvas.beginPath()
            first = item.points[0]
            path.moveTo(first[0], flip(first[1]))
            for point in item.points[1:]:
                path.lineTo(point[0], flip(point[1]))
            if item.close:
                path.close()
            if item.fill:
                canvas.setFillColor(HexColor(item.fill))
            if item.stroke:
                canvas.setStrokeColor(HexColor(item.stroke))
                canvas.setLineWidth(item.stroke_width)
            canvas.drawPath(path, stroke=1 if item.stroke else 0,
                            fill=1 if item.fill else 0)
            canvas.restoreState()
        elif isinstance(item, Text):
            canvas.saveState()
            canvas.setFillColor(HexColor(item.fill))
            font = "Helvetica-Bold" if item.bold else "Helvetica"
            canvas.setFont(font, item.size)
            if item.rotate:
                canvas.translate(item.x, flip(item.y))
                canvas.rotate(item.rotate)
                _draw_anchored(canvas, 0.0, 0.0, item)
            else:
                _draw_anchored(canvas, item.x, flip(item.y), item)
            canvas.restoreState()


def _draw_anchored(canvas: Any, x: float, y: float, item: Text) -> None:
    if item.anchor == "middle":
        canvas.drawCentredString(x, y, item.text)
    elif item.anchor == "end":
        canvas.drawRightString(x, y, item.text)
    else:
        canvas.drawString(x, y, item.text)


# ---------------------------------------------------------------------------
# Document assembly
# ---------------------------------------------------------------------------

def _styles(rl: Dict[str, Any]) -> Dict[str, Any]:
    ParagraphStyle = rl["ParagraphStyle"]
    base = ParagraphStyle("body", fontName="Helvetica", fontSize=8.0, leading=10.2,
                          textColor=rl["colors"].HexColor(INK))
    return {
        "title": ParagraphStyle("title", parent=base, fontName="Helvetica-Bold",
                                fontSize=19, leading=22, spaceAfter=2),
        "subtitle": ParagraphStyle("subtitle", parent=base, fontSize=8.5, leading=11,
                                   textColor=rl["colors"].HexColor(MUTED)),
        "h2": ParagraphStyle("h2", parent=base, fontName="Helvetica-Bold", fontSize=11,
                             leading=13, spaceBefore=10, spaceAfter=4),
        "h3": ParagraphStyle("h3", parent=base, fontName="Helvetica-Bold", fontSize=8,
                             leading=10, spaceBefore=6, spaceAfter=2,
                             textColor=rl["colors"].HexColor(MUTED)),
        "body": base,
        "small": ParagraphStyle("small", parent=base, fontSize=6.8, leading=8.4,
                                textColor=rl["colors"].HexColor(MUTED)),
        "cell": ParagraphStyle("cell", parent=base, fontSize=6.6, leading=8.0),
        "cellbold": ParagraphStyle("cellbold", parent=base, fontName="Helvetica-Bold",
                                   fontSize=6.6, leading=8.0),
        "headline": ParagraphStyle("headline", parent=base, fontSize=10.5, leading=13.5,
                                   fontName="Helvetica-Bold"),
    }


def _xml(text: Any) -> str:
    """Escape for reportlab's mini-markup, which is XML."""
    return (str("" if text is None else text)
            .replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _pairs_table(rl: Dict[str, Any], styles: Dict[str, Any],
                 pairs: Sequence[Tuple[str, str]], label_width: float = 118.0,
                 width: float = CONTENT_WIDTH) -> Any:
    Paragraph, Table, TableStyle = rl["Paragraph"], rl["Table"], rl["TableStyle"]
    colors = rl["colors"]
    data = [[Paragraph("<b>{0}</b>".format(_xml(label)), styles["cell"]),
             Paragraph(_xml(value), styles["cell"])] for label, value in pairs]
    if not data:
        data = [[Paragraph("nothing recorded", styles["cell"]), Paragraph("", styles["cell"])]]
    table = Table(data, colWidths=[label_width, width - label_width], hAlign="LEFT")
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 1.6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("LINEBELOW", (0, 0), (-1, -2), 0.25, colors.HexColor(LINE)),
    ]))
    return table


def _row_table(rl: Dict[str, Any], styles: Dict[str, Any], rows: Sequence[Dict[str, Any]],
               columns: Sequence[Tuple[str, str, float]],
               status_column: str = "") -> List[Any]:
    """A generic annex table: (key, header, width) triples, capped at a stated size."""
    Paragraph, Table, TableStyle = rl["Paragraph"], rl["Table"], rl["TableStyle"]
    colors = rl["colors"]
    flowables: List[Any] = []
    if not rows:
        flowables.append(Paragraph("nothing to show in this table", styles["small"]))
        return flowables
    shown = list(rows[:MAX_ANNEX_ROWS])
    header = [Paragraph("<b>{0}</b>".format(_xml(title)), styles["cell"])
              for _key, title, _w in columns]
    data = [header]
    for row in shown:
        data.append([Paragraph(_xml(_cell_text(row.get(key))), styles["cell"])
                     for key, _title, _w in columns])
    table = Table(data, colWidths=[w for _k, _t, w in columns], repeatRows=1, hAlign="LEFT")
    style = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(PANEL)),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor(LINE)),
        ("GRID", (0, 0), (-1, -1), 0.2, colors.HexColor("#e3e6ec")),
        ("TOPPADDING", (0, 0), (-1, -1), 1.4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.4),
        ("LEFTPADDING", (0, 0), (-1, -1), 2.5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2.5),
    ]
    if status_column:
        index = [k for k, _t, _w in columns].index(status_column)
        for row_index, row in enumerate(shown, start=1):
            value = str(row.get(status_column) or "")
            palette = (NOT_MEASURED_STYLE if value == "not measured"
                       else STATUS_STYLE.get(value))
            if palette:
                style.append(("BACKGROUND", (index, row_index), (index, row_index),
                              colors.HexColor(palette["fill"])))
                style.append(("TEXTCOLOR", (index, row_index), (index, row_index),
                              colors.HexColor(palette["ink"])))
    table.setStyle(TableStyle(style))
    flowables.append(table)
    if len(rows) > MAX_ANNEX_ROWS:
        flowables.append(Paragraph(
            "{0} of {1} rows are shown here; the complete table is in the TSV "
            "artefact written beside this report.".format(MAX_ANNEX_ROWS, len(rows)),
            styles["small"]))
    return flowables


def _cell_text(value: Any) -> str:
    if value is None:
        return T.TSV_NA
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return "{0:.4f}".format(value).rstrip("0").rstrip(".")
    if isinstance(value, (list, tuple)):
        return "; ".join(_cell_text(v) for v in value)
    return str(value)


def _drug_table(rl: Dict[str, Any], styles: Dict[str, Any], result: SampleResult) -> Any:
    """Page 1's drug-by-drug table: call, WHO grade, catalogue agreement, flags."""
    Paragraph, Table, TableStyle = rl["Paragraph"], rl["Table"], rl["TableStyle"]
    colors = rl["colors"]
    header = ["Drug", "Determinant", "Mjolnir call", "WHO grade"]
    header += [c.replace("WHO v2", "WHO") for c in CATALOGUES]
    header += ["Flags"]
    widths = [66.0, 112.0, 104.0, 74.0, 24.0, 30.0, 24.0, 89.0]
    data = [[Paragraph("<b>{0}</b>".format(_xml(h)), styles["cell"]) for h in header]]

    rows = T.drug_rows(result)
    calls = sorted(result.drugs, key=lambda d: T.drug_order(d.drug))
    style_commands = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(PANEL)),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#e3e6ec")),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, colors.HexColor(LINE)),
        ("TOPPADDING", (0, 0), (-1, -1), 2.0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.0),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("ALIGN", (4, 0), (6, -1), "CENTER"),
    ]
    for index, (row, call) in enumerate(zip(rows, calls), start=1):
        determinant = "; ".join(row["supporting_variants"]) or "none matched"
        flags = "; ".join(T.drug_flags(call))
        cells = [
            Paragraph('<b>{0}</b><br/><font size="5.6" color="{1}">{2}</font>'.format(
                _xml(row["drug_code"]), MUTED, _xml(row["drug"])), styles["cell"]),
            Paragraph(_xml(determinant), styles["cell"]),
            Paragraph("<b>{0}</b> {1}".format(_xml(row["call_glyph"]),
                                              _xml(row["call_label"])), styles["cell"]),
            Paragraph(_xml(row["who_grade"] or "not graded by WHO"), styles["cell"]),
        ]
        for catalogue in CATALOGUES:
            key = T.catalogue_key(catalogue)
            cells.append(Paragraph("<b>{0}</b>".format(
                _xml(T.CALL_GLYPH.get(row[key + "_call"], "?"))), styles["cellbold"]))
        cells.append(Paragraph(_xml(flags), styles["cell"]))
        data.append(cells)

        style = call_style(row["call"])
        style_commands.append(("BACKGROUND", (2, index), (2, index),
                               colors.HexColor(style["fill"])))
        style_commands.append(("TEXTCOLOR", (2, index), (2, index),
                               colors.HexColor(style["ink"])))
        for offset, catalogue in enumerate(CATALOGUES):
            key = T.catalogue_key(catalogue)
            cat_style = call_style(row[key + "_call"])
            column = 4 + offset
            style_commands.append(("BACKGROUND", (column, index), (column, index),
                                   colors.HexColor(cat_style["fill"])))
            style_commands.append(("TEXTCOLOR", (column, index), (column, index),
                                   colors.HexColor(cat_style["ink"])))
        if row["disagreement"]:
            # The design requires a drug carrying a catalogue disagreement to be
            # visible on page 1, not merely listed in an annex.
            style_commands.append(("LINEBEFORE", (0, index), (0, index), 2.2,
                                   colors.HexColor(ACCENT)))
            style_commands.append(("TEXTCOLOR", (7, index), (7, index),
                                   colors.HexColor(ACCENT)))
    if not rows:
        data.append([Paragraph("no drug was evaluated for this sample", styles["cell"])]
                    + [Paragraph("", styles["cell"])] * (len(header) - 1))
    table = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    table.setStyle(TableStyle(style_commands))
    return table


def _status_chip(rl: Dict[str, Any], styles: Dict[str, Any], label: str, value: str,
                 palette: Dict[str, str], width: float = CONTENT_WIDTH) -> Any:
    Paragraph, Table, TableStyle = rl["Paragraph"], rl["Table"], rl["TableStyle"]
    colors = rl["colors"]
    table = Table([[Paragraph(
        '<font size="7" color="{0}">{1}</font><br/>'
        '<b><font size="12">{2}</font></b>'.format(
            palette["ink"], _xml(label), _xml(value)), styles["cell"])]],
        colWidths=[width], hAlign="LEFT")
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(palette["fill"])),
        ("BOX", (0, 0), (-1, -1), 0.8, colors.HexColor(palette["border"])),
        ("TEXTCOLOR", (0, 0), (-1, -1), colors.HexColor(palette["ink"])),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
    ]))
    return table


_VALIDITY_PALETTE = {
    VALIDITY_VALID: STATUS_STYLE[STATUS_PASS],
    VALIDITY_SUSPECT: STATUS_STYLE[STATUS_WARN],
    VALIDITY_INVALID: STATUS_STYLE[STATUS_FAIL],
    VALIDITY_NOT_ASSESSED: {"fill": NOT_MEASURED_STYLE["fill"],
                            "ink": NOT_MEASURED_STYLE["ink"],
                            "border": NOT_MEASURED_STYLE["border"]},
}


def _section_identity(rl, styles, result: SampleResult) -> List[Any]:
    Paragraph, Spacer = rl["Paragraph"], rl["Spacer"]
    headline, provenance = T.headline_sentence(result)
    out: List[Any] = [
        Paragraph("Mjolnir report", styles["title"]),
        Paragraph("{0} | {1} | Mjolnir {2} | profile {3}".format(
            _xml(result.sample_id), _xml(result.platform),
            _xml(result.mjolnir_version or __version__), _xml(result.profile)),
            styles["subtitle"]),
        Spacer(1, 8),
        _status_chip(rl, styles, "Sample validity",
                     result.contamination.verdict,
                     _VALIDITY_PALETTE.get(result.contamination.verdict,
                                           STATUS_STYLE[STATUS_WARN])),
        Spacer(1, 4),
        Paragraph(_xml(result.contamination.verdict_reason
                       or "no reason for this verdict was recorded"), styles["small"]),
        Spacer(1, 8),
        Paragraph("Reading", styles["h3"]),
        Paragraph(_xml(headline), styles["headline"]),
        Paragraph("Source of this sentence: {0}. Pass, warn and fail are computed in "
                  "Python from the thresholds in the annex before any model is called."
                  .format(_xml(provenance)), styles["small"]),
        Spacer(1, 8),
        Paragraph("Sample", styles["h3"]),
        _pairs_table(rl, styles, T.identity_pairs(result)),
        Spacer(1, 6),
        Paragraph("Species", styles["h3"]),
        _pairs_table(rl, styles, T.species_pairs(result)),
    ]
    return out


def _section_drugs(rl, styles, result: SampleResult, scenes: Dict[str, Scene]) -> List[Any]:
    Paragraph, Spacer = rl["Paragraph"], rl["Spacer"]
    flowable = _scene_flowable_class(rl)
    out: List[Any] = [
        Paragraph("Drug-by-drug resistance", styles["h2"]),
        _drug_table(rl, styles, result),
        Spacer(1, 5),
    ]
    grid = scenes.get("drug-grid")
    if grid is not None:
        out.append(flowable(grid))
        out.append(Paragraph(_xml(grid.caption), styles["small"]))
    out.append(Spacer(1, 4))
    out.append(Paragraph("How to read a call", styles["h3"]))
    for _call, text in T.CALL_LEGEND:
        out.append(Paragraph(_xml(text), styles["small"]))
    caveats = T.caveat_lines(result)
    if caveats:
        out.append(Spacer(1, 4))
        out.append(Paragraph("Platform and sample caveats", styles["h3"]))
        for caveat in caveats:
            out.append(Paragraph("- " + _xml(caveat), styles["small"]))
    return out


def _section_typing(rl, styles, result: SampleResult, scenes: Dict[str, Scene]) -> List[Any]:
    Paragraph, Spacer = rl["Paragraph"], rl["Spacer"]
    flowable = _scene_flowable_class(rl)
    out: List[Any] = [
        Paragraph("Lineage and barcode evidence", styles["h2"]),
        _pairs_table(rl, styles, T.lineage_pairs(result)),
        Spacer(1, 8),
        Paragraph("Contamination and sample validity", styles["h2"]),
        _pairs_table(rl, styles, T.validity_pairs(result)),
        Spacer(1, 5),
    ]
    out.extend(_row_table(rl, styles, T.contamination_rows(result), (
        ("check", "metric", 132.0), ("value", "value", 52.0),
        ("comparison", "vs", 22.0), ("threshold", "threshold", 50.0),
        ("status", "status", 46.0), ("source", "source", 221.0),
    ), status_column="status"))
    coverage = scenes.get("coverage")
    if coverage is not None:
        out.append(Spacer(1, 8))
        out.append(Paragraph("Coverage and depth", styles["h2"]))
        out.append(flowable(coverage))
        out.append(Paragraph(_xml(coverage.caption), styles["small"]))
    return out


def _section_variants(rl, styles, result: SampleResult, scenes: Dict[str, Scene]) -> List[Any]:
    Paragraph, Spacer = rl["Paragraph"], rl["Spacer"]
    flowable = _scene_flowable_class(rl)
    out: List[Any] = [Paragraph("Annex A - variants", styles["h2"])]
    plot = scenes.get("allele-fraction")
    if plot is not None:
        out.append(flowable(plot))
        out.append(Paragraph(_xml(plot.caption), styles["small"]))
        out.append(Spacer(1, 5))
    columns = [("coordinate", "coordinate", 92.0), ("gene", "gene", 44.0),
               ("hgvs", "HGVS", 84.0), ("variant_type", "type", 34.0),
               ("effect", "effect", 52.0), ("depth", "depth", 28.0),
               ("allele_fraction", "AF", 32.0), ("is_major", "major", 26.0)]
    for catalogue in CATALOGUES:
        columns.append((T.catalogue_key(catalogue) + "_grade",
                        catalogue.replace("WHO v2", "WHO"), 43.66))
    out.extend(_row_table(rl, styles, T.variant_rows(result), columns))
    return out


def _section_disagreements(rl, styles, result: SampleResult) -> List[Any]:
    Paragraph = rl["Paragraph"]
    rows = T.disagreement_rows(result)
    out: List[Any] = [Paragraph("Annex B - catalogue disagreements", styles["h2"])]
    if not rows:
        out.append(Paragraph(
            "No drug in this sample carries a catalogue disagreement. That is a "
            "statement about the three catalogues agreeing, not about the organism.",
            styles["small"]))
        return out
    out.extend(_row_table(rl, styles, rows, (
        ("drug", "drug", 62.0), ("mjolnir_call", "Mjolnir", 62.0),
        ("catalogue", "catalogue", 52.0), ("catalogue_call", "its call", 62.0),
        ("grade", "grade", 78.0), ("disagreement_kind", "kind", 46.0),
        ("note", "note", 161.0),
    )))
    return out


def _section_qc(rl, styles, result: SampleResult) -> List[Any]:
    Paragraph, Spacer = rl["Paragraph"], rl["Spacer"]
    out: List[Any] = [
        Paragraph("Annex C - QC metrics against their thresholds", styles["h2"]),
        Paragraph("Every row states the threshold that was applied and the document it "
                  "came from. A row marked \"not measured\" is a gap in the evidence and "
                  "is never a pass.", styles["small"]),
        Spacer(1, 3),
    ]
    out.extend(_row_table(rl, styles, T.check_rows(T.all_checks(result), result.sample_id), (
        ("category", "panel", 44.0), ("check", "check", 116.0),
        ("value", "value", 46.0), ("comparison", "vs", 20.0),
        ("threshold", "threshold", 44.0), ("status", "status", 44.0),
        ("source", "source", 209.0),
    ), status_column="status"))
    out.append(Spacer(1, 6))
    out.append(Paragraph("Measured without a registered threshold", styles["h3"]))
    out.append(Paragraph(
        "These are reported without a pass or fail because Mjolnir registers no "
        "published bound for them, and a status implies a bound.", styles["small"]))
    out.extend(_row_table(rl, styles, T.observation_rows(result), (
        ("panel", "panel", 62.0), ("observation", "observation", 190.0),
        ("value", "value", 80.0), ("unit", "unit", 52.0), ("note", "note", 139.0),
    )))
    support = T.lineage_support_rows(result)
    keys = [k for k in support[0] if k != "sample"][:6] if support else []
    if keys:
        out.append(Spacer(1, 6))
        out.append(Paragraph("Barcode sites", styles["h3"]))
        width = CONTENT_WIDTH / len(keys)
        out.extend(_row_table(rl, styles, support,
                              [(k, k.replace("_", " "), width) for k in keys]))
    return out


def _section_cohort(rl, styles, cohort: CohortResult, scenes: Dict[str, Scene]) -> List[Any]:
    Paragraph, Spacer = rl["Paragraph"], rl["Spacer"]
    flowable = _scene_flowable_class(rl)
    headline, provenance = T.cohort_headline(cohort)
    out: List[Any] = [
        Paragraph("Annex D - cohort distances and clusters", styles["h2"]),
        Paragraph(_xml(headline), styles["body"]),
        Paragraph("Source of this sentence: {0}.".format(_xml(provenance)), styles["small"]),
        Spacer(1, 4),
        _pairs_table(rl, styles, T.cohort_pairs(cohort)),
        Spacer(1, 6),
    ]
    for name in ("distance-matrix", "dendrogram"):
        scene = scenes.get(name)
        if scene is not None:
            out.append(flowable(scene))
            out.append(Paragraph(_xml(scene.caption), styles["small"]))
            out.append(Spacer(1, 5))
    out.append(Paragraph("Clusters", styles["h3"]))
    out.extend(_row_table(rl, styles, T.cluster_rows(cohort), (
        ("cluster", "cluster", 56.0), ("size", "n", 22.0),
        ("members", "members", 235.0), ("threshold", "threshold", 44.0),
        ("max_distance", "max SNPs", 46.0),
        ("min_shared_callable_sites", "min shared sites", 120.0),
    )))
    out.append(Spacer(1, 5))
    out.append(Paragraph("Pairwise distances with their denominators", styles["h3"]))
    out.append(Paragraph(
        "A distance is only comparable to the published SNP thresholds when the two "
        "samples share enough callable sequence; the floor is {0:,} bp ({1}).".format(
            MIN_SHARED_CALLABLE_SITES, _xml(source_for("min_shared_callable_sites"))),
        styles["small"]))
    out.extend(_row_table(rl, styles, T.distance_rows(cohort), (
        ("sample_a", "sample A", 96.0), ("sample_b", "sample B", 96.0),
        ("snps", "SNPs", 34.0), ("shared_callable_sites", "shared callable", 78.0),
        ("snps_per_mb", "SNPs/Mb", 48.0), ("within_threshold", "within", 38.0),
        ("shared_sites_sufficient", "sites ok", 42.0), ("note", "note", 91.0),
    )))
    return out


def _section_methods(rl, styles, result: Optional[SampleResult],
                     cohort: Optional[CohortResult]) -> List[Any]:
    Paragraph, Spacer = rl["Paragraph"], rl["Spacer"]
    out: List[Any] = [Paragraph("Annex E - methods, versions and checksums", styles["h2"])]
    if result is not None:
        out.append(_pairs_table(rl, styles, T.methods_pairs(result)))
        out.append(Spacer(1, 5))
    tools = dict(result.tool_versions) if result is not None else {}
    databases = list(result.database_versions) if result is not None else []
    if cohort is not None:
        tools.update(cohort.tool_versions)
        databases = databases + list(cohort.database_versions)
    out.append(Paragraph("Databases", styles["h3"]))
    out.append(Paragraph(
        "A catalogue-version mismatch between two installations changes calls, so every "
        "database prints its version and checksum here.", styles["small"]))
    out.extend(_row_table(rl, styles, T.database_rows(databases), (
        ("database", "database", 84.0), ("version", "version", 74.0),
        ("checksum", "checksum", 150.0), ("licence", "licence", 74.0),
        ("citation", "citation", 141.0),
    )))
    out.append(Spacer(1, 5))
    out.append(Paragraph("Tools", styles["h3"]))
    out.extend(_row_table(rl, styles, T.tool_version_rows(tools), (
        ("tool", "tool", 160.0), ("version", "version", 363.0),
    )))
    out.append(Spacer(1, 5))
    out.append(Paragraph("Thresholds whose citation has not been checked", styles["h3"]))
    out.append(Paragraph(
        "An unverified citation is worse than none, because it looks settled. These "
        "numbers are in use and nobody has yet opened the primary document on this "
        "machine.", styles["small"]))
    out.extend(_row_table(rl, styles, T.unverified_rows(), (
        ("threshold", "threshold", 130.0), ("value", "value", 60.0),
        ("source", "source", 200.0), ("note", "note", 133.0),
    )))
    return out


def _story(rl, styles, result: SampleResult, cohort: Optional[CohortResult],
           profile: str) -> List[Any]:
    """Assemble the flowables in profile order.

    ``research`` does not hide anything the clinical profile shows and does not
    show anything it hides. It changes the order only: variants, QC and the
    catalogue evidence come first, and the clinical summary follows.
    """
    PageBreak = rl["PageBreak"]
    scenes = dict(build_scenes(result, cohort))
    identity = _section_identity(rl, styles, result)
    drugs = _section_drugs(rl, styles, result, scenes)
    typing = _section_typing(rl, styles, result, scenes)
    variants = _section_variants(rl, styles, result, scenes)
    disagreements = _section_disagreements(rl, styles, result)
    qc = _section_qc(rl, styles, result)
    methods = _section_methods(rl, styles, result, cohort)
    cohort_section = (_section_cohort(rl, styles, cohort, scenes)
                      if cohort is not None and len(cohort.samples) > 1 else [])

    story: List[Any] = []
    if profile == "research":
        story.extend(identity)
        story.append(PageBreak())
        story.extend(variants)
        story.append(PageBreak())
        story.extend(qc)
        story.append(PageBreak())
        story.extend(disagreements)
        story.extend(drugs)
        story.append(PageBreak())
        story.extend(typing)
    else:
        story.extend(identity)
        story.extend(drugs)
        story.append(PageBreak())
        story.extend(typing)
        story.append(PageBreak())
        story.extend(variants)
        story.append(PageBreak())
        story.extend(disagreements)
        story.extend(qc)
    if cohort_section:
        story.append(PageBreak())
        story.extend(cohort_section)
    story.append(PageBreak())
    story.extend(methods)
    return story


def _page_furniture(title: str, footer: str):
    """A running header and footer, drawn on every page."""
    def draw(canvas: Any, doc: Any) -> None:
        from reportlab.lib.colors import HexColor
        canvas.saveState()
        width, height = doc.pagesize
        canvas.setFont("Helvetica", 6.5)
        canvas.setFillColor(HexColor(MUTED))
        canvas.drawString(36, height - 26, title)
        canvas.drawRightString(width - 36, height - 26, "page {0}".format(doc.page))
        canvas.setStrokeColor(HexColor(LINE))
        canvas.setLineWidth(0.4)
        canvas.line(36, height - 32, width - 36, height - 32)
        canvas.line(36, 30, width - 36, 30)
        canvas.drawString(36, 21, footer)
        canvas.restoreState()
    return draw


def write_pdf(path: Any, result: SampleResult, cohort: Optional[CohortResult] = None,
              profile: str = "") -> Path:
    """Write the single-sample report and hand back the path it wrote.

    *profile* defaults to the profile recorded on the result, so a run does not
    have to say twice which report it asked for.
    """
    rl = _require_reportlab()
    chosen = profile or result.profile or "clinical"
    if chosen not in PROFILES:
        raise MjolnirError("unknown report profile {0!r}; expected one of {1}".format(
            chosen, ", ".join(PROFILES)))
    styles = _styles(rl)
    target = Path(path)
    ensure_dir(target.parent)
    document = rl["SimpleDocTemplate"](
        str(target), pagesize=rl["A4"], leftMargin=36, rightMargin=36,
        topMargin=44, bottomMargin=40,
        title="Mjolnir report - {0}".format(result.sample_id),
        author="Mjolnir {0}".format(result.mjolnir_version or __version__),
        subject="Mycobacterial resistance, lineage, species and contamination report")
    furniture = _page_furniture(
        "Mjolnir {0} | {1} | {2}".format(result.mjolnir_version or __version__,
                                         result.sample_id, result.platform),
        "Verdicts are rule-derived from the thresholds in Annex C; prose is labelled "
        "with its source. Absence of a call is absence, never susceptibility.")
    document.build(_story(rl, styles, result, cohort, chosen),
                   onFirstPage=furniture, onLaterPages=furniture)
    LOG.info("wrote PDF report: %s", target)
    return target


def write_cohort_pdf(path: Any, cohort: CohortResult,
                     results: Sequence[SampleResult] = ()) -> Path:
    """Write the cohort report: the drug grid across samples, distances, clusters."""
    rl = _require_reportlab()
    styles = _styles(rl)
    target = Path(path)
    ensure_dir(target.parent)
    scenes = dict(build_cohort_scenes(cohort, results))
    Paragraph, Spacer = rl["Paragraph"], rl["Spacer"]
    flowable = _scene_flowable_class(rl)
    headline, provenance = T.cohort_headline(cohort)

    story: List[Any] = [
        Paragraph("Mjolnir cohort report", styles["title"]),
        Paragraph("{0} samples | Mjolnir {1}".format(
            len(cohort.samples), _xml(cohort.mjolnir_version or __version__)),
            styles["subtitle"]),
        Spacer(1, 8),
        Paragraph(_xml(headline), styles["headline"]),
        Paragraph("Source of this sentence: {0}.".format(_xml(provenance)), styles["small"]),
        Spacer(1, 8),
    ]
    grid = scenes.get("cohort-drug-grid")
    if grid is not None:
        story.append(Paragraph("Resistance across the cohort", styles["h2"]))
        story.append(flowable(grid))
        story.append(Paragraph(_xml(grid.caption), styles["small"]))
        story.append(rl["PageBreak"]())
    story.extend(_section_cohort(rl, styles, cohort, scenes))
    if results:
        story.append(rl["PageBreak"]())
        story.append(Paragraph("Samples", styles["h2"]))
        story.extend(_row_table(rl, styles, [r.summary_row() for r in results], (
            ("sample", "sample", 86.0), ("platform", "platform", 44.0),
            ("species", "species", 108.0), ("lineage", "lineage", 72.0),
            ("mean_depth", "depth", 36.0), ("sample_validity", "validity", 52.0),
            ("mixture_class", "mixture", 60.0), ("disagreements", "disagreements", 65.0),
        )))
    story.append(rl["PageBreak"]())
    story.extend(_section_methods(rl, styles, None, cohort))

    document = rl["SimpleDocTemplate"](
        str(target), pagesize=rl["A4"], leftMargin=36, rightMargin=36,
        topMargin=44, bottomMargin=40, title="Mjolnir cohort report")
    furniture = _page_furniture(
        "Mjolnir {0} | cohort of {1}".format(cohort.mjolnir_version or __version__,
                                             len(cohort.samples)),
        "Distances are masked and carry their shared-callable denominator. An "
        "uncompared pair is absent, never zero.")
    document.build(story, onFirstPage=furniture, onLaterPages=furniture)
    LOG.info("wrote cohort PDF report: %s", target)
    return target


# ---------------------------------------------------------------------------
# Standalone figures, through matplotlib
# ---------------------------------------------------------------------------

def _require_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot
        from matplotlib.patches import Circle as MplCircle
        from matplotlib.patches import Polygon as MplPolygon
        from matplotlib.patches import Rectangle as MplRectangle
    except ImportError as exc:
        raise MjolnirError(
            "matplotlib is required to write standalone figures but is not installed "
            "({0}).\n  conda install -c conda-forge matplotlib\n"
            "  (or: pip install 'mjolnir-myco[report]')".format(exc)
        )
    return pyplot, MplRectangle, MplCircle, MplPolygon


def scene_to_matplotlib(scene: Scene, path: Any, dpi: int = 300) -> Path:
    """Replay a scene through matplotlib to a standalone figure file.

    The same objects the PDF draws, on a matplotlib canvas, so a figure lifted
    into a manuscript or a slide is the same figure the report shows. ``.svg``
    and ``.pdf`` come out as true vectors; ``.png`` is there for the places that
    still refuse anything else.
    """
    pyplot, MplRectangle, MplCircle, MplPolygon = _require_matplotlib()
    target = Path(path)
    ensure_dir(target.parent)
    figure = pyplot.figure(figsize=(scene.width / 72.0, scene.height / 72.0), dpi=dpi)
    axes = figure.add_axes((0.0, 0.0, 1.0, 1.0))
    axes.set_xlim(0, scene.width)
    # Inverted y so the scene's top-left origin needs no conversion at all.
    axes.set_ylim(scene.height, 0)
    axes.axis("off")
    for item in scene.items:
        if isinstance(item, Rect):
            axes.add_patch(MplRectangle(
                (item.x, item.y), item.w, item.h,
                facecolor=item.fill or "none", edgecolor=item.stroke or "none",
                linewidth=item.stroke_width,
                linestyle=(0, item.dash) if item.dash else "solid"))
        elif isinstance(item, Line):
            axes.plot([item.x1, item.x2], [item.y1, item.y2], color=item.stroke,
                      linewidth=item.width,
                      linestyle=(0, item.dash) if item.dash else "solid",
                      solid_capstyle="butt")
        elif isinstance(item, Circle):
            axes.add_patch(MplCircle((item.cx, item.cy), item.r,
                                     facecolor=item.fill or "none",
                                     edgecolor=item.stroke or "none",
                                     linewidth=item.stroke_width))
        elif isinstance(item, Poly):
            axes.add_patch(MplPolygon(item.points, closed=item.close,
                                      facecolor=item.fill or "none",
                                      edgecolor=item.stroke or "none",
                                      linewidth=item.stroke_width))
        elif isinstance(item, Text):
            axes.text(item.x, item.y, item.text, fontsize=item.size, color=item.fill,
                      ha={"start": "left", "middle": "center", "end": "right"}[item.anchor],
                      va="baseline", rotation=item.rotate,
                      fontweight="bold" if item.bold else "normal")
    figure.savefig(str(target), dpi=dpi, facecolor=PAPER)
    pyplot.close(figure)
    LOG.debug("wrote figure %s", target)
    return target


def write_figures(out_dir: Any, result: SampleResult,
                  cohort: Optional[CohortResult] = None,
                  suffix: str = "svg", dpi: int = 300) -> List[Path]:
    """Write every figure as a standalone file, for slides and manuscripts."""
    directory = ensure_dir(out_dir)
    written: List[Path] = []
    for name, scene in build_scenes(result, cohort):
        written.append(scene_to_matplotlib(
            scene, directory / "{0}.{1}.{2}".format(result.sample_id, name, suffix), dpi))
    LOG.info("wrote %d figure(s) to %s", len(written), directory)
    return written


def generated_stamp() -> str:
    """A human timestamp for the artefacts that carry one."""
    return datetime.now().strftime("%Y-%m-%d %H:%M")
