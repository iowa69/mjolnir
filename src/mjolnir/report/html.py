"""The same report as a single HTML file, with no network at view time.

Everything is inlined — stylesheet, script, and every figure as SVG written out
element by element. A report has to open on a laptop with no network, in a
hospital, years from now, and say exactly what it said the day it was written.
A CDN link is a promise somebody else has to keep.

The content is the PDF's content. Both renderers read the same rows from
:mod:`mjolnir.report.tables` and draw the same :class:`~mjolnir.report.pdf.Scene`
objects, so there is no second implementation of a drug grid to drift out of
step with the first — this file is the SVG transcription of a scene, and
nothing more.

Two deliberate departures from the PDF, both for the screen. The page chrome
follows the reader's light or dark preference, while every figure keeps a white
ground: a chart is a printed insert, and re-tinting a clinical colour scale to
suit a dark theme would change what the colours mean. And the long annex tables
are sortable and filterable, because a variant table is read by searching it.
"""

from __future__ import annotations

import html as html_escape
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .. import __version__
from ..config import CATALOGUES, MIN_SHARED_CALLABLE_SITES, PROFILES, source_for
from ..records import CohortResult, SampleResult, STATUS_FAIL, STATUS_PASS, STATUS_WARN
from ..utils import LOG, MjolnirError, ensure_dir
from . import tables as T
from .pdf import (
    Circle,
    Line,
    NOT_MEASURED_STYLE,
    Poly,
    Rect,
    STATUS_STYLE,
    Scene,
    Text,
    build_cohort_scenes,
    build_scenes,
    call_style,
    generated_stamp,
)

_CSS = """
:root{
  --bg:#ffffff; --fg:#16181d; --muted:#5f6672; --line:#e0e4ea; --panel:#f5f6f9;
  --accent:#1f4e9c; --paper:#ffffff; --shadow:rgba(20,24,32,.08);
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#111318; --fg:#e8eaee; --muted:#9aa3b2; --line:#272b34; --panel:#181b22;
    --accent:#7aa6e8; --shadow:rgba(0,0,0,.4);
  }
}
:root[data-theme="dark"]{
  --bg:#111318; --fg:#e8eaee; --muted:#9aa3b2; --line:#272b34; --panel:#181b22;
  --accent:#7aa6e8; --shadow:rgba(0,0,0,.4);
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font:14px/1.55 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:26px 20px 80px}
header.top{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap}
h1{font-size:25px;margin:0;letter-spacing:-.02em}
h2{font-size:16px;margin:30px 0 8px;letter-spacing:-.01em;
  border-bottom:1px solid var(--line);padding-bottom:5px}
h3{font-size:11px;margin:18px 0 6px;color:var(--muted);font-weight:700;
  text-transform:uppercase;letter-spacing:.07em}
.sub{color:var(--muted);font-size:12.5px}
.verdict{display:flex;gap:14px;align-items:center;border:1px solid;border-radius:10px;
  padding:12px 16px;margin:16px 0 6px}
.verdict .big{font-size:22px;font-weight:700;letter-spacing:-.02em}
.verdict .why{font-size:12.5px}
.headline{font-size:17px;line-height:1.45;font-weight:600;margin:14px 0 4px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(132px,1fr));gap:10px;margin:16px 0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:10px 12px}
.card .n{font-size:21px;font-weight:650;letter-spacing:-.02em}
.card .l{color:var(--muted);font-size:11.5px;margin-top:2px}
dl.pairs{display:grid;grid-template-columns:190px 1fr;gap:0;margin:6px 0 0;
  border-top:1px solid var(--line)}
dl.pairs dt{font-weight:600;padding:5px 8px 5px 0;border-bottom:1px solid var(--line);font-size:12.5px}
dl.pairs dd{margin:0;padding:5px 0;border-bottom:1px solid var(--line);font-size:12.5px;
  color:var(--fg)}
.scroll{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--panel)}
table{border-collapse:collapse;width:100%;font-size:12px}
th,td{padding:5px 8px;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}
thead th{position:sticky;top:0;background:var(--panel);z-index:2;font-weight:650;
  border-bottom:2px solid var(--line);cursor:pointer;user-select:none;white-space:nowrap}
thead th:hover{color:var(--accent)}
tbody tr:hover{background:var(--shadow)}
td.nowrap{white-space:nowrap}
tr.flagged td:first-child{box-shadow:inset 3px 0 0 var(--accent)}
.callcell{display:inline-block;min-width:26px;text-align:center;padding:1px 6px;
  border-radius:4px;font-weight:700;font-size:11.5px;border:1px solid}
.callcell.hatched{background-image:repeating-linear-gradient(45deg,
  rgba(0,0,0,.32) 0 1px,transparent 1px 4px)}
.calltext{font-size:12px}
.pill{display:inline-block;padding:1px 7px;border-radius:999px;font-size:11px;font-weight:650;
  border:1px solid}
.na{color:var(--muted);font-style:italic}
.figure{background:var(--paper);border:1px solid var(--line);border-radius:10px;
  padding:12px;margin:10px 0;overflow-x:auto}
.figure svg{display:block;max-width:100%;height:auto}
.figcap{color:var(--muted);font-size:11.5px;margin-top:6px}
.legend{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:4px 16px;
  margin:8px 0;font-size:12px}
.controls{display:flex;gap:8px;margin:8px 0;flex-wrap:wrap}
input[type=search]{padding:5px 9px;border:1px solid var(--line);border-radius:8px;
  background:var(--bg);color:var(--fg);min-width:230px;font-size:12.5px}
button.tgl{padding:5px 11px;border:1px solid var(--line);border-radius:8px;background:var(--bg);
  color:var(--fg);cursor:pointer;font-size:12.5px}
button.tgl:hover{border-color:var(--accent);color:var(--accent)}
ul.caveats{margin:6px 0;padding-left:18px;font-size:12.5px}
ul.caveats li{margin:3px 0}
footer{margin-top:44px;color:var(--muted);font-size:11.5px;border-top:1px solid var(--line);
  padding-top:12px}
details{margin:8px 0}
summary{cursor:pointer;color:var(--accent);font-size:12.5px}
pre{white-space:pre-wrap;word-break:break-word;font-size:11px;color:var(--muted)}
.empty{color:var(--muted);padding:12px;font-style:italic;font-size:12.5px}
@media print{
  body{background:#fff;color:#000}
  .controls,button.tgl{display:none}
  h2{page-break-after:avoid}
  .figure,table{page-break-inside:avoid}
}
"""

_JS = """
(function(){
  var root=document.documentElement;
  var btn=document.getElementById('theme');
  if(btn){btn.addEventListener('click',function(){
    var now=root.getAttribute('data-theme');
    var next=now==='dark'?'light':'dark';
    root.setAttribute('data-theme',next);
    btn.textContent=next==='dark'?'light theme':'dark theme';
  });}
  document.querySelectorAll('table.sortable thead th').forEach(function(th){
    th.addEventListener('click',function(){
      // cellIndex is the column within this table; a document-wide counter would
      // index past the end of every table after the first.
      var idx=th.cellIndex, table=th.closest('table');
      if(!table||!table.tBodies.length) return;
      var tb=table.tBodies[0], rows=Array.prototype.slice.call(tb.rows);
      var asc=!(th.dataset.asc==='1');
      th.closest('tr').querySelectorAll('th').forEach(function(o){o.dataset.asc='';});
      th.dataset.asc=asc?'1':'0';
      rows.sort(function(a,b){
        var x=a.cells[idx]?a.cells[idx].innerText.trim():'';
        var y=b.cells[idx]?b.cells[idx].innerText.trim():'';
        // NA sorts last in both directions: it is an absence, not a small number.
        if(x==='NA'&&y!=='NA') return 1;
        if(y==='NA'&&x!=='NA') return -1;
        var nx=parseFloat(x),ny=parseFloat(y);
        var both=!isNaN(nx)&&!isNaN(ny)&&x!==''&&y!=='';
        var c=both?(nx-ny):x.localeCompare(y,undefined,{numeric:true});
        return asc?c:-c;
      });
      rows.forEach(function(r){tb.appendChild(r);});
    });
  });
  document.querySelectorAll('input[data-filter]').forEach(function(inp){
    var table=document.getElementById(inp.dataset.filter);
    if(!table||!table.tBodies.length){inp.disabled=true;return;}
    inp.addEventListener('input',function(){
      var q=inp.value.toLowerCase();
      Array.prototype.slice.call(table.tBodies[0].rows).forEach(function(r){
        r.style.display=r.innerText.toLowerCase().indexOf(q)>-1?'':'none';
      });
    });
  });
})();
"""


def _e(value: Any) -> str:
    return html_escape.escape("" if value is None else str(value), quote=True)


# ---------------------------------------------------------------------------
# Scene -> SVG
# ---------------------------------------------------------------------------

def scene_to_svg(scene: Scene, title: str = "") -> str:
    """Render a scene as inline SVG.

    A white ground rectangle is drawn first and deliberately: the figures encode
    clinical meaning in colour, and letting a dark page theme show through would
    change every fill's contrast and with it the greyscale ordering the palette
    was built around.
    """
    parts: List[str] = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {0:g} {1:g}" '
        'width="{0:g}" height="{1:g}" role="img" '
        'font-family="ui-sans-serif,system-ui,Helvetica,Arial,sans-serif">'.format(
            scene.width, scene.height)
    ]
    if title or scene.title:
        parts.append("<title>{0}</title>".format(_e(title or scene.title)))
    parts.append('<rect x="0" y="0" width="{0:g}" height="{1:g}" fill="#ffffff"/>'.format(
        scene.width, scene.height))
    for item in scene.items:
        parts.append(_svg_item(item))
    parts.append("</svg>")
    return "".join(parts)


def _paint(fill: Optional[str], stroke: Optional[str], width: float,
           dash: Optional[Tuple[float, float]] = None) -> str:
    out = ' fill="{0}"'.format(fill if fill else "none")
    if stroke:
        out += ' stroke="{0}" stroke-width="{1:g}"'.format(stroke, width)
    else:
        out += ' stroke="none"'
    if dash:
        out += ' stroke-dasharray="{0}"'.format(" ".join("{0:g}".format(d) for d in dash))
    return out


def _svg_item(item: Any) -> str:
    if isinstance(item, Rect):
        body = '<rect x="{0:g}" y="{1:g}" width="{2:g}" height="{3:g}"{4}{5}>'.format(
            item.x, item.y, max(0.0, item.w), max(0.0, item.h),
            ' rx="{0:g}"'.format(item.radius) if item.radius else "",
            _paint(item.fill, item.stroke, item.stroke_width, item.dash))
        if item.title:
            body += "<title>{0}</title>".format(_e(item.title))
        return body + "</rect>"
    if isinstance(item, Line):
        return '<line x1="{0:g}" y1="{1:g}" x2="{2:g}" y2="{3:g}"{4}/>'.format(
            item.x1, item.y1, item.x2, item.y2,
            _paint(None, item.stroke, item.width, item.dash))
    if isinstance(item, Circle):
        body = '<circle cx="{0:g}" cy="{1:g}" r="{2:g}"{3}>'.format(
            item.cx, item.cy, item.r,
            _paint(item.fill, item.stroke, item.stroke_width))
        if item.title:
            body += "<title>{0}</title>".format(_e(item.title))
        return body + "</circle>"
    if isinstance(item, Poly):
        points = " ".join("{0:g},{1:g}".format(x, y) for x, y in item.points)
        tag = "polygon" if item.close else "polyline"
        return '<{0} points="{1}"{2}/>'.format(
            tag, points, _paint(item.fill, item.stroke, item.stroke_width))
    if isinstance(item, Text):
        # SVG rotates clockwise; the scene's angles are counter-clockwise, as in
        # reportlab, so the sign flips here and only here.
        transform = ('' if not item.rotate else
                     ' transform="rotate({0:g} {1:g} {2:g})"'.format(
                         -item.rotate, item.x, item.y))
        return ('<text x="{0:g}" y="{1:g}" font-size="{2:g}" fill="{3}" '
                'text-anchor="{4}"{5}{6}>{7}</text>').format(
            item.x, item.y, item.size, item.fill, item.anchor,
            ' font-weight="700"' if item.bold else "", transform, _e(item.text))
    return ""


def _figure(scene: Scene) -> str:
    if scene is None:
        return ""
    caption = ('<div class="figcap">{0}</div>'.format(_e(scene.caption))
               if scene.caption else "")
    return '<div class="figure">{0}{1}</div>'.format(scene_to_svg(scene), caption)


# ---------------------------------------------------------------------------
# HTML fragments
# ---------------------------------------------------------------------------

def _call_cell(call: str, glyph: str, title: str = "") -> str:
    style = call_style(call)
    classes = "callcell" + (" hatched" if style["hatch"] == "diag" else "")
    return ('<span class="{0}" style="background:{1};color:{2};border-color:{3}" '
            'title="{4}">{5}</span>').format(
        classes, style["fill"], style["ink"], style["border"],
        _e(title or T.call_label(call)), _e(glyph))


def _status_pill(status: str) -> str:
    if status == "not measured":
        palette = NOT_MEASURED_STYLE
    else:
        palette = STATUS_STYLE.get(status)
    if palette is None:
        return _e(status)
    return ('<span class="pill" style="background:{0};color:{1};border-color:{2}">'
            "{3}</span>").format(palette["fill"], palette["ink"], palette["border"],
                                 _e(status))


def _pairs(pairs: Sequence[Tuple[str, str]]) -> str:
    if not pairs:
        return '<div class="empty">nothing recorded</div>'
    body = "".join("<dt>{0}</dt><dd>{1}</dd>".format(_e(label), _e(value))
                   for label, value in pairs)
    return '<dl class="pairs">{0}</dl>'.format(body)


def _cell_html(value: Any, column: str = "", status_column: str = "") -> str:
    if column and column == status_column:
        return _status_pill(str(value or ""))
    if value is None:
        return '<span class="na">{0}</span>'.format(T.TSV_NA)
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        if value.is_integer():
            return _e(int(value))
        return _e("{0:.4f}".format(value).rstrip("0").rstrip("."))
    if isinstance(value, (list, tuple)):
        return _e("; ".join(str(v) for v in value)) if value else ""
    return _e(value)


def _table(rows: Sequence[Dict[str, Any]], columns: Sequence[Tuple[str, str]],
           table_id: str, status_column: str = "", flag_key: str = "",
           searchable: bool = False, empty: str = "nothing to show in this table") -> str:
    """A sortable annex table. Every table in the report goes through here."""
    if not rows:
        return '<div class="empty">{0}</div>'.format(_e(empty))
    out: List[str] = []
    if searchable:
        out.append('<div class="controls"><input type="search" data-filter="{0}" '
                   'placeholder="filter rows..."></div>'.format(_e(table_id)))
    out.append('<div class="scroll"><table id="{0}" class="sortable"><thead><tr>'.format(
        _e(table_id)))
    for _key, title in columns:
        out.append("<th>{0}</th>".format(_e(title)))
    out.append("</tr></thead><tbody>")
    for row in rows:
        flagged = bool(flag_key and row.get(flag_key))
        out.append('<tr class="flagged">' if flagged else "<tr>")
        for key, _title in columns:
            out.append("<td>{0}</td>".format(_cell_html(row.get(key), key, status_column)))
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def _card(value: Any, label: str) -> str:
    return '<div class="card"><div class="n">{0}</div><div class="l">{1}</div></div>'.format(
        _e(value), _e(label))


def _caveat_list(lines: Sequence[str]) -> str:
    if not lines:
        return ""
    return '<ul class="caveats">{0}</ul>'.format(
        "".join("<li>{0}</li>".format(_e(line)) for line in lines))


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

_VALIDITY_PALETTE = {
    "valid": STATUS_STYLE[STATUS_PASS],
    "suspect": STATUS_STYLE[STATUS_WARN],
    "invalid": STATUS_STYLE[STATUS_FAIL],
    "not-assessed": NOT_MEASURED_STYLE,
}


def _section_identity(result: SampleResult) -> str:
    headline, provenance = T.headline_sentence(result)
    palette = _VALIDITY_PALETTE.get(result.contamination.verdict, STATUS_STYLE[STATUS_WARN])
    resistant = result.resistant_drugs()
    no_determinant = [d for d in result.drugs if d.call == "no-call"]
    cards = [
        _card(len(resistant), "drugs with a determinant"),
        _card(len(no_determinant), "drugs with no determinant"),
        _card(len(result.disagreements()), "catalogue disagreements"),
        _card(len(result.variants), "variants called"),
        _card(T.fmt_number(result.qc.mean_depth, 1, "x", na="NA"), "mean depth"),
        _card(_e(result.lineage.display), "lineage"),
        _card(len(result.unmeasured()), "metrics not measured"),
    ]
    return "".join([
        '<div class="verdict" style="background:{0};border-color:{1};color:{2}">'.format(
            palette["fill"], palette["border"], palette["ink"]),
        '<div><div class="sub" style="color:{0}">Sample validity</div>'.format(palette["ink"]),
        '<div class="big">{0}</div></div>'.format(_e(result.contamination.verdict)),
        '<div class="why">{0}</div></div>'.format(
            _e(result.contamination.verdict_reason
               or "no reason for this verdict was recorded")),
        '<div class="headline">{0}</div>'.format(_e(headline)),
        '<div class="sub">Source of this sentence: {0}. Pass, warn and fail are computed '
        "in Python from the thresholds in the QC annex before any model is called."
        "</div>".format(_e(provenance)),
        '<div class="cards">{0}</div>'.format("".join(cards)),
        "<h2>Sample</h2>", _pairs(T.identity_pairs(result)),
        "<h2>Species</h2>", _pairs(T.species_pairs(result)),
    ])


def _section_drugs(result: SampleResult, scenes: Dict[str, Scene]) -> str:
    rows = T.drug_rows(result)
    calls = sorted(result.drugs, key=lambda d: T.drug_order(d.drug))
    out = ["<h2>Drug-by-drug resistance</h2>"]
    if not rows:
        out.append('<div class="empty">no drug was evaluated for this sample</div>')
    else:
        out.append('<div class="scroll"><table id="tbl-drugs" class="sortable"><thead><tr>')
        headers = ["Drug", "Determinant", "Mjolnir call", "WHO grade"]
        headers += [c.replace("WHO v2", "WHO") for c in CATALOGUES]
        headers += ["Flags"]
        out.extend("<th>{0}</th>".format(_e(h)) for h in headers)
        out.append("</tr></thead><tbody>")
        for row, call in zip(rows, calls):
            flags = T.drug_flags(call)
            out.append('<tr class="flagged">' if row["disagreement"] else "<tr>")
            out.append('<td class="nowrap"><b>{0}</b><br><span class="sub">{1}</span></td>'
                       .format(_e(row["drug_code"]), _e(row["drug"])))
            out.append("<td>{0}</td>".format(
                _e("; ".join(row["supporting_variants"]) or "none matched")))
            out.append('<td class="calltext">{0} {1}</td>'.format(
                _call_cell(row["call"], row["call_glyph"], row["call_label"]),
                _e(row["call_label"])))
            out.append("<td>{0}</td>".format(_e(row["who_grade"] or "not graded by WHO")))
            for catalogue in CATALOGUES:
                key = T.catalogue_key(catalogue)
                cat_call = row[key + "_call"]
                out.append('<td class="nowrap">{0}</td>'.format(_call_cell(
                    cat_call, T.CALL_GLYPH.get(cat_call, "?"),
                    "{0}: {1}{2}".format(catalogue, T.call_label(cat_call),
                                         " [{0}]".format(row[key + "_grade"])
                                         if row[key + "_grade"] else ""))))
            out.append("<td>{0}</td>".format(_e("; ".join(flags))))
            out.append("</tr>")
        out.append("</tbody></table></div>")
    out.append(_figure(scenes.get("drug-grid")))
    out.append("<h3>How to read a call</h3>")
    out.append('<div class="legend">')
    for call, text in T.CALL_LEGEND:
        out.append("<div>{0} {1}</div>".format(
            _call_cell(call, T.CALL_GLYPH.get(call, "?")), _e(text)))
    out.append("</div>")
    caveats = T.caveat_lines(result)
    if caveats:
        out.append("<h3>Platform and sample caveats</h3>")
        out.append(_caveat_list(caveats))
    return "".join(out)


def _section_typing(result: SampleResult, scenes: Dict[str, Scene]) -> str:
    return "".join([
        "<h2>Lineage and barcode evidence</h2>", _pairs(T.lineage_pairs(result)),
        "<h2>Contamination and sample validity</h2>", _pairs(T.validity_pairs(result)),
        "<h3>Contamination metrics against their thresholds</h3>",
        _table(T.contamination_rows(result),
               [("check", "metric"), ("value", "value"), ("comparison", "vs"),
                ("threshold", "threshold"), ("status", "status"), ("source", "source"),
                ("reading", "reading")],
               "tbl-contamination", status_column="status"),
        "<h2>Coverage and depth</h2>", _figure(scenes.get("coverage")),
    ])


def _section_variants(result: SampleResult, scenes: Dict[str, Scene]) -> str:
    columns = [("coordinate", "coordinate"), ("gene", "gene"), ("hgvs", "HGVS"),
               ("variant_type", "type"), ("effect", "effect"), ("depth", "depth"),
               ("alt_reads", "alt reads"), ("allele_fraction", "allele fraction"),
               ("is_major", "major"), ("masked", "masked"), ("filters", "filters")]
    for catalogue in CATALOGUES:
        columns.append((T.catalogue_key(catalogue) + "_grade",
                        "{0} grade".format(catalogue.replace("WHO v2", "WHO"))))
    return "".join([
        "<h2>Annex A - variants</h2>", _figure(scenes.get("allele-fraction")),
        _table(T.variant_rows(result), columns, "tbl-variants", searchable=True,
               empty="no variants were called for this sample"),
    ])


def _section_disagreements(result: SampleResult) -> str:
    rows = T.disagreement_rows(result)
    if not rows:
        return ("<h2>Annex B - catalogue disagreements</h2>"
                '<div class="empty">No drug in this sample carries a catalogue '
                "disagreement. That is a statement about the three catalogues agreeing, "
                "not about the organism.</div>")
    return "".join([
        "<h2>Annex B - catalogue disagreements</h2>",
        _table(rows, [("drug", "drug"), ("mjolnir_call", "Mjolnir"),
                      ("catalogue", "catalogue"), ("catalogue_call", "its call"),
                      ("grade", "grade"), ("disagreement_kind", "kind"),
                      ("evidence", "evidence"), ("note", "note")],
               "tbl-disagreements"),
    ])


def _section_qc(result: SampleResult) -> str:
    out = [
        "<h2>Annex C - QC metrics against their thresholds</h2>",
        '<div class="sub">Every row states the threshold that was applied and the '
        'document it came from. A row marked "not measured" is a gap in the evidence '
        "and is never a pass.</div>",
        _table(T.check_rows(T.all_checks(result), result.sample_id),
               [("category", "panel"), ("check", "check"), ("value", "value"),
                ("comparison", "vs"), ("threshold", "threshold"), ("unit", "unit"),
                ("status", "status"), ("source", "source"), ("reading", "reading")],
               "tbl-checks", status_column="status", searchable=True),
        "<h3>Measured without a registered threshold</h3>",
        '<div class="sub">Reported without a pass or fail because Mjolnir registers no '
        "published bound for them, and a status implies a bound.</div>",
        _table(T.observation_rows(result),
               [("panel", "panel"), ("observation", "observation"), ("value", "value"),
                ("unit", "unit"), ("note", "note")], "tbl-observations"),
    ]
    support = T.lineage_support_rows(result)
    if support:
        keys = [k for k in support[0] if k != "sample"]
        out.append("<h3>Barcode sites</h3>")
        out.append(_table(support, [(k, k.replace("_", " ")) for k in keys],
                          "tbl-barcode", searchable=True))
    return "".join(out)


def _section_cohort(cohort: CohortResult, scenes: Dict[str, Scene]) -> str:
    headline, provenance = T.cohort_headline(cohort)
    return "".join([
        "<h2>Annex D - cohort distances and clusters</h2>",
        '<div class="headline">{0}</div>'.format(_e(headline)),
        '<div class="sub">Source of this sentence: {0}.</div>'.format(_e(provenance)),
        _pairs(T.cohort_pairs(cohort)),
        _figure(scenes.get("distance-matrix")),
        _figure(scenes.get("dendrogram")),
        "<h3>Clusters</h3>",
        _table(T.cluster_rows(cohort),
               [("cluster", "cluster"), ("size", "n"), ("members", "members"),
                ("threshold", "threshold"), ("max_distance", "max SNPs"),
                ("min_shared_callable_sites", "min shared callable sites"),
                ("threshold_basis", "threshold basis")], "tbl-clusters"),
        "<h3>Pairwise distances with their denominators</h3>",
        '<div class="sub">A distance is only comparable to the published SNP thresholds '
        "when the two samples share enough callable sequence; the floor is {0:,} bp "
        "({1}).</div>".format(MIN_SHARED_CALLABLE_SITES,
                              _e(source_for("min_shared_callable_sites"))),
        _table(T.distance_rows(cohort),
               [("sample_a", "sample A"), ("sample_b", "sample B"), ("snps", "SNPs"),
                ("shared_callable_sites", "shared callable"), ("snps_per_mb", "SNPs/Mb"),
                ("masked_sites", "masked"), ("within_threshold", "within threshold"),
                ("shared_sites_sufficient", "sites sufficient"),
                ("same_cluster", "same cluster"), ("note", "note")],
               "tbl-distances", searchable=True),
    ])


def _section_methods(result: Optional[SampleResult],
                     cohort: Optional[CohortResult]) -> str:
    tools = dict(result.tool_versions) if result is not None else {}
    databases = list(result.database_versions) if result is not None else []
    if cohort is not None:
        tools.update(cohort.tool_versions)
        databases = databases + list(cohort.database_versions)
    out = ["<h2>Annex E - methods, versions and checksums</h2>"]
    if result is not None:
        out.append(_pairs(T.methods_pairs(result)))
    out.extend([
        "<h3>Databases</h3>",
        '<div class="sub">A catalogue-version mismatch between two installations changes '
        "calls, so every database prints its version and checksum here.</div>",
        _table(T.database_rows(databases),
               [("database", "database"), ("version", "version"), ("checksum", "checksum"),
                ("licence", "licence"), ("citation", "citation"), ("url", "url"),
                ("fetched", "fetched")], "tbl-databases"),
        "<h3>Tools</h3>",
        _table(T.tool_version_rows(tools), [("tool", "tool"), ("version", "version")],
               "tbl-tools"),
        "<h3>Thresholds whose citation has not been checked</h3>",
        '<div class="sub">An unverified citation is worse than none, because it looks '
        "settled. These numbers are in use and nobody has yet opened the primary "
        "document on this machine.</div>",
        _table(T.unverified_rows(),
               [("threshold", "threshold"), ("value", "value"), ("source", "source"),
                ("note", "note")], "tbl-unverified"),
        "<details><summary>Every registered threshold and its source</summary>",
        _table(T.threshold_rows(),
               [("threshold", "threshold"), ("value", "value"), ("unit", "unit"),
                ("citation_verified", "verified"), ("source", "source"), ("note", "note")],
               "tbl-thresholds", searchable=True),
        "</details>",
    ])
    return "".join(out)


# ---------------------------------------------------------------------------
# Documents
# ---------------------------------------------------------------------------

def _document(title: str, subtitle: str, body: str, footer: str) -> str:
    return """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{css}</style></head>
<body><div class="wrap">
<header class="top"><h1>{h1}</h1><span class="sub">{sub}</span>
<button class="tgl" id="theme" type="button">dark theme</button></header>
{body}
<footer>{footer}</footer>
</div><script>{js}</script></body></html>""".format(
        title=_e(title), css=_CSS, h1=_e(title), sub=subtitle, body=body,
        footer=footer, js=_JS)


def render_html(result: SampleResult, cohort: Optional[CohortResult] = None,
                profile: str = "", generated: str = "") -> str:
    """The complete single-sample report as one self-contained HTML string."""
    chosen = profile or result.profile or "clinical"
    if chosen not in PROFILES:
        raise MjolnirError("unknown report profile {0!r}; expected one of {1}".format(
            chosen, ", ".join(PROFILES)))
    scenes = dict(build_scenes(result, cohort))
    identity = _section_identity(result)
    drugs = _section_drugs(result, scenes)
    typing = _section_typing(result, scenes)
    variants = _section_variants(result, scenes)
    disagreements = _section_disagreements(result)
    qc = _section_qc(result)
    methods = _section_methods(result, cohort)
    cohort_html = (_section_cohort(cohort, scenes)
                   if cohort is not None and len(cohort.samples) > 1 else "")

    if chosen == "research":
        body = identity + variants + qc + disagreements + drugs + typing
    else:
        body = identity + drugs + typing + variants + disagreements + qc
    body += cohort_html + methods

    subtitle = "{0} | {1} | Mjolnir {2} | profile {3}{4}".format(
        _e(result.sample_id), _e(result.platform),
        _e(result.mjolnir_version or __version__), _e(chosen),
        " | " + _e(generated) if generated else "")
    footer = ("Generated by Mjolnir {0}. Verdicts are rule-derived from the thresholds in "
              "Annex C; prose is labelled with its source. Absence of a call is absence, "
              "never susceptibility. This file is self-contained and needs no "
              "network.".format(_e(result.mjolnir_version or __version__)))
    return _document("Mjolnir report - {0}".format(result.sample_id), subtitle, body, footer)


def render_cohort_html(cohort: CohortResult, results: Sequence[SampleResult] = (),
                       generated: str = "") -> str:
    """The cohort report: the drug grid across samples, distances, clusters."""
    scenes = dict(build_cohort_scenes(cohort, results))
    headline, provenance = T.cohort_headline(cohort)
    parts = [
        '<div class="headline">{0}</div>'.format(_e(headline)),
        '<div class="sub">Source of this sentence: {0}.</div>'.format(_e(provenance)),
        '<div class="cards">{0}</div>'.format("".join([
            _card(len(cohort.samples), "samples"),
            _card(len(cohort.clusters), "clusters"),
            _card(cohort.threshold if cohort.threshold is not None else "NA",
                  "SNP threshold"),
            _card(T.fmt_number(cohort.masked_sites, na="NA"), "masked positions"),
            _card(len([p for p in cohort.pairs if p.snps is None]), "pairs not compared"),
        ])),
    ]
    if results:
        parts.append("<h2>Resistance across the cohort</h2>")
        parts.append(_figure(scenes.get("cohort-drug-grid")))
        parts.append('<div class="legend">')
        for call, text in T.CALL_LEGEND:
            parts.append("<div>{0} {1}</div>".format(
                _call_cell(call, T.CALL_GLYPH.get(call, "?")), _e(text)))
        parts.append("</div>")
        parts.append("<h2>Samples</h2>")
        parts.append(_table([r.summary_row() for r in results],
                            [("sample", "sample"), ("platform", "platform"),
                             ("species", "species"), ("lineage", "lineage"),
                             ("mean_depth", "mean depth"),
                             ("breadth_min_depth", "breadth"),
                             ("sample_validity", "validity"),
                             ("mixture_class", "mixture"),
                             ("n_variants", "variants"),
                             ("disagreements", "disagreements")],
                            "tbl-cohort-samples", searchable=True))
    parts.append(_section_cohort(cohort, scenes))
    parts.append(_section_methods(None, cohort))
    subtitle = "{0} samples | Mjolnir {1}{2}".format(
        len(cohort.samples), _e(cohort.mjolnir_version or __version__),
        " | " + _e(generated) if generated else "")
    footer = ("Generated by Mjolnir {0}. Distances are masked and carry their "
              "shared-callable denominator; an uncompared pair is absent, never "
              "zero.".format(_e(cohort.mjolnir_version or __version__)))
    return _document("Mjolnir cohort report", subtitle, "".join(parts), footer)


def write_html(path: Any, result: SampleResult, cohort: Optional[CohortResult] = None,
               profile: str = "", generated: str = "") -> Path:
    target = Path(path)
    ensure_dir(target.parent)
    target.write_text(render_html(result, cohort, profile,
                                  generated or generated_stamp()), encoding="utf-8")
    LOG.info("wrote HTML report: %s", target)
    return target


def write_cohort_html(path: Any, cohort: CohortResult,
                      results: Sequence[SampleResult] = (),
                      generated: str = "") -> Path:
    target = Path(path)
    ensure_dir(target.parent)
    target.write_text(render_cohort_html(cohort, results,
                                         generated or generated_stamp()),
                      encoding="utf-8")
    LOG.info("wrote cohort HTML report: %s", target)
    return target


def write_scene_svgs(out_dir: Any, result: SampleResult,
                     cohort: Optional[CohortResult] = None) -> List[Path]:
    """Each figure as a standalone SVG file, without needing matplotlib.

    :func:`mjolnir.report.pdf.write_figures` goes through matplotlib and can also
    emit PDF and PNG; this is the dependency-free path when SVG is all that is
    wanted.
    """
    directory = ensure_dir(out_dir)
    written: List[Path] = []
    for name, scene in build_scenes(result, cohort):
        target = directory / "{0}.{1}.svg".format(result.sample_id, name)
        target.write_text(scene_to_svg(scene), encoding="utf-8")
        written.append(target)
    LOG.info("wrote %d SVG figure(s) to %s", len(written), directory)
    return written


def render_json_block(payload: Any) -> str:
    """Embed a JSON payload in a collapsible block, escaped for HTML."""
    return "<details><summary>Machine-readable payload</summary><pre>{0}</pre></details>".format(
        _e(json.dumps(payload, indent=2, default=str)))
