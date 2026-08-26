## Engine-agnostic builder for FLIPPer's per-run "<file>_candidate_report.html" - one self-contained
## HTML report, one collapsible "card" per final candidate, combining everything about that
## candidate that's otherwise scattered across separate files: its full sequence and
## physicochemical properties, its repeat-region stats, its repeat copies aligned against each
## other, and an interactive disorder profile with each repeat copy plotted as a block against the
## same sequence coordinates - so a candidate can be assessed at a glance.
##
## Shared between FLIPPer's repeat-detection engines (XSTREAM, DetectRepeats) - each engine is
## responsible for producing the normalized inputs below (a DataFrame in the schema documented on
## build_candidate_report, a sequence lookup dict, and an alignment_provider callable), not for any
## of the HTML/CSS/JS in this module. That's what keeps the two engines' reports visually and
## behaviourally identical instead of two independently-maintained copies of the same template.

## amino acid background colors for the alignment view below - a standard physicochemical
## grouping (hydrophobic/aromatic/polar/negative/positive), the same kind of scheme alignment
## viewers like Clustal or DECIPHER's own BrowseSeqs use. Non-standard characters (e.g. DECIPHER's
## ConsensusSequence uses "+"/"." for weakly-conserved and gap-heavy columns, "X" for no
## consensus) aren't in this dict, so alignment_html_block renders them the same faint, unfilled
## way as a '-' gap (the '.weak' CSS class) instead of giving them a background colour - a solid
## coloured (or same-as-background) square either reads as "another residue type" or, for the old
## surface-2-on-dark-text combination, was invisible outright in dark mode.
REPEAT_AA_COLORS = {
    'A': '#8CFF8C', 'V': '#8CFF8C', 'L': '#8CFF8C', 'I': '#8CFF8C',
    'M': '#8CFF8C', 'P': '#8CFF8C', 'G': '#8CFF8C',
    'F': '#FFC800', 'W': '#FFC800', 'Y': '#FFC800',
    'S': '#79C7FF', 'T': '#79C7FF', 'N': '#79C7FF', 'Q': '#79C7FF', 'C': '#79C7FF',
    'D': '#FF7A7A', 'E': '#FF7A7A',
    'K': '#C9A6FF', 'R': '#C9A6FF', 'H': '#C9A6FF',
}

## copy-index colours for the repeat-unit alignment's dot indicators and their matching blocks
## in the disorder chart above - golden-angle hue stepping gives each consecutive index a visually
## distinct hue without needing to know the total copy count up front, and stays stable/consistent
## between the two views because both call this same function with the same (region-local) index
def _repeat_copy_color(i):
    hue = (i * 137.508) % 360
    return "hsl({:.0f}, 62%, 47%)".format(hue)

## Page shell for build_candidate_report(). Plain string with %%TOKEN%% placeholders (filled in
## by str.replace(), not str.format()) - the <script> below has many literal { } braces that would
## otherwise all need doubling to survive .format()'s escaping.
##
## Layout notes:
## - the disorder chart is a hand-rolled inline SVG (viewBox + width:100%), not a raster image -
##   it scales to fill the full card width with no leftover margin, and lets the hover tooltip
##   below trace an exact residue position/letter/score, neither of which a static plot supports.
## - each detected repeat region gets its own row of blocks under the disorder curve, positioned
##   in the same sequence coordinates as the curve, so it's visually obvious which stretch of the
##   disorder profile each repeat copy corresponds to.
_CANDIDATE_REPORT_TEMPLATE = """<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>FLIPPer Candidate Report</title>
<style>
:root {
  --bg: #f2f5f7;
  --surface: #ffffff;
  --surface-2: #eef2f4;
  --border: #d7dfe4;
  --text: #16212b;
  --text-muted: #5b6b76;
  --accent: #0f6e64;
  --accent-text: #ffffff;
  --accent-soft: #d8f0ec;
  --plddt: #c07f14;
  --shadow: rgba(20, 30, 35, 0.06);
  --font-sans: ui-sans-serif, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  --font-mono: ui-monospace, "SF Mono", "Cascadia Code", "Roboto Mono", Consolas, "Liberation Mono", monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #0e1620;
    --surface: #16212c;
    --surface-2: #1c2833;
    --border: #2b3946;
    --text: #e7edf1;
    --text-muted: #93a4b0;
    --accent: #34d0bd;
    --accent-text: #062420;
    --accent-soft: #113632;
    --plddt: #e8ab4a;
    --shadow: rgba(0, 0, 0, 0.4);
  }
}
:root[data-theme="dark"] {
  --bg: #0e1620;
  --surface: #16212c;
  --surface-2: #1c2833;
  --border: #2b3946;
  --text: #e7edf1;
  --text-muted: #93a4b0;
  --accent: #34d0bd;
  --accent-text: #062420;
  --accent-soft: #113632;
  --shadow: rgba(0, 0, 0, 0.4);
}
* { box-sizing: border-box; }
body {
  font-family: var(--font-sans);
  background: var(--bg);
  color: var(--text);
  margin: 0;
  padding: 2.5rem clamp(1rem, 4vw, 3rem);
}
.page-header { display: flex; flex-direction: column; gap: 0.35rem; margin: 0 auto 2rem; max-width: 1400px; width: 100%; }
.page-header-top { display: flex; align-items: baseline; justify-content: space-between; gap: 1rem; flex-wrap: wrap; }
.page-header h1 { font-size: 1.5rem; margin: 0; text-wrap: balance; letter-spacing: -0.01em; }
.page-header .summary { color: var(--text-muted); font-size: 0.9rem; }
.toggle-all-btn {
  font-family: var(--font-sans); font-size: 0.78rem; font-weight: 600;
  background: var(--surface-2); color: var(--text); border: 1px solid var(--border);
  border-radius: 8px; padding: 0.35rem 0.75rem; cursor: pointer; white-space: nowrap;
}
.toggle-all-btn:hover { background: var(--accent-soft); color: var(--accent); }
.page-header-actions { display: flex; align-items: center; gap: 0.5rem; }
.sort-select {
  font-family: var(--font-sans); font-size: 0.78rem; font-weight: 600;
  background: var(--surface-2); color: var(--text); border: 1px solid var(--border);
  border-radius: 8px; padding: 0.35rem 0.6rem; cursor: pointer;
}
.cards { display: flex; flex-direction: column; gap: 1.25rem; margin: 0 auto; max-width: 1400px; width: 100%; }
.card {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: 12px;
  padding: 1.4rem 1.6rem;
  box-shadow: 0 1px 3px var(--shadow);
}
.card-head { display: flex; flex-direction: column; gap: 0.7rem; margin-bottom: 0.9rem; cursor: pointer; list-style: none; }
.card-head::-webkit-details-marker { display: none; }
.card-head-row { display: flex; align-items: baseline; justify-content: space-between; gap: 1rem; }
.card-head h2 { margin: 0; font-family: var(--font-mono); font-size: 1.05rem; font-weight: 600; }
.card-head h2::before { content: '\\25b8'; display: inline-block; margin-right: 0.55rem; color: var(--text-muted); font-size: 0.7rem; transition: transform 0.15s ease; }
details.card[open] > .card-head h2::before { transform: rotate(90deg); }
/* only shown while the card is collapsed - lets a closed card still be scanned for its
   disorder shape and repeat-block layout without opening it, since <details> otherwise hides
   everything but the summary when closed */
.card-preview { display: none; }
details.card:not([open]) .card-preview { display: block; }
.preview-svg { display: block; width: 100%; height: 3.2rem; }
.preview-area { fill: var(--accent); opacity: 0.14; stroke: none; }
.preview-line { fill: none; stroke: var(--accent); stroke-width: 2; vector-effect: non-scaling-stroke; }
.preview-block { opacity: 0.8; stroke: var(--surface); stroke-width: 1; vector-effect: non-scaling-stroke; }
.pill { font-size: 0.75rem; font-weight: 600; padding: 0.2rem 0.6rem; border-radius: 999px; white-space: nowrap; font-variant-numeric: tabular-nums; cursor: help; }
.pill-high { background: var(--accent-soft); color: var(--accent); }
.pill-mid { background: var(--surface-2); color: var(--text-muted); }
.pill-low { background: var(--surface-2); color: var(--text-muted); }
h3 { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.06em; color: var(--text-muted); font-weight: 600; margin: 1.1rem 0 0.5rem; }
h3 .dim { text-transform: none; letter-spacing: normal; font-weight: 400; }
h4 { font-size: 0.78rem; color: var(--text-muted); font-weight: 600; margin: 0.7rem 0 0.35rem; }
.region-chips { display: flex; flex-wrap: wrap; gap: 0.5rem; }
.region-chip {
  display: flex; align-items: center; gap: 0.5rem;
  background: var(--surface-2); border: 1px solid var(--border);
  border-radius: 8px; padding: 0.3rem 0.7rem; font-size: 0.78rem; color: var(--text-muted);
  font-variant-numeric: tabular-nums;
}
.region-chip-label { font-weight: 600; color: var(--accent); }
.region-chip b { color: var(--text); font-weight: 600; }
pre.sequence {
  font-family: var(--font-mono); font-size: 0.78rem; line-height: 1.4;
  white-space: pre-wrap; word-break: break-all;
  background: var(--surface-2); border: 1px solid var(--border); border-radius: 8px;
  padding: 0.7rem 0.85rem; max-height: 7rem; overflow-y: auto; color: var(--text);
}
pre.sequence .seq-res { border-radius: 2px; color: #16212b; }
/* two-tone ring (surface, then text colour) instead of a single outline - a plain outline in
   one colour disappears against residue backgrounds close to that colour (pale greens/yellows
   in dark mode, the grey alignment gutter, etc); the inner surface-coloured ring guarantees
   separation from whatever colour is directly underneath, whatever that colour is */
.seq-res.active { outline: none; box-shadow: 0 0 0 1.5px var(--surface), 0 0 0 3.5px var(--text); }
.section-head { display: flex; align-items: baseline; justify-content: space-between; gap: 1rem; flex-wrap: wrap; margin: 1.1rem 0 0.5rem; }
.section-head h3 { margin: 0; }
.blast-btn {
  font-family: var(--font-sans); font-size: 0.78rem; font-weight: 600;
  background: var(--surface-2); color: var(--text); border: 1px solid var(--border);
  border-radius: 8px; padding: 0.3rem 0.7rem; white-space: nowrap; text-decoration: none;
}
.blast-btn:hover { background: var(--accent-soft); color: var(--accent); }
.chart-legend { display: flex; gap: 1.1rem; font-size: 0.72rem; color: var(--text-muted); margin: 0.15rem 0 0.5rem; }
.section-head .chart-legend { margin: 0; }
.chart-legend .legend-dot { display: inline-block; width: 0.55rem; height: 0.55rem; border-radius: 50%; margin-right: 0.35rem; vertical-align: middle; }
.chart-legend .legend-dot.disorder { background: var(--accent); }
.chart-legend .legend-dot.plddt { background: var(--plddt); }
.chart { position: relative; margin-top: 0.25rem; width: 100%; }
.chart-svg { display: block; width: 100%; height: auto; }
.chart-area { fill: var(--accent); opacity: 0.14; stroke: none; }
.chart-line { fill: none; stroke: var(--accent); stroke-width: 1.6; }
.chart-line-plddt { fill: none; stroke: var(--plddt); stroke-width: 1.4; stroke-dasharray: 4 2; opacity: 0.9; }
.chart-gridline { stroke: var(--border); stroke-width: 1; }
.chart-threshold { stroke: var(--text-muted); stroke-width: 1; stroke-dasharray: 3 2; }
.chart-axis-text { font-family: var(--font-mono); font-size: 7.5px; fill: var(--text-muted); }
.chart-block { opacity: 0.8; stroke: var(--surface); stroke-width: 1; }
.chart-guide { stroke: var(--text); stroke-width: 1; opacity: 0; pointer-events: none; }
.chart-dot { fill: var(--accent); stroke: var(--surface); stroke-width: 1.5; opacity: 0; pointer-events: none; }
.chart-dot-plddt { fill: var(--plddt); }
.chart-capture { fill: transparent; cursor: crosshair; }
.chart-tooltip {
  position: absolute; top: -0.35rem; transform: translate(-50%, -100%);
  background: var(--text); color: var(--surface);
  font-family: var(--font-mono); font-size: 0.72rem;
  padding: 0.28rem 0.55rem; border-radius: 6px; white-space: nowrap;
  opacity: 0; pointer-events: none; transition: opacity 0.06s ease;
  left: 0;
}
.alignment { overflow: auto; max-height: 14rem; margin-top: 0.6rem; font-family: var(--font-mono); font-size: 0.78rem; border: 1px solid var(--border); border-radius: 8px; padding: 0.55rem 0.7rem; background: transparent; }
.region-block + .region-block { margin-top: 1.1rem; }
.align-row { white-space: nowrap; line-height: 1.4; }
.align-name { display: inline-block; color: var(--text-muted); font-variant-numeric: tabular-nums; }
.align-seq span { display: inline-block; width: 0.85em; text-align: center; color: #16212b; border-radius: 2px; }
.align-seq span.weak { background: transparent !important; color: var(--text-muted); opacity: 0.55; font-weight: 300; }
.copy-dot { display: inline-block; width: 0.55rem; height: 0.55rem; border-radius: 50%; margin-right: 0.4rem; vertical-align: middle; box-shadow: 0 0 0 1px var(--surface); }
.align-row.consensus { font-weight: 700; border-top: 1px solid var(--border); margin-top: 0.2rem; padding-top: 0.2rem; }
.align-row.consensus .align-name { color: var(--text); }
.align-row.disorder-row { font-weight: 400; font-style: italic; }
.align-row.disorder-row .align-name { color: var(--text-muted); }
.align-row.disorder-row .align-seq span { width: 0.85em; height: 0.9em; vertical-align: middle; }
.muted { color: var(--text-muted); font-size: 0.85rem; }
</style>
</head><body>
<div class="page-header">
  <div class="page-header-top">
    <h1>FLIPPer Candidate Report</h1>
    <div class="page-header-actions">
      <select id="sort-select" class="sort-select">
        <option value="score-desc">Sort: score (high &rarr; low)</option>
        <option value="score-asc">Sort: score (low &rarr; high)</option>
        <option value="id-asc">Sort: ID (A &rarr; Z)</option>
        <option value="id-desc">Sort: ID (Z &rarr; A)</option>
      </select>
      <button type="button" id="toggle-all" class="toggle-all-btn">Collapse all</button>
    </div>
  </div>
  <div class="summary">%%SUMMARY%%</div>
</div>
<div class="cards">
%%CARDS%%
</div>
<script>
function initChart(container) {
  var data = JSON.parse(container.getAttribute('data-chart'));
  var scores = data.scores, seq = data.seq, blocks = data.blocks, plddt = data.plddt || [];
  var hasPlddt = plddt.length === scores.length;
  var n = scores.length;
  if (n === 0) return;
  var W = 1000;
  var padL = 34, padR = 10;
  var padT = 10, plotH = 96;
  var blockGap = 12, blockH = 16, blockRowGap = 4;
  var numRows = blocks.length;
  var blocksH = numRows > 0 ? numRows * blockH + (numRows - 1) * blockRowGap : 0;
  var axisGap = 8, axisH = 16, padB = 6;
  var blockTop = padT + plotH + blockGap;
  var axisY = blockTop + blocksH + axisGap;
  var H = axisY + axisH + padB;

  function xAt(i) { return padL + (n <= 1 ? 0 : (i / (n - 1)) * (W - padL - padR)); }
  function yAt(s) { return padT + (1 - s) * plotH; }

  var svgNS = 'http://www.w3.org/2000/svg';
  var svg = document.createElementNS(svgNS, 'svg');
  svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
  svg.setAttribute('class', 'chart-svg');

  [0, 0.5, 1].forEach(function (v) {
    var y = yAt(v);
    var gline = document.createElementNS(svgNS, 'line');
    gline.setAttribute('x1', padL); gline.setAttribute('x2', W - padR);
    gline.setAttribute('y1', y); gline.setAttribute('y2', y);
    gline.setAttribute('class', v === 0.5 ? 'chart-threshold' : 'chart-gridline');
    svg.appendChild(gline);
    var label = document.createElementNS(svgNS, 'text');
    label.setAttribute('x', padL - 6); label.setAttribute('y', y + 3);
    label.setAttribute('class', 'chart-axis-text'); label.setAttribute('text-anchor', 'end');
    label.textContent = v;
    svg.appendChild(label);
  });

  var areaD = 'M ' + xAt(0) + ' ' + (padT + plotH);
  for (var i = 0; i < n; i++) areaD += ' L ' + xAt(i) + ' ' + yAt(scores[i]);
  areaD += ' L ' + xAt(n - 1) + ' ' + (padT + plotH) + ' Z';
  var area = document.createElementNS(svgNS, 'path');
  area.setAttribute('d', areaD); area.setAttribute('class', 'chart-area');
  svg.appendChild(area);

  var lineD = '';
  for (var i = 0; i < n; i++) lineD += (i === 0 ? 'M ' : ' L ') + xAt(i) + ' ' + yAt(scores[i]);
  var linePath = document.createElementNS(svgNS, 'path');
  linePath.setAttribute('d', lineD); linePath.setAttribute('class', 'chart-line');
  svg.appendChild(linePath);

  if (hasPlddt) {
    var plddtD = '';
    for (var i = 0; i < n; i++) plddtD += (i === 0 ? 'M ' : ' L ') + xAt(i) + ' ' + yAt(plddt[i]);
    var plddtPath = document.createElementNS(svgNS, 'path');
    plddtPath.setAttribute('d', plddtD); plddtPath.setAttribute('class', 'chart-line-plddt');
    svg.appendChild(plddtPath);
  }

  blocks.forEach(function (row, rowIdx) {
    var rowY = blockTop + rowIdx * (blockH + blockRowGap);
    row.forEach(function (triple, i) {
      var x0 = xAt(triple[0] - 1), x1 = xAt(triple[1] - 1);
      var rect = document.createElementNS(svgNS, 'rect');
      rect.setAttribute('x', x0); rect.setAttribute('y', rowY);
      rect.setAttribute('width', Math.max(1.5, x1 - x0)); rect.setAttribute('height', blockH);
      rect.setAttribute('rx', 3);
      rect.setAttribute('class', 'chart-block');
      // colour computed server-side (build_candidate_report) from the same (copy index, region)
      // pair used for that copy's dot in the alignment view below, so the two stay in sync
      rect.style.fill = triple[2];
      svg.appendChild(rect);
    });
  });

  var tickCount = Math.min(5, n - 1);
  for (var t = 0; t <= tickCount; t++) {
    var idx = tickCount === 0 ? 0 : Math.round((t / tickCount) * (n - 1));
    var x = xAt(idx);
    var tick = document.createElementNS(svgNS, 'text');
    tick.setAttribute('x', x); tick.setAttribute('y', axisY + 11);
    tick.setAttribute('class', 'chart-axis-text');
    tick.setAttribute('text-anchor', t === 0 ? 'start' : (t === tickCount ? 'end' : 'middle'));
    tick.textContent = idx + 1;
    svg.appendChild(tick);
  }

  var guide = document.createElementNS(svgNS, 'line');
  guide.setAttribute('y1', padT); guide.setAttribute('y2', blockTop + blocksH);
  guide.setAttribute('class', 'chart-guide');
  svg.appendChild(guide);
  var dot = document.createElementNS(svgNS, 'circle');
  dot.setAttribute('r', 3.2); dot.setAttribute('class', 'chart-dot');
  svg.appendChild(dot);
  var plddtDot = null;
  if (hasPlddt) {
    plddtDot = document.createElementNS(svgNS, 'circle');
    plddtDot.setAttribute('r', 3.2); plddtDot.setAttribute('class', 'chart-dot chart-dot-plddt');
    svg.appendChild(plddtDot);
  }

  var capture = document.createElementNS(svgNS, 'rect');
  capture.setAttribute('x', 0); capture.setAttribute('y', 0);
  capture.setAttribute('width', W); capture.setAttribute('height', H);
  capture.setAttribute('class', 'chart-capture');
  svg.appendChild(capture);

  container.appendChild(svg);
  var tooltip = container.querySelector('.chart-tooltip');

  var card = container.closest('.card');
  var seqSpans = card ? card.querySelectorAll('.sequence .seq-res') : [];
  // the alignment view's per-copy residues carry the real sequence position they came from
  // (data-pos, set in build_candidate_report) so they can be matched up with the chart/full
  // sequence by position, the same way seqSpans is matched by array index
  var alignSpans = card ? card.querySelectorAll('.alignment .seq-res[data-pos]') : [];
  var alignByPos = {};
  alignSpans.forEach(function (el) {
    var p = el.getAttribute('data-pos');
    (alignByPos[p] = alignByPos[p] || []).push(el);
  });
  var activeEls = [];
  function highlightResidue(idx) {
    clearHighlight();
    var seqSpan = seqSpans[idx];
    if (seqSpan) {
      seqSpan.classList.add('active');
      activeEls.push(seqSpan);
      seqSpan.scrollIntoView({ block: 'nearest', inline: 'nearest' });
    }
    var matches = alignByPos[idx + 1];
    if (matches) {
      matches.forEach(function (el) { el.classList.add('active'); activeEls.push(el); });
    }
  }
  function clearHighlight() {
    activeEls.forEach(function (el) { el.classList.remove('active'); });
    activeEls = [];
  }

  function showPosition(idx) {
    guide.setAttribute('x1', xAt(idx)); guide.setAttribute('x2', xAt(idx));
    guide.style.opacity = 1;
    dot.setAttribute('cx', xAt(idx)); dot.setAttribute('cy', yAt(scores[idx]));
    dot.style.opacity = 1;
    if (plddtDot) {
      plddtDot.setAttribute('cx', xAt(idx)); plddtDot.setAttribute('cy', yAt(plddt[idx]));
      plddtDot.style.opacity = 1;
    }
    tooltip.style.opacity = 1;
    tooltip.style.left = ((xAt(idx) / W) * 100) + '%';
    var text = 'Position ' + (idx + 1) + ' (' + (seq[idx] || '?') + ') \\u00b7 disorder ' + scores[idx].toFixed(2);
    if (hasPlddt) text += ' \\u00b7 pLDDT ' + Math.round(plddt[idx] * 100);
    tooltip.textContent = text;
    highlightResidue(idx);
  }
  function hidePosition() {
    guide.style.opacity = 0; dot.style.opacity = 0; tooltip.style.opacity = 0;
    if (plddtDot) plddtDot.style.opacity = 0;
    clearHighlight();
  }

  capture.addEventListener('mousemove', function (evt) {
    var rect = svg.getBoundingClientRect();
    var xFrac = (evt.clientX - rect.left) / rect.width;
    var xSvg = xFrac * W;
    var idx = Math.round(((xSvg - padL) / (W - padL - padR)) * (n - 1));
    if (idx < 0) idx = 0;
    if (idx > n - 1) idx = n - 1;
    showPosition(idx);
  });
  capture.addEventListener('mouseleave', hidePosition);

  alignSpans.forEach(function (el) {
    el.addEventListener('mouseenter', function () {
      var idx = parseInt(el.getAttribute('data-pos'), 10) - 1;
      if (idx >= 0 && idx < n) showPosition(idx);
    });
    el.addEventListener('mouseleave', hidePosition);
  });
}
document.querySelectorAll('.chart').forEach(initChart);
var toggleAllBtn = document.getElementById('toggle-all');
if (toggleAllBtn) {
  toggleAllBtn.addEventListener('click', function () {
    var cards = document.querySelectorAll('details.card');
    var anyOpen = Array.prototype.some.call(cards, function (c) { return c.open; });
    cards.forEach(function (c) { c.open = !anyOpen; });
    toggleAllBtn.textContent = anyOpen ? 'Expand all' : 'Collapse all';
  });
}
var sortSelect = document.getElementById('sort-select');
var cardsContainer = document.querySelector('.cards');
if (sortSelect && cardsContainer) {
  sortSelect.addEventListener('change', function () {
    var mode = sortSelect.value;
    var cards = Array.prototype.slice.call(cardsContainer.querySelectorAll('details.card'));
    cards.sort(function (a, b) {
      if (mode === 'score-desc') return parseFloat(b.dataset.score) - parseFloat(a.dataset.score);
      if (mode === 'score-asc') return parseFloat(a.dataset.score) - parseFloat(b.dataset.score);
      if (mode === 'id-asc') return a.dataset.seqId.localeCompare(b.dataset.seqId);
      if (mode === 'id-desc') return b.dataset.seqId.localeCompare(a.dataset.seqId);
      return 0;
    });
    // re-appending an already-attached node moves it rather than duplicating it, so this
    // reorders the existing cards (and their live chart state) in place
    cards.forEach(function (c) { cardsContainer.appendChild(c); });
  });
}
</script>
</body></html>"""

## module to build a single self-contained HTML report, one "card" per final candidate, combining
## everything about that candidate that's otherwise scattered across separate files: its full
## sequence and physicochemical properties, its repeat-region stats, its repeat copies aligned
## against each other, and an interactive disorder profile with each repeat copy plotted as a
## block against the same sequence coordinates - so a candidate can be assessed at a glance
## instead of cross-referencing several separate files.
##
## report_df: one row per detected repeat region, columns ID, Begin, End, Period, Copies, Score,
## RepeatIndex, Coverage (Score is assumed higher-is-better - a caller whose native metric is
## lower-is-better, like XSTREAM's ConsensusError, must invert it before calling this).
## Optional UnitLefts/UnitRights columns (semicolon-joined 1-based inclusive start/end per repeat
## copy) draw the per-copy blocks in the disorder chart; omitted/blank rows just skip that region's
## blocks.
##
## sequences: {seq_id: full amino acid sequence} for every ID in report_df.
##
## alignment_provider(seq_id, repeat_index) -> list[(name, seq)] | None: returns that region's
## repeat-copy alignment (one entry per copy, all the same length, plus optionally one named
## 'consensus') for the alignment view, or None/empty if unavailable. Each engine supplies its own -
## DetectRepeats reads a per-hit alignment FASTA from disk, XSTREAM reconstructs it in memory from
## its own HTML report.
##
## score_meta: optional {"label": str, "explanation": str, "format": "{:.1f}"-style str} controlling
## how the Score column is labelled/formatted in the score pill and summary line.
def build_candidate_report(report_df, sequences, alignment_provider, out_html, score_meta=None):
    import json
    import html as html_module
    import pandas as pd
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    import metapredict as meta

    score_meta = score_meta or {}
    score_label = score_meta.get("label", "score")
    score_format = score_meta.get("format", "{:.1f}")
    score_explanation = score_meta.get("explanation", "")

    if report_df is None:
        report_df = pd.DataFrame()

    ## per-column average of the full protein's own metapredict disorder score across each
    ## repeat copy's real residues (not a prediction on the consensus string itself, which
    ## wouldn't correspond to any actual sequence) - walks each copy's aligned string against
    ## its real start position (lefts[i]) to map alignment columns back to disorder_scores
    def region_disorder_row(copies, lefts, disorder_scores):
        if not copies or not disorder_scores:
            return None
        aln_len = len(copies[0][1])
        sums = [0.0] * aln_len
        counts = [0] * aln_len
        for (_, aseq), left in zip(copies, lefts):
            pos = left
            for col, ch in enumerate(aseq):
                if ch != '-':
                    if 0 <= pos - 1 < len(disorder_scores):
                        sums[col] += disorder_scores[pos - 1]
                        counts[col] += 1
                    pos += 1
        return [(sums[i] / counts[i]) if counts[i] else None for i in range(aln_len)]

    def alignment_html_block(aligned_records, lefts, rights, disorder_scores):
        copies = [(name, seq) for name, seq in aligned_records if name != 'consensus']
        consensus = next((seq for name, seq in aligned_records if name == 'consensus'), None)
        if not copies:
            return "<p class='muted'>No repeat-unit alignment available (fewer than 2 copies).</p>"
        labels = ["{}-{}".format(l, r) for l, r in zip(lefts, rights)]
        ## +3 (not +1) reserves room for the copy-dot indicator prefixed onto each label below,
        ## so the position numbers still line up in a column even though copy rows now start
        ## with a dot that "consensus"/"disorder" don't have
        name_width = max([len(l) for l in labels] + [len('consensus'), len('disorder')]) + 3
        blank_dot = "<i class='copy-dot' style='background:transparent;box-shadow:none;'></i>"
        row_html = []
        ## each non-gap residue gets a data-pos (its real position in the full sequence, same
        ## coordinate space as the sequence box and disorder chart) so the hover-tracking in the
        ## report's <script> can highlight it from either place, the same way it already does
        ## for the full sequence box
        for i, ((name, seq), label, left) in enumerate(zip(copies, labels, lefts)):
            dot_color = _repeat_copy_color(i)
            pos = left
            span_parts = []
            for ch in seq:
                if ch == '-':
                    span_parts.append("<span class='weak'>-</span>")
                else:
                    span_parts.append(
                        "<span class='seq-res' data-pos='{pos}' style='background:{bg}' title='{pos}: {aa}'>{aa}</span>".format(
                            pos=pos, bg=REPEAT_AA_COLORS.get(ch, '#eee'), aa=html_module.escape(ch)))
                    pos += 1
            spans = "".join(span_parts)
            ## same colour as this copy's block in the disorder chart above (both come from
            ## _repeat_copy_color(i) for this region), so the two views can be matched by eye
            row_html.append(
                "<div class='align-row'><span class='align-name' style='width:{w}ch'>"
                "<i class='copy-dot' style='background:{dot}'></i>{label}</span>"
                "<span class='align-seq'>{spans}</span></div>".format(
                    w=name_width, dot=dot_color, label=html_module.escape(label), spans=spans))
        if consensus:
            def _cons_span(c):
                ## '-' (true gap) and other engines' ambiguous/weak-consensus symbols (e.g.
                ## DECIPHER's "+", ".", "X") both mean "no confident single residue here", and
                ## both get the same faint rendering rather than a background colour
                if c not in REPEAT_AA_COLORS:
                    return "<span class='weak'>{}</span>".format(html_module.escape(c))
                return "<span style='background:{}'>{}</span>".format(
                    REPEAT_AA_COLORS[c], html_module.escape(c))
            cons_spans = "".join(_cons_span(c) for c in consensus)
            row_html.append("<div class='align-row consensus'><span class='align-name' style='width:{w}ch'>{blank}consensus</span><span class='align-seq'>{spans}</span></div>".format(
                w=name_width, blank=blank_dot, spans=cons_spans))
            disorder_row = region_disorder_row(copies, lefts, disorder_scores)
            if disorder_row is not None:
                ## coloured relative to THIS region's own min-max range, not the fixed 0-1
                ## disorder scale - these candidates are selected for being highly disordered,
                ## so their copies often all sit in a narrow high band (e.g. 0.8-0.95), which
                ## made a fixed-scale bar look almost uniformly dark and hid the actual shape
                valid = [v for v in disorder_row if v is not None]
                lo, hi = (min(valid), max(valid)) if valid else (0.0, 0.0)
                span_range = hi - lo
                dis_spans = "".join(
                    "<span style='background:{bg}' title='avg disorder {v}'></span>".format(
                        bg=(
                            "color-mix(in srgb, var(--accent) {}%, var(--surface-2))".format(
                                round(((v - lo) / span_range if span_range > 1e-9 else 0.5) * 100))
                            if v is not None else "transparent"
                        ),
                        v="{:.2f}".format(v) if v is not None else "n/a")
                    for v in disorder_row
                )
                row_html.append(
                    "<div class='align-row disorder-row'><span class='align-name' style='width:{w}ch' "
                    "title='colour scaled to this region&#39;s own range: {lo:.2f}-{hi:.2f}'>{blank}disorder</span>"
                    "<span class='align-seq'>{spans}</span></div>".format(w=name_width, lo=lo, hi=hi, blank=blank_dot, spans=dis_spans))
        return "<div class='alignment'>" + "".join(row_html) + "</div>"

    def score_pill_class(score, max_score):
        if max_score <= 0:
            return "mid"
        frac = score / max_score
        if frac >= 0.66:
            return "high"
        if frac >= 0.33:
            return "mid"
        return "low"

    def properties_chip_row(seq):
        if not seq:
            return ""
        analysis = ProteinAnalysis(seq)
        return (
            "<div class='region-chip'><span class='region-chip-label'>pI</span><span class='region-chip-stat'><b>{pi:.2f}</b></span></div>"
            "<div class='region-chip'><span class='region-chip-label'>Aromaticity</span><span class='region-chip-stat'><b>{arom:.1%}</b></span></div>"
            "<div class='region-chip'><span class='region-chip-label'>Instability</span><span class='region-chip-stat'><b>{inst:.0f}</b></span></div>"
            "<div class='region-chip'><span class='region-chip-label'>GRAVY</span><span class='region-chip-stat'><b>{gravy:.2f}</b></span></div>"
        ).format(
            pi=analysis.isoelectric_point(),
            arom=analysis.aromaticity(),
            inst=analysis.instability_index(),
            gravy=analysis.gravy(),
        )

    ## colours each residue by amino acid physicochemical class, the same REPEAT_AA_COLORS
    ## scheme used in the repeat-unit alignment below, so a residue's colour means the same thing
    ## in both places instead of the sequence box using a separate disorder-based encoding
    def sequence_html(seq):
        if not seq:
            return "<pre class='sequence'>(sequence not found)</pre>"
        spans = []
        for i, aa in enumerate(seq):
            color = REPEAT_AA_COLORS.get(aa, '#eee')
            spans.append(
                "<span class='seq-res' data-pos='{pos}' style='background:{color}' title='{pos}: {aa}'>{aa}</span>".format(
                    color=color, pos=i + 1, aa=html_module.escape(aa)))
        return "<pre class='sequence'>" + "".join(spans) + "</pre>"

    ## NCBI's interactive BLAST page accepts a query pre-filled straight from the URL (the same
    ## QUERY GET param pattern UniProt/other sequence databases use for their own "BLAST this"
    ## links) - no server-side submission or API key needed, it just opens the page with the
    ## sequence already in the search box for the user to launch themselves.
    def blast_button(seq):
        if not seq:
            return ""
        from urllib.parse import quote
        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&QUERY=" + quote(seq)
        return (
            "<a class='blast-btn' href='{url}' target='_blank' rel='noopener noreferrer' "
            "title='Open this sequence in NCBI BLASTp'>BLASTp &#8599;</a>"
        ).format(url=url)

    def chart_html(seq, disorder_scores, plddt_scores, unit_rows):
        if not seq:
            return "<p class='muted'>Sequence not found.</p>"
        payload = {
            "scores": disorder_scores,
            "plddt": plddt_scores,
            "seq": seq,
            "blocks": unit_rows,
        }
        return (
            "<div class='chart' data-chart='{data}'>"
            "<div class='chart-tooltip'></div>"
            "</div>"
        ).format(data=html_module.escape(json.dumps(payload), quote=True))

    ## static (non-interactive) rendering of the same disorder curve + repeat blocks as chart_html,
    ## sized to sit inside a collapsed card's <summary> - a card that's closed hides everything
    ## except its summary, so this is what makes a collapsed card scannable at a glance instead of
    ## needing to be opened first. Built server-side rather than reusing initChart() in JS because
    ## it needs no hover/tooltip behaviour, just the shape.
    def preview_chart_svg(disorder_scores, unit_rows):
        n = len(disorder_scores)
        if n == 0:
            return ""
        W = 1000
        plotH = 34
        blockH = 6
        blockRowGap = 3
        blockGap = 6
        numRows = len(unit_rows)
        blocksH = numRows * blockH + max(0, numRows - 1) * blockRowGap
        blockTop = plotH + (blockGap if numRows else 0)
        H = blockTop + blocksH

        def xAt(i):
            return 0.0 if n <= 1 else (i / (n - 1)) * W

        def yAt(s):
            return (1 - s) * plotH

        line_d = " ".join(
            "{} {:.1f} {:.1f}".format("M" if i == 0 else "L", xAt(i), yAt(s))
            for i, s in enumerate(disorder_scores)
        )
        area_d = "M {:.1f} {:.1f} L {} L {:.1f} {:.1f} Z".format(
            xAt(0), plotH, line_d[2:], xAt(n - 1), plotH)

        blocks_svg = []
        for row_idx, row in enumerate(unit_rows):
            row_y = blockTop + row_idx * (blockH + blockRowGap)
            for (left, right, color) in row:
                x0, x1 = xAt(left - 1), xAt(right - 1)
                blocks_svg.append(
                    "<rect x='{:.1f}' y='{:.1f}' width='{:.1f}' height='{}' rx='1.5' class='preview-block' style='fill:{}'></rect>".format(
                        x0, row_y, max(1.5, x1 - x0), blockH, color))

        return (
            "<svg class='preview-svg' viewBox='0 0 {W} {H}' preserveAspectRatio='none'>"
            "<path d='{area_d}' class='preview-area'></path>"
            "<path d='{line_d}' class='preview-line'></path>"
            "{blocks}"
            "</svg>"
        ).format(W=W, H=H, area_d=area_d, line_d=line_d, blocks="".join(blocks_svg))

    cards = []
    all_scores = report_df['Score'].tolist() if not report_df.empty else []
    max_score = max(all_scores) if all_scores else 0
    if not report_df.empty:
        ## show the candidates the engine was most confident about first
        order = report_df.groupby('ID')['Score'].max().sort_values(ascending=False).index
        for seq_id in order:
            group = report_df[report_df['ID'] == seq_id]
            seq = sequences.get(seq_id, "")
            top_score = group['Score'].max()
            disorder_scores = [round(float(s), 4) for s in meta.predict_disorder(seq)] if seq else []
            plddt_scores = [round(float(p), 4) for p in meta.predict_pLDDT(seq, return_decimals=True)] if seq else []
            region_blocks = []
            unit_rows = []
            for _, row in group.iterrows():
                repeat_index = int(row['RepeatIndex'])
                coverage = row['Coverage'] if 'Coverage' in row else float('nan')
                ## the region's own stats (position/period/copies/coverage) sit right above its
                ## alignment rather than bunched with the whole-sequence properties at the top of
                ## the card, so they read together with the repeat they describe
                region_chip = (
                    "<div class='region-chips'><div class='region-chip'><span class='region-chip-label'>Region {idx}</span>"
                    "<span class='region-chip-stat'><b>{begin:.0f}</b>-<b>{end:.0f}</b> aa</span>"
                    "<span class='region-chip-stat'>{period:.0f} aa &times; <b>{copies:.0f}</b> copies</span>"
                    "<span class='region-chip-stat'>{coverage:.0%} coverage</span></div></div>"
                ).format(idx=repeat_index, begin=row['Begin'], end=row['End'],
                         period=row['Period'], copies=row['Copies'], coverage=coverage)
                lefts, rights = [], []
                if 'UnitLefts' in row and 'UnitRights' in row and pd.notna(row['UnitLefts']):
                    lefts = [int(x) for x in str(row['UnitLefts']).split(';')]
                    rights = [int(x) for x in str(row['UnitRights']).split(';')]
                    ## same per-index colour as this copy's dot in the alignment view, so the
                    ## chart block and its alignment row can be matched by eye
                    unit_rows.append([(l, r, _repeat_copy_color(i)) for i, (l, r) in enumerate(zip(lefts, rights))])
                aligned_records = alignment_provider(seq_id, repeat_index)
                if aligned_records:
                    block = alignment_html_block(aligned_records, lefts, rights, disorder_scores)
                else:
                    block = "<p class='muted'>No alignment available.</p>"
                region_blocks.append("<div class='region-block'>" + region_chip + block + "</div>")

            cards.append("""
<details class="card" data-score="{score_attr}" data-seq-id="{seq_id_attr}" open>
  <summary class="card-head">
    <div class="card-head-row">
      <h2>{seq_id}</h2>
      <span class="pill pill-{pill_class}" title="{score_explanation}">{score_label} {score_display}</span>
    </div>
    <div class="card-preview">{preview}</div>
  </summary>
  <h4>Whole-sequence properties</h4>
  <div class="region-chips">{properties}</div>
  <div class="section-head">
    <h3>Sequence <span class="dim">&middot; {seq_len} aa &middot; coloured by residue type</span></h3>
    {blast_button}
  </div>
  {seq}
  <div class="section-head">
    <h3>Disorder profile &amp; repeat positions <span class="dim">&middot; hover to trace position</span></h3>
    <div class="chart-legend">
      <span><i class="legend-dot disorder"></i>Disorder</span>
      <span><i class="legend-dot plddt"></i>Predicted AF2 pLDDT</span>
    </div>
  </div>
  {chart}
  <h3>Repeat-unit alignment <span class="dim">&middot; copies aligned against consensus; bottom row is the mean per-residue disorder across copies, scaled to each region's own range, darker = more disordered; each copy's dot matches its block in the diagram above</span></h3>
  {alignment_blocks}
</details>""".format(
                seq_id=html_module.escape(seq_id),
                seq_id_attr=html_module.escape(seq_id, quote=True),
                score_attr="{:.4f}".format(float(top_score)),
                pill_class=score_pill_class(top_score, max_score),
                score_label=html_module.escape(score_label),
                score_display=html_module.escape(score_format.format(top_score)),
                score_explanation=html_module.escape(score_explanation, quote=True),
                properties=properties_chip_row(seq),
                seq_len=len(seq),
                seq=sequence_html(seq),
                blast_button=blast_button(seq),
                chart=chart_html(seq, disorder_scores, plddt_scores, unit_rows),
                preview=preview_chart_svg(disorder_scores, unit_rows),
                alignment_blocks="".join(region_blocks) if region_blocks else "<p class='muted'>No repeat regions.</p>",
            ))

    summary = (
        "{n} candidate{plural}, {label} {lo}-{hi}".format(
            n=len(cards), plural="" if len(cards) == 1 else "s", label=score_label,
            lo=score_format.format(min(all_scores)), hi=score_format.format(max(all_scores)))
        if all_scores else "no candidates passed filtering"
    )

    ## built with %%TOKEN%% placeholders + str.replace(), not str.format() - the inline <script>
    ## below is full of literal { } braces (JS objects, functions) that would otherwise all need
    ## doubling to survive .format()'s escaping, which is unmaintainable at this size
    page = _CANDIDATE_REPORT_TEMPLATE
    page = page.replace("%%SUMMARY%%", html_module.escape(summary))
    page = page.replace("%%CARDS%%", "\n".join(cards))

    with open(out_html, "w") as f:
        f.write(page)
