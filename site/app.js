"use strict";
/*
 * ICP Immune Map - dashboard script (no build step: edit and reload).
 *
 * Reads data/datasets.json (the list of datasets) and data/<dataset>.json, both written by
 * analysis/step4_export_site_data.py, and draws every chart and table with Plotly.
 *
 * Sections: helpers · start-up & loading · header & tiles · immune processes (ImmPort) ·
 *           immune cell types (atlas) · filters · arrow chart · gene table · gene detail ·
 *           heatmap · quality checks · chart/table toggle
 */

const FONT = 'system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif';
const PLOT_CONFIG = { displayModeBar: false, responsive: true };
const TOP_ARROWS = 30;        // genes shown in the arrow chart
const TOP_HEATMAP = 50;       // genes shown in the heatmap
const TOP_IMMUNE_LISTS = 15;  // immune lists shown in the immune-process chart
const ROW_PX = 26;            // height of one row in the bar, arrow and cell-type charts
const TABLE_ROWS = 25;        // rows shown in the gene table before "Show all"

const state = {
  data: null,                 // the dataset JSON
  lists: new Map(),           // ImmPort list id -> its summary (name, link, counts)
  cellTypes: new Map(),       // atlas population id -> its scores (only populations that could be scored)
  filters: { search: "", direction: "all", immuneOnly: false, list: "" },
  sort: { key: "pvalue", ascending: true },
  showAllRows: false,         // gene table: all rows, or only the first TABLE_ROWS
  selected: null,             // id of the gene shown in the gene panel
  cellSection: null,          // atlas section shown in the cell-type chart
  selectedCell: null,         // atlas population shown in the cell-type panel
};

/* ------------------------------------------------------------------ helpers */

const $ = (id) => document.getElementById(id);

/** Create an element. Strings become text nodes (never HTML), so data values cannot inject markup. */
function h(tag, props = {}, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(props)) {
    if (value == null || value === false) continue;
    if (key === "class") node.className = value;
    else if (key.startsWith("on")) node.addEventListener(key.slice(2), value);
    else node.setAttribute(key, value === true ? "" : String(value));
  }
  for (const child of children.flat()) {
    if (child == null || child === false) continue;
    node.append(child instanceof Node ? child : document.createTextNode(String(child)));
  }
  return node;
}

/** Value of a color token defined in style.css. */
const token = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

const fmt = {
  int: (n) => Math.round(n).toLocaleString("en-US"),
  count: (n) => (n == null ? "–" : n >= 1e6 ? `${(n / 1e6).toFixed(1)}M` : n >= 1e4 ? `${Math.round(n / 1e3)}K`
    : n >= 1e3 ? `${(n / 1e3).toFixed(1)}K` : n >= 10 ? String(Math.round(n)) : n.toFixed(1)),
  p: (p) => (p == null ? "–" : p < 0.001 ? p.toExponential(1) : p.toFixed(3)),
  log2fc: (x) => `${x > 0 ? "+" : "−"}${Math.abs(x).toFixed(2)}`,
  fold: (x) => `${(2 ** Math.abs(x)).toFixed(1)}× ${x > 0 ? "higher" : "lower"}`,
  sd: (x) => `${x > 0 ? "+" : "−"}${Math.abs(x).toFixed(2)} SD`,
  number: (x) => String(Number(x.toFixed(2))),
};

const byPvalue = (a, b) => a.pvalue - b.pvalue;
/** Words for the dataset's measurements: RNA-seq read counts or microarray intensities. */
const units = () => (state.data.measure === "intensity"
  ? { axis: "Normalized intensity", value: "intensity (normalized)", mean: "Mean intensity", source: "on the array",
      detail: "Normalized microarray intensities (quantile-normalized across arrays)." }
  : { axis: "Normalized read count", value: "reads (normalized)", mean: "Mean count", source: "in the file",
      detail: "Normalized read counts (adjusted for sequencing depth)." });
const adjustedForSex = () => state.data.design.includes("fetal_sex");
const isImmune = (gene) => gene.immune;
const geneLabel = (gene) => gene.symbol + (isImmune(gene) ? " ◆" : "") + (gene.raw_disagrees ? " ⚠" : "");
const shorten = (text, max) => (text.length > max ? `${text.slice(0, max - 1)}…` : text);
const sexName = (sex) => ({ F: "female", M: "male" }[sex] ?? "unknown");
/** Atlas IDs sort by number: 1.2 before 1.10. */
const byPopulationId = (a, b) => a.localeCompare(b, undefined, { numeric: true });
/** Plotly hover labels understand a little HTML, so escape names taken from the data. */
const plotlyText = (text) => String(text).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;");

function groupColor(group) {
  const d = state.data;
  return group === d.reference_group ? token("--healthy") : group === d.compare_group ? token("--icp") : token("--ink-muted");
}

function median(values) {
  const sorted = [...values].sort((a, b) => a - b);
  const mid = sorted.length >> 1;
  return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
}

function countBy(items, key) {
  const counts = new Map();
  for (const item of items) counts.set(key(item), (counts.get(key(item)) ?? 0) + 1);
  return counts;
}

function chromosomeRank(chromosome) {
  if (chromosome == null) return 99;
  const n = Number(chromosome);
  return Number.isInteger(n) ? n : ({ X: 23, Y: 24, MT: 25 }[chromosome] ?? 26);
}

/** Plotly axis styled with the page's tokens: hairline grid, muted tick labels. */
function axis(title = "", extra = {}) {
  return {
    title: { text: title, font: { color: token("--ink-2") }, standoff: 10 },
    gridcolor: token("--grid"),
    linecolor: token("--axis"),
    zerolinecolor: token("--axis"),
    tickfont: { color: token("--ink-muted") },
    automargin: true,
    fixedrange: true,
    ...extra,
  };
}

function layout(extra = {}) {
  return {
    paper_bgcolor: token("--surface"),
    plot_bgcolor: token("--surface"),
    font: { family: FONT, size: 12, color: token("--ink-2") },
    margin: { l: 10, r: 20, t: 36, b: 40 },
    hoverlabel: { bgcolor: token("--surface"), bordercolor: token("--axis"),
                  font: { family: FONT, size: 12, color: token("--ink-1") } },
    showlegend: false,
    legend: { orientation: "h", x: 0, xanchor: "left", y: 1, yanchor: "bottom",
              font: { color: token("--ink-2") }, bgcolor: "rgba(0,0,0,0)" },
    ...extra,
  };
}

/** Draw (or redraw) a chart, and attach its click handler once. */
function draw(id, traces, chartLayout, onClick) {
  const div = $(id);
  Plotly.react(div, traces, chartLayout, PLOT_CONFIG);
  if (onClick && !div.dataset.clickBound) {
    div.on("plotly_click", onClick);
    div.dataset.clickBound = "yes";
  }
}

function clearChart(id) {
  Plotly.purge($(id));
  delete $(id).dataset.clickBound;
}

/** A plain table: columns are [heading, value function, numeric?]. */
function simpleTable(columns, rows) {
  return h("table", {},
    h("thead", {}, h("tr", {}, columns.map(([heading, , numeric]) => h("th", { scope: "col", class: numeric ? "num" : null }, heading)))),
    h("tbody", {}, rows.map((row) => h("tr", {}, columns.map(([, value, numeric]) => h("td", { class: numeric ? "num" : null }, value(row)))))));
}

/** One dot per patient in a column per group; shape = fetal sex; a line marks each group's median. */
function drawPatientDots(id, values, { yTitle, format, unit, zeroLine = false }) {
  const d = state.data;
  const groups = [...new Set(d.samples.map((s) => s.group))];
  const traces = groups.map((group, column) => {
    const members = d.samples.map((s, i) => ({ s, i })).filter(({ s }) => s.group === group);
    const spread = Math.min(0.08, 0.5 / members.length);  // spread the dots so they don't overlap
    return {
      type: "scatter", mode: "markers", name: group,
      x: members.map((_, k) => column + (k - (members.length - 1) / 2) * spread),
      y: members.map(({ i }) => values[i]),
      customdata: members.map(({ s, i }) => [plotlyText(s.id), format(values[i]), sexName(s.fetal_sex)]),
      hovertemplate: `<b>%{customdata[0]}</b><br>%{customdata[1]} ${unit}<br>fetal sex: %{customdata[2]}<extra></extra>`,
      marker: { color: groupColor(group), size: 11, symbol: members.map(({ s }) => (s.fetal_sex === "M" ? "diamond" : "circle")),
                line: { color: token("--surface"), width: 2 } },
    };
  });
  const medians = groups.map((group, column) => {
    const y = median(d.samples.flatMap((s, i) => (s.group === group ? [values[i]] : [])));
    return { type: "line", xref: "x", yref: "y", x0: column - 0.3, x1: column + 0.3, y0: y, y1: y,
             line: { color: token("--ink-2"), width: 2 } };
  });
  draw(id, traces, layout({
    height: 260,
    margin: { l: 10, r: 10, t: 10, b: 30 },
    shapes: medians,
    xaxis: axis("", { tickvals: groups.map((_, i) => i), ticktext: groups, range: [-0.6, groups.length - 0.4],
                      showgrid: false, zeroline: false, tickfont: { color: token("--ink-2"), size: 13 } }),
    yaxis: axis(yTitle, zeroLine ? { zeroline: true } : { rangemode: "tozero" }),
  }));
}

/* ------------------------------------------------------- start-up & loading */

async function fetchJSON(path) {
  const response = await fetch(path, { cache: "no-cache" });
  if (!response.ok) throw new Error(`${path}: HTTP ${response.status}`);
  return response.json();
}

function showError(message, error) {
  const box = $("load-error");
  box.textContent = message + (error ? ` (${error.message})` : "");
  box.hidden = false;
}

async function start() {
  let datasets;
  try {
    datasets = await fetchJSON("data/datasets.json");
  } catch (error) {
    showError("Could not load data/datasets.json. Run ./run_all.sh to create it, then serve the site folder "
              + "over HTTP: uv run python -m http.server --directory site 8000", error);
    return;
  }
  if (!datasets.length) return showError("data/datasets.json lists no datasets.");

  const select = $("dataset-select");
  select.replaceChildren(...datasets.map((d) => h("option", { value: d.id }, d.tissue ? `${d.id} · ${d.tissue}` : d.id)));
  const fromLink = () => datasets.find((d) => d.id === decodeURIComponent(location.hash.slice(1)));
  select.addEventListener("change", () => { location.hash = select.value; });
  window.addEventListener("hashchange", () => { if (fromLink()) loadDataset(fromLink()); });

  setUpFilters();
  setUpTableToggles();
  $("cells-section").addEventListener("change", (e) => { state.cellSection = e.target.value; state.selectedCell = null; renderCellTypes(); });
  $("qc").addEventListener("toggle", () => { if ($("qc").open) renderQualityChecks(); });

  await loadDataset(fromLink() ?? datasets[0]);
}

async function loadDataset(entry) {
  $("dataset-select").value = entry.id;
  try {
    state.data = await fetchJSON(`data/${entry.file}`);
  } catch (error) {
    showError(`Could not load data/${entry.file}.`, error);
    return;
  }
  $("load-error").hidden = true;
  state.lists = new Map(state.data.immune_lists.map((list) => [list.id, list]));
  state.cellTypes = new Map(state.data.cell_types.map((cell) => [cell.id, cell]));
  state.selected = null;
  state.cellSection = null;
  state.selectedCell = null;
  resetFilters();
  renderAll();
}

function renderAll() {
  if (!state.data) return;
  renderHeader();
  renderTiles();
  renderCellMarkers();
  renderImmuneMap();
  renderCellTypes();
  renderFilterOptions();
  renderFiltered();
  if ($("qc").open) renderQualityChecks();
}

/* ----------------------------------------------------------- header & tiles */

function geoLink() {
  const d = state.data;
  if (!d.geo_accession) return d.id;
  return h("a", { href: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${encodeURIComponent(d.geo_accession)}` },
    `GEO ${d.geo_accession}`);
}

/** On GitHub Pages (user.github.io/repo/), link back to the repository. */
function githubRepoUrl() {
  const user = location.hostname.match(/^([^.]+)\.github\.io$/)?.[1];
  if (!user) return null;
  const repo = location.pathname.split("/").filter(Boolean)[0];
  return `https://github.com/${user}/${repo ?? `${user}.github.io`}`;
}

function renderHeader() {
  const d = state.data;
  const ref = d.reference_group, cmp = d.compare_group;
  const sizes = countBy(d.samples, (s) => s.group);
  $("dataset-facts").replaceChildren(
    geoLink(), ` · ${d.title} · ${sizes.get(ref)} ${ref} and ${sizes.get(cmp)} ${cmp} samples`,
    ` · ${d.method} (${d.software}), model ${d.design}`,
    d.chr_y_genes_excluded ? ` · fetal sex not adjusted for; ${d.chr_y_genes_excluded} chrY genes left out` : "");

  const c = d.cutoffs;
  const pRule = c.padj != null ? `FDR-adjusted p-value is below ${c.padj}` : `p-value is below ${c.pvalue}`;
  const note = [`A gene counts as significant when its ${pRule} and its expression is at least ${fmt.number(2 ** c.abs_log2fc)}-fold `
    + `higher or lower in ${cmp} than in ${ref} (|log2 fold change| ≥ ${c.abs_log2fc}).`];
  if (c.padj == null) {
    note.push(" These p-values are not corrected for testing thousands of genes at once, so some of these genes are "
      + "expected to be chance findings. The FDR column gives the corrected value.");
  }
  $("cutoff-note").replaceChildren(...note);
  $("cutoff-note").hidden = false;

  $("footer-geo").replaceChildren(geoLink());
  $("footer-method").replaceChildren(...(d.measure === "intensity"
    ? ["analysed with a limma-style moderated t-test (Smyth, ", h("em", {}, "Stat Appl Genet Mol Biol"), " 2004)"]
    : ["analysed with DESeq2 (Love, Huber & Anders, ", h("em", {}, "Genome Biology"), " 2014) via ",
       h("a", { href: "https://github.com/owkin/PyDESeq2" }, "pydeseq2")]));
  $("footer-immport").hidden = !d.immport;
  $("table-immune-note").textContent = d.immport
    ? "◆ = immune-related: a marker gene of an immune cell population in the Human Immune Cell Atlas, or in an ImmPort immune gene list."
    : "◆ = immune-related: a marker gene of an immune cell population in the Human Immune Cell Atlas.";
  const repo = githubRepoUrl();
  if (repo) {
    $("footer-atlas").replaceChildren(h("a", { href: `${repo}/blob/main/resources/human_immune_cell_atlas.txt` }, "Human Immune Cell Atlas"),
      " (", h("a", { href: `${repo}/blob/main/resources/immune_cell_markers.tsv` }, "marker table"), ")");
  }
  $("footer-meta").replaceChildren(`Data exported ${d.generated}.`, repo ? " " : "", repo ? h("a", { href: repo }, "View the code on GitHub") : "");
}

function directionSplit(higher, lower, cmp) {
  return [h("span", { class: "swatch higher" }), `${higher} higher in ${cmp}`, "  ·  ",
          h("span", { class: "swatch lower" }), `${lower} lower`];
}

function renderTiles() {
  const d = state.data, cmp = d.compare_group;
  const higher = d.genes.filter((g) => g.log2fc > 0).length;
  const immune = d.genes.filter(isImmune);
  const immuneHigher = immune.filter((g) => g.log2fc > 0).length;
  $("tile-tested").textContent = fmt.int(d.genes_tested);
  $("tile-tested-sub").textContent = `of ${fmt.int(d.genes_in_file)} ${units().source} · ${fmt.int(d.immune_genes_tested)} immune-related`;
  $("tile-significant").textContent = fmt.int(d.genes.length);
  $("tile-significant-sub").replaceChildren(...directionSplit(higher, d.genes.length - higher, cmp));
  $("tile-immune-label").textContent = d.immport ? "Immune-related significant genes" : "Significant immune cell marker genes";
  $("tile-immune").textContent = fmt.int(immune.length);
  $("tile-immune-sub").replaceChildren(...directionSplit(immuneHigher, immune.length - immuneHigher, cmp));

  // the major-lineage population (first atlas section) with the smallest p-value
  const lead = d.cell_types.filter((c) => c.section === d.atlas_sections[0]).sort(byPvalue)[0];
  $("tile-cell").textContent = lead ? lead.name : "–";
  $("tile-cell-sub").replaceChildren(...(lead
    ? [h("span", { class: `swatch ${lead.effect > 0 ? "higher" : "lower"}` }),
       `${lead.effect > 0 ? "higher" : "lower"} in ${cmp} (${fmt.sd(lead.effect)}) · p = ${fmt.p(lead.pvalue)} · FDR ${fmt.p(lead.fdr)}`]
    : ["No atlas population could be scored"]));
}

/* -------------------------------------- immune cells and their changed marker genes */

/**
 * Grid: rows are atlas immune cell populations, columns are significant genes that are markers of
 * at least one of them. A colored square means "this gene is a marker of this cell population"; its
 * color and arrow show how the gene changed in ICP (darker = bigger change).
 */
function renderCellMarkers() {
  const d = state.data, cmp = d.compare_group;
  const cells = d.cell_marker_changes;
  const byId = new Map(d.genes.map((g) => [g.id, g]));
  const change = (id) => byId.get(id).log2fc;
  const sectionNumber = (c) => c.section.split(".")[0];

  $("markers-empty").hidden = cells.length > 0;
  $("markers-caption").textContent = "Each row is an immune cell population from the Human Immune Cell Atlas, each column one "
    + "of the significant genes that is a marker of it. A square means the gene is a marker of that population: ▲ orange = "
    + `higher in ${cmp}, ▼ blue = lower, darker = bigger change (hover for the fold change and p-value). Populations share `
    + "markers, e.g. PRF1 marks NK cells and cytotoxic T cells, so one gene can appear in several rows. § is the atlas "
    + "section. Click a square to show that population's genes below.";
  $("markers-table").replaceChildren(simpleTable([
    ["Immune cell population", (c) => c.name],
    ["Atlas section", (c) => c.section],
    ["Changed marker genes", (c) => c.genes.map((id) => `${byId.get(id).symbol} ${change(id) > 0 ? "▲" : "▼"} ${fmt.fold(change(id))}`).join(", ")],
    ["Markers tested", (c) => `${c.markers_tested} of ${c.markers_in_panel}`, true],
  ], cells));
  if (!cells.length) return clearChart("markers-chart");

  // columns: marker genes, most increased first; rows: populations in atlas order, top to bottom
  const genes = [...new Set(cells.flatMap((c) => c.genes))].sort((a, b) => change(b) - change(a));
  const rows = [...cells].reverse();  // Plotly draws the first row at the bottom
  const columnsX = genes.map((_, i) => i), rowsY = rows.map((_, i) => i);
  const reach = Math.max(1, ...genes.map((id) => Math.abs(change(id))));
  const scale = [[0, token("--div-low-3")], [1 / 6, token("--div-low-2")], [2 / 6, token("--div-low-1")], [0.5, token("--div-mid")],
                 [4 / 6, token("--div-high-1")], [5 / 6, token("--div-high-2")], [1, token("--div-high-3")]];
  const marks = rows.flatMap((c, y) => c.genes.map((id) => ({ x: genes.indexOf(id), y, id })));
  const height = rows.length * 26 + 90;
  // a thin line between atlas sections
  const dividers = rows.slice(1).flatMap((c, i) => (sectionNumber(c) === sectionNumber(rows[i]) ? [] : [{
    type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: i + 0.5, y1: i + 0.5, line: { color: token("--axis"), width: 1 } }]));
  $("markers-chart").style.maxWidth = `${420 + 64 * genes.length}px`;  // keep squares compact when there are few genes

  draw("markers-chart", [
    {  // light background squares, so the grid reads as a grid
      type: "heatmap", x: columnsX, y: rowsY, z: rows.map(() => genes.map(() => 0)), showscale: false, hoverinfo: "skip",
      colorscale: [[0, token("--page")], [1, token("--page")]], xgap: 3, ygap: 3,
    },
    {  // the changes
      type: "heatmap", x: columnsX, y: rowsY, xgap: 3, ygap: 3, hoverongaps: false,
      z: rows.map((c) => genes.map((id) => (c.genes.includes(id) ? change(id) : null))),
      zmin: -reach, zmax: reach, zmid: 0, colorscale: scale,
      customdata: rows.map((c) => genes.map((id) => [plotlyText(c.name), plotlyText(c.section), plotlyText(byId.get(id).symbol),
                                                     fmt.fold(change(id)), fmt.p(byId.get(id).pvalue), c.id])),
      hovertemplate: `<b>%{customdata[0]}</b><br>%{customdata[1]}<br>marker %{customdata[2]}: %{customdata[3]} in ${plotlyText(cmp)}`
        + " (p = %{customdata[4]})<extra></extra>",
      colorbar: { title: { text: "log2 fold change", side: "right", font: { color: token("--ink-2") } }, thickness: 10,
                  outlinewidth: 0, tickfont: { color: token("--ink-muted") }, len: Math.min(1, 240 / height) },
    },
    {  // ▲/▼ inside each square: white on dark squares, dark on light ones
      type: "scatter", mode: "text", x: marks.map((m) => m.x), y: marks.map((m) => m.y), hoverinfo: "skip",
      text: marks.map((m) => (change(m.id) > 0 ? "▲" : "▼")),
      textfont: { size: 11, color: marks.map((m) => (Math.abs(change(m.id)) / reach > 0.45 ? "#ffffff" : token("--ink-1"))) },
    },
  ], layout({
    height,
    margin: { l: 10, r: 10, t: 10, b: 10 },
    shapes: dividers,
    xaxis: axis("", { tickvals: columnsX, ticktext: genes.map((id) => byId.get(id).symbol + (byId.get(id).raw_disagrees ? " ⚠" : "")),
                      side: "top", showgrid: false, zeroline: false, tickfont: { color: token("--ink-1"), size: 12 } }),
    yaxis: axis("", { tickvals: rowsY, showgrid: false, zeroline: false, tickfont: { color: token("--ink-2") },
                      ticktext: rows.map((c) => `${plotlyText(shorten(c.name, 36))} <span style="color:${token("--ink-muted")}">§${sectionNumber(c)}</span>`) }),
  }), (event) => {
    const id = event.points?.find((p) => p.customdata)?.customdata?.[5];
    if (!id) return;
    state.filters.list = id;
    $("filter-list").value = id;
    renderFiltered();
    document.querySelector(".filters").scrollIntoView({ behavior: "smooth", block: "start" });
  });
}

/* ------------------------------------------- ImmPort lists (only with use_immport) */

function renderImmuneMap() {
  const d = state.data, cmp = d.compare_group;
  $("immport-card").hidden = !d.immport;
  if (!d.immport) return clearChart("immune-map-chart");
  const lists = d.immune_lists.slice(0, TOP_IMMUNE_LISTS).reverse();  // Plotly draws the first row at the bottom
  // genes counted in a list whose plain averages point the other way than the model (⚠)
  const flagged = (l) => d.genes.filter((g) => g.raw_disagrees && g.immune_lists.includes(l.id)).map((g) => g.symbol);
  $("immune-map-empty").hidden = lists.length > 0;
  $("immune-map-table").replaceChildren(simpleTable([
    ["Immune gene list", (l) => h("a", { href: l.url, target: "_blank", rel: "noopener" }, l.name)],
    ["Source", (l) => l.source],
    [`Higher in ${cmp}`, (l) => l.higher, true],
    [`Lower in ${cmp}`, (l) => l.lower, true],
    ["List genes tested", (l) => `${l.genes_tested} of ${l.genes_in_list}`, true],
    ["⚠ genes", (l) => flagged(l).join(", ") || "–"],
  ], d.immune_lists));
  if (!lists.length) return clearChart("immune-map-chart");

  const rows = lists.map((_, i) => i);
  const custom = lists.map((l) => {
    const odd = flagged(l);
    return [plotlyText(l.name), l.higher, l.lower, l.genes_tested, l.genes_in_list, l.id,
            odd.length ? `<br>⚠ ${plotlyText(odd.join(", "))}: plain averages point the other way (see gene panel)` : ""];
  });
  const hover = `<b>%{customdata[0]}</b><br>%{customdata[1]} higher, %{customdata[2]} lower in ${plotlyText(cmp)}`
    + "<br>%{customdata[3]} of the list's %{customdata[4]} genes were tested%{customdata[6]}<extra></extra>";
  const bars = (name, values, color) => ({
    type: "bar", orientation: "h", name, y: rows, x: values, customdata: custom, hovertemplate: hover,
    marker: { color, cornerradius: 4 }, text: values.map((v) => (v ? Math.abs(v) : "")),
    textposition: "outside", cliponaxis: false, textfont: { color: token("--ink-2") },
  });
  const most = Math.max(1, ...lists.map((l) => Math.max(l.higher, l.lower)));
  const step = [1, 2, 5, 10, 20, 50, 100].find((s) => most / s <= 4) ?? Math.ceil(most / 4);
  const ticks = [];
  for (let t = -Math.floor(most / step) * step; t <= most; t += step) ticks.push(t);
  const labelLength = $("immune-map-chart").clientWidth < 600 ? 24 : 46;

  draw("immune-map-chart", [
    bars(`Lower in ${cmp}`, lists.map((l) => -l.lower), token("--healthy")),
    bars(`Higher in ${cmp}`, lists.map((l) => l.higher), token("--icp")),
  ], layout({
    barmode: "relative",
    bargap: 0.35,
    height: lists.length * ROW_PX + 100,
    showlegend: true,
    xaxis: axis("Significant genes", { range: [-most * 1.3, most * 1.3], tickvals: ticks,
                                        ticktext: ticks.map((t) => String(Math.abs(t))), zeroline: true }),
    yaxis: axis("", { tickvals: rows, ticktext: lists.map((l) => shorten(l.name, labelLength) + (flagged(l).length ? " ⚠" : "")),
                      showgrid: false, tickfont: { color: token("--ink-2") } }),
  }), (event) => {
    const listId = event.points?.[0]?.customdata?.[5];
    if (!listId) return;
    state.filters.list = listId;
    $("filter-list").value = listId;
    renderFiltered();
    document.querySelector(".filters").scrollIntoView({ behavior: "smooth", block: "start" });
  });
}

/* --------------------------------------------------- immune cell types (atlas) */

/**
 * Forest plot of the atlas cell-type scores for one atlas section: a dot at the difference
 * between groups (in SD of the score), a line for its 95% confidence interval.
 */
function renderCellTypes() {
  const d = state.data, ref = d.reference_group, cmp = d.compare_group;
  const scored = d.cell_types;
  const lowestFdr = Math.min(...scored.map((c) => c.fdr ?? 1));
  const passing = scored.filter((c) => c.fdr != null && c.fdr < 0.05).length;
  $("cells-caption").textContent = "Immune cell populations and states from the Human Immune Cell Atlas. Each patient gets a "
    + "score per population: the average of its marker genes, each scaled across patients (using all "
    + `${fmt.int(d.genes_tested)} tested genes, not only the significant ones). The dot shows how much higher (orange, right) `
    + `or lower (blue, left) the score is in ${cmp} than in ${ref}, in standard deviations (SD), adjusted for fetal sex; the `
    + "line is the 95% confidence interval and filled dots have p < 0.05. "
    + (passing ? `${passing} of ${scored.length} populations pass FDR < 0.05. `
      : `None of the ${scored.length} scored populations passes FDR < 0.05 (lowest ${fmt.p(lowestFdr)}), so treat them as leads. `)
    + "Click a population to see each patient.";

  const sections = d.atlas_sections.filter((s) => scored.some((c) => c.section === s));
  if (!sections.includes(state.cellSection)) state.cellSection = sections[0] ?? null;
  $("cells-section").replaceChildren(...sections.map((s) => h("option", { value: s }, s)));
  $("cells-section").value = state.cellSection ?? "";
  renderCellTable(scored);

  const rows = scored.filter((c) => c.section === state.cellSection).sort((a, b) => a.effect - b.effect);
  if (!rows.some((c) => c.id === state.selectedCell)) state.selectedCell = [...rows].sort(byPvalue)[0]?.id ?? null;
  $("cells-empty").hidden = rows.length > 0;
  if (!rows.length) {
    clearChart("cells-chart");
    renderCellDetail();
    return;
  }

  const y = rows.map((_, i) => i);
  const colors = rows.map((c) => (c.effect > 0 ? token("--icp") : token("--healthy")));
  const reach = Math.max(0.5, ...rows.flatMap((c) => c.ci.map(Math.abs)));
  const label = (c) => {
    const name = plotlyText(shorten(c.name, 34));
    return c.id === state.selectedCell ? `<b>${name}</b>` : name;
  };
  draw("cells-chart", [{
    type: "scatter", mode: "markers", x: rows.map((c) => c.effect), y,
    error_x: { type: "data", symmetric: false, array: rows.map((c) => c.ci[1] - c.effect),
               arrayminus: rows.map((c) => c.effect - c.ci[0]), color: token("--ink-muted"), thickness: 2, width: 0 },
    // filled dot = p < 0.05, open dot otherwise
    marker: { size: 11, color: rows.map((c, i) => (c.pvalue < 0.05 ? colors[i] : token("--surface"))), line: { color: colors, width: 2 } },
    customdata: rows.map((c) => [plotlyText(c.name), fmt.sd(c.effect), `${fmt.sd(c.ci[0])} to ${fmt.sd(c.ci[1])}`, fmt.p(c.pvalue),
                                 fmt.p(c.fdr), c.markers.length, c.markers.length + c.markers_missing.length, c.id]),
    hovertemplate: `<b>%{customdata[0]}</b><br>%{customdata[1]} in ${plotlyText(cmp)} (95% CI %{customdata[2]})`
      + "<br>p = %{customdata[3]}, FDR = %{customdata[4]}<br>%{customdata[5]} of %{customdata[6]} markers tested<extra></extra>",
    cliponaxis: false,
  }], layout({
    height: rows.length * ROW_PX + 80,
    margin: { l: 10, r: 20, t: 10, b: 40 },
    xaxis: axis(`Score difference, ${cmp} vs ${ref} (SD)`, { range: [-reach * 1.1, reach * 1.1], zeroline: true }),
    yaxis: axis("", { tickvals: y, ticktext: rows.map(label), range: [-0.7, rows.length - 0.3],
                      tickfont: { color: token("--ink-2") } }),
  }), (event) => {
    const id = event.points?.[0]?.customdata?.[7];
    if (id) selectCellType(id);
  });
  renderCellDetail();
}

function selectCellType(id) {
  state.selectedCell = id;
  renderCellTypes();  // redraws the chart (to highlight the label) and the panel
}

/** Jump to a population from elsewhere on the page (e.g. a gene's atlas chips). */
function showCellType(id) {
  const population = state.data.atlas_populations[id];
  if (!population || !state.cellTypes.has(id)) return;
  state.cellSection = population.section;
  selectCellType(id);
  $("cells-card").scrollIntoView({ behavior: "smooth", block: "start" });
}

function renderCellDetail() {
  const d = state.data, cmp = d.compare_group;
  const cell = state.cellTypes.get(state.selectedCell);
  $("cell-detail").hidden = !cell;
  if (!cell) return;

  $("cell-name").textContent = cell.name + (cell.kind === "functional state" ? " (functional state)" : "");
  $("cell-summary").textContent = `${cell.section} · score ${fmt.sd(cell.effect)} in ${cmp} `
    + `(95% CI ${fmt.sd(cell.ci[0])} to ${fmt.sd(cell.ci[1])}) · p = ${fmt.p(cell.pvalue)} · FDR = ${fmt.p(cell.fdr)}`;
  drawPatientDots("cell-chart", cell.scores, { yTitle: "Marker score (SD)", format: (v) => v.toFixed(2), unit: "SD", zeroLine: true });

  const check = d.composition_check;
  $("cell-composition").textContent = check && cell.composition_r != null
    ? `Tissue check: across patients this score correlates with ${check.name.toLowerCase()} content `
      + `(${check.markers.join(", ")}; r = ${cell.composition_r.toFixed(2)}). Adjusted for it, the difference is `
      + `${fmt.sd(cell.effect_composition_adjusted)} (p = ${fmt.p(cell.pvalue_composition_adjusted)}).`
    : "";

  $("cell-markers").replaceChildren(
    ...cell.markers.map((m) => h("li", { class: "chip", title: `log2 fold change ${fmt.log2fc(m.log2fc)}, p = ${fmt.p(m.pvalue)}` },
      h("span", { class: `swatch ${m.log2fc > 0 ? "higher" : "lower"}` }), `${m.symbol} `,
      h("span", { class: "muted" }, fmt.fold(m.log2fc) + (m.pvalue < 0.05 ? " · p < 0.05" : "")))),
    ...cell.markers_missing.map((symbol) => h("li", { class: "chip muted" }, `${symbol} · not tested`)));
}

function renderCellTable(cells) {
  const check = state.data.composition_check;
  const columns = [
    ["§", (c) => c.section.split(".")[0], true],
    ["Population", (c) => c.name],
    ["Markers tested", (c) => `${c.markers.length} of ${c.markers.length + c.markers_missing.length}`, true],
    ["Difference", (c) => fmt.sd(c.effect), true],
    ["95% CI", (c) => `${c.ci[0].toFixed(2)} to ${c.ci[1].toFixed(2)}`, true],
    ["p-value", (c) => fmt.p(c.pvalue), true],
    ["FDR", (c) => fmt.p(c.fdr), true],
  ];
  if (check) {
    columns.push([`r with ${check.name}`, (c) => c.composition_r?.toFixed(2) ?? "–", true],
                 [`p adjusted for ${check.name}`, (c) => fmt.p(c.pvalue_composition_adjusted), true]);
  }
  const file = `${state.data.id}_cell_type_scores.csv`;
  $("cells-table").replaceChildren(simpleTable(columns, [...cells].sort(byPvalue)),
    h("p", {}, h("a", { href: `data/${file}`, download: file }, "Download all cell-type scores (CSV)")));
}

/* ------------------------------------------------------------------ filters */

function setUpFilters() {
  $("filter-search").addEventListener("input", (e) => { state.filters.search = e.target.value; renderFiltered(); });
  $("filter-direction").addEventListener("change", (e) => { state.filters.direction = e.target.value; renderFiltered(); });
  $("filter-immune").addEventListener("change", (e) => { state.filters.immuneOnly = e.target.checked; renderFiltered(); });
  $("filter-list").addEventListener("change", (e) => { state.filters.list = e.target.value; renderFiltered(); });
  $("filter-reset").addEventListener("click", () => { resetFilters(); renderFiltered(); });
}

function resetFilters() {
  state.filters = { search: "", direction: "all", immuneOnly: false, list: "" };
  $("filter-search").value = "";
  $("filter-direction").value = "all";
  $("filter-immune").checked = false;
  $("filter-list").value = "";
}

function renderFilterOptions() {
  const d = state.data;
  $("option-up").textContent = `Higher in ${d.compare_group}`;
  $("option-down").textContent = `Lower in ${d.compare_group}`;
  $("filter-list").replaceChildren(
    h("option", { value: "" }, "Any"),
    h("optgroup", { label: "Immune cell populations (atlas)" },
      d.cell_marker_changes.map((c) => h("option", { value: c.id }, `${shorten(c.name, 42)} (${c.genes.length})`))),
    d.immport ? h("optgroup", { label: "ImmPort immune gene lists" },
      d.immune_lists.map((l) => h("option", { value: l.id }, `${shorten(l.name, 42)} (${l.higher + l.lower})`))) : null);
  $("filter-list").value = state.filters.list;
}

function filteredGenes() {
  const { search, direction, immuneOnly, list } = state.filters;
  const query = search.trim().toLowerCase();
  return state.data.genes.filter((g) =>
    (direction === "all" || (direction === "up" ? g.log2fc > 0 : g.log2fc < 0))
    && (!immuneOnly || isImmune(g))
    && (!list || g.immune_lists.includes(list) || g.atlas_markers.includes(list))
    && (!query || [g.symbol, g.current_symbol, g.description, g.id].some((v) => v && v.toLowerCase().includes(query))));
}

/** Redraw everything below the filter row. */
function renderFiltered() {
  const genes = filteredGenes();
  $("filter-count").textContent = `${genes.length} of ${state.data.genes.length} genes`;
  if (!genes.some((g) => g.id === state.selected)) state.selected = [...genes].sort(byPvalue)[0]?.id ?? null;
  renderArrows(genes);
  renderTable(genes);
  renderDetail();
  renderHeatmap(genes);
}

/* -------------------------------------------------------------- arrow chart */

/**
 * One row per gene: a dot at the gene's level in the reference group (Healthy) and an arrow to
 * its level in the compare group (ICP). Right-pointing arrows = higher in ICP. The levels are the
 * model's estimates (other terms, such as fetal sex, held at the cohort average), so every arrow matches the
 * fold change in the table.
 */
function renderArrows(genes) {
  const d = state.data, ref = d.reference_group, cmp = d.compare_group;
  // the most significant genes, ordered so the biggest increase is at the top
  const top = [...genes].sort(byPvalue).slice(0, TOP_ARROWS).sort((a, b) => a.log2fc - b.log2fc);
  $("arrows-empty").hidden = top.length > 0;
  $("arrows-caption").textContent = `Each row is a gene. The blue dot is its expression level in ${ref} patients and the `
    + `arrow points to its level in ${cmp} patients (orange): right means higher in ${cmp}, left means lower. Levels are `
    + `the model's estimates${adjustedForSex() ? " with fetal sex held at the cohort average" : ""}, so each arrow is the fold change in the `
    + `table; click a gene to see every patient's actual value. Showing the ${top.length} most significant of `
    + `${genes.length} genes; ◆ marks immune-related genes, ⚠ genes whose plain averages point the other way.`;
  if (!top.length) return clearChart("arrows-chart");

  const rows = top.map((_, i) => i);
  const atLeast = (v) => Math.max(v, 0.5);  // a log axis cannot show 0; a level of 0 is drawn at 0.5
  const before = top.map((g) => g.group_estimates[ref]);
  const after = top.map((g) => g.group_estimates[cmp]);
  const dots = (group, values, color) => ({
    type: "scatter", mode: "markers", name: `${group} level`, y: rows, x: values.map(atLeast),
    customdata: top.map((g, i) => [plotlyText(g.symbol), fmt.count(values[i]), fmt.fold(g.log2fc), fmt.p(g.pvalue), g.id]),
    hovertemplate: `<b>%{customdata[0]}</b><br>${plotlyText(group)} level: %{customdata[1]} (model estimate)`
      + `<br>%{customdata[2]} in ${plotlyText(cmp)} (p = %{customdata[3]})<extra></extra>`,
    marker: { color, size: 11, line: { color: token("--surface"), width: 2 } },
    cliponaxis: false,
  });
  // On a log axis, Plotly annotation positions are given as log10 of the value.
  const arrows = top.map((g, i) => ({
    xref: "x", yref: "y", axref: "x", ayref: "y",
    x: Math.log10(atLeast(after[i])), y: i, ax: Math.log10(atLeast(before[i])), ay: i,
    showarrow: true, arrowhead: 2, arrowsize: 1, arrowwidth: 2, arrowcolor: token("--ink-2"),
    standoff: 8, startstandoff: 8, text: "",
  }));

  draw("arrows-chart", [dots(ref, before, token("--healthy")), dots(cmp, after, token("--icp"))], layout({
    height: top.length * ROW_PX + 100,
    showlegend: true,
    annotations: arrows,
    xaxis: axis(`${units().axis} (log scale)`, { type: "log", dtick: 1, tickformat: "~s" }),
    yaxis: axis("", { tickvals: rows, ticktext: top.map(geneLabel), range: [-0.7, top.length - 0.3],
                      tickfont: { color: token("--ink-2") } }),
  }), (event) => {
    const id = event.points?.[0]?.customdata?.[4];
    if (id) selectGene(id);
  });
}

/* --------------------------------------------------------------- gene table */

const COLUMNS = [
  { key: "symbol", label: "Gene", value: (g) => g.symbol.toLowerCase() },
  { key: "log2fc", label: "Change", numeric: true, descendingFirst: true, value: (g) => g.log2fc },
  { key: "pvalue", label: "p-value", numeric: true, value: (g) => g.pvalue },
  { key: "padj", label: "FDR", numeric: true, value: (g) => g.padj ?? Infinity },
  { key: "base_mean", label: "Mean", numeric: true, descendingFirst: true, value: (g) => g.base_mean },
  { key: "chromosome", label: "Chr", value: (g) => chromosomeRank(g.chromosome) },
  { key: "immune", label: "Immune", numeric: true, descendingFirst: true,
    value: (g) => (isImmune(g) ? g.immune_lists.length + g.atlas_markers.length : -1) },
];

function sortBy(key) {
  const column = COLUMNS.find((c) => c.key === key);
  state.sort = state.sort.key === key
    ? { key, ascending: !state.sort.ascending }
    : { key, ascending: !column.descendingFirst };
  renderTable(filteredGenes());
}

function renderTable(genes) {
  const cmp = state.data.compare_group;
  document.querySelector("#gene-table thead tr").replaceChildren(...COLUMNS.map((c) => {
    const sorted = state.sort.key === c.key ? (state.sort.ascending ? "ascending" : "descending") : "none";
    return h("th", { scope: "col", class: c.numeric ? "num" : null, "aria-sort": sorted },
      h("button", { type: "button", onclick: () => sortBy(c.key) },
        c.key === "log2fc" ? `Change in ${cmp}` : c.key === "base_mean" ? units().mean : c.label));
  }));

  const column = COLUMNS.find((c) => c.key === state.sort.key);
  const direction = state.sort.ascending ? 1 : -1;
  const sorted = [...genes].sort((a, b) => {
    const x = column.value(a), y = column.value(b);
    return (x < y ? -1 : x > y ? 1 : 0) * direction;
  });
  const shown = state.showAllRows ? sorted : sorted.slice(0, TABLE_ROWS);
  document.querySelector("#gene-table tbody").replaceChildren(...shown.map(geneRow));
  $("table-empty").hidden = genes.length > 0;
  const more = $("table-more");
  more.hidden = sorted.length <= TABLE_ROWS;
  more.textContent = state.showAllRows ? `Show the first ${TABLE_ROWS} only` : `Show all ${sorted.length} genes`;
  more.onclick = () => { state.showAllRows = !state.showAllRows; renderTable(filteredGenes()); };
  const csv = $("csv-link");
  csv.href = `data/${state.data.id}_significant_genes.csv`;
  csv.setAttribute("download", `${state.data.id}_significant_genes.csv`);
}

function geneRow(g) {
  const choose = () => selectGene(g.id);
  const annotations = g.immune_lists.length + g.atlas_markers.length;
  return h("tr", {
    "data-id": g.id, tabindex: 0, "aria-selected": String(g.id === state.selected), onclick: choose,
    onkeydown: (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); choose(); } },
  },
  h("td", {}, h("div", { class: "gene-name" }, g.symbol), g.description ? h("div", { class: "gene-desc" }, g.description) : null),
  h("td", { class: "num" }, h("span", { class: `swatch ${g.log2fc > 0 ? "higher" : "lower"}` }), fmt.fold(g.log2fc),
    g.raw_disagrees ? h("span", { title: "Plain averages point the other way; see the gene panel" }, " ⚠") : null,
    h("div", { class: "muted" }, `log2 ${fmt.log2fc(g.log2fc)}`)),
  h("td", { class: "num" }, fmt.p(g.pvalue)),
  h("td", { class: "num" }, fmt.p(g.padj)),
  h("td", { class: "num" }, fmt.count(g.base_mean)),
  h("td", {}, g.chromosome ?? "–"),
  h("td", { class: "num" }, isImmune(g)
    ? h("span", { class: "badge", title: `${g.immune_lists.length} ImmPort lists, ${g.atlas_markers.length} atlas populations` }, `◆ ${annotations}`)
    : h("span", { class: "muted" }, "–")));
}

function selectGene(id) {
  state.selected = id;
  for (const row of document.querySelectorAll("#gene-table tbody tr")) {
    row.setAttribute("aria-selected", String(row.dataset.id === id));
  }
  renderDetail();
}

/* -------------------------------------------------------------- gene detail */

function renderDetail() {
  const d = state.data, ref = d.reference_group, cmp = d.compare_group;
  const gene = d.genes.find((g) => g.id === state.selected);
  $("detail-placeholder").hidden = Boolean(gene);
  $("detail-body").hidden = !gene;
  if (!gene) return;

  $("detail-symbol").textContent = gene.symbol;
  $("detail-renamed").textContent = gene.current_symbol ? `now named ${gene.current_symbol}` : "";
  $("detail-immune").textContent = isImmune(gene) ? "◆ Immune-related" : "Not immune-annotated";
  $("detail-description").textContent = gene.description ?? "";
  const stats = [
    [`Change in ${cmp}`, fmt.fold(gene.log2fc)],
    ["log2 fold change", `${fmt.log2fc(gene.log2fc)} ± ${gene.lfc_se?.toFixed(2) ?? "–"}`],
    ["p-value", fmt.p(gene.pvalue)],
    ["FDR (padj)", fmt.p(gene.padj)],
    [`${ref} level (model)`, fmt.count(gene.group_estimates[ref])],
    [`${cmp} level (model)`, fmt.count(gene.group_estimates[cmp])],
    ["Chromosome", gene.chromosome ?? "–"],
    ["Type", gene.biotype?.replaceAll("_", " ") ?? "–"],
  ];
  $("detail-stats").replaceChildren(...stats.map(([term, value]) => h("div", {}, h("dt", {}, term), h("dd", {}, value))));

  // the change within each fetal sex, from plain averages
  const bySex = [["F", "female"], ["M", "male"]].filter(([sex]) => gene.sex_log2fc?.[sex] != null).map(([sex, label]) => {
    const n = d.sex_counts?.[sex] ?? {};
    return `${label} placentas (${n[ref] ?? 0} ${ref} vs ${n[cmp] ?? 0} ${cmp}): ${fmt.fold(gene.sex_log2fc[sex])}`;
  });
  $("detail-sex").textContent = bySex.length ? `By fetal sex, plain averages: ${bySex.join(" · ")}.` : "";
  $("detail-warning").hidden = !gene.raw_disagrees;
  $("detail-warning").replaceChildren(...(gene.raw_disagrees ? [
    h("strong", {}, "⚠ Plain averages point the other way. "),
    `Averaged over all patients, ${gene.symbol} reads ${fmt.count(gene.plain_means[cmp])} in ${cmp} vs `
    + `${fmt.count(gene.plain_means[ref])} in ${ref}, which looks ${gene.plain_log2fc > 0 ? "higher" : "lower"} in ${cmp}. `
    + `The change differs between female and male placentas and the groups have different fetal-sex mixes, so the `
    + (adjustedForSex()
      ? `model compares like with like by adjusting for fetal sex: ${fmt.fold(gene.log2fc)} in ${cmp}. Treat this gene with caution.`
      : `model compares typical (log-scale) values, which one very high sample sways less: ${fmt.fold(gene.log2fc)} in ${cmp}. `
        + "Treat this gene with caution."),
  ] : []));

  $("detail-immport").hidden = !d.immport;
  $("detail-lists").replaceChildren(...(gene.immune_lists.length
    ? gene.immune_lists.map((id) => {
      const list = state.lists.get(id) ?? { name: id, url: null, source: "" };
      return h("li", {}, h("a", { href: list.url, target: "_blank", rel: "noopener" }, list.name, h("span", { class: "source" }, list.source)));
    })
    : [h("li", { class: "muted" }, "None")]));

  $("detail-cells").replaceChildren(...(gene.atlas_markers.length
    ? gene.atlas_markers.map((id) => {
      const population = d.atlas_populations[id] ?? { name: id, section: "", kind: "" };
      const scored = state.cellTypes.has(id);
      return h("li", {}, h("button", {
        type: "button", class: "chip-button", disabled: !scored, onclick: () => showCellType(id),
        title: scored ? "Show this population's scores" : "Not scored: fewer than 2 of its markers were tested",
      }, population.name, h("span", { class: "source" },
        population.kind === "functional state" ? "state" : `§${population.section.split(".")[0]}`)));
    })
    : [h("li", { class: "muted" }, "None")]));

  const symbol = encodeURIComponent(gene.current_symbol ?? gene.symbol);
  const links = [["GeneCards", `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${symbol}`],
                 ["NCBI Gene", `https://www.ncbi.nlm.nih.gov/gene/?term=${symbol}%5Bsym%5D+AND+human%5Borgn%5D`]];
  if (gene.id.startsWith("ENSG")) links.push(["Ensembl", `https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${encodeURIComponent(gene.id.split(".")[0])}`]);
  $("detail-links").replaceChildren("Look up: ", ...links.flatMap(([name, url], i) =>
    [i ? " · " : "", h("a", { href: url, target: "_blank", rel: "noopener" }, name)]));

  $("detail-measure").textContent = units().detail;
  drawPatientDots("detail-chart", gene.normalized, { yTitle: units().axis, format: fmt.count, unit: units().value });
  $("detail-table").replaceChildren(simpleTable([
    ["Sample", (r) => r.s.id], ["Group", (r) => r.s.group], ["Fetal sex", (r) => sexName(r.s.fetal_sex)],
    ["Normalized count", (r) => fmt.count(gene.normalized[r.i]), true],
  ], d.samples.map((s, i) => ({ s, i }))));
}

/* ------------------------------------------------------------------ heatmap */

function renderHeatmap(genes) {
  const d = state.data, samples = d.samples;
  const top = [...genes].sort(byPvalue).slice(0, TOP_HEATMAP).sort((a, b) => b.log2fc - a.log2fc);
  $("heatmap-empty").hidden = top.length > 0;
  $("heatmap-caption").textContent = "Each row is a gene, each column a patient. Color shows whether the patient's expression is "
    + "above (orange) or below (blue) that gene's average across all patients, in standard deviations (z-score of "
    + `variance-stabilized counts). Showing the ${top.length} most significant of ${genes.length} genes, highest increase in `
    + `${d.compare_group} first.`;
  if (!top.length) return clearChart("heatmap-chart");

  const rows = top.map((_, i) => i);
  const scale = [[0, token("--div-low-3")], [1 / 6, token("--div-low-2")], [2 / 6, token("--div-low-1")], [0.5, token("--div-mid")],
                 [4 / 6, token("--div-high-1")], [5 / 6, token("--div-high-2")], [1, token("--div-high-3")]];
  const height = top.length * 18 + 130;
  // a line between the groups, and the group names above the columns
  const groups = [...new Set(samples.map((s) => s.group))];
  const firstColumn = (group) => samples.findIndex((s) => s.group === group);
  const lastColumn = (group) => samples.findLastIndex((s) => s.group === group);
  const dividers = groups.slice(1).map((group) => ({
    type: "line", xref: "x", yref: "paper", x0: firstColumn(group) - 0.5, x1: firstColumn(group) - 0.5, y0: 0, y1: 1,
    line: { color: token("--surface"), width: 6 } }));
  const groupNames = groups.map((group) => ({
    xref: "x", yref: "paper", x: (firstColumn(group) + lastColumn(group)) / 2, y: 1, yanchor: "bottom", showarrow: false,
    text: plotlyText(group), font: { color: token("--ink-1"), size: 13 } }));

  draw("heatmap-chart", [{
    type: "heatmap", z: top.map((g) => g.z), x: samples.map((_, i) => i), y: rows,
    zmin: -2.5, zmax: 2.5, zmid: 0, colorscale: scale, xgap: 2, ygap: 2,
    customdata: top.map((g) => samples.map((s, i) => [plotlyText(g.symbol), plotlyText(s.id), plotlyText(s.group), fmt.count(g.normalized[i])])),
    hovertemplate: `<b>%{customdata[0]}</b> · %{customdata[1]} (%{customdata[2]})<br>%{customdata[3]} ${units().value}, z = %{z:.2f}<extra></extra>`,
    colorbar: { title: { text: "z-score", side: "right", font: { color: token("--ink-2") } }, thickness: 10,
                outlinewidth: 0, tickfont: { color: token("--ink-muted") }, len: Math.min(1, 320 / height) },
  }], layout({
    height,
    margin: { l: 10, r: 10, t: 30, b: 10 },
    shapes: dividers,
    annotations: groupNames,
    xaxis: axis("", { tickvals: samples.map((_, i) => i), ticktext: samples.map((s) => s.id), showgrid: false, zeroline: false }),
    yaxis: axis("", { tickvals: rows, ticktext: top.map(geneLabel), autorange: "reversed", showgrid: false,
                      tickfont: { color: token("--ink-2") } }),
  }), (event) => {
    const point = event.points?.[0];
    if (point) selectGene(top[point.y].id);
  });
}

/* ----------------------------------------------------------- quality checks */

function renderQualityChecks() {
  const d = state.data;
  const sizes = countBy(d.samples, (s) => s.group);
  const [pc1, pc2] = d.pca_variance.map((v) => Math.round(v * 100));
  $("qc-caption").textContent = "Principal component analysis of all patients, using the 500 most variable of the "
    + `${fmt.int(d.genes_tested)} tested genes (not only the significant ones). Patients close together have similar overall `
    + `expression; if the groups do not separate, ${d.compare_group} has only a subtle overall effect. With `
    + `${[...sizes.values()].join(" and ")} patients per group, a single unusual sample can sway the results. `
    + "● female fetus, ◆ male fetus.";

  const traces = [...sizes.keys()].map((group) => {
    const members = d.samples.filter((s) => s.group === group);
    return {
      type: "scatter", mode: "markers", name: group, x: members.map((s) => s.pc1), y: members.map((s) => s.pc2),
      customdata: members.map((s) => [plotlyText(s.id), sexName(s.fetal_sex)]),
      hovertemplate: "<b>%{customdata[0]}</b><br>fetal sex: %{customdata[1]}<extra></extra>",
      marker: { color: groupColor(group), size: 12, symbol: members.map((s) => (s.fetal_sex === "M" ? "diamond" : "circle")),
                line: { color: token("--surface"), width: 2 } },
    };
  });
  draw("pca-chart", traces, layout({
    height: 360, showlegend: true,
    xaxis: axis(`PC1 (${pc1}% of variation)`, { zeroline: false }),
    yaxis: axis(`PC2 (${pc2}% of variation)`, { zeroline: false }),
  }));
  $("samples-table").replaceChildren(simpleTable([
    ["Sample", (s) => s.id], ["Group", (s) => s.group], ["Fetal sex", (s) => sexName(s.fetal_sex)],
    ["Sequenced reads", (s) => (s.library_size == null ? "–" : fmt.int(s.library_size)), true], ["PC1", (s) => s.pc1?.toFixed(1), true], ["PC2", (s) => s.pc2?.toFixed(1), true],
  ], d.samples));
}

/* ------------------------------------------------------- chart/table toggle */

function setUpTableToggles() {
  for (const button of document.querySelectorAll(".table-toggle")) {
    button.addEventListener("click", () => {
      const view = $(button.dataset.chart), table = $(button.dataset.table);
      const showTable = table.hidden;
      table.hidden = !showTable;
      view.hidden = showTable;
      button.textContent = showTable ? "Show as chart" : "Show as table";
      if (!showTable) {
        // charts drawn while hidden need to re-measure their width
        const plots = view.classList.contains("js-plotly-plot") ? [view] : [...view.querySelectorAll(".js-plotly-plot")];
        plots.forEach((plot) => Plotly.Plots.resize(plot));
      }
    });
  }
}

start();
