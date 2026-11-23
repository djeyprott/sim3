// The viewer polls the local simulator for its state and new events.

const $ = (id) => document.getElementById(id);
let state = null;
let nextEventId = 0;

// Simulated days per real second.
const SPEED_NAMES = {
  1: "1 day / sec",
  30: "1 month / sec",
  365: "1 year / sec",
  3652: "10 years / sec",
  36524: "100 years / sec",
  0: "Maximum",
};

// How each kind of event reads in the chronicle. Later phases add many more.
const EVENT_TEXT = {
  simulation_started: () => "The simulation began.",
};

async function api(path, body) {
  const opts = body
    ? { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) }
    : {};
  const res = await fetch(path, opts);
  if (!res.ok) throw new Error(`${path}: ${res.status}`);
  return res.json();
}

function render() {
  if (!state) return;
  $("date").textContent = state.date_label;
  const years = state.years_elapsed;
  $("elapsed").textContent = `${years.toLocaleString()} ${years === 1 ? "year" : "years"} simulated`;
  const ys = state.measured_years_per_second;
  $("rate").textContent = state.finished ? "run complete" : ys > 0.05 ? `running at ${ys >= 10 ? Math.round(ys) : ys.toFixed(1)} years/sec` : "";
  $("toggle").textContent = state.paused ? "Play" : "Pause";

  const select = $("speed");
  if (select.options.length === 0) {
    for (const s of state.speeds) {
      const opt = document.createElement("option");
      opt.value = s;
      opt.textContent = SPEED_NAMES[s] ?? `${s} days / sec`;
      select.appendChild(opt);
    }
  }
  if (document.activeElement !== select) select.value = state.days_per_second;
  if (document.activeElement !== $("run")) $("run").value = state.scenario;
}

async function poll() {
  try {
    state = await api("/api/state");
    $("connection").textContent = "connected";
    $("connection").className = "status ok";
    render();
    if (state.event_count > nextEventId) {
      const events = await api(`/api/events?from=${nextEventId}&limit=200`);
      const list = $("events");
      for (const e of events) {
        const li = document.createElement("li");
        const when = document.createElement("span");
        when.className = "when";
        when.textContent = e.date_label;
        li.appendChild(when);
        li.appendChild(document.createTextNode((EVENT_TEXT[e.type] ?? (() => e.type))(e)));
        list.prepend(li);
        nextEventId = e.id + 1;
      }
    }
  } catch {
    $("connection").textContent = "not connected — is sim3 running?";
    $("connection").className = "status bad";
  }
}

$("toggle").addEventListener("click", async () => {
  await api("/api/control", { action: state?.paused ? "resume" : "pause" });
  poll();
});
for (const unit of ["day", "month", "year"]) {
  $(`step-${unit}`).addEventListener("click", async () => {
    await api("/api/control", { action: "step", unit, count: 1 });
    setTimeout(poll, 60);
  });
}
$("run").addEventListener("change", async (ev) => {
  const label = ev.target.selectedOptions[0].textContent;
  if (!confirm(`Start a new run: ${label}? The current run will be discarded.`)) {
    ev.target.value = state.scenario;
    return;
  }
  await api("/api/run", { scenario: ev.target.value });
  nextEventId = 0;
  $("events").innerHTML = "";
  $("cell-chart").hidden = true;
  poll();
});
$("speed").addEventListener("change", async (ev) => {
  await api("/api/control", { action: "speed", days_per_second: Number(ev.target.value) });
  poll();
});

poll();
setInterval(poll, 250);

// ---- map ----

const fmt = (m) => `${Math.round(m).toLocaleString()} m`;
const fmtC = (t) => `${t >= 0 ? "" : "−"}${Math.abs(t).toFixed(1)} °C`;
const SURFACE = ["Sea", "Land", "Below sea level, cut off from the ocean"];
const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
const MONTH_NAMES = ["January", "February", "March", "April", "May", "June", "July", "August", "September",
  "October", "November", "December"];

function coords(h) {
  const ns = h.lat >= 0 ? "N" : "S", ew = h.lon >= 0 ? "E" : "W";
  return `${Math.abs(h.lat).toFixed(2)}°${ns}, ${Math.abs(h.lon).toFixed(2)}°${ew}`;
}

function hoverHtml(h) {
  const rows = [
    `<div class="coords">${coords(h)} · ${h.epochLabel}</div>`,
    `<div>${h.ice > 0 && h.kind !== 0 ? "Ice" : SURFACE[h.kind]}</div>`,
    h.kind === 0
      ? `<div>Average depth ${fmt(-h.elev)} (deepest ${fmt(-h.elevMin)})</div>`
      : `<div>Average height ${fmt(h.elev)} (lowest ${fmt(h.elevMin)}, highest ${fmt(h.elevMax)})</div>`,
    `<div>Land: ${Math.round(h.landShare * 100)}% of the cell</div>`,
  ];
  if (h.ice > 0) rows.push(`<div>Ice thickness ${fmt(h.ice)}</div>`);
  if (h.temp !== null && h.temp !== undefined) {
    rows.push(`<div>${h.tempTitle}: ${h.tempMode === "diff" ? (h.temp >= 0 ? "+" : "−") + Math.abs(h.temp).toFixed(1) + " °C" : fmtC(h.temp)}</div>`);
  }
  rows.push(`<div class="coords">Click for this place's climate</div>`);
  return rows.join("");
}

function optionButtons(container, items, current, onPick) {
  container.innerHTML = "";
  for (const it of items) {
    const b = document.createElement("button");
    b.type = "button";
    b.textContent = it.label;
    b.dataset.id = it.id;
    b.setAttribute("aria-pressed", String(it.id === current));
    b.addEventListener("click", () => {
      for (const o of container.querySelectorAll("button")) o.setAttribute("aria-pressed", "false");
      b.setAttribute("aria-pressed", "true");
      onPick(it.id);
    });
    container.appendChild(b);
  }
}

// ---- small SVG line charts ----

function lineChart({ series, xLabels, xTicks, width = 340, height = 170, yUnit = "°C" }) {
  const pad = { l: 34, r: 8, t: 8, b: 20 };
  const all = series.flatMap((s) => s.values.map((p) => p[1])).filter(Number.isFinite);
  if (!all.length) return "<p>No data yet.</p>";
  let lo = Math.min(...all), hi = Math.max(...all);
  if (hi - lo < 2) { lo -= 1; hi += 1; }
  const step = [1, 2, 5, 10, 20][[1, 2, 5, 10, 20].findIndex((s) => (hi - lo) / s <= 6)] ?? 50;
  lo = Math.floor(lo / step) * step; hi = Math.ceil(hi / step) * step;
  const xs = series.flatMap((s) => s.values.map((p) => p[0]));
  const x0 = Math.min(...xs), x1 = Math.max(...xs) === x0 ? x0 + 1 : Math.max(...xs);
  const X = (x) => pad.l + (x - x0) / (x1 - x0) * (width - pad.l - pad.r);
  const Y = (y) => pad.t + (hi - y) / (hi - lo) * (height - pad.t - pad.b);
  let out = `<svg viewBox="0 0 ${width} ${height}" role="img">`;
  for (let y = lo; y <= hi + 1e-9; y += step) {
    out += `<line x1="${pad.l}" x2="${width - pad.r}" y1="${Y(y)}" y2="${Y(y)}" stroke="#3a342a" stroke-width="1"/>`;
    out += `<text x="${pad.l - 4}" y="${Y(y) + 3}" text-anchor="end">${y}${y === hi ? " " + yUnit : ""}</text>`;
  }
  for (const [x, label] of xTicks) out += `<text x="${X(x)}" y="${height - 6}" text-anchor="middle">${label}</text>`;
  for (const s of series) {
    const pts = s.values.filter((p) => Number.isFinite(p[1])).map((p) => `${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join(" ");
    out += `<polyline fill="none" stroke="${s.color}" stroke-width="2" ${s.dash ? 'stroke-dasharray="5 4"' : ""} points="${pts}"/>`;
  }
  return out + "</svg>";
}

async function showCellChart(h) {
  const d = await api(`/api/climate/cell?cell=${h.cell}`);
  const months = (arr) => arr.map((v, i) => [i, v]);
  const series = [];
  const key = [];
  if (d.window) {
    series.push({ values: months(d.window), color: "#d9a441" });
    key.push(`<span class="model">Model, 1981–2010</span>`);
  } else {
    series.push({ values: months(d.normal), color: "#d9a441" });
    key.push(`<span class="model">Recent ${Math.min(d.years, 30)} years (model)</span>`);
  }
  if (d.observed) {
    series.push({ values: months(d.observed), color: "#8fb3d9", dash: true });
    key.push(`<span class="obs">Observed, 1981–2010</span>`);
  }
  const box = $("cell-chart");
  box.innerHTML = `<div class="head"><strong>Climate at ${coords(h)}</strong><button type="button" id="cell-close">×</button></div>` +
    `<div class="coords">Average temperature by month · height ${fmt(Math.max(0, h.elev))}</div>` +
    lineChart({ series, xTicks: MONTHS.map((m, i) => [i, m[0]]) }) + `<div class="key">${key.join("")}</div>`;
  box.hidden = false;
  $("cell-close").addEventListener("click", () => { box.hidden = true; });
}

async function refreshWorldChart() {
  if ($("world-chart").hidden) return;
  const d = await api("/api/climate/global");
  const pts = d.points;
  const yr = (y) => (y <= 0 ? `${1 - Math.round(y)} BCE` : `${Math.round(y)} CE`);
  const ticks = pts.length ? [[pts[0][0], yr(pts[0][0])], [pts[pts.length - 1][0], yr(pts[pts.length - 1][0])]] : [];
  $("world-chart").innerHTML = `<div class="head"><strong>World average temperature</strong></div>` +
    `<div class="coords">Yearly average of the whole surface, land and sea</div>` +
    lineChart({ series: [{ values: pts, color: "#d9a441" }], xTicks: ticks, width: 400, height: 160 });
}

function checkHtml(d) {
  const [a, b] = d.window;
  if (!d.scores) {
    const done = state?.check_years ?? 0;
    return `<p>The model settles from 1931 to 1980, then averages ${a}–${b} and compares with observations.</p>` +
      `<p>Averaged so far: ${done} of ${b - a + 1} years. Press Play (fastest speed is fine).</p>`;
  }
  const s = d.scores, pm = d.pass_marks;
  const row = (label, value, ok) => `<tr><td>${label}</td><td class="num">${value}</td>` +
    `<td class="num ${ok === undefined ? "" : ok ? "pass" : "fail"}">${ok === undefined ? "" : ok ? "pass" : "miss"}</td></tr>`;
  return `<p>Averaged ${s.years_averaged} of ${b - a + 1} years (${a}–${b}).</p><table>` +
    row(`World average (observed ${d.observed_global_mean_c} ± ${d.observed_global_uncertainty_c} °C)`, fmtC(s.global_mean_model_c),
      s.pass_global_mean) +
    row(`Typical land error, yearly average (mark ${pm.annual_typical_error_c} °C)`, `${s.annual_typical_error_c.toFixed(1)} °C`,
      s.pass_annual) +
    row(`Typical land error, summer–winter swing (mark ${pm.seasonal_swing_error_c} °C)`, `${s.seasonal_swing_error_c.toFixed(1)} °C`,
      s.pass_seasonal) +
    row("Land on average too warm (+) or cold (−)", `${s.land_bias_c >= 0 ? "+" : "−"}${Math.abs(s.land_bias_c).toFixed(1)} °C`) +
    `</table><p class="note">At this step only the world average is expected to pass. The land marks need rain, winds and sea ice (Steps 3–4). ` +
    `Compare maps with the Temperature layer's "Observed" and "Model − observed" views.</p>`;
}

async function refreshCheck() {
  const box = $("check");
  if (!state || state.scenario !== "present_day_check") { box.hidden = true; return; }
  box.hidden = false;
  $("check-body").innerHTML = checkHtml(await api("/api/climate/check"));
}

const temp = { view: "normal", month: null };
let lastTempTick = -1;

function tempOptions(map) {
  const box = $("temp-options");
  box.hidden = map.layer() !== "temperature";
  if (box.hidden) return;
  const views = [{ id: "current", label: "This month" }, { id: "normal", label: "Recent average" }];
  if (state?.scenario === "present_day_check") {
    views.push({ id: "observed", label: "Observed" }, { id: "difference", label: "Model − observed" });
  }
  if (!views.some((v) => v.id === temp.view)) temp.view = "normal";
  optionButtons(box, views, temp.view, (id) => { temp.view = id; tempOptions(map); refreshTemperature(map, true); });
  if (temp.view !== "current") {
    const sel = document.createElement("select");
    sel.innerHTML = `<option value="">Whole year</option>` + MONTH_NAMES.map((m, i) => `<option value="${i}">${m}</option>`).join("");
    sel.value = temp.month ?? "";
    sel.addEventListener("change", () => { temp.month = sel.value === "" ? null : Number(sel.value); refreshTemperature(map, true); });
    box.appendChild(sel);
  }
}

async function refreshTemperature(map, force) {
  if (map.layer() !== "temperature" || !state) return;
  const fixed = temp.view === "observed";
  if (!force && (fixed || state.tick === lastTempTick)) return;
  lastTempTick = state.tick;
  const q = `/api/climate/map?view=${temp.view}${temp.view !== "current" && temp.month !== null ? `&month=${temp.month}` : ""}`;
  const buf = await fetch(q).then((r) => r.arrayBuffer());
  const when = temp.view === "current" ? state.date_label.replace(/^\d+ /, "")
    : temp.month === null ? "whole year" : MONTH_NAMES[temp.month];
  const titles = {
    current: `Temperature, ${when} (°C)`,
    normal: `Average temperature, ${when}, recent years (°C)`,
    observed: `Observed average, ${when}, 1981–2010 (°C)`,
    difference: `Model minus observed, ${when} (°C)`,
  };
  map.setTemperature(new Int16Array(buf), temp.view === "difference" ? "diff" : "temp", titles[temp.view]);
  $("legend").innerHTML = map.legendHtml();
}

WorldMap.create($("map"), (h) => {
  const box = $("hover");
  box.hidden = !h;
  if (h) box.innerHTML = hoverHtml(h);
}, (h) => showCellChart(h).catch(() => {})).then((map) => {
  $("map-loading").remove();
  const showLegend = () => { $("legend").innerHTML = map.legendHtml(); };
  const epochButtons = () => optionButtons($("epochs"), map.epochs.map((e) => ({ id: e.index, label: e.label })), map.epoch(),
    (i) => map.setEpoch(i));
  epochButtons();
  optionButtons($("layers"), map.layers.map((l) => ({ id: l.id, label: l.name })), "terrain", (id) => {
    map.setLayer(id);
    tempOptions(map);
    showLegend();
    refreshTemperature(map, true);
  });
  showLegend();
  $("world-toggle").addEventListener("click", (e) => {
    const box = $("world-chart");
    box.hidden = !box.hidden;
    e.target.setAttribute("aria-pressed", String(!box.hidden));
    refreshWorldChart();
  });
  let lastScenario = null;
  setInterval(() => {
    if (!state) return;
    if (state.scenario !== lastScenario) {
      lastScenario = state.scenario;
      map.setEpoch(state.epoch);
      epochButtons();
      tempOptions(map);
      refreshTemperature(map, true);
    }
    refreshTemperature(map, false);
  }, 1000);
  setInterval(() => { refreshWorldChart().catch(() => {}); refreshCheck().catch(() => {}); }, 2000);
}).catch((err) => {
  $("map-loading").textContent = `The map could not be drawn: ${err.message}`;
});
