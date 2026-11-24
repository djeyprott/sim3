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
let ZONE_NAMES = [];
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
  if (h.field !== null && h.field !== undefined) {
    const v = h.field;
    const text = {
      temp: () => fmtC(v / 100),
      diff: () => `${v >= 0 ? "+" : "−"}${Math.abs(v / 100).toFixed(1)} °C`,
      rain: () => `${Math.round(v).toLocaleString()} mm`,
      ratio: () => { const r = Math.pow(10, v / 100); return r >= 1 ? `×${r.toFixed(1)}` : `÷${(1 / r).toFixed(1)}`; },
      snow: () => `${Math.round(v).toLocaleString()} mm of water`,
      seaice: () => (v > 0 ? `${(v / 100).toFixed(1)} m thick` : "open sea"),
      zone: () => { const z = ZONE_NAMES[v]; return z ? `${z.code} · ${z.name}` : "—"; },
    }[h.fieldMode] ?? (() => String(v));
    rows.push(`<div>${h.fieldTitle}: ${text()}</div>`);
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

// Monthly rainfall as bars (model) with observed marks, sharing the
// month axis of the temperature chart above it.
function rainBars({ model, observed, width = 340, height = 110 }) {
  const pad = { l: 34, r: 8, t: 8, b: 18 };
  const all = [...model, ...(observed ?? [])].filter(Number.isFinite);
  let hi = Math.max(10, ...all);
  const step = [10, 20, 50, 100, 200, 500][[10, 20, 50, 100, 200, 500].findIndex((s) => hi / s <= 4)] ?? 1000;
  hi = Math.ceil(hi / step) * step;
  const X = (i) => pad.l + i / 11 * (width - pad.l - pad.r);
  const Y = (v) => pad.t + (hi - v) / hi * (height - pad.t - pad.b);
  const bw = (width - pad.l - pad.r) / 11 * 0.55;
  let out = `<svg viewBox="0 0 ${width} ${height}" role="img">`;
  for (let v = 0; v <= hi + 1e-9; v += step) {
    out += `<line x1="${pad.l}" x2="${width - pad.r}" y1="${Y(v)}" y2="${Y(v)}" stroke="#3a342a" stroke-width="1"/>`;
    out += `<text x="${pad.l - 4}" y="${Y(v) + 3}" text-anchor="end">${v}${v === hi ? " mm" : ""}</text>`;
  }
  model.forEach((v, i) => {
    out += `<rect x="${X(i) - bw / 2}" y="${Y(v)}" width="${bw}" height="${Y(0) - Y(v)}" fill="#5b9bc9"/>`;
  });
  if (observed) {
    observed.forEach((v, i) => {
      out += `<line x1="${X(i) - bw * 0.8}" x2="${X(i) + bw * 0.8}" y1="${Y(v)}" y2="${Y(v)}" stroke="#e8e2d0" stroke-width="2"/>`;
    });
  }
  MONTHS.forEach((m, i) => { out += `<text x="${X(i)}" y="${height - 5}" text-anchor="middle">${m[0]}</text>`; });
  return out + "</svg>";
}

const zoneText = (code) => {
  const z = ZONE_NAMES.find((z) => z.code === code);
  return z ? `${z.code} · ${z.name}` : "—";
};

async function showCellChart(h) {
  const d = await api(`/api/climate/cell?cell=${h.cell}`);
  const months = (arr) => arr.map((v, i) => [i, v]);
  const series = [];
  const key = [];
  const useWindow = !!d.window;
  if (useWindow) {
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
  const rain = useWindow ? d.rain_window : d.rain_normal;
  const total = (a) => Math.round(a.reduce((x, y) => x + y, 0)).toLocaleString();
  const zone = useWindow ? d.zone_window : d.zone;
  const box = $("cell-chart");
  box.innerHTML = `<div class="head"><strong>Climate at ${coords(h)}</strong><button type="button" id="cell-close">×</button></div>` +
    `<div class="coords">Average temperature by month · height ${fmt(Math.max(0, h.elev))}</div>` +
    lineChart({ series, xTicks: MONTHS.map((m, i) => [i, m[0]]) }) + `<div class="key">${key.join("")}</div>` +
    `<div class="coords">Rain and snow by month: model ${total(rain)} mm a year` +
    (d.rain_observed ? `, observed ${total(d.rain_observed)} mm (white marks)` : "") + `</div>` +
    rainBars({ model: rain, observed: d.rain_observed }) +
    `<div class="coords">Climate zone: ${zoneText(zone)}` +
    (d.zone_observed ? ` · observed ${zoneText(d.zone_observed)}` : "") + `</div>`;
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
    `</table>` + rainCheckHtml(d) + oceanCheckHtml(d) +
    `<p class="note">Compare maps with each layer's "Observed" and "Model − observed" views.</p>`;
}

function rainCheckHtml(d) {
  const r = d.rain, pm = d.pass_marks;
  if (!r) return "";
  const row = (label, value, ok) => `<tr><td>${label}</td><td class="num">${value}</td>` +
    `<td class="num ${ok === undefined ? "" : ok ? "pass" : "fail"}">${ok === undefined ? "" : ok ? "pass" : "miss"}</td></tr>`;
  const pct = (x) => `${Math.round(100 * x)}%`;
  let out = `<h3>Rain and climate zones</h3><table>` +
    row(`Pattern of yearly rain matches observed (correlation, mark ${pm.rain_correlation})`, r.correlation.toFixed(2), r.pass_correlation) +
    row(`Land in the right climate zone (mark ${pct(pm.zone_agreement)})`, pct(r.zone_agreement), r.pass_zones) +
    row("Average yearly rain on land, model / observed", `${Math.round(r.land_mean_model_mm)} / ${Math.round(r.land_mean_observed_mm)} mm`) +
    `</table><table class="zones"><tr><th>Zone</th><th class="num">Model</th><th class="num">Observed</th></tr>` +
    d.zones.map((z, i) => `<tr><td>${zoneText(z)}</td><td class="num">${pct(r.zone_share_model[i])}</td>` +
      `<td class="num">${pct(r.zone_share_observed[i])}</td></tr>`).join("") +
    `</table><table>`;
  const what = { dry: "dry (under 250 mm)", wet: "wet (over 1,500 mm)", monsoon: "most rain June–September, 500 mm+" };
  for (const l of r.landmarks) {
    const v = l.expect === "monsoon" ? `${pct(l.model)} (obs ${pct(l.observed)})`
      : `${Math.round(l.model).toLocaleString()} mm (obs ${Math.round(l.observed).toLocaleString()})`;
    out += row(`${l.name}: ${what[l.expect]}`, v, l.pass);
  }
  return out + `</table>`;
}

function oceanCheckHtml(d) {
  const o = d.ocean;
  if (!o) return "";
  const row = (label, value, ok) => `<tr><td>${label}</td><td class="num">${value}</td>` +
    `<td class="num ${ok ? "pass" : "fail"}">${ok ? "pass" : "miss"}</td></tr>`;
  let out = `<h3>Oceans and sea ice</h3><table>`;
  for (const i of o.sea_ice) {
    out += row(`Sea ice, ${i.name} (observed ${i.observed_million_km2} million km², within 25%)`,
      `${i.model_million_km2.toFixed(1)} million km²`, i.pass);
  }
  for (const c of o.contrasts) {
    out += row(`${c.name} (observed by ${c.observed_c.toFixed(0)} °C)`, `by ${c.model_c.toFixed(1)} °C`, c.pass);
  }
  return out + `</table>`;
}

async function refreshCheck() {
  const box = $("check");
  if (!state || state.scenario !== "present_day_check") { box.hidden = true; return; }
  box.hidden = false;
  $("check-body").innerHTML = checkHtml(await api("/api/climate/check"));
}

// The climate layers and what each can show.
const FIELDS = {
  temperature: { views: ["current", "normal", "observed", "difference"] },
  rain: { views: ["current", "normal", "observed", "difference"] },
  snow: { views: ["current", "normal"] },
  seaice: { views: ["current"] },
  zone: { views: ["normal", "observed"] },
};
const VIEW_LABELS = { current: "This month", normal: "Recent average", observed: "Observed", difference: "Model − observed" };
const temp = { view: "normal", month: null };
let lastTempTick = -1;

function tempOptions(map) {
  const box = $("temp-options");
  const f = FIELDS[map.layer()];
  box.hidden = !f;
  if (box.hidden) return;
  const checkOnly = ["observed", "difference"];
  const views = f.views
    .filter((v) => state?.scenario === "present_day_check" || !checkOnly.includes(v))
    .map((id) => ({ id, label: map.layer() === "zone" && id === "normal" ? "Model" : VIEW_LABELS[id] }));
  if (!views.some((v) => v.id === temp.view)) temp.view = "normal";
  optionButtons(box, views, temp.view, (id) => { temp.view = id; tempOptions(map); refreshTemperature(map, true); });
  if (temp.view !== "current" && map.layer() !== "zone") {
    const sel = document.createElement("select");
    sel.innerHTML = `<option value="">Whole year</option>` + MONTH_NAMES.map((m, i) => `<option value="${i}">${m}</option>`).join("");
    sel.value = temp.month ?? "";
    sel.addEventListener("change", () => { temp.month = sel.value === "" ? null : Number(sel.value); refreshTemperature(map, true); });
    box.appendChild(sel);
  }
}

async function refreshTemperature(map, force) {
  const layer = map.layer();
  if (!FIELDS[layer] || !state) return;
  const fixed = temp.view === "observed";
  if (!force && (fixed || state.tick === lastTempTick)) return;
  lastTempTick = state.tick;
  const month = layer !== "zone" && temp.view !== "current" ? temp.month : null;
  const q = `/api/climate/map?field=${layer}&view=${temp.view}${month !== null ? `&month=${month}` : ""}`;
  const buf = await fetch(q).then((r) => r.arrayBuffer());
  if (map.layer() !== layer) return;
  const when = temp.view === "current" ? state.date_label.replace(/^\d+ /, "")
    : month === null ? "whole year" : MONTH_NAMES[month];
  const recent = state.scenario === "present_day_check" ? "1981–2010" : "recent years";
  let mode, title;
  if (layer === "temperature") {
    mode = temp.view === "difference" ? "diff" : "temp";
    title = {
      current: `Temperature, ${when} (°C)`,
      normal: `Average temperature, ${when}, recent years (°C)`,
      observed: `Observed average, ${when}, 1981–2010 (°C)`,
      difference: `Model minus observed, ${when} (°C)`,
    }[temp.view];
  } else if (layer === "rain") {
    mode = temp.view === "difference" ? "ratio" : "rain";
    title = {
      current: `Rain and snow, ${when} (mm)`,
      normal: `Rain and snow, ${when}, recent years (mm)`,
      observed: `Observed rain and snow, ${when}, 1981–2010 (mm)`,
      difference: `Model rain as a multiple of observed, ${when}`,
    }[temp.view];
  } else if (layer === "seaice") {
    mode = "seaice";
    title = `Sea ice thickness now, ${state.date_label.replace(/^\d+ /, "")} (cm)`;
  } else if (layer === "snow") {
    mode = "snow";
    title = temp.view === "current" ? `Lying snow now (mm of water)`
      : month === null ? "Deepest snow of the year, recent years (mm of water)" : `Lying snow at the end of ${when}, recent years (mm of water)`;
  } else {
    mode = "zone";
    title = temp.view === "observed" ? "Observed climate zones, 1981–2010" : `Climate zones, model, ${recent}`;
  }
  const monthly = layer === "rain" && (temp.view === "current" || month !== null);
  map.setField(layer, new Int16Array(buf), mode, title, monthly);
  $("legend").innerHTML = map.legendHtml();
}

// Arrows over the map: this month's winds or this year's ocean currents.
const winds = { on: null, tick: -1 };
async function refreshWinds(map, force) {
  if (!winds.on) { map.setArrows(null); return; }
  if (!state || (!force && state.tick === winds.tick)) return;
  winds.tick = state.tick;
  if (winds.on === "winds") {
    const d = await api("/api/climate/winds");
    if (winds.on === "winds") map.setArrows(d.winds, 2.2, "rgba(20, 18, 14, 0.85)");
  } else {
    const d = await api("/api/climate/currents");
    if (winds.on === "currents") map.setArrows(d.currents, 100, "rgba(10, 30, 70, 0.9)", true);
  }
}

WorldMap.create($("map"), (h) => {
  const box = $("hover");
  box.hidden = !h;
  if (h) box.innerHTML = hoverHtml(h);
}, (h) => showCellChart(h).catch(() => {})).then((map) => {
  $("map-loading").remove();
  ZONE_NAMES = map.zones;
  for (const [id, kind] of [["winds-toggle", "winds"], ["currents-toggle", "currents"]]) {
    $(id).addEventListener("click", () => {
      winds.on = winds.on === kind ? null : kind;
      $("winds-toggle").setAttribute("aria-pressed", String(winds.on === "winds"));
      $("currents-toggle").setAttribute("aria-pressed", String(winds.on === "currents"));
      refreshWinds(map, true).catch(() => {});
    });
  }
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
    refreshWinds(map, false).catch(() => {});
  }, 1000);
  setInterval(() => { refreshWorldChart().catch(() => {}); refreshCheck().catch(() => {}); }, 2000);
}).catch((err) => {
  $("map-loading").textContent = `The map could not be drawn: ${err.message}`;
});
