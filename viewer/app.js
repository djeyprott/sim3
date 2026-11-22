// The viewer polls the local simulator for its state and new events.

const $ = (id) => document.getElementById(id);
let state = null;
let nextEventId = 0;

const SPEED_NAMES = {
  1: "1 month / sec",
  12: "1 year / sec",
  120: "10 years / sec",
  1200: "100 years / sec",
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
  $("year").textContent = state.year_label;
  $("date").textContent = state.date_label;
  const years = Math.floor(state.tick / 12);
  $("elapsed").textContent = `${years.toLocaleString()} years simulated`;
  $("toggle").textContent = state.paused ? "Play" : "Pause";

  const select = $("speed");
  if (select.options.length === 0) {
    for (const s of state.speeds) {
      const opt = document.createElement("option");
      opt.value = s;
      opt.textContent = SPEED_NAMES[s] ?? `${s} months / sec`;
      select.appendChild(opt);
    }
  }
  if (document.activeElement !== select) select.value = state.months_per_second;
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
$("step").addEventListener("click", async () => {
  await api("/api/control", { action: "step", months: 1 });
  setTimeout(poll, 60);
});
$("step-year").addEventListener("click", async () => {
  await api("/api/control", { action: "step", months: 12 });
  setTimeout(poll, 60);
});
$("speed").addEventListener("change", async (ev) => {
  await api("/api/control", { action: "speed", months_per_second: Number(ev.target.value) });
  poll();
});

poll();
setInterval(poll, 250);
