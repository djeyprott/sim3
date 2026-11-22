// World map: draws the 655,362 cells with WebGL2.
//
// The server sends a lookup image saying which cell each map pixel belongs
// to (equirectangular, 4096 x 2048) and per-cell arrays for each moment.
// This file colours every cell for the chosen layer, and the graphics card
// paints each screen pixel with the colour of its cell.

const WorldMap = (() => {
  const COLOR_TEX_W = 1024;

  // ---- colour scales (display only) ----
  const OCEAN_STOPS = [ // depth below sea level (m) -> colour
    [0, [96, 150, 186]], [200, [70, 120, 165]], [2000, [44, 86, 132]], [5000, [27, 55, 96]], [8000, [16, 34, 66]],
  ];
  const LAND_STOPS = [ // height above sea level (m) -> colour
    [0, [94, 128, 76]], [300, [130, 150, 90]], [800, [178, 166, 112]], [1500, [160, 128, 92]],
    [2500, [134, 110, 92]], [4000, [168, 160, 152]], [6000, [230, 228, 224]],
  ];
  const ICE_STOPS = [ // ice thickness (m) -> colour
    [0, [226, 236, 242]], [1000, [188, 214, 232]], [2500, [132, 176, 214]], [4000, [80, 132, 190]],
  ];
  const BASIN = [104, 134, 128];
  const FLAT = { ocean: [43, 74, 107], land: [139, 139, 90], ice: [232, 238, 243], basin: BASIN };

  function ramp(stops, v) {
    if (v <= stops[0][0]) return stops[0][1];
    for (let i = 1; i < stops.length; i++) {
      if (v <= stops[i][0]) {
        const [v0, c0] = stops[i - 1], [v1, c1] = stops[i];
        const t = (v - v0) / (v1 - v0);
        return [0, 1, 2].map((k) => c0[k] + t * (c1[k] - c0[k]));
      }
    }
    return stops[stops.length - 1][1];
  }

  // Shaded relief: light from the north-west, 45° up. Slopes are exaggerated
  // 20x so 28 km cells still show relief (display only).
  const LIGHT = (() => { const l = [-0.5, 0.5, Math.SQRT1_2]; return l; })();
  function shade(se, sn) {
    const k = 20e-4; // slopes arrive as m/m x 1e4
    const nx = -se * k, ny = -sn * k, nz = 1;
    const len = Math.hypot(nx, ny, nz);
    const d = (nx * LIGHT[0] + ny * LIGHT[1] + nz * LIGHT[2]) / len / LIGHT[2];
    return Math.min(1.35, Math.max(0.55, d));
  }

  const LAYERS = {
    terrain: {
      name: "Terrain",
      color(c, i) {
        const kind = c.kind[i];
        if (c.ice[i] > 0 && kind !== 0) return ramp(ICE_STOPS, c.ice[i]).map((x) => x * shade(c.se[i], c.sn[i]));
        if (kind === 0) return ramp(OCEAN_STOPS, -c.elev[i]);
        const base = kind === 2 ? BASIN : ramp(LAND_STOPS, c.elev[i]);
        return base.map((x) => x * shade(c.se[i], c.sn[i]));
      },
      legend: () => rampLegend("Height (m)", LAND_STOPS, ["0", "1,500", "4,000", "6,000"]) +
        rampLegend("Sea depth (m)", OCEAN_STOPS, ["0", "2,000", "5,000", "8,000"]) +
        swatch(ICE_STOPS[0][1], "Ice") + swatch(BASIN, "Below sea level, cut off from the ocean"),
    },
    surface: {
      name: "Land, sea, ice",
      color(c, i) {
        if (c.ice[i] > 0 && c.kind[i] !== 0) return FLAT.ice;
        return [FLAT.ocean, FLAT.land, FLAT.basin][c.kind[i]];
      },
      legend: () => swatch(FLAT.land, "Land") + swatch(FLAT.ocean, "Sea") + swatch(FLAT.ice, "Ice") +
        swatch(FLAT.basin, "Below sea level, cut off from the ocean"),
    },
    ice: {
      name: "Ice thickness",
      color(c, i) {
        if (c.ice[i] > 0) return ramp(ICE_STOPS, c.ice[i]);
        return c.kind[i] === 0 ? [30, 40, 52] : [70, 68, 62];
      },
      legend: () => rampLegend("Ice thickness (m)", ICE_STOPS, ["0", "1,000", "2,500", "4,000"]),
    },
  };

  function css(c) { return `rgb(${c.map(Math.round).join(",")})`; }
  function swatch(c, label) {
    return `<div class="row"><span class="swatch" style="background:${css(c)}"></span>${label}</div>`;
  }
  function rampLegend(title, stops, ticks) {
    const max = stops[stops.length - 1][0];
    const grad = stops.map(([v, c]) => `${css(c)} ${(v / max) * 100}%`).join(",");
    return `<div>${title}</div><div class="ramp" style="background:linear-gradient(to right,${grad})"></div>` +
      `<div class="ticks">${ticks.map((t) => `<span>${t}</span>`).join("")}</div>`;
  }

  // ---- WebGL ----
  const VS = `#version 300 es
    in vec2 pos; void main() { gl_Position = vec4(pos, 0.0, 1.0); }`;
  const FS = `#version 300 es
    precision highp float; precision highp usampler2D;
    uniform usampler2D cellIndex; uniform sampler2D cellColor;
    uniform vec2 resolution; uniform vec2 center; uniform float scale;
    uniform ivec2 rasterSize; uniform int colorWidth;
    out vec4 outColor;
    void main() {
      vec2 p = gl_FragCoord.xy;
      float lon = center.x + (p.x - resolution.x * 0.5) / scale;
      float lat = center.y + (p.y - resolution.y * 0.5) / scale;
      if (lat > 90.0 || lat < -90.0) { outColor = vec4(0.05, 0.047, 0.035, 1.0); return; }
      lon = mod(lon + 180.0, 360.0);
      int x = min(int(lon / 360.0 * float(rasterSize.x)), rasterSize.x - 1);
      int y = min(int((90.0 - lat) / 180.0 * float(rasterSize.y)), rasterSize.y - 1);
      uint idx = texelFetch(cellIndex, ivec2(x, y), 0).r;
      int i = int(idx);
      outColor = texelFetch(cellColor, ivec2(i % colorWidth, i / colorWidth), 0);
    }`;

  function compile(gl, type, src) {
    const s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
    return s;
  }

  async function create(canvas, onHover) {
    const gl = canvas.getContext("webgl2", { antialias: false });
    if (!gl) throw new Error("This browser does not support WebGL2.");

    const [info, rasterBuf] = await Promise.all([
      fetch("/api/map/info").then((r) => r.json()),
      fetch("/api/map/raster").then((r) => r.arrayBuffer()),
    ]);
    const n = info.cells;
    const head = new Uint32Array(rasterBuf, 0, 2);
    const rw = head[0], rh = head[1];
    const raster = new Uint32Array(rasterBuf, 8, rw * rh);

    const epochs = await Promise.all(info.epochs.map(async (e) => {
      const b = await fetch(`/api/map/cells?epoch=${e.index}`).then((r) => r.arrayBuffer());
      return {
        ...e,
        elev: new Int16Array(b, 0, n), elevMin: new Int16Array(b, 2 * n, n), elevMax: new Int16Array(b, 4 * n, n),
        ice: new Uint16Array(b, 6 * n, n), se: new Int16Array(b, 8 * n, n), sn: new Int16Array(b, 10 * n, n),
        land: new Uint8Array(b, 12 * n, n), kind: new Uint8Array(b, 13 * n, n),
      };
    }));

    const prog = gl.createProgram();
    gl.attachShader(prog, compile(gl, gl.VERTEX_SHADER, VS));
    gl.attachShader(prog, compile(gl, gl.FRAGMENT_SHADER, FS));
    gl.linkProgram(prog);
    gl.useProgram(prog);
    const buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), gl.STATIC_DRAW);
    const loc = gl.getAttribLocation(prog, "pos");
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);

    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    const idxTex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, idxTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32UI, rw, rh, 0, gl.RED_INTEGER, gl.UNSIGNED_INT, raster);
    for (const p of [gl.TEXTURE_MIN_FILTER, gl.TEXTURE_MAG_FILTER]) gl.texParameteri(gl.TEXTURE_2D, p, gl.NEAREST);

    const colorH = Math.ceil(n / COLOR_TEX_W);
    const colors = new Uint8ClampedArray(COLOR_TEX_W * colorH * 4);
    const colTex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, colTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, COLOR_TEX_W, colorH, 0, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array(colors.buffer));
    for (const p of [gl.TEXTURE_MIN_FILTER, gl.TEXTURE_MAG_FILTER]) gl.texParameteri(gl.TEXTURE_2D, p, gl.NEAREST);

    const u = (name) => gl.getUniformLocation(prog, name);
    gl.uniform1i(u("cellIndex"), 0);
    gl.uniform1i(u("cellColor"), 1);
    gl.uniform2i(u("rasterSize"), rw, rh);
    gl.uniform1i(u("colorWidth"), COLOR_TEX_W);

    const state = { epoch: 0, layer: "terrain", center: [10, 20], scale: 0 };

    function recolor() {
      const c = epochs[state.epoch], layer = LAYERS[state.layer];
      for (let i = 0; i < n; i++) {
        const col = layer.color(c, i);
        colors[i * 4] = col[0]; colors[i * 4 + 1] = col[1]; colors[i * 4 + 2] = col[2]; colors[i * 4 + 3] = 255;
      }
      gl.activeTexture(gl.TEXTURE1);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, COLOR_TEX_W, colorH, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array(colors.buffer));
      draw();
    }

    function minScale() { return canvas.height / 180; }
    function resize() {
      const dpr = window.devicePixelRatio || 1;
      canvas.width = Math.round(canvas.clientWidth * dpr);
      canvas.height = Math.round(canvas.clientHeight * dpr);
      if (!state.scale) state.scale = Math.max(minScale(), canvas.width / 360);
      state.scale = Math.max(state.scale, minScale());
      draw();
    }
    function clampCenter() {
      const halfLat = canvas.height / 2 / state.scale;
      state.center[1] = Math.max(-90 + halfLat, Math.min(90 - halfLat, state.center[1]));
      if (halfLat >= 90) state.center[1] = 0;
    }
    function draw() {
      clampCenter();
      gl.viewport(0, 0, canvas.width, canvas.height);
      gl.uniform2f(u("resolution"), canvas.width, canvas.height);
      gl.uniform2f(u("center"), state.center[0], state.center[1]);
      gl.uniform1f(u("scale"), state.scale);
      gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    }

    function lonLatAt(clientX, clientY) {
      const r = canvas.getBoundingClientRect(), dpr = canvas.width / r.width;
      const px = (clientX - r.left) * dpr, py = (r.bottom - clientY) * dpr;
      const lon = state.center[0] + (px - canvas.width / 2) / state.scale;
      const lat = state.center[1] + (py - canvas.height / 2) / state.scale;
      return [((lon + 180) % 360 + 360) % 360 - 180, lat];
    }
    function cellAt(lon, lat) {
      if (lat > 90 || lat < -90) return -1;
      const x = Math.min(rw - 1, Math.floor((lon + 180) / 360 * rw));
      const y = Math.min(rh - 1, Math.floor((90 - lat) / 180 * rh));
      return raster[y * rw + x];
    }

    let drag = null;
    canvas.addEventListener("mousedown", (e) => { drag = [e.clientX, e.clientY]; canvas.classList.add("dragging"); });
    window.addEventListener("mouseup", () => { drag = null; canvas.classList.remove("dragging"); });
    window.addEventListener("mousemove", (e) => {
      if (drag) {
        const dpr = canvas.width / canvas.getBoundingClientRect().width;
        state.center[0] -= (e.clientX - drag[0]) * dpr / state.scale;
        state.center[1] += (e.clientY - drag[1]) * dpr / state.scale;
        drag = [e.clientX, e.clientY];
        draw();
      }
      if (e.target === canvas) {
        const [lon, lat] = lonLatAt(e.clientX, e.clientY);
        const i = cellAt(lon, lat);
        onHover(i < 0 ? null : describe(i, lon, lat));
      }
    });
    canvas.addEventListener("mouseleave", () => onHover(null));
    canvas.addEventListener("wheel", (e) => {
      e.preventDefault();
      const [lon0, lat0] = lonLatAt(e.clientX, e.clientY);
      const factor = Math.exp(-e.deltaY * 0.0015);
      state.scale = Math.min(400, Math.max(minScale(), state.scale * factor));
      const [lon1, lat1] = lonLatAt(e.clientX, e.clientY);
      let dlon = lon0 - lon1;
      if (dlon > 180) dlon -= 360; else if (dlon < -180) dlon += 360;
      state.center[0] += dlon;
      state.center[1] += lat0 - lat1;
      draw();
    }, { passive: false });
    window.addEventListener("resize", resize);

    function describe(i, lon, lat) {
      const c = epochs[state.epoch];
      return {
        lat, lon, cell: i, kind: c.kind[i], elev: c.elev[i], elevMin: c.elevMin[i], elevMax: c.elevMax[i],
        landShare: c.land[i] / 255, ice: c.ice[i], epochLabel: c.label,
      };
    }

    // A view can be opened directly with a link ending in #<lat>,<lon>,<zoom>,
    // e.g. #54,4,8 for the North Sea.
    function viewFromHash() {
      const m = location.hash.slice(1).split(",").map(Number);
      if (m.length >= 2 && m.every(Number.isFinite)) {
        state.center = [m[1], m[0]];
        if (m[2] > 0) state.scale = minScale() * m[2];
      }
    }

    resize();
    viewFromHash();
    window.addEventListener("hashchange", () => { viewFromHash(); draw(); });
    recolor();

    return {
      epochs: info.epochs,
      layers: Object.entries(LAYERS).map(([id, l]) => ({ id, name: l.name })),
      setEpoch(i) { state.epoch = i; recolor(); },
      setLayer(id) { state.layer = id; recolor(); },
      legendHtml() { return LAYERS[state.layer].legend(); },
    };
  }

  return { create };
})();
