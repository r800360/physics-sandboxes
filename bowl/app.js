import { defaultParams, simulate, world } from "./core.mjs";

const view = document.getElementById("view");
const chart = document.getElementById("chart");
const vctx = view.getContext("2d");
const cctx = chart.getContext("2d");
const video = document.querySelector("video");
video.addEventListener("error", () => video.remove());

const state = { playing: true, frame: 0, rk: null, sy: null, p: null };

function readParams() {
  const p = defaultParams();
  p.alpha = Number(document.getElementById("alpha").value);
  p.k = 10 ** Number(document.getElementById("k").value);
  p.dt = Number(document.getElementById("dt").value);
  p.T = 8;
  document.getElementById("aVal").textContent = p.alpha.toFixed(2);
  document.getElementById("kVal").textContent = p.k.toExponential(1);
  document.getElementById("dtVal").textContent = `${p.dt.toFixed(3)} s`;
  document.getElementById("xVal").textContent = Number(document.getElementById("x0").value).toFixed(2);
  return p;
}

function drift(frames) {
  let m = 0;
  let pen = 0;
  for (const f of frames) {
    m = Math.max(m, Math.abs(f.e - frames[0].e));
    pen = Math.max(pen, f.pen);
  }
  return { m, pen };
}

function rebuild() {
  const p = readParams();
  const x0 = Number(document.getElementById("x0").value);
  const q0 = [x0, 1.8, 0.45];
  const v0 = [0, 0, 0];
  state.p = p;
  state.rk = simulate("rk4", p, q0, v0);
  state.sy = simulate("symplectic", p, q0, v0);
  state.frame = 0;
  const a = drift(state.rk);
  const b = drift(state.sy);
  document.getElementById("erk").textContent = a.m.toFixed(3);
  document.getElementById("esy").textContent = b.m.toFixed(2);
  document.getElementById("pen").textContent = a.pen.toFixed(3);
}

function resize(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  if (canvas.width !== Math.floor(w * dpr) || canvas.height !== Math.floor(h * dpr)) {
    canvas.width = Math.floor(w * dpr);
    canvas.height = Math.floor(h * dpr);
  }
  return { w, h, dpr };
}

function camera(w, h) {
  const x0 = -3.3;
  const x1 = 3.3;
  const y0 = -0.45;
  const y1 = 3.55;
  const s = Math.min((w - 28) / (x1 - x0), (h - 28) / (y1 - y0));
  const ox = (w - s * (x1 - x0)) / 2;
  const oy = (h - s * (y1 - y0)) / 2;
  return { x0, y0, y1, s, ox, oy, h };
}

function mapPoint(x, y, cam) {
  return {
    x: cam.ox + (x - cam.x0) * cam.s,
    y: cam.h - cam.oy - (y - cam.y0) * cam.s,
  };
}

function drawBody(q, color, width, cam) {
  const p = state.p;
  const ends = [[-p.l1, 0], [p.l2, 0]].map(([rx, ry]) => {
    const X = world(q, rx, ry);
    return mapPoint(X[0], X[1], cam);
  });
  vctx.strokeStyle = color;
  vctx.lineWidth = width;
  vctx.beginPath();
  vctx.moveTo(ends[0].x, ends[0].y);
  vctx.lineTo(ends[1].x, ends[1].y);
  vctx.stroke();
  const radii = [p.r1, p.r2];
  const centers = [[-p.l1, 0], [p.l2, 0]];
  centers.forEach(([rx, ry], k) => {
    const X = world(q, rx, ry);
    const o = mapPoint(X[0], X[1], cam);
    const edge = mapPoint(X[0] + radii[k], X[1], cam);
    vctx.beginPath();
    vctx.arc(o.x, o.y, Math.abs(edge.x - o.x), 0, Math.PI * 2);
    vctx.stroke();
  });
}

function draw() {
  if (!state.rk) return;
  const { w, h, dpr } = resize(view);
  vctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  vctx.clearRect(0, 0, w, h);
  vctx.strokeStyle = "#1b1814";
  vctx.lineWidth = 2;
  const cam = camera(w, h);
  vctx.beginPath();
  for (let i = 0; i <= 80; i++) {
    const x = -3.2 + (6.4 * i) / 80;
    const y = state.p.alpha * x * x;
    const pt = mapPoint(x, y, cam);
    if (i === 0) vctx.moveTo(pt.x, pt.y);
    else vctx.lineTo(pt.x, pt.y);
  }
  vctx.stroke();
  const i = Math.min(state.frame, state.rk.length - 1);
  vctx.strokeStyle = "rgba(184,92,20,0.35)";
  vctx.lineWidth = 1.2;
  vctx.beginPath();
  state.sy.slice(0, i + 1).forEach((f, k) => {
    const pt = mapPoint(f.q[0], f.q[1], cam);
    if (k === 0) vctx.moveTo(pt.x, pt.y);
    else vctx.lineTo(pt.x, pt.y);
  });
  vctx.stroke();
  vctx.strokeStyle = "rgba(15,110,104,0.8)";
  vctx.beginPath();
  state.rk.slice(0, i + 1).forEach((f, k) => {
    const pt = mapPoint(f.q[0], f.q[1], cam);
    if (k === 0) vctx.moveTo(pt.x, pt.y);
    else vctx.lineTo(pt.x, pt.y);
  });
  vctx.stroke();
  drawBody(state.sy[i].q, "rgba(184,92,20,0.55)", 1.5, cam);
  drawBody(state.rk[i].q, "#0f6e68", 2.4, cam);
  document.getElementById("tt").textContent = `${state.rk[i].t.toFixed(2)} s`;
  drawChart(i);
}

function drawChart(i) {
  const { w, h, dpr } = resize(chart);
  cctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  cctx.clearRect(0, 0, w, h);
  const series = [state.rk, state.sy];
  const colors = ["#0f6e68", "#b85c14"];
  let max = 1e-6;
  for (const run of series) {
    for (const f of run) max = Math.max(max, Math.abs(f.e - run[0].e));
  }
  series.forEach((run, s) => {
    cctx.strokeStyle = colors[s];
    cctx.lineWidth = 1.5;
    cctx.beginPath();
    run.forEach((f, k) => {
      const x = (k / (run.length - 1)) * (w - 16) + 8;
      const y = h / 2 - ((f.e - run[0].e) / max) * (h * 0.38);
      if (k === 0) cctx.moveTo(x, y);
      else cctx.lineTo(x, y);
    });
    cctx.stroke();
  });
  cctx.fillStyle = "#665e54";
  cctx.font = "12px Segoe UI";
  cctx.fillText("energy − E₀    teal RK4    amber symplectic", 10, 16);
}

let acc = 0;
let last = performance.now();
function loop(now) {
  const elapsed = Math.min(0.05, (now - last) / 1000);
  last = now;
  if (state.playing && state.rk) {
    acc += elapsed;
    const hop = Math.min(8, Math.floor(acc / (state.p.dt * 5)));
    if (hop > 0) {
      state.frame = Math.min(state.rk.length - 1, state.frame + hop);
      acc -= hop * state.p.dt * 5;
      if (state.frame >= state.rk.length - 1) state.playing = false;
    }
  }
  draw();
  requestAnimationFrame(loop);
}

document.getElementById("play").onclick = () => {
  state.playing = !state.playing;
  if (state.rk && state.frame >= state.rk.length - 1) {
    state.frame = 0;
    state.playing = true;
  }
  document.getElementById("play").textContent = state.playing ? "Pause" : "Play";
};
document.getElementById("reset").onclick = () => {
  document.getElementById("alpha").value = "0.35";
  document.getElementById("k").value = "4";
  document.getElementById("dt").value = "0.002";
  document.getElementById("x0").value = "0.8";
  rebuild();
  state.playing = true;
};
let pending = 0;
for (const id of ["alpha", "k", "dt", "x0"]) {
  document.getElementById(id).oninput = () => {
    clearTimeout(pending);
    pending = setTimeout(rebuild, 80);
  };
}
window.addEventListener("keydown", (event) => {
  if (event.code === "Space") {
    event.preventDefault();
    document.getElementById("play").click();
  }
});

rebuild();
requestAnimationFrame(loop);
