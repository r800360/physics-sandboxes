import { SCENE_PRESETS, createSheet, sheetEnergy, stepSheet, vonMises } from "./core.mjs";

const view = document.getElementById("view");
const ctx = view.getContext("2d");
const video = document.querySelector("video");
video.addEventListener("error", () => video.remove());

const names = ["fall", "drop_spin", "pinned", "shake", "stretch", "twist"];
const labels = {
  fall: "Fall",
  drop_spin: "Spin",
  pinned: "Pin",
  shake: "Shake",
  stretch: "Stretch",
  twist: "Twist",
};

const state = { scene: "fall", sim: null, playing: true, clim: 0.05 };
const scenes = document.getElementById("scenes");
for (const name of names) {
  const button = document.createElement("button");
  button.textContent = labels[name];
  button.dataset.scene = name;
  button.onclick = () => start(name);
  scenes.appendChild(button);
}

function start(name) {
  state.scene = name;
  const preset = SCENE_PRESETS[name];
  state.sim = createSheet({
    nx: 15,
    ny: 7,
    kind: preset.kind,
    cy: preset.cy,
    young: preset.young,
    damping: preset.damping,
    dt: 0.01,
    newtonIters: 5,
    kc: 4e5,
  });
  state.clim = 0.04;
  state.playing = true;
  document.getElementById("play").textContent = "Pause";
  for (const button of scenes.children) {
    button.classList.toggle("on", button.dataset.scene === name);
  }
}

function resize() {
  const dpr = window.devicePixelRatio || 1;
  const w = view.clientWidth;
  const h = view.clientHeight;
  if (view.width !== Math.floor(w * dpr) || view.height !== Math.floor(h * dpr)) {
    view.width = Math.floor(w * dpr);
    view.height = Math.floor(h * dpr);
  }
  return { w, h, dpr };
}

const MAGMA = [
  [244, 236, 220],
  [236, 196, 150],
  [220, 120, 70],
  [180, 50, 45],
  [90, 20, 30],
];

function magma(t) {
  const x = Math.min(1, Math.max(0, t)) * (MAGMA.length - 1);
  const i = Math.min(MAGMA.length - 2, Math.floor(x));
  const f = x - i;
  const c = MAGMA[i].map((v, k) => Math.round(v + (MAGMA[i + 1][k] - v) * f));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}

function draw() {
  const sim = state.sim;
  if (!sim) return;
  const { w, h, dpr } = resize();
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  const strains = vonMises(sim.x, sim.elems);
  let maxS = 0;
  const sorted = Array.from(strains).sort((a, b) => a - b);
  const p98 = sorted[Math.floor(0.98 * (sorted.length - 1))] || 0;
  state.clim = Math.max(p98, state.clim * 0.992, 0.02);
  for (const s of strains) maxS = Math.max(maxS, s);
  let minX = Infinity;
  let maxX = -Infinity;
  let maxY = 0.2;
  for (const p of sim.rest) {
    minX = Math.min(minX, p[0]);
    maxX = Math.max(maxX, p[0]);
    maxY = Math.max(maxY, p[1]);
  }
  minX -= 0.9;
  maxX += 1.3;
  const minY = -0.25;
  maxY += 0.55;
  const project = (x, y) => ({
    x: ((x - minX) / (maxX - minX)) * (w - 28) + 14,
    y: h - 16 - ((y - minY) / (maxY - minY)) * (h - 32),
  });
  const floor0 = project(minX, 0);
  const floor1 = project(maxX, 0);
  ctx.fillStyle = "rgba(28,24,20,0.05)";
  ctx.fillRect(floor0.x, floor0.y, floor1.x - floor0.x, h - floor0.y);
  ctx.strokeStyle = "#1b1814";
  ctx.lineWidth = 1.4;
  ctx.beginPath();
  ctx.moveTo(floor0.x, floor0.y);
  ctx.lineTo(floor1.x, floor1.y);
  ctx.stroke();
  ctx.strokeStyle = "rgba(28,24,20,0.18)";
  ctx.lineWidth = 0.6;
  const seen = new Set();
  for (const tri of sim.triangles) {
    for (let e = 0; e < 3; e++) {
      const a = tri[e];
      const b = tri[(e + 1) % 3];
      const key = a < b ? `${a}-${b}` : `${b}-${a}`;
      if (seen.has(key)) continue;
      seen.add(key);
      const p = project(sim.rest[a][0], sim.rest[a][1]);
      const q = project(sim.rest[b][0], sim.rest[b][1]);
      ctx.beginPath();
      ctx.moveTo(p.x, p.y);
      ctx.lineTo(q.x, q.y);
      ctx.stroke();
    }
  }
  sim.triangles.forEach((tri, idx) => {
    ctx.beginPath();
    tri.forEach((id, k) => {
      const p = project(sim.x[2 * id], sim.x[2 * id + 1]);
      if (k === 0) ctx.moveTo(p.x, p.y);
      else ctx.lineTo(p.x, p.y);
    });
    ctx.closePath();
    ctx.fillStyle = magma(strains[idx] / state.clim);
    ctx.fill();
    ctx.strokeStyle = "rgba(20,16,12,0.35)";
    ctx.lineWidth = 0.4;
    ctx.stroke();
  });
  for (const id of sim.handles.ids) {
    const p = project(sim.x[2 * id], sim.x[2 * id + 1]);
    ctx.fillStyle = "#fffdf8";
    ctx.strokeStyle = "#1b1814";
    ctx.lineWidth = 1.2;
    ctx.fillRect(p.x - 4, p.y - 4, 8, 8);
    ctx.strokeRect(p.x - 4, p.y - 4, 8, 8);
  }
  const info = sheetEnergy(sim);
  document.getElementById("tt").textContent = sim.t.toFixed(2);
  document.getElementById("st").textContent = maxS.toFixed(3);
  document.getElementById("pen").textContent = info.penetration.toExponential(1);
  document.getElementById("res").textContent = sim.residual.toExponential(1);
}

let last = performance.now();
let acc = 0;
function loop(now) {
  const elapsed = Math.min(0.05, (now - last) / 1000);
  last = now;
  if (state.playing && state.sim) {
    acc += elapsed;
    let guard = 0;
    while (acc >= state.sim.dt && guard < 1) {
      stepSheet(state.sim);
      acc -= state.sim.dt;
      guard += 1;
    }
  }
  draw();
  requestAnimationFrame(loop);
}

document.getElementById("play").onclick = () => {
  state.playing = !state.playing;
  document.getElementById("play").textContent = state.playing ? "Pause" : "Play";
};
document.getElementById("reset").onclick = () => start(state.scene);
window.addEventListener("keydown", (event) => {
  if (event.code === "Space") {
    event.preventDefault();
    document.getElementById("play").click();
  }
});

start("fall");
requestAnimationFrame(loop);
