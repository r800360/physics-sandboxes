import { energy, rk4Step, simpleStep } from "./core.mjs";

const view = document.getElementById("view");
const chart = document.getElementById("chart");
const vctx = view.getContext("2d");
const cctx = chart.getContext("2d");
const g = 9.81;
const theta0Default = 1.05;

const state = {
  L: 1.2,
  dt: 0.005,
  theta0: theta0Default,
  playing: true,
  showSimple: true,
  frame: 0,
  rod: null,
  simple: null,
};

const lengthInput = document.getElementById("length");
const dtInput = document.getElementById("dt");
const playBtn = document.getElementById("play");
const video = document.querySelector("video");
video.addEventListener("error", () => video.remove());

function integrate(stepFn) {
  const n = Math.max(2, Math.round(12 / state.dt));
  const theta = new Float64Array(n);
  const omega = new Float64Array(n);
  const eng = new Float64Array(n);
  let y = [state.theta0, 0];
  for (let i = 0; i < n; i++) {
    theta[i] = y[0];
    omega[i] = y[1];
    eng[i] = energy(y[0], y[1], 1, state.L, g);
    y = stepFn(y, state.dt, g, state.L);
  }
  let drift = 0;
  for (let i = 0; i < n; i++) drift = Math.max(drift, Math.abs(eng[i] - eng[0]));
  return { theta, omega, energy: eng, drift, n };
}

function rebuild() {
  state.L = Number(lengthInput.value);
  state.dt = Number(dtInput.value);
  document.getElementById("lVal").textContent = `${state.L.toFixed(2)} m`;
  document.getElementById("dtVal").textContent = `${state.dt.toFixed(3)} s`;
  state.rod = integrate(rk4Step);
  state.simple = integrate(simpleStep);
  state.frame = 0;
  document.getElementById("e0").textContent = state.rod.energy[0].toFixed(3);
  document.getElementById("de").textContent = state.rod.drift.toExponential(2);
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

function draw() {
  const { w, h, dpr } = resize(view);
  vctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  vctx.clearRect(0, 0, w, h);
  const pivot = { x: w * 0.5, y: h * 0.22 };
  const scale = Math.min(w * 0.7, h * 0.62) / state.L;
  const i = Math.min(state.frame, state.rod.n - 1);
  const drawRod = (theta, color, width) => {
    const x = state.L * Math.sin(theta);
    const z = -state.L * Math.cos(theta);
    const tip = { x: pivot.x + x * scale, y: pivot.y - z * scale };
    vctx.strokeStyle = color;
    vctx.lineWidth = width;
    vctx.lineCap = "round";
    vctx.beginPath();
    vctx.moveTo(pivot.x, pivot.y);
    vctx.lineTo(tip.x, tip.y);
    vctx.stroke();
    return tip;
  };
  vctx.setLineDash([4, 6]);
  vctx.strokeStyle = "#c8bfb2";
  vctx.lineWidth = 1;
  vctx.beginPath();
  vctx.moveTo(pivot.x, pivot.y);
  vctx.lineTo(pivot.x, pivot.y + state.L * scale);
  vctx.stroke();
  vctx.setLineDash([]);
  if (state.showSimple) drawRod(state.simple.theta[i], "#b85c14", 3);
  const tip = drawRod(state.rod.theta[i], "#0f6e68", 5);
  const midTheta = state.rod.theta[i];
  const mx = (state.L / 2) * Math.sin(midTheta);
  const mz = -(state.L / 2) * Math.cos(midTheta);
  vctx.fillStyle = "#1b1814";
  vctx.beginPath();
  vctx.arc(pivot.x + mx * scale, pivot.y - mz * scale, 4, 0, Math.PI * 2);
  vctx.fill();
  vctx.beginPath();
  vctx.arc(tip.x, tip.y, 7, 0, Math.PI * 2);
  vctx.fill();
  vctx.fillStyle = "#b85c14";
  vctx.beginPath();
  vctx.arc(pivot.x, pivot.y, 6, 0, Math.PI * 2);
  vctx.fill();
  document.getElementById("th").textContent = state.rod.theta[i].toFixed(3);
  document.getElementById("om").textContent = state.rod.omega[i].toFixed(3);
  drawChart(i);
}

function drawChart(i) {
  const { w, h, dpr } = resize(chart);
  cctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  cctx.clearRect(0, 0, w, h);
  const series = state.rod.energy;
  const e0 = series[0];
  let max = 1e-8;
  for (let k = 0; k < series.length; k++) max = Math.max(max, Math.abs(series[k] - e0));
  max *= 1.3;
  cctx.strokeStyle = "#0f6e68";
  cctx.lineWidth = 1.6;
  cctx.beginPath();
  for (let k = 0; k < series.length; k++) {
    const x = (k / (series.length - 1)) * (w - 16) + 8;
    const y = h / 2 - ((series[k] - e0) / max) * (h * 0.38);
    if (k === 0) cctx.moveTo(x, y);
    else cctx.lineTo(x, y);
  }
  cctx.stroke();
  const cx = (i / (series.length - 1)) * (w - 16) + 8;
  cctx.strokeStyle = "#b85c14";
  cctx.beginPath();
  cctx.moveTo(cx, 8);
  cctx.lineTo(cx, h - 8);
  cctx.stroke();
  cctx.fillStyle = "#665e54";
  cctx.font = "12px Segoe UI";
  cctx.fillText("energy − E₀", 10, 16);
}

let acc = 0;
let last = performance.now();
function loop(now) {
  const elapsed = Math.min(0.05, (now - last) / 1000);
  last = now;
  if (state.playing && state.rod) {
    acc += elapsed;
    const hop = Math.min(40, Math.floor(acc / state.dt));
    if (hop > 0) {
      state.frame = (state.frame + hop) % state.rod.n;
      acc -= hop * state.dt;
    }
  }
  draw();
  requestAnimationFrame(loop);
}

view.addEventListener("pointerdown", (event) => {
  const rect = view.getBoundingClientRect();
  const pivot = { x: rect.width * 0.5, y: rect.height * 0.22 };
  const dx = event.clientX - rect.left - pivot.x;
  const below = event.clientY - rect.top - pivot.y;
  state.theta0 = Math.atan2(dx, below);
  state.playing = false;
  playBtn.textContent = "Play";
  rebuild();
  draw();
});

document.getElementById("reset").onclick = () => {
  state.theta0 = theta0Default;
  lengthInput.value = "1.2";
  dtInput.value = "0.005";
  rebuild();
};
playBtn.onclick = () => {
  state.playing = !state.playing;
  playBtn.textContent = state.playing ? "Pause" : "Play";
};
lengthInput.oninput = rebuild;
dtInput.oninput = rebuild;
document.getElementById("simple").onchange = (event) => {
  state.showSimple = event.target.checked;
};
window.addEventListener("keydown", (event) => {
  if (event.code === "Space") {
    event.preventDefault();
    playBtn.click();
  }
});
window.addEventListener("resize", draw);

rebuild();
requestAnimationFrame(loop);
