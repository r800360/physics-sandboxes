/**
 * Live corotated sheet. Same energy, same incremental potential, same
 * rest-pose Hessian as elastic_fem/. The mesh is coarser so it can run
 * in the browser. The linear solve is dense; the Python run uses CG
 * on the same matrix.
 */

const GRAVITY = -9.81;

function smoothstep(t) {
  const s = Math.min(1, Math.max(0, t));
  return s * s * (3 - 2 * s);
}

function polar(F00, F01, F10, F11) {
  const p = F00 * F00 + F10 * F10;
  const q = F00 * F01 + F10 * F11;
  const r = F01 * F01 + F11 * F11;
  const tr = p + r;
  const det = p * r - q * q;
  const disc = Math.sqrt(Math.max(0, (tr * tr) / 4 - det));
  const s1 = Math.sqrt(Math.max(0, tr / 2 + disc));
  const s2 = Math.sqrt(Math.max(0, tr / 2 - disc));
  const l1 = tr / 2 + disc;
  let v0;
  let v1;
  if (Math.abs(q) > 1e-15) {
    v0 = q;
    v1 = l1 - p;
    const nrm = Math.hypot(v0, v1) || 1;
    v0 /= nrm;
    v1 /= nrm;
  } else if (p >= r) {
    v0 = 1;
    v1 = 0;
  } else {
    v0 = 0;
    v1 = 1;
  }
  const w0 = -v1;
  const w1 = v0;
  const column = (vx, vy, sig) => {
    const ux = F00 * vx + F01 * vy;
    const uy = F10 * vx + F11 * vy;
    if (sig < 1e-14) return [0, 0];
    return [ux / sig, uy / sig];
  };
  let [u0x, u0y] = column(v0, v1, s1);
  let [u1x, u1y] = column(w0, w1, s2);
  if (s2 < 1e-12) {
    const n = Math.hypot(u0x, u0y) || 1;
    u0x /= n;
    u0y /= n;
    u1x = -u0y;
    u1y = u0x;
  }
  let r00 = u0x * v0 + u1x * w0;
  let r01 = u0x * v1 + u1x * w1;
  let r10 = u0y * v0 + u1y * w0;
  let r11 = u0y * v1 + u1y * w1;
  if (r00 * r11 - r01 * r10 < 0) {
    u1x = -u1x;
    u1y = -u1y;
    r00 = u0x * v0 + u1x * w0;
    r01 = u0x * v1 + u1x * w1;
    r10 = u0y * v0 + u1y * w0;
    r11 = u0y * v1 + u1y * w1;
  }
  return [r00, r01, r10, r11];
}

function meshOf(nx, ny, width, height, cx, cy) {
  const V = [];
  for (let j = 0; j < ny; j++) {
    const y = -height / 2 + (height * j) / (ny - 1) + cy;
    for (let i = 0; i < nx; i++) {
      const x = -width / 2 + (width * i) / (nx - 1) + cx;
      V.push([x, y]);
    }
  }
  const T = [];
  const vid = (i, j) => j * nx + i;
  for (let j = 0; j < ny - 1; j++) {
    for (let i = 0; i < nx - 1; i++) {
      const a = vid(i, j);
      const b = vid(i + 1, j);
      const c = vid(i, j + 1);
      const d = vid(i + 1, j + 1);
      if ((i + j) % 2 === 0) {
        T.push([a, b, d]);
        T.push([a, d, c]);
      } else {
        T.push([a, b, c]);
        T.push([b, d, c]);
      }
    }
  }
  return { V, T, nx, ny };
}

function buildElements(V, T) {
  return T.map((tri) => {
    const [x0, y0] = V[tri[0]];
    const [x1, y1] = V[tri[1]];
    const [x2, y2] = V[tri[2]];
    const a = x1 - x0;
    const b = x2 - x0;
    const c = y1 - y0;
    const d = y2 - y0;
    const det = a * d - b * c;
    const inv = 1 / det;
    return {
      tri,
      area: 0.5 * det,
      DmInv: [d * inv, -b * inv, -c * inv, a * inv],
    };
  });
}

function assembleK0(V, elems, mu, lam) {
  const n = V.length;
  const N = 2 * n;
  const K = new Float64Array(N * N);
  const D = [lam + 2 * mu, lam, 0, lam, lam + 2 * mu, 0, 0, 0, mu];
  for (const e of elems) {
    const [i0, i1, i2] = e.tri;
    const [x1, y1] = V[i0];
    const [x2, y2] = V[i1];
    const [x3, y3] = V[i2];
    const signed2A = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    const area = 0.5 * signed2A;
    const b1 = y2 - y3;
    const b2 = y3 - y1;
    const b3 = y1 - y2;
    const c1 = x3 - x2;
    const c2 = x1 - x3;
    const c3 = x2 - x1;
    const s = 1 / (2 * area);
    const B = [
      s * b1, 0, s * b2, 0, s * b3, 0,
      0, s * c1, 0, s * c2, 0, s * c3,
      s * c1, s * b1, s * c2, s * b2, s * c3, s * b3,
    ];
    const DB = new Float64Array(18);
    for (let row = 0; row < 3; row++) {
      for (let col = 0; col < 6; col++) {
        let sum = 0;
        for (let k = 0; k < 3; k++) sum += D[row * 3 + k] * B[k * 6 + col];
        DB[row * 6 + col] = sum;
      }
    }
    const dofs = [2 * i0, 2 * i0 + 1, 2 * i1, 2 * i1 + 1, 2 * i2, 2 * i2 + 1];
    for (let aa = 0; aa < 6; aa++) {
      for (let bb = 0; bb < 6; bb++) {
        let sum = 0;
        for (let k = 0; k < 3; k++) sum += B[k * 6 + aa] * DB[k * 6 + bb];
        K[dofs[aa] * N + dofs[bb]] += area * sum;
      }
    }
  }
  return K;
}

function deformation(x, e) {
  const [i0, i1, i2] = e.tri;
  const Ds00 = x[2 * i1] - x[2 * i0];
  const Ds01 = x[2 * i2] - x[2 * i0];
  const Ds10 = x[2 * i1 + 1] - x[2 * i0 + 1];
  const Ds11 = x[2 * i2 + 1] - x[2 * i0 + 1];
  const [m00, m01, m10, m11] = e.DmInv;
  return [
    Ds00 * m00 + Ds01 * m10,
    Ds00 * m01 + Ds01 * m11,
    Ds10 * m00 + Ds11 * m10,
    Ds10 * m01 + Ds11 * m11,
  ];
}

export function vonMises(x, elems) {
  const out = new Float64Array(elems.length);
  elems.forEach((e, idx) => {
    const [F00, F01, F10, F11] = deformation(x, e);
    const c00 = F00 * F00 + F10 * F10;
    const c01 = F00 * F01 + F10 * F11;
    const c11 = F01 * F01 + F11 * F11;
    const exx = 0.5 * (c00 - 1);
    const eyy = 0.5 * (c11 - 1);
    const exy = 0.5 * c01;
    out[idx] = Math.sqrt(Math.max(0, exx * exx - exx * eyy + eyy * eyy + 3 * exy * exy));
  });
  return out;
}

function elastic(x, elems, mu, lam) {
  const grad = new Float64Array(x.length);
  let energy = 0;
  for (const e of elems) {
    const [F00, F01, F10, F11] = deformation(x, e);
    const [R00, R01, R10, R11] = polar(F00, F01, F10, F11);
    const J = F00 * F11 - F01 * F10;
    const Js = Math.max(J, 1e-6);
    energy += e.area * (mu * ((F00 - R00) ** 2 + (F01 - R01) ** 2 + (F10 - R10) ** 2 + (F11 - R11) ** 2) + 0.5 * lam * (Js - 1) ** 2);
    let I00 = 1;
    let I01 = 0;
    let I10 = 0;
    let I11 = 1;
    if (Math.abs(J) > 1e-8) {
      const invJ = 1 / J;
      I00 = F11 * invJ;
      I01 = -F10 * invJ;
      I10 = -F01 * invJ;
      I11 = F00 * invJ;
    }
    const P00 = 2 * mu * (F00 - R00) + lam * (Js - 1) * Js * I00;
    const P01 = 2 * mu * (F01 - R01) + lam * (Js - 1) * Js * I01;
    const P10 = 2 * mu * (F10 - R10) + lam * (Js - 1) * Js * I10;
    const P11 = 2 * mu * (F11 - R11) + lam * (Js - 1) * Js * I11;
    const [m00, m01, m10, m11] = e.DmInv;
    const G00 = e.area * (P00 * m00 + P01 * m01);
    const G01 = e.area * (P00 * m10 + P01 * m11);
    const G10 = e.area * (P10 * m00 + P11 * m01);
    const G11 = e.area * (P10 * m10 + P11 * m11);
    const [i0, i1, i2] = e.tri;
    grad[2 * i1] += G00;
    grad[2 * i1 + 1] += G10;
    grad[2 * i2] += G01;
    grad[2 * i2 + 1] += G11;
    grad[2 * i0] += -G00 - G01;
    grad[2 * i0 + 1] += -G10 - G11;
  }
  return { energy, grad };
}

function contactTerms(x, k, floor) {
  const n = x.length / 2;
  const grad = new Float64Array(x.length);
  const diag = new Float64Array(x.length);
  let energy = 0;
  let maxPen = 0;
  for (let i = 0; i < n; i++) {
    const pen = floor - x[2 * i + 1];
    if (pen > maxPen) maxPen = pen;
    if (pen <= 0) continue;
    energy += 0.5 * k * pen * pen;
    grad[2 * i + 1] += k * (x[2 * i + 1] - floor);
    diag[2 * i + 1] = k;
  }
  return { energy, grad, diag, maxPen: Math.max(0, maxPen) };
}

function solveDense(A, b) {
  const n = b.length;
  const M = Float64Array.from(A);
  const x = Float64Array.from(b);
  for (let k = 0; k < n; k++) {
    let piv = k;
    let best = Math.abs(M[k * n + k]);
    for (let i = k + 1; i < n; i++) {
      const v = Math.abs(M[i * n + k]);
      if (v > best) {
        best = v;
        piv = i;
      }
    }
    if (piv !== k) {
      for (let j = k; j < n; j++) {
        const tmp = M[k * n + j];
        M[k * n + j] = M[piv * n + j];
        M[piv * n + j] = tmp;
      }
      const tb = x[k];
      x[k] = x[piv];
      x[piv] = tb;
    }
    const diag = M[k * n + k] || 1e-18;
    for (let i = k + 1; i < n; i++) {
      const f = M[i * n + k] / diag;
      for (let j = k + 1; j < n; j++) M[i * n + j] -= f * M[k * n + j];
      x[i] -= f * x[k];
    }
  }
  for (let i = n - 1; i >= 0; i--) {
    let s = x[i];
    for (let j = i + 1; j < n; j++) s -= M[i * n + j] * x[j];
    x[i] = s / (M[i * n + i] || 1e-18);
  }
  return x;
}

function topIds(V, nx) {
  let ymax = -Infinity;
  for (const p of V) ymax = Math.max(ymax, p[1]);
  const ids = [];
  for (let i = 0; i < V.length; i++) if (V[i][1] >= ymax - 1e-8) ids.push(i);
  ids.sort((a, b) => V[a][0] - V[b][0]);
  return { ids, nx };
}

function leftIds(V) {
  let xmin = Infinity;
  for (const p of V) xmin = Math.min(xmin, p[0]);
  const ids = [];
  for (let i = 0; i < V.length; i++) if (V[i][0] <= xmin + 1e-8) ids.push(i);
  return ids;
}

export const SCENE_PRESETS = {
  fall: { cy: 1.75, young: 7500, damping: 0.996, kind: "fall" },
  drop_spin: { cy: 1.95, young: 7000, damping: 0.996, kind: "spin" },
  pinned: { cy: 1.25, young: 4500, damping: 0.9985, kind: "pinned" },
  shake: { cy: 1.35, young: 5500, damping: 0.995, kind: "shake" },
  stretch: { cy: 1.2, young: 6500, damping: 0.997, kind: "stretch" },
  twist: { cy: 1.35, young: 5000, damping: 0.995, kind: "twist" },
};

export function createSheet(opts) {
  const nx = opts.nx ?? 11;
  const ny = opts.ny ?? 5;
  const young = opts.young ?? 8000;
  const poisson = opts.poisson ?? 0.35;
  const density = opts.density ?? 1;
  const width = opts.width ?? 2.4;
  const height = opts.height ?? 1;
  const mu = young / (2 * (1 + poisson));
  const lam = (young * poisson) / ((1 + poisson) * (1 - 2 * poisson));
  const built = meshOf(nx, ny, width, height, opts.cx ?? 0, opts.cy ?? 1);
  const rest = built.V.map((p) => p.slice());
  const V0 = built.V.map((p) => p.slice());
  if (opts.kind === "pinned") {
    let xmin = Infinity;
    let xmax = -Infinity;
    for (const p of rest) {
      xmin = Math.min(xmin, p[0]);
      xmax = Math.max(xmax, p[0]);
    }
    for (const p of V0) {
      const xn = (p[0] - xmin) / Math.max(xmax - xmin, 1e-12);
      const bend = xn * xn * (3 - 2 * xn);
      p[1] += 0.32 * bend;
    }
  }
  const elems = buildElements(rest, built.T);
  const K0 = assembleK0(rest, elems, mu, lam);
  const n = rest.length;
  const mass = new Float64Array(n);
  for (const e of elems) {
    const share = (density * e.area) / 3;
    for (const vi of e.tri) mass[vi] += share;
  }
  const massDiag = new Float64Array(2 * n);
  for (let i = 0; i < n; i++) {
    const m = Math.max(mass[i], 1e-10);
    massDiag[2 * i] = m;
    massDiag[2 * i + 1] = m;
  }
  const handles = chooseHandles(rest, nx, opts.kind);
  const x = new Float64Array(2 * n);
  const v = new Float64Array(2 * n);
  for (let i = 0; i < n; i++) {
    x[2 * i] = V0[i][0];
    x[2 * i + 1] = V0[i][1];
  }
  if (opts.kind === "fall") {
    for (let i = 0; i < n; i++) v[2 * i + 1] = -0.05;
  } else if (opts.kind === "spin") {
    let cx = 0;
    let cy = 0;
    for (const p of V0) {
      cx += p[0];
      cy += p[1];
    }
    cx /= n;
    cy /= n;
    const omega = -3;
    for (let i = 0; i < n; i++) {
      const rx = V0[i][0] - cx;
      const ry = V0[i][1] - cy;
      v[2 * i] = -omega * ry + 0.12;
      v[2 * i + 1] = omega * rx - 0.08;
    }
  } else if (opts.kind === "pinned") {
    let xmin = Infinity;
    let xmax = -Infinity;
    for (const p of V0) {
      xmin = Math.min(xmin, p[0]);
      xmax = Math.max(xmax, p[0]);
    }
    for (let i = 0; i < n; i++) {
      const w = (V0[i][0] - xmin) / Math.max(xmax - xmin, 1e-12);
      v[2 * i] = 0.05 * w;
      v[2 * i + 1] = -0.35 * w;
    }
  }
  const sim = {
    rest, elems, triangles: built.T, K0, massDiag, mu, lam, n,
    handles, x, v, t: 0,
    dt: opts.dt ?? 0.01,
    damping: opts.damping ?? 0.996,
    kc: opts.kc ?? 4e5,
    floor: 0,
    newtonIters: opts.newtonIters ?? 4,
    residual: 0,
    iters: 0,
  };
  applyHandles(sim, 0);
  for (const id of handles.ids) {
    v[2 * id] = 0;
    v[2 * id + 1] = 0;
  }
  return sim;
}

function chooseHandles(rest, nx, kind) {
  if (kind === "pinned") return { kind, ids: leftIds(rest), rest };
  if (kind === "shake") {
    const { ids } = topIds(rest);
    return { kind, ids: ids.slice(0, Math.max(4, Math.floor(nx / 4))), rest, amp: 0.32, freq: 1.1 };
  }
  if (kind === "stretch") {
    const { ids } = topIds(rest);
    const count = Math.max(4, Math.floor(nx / 5));
    return {
      kind, rest, amp: 0.6, freq: 0.65,
      ids: ids.slice(0, count).concat(ids.slice(-count)),
      nLeft: count,
    };
  }
  if (kind === "twist") {
    const { ids } = topIds(rest);
    const count = Math.max(4, Math.floor(nx / 5));
    return {
      kind, rest, amp: 0.3, freq: 0.85,
      ids: ids.slice(0, count).concat(ids.slice(-count)),
      nLeft: count,
    };
  }
  return { kind: "free", ids: [], rest };
}

function handlePositions(sim, t) {
  const h = sim.handles;
  const out = [];
  if (!h.ids.length) return out;
  if (h.kind === "pinned") {
    for (const id of h.ids) out.push([id, h.rest[id][0], h.rest[id][1]]);
    return out;
  }
  if (h.kind === "shake") {
    const ramp = smoothstep(t / 0.6);
    const phase = 2 * Math.PI * h.freq * t;
    for (const id of h.ids) {
      out.push([
        id,
        h.rest[id][0] + ramp * h.amp * Math.sin(phase),
        h.rest[id][1] + ramp * 0.35 * h.amp * Math.cos(1.3 * phase),
      ]);
    }
    return out;
  }
  if (h.kind === "stretch") {
    const phase = 2 * Math.PI * h.freq * t;
    const s = smoothstep(t / 1.4);
    h.ids.forEach((id, k) => {
      let x = h.rest[id][0];
      let y = h.rest[id][1];
      if (k >= h.nLeft) {
        x += h.amp * s;
        y += 0.08 * h.amp * Math.sin(phase) * s;
      }
      out.push([id, x, y]);
    });
    return out;
  }
  if (h.kind === "twist") {
    const ramp = smoothstep(t / 0.6);
    const phase = 2 * Math.PI * h.freq * t;
    h.ids.forEach((id, k) => {
      const sign = k < h.nLeft ? 1 : -1;
      out.push([
        id,
        h.rest[id][0] + sign * ramp * h.amp * Math.sin(phase),
        h.rest[id][1] + sign * ramp * 0.32 * h.amp * Math.cos(phase),
      ]);
    });
  }
  return out;
}

function applyHandles(sim, t) {
  for (const [id, px, py] of handlePositions(sim, t)) {
    sim.x[2 * id] = px;
    sim.x[2 * id + 1] = py;
  }
}

function freeDofs(sim) {
  const marked = new Uint8Array(sim.n * 2);
  for (const id of sim.handles.ids) {
    marked[2 * id] = 1;
    marked[2 * id + 1] = 1;
  }
  const fd = [];
  for (let i = 0; i < marked.length; i++) if (!marked[i]) fd.push(i);
  return fd;
}

export function stepSheet(sim) {
  const dt = sim.dt;
  const N = sim.n * 2;
  const y = new Float64Array(N);
  for (let i = 0; i < sim.n; i++) {
    y[2 * i] = sim.x[2 * i] + dt * sim.v[2 * i];
    y[2 * i + 1] = sim.x[2 * i + 1] + dt * sim.v[2 * i + 1] + dt * dt * GRAVITY;
  }
  const xOld = Float64Array.from(sim.x);
  const targetT = sim.t + dt;
  for (const [id, px, py] of handlePositions(sim, targetT)) {
    y[2 * id] = px;
    y[2 * id + 1] = py;
  }
  const xk = Float64Array.from(y);
  const fd = freeDofs(sim);
  const massOver = sim.massDiag.map((m) => m / (dt * dt));
  let residual = 0;
  let iters = 0;
  const potential = (x) => {
    let val = 0;
    for (let i = 0; i < N; i++) {
      const d = x[i] - y[i];
      val += 0.5 * massOver[i] * d * d;
    }
    val += elastic(x, sim.elems, sim.mu, sim.lam).energy;
    val += contactTerms(x, sim.kc, sim.floor).energy;
    return val;
  };
  for (let it = 0; it < sim.newtonIters; it++) {
    iters = it + 1;
    const eg = elastic(xk, sim.elems, sim.mu, sim.lam);
    const ct = contactTerms(xk, sim.kc, sim.floor);
    const grad = new Float64Array(N);
    for (let i = 0; i < N; i++) grad[i] = massOver[i] * (xk[i] - y[i]) + eg.grad[i] + ct.grad[i];
    if (fd.length === 0) {
      residual = 0;
      break;
    }
    let gnorm = 0;
    for (const i of fd) gnorm += grad[i] * grad[i];
    residual = Math.sqrt(gnorm) / fd.length;
    if (residual < 1e-7) break;
    const m = fd.length;
    const A = new Float64Array(m * m);
    const rhs = new Float64Array(m);
    for (let a = 0; a < m; a++) {
      rhs[a] = -grad[fd[a]];
      for (let b = 0; b < m; b++) {
        A[a * m + b] = sim.K0[fd[a] * N + fd[b]];
      }
      A[a * m + a] += massOver[fd[a]] + ct.diag[fd[a]] + 1e-8;
    }
    const delta = solveDense(A, rhs);
    let oldPhi = potential(xk);
    let directional = 0;
    for (let a = 0; a < m; a++) directional += grad[fd[a]] * delta[a];
    let alpha = 1;
    let accepted = false;
    for (let ls = 0; ls < 18; ls++) {
      const trial = Float64Array.from(xk);
      for (let a = 0; a < m; a++) trial[fd[a]] += alpha * delta[a];
      for (const [id, px, py] of handlePositions(sim, targetT)) {
        trial[2 * id] = px;
        trial[2 * id + 1] = py;
      }
      const trialPhi = potential(trial);
      if (Number.isFinite(trialPhi) && trialPhi <= oldPhi + 1e-4 * alpha * directional) {
        xk.set(trial);
        accepted = true;
        break;
      }
      alpha *= 0.5;
    }
    if (!accepted) {
      for (let a = 0; a < m; a++) xk[fd[a]] += 0.01 * delta[a];
      for (const [id, px, py] of handlePositions(sim, targetT)) {
        xk[2 * id] = px;
        xk[2 * id + 1] = py;
      }
    }
  }
  for (let i = 0; i < N; i++) sim.v[i] = sim.damping * ((xk[i] - xOld[i]) / dt);
  for (const [id, px, py] of handlePositions(sim, targetT)) {
    sim.v[2 * id] = (px - xOld[2 * id]) / dt;
    sim.v[2 * id + 1] = (py - xOld[2 * id + 1]) / dt;
    sim.x[2 * id] = px;
    sim.x[2 * id + 1] = py;
  }
  for (let i = 0; i < N; i++) sim.x[i] = xk[i];
  // Re-apply after the bulk copy so handles win.
  for (const [id, px, py] of handlePositions(sim, targetT)) {
    sim.x[2 * id] = px;
    sim.x[2 * id + 1] = py;
  }
  sim.t += dt;
  sim.residual = residual;
  sim.iters = iters;
  return sim;
}

export function sheetEnergy(sim) {
  const eg = elastic(sim.x, sim.elems, sim.mu, sim.lam);
  const ct = contactTerms(sim.x, sim.kc, sim.floor);
  let kinetic = 0;
  for (let i = 0; i < sim.x.length; i++) kinetic += 0.5 * sim.massDiag[i] * sim.v[i] * sim.v[i];
  return { elastic: eg.energy, contact: ct.energy, kinetic, penetration: ct.maxPen };
}
