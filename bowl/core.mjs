/**
 * Asymmetric dumbbell in y = alpha x^2.
 * q = [x, y, theta], v = [xdot, ydot, thetadot].
 * Matches the Python integrator, including sample order.
 */

function rot(theta, rx, ry) {
  const c = Math.cos(theta);
  const s = Math.sin(theta);
  return [c * rx - s * ry, s * rx + c * ry];
}

function dRot(theta, rx, ry) {
  const c = Math.cos(theta);
  const s = Math.sin(theta);
  return [-s * rx - c * ry, c * rx - s * ry];
}

export function contactSamples(p) {
  const pts = [];
  for (let i = 0; i < p.nDisk; i++) {
    const a = (i * 2 * Math.PI) / p.nDisk;
    pts.push([-p.l1 + p.r1 * Math.cos(a), p.r1 * Math.sin(a)]);
  }
  for (let i = 0; i < p.nDisk; i++) {
    const a = (i * 2 * Math.PI) / p.nDisk;
    pts.push([p.l2 + p.r2 * Math.cos(a), p.r2 * Math.sin(a)]);
  }
  for (let i = 0; i < 9; i++) {
    const x = -p.l1 + ((p.l1 + p.l2) * i) / 8;
    pts.push([x, 0]);
  }
  return pts;
}

export function defaultParams() {
  return {
    m1: 1, m2: 2, l1: 0.8, l2: 0.4, g: 9.81,
    alpha: 0.35, k: 1e4, r1: 0.18, r2: 0.26, nDisk: 24,
    dt: 0.002, T: 12,
  };
}

function massInertia(p) {
  const M = p.m1 + p.m2;
  const I = p.m1 * p.l1 * p.l1 + p.m2 * p.l2 * p.l2;
  return { M, I };
}

function gapAndGrad(q, r, p) {
  const X = world(q, r[0], r[1]);
  const slope = 2 * p.alpha * X[0];
  const dth = dRot(q[2], r[0], r[1]);
  const g = p.alpha * X[0] * X[0] - X[1];
  // ∇g = b'(x) * dX/dq_x-row minus dX/dq_y-row
  return {
    g,
    grad: [slope * 1 - 0, slope * 0 - 1, slope * dth[0] - dth[1]],
  };
}

export function world(q, rx, ry) {
  const d = rot(q[2], rx, ry);
  return [q[0] + d[0], q[1] + d[1]];
}

export function gradPotential(q, samples, p) {
  const { M } = massInertia(p);
  const grad = [0, M * p.g, 0];
  for (const r of samples) {
    const { g, grad: gg } = gapAndGrad(q, r, p);
    if (g > 0) {
      const s = p.k * g;
      grad[0] += s * gg[0];
      grad[1] += s * gg[1];
      grad[2] += s * gg[2];
    }
  }
  return grad;
}

function acceleration(q, samples, p) {
  const { M, I } = massInertia(p);
  const g = gradPotential(q, samples, p);
  return [-g[0] / M, -g[1] / M, -g[2] / I];
}

function rhs(z, samples, p) {
  const q = z.slice(0, 3);
  const v = z.slice(3);
  const a = acceleration(q, samples, p);
  return [v[0], v[1], v[2], a[0], a[1], a[2]];
}

function add(z, s, k) {
  return z.map((v, i) => v + s * k[i]);
}

export function rk4Step(q, v, dt, samples, p) {
  let z = [q[0], q[1], q[2], v[0], v[1], v[2]];
  const k1 = rhs(z, samples, p);
  const k2 = rhs(add(z, 0.5 * dt, k1), samples, p);
  const k3 = rhs(add(z, 0.5 * dt, k2), samples, p);
  const k4 = rhs(add(z, dt, k3), samples, p);
  z = z.map((v, i) => v + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
  return [z.slice(0, 3), z.slice(3)];
}

export function symplecticStep(q, v, dt, samples, p) {
  const a = acceleration(q, samples, p);
  const vn = [v[0] + dt * a[0], v[1] + dt * a[1], v[2] + dt * a[2]];
  const qn = [q[0] + dt * vn[0], q[1] + dt * vn[1], q[2] + dt * vn[2]];
  return [qn, vn];
}

export function totalEnergy(q, v, samples, p) {
  const { M, I } = massInertia(p);
  let U = M * p.g * q[1];
  for (const r of samples) {
    const { g } = gapAndGrad(q, r, p);
    if (g > 0) U += 0.5 * p.k * g * g;
  }
  const K = 0.5 * M * (v[0] * v[0] + v[1] * v[1]) + 0.5 * I * v[2] * v[2];
  return K + U;
}

export function maxPenetration(q, samples, p) {
  let m = 0;
  for (const r of samples) m = Math.max(m, gapAndGrad(q, r, p).g);
  return Math.max(0, m);
}

export function simulate(method, p, q0, v0) {
  const samples = contactSamples(p);
  const n = Math.round(p.T / p.dt);
  const keep = 5;
  const frames = [];
  let q = q0.slice();
  let v = v0.slice();
  const push = (step) => {
    frames.push({
      t: step * p.dt,
      q: q.slice(),
      e: totalEnergy(q, v, samples, p),
      pen: maxPenetration(q, samples, p),
    });
  };
  push(0);
  const stepFn = method === "rk4" ? rk4Step : symplecticStep;
  for (let s = 1; s <= n; s++) {
    [q, v] = stepFn(q, v, p.dt, samples, p);
    if (s % keep === 0) push(s);
  }
  return frames;
}
