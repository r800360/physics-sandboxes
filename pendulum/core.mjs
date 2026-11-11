/** Same RK4 step as pendulum.py. Theta is from the downward vertical. */

export function rk4Step(y, dt, g, L) {
  const f = (yy) => [yy[1], -(3 * g / (2 * L)) * Math.sin(yy[0])];
  const k1 = f(y);
  const k2 = f([y[0] + 0.5 * dt * k1[0], y[1] + 0.5 * dt * k1[1]]);
  const k3 = f([y[0] + 0.5 * dt * k2[0], y[1] + 0.5 * dt * k2[1]]);
  const k4 = f([y[0] + dt * k3[0], y[1] + dt * k3[1]]);
  return [
    y[0] + (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
    y[1] + (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]),
  ];
}

export function simpleStep(y, dt, g, L) {
  const f = (yy) => [yy[1], -(g / L) * Math.sin(yy[0])];
  const k1 = f(y);
  const k2 = f([y[0] + 0.5 * dt * k1[0], y[1] + 0.5 * dt * k1[1]]);
  const k3 = f([y[0] + 0.5 * dt * k2[0], y[1] + 0.5 * dt * k2[1]]);
  const k4 = f([y[0] + dt * k3[0], y[1] + dt * k3[1]]);
  return [
    y[0] + (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
    y[1] + (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]),
  ];
}

export function energy(theta, omega, m, L, g) {
  const inertia = (m * L * L) / 3;
  const kinetic = 0.5 * inertia * omega * omega;
  const potential = m * g * (L / 2) * (1 - Math.cos(theta));
  return kinetic + potential;
}

export function simulate(theta0, omega0, L, g, T, dt) {
  const n = Math.round(T / dt);
  const theta = new Float64Array(n + 1);
  const omega = new Float64Array(n + 1);
  const eng = new Float64Array(n + 1);
  let y = [theta0, omega0];
  theta[0] = y[0];
  omega[0] = y[1];
  eng[0] = energy(y[0], y[1], 1, L, g);
  for (let i = 0; i < n; i++) {
    y = rk4Step(y, dt, g, L);
    theta[i + 1] = y[0];
    omega[i + 1] = y[1];
    eng[i + 1] = energy(y[0], y[1], 1, L, g);
  }
  return { theta, omega, energy: eng, dt, n };
}
