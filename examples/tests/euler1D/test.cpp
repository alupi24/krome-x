// THIS TEST EVOLVES 1D HYDRODYNAMICAL EQUATIONS USING HLL+RK2
// CLASSIC SOD SHOCK TEST, 0.125 density factor, 0.1 pressure factor
// NORMALIZED TO CGS UNITS
// (KROME-X port of legacy tests/euler1D: after every hydro step the chemistry
//  and cooling of all the cells are solved in one solve_cells_x call)
// Cells are 0-based here (legacy cell i is index i-1).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

namespace {

// the legacy constants are default-real (single precision) literals: kept
// as float literals, so the grid and the initial state are the same
const double pc2cm = 3.08567758128e18f;  // cm / pc
const double kboltzmann = 1.380649e-16f;  // erg/K
const double pmass = 1.67262192e-24f;     // g

const int ncell = 1024;        // number of grid cells
const int nall = ncell + 2;    // number of grid cells including ghosts
const int idx_rho = 0;         // density in the hydro variables, g/cm3
const int idx_rvx = 1;         // momentum, g/cm2/s
const int idx_energy = 2;      // energy density, erg/cm3
const int nhydro = 3;          // number of hydro variables
const int nscalar = NSPEC;     // advected scalars (the chemical species)
const int nvars = nhydro + nscalar;  // total number of hydro variables
const double gamma_ad = 5.f / 3.f;   // adiabatic index
const double mu = 1.25;        // mean molecular weight (fixed)
const double dx = pc2cm / ncell;  // grid spacing, cm

// the state: u[cell * nvars + variable]
using State = std::vector<double>;

// ************************
// RHS of hydro equations
State rhs(const State& u) {
  State ug(nall * nvars), fg(nall * nvars), dy(ncell * nvars);
  std::vector<double> p(nall), vx(nall), cs(nall), lp(nall), lm(nall);
  std::vector<double> fp(nvars), fm(nvars);

  // set ghost cells
  for (int i = 0; i < ncell; ++i)
    for (int k = 0; k < nvars; ++k) ug[(i + 1) * nvars + k] = u[i * nvars + k];
  for (int k = 0; k < nvars; ++k) {
    ug[k] = u[k];
    ug[(nall - 1) * nvars + k] = u[(ncell - 1) * nvars + k];
  }

  for (int i = 0; i < nall; ++i) {
    const double* g = &ug[i * nvars];
    double* f = &fg[i * nvars];
    // velocity and pressure for later use
    vx[i] = g[idx_rvx] / g[idx_rho];
    p[i] = (g[idx_energy] - g[idx_rho] * vx[i] * vx[i] / 2.) * (gamma_ad - 1.0);

    // fluxes HD from conservation laws
    f[idx_rho] = g[idx_rvx];
    f[idx_rvx] = g[idx_rvx] * vx[i] + p[i];
    f[idx_energy] = (g[idx_energy] + p[i]) * vx[i];

    // chemical species advection, g/cm3
    for (int k = nhydro; k < nvars; ++k) f[k] = vx[i] * g[k];

    // speed of sound, cm/s
    cs[i] = std::sqrt(gamma_ad * p[i] / g[idx_rho]);

    // lambdas
    lp[i] = vx[i] + cs[i];
    lm[i] = vx[i] - cs[i];
  }

  // compute left and right fluxes, then the spatial gradient (RHS)
  for (int i = 1; i < nall - 1; ++i) {
    const double* gm = &ug[(i - 1) * nvars];
    const double* g0 = &ug[i * nvars];
    const double* gp = &ug[(i + 1) * nvars];
    const double* fgm = &fg[(i - 1) * nvars];
    const double* fg0 = &fg[i * nvars];
    const double* fgp = &fg[(i + 1) * nvars];

    double ap = std::max(0.0, std::max(lp[i], lp[i + 1]));
    double am = std::max(0.0, std::max(-lm[i], -lm[i + 1]));
    for (int k = 0; k < nvars; ++k)
      fp[k] = (ap * fg0[k] + am * fgp[k] - ap * am * (gp[k] - g0[k])) /
              (ap + am);

    ap = std::max(0.0, std::max(lp[i], lp[i - 1]));
    am = std::max(0.0, std::max(-lm[i], -lm[i - 1]));
    for (int k = 0; k < nvars; ++k)
      fm[k] = (ap * fgm[k] + am * fg0[k] - ap * am * (g0[k] - gm[k])) /
              (ap + am);

    for (int k = 0; k < nvars; ++k)
      dy[(i - 1) * nvars + k] = -(fp[k] - fm[k]) / dx;
  }
  return dy;
}

// ********************
// use RK2 scheme to advance the solution; xall[cell * NSPEC + species] in
// g/cm3
void solve_hydro(std::vector<double>& rho, std::vector<double>& vx,
                 std::vector<double>& tgas, std::vector<double>& xall,
                 double tend) {
  const double cmax = 1e-1;          // Courant factor
  const double dtmax = dx / 1e6;     // largest time-step allowed, s
  State y0(ncell * nvars), u1(ncell * nvars);
  std::vector<double> cs(ncell);

  // KROME: one parameter set and the solver options for every cell
  // (thermal evolution on: -cooling=H2)
  Params par;
  SolverOpts o;
  o.evolve_T = true;
  o.renormalize_x = true;  // mass conservation, as the legacy krome call
  std::vector<double> x(ncell * NSPEC), rhoc(ncell), T(ncell), dtc(ncell);

  // copy initial conditions for all cells from the arguments
  for (int i = 0; i < ncell; ++i) {
    const double p = rho[i] * kboltzmann * tgas[i] / mu / pmass;  // erg/cm3
    cs[i] = std::sqrt(gamma_ad * p / rho[i]);  // speed of sound, cm/s
    double* u = &y0[i * nvars];
    u[idx_rho] = rho[i];                   // density, g/cm3
    u[idx_rvx] = rho[i] * vx[i];           // momentum, g/cm2/s
    u[idx_energy] = p / (gamma_ad - 1.0) + rho[i] * vx[i] * vx[i] / 2.;
    for (int k = 0; k < nscalar; ++k) u[nhydro + k] = xall[i * NSPEC + k];
  }

  double t = 0.0;
  for (;;) {
    // compute time-step (cs: the initial sound speeds, as legacy)
    double vmax = 0.0, vmin = 0.0;
    for (int i = 0; i < ncell; ++i) {
      const double v = std::abs(y0[i * nvars + idx_rvx] / y0[i * nvars + idx_rho]);
      vmax = std::max(vmax, v + cs[i]);
      vmin = std::max(vmin, v - cs[i]);
    }
    vmax = std::max(vmax + 1e-40, vmin + 1e-40);
    double dt = std::min(dx / vmax * cmax, dtmax);
    dt = std::min(dt, tend - t);

    // advance RK2
    const State k1 = rhs(y0);
    for (int j = 0; j < ncell * nvars; ++j) u1[j] = y0[j] + 0.5 * dt * k1[j];
    const State k2 = rhs(u1);
    for (int j = 0; j < ncell * nvars; ++j) y0[j] += dt * k2[j];

    // KROME: every cell in one call (mass fractions, legacy -useX)
    for (int i = 0; i < ncell; ++i) {
      const double* u = &y0[i * nvars];
      rhoc[i] = u[idx_rho];
      for (int k = 0; k < NSPEC; ++k)
        x[i * NSPEC + k] = u[nhydro + k] / u[idx_rho];  // mass fractions
      const double p1 = (u[idx_energy] - u[idx_rvx] * u[idx_rvx] / u[idx_rho] / 2.) *
                        (gamma_ad - 1.0);  // pressure, erg/cm3
      T[i] = p1 / u[idx_rho] / kboltzmann * mu * pmass;  // tgas, K
      dtc[i] = dt;
    }
    krome_test::check(solve_cells_x(ncell, x.data(), rhoc.data(), T.data(),
                                    dtc.data(), par, o),
                      "solve_cells_x");
    for (int i = 0; i < ncell; ++i) {
      double* u = &y0[i * nvars];
      for (int k = 0; k < NSPEC; ++k)  // update abundances, g/cm3
        u[nhydro + k] = x[i * NSPEC + k] * u[idx_rho];
      const double p1 = T[i] / mu / pmass * kboltzmann * u[idx_rho];
      // update energy density, erg/cm3
      u[idx_energy] = p1 / (gamma_ad - 1.0) +
                      u[idx_rvx] * u[idx_rvx] / u[idx_rho] / 2.;
    }

    // advance time
    t += dt;
    std::printf("%20.12E\n", t / tend);
    if (t >= tend) break;
  }

  // copy final solutions to the arguments
  for (int i = 0; i < ncell; ++i) {
    const double* u = &y0[i * nvars];
    rho[i] = u[idx_rho];
    vx[i] = u[idx_rvx] / rho[i];
    const double p = (u[idx_energy] - rho[i] * vx[i] * vx[i] / 2.) * (gamma_ad - 1.0);
    tgas[i] = p / rho[i] / kboltzmann * mu * pmass;
    for (int k = 0; k < nscalar; ++k) xall[i * NSPEC + k] = u[nhydro + k];
  }
}

}  // namespace

// ****************
// driver code
int main() {
  const double spy = 3600.f * 24.f * 365;  // s / yr
  std::vector<double> rho(ncell), vx(ncell), tgas(ncell), xall(ncell * NSPEC);

  // density left and right conditions, g/cm3
  const double rhoL = 1e-20;
  const double rhoR = 0.125 * rhoL;
  // temperature left and right conditions (to have pressure factor 0.1), K
  const double TgasL = 1.4e3;
  const double TgasR = 0.1f * rhoL * TgasL / rhoR;

  // final integration time, s
  const double tend = 2.45e4 * spy;

  // set initial conditions, gas initially at rest (cm/s)
  for (int i = 0; i < ncell; ++i) {
    vx[i] = 0.0;
    rho[i] = i < ncell / 2 ? rhoL : rhoR;
    tgas[i] = i < ncell / 2 ? TgasL : TgasR;
  }

  // set initial chemistry conditions (g/cm3)
  for (int i = 0; i < ncell; ++i) {
    double* x = &xall[i * NSPEC];
    for (int k = 0; k < NSPEC; ++k) x[k] = 1e-40;
    x[IDX_H] = 0.75615f * rho[i];
    x[IDX_E] = 4.4983e-8 * rho[i];
    x[IDX_Hj] = 8.1967e-5 * rho[i];
    x[IDX_He] = 0.24375f * rho[i];
    x[IDX_H2] = 1.5123e-6 * rho[i];
  }

  // do the shock test
  solve_hydro(rho, vx, tgas, xall, tend);

  // write output to file
  std::FILE* f = std::fopen("profiles_final.dat", "w");
  std::fprintf(f, "#x rho vx Tgas H2 H- e-\n");
  for (int i = 0; i < ncell; ++i)
    std::fprintf(f, "%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E\n",
                 (i + 1) * dx / pc2cm, rho[i], vx[i], tgas[i],
                 xall[i * NSPEC + IDX_H2], xall[i * NSPEC + IDX_Hk],
                 xall[i * NSPEC + IDX_E]);
  std::fclose(f);

  std::printf(" To plot: python3 plot.py\n done, bye!\n");
  return 0;
}
