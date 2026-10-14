// ****************SEDOV-LIKE SHOCK 1D LAGRANGIAN TEST***************
// REFER TO: Bodenheimer P., Laughlin G., Rozyczka M., Yorke H.
//  2006, Numerical Methods in Astrophysics:
//  An Introduction. Series in Astronomy and Astrophysics,
//  Taylor&Francis
// (KROME-X port of legacy tests/shock1D: the chemistry of all the shells is
//  solved in one solve_cells_x call per hydro step)
// *******************************************************************
//
// This Lagrangian code follows the adiabatic expansion of a hot spherical
// bubble in a uniform ambient medium, in spherical geometry, with
// primordial chemistry (no cooling: the chemistry runs at fixed Tgas).
// Shells are 0-based here (legacy shell ix is index ix-1).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

namespace {

const int nx = 100, il = 20;  // shells, shells of the inner ejecta
const double pi = 3.1415, pi43 = 4.0 / 3.0 * pi;
const double spy = 365. * 24. * 3600.;  // s / yr (the legacy driver's)
// physical constants (kromex/constants.py)
const double kb = constants::boltzmann_erg;  // erg / K
const double pmass = constants::p_mass;      // g
const double gravity = constants::gravity;   // cm3 / g / s2

// the hydro state (legacy module commons)
struct Grid {
  double cflfactor = 0, dfactor = 0, gamma = 0, dt12 = 0, dt = 0, t = 0;
  double fmratio = 0, tmax = 0, radius = 0, rhogas = 0;
  std::vector<double> Tgas, dm12, rho, eps, p, u, fm, r, v, at12, rold, a,
      dr12, dm, w, ak12, massr;
  std::vector<double> x;  // mass fractions x[ix * NSPEC + i]
  Grid()
      : Tgas(nx), dm12(nx), rho(nx), eps(nx), p(nx), u(nx), fm(nx), r(nx),
        v(nx), at12(nx), rold(nx), a(nx), dr12(nx), dm(nx), w(nx), ak12(nx),
        massr(nx), x(nx * NSPEC) {}
};

// *****INITIALIZATION OF THE CHEMICAL ABUNDANCES (mass fractions)*******
void initchem(Grid& g) {
  for (int ix = 0; ix < nx; ++ix) {
    double* x = &g.x[ix * NSPEC];
    for (int i = 0; i < NSPEC; ++i) x[i] = 0.0;
    x[IDX_H] = 0.9225;  // H
    x[IDX_E] = 1e-4;    // E
    x[IDX_Hj] = 1e-4;   // H+
    x[IDX_D] = 1e-20;   // D
    x[IDX_Dj] = 1e-20;  // D+
    x[IDX_He] = 0.0972; // He
    x[IDX_Hej] = 1e-20; // He+
    x[IDX_H2j] = 1e-20; // H2+
    x[IDX_H2] = 1e-5;   // H2
    x[IDX_HD] = 1e-8;   // HD
    x[IDX_Hk] = 1e-20;  // H-
    x[IDX_Hejj] = 1e-20;  // He++
    // normalize
    double s = 0.0;
    for (int i = 0; i < NSPEC; ++i) s += x[i];
    for (int i = 0; i < NSPEC; ++i) x[i] /= s;
  }
}

// ******TIME STEP EVALUATION********
void tmsc(Grid& g) {
  double dtc = 1e99;
  for (int ix = 0; ix < nx - 1; ++ix)
    dtc = std::min(dtc, g.dr12[ix] / (std::abs(g.u[ix]) +
                                      std::sqrt(g.gamma * g.eps[ix])));

  // courant
  dtc = g.cflfactor * dtc;
  if (g.t + dtc > g.tmax) dtc = g.tmax - g.t;

  // diffusion limit
  double dtd = 1e-99;
  for (int ix = 0; ix < nx - 1; ++ix)
    dtd = std::max(dtd, std::abs(g.at12[ix + 1] * g.u[ix + 1] -
                                 g.at12[ix] * g.u[ix]) /
                            (g.v[ix + 1] - g.v[ix]));

  dtd = 0.5 / dtd / g.dfactor;
  dtc = std::min(dtc, dtd);
  g.dt = 0.5 * (g.dt12 + dtc);
  g.dt12 = dtc;
}

// *****HYDRODINAMIC INITIAL CONDITIONS********
void inicond(Grid& g, std::FILE* f22) {
  // Courant factor, must be smaller than 1
  g.cflfactor = 0.1;
  g.fmratio = 0.1;

  // diffusion factor, must be greater that 1
  g.dfactor = 1.5 * 1.5;

  const double rSN = g.radius;  // cm

  // inner: high pressure ejecta
  const double rhoej =
      (std::pow(double(nx) / il, 3) - 1.0) * g.fmratio * g.rhogas;
  const double dmej = pi43 * std::pow(double(il) / nx * rSN, 3) * rhoej / il;
  double efac = 1e6 * kb / pmass;

  std::printf(" ***********************\n inner: high pressure ejecta\n");
  std::printf(" rho (g/cm3): %E\n dm (g): %E\n eps (erg/g): %E\n", rhoej,
              dmej, efac);

  for (int ix = 0; ix < il; ++ix) {
    g.dm12[ix] = dmej;
    g.rho[ix] = rhoej;
    g.eps[ix] = efac;
    g.p[ix] = (g.gamma - 1.0) * g.rho[ix] * g.eps[ix];
    g.u[ix] = 0.0;
  }

  // outer: low pressure ambient medium
  const double rhoamb = g.rhogas;
  const double dmamb = pi43 * std::pow(rSN, 3) *
                       (1.0 - std::pow(double(il) / nx, 3)) / (nx - il) *
                       rhoamb;
  efac = 1e1 * kb / pmass;
  std::printf(" ***********************\n outer: low pressure ambient medium\n");
  std::printf(" rho (g/cm3): %E\n dm (g): %E\n eps (erg/g): %E\n", rhoamb,
              dmamb, efac);

  for (int ix = il; ix < nx; ++ix) {
    g.dm12[ix] = dmamb;
    g.rho[ix] = rhoamb;
    g.eps[ix] = efac;
    g.p[ix] = (g.gamma - 1.0) * g.rho[ix] * g.eps[ix];
    g.u[ix] = 0.0;
  }

  // initialize integral mass
  g.fm[0] = g.dm12[0];
  for (int ix = 1; ix < nx; ++ix) g.fm[ix] = g.fm[ix - 1] + g.dm12[ix];
  // initialize shell masses
  for (int ix = 1; ix < nx; ++ix) g.dm[ix] = 0.5 * (g.dm12[ix] + g.dm12[ix - 1]);

  // initialize shell radius, volume, surface
  g.r[0] = 0.0;
  g.v[0] = 0.0;
  for (int ix = 1; ix < nx; ++ix) {
    g.v[ix] = g.v[ix - 1] + g.dm12[ix - 1] / g.rho[ix - 1];
    // legacy (1./3.) is a single-precision exponent
    g.r[ix] = std::pow(g.v[ix] / pi43, double(1.f / 3.f));
    g.a[ix] = 4.0 * pi * g.r[ix] * g.r[ix];
  }
  // initialize shell sizes
  for (int ix = 0; ix < nx - 1; ++ix) g.dr12[ix] = g.r[ix + 1] - g.r[ix];

  // dump data
  double mtot = 0.0;
  for (double m : g.dm) mtot += m;
  std::printf(" ***********************\n Rtot: %E\n Vtot: %E\n Mtot: %E\n",
              g.r[nx - 1], g.v[nx - 1], mtot);

  std::fprintf(f22, "#radius dr12 dm12 eps rho\n");
  for (int ix = 0; ix < nx; ++ix)
    std::fprintf(f22, "%17.8E%17.8E%17.8E%17.8E%17.8E\n", g.r[ix], g.dr12[ix],
                 g.dm12[ix], g.eps[ix], g.rho[ix]);
  g.t = 0.0;
  std::printf(" ***********************\n");
}

// one hydro row: idx time radius rho velocity pressure Tgas
void hydro_row(std::FILE* f, const Grid& g, int ix) {
  std::fprintf(f, "%8d%17.8E%17.8E%17.8E%17.8E%17.8E%17.8E\n", ix + 1, g.t,
               g.r[ix], g.rho[ix], std::abs(g.u[ix]), g.p[ix], g.Tgas[ix]);
}

// one chemistry row: idx, the lead values, the mass fractions
void chem_row(std::FILE* f, int ix, std::initializer_list<double> lead,
              const double* x) {
  std::fprintf(f, "%8d", ix + 1);
  krome_test::write_row(f, lead, x);
}

}  // namespace

// *****************START PROGRAM SEDOV****************************
int main() {
  Grid g;
  SolverOpts o;  // solver options: defaults (isothermal)
  o.renormalize_x = true;  // mass conservation, as the legacy krome call
  Params par;    // one parameter set for every shell (network defaults)

  g.gamma = 5.0 / 3.0;          // adiabatic index
  const int nsteps = int(1e4);  // maximum number of time-steps
  const int ndump = 10;         // dump interval
  for (double& T : g.Tgas) T = 1e3;  // default gas temp (K)
  g.tmax = spy * 5e3;           // maximum integration time (s)
  const double q = 2.0;         // artificial viscosity parameter
  g.radius = 3.08568025e18 * 1.0;  // radius of the box (cm)
  g.rhogas = 1e-24;             // gas density (g/cm3)

  // write input to check
  std::printf("%20s%11d\n", "nsteps:", nsteps);
  std::printf("%20s%11.3E\n", "tmax (s):", g.tmax);
  std::printf("%20s%11.3E\n", "Tgas (K):", 1e3);
  std::printf("%20s%11.3E\n", "q:", q);
  std::printf("%20s%11.3E\n", "radius (cm):", g.radius);
  std::printf("%20s%11.3E\n", "density (g/cm3):", g.rhogas);
  std::printf("%20s%11.3E\n", "density (1/cm3):", g.rhogas / pmass);

  std::FILE* f22 = std::fopen("initial.dat", "w");
  std::FILE* f23 = std::fopen("hydro.dat", "w");
  std::FILE* f33 = std::fopen("chemistry.dat", "w");
  inicond(g, f22);  // set initial conditions
  initchem(g);      // init chemistry, abundances in mass fraction
  std::fclose(f22);

  // output headers of the dumps
  const std::string names = names_header();
  std::fprintf(f23, "#idx time radius rho velocity pressure Tgas\n");
  std::fprintf(f33, "#idx time %s\n", names.c_str());

  std::vector<double> dtcell(nx);

  // begin loop over time steps MAIN LOOP
  for (int istep = 1; istep <= nsteps; ++istep) {
    // determine time step
    tmsc(g);
    const double dt = g.dt, dt12 = g.dt12;

    // update velocities (eqn 6.42)
    for (int ix = 1; ix < nx; ++ix) {
      g.u[ix] = g.u[ix] - g.a[ix] * (g.p[ix] - g.p[ix - 1]) * dt / g.dm[ix] -
                0.5 * (g.w[ix] * (3.0 * g.ak12[ix] - g.a[ix]) -
                       g.w[ix - 1] * (3.0 * g.ak12[ix - 1] - g.a[ix])) *
                    dt / g.dm[ix];
      g.u[ix] = g.u[ix] - gravity * g.massr[ix] / (g.r[ix] * g.r[ix]) * dt;
    }

    g.u[0] = 0.0;  // velocity of the first shell

    // update radii, surfaces and volumes (eqn 6.41_2)
    g.rold = g.r;
    for (int ix = 0; ix < nx; ++ix) g.r[ix] = g.rold[ix] + g.u[ix] * dt12;

    // control shell position
    if (*std::min_element(g.r.begin(), g.r.end()) < 0.0) {
      std::printf("r(ix)<0\n");
      return 1;
    }

    // update shell radius
    for (int ix = 0; ix < nx - 1; ++ix) g.dr12[ix] = g.r[ix + 1] - g.r[ix];

    // update surface and volume
    for (int ix = 0; ix < nx; ++ix) {
      const double rh = 0.5 * (g.r[ix] + g.rold[ix]);
      g.at12[ix] = 4.0 * pi * rh * rh;
      g.a[ix] = 4.0 * pi * g.r[ix] * g.r[ix];
      g.v[ix] = pi43 * std::pow(g.r[ix], 3);
    }

    for (int ix = 0; ix < nx - 1; ++ix)
      g.ak12[ix] = 0.5 * (g.at12[ix + 1] + g.at12[ix]);

    // update densities (6.41_3)
    for (int ix = 0; ix < nx - 1; ++ix)
      g.rho[ix] = g.dm12[ix] / (g.v[ix + 1] - g.v[ix]);
    g.rho[nx - 1] = g.rho[nx - 2];

    // update shells integral mass
    g.massr[0] = 0.0;
    for (int ix = 1; ix < nx; ++ix)
      g.massr[ix] = g.massr[ix - 1] +
                    pi43 * (std::pow(g.r[ix], 3) - std::pow(g.r[ix - 1], 3)) *
                        g.rho[ix];

    // artificial viscosity
    for (int ix = 0; ix < nx - 1; ++ix)
      g.w[ix] = -q * q * g.rho[ix] * std::abs(g.u[ix + 1] - g.u[ix]) *
                (g.u[ix + 1] * (1.0 - g.at12[ix + 1] / 3.0 / g.ak12[ix]) -
                 g.u[ix] * (1.0 - g.at12[ix] / 3.0 / g.ak12[ix]));

    // no viscosity for postive velocity divergence
    for (int ix = 0; ix < nx - 1; ++ix)
      if (g.u[ix + 1] > g.u[ix]) g.w[ix] = 0.0;

    // update internal energies and pressures (6.43)
    for (int ix = 0; ix < nx - 1; ++ix) {
      const double aux = g.eps[ix] - g.p[ix] *
                                         (g.at12[ix + 1] * g.u[ix + 1] -
                                          g.at12[ix] * g.u[ix]) *
                                         dt12 / g.dm12[ix];
      g.p[ix] = 0.5 * (g.p[ix] + (g.gamma - 1.0) * g.rho[ix] * aux);
    }

    // internal energies
    for (int ix = 0; ix < nx - 1; ++ix)
      g.eps[ix] = g.eps[ix] - g.p[ix] *
                                  (g.at12[ix + 1] * g.u[ix + 1] -
                                   g.at12[ix] * g.u[ix]) *
                                  dt12 / g.dm12[ix];

    // contribution from artificial viscosity
    for (int ix = 0; ix < nx - 1; ++ix) {
      g.eps[ix] = g.eps[ix] - 0.5 * g.w[ix] * dt12 / g.dm12[ix] *
                                  (g.u[ix + 1] * (3.0 * g.ak12[ix] - g.at12[ix + 1]) -
                                   g.u[ix] * (3.0 * g.ak12[ix] - g.at12[ix]));
      g.Tgas[ix] = (g.gamma - 1.0) * g.eps[ix] / kb * pmass;
    }

    // DO CHEMISTRY: ONE KROME CALL FOR ALL THE SHELLS
    // x in mass fraction, density in g/cm3, the same dt for every shell
    for (int ix = 0; ix < nx; ++ix) dtcell[ix] = dt;
    krome_test::check(solve_cells_x(nx, g.x.data(), g.rho.data(),
                                    g.Tgas.data(), dtcell.data(), par, o),
                      "solve_cells_x");

    // update internal energy with new temperature
    for (int ix = 0; ix < nx; ++ix)
      g.eps[ix] = g.Tgas[ix] / (g.gamma - 1.0) / pmass * kb;

    // update pressure
    for (int ix = 0; ix < nx - 1; ++ix)
      g.p[ix] = (g.gamma - 1.0) * g.rho[ix] * g.eps[ix];

    // update last shell
    g.p[nx - 1] = g.p[nx - 2];
    g.eps[nx - 1] = g.eps[nx - 2];

    g.t += dt12;  // update time

    // dump hydrodinamics and chemistry
    if (istep % ndump == 0) {
      std::printf("step:%6d; t:%10.2E; %6.1f%%\n", istep, g.t / spy,
                  g.t / g.tmax * 100);
      for (int ix = 0; ix < nx; ++ix) {
        hydro_row(f23, g, ix);
        chem_row(f33, ix, {g.t}, &g.x[ix * NSPEC]);
      }
      std::fprintf(f23, "\n");
      std::fprintf(f33, "\n");
    }

    // end loop over time steps when overshoot tmax
    if (g.t >= g.tmax) break;
  }
  std::fclose(f23);
  std::fclose(f33);

  // dump final conditions and chemistry (mass fractions)
  std::FILE* f24 = std::fopen("hydro_final.dat", "w");
  std::FILE* f34 = std::fopen("chemistry_final.dat", "w");
  std::fprintf(f24, "#idx time radius rho velocity pressure Tgas\n");
  std::fprintf(f34, "#idx time radius %s\n", names.c_str());
  for (int ix = 0; ix < nx - 1; ++ix) {
    hydro_row(f24, g, ix);
    chem_row(f34, ix, {g.t, g.r[ix]}, &g.x[ix * NSPEC]);
  }
  std::fclose(f24);
  std::fclose(f34);

  std::printf(" finish!\nTo plot: python3 plot.py\n");
  return 0;
}
