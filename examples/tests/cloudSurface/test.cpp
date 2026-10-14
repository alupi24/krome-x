// The following test provides a benchmark for gas-grain and grain-surface
// chemistry. Initial conditions from Semenov+2010 (KROME-X port of legacy
// tests/cloudSurface)
#include <cmath>
#include <cstdio>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const double pi = 3.14159265358979;
  const double spy = 3600. * 24. * 365.;  // seconds per year
  double Tgas = 1e1;                      // gas temperature (K)
  const double xH = 2e4;         // Hydrogen density xH = n(H) + 2 n(H2)
  const double CR_ion = 1.3e-17;  // cosmic rays ionisation rate in s^-1
  const double Av_opacity = 1e1;  // extinction coefficient
  const double d2g = 0.01;        // dust to gas ratio
  const double rho0 = 3.0;        // specific density for grain silicates
  const double mmw = 1.43;        // mean molecular weight
  const double adust = 1e-5;      // 0.1 micron average size for grains

  // runtime parameters (legacy krome_set_user_*): opacity, CR rate,
  // Tdust, grain size and its square; xdust below
  Params p;
  p.Av = Av_opacity;
  p.crflux = CR_ion;
  p.Tdust = Tgas;
  p.gsize = adust;
  p.gsize2 = adust * adust;

  // solver options: fixed Tgas (legacy: no thermal processes), the
  // test's -RTOL=1e-6 is the project default; a large step budget for the
  // 1e9 yr steps
  SolverOpts o;
  o.max_steps = 1000000;

  double x[NSPEC];
  for (double& v : x) v = 1e-40;
  // initial densities (Semenov2010)
  x[IDX_He] = 9e-2 * xH;
  x[IDX_H2] = 0.5 * xH;
  x[IDX_Cj] = 1.2e-4 * xH;
  x[IDX_N] = 7.6e-5 * xH;
  x[IDX_O] = 2.56e-4 * xH;
  x[IDX_Sj] = 8.0e-8 * xH;
  x[IDX_Sij] = 8.0e-9 * xH;
  x[IDX_Fej] = 3e-9 * xH;
  x[IDX_Naj] = 2e-9 * xH;
  x[IDX_Mgj] = 7.0e-9 * xH;
  x[IDX_Pj] = 2.0e-10 * xH;
  x[IDX_Clj] = 1.0e-9 * xH;

  // calculate electrons (neutral cloud)
  x[IDX_E] = get_electrons_neutral(x);

  // dust number density (From Semenov 2010), legacy user_xdust
  const double rhogas = get_rho(x);
  const double dmass = 4. / 3. * pi * rho0 * std::pow(adust, 3);
  const double xdust = d2g * rhogas / (dmass * mmw);
  p.xdust = xdust;
  x[IDX_GRAIN0] = xdust;

  double dt = 1.0 * spy;  // time-step (s)
  double t = 0.0;         // initial time (s)

  // output header
  std::FILE* f = std::fopen("abundances.dat", "w");
  std::fprintf(f, "#time %s\n", names_header().c_str());

  for (;;) {
    std::printf("%10s%11.3E%3s\n", "time:", t / spy, "yr");
    krome_test::check(solve(x, Tgas, p, dt, o), "solve");  // call KROME
    t += dt;                        // increase time
    dt = std::fmax(1e2, t / 3.0);  // increase time-step
    krome_test::write_row(f, {t / spy}, x, xH);
    if (t > 1e9 * spy) break;  // exit when overshoot 1e9 years
  }
  std::fclose(f);

  std::printf("Output dump in abundances.dat\nTo plot: python3 plot.py\n"
              "That's all! have a nice day!\n");
  return 0;
}
