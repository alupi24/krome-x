// This test is designed to check if CO, OH, H2O, and HCN cooling produce
// reasonable temperature evolutions. WARNING: do not use the network for
// production (KROME-X port of legacy tests/cooling_molecules)
#include <cstdio>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const double spy = 365.0 * 24.0 * 3600.0;  // legacy krome_seconds_per_year

  // solver options: thermal evolution on
  SolverOpts o;
  o.evolve_T = true;

  // initial conditions
  const double ntot = 1e6;      // total density, 1/cm3
  double Tgas = 3e2;            // temperature, K
  const double tend = 1e8 * spy;  // integration time, s

  // species initialization, 1/cm3
  double x[NSPEC];
  for (double& v : x) v = 1e-40;
  x[IDX_H] = ntot;
  x[IDX_H2] = 1e-6 * ntot;
  x[IDX_CN] = 7e-5 * ntot;
  x[IDX_HCO] = 1e-4 * ntot;
  x[IDX_O] = 1e-4 * ntot;

  // set cosmic rays (legacy krome_set_user_crate)
  Params p;
  p.crate = 1.3e-17;

  // open file to write output
  std::FILE* f = std::fopen("data.out", "w");
  std::fprintf(f, "#time Tgas %s\n", names_header().c_str());

  double t = 0.0, dt = spy;

  // pseudo time-dependent loop
  for (;;) {
    dt *= 1.1;
    // solve the chemistry
    krome_test::check(solve(x, Tgas, p, dt, o), "solve");
    t += dt;
    krome_test::write_row(f, {t / spy, Tgas}, x, ntot);
    if (t >= tend) break;
  }
  std::fclose(f);

  std::printf("To plot: python3 plot.py\nThat's all! have a nice day!\n");
  return 0;
}
