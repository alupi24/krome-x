// This test is based on the Zel'dovich nitric oxide mechanism involving
// two reactions and its reverse as discussed in Al-Khateeb+2009. The aim
// of the test is to evaluate the inverse reaction rate coefficients
// applying the thermodynamic method discussed in Grassi+2014 and employed
// in KROME. Initial conditions are n = 1e-3 mols, Volume = 1e3 cm3 and
// T = 4000 K. (KROME-X port of legacy tests/reverse)
#include <cstdio>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const double N_avogadro = constants::N_avogadro;
  const double volume = 1e3;            // cm3
  double x[NSPEC];
  for (double& v : x) v = 1e-3 * N_avogadro / volume;  // 1e-3 mol, cm-3
  const double tend = 1e-4;  // max time, s

  double Tgas = 4e3;  // gas temperature (K)
  double dt = 1e-10;  // initial time-step (s)
  double t = 0.0;     // time (s)
  Params p;           // runtime parameters: defaults
  SolverOpts o;       // solver options: defaults (isothermal)

  std::FILE* f = std::fopen("abundances.dat", "w");
  std::fprintf(f, "#time %s\n", names_header().c_str());
  std::printf("solving...\n");
  for (;;) {
    // call the solver
    krome_test::check(solve(x, Tgas, p, dt, o), "solve");
    dt *= 1.5;
    t += dt;  // as legacy: the time written is t + the next dt
    std::fprintf(f, "%12.3E", t);
    for (int i = 0; i < NSPEC; ++i)
      std::fprintf(f, "%12.3E", x[i] / N_avogadro * volume);  // mol
    std::fprintf(f, "\n");
    if (t > tend) break;
  }
  std::fclose(f);

  std::printf("done!\nTest OK!\nTo plot: python3 plot.py\n");
  return 0;
}
