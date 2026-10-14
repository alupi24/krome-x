// HELLO KROME test evolves the network
//   FK1 -> FK2 -> FK3
// where FK* are fake species (KROME-X port of legacy tests/hello)
#include <cstdio>
#include "krome.hpp"            // the C++ API: solve(), queries, constants
#include "krome_test_utils.hpp" // check, write_row
using namespace krome;

int main() {
  double n[NSPEC];
  for (double& v : n) v = 0.0;  // default abundances (number density)
  n[IDX_FK1] = 1.0;             // FK1 initial abundance (number density)

  double Tgas = 1e2;  // gas temperature, not used (K)
  double dt = 1e-5;   // time-step (arbitrary)
  double t = 0.0;     // time (arbitrary)
  Params p;           // runtime parameters (none for this network)
  SolverOpts o;       // solver options: defaults (isothermal)

  std::FILE* f = std::fopen("abundances.dat", "w");
  std::fprintf(f, "#time %s\n", names_header().c_str());
  for (;;) {
    dt *= 1.1;  // increase timestep
    krome_test::check(solve(n, Tgas, p, dt, o), "solve");
    krome_test::write_row(f, {t}, n);
    t += dt;           // increase time
    if (t > 5.0) break;  // break loop after 5 time units
  }
  std::fclose(f);

  std::printf("Test OK!\nTo plot: python3 plot.py\n");
  return 0;
}
