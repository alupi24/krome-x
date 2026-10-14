// Test based on [R2008] Reinhardt, V. et al. JPhysChemA 112, 1712 (2008)
// Slow manifold trajectories example.
// This program evolves the chemical network for jmax*kmax initial
// conditions of H2O and H2 (KROME-X port of legacy tests/slowmanifold)
#include <cstdio>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  // grid points
  const int jmax = 14, kmax = 7;
  double x[NSPEC], xi[NSPEC], xold[NSPEC];

  double Tgas = 1.2e3;  // gas Temperature (K), dummy
  Params p;             // runtime parameters: defaults
  SolverOpts o;         // solver options: defaults (isothermal)
  std::printf("solving...\n");

  std::FILE* f66 = std::fopen("abundances.dat", "w");
  std::FILE* f77 = std::fopen("initial.dat", "w");
  std::fprintf(f66, "#time %s\n", names_header().c_str());
  std::fprintf(f77, "#time %s\n", names_header().c_str());

  // loop on H2O points
  for (int j = 1; j <= jmax; ++j) {
    // loop on H2 points
    for (int k = 1; k <= kmax; ++k) {
      // set initial conditions
      for (double& v : xi) v = 0.0;
      xi[IDX_OH] = 0.0;  // OH
      xi[IDX_O] = 0.0;   // O
      xi[IDX_H2O] = (j - 1) * (0.65 - 0.05) / (jmax - 1) + 0.05;  // H2O
      xi[IDX_H2] = (k - 1) * (0.9 - 0.3) / (kmax - 1) + 0.3;      // H2
      // O2 and H from Eqn.(16) of [R2008]
      xi[IDX_O2] = (1.0 - xi[IDX_H2O] - xi[IDX_OH] - xi[IDX_O]) / 2.0;
      xi[IDX_H] = 2.0 - xi[IDX_OH] - 2.0 * xi[IDX_H2O] - 2.0 * xi[IDX_H2];
      // use only initial values on the lower part of the H2-H2O space
      if (xi[IDX_H2] > .96 - xi[IDX_H2O]) continue;
      double dt = 1e-10;  // initial time-step (s)
      double t = 0.0;     // initial time (s)
      for (int i = 0; i < NSPEC; ++i) x[i] = xi[i];  // init abundances
      krome_test::write_row(f77, {t}, x);  // initial conditions
      // loop on time
      for (;;) {
        for (int i = 0; i < NSPEC; ++i) xold[i] = x[i];  // store fractions
        dt *= 1.2;  // increase time-step
        krome_test::check(solve(x, Tgas, p, dt, o), "solve");  // call KROME
        t += dt;    // increase time
        krome_test::write_row(f66, {t}, x);  // dump evolution
        double d2 = 0.0;
        for (int i = 0; i < NSPEC; ++i) d2 += (xold[i] - x[i]) * (xold[i] - x[i]);
        if (d2 / NSPEC < 1e-20) break;  // check the steady-state
      }
      std::fprintf(f66, "\n\n");  // blank lines
    }
  }
  std::fclose(f66);
  std::fclose(f77);

  // print final abundances to screen
  std::printf("Final equilibrium abundances\n");
  for (int i = 0; i < NSPEC; ++i)
    std::printf("%-4s%11.3E\n", species_names[i], x[i]);

  std::printf("done!\nOutput written in abundances.dat\n"
              "To plot: python3 plot.py\n"
              "that's all! have a nice day!\n");
  return 0;
}
