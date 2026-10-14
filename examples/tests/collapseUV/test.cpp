// This is a simple one-zone collapse test following the chemical and
// thermal evolution of a primordial cloud. The dynamics is described by
// the Larson-Penston-type similar solution and includes cooling, heating
// processes and a UV background. For additional details look also to
// Omukai 2005. (KROME-X port of legacy tests/collapseUV)
#include <cstdio>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const int rstep = 500000;
  const double j21s[] = {0.0, 1.0, 1e2, 1e5};  // preset J21 values

  // solver options: thermal evolution on; legacy -conserve (element and
  // charge conservation) is the runtime option conserve = 2
  SolverOpts o;
  o.evolve_T = true;
  o.conserve = 2;

  // runtime parameters: legacy krome_set_zredshift(15) set only the
  // redshift (Compton cooling); the legacy CMB temperature (CIE cooling,
  // lower temperature cap) stayed at its default 2.73 K. tff is set at
  // every step
  Params p;
  p.zredshift = 15.0;
  p.Tcmb = 2.73;

  // output (legacy wrote the stored results at the end)
  std::FILE* f22 = std::fopen("abundances.dat", "w");
  std::fprintf(f22, "#J21 ntot Tgas H2 H\n");

  for (int j = 0; j < 4; ++j) {
    const double j21 = j21s[j];
    const double ntot = 1.0;  // total density in 1/cm3
    double Tgas = 1e2;        // temperature in kelvin

    p.J21 = j21;  // legacy krome_set_user_J21

    // species initialization in 1/cm3
    double x[NSPEC];
    for (double& v : x) v = 1e-40;

    x[IDX_H] = ntot;            // H
    x[IDX_H2] = 1e-6 * ntot;    // H2
    x[IDX_E] = 1e-4 * ntot;     // E
    x[IDX_Hj] = 1e-4 * ntot;    // H+
    x[IDX_He] = 0.0775 * ntot;  // He

    double dd = ntot;

    std::printf(" solving for J21= %g\n%5s%11s%11s\n", j21, "step",
                "n(cm-3)", "Tgas(K)");

    // loop over the hydro time-step
    for (int i = 1; i <= rstep; ++i) {
      const double dd1 = dd;

      // free-fall time (legacy user_tff)
      const double tff = get_free_fall_time(x);
      p.tff = tff;
      const double dtH = 0.01 * tff;  // TIME-STEP
      dd += dd / tff * dtH;           // UPDATE DENSITY

      for (double& v : x) v *= dd / dd1;

      if (dd > 1e16) break;

      // solve the chemistry
      krome_test::check(solve(x, Tgas, p, dtH, o), "solve");

      std::fprintf(f22, "%17.8E%17.8E%17.8E%17.8E%17.8E\n", j21, dd, Tgas,
                   x[IDX_H2] / dd, x[IDX_H] / dd);
      if (i % 100 == 0) std::printf("%5d%5d%11.3E%11.3E\n", j + 1, i, dd,
                                   Tgas);
    }
    std::fprintf(f22, "\n");
    std::printf("\n");
  }
  std::fclose(f22);

  std::printf("To plot: python3 plot.py\nThat's all! have a nice day!\n");
  return 0;
}
