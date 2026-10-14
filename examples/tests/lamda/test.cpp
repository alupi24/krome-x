// This test uses the data from the LAMDA database to compute the cooling
// from the O2 molecule. It consists in a set of runs with different
// densities and temperatures. The chemistry is solved to the equilibrium,
// then the cooling is computed. The cooling function is not a table, but
// is computed on the fly using the atomic data provided by LAMDA.
// http://home.strw.leidenuniv.nl/~moldata/ (KROME-X port of legacy tests/lamda)
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const double spy = constants::seconds_per_year;
  const int imax = 10, jmax = 10;

  // solver options: defaults, Tgas constant (legacy krome_thermo_OFF)
  SolverOpts o;
  o.evolve_T = false;

  // the O2 cooling term of get_thermal (legacy krome_coolingO2)
  int ncool = 0, iO2 = -1;
  for (; cooling_name(ncool) != nullptr; ++ncool)
    if (std::strcmp(cooling_name(ncool), "O2") == 0) iO2 = ncool;
  if (iO2 < 0) {
    std::printf("no O2 cooling term\n");
    return 1;
  }
  std::vector<double> cools(ncool);

  std::FILE* f44 = std::fopen("cooling.dat", "w");
  std::FILE* f55 = std::fopen("equilibrium.dat", "w");
  std::FILE* f66 = std::fopen("abundances.dat", "w");
  std::fprintf(f44, "#ntot Tgas cooling\n");
  std::fprintf(f55, "#Tgas %s\n", names_header().c_str());
  std::fprintf(f66, "#time Tgas %s\n", names_header().c_str());

  std::printf("%8s%12s\n", "%", "ntot");
  for (int j = 1; j <= jmax; ++j) {
    const double ntot =
        std::pow(10.0, (j - 1) * (9.0 - 3.0) / (jmax - 1) + 3.0);
    std::printf("%8.1f%12.3E\n", j * 1e2 / jmax, ntot);
    for (int i = 1; i <= imax; ++i) {
      double Tgas = std::pow(10.0, (i - 1) * (4.0 - 1.0) / (imax - 1) + 1.0);

      double x[NSPEC];
      for (double& v : x) v = 1e-40;

      x[IDX_H] = ntot;             // H
      x[IDX_H2] = 1e-6 * ntot;     // H2
      x[IDX_Hj] = 1e-4 * ntot;     // H+
      x[IDX_He] = 0.0775 * ntot;   // He

      Params p;
      // set cosmic rays (legacy krome_set_user_crate)
      p.crate = 1.3e-17;

      // rescale metallicity for neutral metals (C,Si,O)
      scale_Z(x, 1e-1);

      // legacy krome_set_user_Av
      p.Av = std::pow(ntot * 1e-3, 2.0 / 3.0);

      x[IDX_Cj] = x[IDX_C];  // carbon is fully ionized
      x[IDX_C] = 1e-40;

      x[IDX_E] = get_electrons_neutral(x);

      double t = 0.0, dt = spy;
      for (;;) {
        dt *= 1.1;
        krome_test::check(solve(x, Tgas, p, dt, o), "solve");
        t += dt;
        krome_test::write_row(f66, {t / spy, Tgas}, x);
        if (t > 1e8 * spy) break;
      }
      krome_test::check(get_thermal(x, Tgas, p, o, nullptr, nullptr,
                                    nullptr, cools.data()),
                        "get_thermal");
      std::fprintf(f44, "%17.8E%17.8E%17.8E\n", ntot, Tgas,
                   cools[iO2] / ntot);
      krome_test::write_row(f55, {Tgas}, x);
    }
    std::fprintf(f44, "\n");
  }
  std::fclose(f44);
  std::fclose(f55);
  std::fclose(f66);

  std::printf("Test OK!\nTo plot: python3 plot.py\n");
  return 0;
}
