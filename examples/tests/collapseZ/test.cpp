// This is a simple one-zone collapse test following the chemical and
// thermal evolution of a metal-enriched cloud. The dynamics is described
// by the Larson-Penston-type similar solution and includes cooling and
// heating processes. For additional details look also to Omukai 2005.
// (KROME-X port of legacy tests/collapseZ)
#include <cstdio>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const int rstep = 500000;
  const double zs[] = {-99.0, -4.0, -3.0, -2.0, -1.0};  // log Z/Zsun

  // solver options: thermal evolution on
  SolverOpts o;
  o.evolve_T = true;

  // runtime parameters: legacy krome_set_zredshift(15) set only the
  // redshift (Compton cooling); the legacy CMB temperature, the lower
  // temperature cap of the rates and cooling, stayed at its default
  // 2.73 K. tff is set at every step
  Params p;
  p.zredshift = 15.0;
  p.Tcmb = 2.73;

  // output header (legacy wrote the stored results at the end)
  std::FILE* f22 = std::fopen("abundances.dat", "w");
  std::fprintf(f22, "#Z rho Tgas %s\n", names_header().c_str());

  for (int jz = 0; jz < 5; ++jz) {
    // INITIAL CONDITIONS
    const double ntot = 0.1;  // total density in 1/cm3
    double Tgas = 3e2;        // temperature in kelvin

    // species initialization in 1/cm3
    double x[NSPEC];
    for (double& v : x) v = 1e-40;

    x[IDX_H] = ntot;             // H
    x[IDX_H2] = 1e-6 * ntot;     // H2
    x[IDX_E] = 1e-4 * ntot;      // E
    x[IDX_Hj] = 1e-4 * ntot;     // H+
    x[IDX_He] = 0.0775 * ntot;   // He

    // rescale metallicity for neutral metals (C,Fe,Si,O)
    scale_Z(x, zs[jz]);

    x[IDX_Cj] = x[IDX_C];  // carbon is fully ionized
    x[IDX_C] = 1e-40;

    // list abundances
    print_info(x, Tgas);

    double dd = ntot;

    std::printf(" solving...\n%5s%5s%11s%11s\n", "z", "step", "n(cm-3)",
                "Tgas(K)");

    // loop over the hydro time-step
    for (int i = 1; i <= rstep; ++i) {
      const double dd1 = dd;

      // free-fall time (legacy user_tff)
      const double tff = get_free_fall_time(x);
      p.tff = tff;
      const double dtH = 0.01 * tff;  // TIME-STEP
      dd += dd / tff * dtH;           // UPDATE DENSITY

      for (double& v : x) v *= dd / dd1;  // rescale abundances

      if (dd > 1e18) break;  // quit after 1e18 1/cm3

      // solve the chemistry
      krome_test::check(solve(x, Tgas, p, dtH, o), "solve");

      // dump Tgas and normalized abundances
      krome_test::write_row(f22, {zs[jz], dd, Tgas}, x, dd);
      if (i % 100 == 0)  // print every 100 steps
        std::printf("%5d%5d%11.3E%11.3E\n", jz + 1, i, dd, Tgas);
    }
    std::fprintf(f22, "\n");  // blank line between the metallicities
  }
  std::fclose(f22);

  std::printf("To plot: python3 plot.py\nThat's all! have a nice day!\n");
  return 0;
}
