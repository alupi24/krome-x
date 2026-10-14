// This is a simple one-zone collapse test following the chemical and
// thermal evolution of a metal-enriched cloud. The dynamics is described
// by the Larson-Penston-type similar solution and includes cooling and
// heating processes. For additional details look also to Omukai 2005.
// (KROME-X port of legacy tests/collapseCO: solar metallicity, CO and H2O
// freeze-out and evaporation on dust, tabulated dust temperature,
// dust cooling and H2 formation on dust)
#include <cmath>
#include <cstdio>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const int rstep = 500000;

  // solver options: thermal evolution on (near n = 1e12 cm^-3 the gas
  // drops from ~1e9 K to ~12 K within one step: the restarts shrink the
  // first step)
  SolverOpts o;
  o.evolve_T = true;

  // runtime parameters: the legacy defaults the test relied on (CMB
  // temperature 2.73 K, the Draine field G0 = 1.69 Habing of the legacy
  // PHOTODUST heating), the cosmic rays (legacy krome_set_user_crate);
  // Av and tff are set at every step
  Params p;
  p.Tcmb = 2.73;
  p.G0 = 1.69;
  p.crate = 1.3e-17;

  // INITIAL CONDITIONS
  const double ntot = 0.1;  // total density in 1/cm3
  double Tgas = 3e2;        // temperature in kelvin

  // species initialization in 1/cm3
  double x[NSPEC];
  for (double& v : x) v = 1e-40;

  x[IDX_H] = ntot;            // H
  x[IDX_H2] = 1e-6 * ntot;    // H2
  x[IDX_Hj] = 1e-4 * ntot;    // H+
  x[IDX_He] = 0.0775 * ntot;  // He

  // rescale metallicity for neutral metals (C,Fe,Si,O)
  scale_Z(x, 0.0);

  x[IDX_Cj] = x[IDX_C];  // carbon is fully ionized
  x[IDX_C] = 1e-40;

  x[IDX_E] = get_electrons_neutral(x);

  double dd = ntot;
  std::printf("%5s%5s%11s%11s\n", "Zstep", "step", "n(cm-3)", "Tgas(K)");

  // output header (legacy abundances.dat had none)
  std::FILE* f22 = std::fopen("abundances.dat", "w");
  std::fprintf(f22, "#ntot Tgas Tdust %s\n",
               names_header().c_str());

  // loop over the hydro time-step
  for (int i = 1; i <= rstep; ++i) {
    const double dd1 = dd;

    // free-fall time (legacy user_tff)
    const double tff = get_free_fall_time(x);
    p.tff = tff;
    const double dtH = 0.01 * tff;  // TIME-STEP
    dd += dd / tff * dtH;           // UPDATE DENSITY

    for (double& v : x) v *= dd / dd1;  // rescale abundances

    x[IDX_E] = get_electrons_neutral(x);

    if (dd > 1e16) break;

    // local approximation for Av (legacy krome_set_user_Av)
    p.Av = std::pow(dd * 1e-3, 2.0 / 3.0);

    // solve the chemistry
    krome_test::check(solve(x, Tgas, p, dtH, o), "solve");

    // the tabulated dust temperature (legacy krome_get_table_Tdust)
    double Td = 0.0;
    krome_test::check(get_Tdust(x, Tgas, p, o, &Td), "get_Tdust");
    krome_test::write_row(f22, {dd, Tgas, Td}, x, dd);

    // print every 100 steps
    if (i % 100 == 0) std::printf("%5d%11.3E%11.3E\n", i, dd, Tgas);
  }
  std::fclose(f22);

  std::printf("To plot: python3 plot.py\nThat's all! have a nice day!\n");
  return 0;
}
