// This is a simple one-zone collapse test following the chemical and
// thermal evolution of a primordial cloud. The dynamics is described by
// the Larson-Penston-type similar solution and includes cooling and
// heating processes. For additional details look also to Omukai 2000 and
// the KROME paper. (KROME-X port of legacy tests/collapse)
#include <cstdio>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const int rstep = 500000;
  const double ntot = 1.0;  // total density, cm-3
  double Tgas = 1e2;        // temperature, K

  // solver options: thermal evolution on
  SolverOpts o;
  o.evolve_T = true;

  // runtime parameters: the legacy default Tcmb = 2.73 K (the legacy
  // driver's krome_redshift = 15 is an unused user variable: legacy ran
  // at z = 0); tff is set at every step
  Params p;
  p.Tcmb = 2.73;

  // species default, cm-3
  double x[NSPEC];
  for (double& v : x) v = 1e-40;

  // set individual species
  x[IDX_H] = ntot;
  x[IDX_H2] = 1e-6 * ntot;
  x[IDX_E] = 1e-4 * ntot;
  x[IDX_Hj] = 1e-4 * ntot;
  x[IDX_He] = 0.0775 * ntot;

  double dd = ntot;  // initial density

  // explore data, and the output with its header
  std::FILE* fexp = std::fopen("explore.dat", "w");
  std::FILE* f22 = std::fopen("abundances.dat", "w");
  std::fprintf(f22, "#ntot Tgas %s\n", names_header().c_str());

  std::printf("solving...\n%5s%11s%11s\n", "step", "n(cm-3)", "Tgas(K)");

  // loop on density steps
  for (int i = 1; i <= rstep; ++i) {
    const double dd1 = dd;  // old density

    // free-fall time, s: the compression timescale of COMPRESS heating
    // (legacy user_tff)
    const double tff = get_free_fall_time(x);
    p.tff = tff;
    const double dtH = 0.01 * tff;  // time-step
    dd += dd / tff * dtH;           // density increase

    // rescale density
    for (double& v : x) v *= dd / dd1;

    if (dd > 1e18) break;  // max density reached

    // solve the chemistry
    krome_test::check(solve(x, Tgas, p, dtH, o), "solve");

    // output
    krome_test::write_row(f22, {dd, Tgas}, x, dd);
    if (i % 100 == 0) {
      std::printf("%5d%11.3E%11.3E\n", i, dd, Tgas);
      print_best_flux(x, Tgas, p, 5, o);
      explore_flux(fexp, x, Tgas, p, dd, o);
    }
  }
  std::fclose(fexp);
  std::fclose(f22);

  std::printf("To plot: python3 plot.py\nThat's all! have a nice day!\n");
  return 0;
}
