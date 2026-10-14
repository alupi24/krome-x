// The Fortran-interface test, C++ side: the collapse of test.f90 through
// the C++ API (same loop, parameters and output columns, into
// abundances_cpp.dat)
#include <cstdio>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  SolverOpts o;
  o.evolve_T = true;
  Params p;
  p.Tcmb = 2.73;  // the legacy default (krome_redshift was unused)
  const double ntot = 1.0;
  double Tgas = 1e2, x[NSPEC];
  for (double& v : x) v = 1e-40;
  x[IDX_H] = ntot;
  x[IDX_H2] = 1e-6 * ntot;
  x[IDX_E] = 1e-4 * ntot;
  x[IDX_Hj] = 1e-4 * ntot;
  x[IDX_He] = 0.0775 * ntot;

  std::FILE* f = std::fopen("abundances_cpp.dat", "w");
  std::fprintf(f, "#ntot Tgas %s\n", names_header().c_str());
  double dd = ntot;
  for (int i = 1; i <= 500000; ++i) {
    const double dd1 = dd;
    const double tff = get_free_fall_time(x);
    p.tff = tff;
    const double dtH = 0.01 * tff;
    dd += dd / tff * dtH;
    for (double& v : x) v *= dd / dd1;
    if (dd > 1e18) break;
    krome_test::check(solve(x, Tgas, p, dtH, o), "solve");
    krome_test::write_row(f, {dd, Tgas}, x, dd);
  }
  std::fclose(f);
  std::printf("wrote abundances_cpp.dat; now: python3 check_diffs.py\n");
  return 0;
}
