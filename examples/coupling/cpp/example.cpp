// KROME-X from C++ (the header-only API): the cell of the quick start
// (README), primordial gas at nH = 1e4 cm^-3 and 100 K, followed in time
// to 1e6 yr; the C, Fortran and Python examples solve the same cell
// (README.md)
#include <cstdio>
#include "krome.hpp"
using namespace krome;

int main() {
  const double spy = 3.15e7, nH = 1e4;
  SolverOpts o;
  o.evolve_T = true;
  Params p;
  p.tff = 1e14;  // the free-fall time of the COMPRESS heating
  double x[NSPEC] = {0.0};
  x[IDX_H] = nH;
  x[IDX_Hj] = 1.0;
  x[IDX_E] = 1.0;
  x[IDX_He] = 775.0;
  double Tgas = 100.0;

  std::FILE* f = std::fopen("abundances.dat", "w");
  std::fprintf(f, "#t_yr Tgas %s\n", names_header().c_str());
  double t = 0.0, dt = spy;
  while (t < 1e6 * spy) {
    const int flag = solve(x, Tgas, p, dt, o);
    if (flag != 0) {
      std::printf("KROME failed: %s\n", flag_message(flag));
      return 1;
    }
    t += dt;
    dt *= 1.1;
    std::fprintf(f, "%17.8E%17.8E", t / spy, Tgas);
    for (int i = 0; i < NSPEC; ++i) std::fprintf(f, "%17.8E", x[i] / nH);
    std::fprintf(f, "\n");
  }
  std::fclose(f);
  std::printf("wrote abundances.dat\n");
  return 0;
}
