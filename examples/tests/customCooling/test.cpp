// this test shows how to use custom cooling functions
// check the network networks/customCool.knet for details
// (KROME-X port of legacy tests/customCooling)
#include <cmath>
#include <cstdio>
#include <vector>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const int jmax = 30;
  Params p;      // runtime parameters: defaults
  SolverOpts o;  // options of the queries: defaults

  // the configured cooling terms (legacy krome_get_cooling_array had
  // every legacy slot, CUSTOM among them; KROME-X only the configured ones,
  // here the network's custom cooling)
  int ncool = 0;
  while (cooling_name(ncool) != nullptr) ++ncool;
  std::vector<double> cools(ncool > 0 ? ncool : 1);

  std::FILE* f = std::fopen("cooling.dat", "w");
  std::fprintf(f, "#Tgas");
  for (int i = 0; i < ncool; ++i) std::fprintf(f, " %s", cooling_name(i));
  std::fprintf(f, "\n");

  // loop on Tgas points
  for (int j = 1; j <= jmax; ++j) {
    // compute Tgas
    const double Tgas =
        std::pow(10.0, (j - 1) * (5.0 - 1.0) / (jmax - 1) + 1.0);

    double x[NSPEC];
    for (double& v : x) v = 0.0;
    x[IDX_H2] = 1.0;
    x[IDX_H2O] = 1.0;
    x[IDX_H] = 1.0;
    double heat, cool;
    krome_test::check(get_thermal(x, Tgas, p, o, &heat, &cool, nullptr,
                                  cools.data()), "get_thermal");
    std::fprintf(f, "%17.8E", Tgas);
    for (int i = 0; i < ncool; ++i) std::fprintf(f, "%17.8E", cools[i]);
    std::fprintf(f, "\n");
  }
  std::fclose(f);

  std::printf("test DONE!\nTo plot: python3 plot.py\n");
  return 0;
}
