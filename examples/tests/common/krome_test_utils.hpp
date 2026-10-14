// helpers of the test drivers (python3 -m kromex --test), on the KROME-X C++ API; the
// legacy krome_user conveniences (names_header, print_best_flux,
// explore_flux, print_info, scale_Z) are part of the API. Include after
// krome.hpp.
#pragma once
#include <cstdio>
#include <cstdlib>
#include <initializer_list>

namespace krome_test {
using namespace krome;

// stop with the solver's message on a failed solve
inline void check(int flag, const char* where) {
  if (flag == 0) return;
  std::printf("KROME failed (%s): flag %d %s\n", where, flag,
              flag_message(flag));
  std::exit(1);
}

// one output row: the leading values, then the species (divided by norm)
inline void write_row(std::FILE* f, std::initializer_list<double> lead,
                      const double* n, double norm = 1.0) {
  for (double v : lead) std::fprintf(f, "%17.8E", v);
  for (int i = 0; i < NSPEC; ++i) std::fprintf(f, "%17.8E", n[i] / norm);
  std::fprintf(f, "\n");
}

}  // namespace krome_test
