// This is a test file to understand how automatic reactions and
// bin-based photochemistry work. This is a simple network with only
//   C+ + e -> C
//   C + e -> C+ + e + e
//   C + photon -> C+ + e
// the evolution is computed for different j21 values with constant Tgas
// and density. See networks/auto.knet. (KROME-X port of legacy tests/auto)
#include <cstdio>
#include <vector>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const double spy = 365.0 * 24.0 * 3600.0;  // seconds per year
  const double j21s[3] = {1e-3, 1e-1, 1e0};  // list of j21 values

  // photon bins (legacy -photoBins=10): 10 log-spaced bins between 1e1
  // and 1e2 eV, and the network's photo cross sections on them
  const int B = 10;
  Radiation rad = Radiation::log_bins(B, 1e1, 1e2);
  PhotoModule photo;
  photo.init(rad);
  std::vector<double> phi(B);

  Params p;      // runtime parameters: defaults
  SolverOpts o;  // solver options: defaults (isothermal), the photon field
  o.photo = &photo;
  o.Phi = phi.data();

  // header
  std::FILE* f = std::fopen("abundances.dat", "w");
  std::fprintf(f, "#j21 time %s\n", names_header().c_str());
  // loop on j21 values
  for (const double j21 : j21s) {
    // radiation field between 1e1 and 1e2 eV, 1/E profile (J21 power
    // law, alpha = 1) scaled according to j21
    spectra::j21(B, rad.Elo.data(), rad.Ehi.data(), j21, 1.0, phi.data());
    std::printf("running with j21= %g\n", j21);
    // init abundances
    double x[NSPEC];
    for (double& v : x) v = 1e-40;  // default abundances
    x[IDX_C] = .5;
    x[IDX_Cj] = .5;
    x[IDX_E] = get_electrons_neutral(x);

    double Tgas = 1.2e4;  // gas temperature (K)
    double dt = spy;      // initial time-step (s)
    double t = 0.0;       // start time (s)

    for (;;) {
      dt *= 1.1;  // increase time-step
      t += dt;    // advance time
      krome_test::check(solve(x, Tgas, p, dt, o), "solve");  // call KROME
      krome_test::write_row(f, {j21, t / spy}, x);           // dump
      if (t > 1e8 * spy) break;  // exit when 1e8 years is reached
    }
    std::fprintf(f, "\n");
  }
  std::fclose(f);

  std::printf("done!\nTo plot: python3 plot.py\n");
  return 0;
}
