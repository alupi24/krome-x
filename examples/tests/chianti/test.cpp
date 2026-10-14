// This is a test to plot cooling function using e.g. the information from
// the Chianti database. http://www.chiantidatabase.org
// It computes a cooling curve similar to the one shown in Gnat+Ferland
// 2012 ApJS 199 20 (Fig.2)
// http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1111.6980
// (KROME-X port of legacy tests/chianti)
#include <cmath>
#include <cstdio>
#include <vector>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const double spy = constants::seconds_per_year;
  const double j21s[3] = {0.0, 1.0, 1e2};  // list of j21 values
  const int kmax = 30, imaxx = 5000;

  // photon bins (legacy -photoBins=10): 10 log-spaced bins between 1e1
  // and 1e2 eV, and the network's photo cross sections on them
  const int B = 10;
  Radiation rad = Radiation::log_bins(B, 1e1, 1e2);
  PhotoModule photo;
  photo.init(rad);
  std::vector<double> phi(B);

  Params p;      // runtime parameters: defaults
  SolverOpts o;  // defaults: Tgas constant (legacy krome_thermo_off)
  o.evolve_T = false;
  o.photo = &photo;
  o.Phi = phi.data();

  std::FILE* f55 = std::fopen("cooling.dat", "w");
  std::FILE* f66 = std::fopen("abundances.dat", "w");
  std::FILE* f77 = std::fopen("equilibrium.dat", "w");
  std::fprintf(f66, "#j time Tgas %s\n", names_header().c_str());
  std::fprintf(f77, "#j Tgas %s\n", names_header().c_str());

  // loop on j21 values
  for (int j = 1; j <= 3; ++j) {
    // radiation field between 1e1 and 1e2 eV, 1/E profile (legacy
    // krome_set_photoBin_J21log), scaled according to j21 (legacy
    // krome_photoBin_scale)
    const double j21 = j21s[j - 1];
    spectra::j21(B, rad.Elo.data(), rad.Ehi.data(), j21, 1.0, phi.data());
    std::printf(" running with j21= %g\n", j21);

    std::fprintf(f55, "#J21 Tgas cooling\n");
    // loop on temperature values
    for (int k = 1; k <= kmax; ++k) {
      // init abundances
      double x[NSPEC];
      for (double& v : x) v = 0.0;  // default abundances
      x[IDX_Cj] = 1.0;
      x[IDX_E] = 1.0;

      // gas temperature (K)
      double Tgas = std::pow(10.0, (k - 1) * (8.0 - 4.0) / (kmax - 1) + 4.0);
      if (k % 10 == 0) std::printf("%7.2f %%\n", (k - 1) * 1e2 / (kmax - 1));

      double dt = spy;  // initial time-step (s)
      double t = 0.0;   // start time (s)

      // equilibrium
      for (int i = 1;; ++i) {
        dt *= 1.1;  // increase time-step
        t += dt;    // advance time
        // electron conservation
        x[IDX_E] = get_electrons_neutral(x);
        krome_test::check(solve(x, Tgas, p, dt, o), "solve");
        if (i <= imaxx) {
          std::fprintf(f66, "%5d", j);
          krome_test::write_row(f66, {t / spy, Tgas}, x);
        }
        if (t > 1e8 * spy) break;  // exit when 1e8 years reached
      }
      std::fprintf(f66, "\n");

      // the abundances of the species and the corresponding cooling
      double cool = 0.0;
      krome_test::check(get_thermal(x, Tgas, p, o, nullptr, &cool),
                        "get_thermal");
      std::fprintf(f77, "%5d", j);
      krome_test::write_row(f77, {Tgas}, x);
      std::fprintf(f55, "%17.8E%17.8E%17.8E\n", j21, Tgas, cool);
    }
    std::fprintf(f55, "\n");
    std::fprintf(f77, "\n");
    std::fprintf(f66, "\n");
  }
  std::fclose(f55);
  std::fclose(f66);
  std::fclose(f77);

  std::printf("done!\nTo plot: python3 plot.py\n");
  return 0;
}
