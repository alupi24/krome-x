// This test consists of a gas with an initial population of dust, both
// silicon and carbon based, with 10 size-bins each. The dust-to-gas mass
// ratio is 1e-2 and the dust distribution follows an MRN power-law. For
// each of 100 gas temperatures between 10 and 1e4 K (Tdust = Tgas) the
// gas and the grains evolve for 1000 yr: the key elements C and Si
// accrete onto the grains, which evaporate at high Tdust. (KROME-X port of
// legacy tests/dust)
#include <cmath>
#include <cstdio>
#include <vector>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const int imax = 100;
  const double spy = constants::seconds_per_year;
  const double xH = 1e6;  // hydrogen density

  // the grain population: 10 bins of silicates and 10 of carbon grains,
  // MRN 5e-7 - 2.5e-5 cm (legacy krome_init_dust_distribution defaults)
  DustBins dust;
  dust.nd = 10;
  dust.init();
  const int nd = dust.nb;

  // number densities (cm-3)
  double xi[NSPEC];
  for (double& v : xi) v = 1e-20;
  xi[IDX_H] = 1.0 * xH;
  xi[IDX_Hj] = 1e-4 * xH;
  xi[IDX_H2] = 1e-5 * xH;
  xi[IDX_C] = 1e-4 * xH;
  xi[IDX_Si] = 1e-4 * xH;
  xi[IDX_E] = get_electrons_neutral(xi);  // globally neutral

  Params p;
  p.d2g = 1e-2;  // dust-to-gas mass ratio of the population
  SolverOpts o;  // thermal evolution off: Tgas and Tdust held
  o.dust = &dust;
  const double dt = 1e3 * spy;  // time-step (s)

  std::FILE* f66 = std::fopen("gas.dat", "w");
  std::FILE* f88 = std::fopen("mass_check.dat", "w");
  std::fprintf(f66, "#Tgas mC mSi\n");
  std::fprintf(f88, "#Tgas mtot\n");
  std::vector<double> data(static_cast<size_t>(imax) * nd), dataT(imax);
  const double mC = species_mass_g[IDX_C], mSi = species_mass_g[IDX_Si];

  // loop on Tgas
  for (int i = 0; i < imax; ++i) {
    double x[NSPEC];
    for (int s = 0; s < NSPEC; ++s) x[s] = xi[s];
    // re-init the dust: the reference sizes, Tdust = Tgas
    double Tgas = std::pow(10.0, i * (4.0 - 1.0) / (imax - 1) + 1.0);
    for (int b = 0; b < nd; ++b) {
      dust.a[b] = dust.asize[b];
      dust.Td[b] = Tgas;
    }
    dust.set_radii(nullptr);
    if ((i + 1) % 10 == 0) std::printf("%g\n", Tgas);

    // the grain counts of this gas (fixed while the radii evolve)
    const double rho = get_rho(x);
    krome_test::check(solve(x, Tgas, p, dt, o), "solve");

    // mass density in the dust bins, with each type's bulk density
    // (legacy used 2.2 g cm^-3 for both)
    double mdust = 0.0;
    for (int b = 0; b < nd; ++b) {
      const double nb = dust.xw[b] * rho * p.d2g;
      data[static_cast<size_t>(i) * nd + b] = 4.0 * constants::pi / 3.0
          * std::pow(dust.a[b], 3) * (dust.isC[b] ? 2.25 : 3.13) * nb;
      mdust += data[static_cast<size_t>(i) * nd + b];
    }
    // mass density in the gas phase, C and Si
    std::fprintf(f66, "%17.8E%17.8E%17.8E\n", Tgas, x[IDX_C] * mC,
                 x[IDX_Si] * mSi);
    std::fprintf(f88, "%17.8E%17.8E\n", Tgas,
                 mdust + x[IDX_C] * mC + x[IDX_Si] * mSi);
    dataT[i] = Tgas;
  }
  std::fclose(f66);
  std::fclose(f88);

  // dust mass density evolution, bin by bin (legacy type index: 0
  // carbon, 1 silicon)
  std::FILE* f77 = std::fopen("dust.dat", "w");
  std::fprintf(f77, "#dustType Tgas rhodust\n");
  for (int b = 0; b < nd; ++b) {
    for (int i = 0; i < imax; ++i)
      std::fprintf(f77, "%4d%17.8E%17.8E\n", dust.isC[b] ? 0 : 1, dataT[i],
                   data[static_cast<size_t>(i) * nd + b]);
    std::fprintf(f77, "\n");
  }
  std::fclose(f77);

  std::printf("To plot: python3 plot.py\nThat's all! have a nice day!\n");
  return 0;
}
