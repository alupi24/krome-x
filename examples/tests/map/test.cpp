// This test represents a series of one-zone models with different photon
// flux (J21) and an initial temperature of 10^4 Kelvin.
// (KROME-X port of legacy tests/map)
#include <cmath>
#include <cstdio>
#include <vector>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const int itmax = 50;  // number of J21 values
  const int irmax = 50;  // number of rho values
  const double pi = krome::constants::pi;
  const double gravity = krome::constants::gravity;

  // photon bins (legacy -photoBins=10): 10 log-spaced bins between 5 and
  // 50 eV, and the network's photo cross sections on them
  const int B = 10;
  Radiation rad = Radiation::log_bins(B, 5.0, 5e1);
  PhotoModule photo;
  photo.init(rad);

  // solver options: thermal evolution on, the photon bins
  SolverOpts o;
  o.evolve_T = true;
  o.photo = &photo;
  Params p;  // runtime parameters: defaults, J21 set below

  const double J21min = std::log10(1e-3);   // min J21 value (log)
  const double J21max = std::log10(1e1);    // max J21 value (log)
  const double rhomin = std::log10(1e-25);  // min rho value (log)
  const double rhomax = std::log10(1e-21);  // max rho value (log)

  // initial mass fractions (legacy -useX)
  double xo[NSPEC];
  for (double& v : xo) v = 1e-20;  // default value
  xo[IDX_H] = 0.9225;   // H
  xo[IDX_Hj] = 1e-4;    // H+
  xo[IDX_He] = 0.0972;  // He
  xo[IDX_H2] = 1e-5;    // H2
  xo[IDX_HD] = 1e-8;    // HD

  // computes electrons (as legacy, from the mass fractions)
  xo[IDX_E] = get_electrons_neutral(xo);

  const double Tgaso = 1e4;  // initial temperature
  double sum = 0.0;
  for (double v : xo) sum += v;
  for (double& v : xo) v /= sum;  // normalize

  // the densities of one J21: one cell each, with its own photon field
  std::vector<double> x(irmax * NSPEC), rho(irmax), Tgas(irmax), dt(irmax),
      phi(irmax * B);
  CellsOpts co;
  co.Phi = phi.data();

  std::FILE* f = std::fopen("map.dat", "w");
  std::fprintf(f, "#rho J21 TgasInit Tgas %s\n",
               names_header().c_str());
  std::printf("%12s\n", "J21");
  // loop on J21s
  for (int it = 1; it <= itmax; ++it) {
    const double J21 =
        std::pow(1e1, (it - 1) * (J21max - J21min) / (itmax - 1) + J21min);

    // the H2 photodissociation rate 9e-13 J21 (legacy user_J21), and the
    // photon field: 1/E profile scaled by J21, the same in every cell
    p.J21 = J21;
    spectra::j21(B, rad.Elo.data(), rad.Ehi.data(), J21, 1.0, phi.data());
    for (int ir = 1; ir < irmax; ++ir)
      for (int b = 0; b < B; ++b) phi[ir * B + b] = phi[b];

    std::printf("%12.3E%12.2f%%\n", J21, it * 1e2 / itmax);
    // loop on densities
    for (int ir = 0; ir < irmax; ++ir) {
      for (int i = 0; i < NSPEC; ++i) x[ir * NSPEC + i] = xo[i];  // initial species
      rho[ir] = std::pow(1e1, ir * (rhomax - rhomin) / (irmax - 1) + rhomin);  // g/cm3
      dt[ir] = std::sqrt(3.0 * pi / 32.0 / gravity / rho[ir]);  // free-fall time
      Tgas[ir] = Tgaso;  // copy initial Tgas
    }
    krome_test::check(solve_cells_x(irmax, x.data(), rho.data(), Tgas.data(),
                                    dt.data(), p, o, co), "solve");  // KROME
    for (int ir = 0; ir < irmax; ++ir)
      krome_test::write_row(f, {rho[ir], J21, Tgaso, Tgas[ir]},
                            &x[ir * NSPEC]);  // dump
    std::fprintf(f, "\n");
  }
  std::fclose(f);
  std::printf("Models done: %d\nEnd of map test!\n"
              "To plot the map: python3 plot.py\n",
              itmax * irmax);
  return 0;
}
