// This test is aimed to reproduce the evolution of the primordial gas
// based on the Galli&Palla 1998 paper (KROME-X port of legacy
// tests/earlyUniverse)
#include <cmath>
#include <cstdio>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const int rstep = 5000;
  // the legacy KROME cosmology: h, Omega0, OmegaB and Hubble0 = 100 h
  // km/s/Mpc [1/s] (legacy took 1 Mpc = 3.08e24 cm)
  const double hubble = 0.704, Omega0 = 1.0, OmegaB = 0.0456;
  const double Hubble0 = 1e2 * hubble * 1e5 / krome::constants::Mpc_to_cm;

  // solver options: thermal evolution on
  SolverOpts o;
  o.evolve_T = true;
  // electrons from charge neutrality after each solve (legacy
  // -computeElectrons): the network's freeze stops dn(e-)/dt at
  // z > 2243, not He's recombination
  o.conserve = 3;

  // runtime parameters: the cosmology of the EXPANSION cooling and of the
  // RECFAST recombination (REQUIRED), and the legacy temperature floor
  // phys_Tcmb = 2.73 K (legacy krome_set_zredshift did not change it);
  // zredshift and Trad are set at every step
  Params p;
  p.hubble0 = Hubble0;
  p.omega0 = Omega0;
  p.Tcmb = 2.73;

  // INITIAL CONDITIONS
  const double Tcmb0 = 2.73;              // K
  const double rhoc = 9.19e-6;            // h2/cm3
  const double zinit = 1e4;               // initial redshift
  double aexp = 1.0 + zinit;              // 1 + z factor
  const double ntot = OmegaB * rhoc;      // co-moving number density, 1/cm3
  double dd = ntot * aexp * aexp * aexp;  // proper initial number density
  double Trad = Tcmb0 * aexp;             // initial radiation temperature, K
  double Tgas = Trad;                     // initial gas temperature, K

  // INTEGRATION PARAMETERS
  const double zend = 1e1;                           // final redshift
  const double dzeta = (zend - zinit) / (rstep - 1.0);  // redshift step

  // species initialization in 1/cm3
  double x[NSPEC];
  for (double& v : x) v = 1e-40;

  // INIT: FULLY IONIZED GAS IN NUMBER DENSITY
  x[IDX_Hj] = 0.924 * dd;   // H+
  x[IDX_Hejj] = 0.076 * dd;  // He++
  x[IDX_Dj] = 4.3e-5 * dd;  // D+
  x[IDX_E] = get_electrons_neutral(x);

  // the initial state (legacy krome_get_info)
  std::printf("%4s %-9s%11s%5s%11s\n", "#", "Name", "m (g)", "Chrg", "x");
  double sum = 0.0;
  for (int i = 0; i < NSPEC; ++i) {
    std::printf("%4d %-9s%11.3E%5d%11.3E\n", i + 1, species_names[i],
                species_mass_g[i], (int)species_charge[i], x[i]);
    sum += x[i];
  }
  std::printf("%30s%11.3E\n%14s%11.3E\n", "sum", sum, "Tgas", Tgas);

  std::printf("solving...\n%7s%11s%11s%11s%11s\n", "step", "z", "n(cm-3)",
              "Tgas(K)", "Trad(K)");

  // header
  std::FILE* f = std::fopen("abundances.dat", "w");
  std::fprintf(f, "#z density Tgas Trad %s\n",
               names_header().c_str());

  // loop over the hydro time-step
  for (int i = 1; i <= rstep; ++i) {
    const double dd1 = dd;
    const double redshift = zinit + dzeta * (i - 1);

    // H and D are computed using the approach of Peebles 1968 book, see
    // eqn.(6.95). note also that for redshift higher than 2.243d3 the
    // reactions involving everything except He are set to zero, by the
    // freeze directive of the network earlyUniverse.knet
    if (redshift > 2.243e3) {
      const double u = std::pow(
          10.0, 20.99 - std::log10(OmegaB * hubble * hubble
                                   * std::pow(1.0 + redshift, 1.5))
                    - 2.5050e4 / (1.0 + redshift));
      x[IDX_H] = 0.924 * dd / (2.0 + u);
      x[IDX_D] = 4.3e-5 * dd / (2.0 + u);
    }

    // set the redshift (legacy krome_set_zredshift)
    p.zredshift = redshift;

    aexp = 1.0 + redshift;  // update 1 + z

    // convert redshift step in timestep, eq. (5) Galli&Palla 1998
    const double dtH =
        -dzeta / (Hubble0 * aexp * aexp * std::sqrt(Omega0 * redshift + 1.0));

    // radiation temperature (legacy krome_set_user_Trad)
    Trad = Tcmb0 * aexp;
    p.Trad = Trad;

    // the new density, species update
    dd = ntot * aexp * aexp * aexp;
    for (double& v : x) v *= dd / dd1;

    // solve the chemistry
    krome_test::check(solve(x, Tgas, p, dtH, o), "solve");

    krome_test::write_row(f, {aexp, dd, Tgas, Trad}, x, dd);
    if (i % 100 == 0)
      std::printf("%7d%11.3E%11.3E%11.3E%11.3E\n", i, aexp, dd, Tgas, Trad);
  }
  std::fclose(f);

  std::printf("To plot: python3 plot.py\nThat's all! have a nice day!\n");
  return 0;
}
