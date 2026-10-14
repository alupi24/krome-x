// KROME-X coupling example: a GIZMO/GADGET-like C++ particle code, on
// the C++ API (header-only krome.hpp + the project's net_<name>_data.o).
//
// The gas particles are an array of structs (position, mass, density,
// specific internal energy, species mass fractions, and their own KROME
// parameters). The active particles are gathered into KROME's layout
// n[npart * NSPEC], integrated with ONE solve_cells call with PER-PARTICLE
// Params and dt (the per-cell path runs the particles in parallel with
// OpenMP when compiled with -fopenmp: make OPENMP=-fopenmp), and scattered
// back.
//
// Operator splitting per step (the same problem as the Fortran, C and
// Python examples):
//   1. "hydro": a uniform compression rho -> 1.02 rho, adiabatic heating
//      u -> u 1.02^(gamma-1)
//   2. chemistry + cooling: Tgas from u with KROME's mu and gamma,
//      solve_cells over dt, then u from the new Tgas with the new mu and
//      gamma
// A particle whose solve fails (co.flags[i] != 0, left untouched by KROME)
// is retried with two, then four substeps.
//
// Output: result.dat, per particle: index, initial nH [cm^-3], final Tgas
// [K], n(H2)/nH, n(e-)/nH.
#include <cmath>
#include <cstdio>
#include <vector>
#include "krome.hpp"
using namespace krome;

constexpr double KB = 1.380649e-16, MP = 1.67262192e-24;
constexpr double DT = 3e11, FCOMP = 1.02;
constexpr int NPART = 16, NSTEP = 20;

struct Particle {
  double pos[3];
  double mass;       // [g]
  double rho;        // [g cm^-3]
  double u;          // specific internal energy [erg g^-1]
  double X[NSPEC];   // species mass fractions
  Params par;        // this particle's KROME parameters
};

// number densities [cm^-3] of a particle
static void densities(const Particle& q, double* n) {
  for (int s = 0; s < NSPEC; ++s)
    n[s] = q.rho * q.X[s] / species_mass_g[s];
}

// Tgas [K] from u; gamma may depend on Tgas (codegen thermal.gamma): two
// fixed-point passes from 1e4 K
static double temperature(const Particle& q, const double* n) {
  double T = 1e4;
  for (int pass = 0; pass < 2; ++pass)
    T = q.u * (get_gamma(n, T, q.par) - 1.0) * get_mu(n) * MP / KB;
  return T;
}

static double energy(const Particle& q, const double* n, double T) {
  return KB * T / ((get_gamma(n, T, q.par) - 1.0) * get_mu(n) * MP);
}

// a failed particle: 2, then 4 substeps of dt (each from the last good
// state)
static int subcycle(double* n, double& T, const Params& p,
                    const SolverOpts& o) {
  for (int nsub = 2; nsub <= 4; nsub += 2) {
    double ns[NSPEC], Ts = T;
    int ierr = 0;
    for (int s = 0; s < NSPEC; ++s) ns[s] = n[s];
    for (int j = 0; j < nsub && !ierr; ++j)
      ierr = solve(ns, Ts, p, DT / nsub, o);
    if (!ierr) {
      for (int s = 0; s < NSPEC; ++s) n[s] = ns[s];
      T = Ts;
      return 0;
    }
  }
  return -1;
}

int main() {
  SolverOpts o;  // solver options, once
  o.evolve_T = true;

  // initial state: nH from 1 to 1e4 cm^-3, 1e4 K, 1% ionized
  std::vector<Particle> part(NPART);
  std::vector<double> nH0(NPART);
  for (int i = 0; i < NPART; ++i) {
    Particle& q = part[i];
    nH0[i] = std::pow(10.0, 4.0 * i / (NPART - 1));
    double n[NSPEC];
    for (int s = 0; s < NSPEC; ++s) n[s] = 1e-20 * nH0[i];
    n[IDX_H] = 0.99 * nH0[i];
    n[IDX_Hj] = 0.01 * nH0[i];
    n[IDX_H2] = 1e-6 * nH0[i];
    n[IDX_He] = 0.0775 * nH0[i];
    n[IDX_E] = 0.01 * nH0[i];
    q.pos[0] = i;
    q.pos[1] = q.pos[2] = 0.0;
    q.rho = get_rho(n);
    q.mass = q.rho;  // a 1 cm^3 particle
    for (int s = 0; s < NSPEC; ++s) q.X[s] = species_mass_g[s] * n[s] / q.rho;
    // q.par: the network defaults; per-particle values go here (e.g. a
    // local radiation field, cosmic-ray rate or metallicity parameter)
    q.u = energy(q, n, 1e4);
  }

  std::vector<double> n(NPART * NSPEC), T(NPART), dt(NPART);
  std::vector<Params> par(NPART);
  std::vector<int> flags(NPART);
  for (int step = 0; step < NSTEP; ++step) {
    // 1. hydro: uniform compression
    for (Particle& q : part) {
      double nq[NSPEC];
      densities(q, nq);
      const double g = get_gamma(nq, temperature(q, nq), q.par);
      q.rho *= FCOMP;
      q.u *= std::pow(FCOMP, g - 1.0);
    }

    // 2. chemistry and cooling: gather, one call, scatter
    for (int i = 0; i < NPART; ++i) {
      densities(part[i], &n[i * NSPEC]);
      T[i] = temperature(part[i], &n[i * NSPEC]);
      dt[i] = DT;  // per particle (individual timesteps)
      par[i] = part[i].par;
    }
    CellsOpts co;
    co.flags = flags.data();
    if (solve_cells(NPART, n.data(), T.data(), dt.data(), par.data(), o, co))
      for (int i = 0; i < NPART; ++i)
        if (flags[i] && subcycle(&n[i * NSPEC], T[i], par[i], o))
          std::fprintf(stderr, "particle %d failed after subcycling: %s\n",
                       i, flag_message(flags[i]));
    for (int i = 0; i < NPART; ++i) {
      Particle& q = part[i];
      const double* ni = &n[i * NSPEC];
      for (int s = 0; s < NSPEC; ++s)
        q.X[s] = species_mass_g[s] * ni[s] / q.rho;
      q.u = energy(q, ni, T[i]);
    }
  }

  std::FILE* f = std::fopen("result.dat", "w");
  for (int i = 0; i < NPART; ++i) {
    double nq[NSPEC];
    densities(part[i], nq);
    const double nH = get_Hnuclei(nq);
    std::fprintf(f, "%4d %19.11E %19.11E %19.11E %19.11E\n", i, nH0[i],
                 temperature(part[i], nq), nq[IDX_H2] / nH, nq[IDX_E] / nH);
  }
  std::fclose(f);
  std::printf("wrote result.dat\n");
  return 0;
}
