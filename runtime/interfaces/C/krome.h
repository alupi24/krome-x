/* KROME-X C interface: the ABI of the libkrome library, for C and C++
 * callers and for the optional language bindings built on it (Fortran:
 * runtime/interfaces/fortran/krome.f90, codegen project.bindings = ["fortran"]; Python:
 * runtime/interfaces/python/krome.py, codegen project.bindings = ["python"]). Network-agnostic: the network-specific constants
 * (KROME_NSPEC, KROME_IDX_<species>, KROME_PAR_<param>) are in the
 * generated krome_net.h; the implementation (krome.cpp) is compiled
 * once per project against the generated network.
 *
 * Array layouts (C row-major == Fortran column-major, so no transposes):
 *   n[ncell][nspec]    number densities [cm^-3]        Fortran n(nspec,ncell)
 *   par[ncell][npar]   runtime parameters, param order  Fortran par(npar,ncell)
 *                      (par_stride = npar), or ONE row broadcast to every
 *                      cell (par_stride = 0); krome_params_default fills
 *                      a row with the network defaults
 *   Phi[ncell][nbins]  photon flux per bin [cm^-2 s^-1], used by the
 *                      photochemistry, photoheating and grains: held fixed
 *                      during the call, or co-evolved and updated in place
 *                      when the network was generated with radiation.evolve
 *
 * Species and parameter indices are 0-based here (the Fortran module
 * provides 1-based ones). Functions returning int return 0 on success.
 */
#ifndef KROME_H
#define KROME_H

/* return flags: 0 success; a negative CVODE flag (-1 too much work, -3
 * error-test failures, -4 convergence failures, ...) is a failed
 * integration, the cell left untouched; the KROME flags below reject the
 * inputs before integrating. Functions other than the solvers return -1
 * for bad arguments. krome_flag_message() describes any flag. */
/* physical constants (kromex/constants.py: CODATA 2018, IAU), the
 * numbers the solver uses (C++ krome::constants, Fortran
 * krome_constants) */
#define KROME_CONST_boltzmann_eV 8.617333262145179e-05  /* eV / K; exact */
#define KROME_CONST_boltzmann_J 1.380649e-23  /* J / K; exact */
#define KROME_CONST_boltzmann_erg 1.380649e-16  /* erg / K; exact */
#define KROME_CONST_planck_eV 4.135667696923859e-15  /* eV s; exact */
#define KROME_CONST_planck_J 6.62607015e-34  /* J s; exact */
#define KROME_CONST_planck_erg 6.62607015e-27  /* erg s; exact */
#define KROME_CONST_hbar_erg 1.0545718176461565e-27  /* erg s; exact */
#define KROME_CONST_hc_k 1.4387768775039336  /* cm K; second radiation constant: E/k [K] = hc_k x E [cm^-1] */
#define KROME_CONST_gravity 6.6743e-08  /* cm^3 g^-1 s^-2; CODATA 2018 */
#define KROME_CONST_e_mass 9.1093837015e-28  /* g; CODATA 2018 */
#define KROME_CONST_p_mass 1.67262192369e-24  /* g; CODATA 2018 */
#define KROME_CONST_n_mass 1.67492749804e-24  /* g; CODATA 2018 */
#define KROME_CONST_amu 1.6605390666e-24  /* g; CODATA 2018 */
#define KROME_CONST_clight 29979245800.0  /* cm s^-1; exact */
#define KROME_CONST_pi 3.141592653589793
#define KROME_CONST_eV_to_erg 1.602176634e-12  /* eV -> erg; exact */
#define KROME_CONST_ry_to_eV 13.605693122994  /* Rydberg -> eV; CODATA 2018 */
#define KROME_CONST_ry_to_erg 2.1798723611035473e-11  /* Rydberg -> erg; CODATA 2018 */
#define KROME_CONST_seconds_per_year 31557600.0  /* yr -> s; Julian year, 365.25 d */
#define KROME_CONST_km_to_cm 100000.0  /* km -> cm */
#define KROME_CONST_pc_to_cm 3.085677581491368e+18  /* pc -> cm; IAU 2015 */
#define KROME_CONST_Mpc_to_cm 3.085677581491368e+24  /* Mpc -> cm; IAU 2015 */
#define KROME_CONST_stefboltz_erg 5.6703744191844286e-05  /* erg s^-1 cm^-2 K^-4; exact */
#define KROME_CONST_N_avogadro 6.02214076e+23  /* mol^-1; exact */
#define KROME_CONST_Rgas_J 8.31446261815324  /* J K^-1 mol^-1; exact */
#define KROME_CONST_Rgas_kJ 0.00831446261815324  /* kJ K^-1 mol^-1; exact */

/* the integrator's failures (CVODE's values) */
#define KROME_FLAG_TOO_MUCH_WORK (-1)           /* max_steps exhausted before t_end */
#define KROME_FLAG_TOO_MUCH_ACC (-2)            /* tolerances too small for machine precision */
#define KROME_FLAG_ERR_FAILURE (-3)             /* repeated error-test failures */
#define KROME_FLAG_CONV_FAILURE (-4)            /* repeated Newton convergence failures */
#define KROME_FLAG_LSETUP_FAIL (-6)             /* Newton matrix setup / factorization failed */
#define KROME_FLAG_LSOLVE_FAIL (-7)             /* linear solve failed */
#define KROME_FLAG_RHSFUNC_FAIL (-8)            /* right-hand side failed unrecoverably */
#define KROME_FLAG_FIRST_RHSFUNC_ERR (-9)       /* right-hand side failed at the initial state */
#define KROME_FLAG_REPTD_RHSFUNC_ERR (-10)      /* repeated recoverable right-hand side errors */
#define KROME_FLAG_UNREC_RHSFUNC_ERR (-11)      /* right-hand side failed unrecoverably in a restart */
#define KROME_FLAG_NLS_FAIL (-16)               /* nonlinear solver failure */
#define KROME_FLAG_ILL_INPUT (-22)              /* illegal input (e.g. error weights not positive) */
#define KROME_FLAG_TOO_CLOSE (-27)              /* t_end too close to the start */
/* KROME's flags */
#define KROME_FLAG_OK 0
#define KROME_FLAG_NO_DUST (-102)       /* DUSTBINS cooling, no dust bins */
#define KROME_FLAG_NO_FIELD (-103)      /* photo reactions, no Phi */
#define KROME_FLAG_UNSUPPORTED (-104)   /* path/backend cannot do this */
#define KROME_FLAG_PARAM_UNSET (-105)   /* a required parameter is NaN */
#define KROME_FLAG_CONSERVE_GAS_TOTAL (-106) /* conserve with dust.gas_total */
#define KROME_FLAG_SURFACE_BINS (-107)  /* surface.bins: dust bins */
#define KROME_FLAG_NO_LAPACK (-108)     /* linsol LAPACK not built */
#define KROME_FLAG_NO_KOKKOS (-109)     /* device before initialization */
#define KROME_FLAG_NONPOSITIVE_T (-110) /* ended at Tgas <= 0: untouched */
#define KROME_FLAG_DEVICE_MEMORY (-111) /* device: no cell fits in memory */

#ifdef __cplusplus
extern "C" {
#endif

/* a one-line description of a solver flag (never NULL) */
const char* krome_flag_message(int flag);

/* network description */
const char* krome_version(void);  /* the KROME-X version, e.g.
                                     2.0.0-dev+g04bb53f */
const char* krome_network_name(void);        /* the generated network */
int krome_nspec(void);
int krome_nrea(void);
int krome_npar(void);
const char* krome_species_name(int i);
int krome_species_index(const char* name);   /* -1 if unknown */
const char* krome_param_name(int i);
int krome_param_index(const char* name);     /* -1 if unknown */
/* one row of npar default values; a REQUIRED parameter (no default
 * exists) is NaN and must be set before solving */
void krome_params_default(double* par);
/* "" for a parameter with a default, else what requires it, e.g.
 * "thermal.heating = ["COMPRESS"]" or "elec_recomb_ST93()" */
const char* krome_param_required_by(int i);

/* solver options: opaque handle, set by name. Keys (value as double;
 * booleans are 0/1):
 *   evolve_T                     thermal evolution on/off (off: Tgas, and
 *                                Tdust state if generated, stay fixed)
 *   jac                          0 auto, 1 analytic, 2 finite differences
 *   linsol                       LU of the Newton systems: 0 auto (the
 *                                sparse LU for sparse networks, else
 *                                LAPACK from 64 equations when built
 *                                with LAPACK), 1 the integrator's dense
 *                                LU, 2 LAPACK (-108 from a solve in a
 *                                build without it), 3 the sparse LU
 *   rtol, atol_n (also the atoms locked in a dust bin, dust.evolve),
 *   atol_T, atol_u
 *   max_steps                    CVODE step budget (exhausting it is final)
 *   max_err_fails, max_conv_fails  CVODE failure limits per step (7, 10)
 *   max_restarts                 CVODE restarts after error-test/convergence
 *                                failures (5), each from a tenth of the
 *                                step size CVODE had reduced to
 *   clip_negative                n = max(n, 0) after a solve (1)
 *   renormalize_x                krome_solve_cells_x: rescale the returned
 *                                mass fractions to their input sum (0;
 *                                legacy -useX did)
 *   error_report                 write a diagnostic report when a solve
 *                                fails (1; legacy KROME_ERROR_REPORT: the
 *                                state, dn/dt, sorted rates and fluxes,
 *                                solver statistics), to the file set by
 *                                krome_opts_set_error_file
 *   conserve                     post-solve conservation: 0 off (default),
 *                                1 electrons restore the input charge, 2
 *                                minimal correction restoring the input
 *                                element totals and charge, 3 electrons
 *                                neutralize the ions (even when frozen;
 *                                legacy -computeElectrons)
 *   path                         0 auto, 1 per-cell, 2 device (the
 *                                Kokkos backend: -104 in a build without
 *                                it)
 *   device_cache                 keep the device buffers between calls
 *                                (1, default; 0: allocate in every call)
 *   device_memory                the device memory a call may use
 *                                [bytes]: larger calls run in chunks
 *                                that fit (0, default: 80 % of the free
 *                                memory on CUDA / HIP; < 0 unbounded);
 *                                when no cell fits, path auto runs on
 *                                the host, path device returns -111
 *   min_device                   auto with a device backend: calls from
 *                                this many cells go to the device (1e4)
 *   min_interleave               the device layout: interleaved across
 *                                cells from this many cells, contiguous
 *                                below (-1, default: contiguous; on
 *                                CUDA interleaved above the cells
 *                                resident at once, 27 648 on an A100)
 * krome_opts_set returns -1 for any other key. Per-cell physical
 * quantities are parameters (par rows, see krome_param_name). */
typedef struct krome_opts krome_opts;
krome_opts* krome_opts_new(void);
void krome_opts_free(krome_opts* h);
int krome_opts_set(krome_opts* h, const char* key, double value);
int krome_opts_get(const krome_opts* h, const char* key, double* value);

/* per-species tolerances (legacy -customATOL / -customRTOL): atol[NSPEC]
 * [cm^-3] / rtol[NSPEC] replacing the scalars atol_n / rtol for the
 * species; a negative entry takes the scalar, NULL the scalar for every
 * species. Copied into the handle, which starts from the project's
 * defaults (codegen solver.atol_file / solver.rtol_file). Set returns -1 for a
 * NaN, or a zero atol. Get writes the tolerance each species gets (the
 * entry, or the scalar). */
int krome_opts_set_atol_species(krome_opts* h, const double* atol);
int krome_opts_set_rtol_species(krome_opts* h, const double* rtol);
int krome_opts_get_atol_species(const krome_opts* h, double* atol);
int krome_opts_get_rtol_species(const krome_opts* h, double* rtol);

/* device backend (Kokkos: libkrome built with the Makefile's KOKKOS=
 * option). krome_has_device: 1 when compiled in. Kokkos must be
 * initialized before a device solve (-109 otherwise): a code that uses
 * Kokkos itself initializes it; any other calls krome_device_initialize
 * once, and krome_device_finalize before exiting. Both return -104 in a
 * build without the device backend. */
int krome_has_device(void);
int krome_device_initialize(void);
int krome_device_finalize(void);

/* photon bins [eV edges] for photo rates, photoheating and grains;
 * nbins bins with edges Elo[b] < Ehi[b]. After this call, Phi arrays
 * passed to krome_solve_cells have nbins entries per cell. */
int krome_opts_set_photo(krome_opts* h, int nbins, const double* Elo,
                         const double* Ehi);
int krome_opts_nbins(const krome_opts* h);
/* Doppler broadening of the tabulated (file) cross sections, applied when
 * the bins are built (legacy -useBroadening, fixed at bin setup): a
 * temperature Tgas [K] and a turbulent velocity vturb [cm s^-1] (0, 0 =
 * none). Rebuilds the bins already set; -1 for negative values. */
int krome_opts_set_photo_broadening(krome_opts* h, double Tgas,
                                    double vturb);

/* binned dust population (MRN bins per grain type, radiative-balance
 * Tdust by bisection, dust-chemistry coupling); nd bins per type.
 * Drives the surface rates always, and the gas-dust cooling when the
 * network was generated with thermal.cooling = ["DUSTBINS"] (required then). */
int krome_opts_set_dustbins(krome_opts* h, int nd, int use_Si, int use_C,
                            double alow, double aup, double phi);
/* the dust bins' temperatures Td[ndust] [K] (ndust = nd x grain types; 0
 * without bins): the starting Td of every cell when krome_solve_cells gets
 * no per-cell Td, and the Td the queries use (as C++ DustBins::Td). Set
 * returns -1 without bins or for a non-positive value. */
int krome_opts_ndust(const krome_opts* h);
int krome_opts_set_dust_Td(krome_opts* h, const double* Td);
int krome_opts_get_dust_Td(const krome_opts* h, double* Td);
/* the bins' grain radii a[ndust] [cm] (--evolve-dust): the starting radii
 * of every cell without per-cell adust (the reference sizes after
 * krome_opts_set_dustbins); the grain counts nd[ndust] [cm^-3] of a cell
 * with gas density rho_gas [g cm^-3] and parameter d2g (nd = the
 * population's grains per dust mass at the reference sizes x rho_gas d2g,
 * fixed while the radii evolve); the grain type of each bin (1 carbon,
 * 0 silicate). -1 without bins (set: for a negative radius). */
int krome_opts_set_dust_a(krome_opts* h, const double* a);
int krome_opts_get_dust_a(const krome_opts* h, double* a);
int krome_opts_get_dust_counts(const krome_opts* h, double rho_gas,
                               double d2g, double* nd);
int krome_opts_get_dust_types(const krome_opts* h, int* isC);

/* the optional per-cell arrays of a multi-cell call (krome_solve_cells,
 * krome_solve_cells_x, krome_equilibrium_cells): pass NULL for none, or
 * a struct whose unused members are NULL (initialize it with {0}). New
 * per-cell fields are added here, so callers keep compiling.
 *   Phi    [ncell][nbins]  photon flux per bin, required with photon bins
 *                          on the handle (-103 otherwise): held fixed, or
 *                          updated in place with radiation.evolve
 *   Td     [ncell][ndust]  dust temperatures [K] with dust bins on the
 *                          handle: each cell starts from its row (with
 *                          dust.evolve_Td its initial state) and the row
 *                          receives its final Td; NULL: every cell starts
 *                          from the handle's (krome_opts_set_dust_Td), the
 *                          results are discarded
 *   adust  [ncell][ndust]  grain radii [cm] with dust.evolve, read and
 *                          updated like Td (0: an evaporated bin); NULL:
 *                          from the handle's radii (krome_opts_set_dust_a)
 * A failed cell's rows are untouched. */
typedef struct {
  double* Phi;
  double* Td;
  double* adust;
} krome_cells;

/* integrate ncell independent cells, each over its own dt[c] [s].
 * n and Tgas (when evolve_T) are updated in place; extra: the optional
 * per-cell arrays (krome_cells, may be NULL). A cell that fails (after
 * the restarts) is left untouched. flags[ncell] (optional, may be NULL)
 * receives per-cell solver flags. Returns 0, or the flag of the first
 * failing cell; -103: photon bins on the handle and no extra->Phi; -105:
 * a REQUIRED parameter is unset (NaN) in the flagged cells, nothing
 * integrated; -107: generated with surface.bins = NB and the handle's dust
 * bins (krome_opts_set_dustbins) are missing or not NB; -108: linsol
 * LAPACK in a build without LAPACK. */
int krome_solve_cells(int ncell, double* n, double* Tgas, const double* dt,
                      const double* par, int par_stride,
                      const krome_cells* extra, int* flags,
                      const krome_opts* h);

/* the same with mass fractions x[ncell][nspec] and the per-cell gas
 * density rho[ncell] [g cm^-3] (legacy -useX): n = rho x / m on the way
 * in, x = m n / rho on the way out (the caller's rho; rescaled to the
 * input sum with the renormalize_x key; see also the conserve option),
 * for the cells that succeed; a failed cell keeps its input x. The
 * solver works on number densities. */
int krome_solve_cells_x(int ncell, double* x, const double* rho,
                        double* Tgas, const double* dt, const double* par,
                        int par_stride, const krome_cells* extra,
                        int* flags, const krome_opts* h);

/* the file of the error report (default "KROME_ERROR_REPORT", replaced at
 * each failed solve; in an MPI run give each rank its own name) */
int krome_opts_set_error_file(krome_opts* h, const char* path);

/* integrate every cell to steady state (legacy krome_equilibrium), with
 * the krome_solve_cells conventions; isothermal unless evolve_T. Handle
 * keys: eq_t0 (first output, 100 yr), eq_tmax (5e8 yr), eq_factor (3),
 * eq_rtol (1e-5), eq_atol (1e-12). flags[c]: 0 converged, 1 eq_tmax
 * reached first (state returned), <0 failure (cell untouched). extra as
 * in krome_solve_cells. Returns
 * the flag of the first cell that did not converge (0 if all did). */
int krome_equilibrium_cells(int ncell, double* n, double* Tgas,
                            const double* par, int par_stride,
                            const krome_cells* extra, int* flags,
                            const krome_opts* h);

/* ---- query / diagnostic functions (legacy krome_get_*) ------------------
 * Same input conventions as krome_solve_cells (n[ncell][nspec], par rows
 * or one broadcast row with par_stride 0, Phi[ncell][nbins] when photon
 * bins are configured on h; h may be NULL for default options). Values are
 * evaluated as the solver sees them (temperature capped to [Tcmb, 1e9 K]
 * for the rates, to [2.73 K, 1e9 K] for the thermal terms and gamma).
 * Functions returning int return 0, -1 (bad arguments), -103 (photo bins
 * configured, no Phi) or -105 (a REQUIRED parameter unset). */
/* rate coefficients k[ncell][nrea] and fluxes k*prod(n) [cm^-3 s^-1] */
int krome_get_rates(int ncell, const double* n, const double* Tgas,
                    const double* par, int par_stride, const double* Phi,
                    const krome_opts* h, double* k);
int krome_get_fluxes(int ncell, const double* n, const double* Tgas,
                     const double* par, int par_stride, const double* Phi,
                     const krome_opts* h, double* fluxes);
/* the dust temperature the rates see, Td[ncell] [K] (legacy
 * krome_get_table_Tdust): the dust.tables value, else the Tdust
 * parameter (the binned-dust mean with dust bins on h); -1 when the
 * network has no dust temperature */
int krome_get_Tdust(int ncell, const double* n, const double* Tgas,
                    const double* par, int par_stride, const double* Phi,
                    const krome_opts* h, double* Td);
/* thermal balance [erg cm^-3 s^-1]: totals heat[ncell], cool[ncell], the
 * individual terms heats[ncell][nheat], cools[ncell][ncool] (names:
 * krome_heating_name / krome_cooling_name) and dTdt[ncell] [K/s]; any
 * output may be NULL. The FLOOR cooling term is the (negative) floor
 * correction; CHEM is the H2 formation heating / dissociation cooling. */
int krome_nheat(void);
int krome_ncool(void);
const char* krome_heating_name(int i);
const char* krome_cooling_name(int i);
int krome_get_thermal(int ncell, const double* n, const double* Tgas,
                      const double* par, int par_stride, const double* Phi,
                      const krome_opts* h, double* heat, double* cool,
                      double* heats, double* cools, double* dTdt);
/* energy released per reaction event q[ncell][nrea] [MeV] (; qeff=, 0
 * elsewhere; legacy krome_get_qeff); thermal.heating = ["NUCLEAR"] adds sum flux*qeff */
int krome_get_qeff(int ncell, const double* n, const double* Tgas,
                   const double* par, int par_stride, double* q);
/* adiabatic index (codegen thermal.gamma) */
int krome_get_gamma(int ncell, const double* n, const double* Tgas,
                    const double* par, int par_stride, double* gamma);
/* per-cell scalars from n[ncell][nspec] into out[ncell]: mean molecular
 * weight, mass density [g cm^-3], total number density, H nuclei, net
 * charge [e cm^-3], the electron density for neutrality (legacy
 * get_electrons), free-fall time [s]; Jeans length [cm] needs Tgas */
int krome_get_mu(int ncell, const double* n, double* out);
int krome_get_rho(int ncell, const double* n, double* out);
int krome_get_ntot(int ncell, const double* n, double* out);
int krome_get_Hnuclei(int ncell, const double* n, double* out);
int krome_get_charge(int ncell, const double* n, double* out);
int krome_get_electrons(int ncell, const double* n, double* out);
int krome_get_free_fall_time(int ncell, const double* n, double* out);
int krome_get_jeans_length(int ncell, const double* n, const double* Tgas,
                           double* out);
/* the free-fall time [s] from the mass density rho[ncell] [g cm^-3]
 * (legacy krome_get_free_fall_time_rho); the column density [cm^-2] of the
 * number densities num[ncell] [cm^-3] in each cell with the network's
 * estimator (radiation.num2col; legacy krome_num2col) */
int krome_get_free_fall_time_rho(int ncell, const double* rho, double* out);
int krome_num2col(int ncell, const double* num, const double* n,
                  const double* Tgas, const double* par, int par_stride,
                  double* out);
/* the line coolants: the coolZ.dat metals of -cooling=Z, then the
 * thermal.coolant_files coolants. Their number, legacy names (CII, FeII, O2), species
 * (0-based index), level count, and levels (energies E[nlev] [K],
 * degeneracies g[nlev], in the file's numbering); per cell each coolant's
 * line cooling cool[ncell][ncoolant] [erg cm^-3 s^-1] (legacy
 * krome_cooling<metal>; the Z and coolant-file cooling terms are their
 * sums) and the level populations of coolant i, pop[ncell][nlev_i]
 * [cm^-3] (legacy krome_popcool_dump) */
int krome_ncoolant(void);
const char* krome_coolant_name(int i);
int krome_coolant_species(int i);
int krome_coolant_nlev(int i);
int krome_coolant_levels(int i, double* E, double* g);
int krome_get_coolant_cooling(int ncell, const double* n, const double* Tgas,
                              const double* par, int par_stride,
                              const double* Phi, const krome_opts* h,
                              double* cool);
int krome_get_coolant_populations(int ncell, const double* n,
                                  const double* Tgas, const double* par,
                                  int par_stride, const double* Phi,
                                  const krome_opts* h, int i, double* pop);
/* photoheating per photo reaction heat[ncell][nphoto] [erg cm^-3 s^-1],
 * summing to the PHOTO heating term (legacy krome_get_photoBin_heats gave
 * it per absorber: divide by n[krome_photo_reactant(r)]); the field over
 * the Draine ISRF of Phi[ncell][nbins] (legacy krome_get_ratioFluxDraine).
 * -1 without photon bins on h */
int krome_get_photo_heating(int ncell, const double* n, const double* Tgas,
                            const double* par, int par_stride,
                            const double* Phi, const krome_opts* h,
                            double* heat);
int krome_get_ratio_flux_draine(int ncell, const double* Phi,
                                const krome_opts* h, double* ratio);
/* element totals tot[ncell][nelem] [atoms cm^-3] (krome_element_name) */
int krome_get_element_totals(int ncell, const double* n, double* tot);
/* number densities <-> mass fractions x = m n / rho; rho[ncell] */
int krome_n2x(int ncell, const double* n, double* x);

/* ---- legacy conveniences (krome_user) ----------------------------------
 * one cell: n[nspec], par one row (NULL: defaults), Phi[nbins] with photon
 * bins on h. The species names separated by spaces (legacy
 * krome_get_names_header); the nbest largest reaction fluxes to stdout
 * (legacy krome_print_best_flux); every reaction's flux, absolute and
 * relative, appended to the file path, labelled by xvar (legacy
 * krome_explore_flux); the species table to stdout (legacy krome_get_info);
 * ncell cells' neutral atomic metals set to the Asplund+09 abundances times
 * 10^Z (legacy krome_scale_Z). Functions returning int return 0, -1 (bad
 * arguments or an unwritable file) or a query flag. */
const char* krome_names_header(void);
int krome_print_best_flux(const double* n, double Tgas, const double* par,
                          const double* Phi, const krome_opts* h,
                          int nbest);
int krome_explore_flux(const char* path, const double* n, double Tgas,
                       const double* par, const double* Phi,
                       const krome_opts* h, double xvar);
void krome_print_info(const double* n, double Tgas);
int krome_scale_Z(int ncell, double* n, double Z);
/* appended to the file path: every coolant's level populations, one row
 * per level (name, level from 1, Tgas, population, total; legacy
 * krome_popcool_dump); one row Tgas, total cooling, every cooling term
 * (legacy krome_dump_cooling); the rate coefficients on npts log
 * temperatures Tmin..Tmax at n = 1e-40, rows Tgas, reaction from 1, k
 * (legacy krome_dump_rates). Every parameter of the row par (NULL: the
 * defaults), 'name: value', to stdout (legacy krome_print_phys_variables) */
int krome_print_coolant_populations(const char* path, const double* n,
                                    double Tgas, const double* par,
                                    const double* Phi, const krome_opts* h);
int krome_dump_cooling(const char* path, const double* n, double Tgas,
                       const double* par, const double* Phi,
                       const krome_opts* h);
int krome_dump_rates(const char* path, double Tmin, double Tmax, int npts,
                     const double* par, const double* Phi,
                     const krome_opts* h);
void krome_print_params(const double* par);
int krome_x2n(int ncell, const double* x, const double* rho, double* n);

/* ---- radiation (legacy krome_set_photoBin_*, krome_get_opacity*) ---------
 * The field is the caller's Phi[nbins] (photon flux per bin [cm^-2 s^-1]);
 * these build binnings and fields for it. Bins: edges Elo[b] < Ehi[b] [eV],
 * Elo[b] > 0. Every spectrum is the bin integral
 *   Phi_b = 4 pi / h  int_bin J(E) d ln E,  J [eV s^-1 cm^-2 Hz^-1 sr^-1]
 * (legacy sampled J at the bin midpoint). Returns 0, or -1 bad arguments. */
/* nbins linear / logarithmic bins in [lower, upper] eV (lower > 0) */
int krome_bins_lin(int nbins, double lower, double upper, double* Elo,
                   double* Ehi);
int krome_bins_log(int nbins, double lower, double upper, double* Elo,
                   double* Ehi);
/* J constant across each bin (legacy photoBinJ) <-> Phi:
 * Phi_b = 4 pi J_b ln(Ehi/Elo) / h */
int krome_J_to_phi(int nbins, const double* Elo, const double* Ehi,
                   const double* J, double* phi);
int krome_phi_to_J(int nbins, const double* Elo, const double* Ehi,
                   const double* phi, double* J);
/* undiluted blackbody at Tbb [K] (multiply by a dilution factor as
 * needed); krome_blackbody_range: the legacy BBlog_auto energy range where
 * B_nu > 1e-6 (-1 when there is none) */
int krome_spectrum_blackbody(int nbins, const double* Elo, const double* Ehi,
                             double Tbb, double* phi);
int krome_blackbody_range(double Tbb, double* lower, double* upper);
/* Draine 1978 ISRF (5-13.6 eV), exact bin integral */
int krome_spectrum_draine(int nbins, const double* Elo, const double* Ehi,
                          double* phi);
/* J = J21 1e-21 erg s^-1 cm^-2 Hz^-1 sr^-1 (13.6 eV/E)^alpha (legacy J21:
 * J21 = 1, alpha = 1), exact bin integral */
int krome_spectrum_j21(int nbins, const double* Elo, const double* Ehi,
                       double J21, double alpha, double* phi);
/* Haardt & Madau 2012 UV background at redshift z in [0, 15.66] (-1
 * outside); zero outside the tabulated 0.1237 eV - 49.97 MeV */
int krome_spectrum_hm2012(int nbins, const double* Elo, const double* Ehi,
                          double z, double* phi);
/* tabulated J(E) at nE ascending energies [eV], linearly (loglog 0) or
 * log-log (loglog 1) interpolated, zero outside the table (legacy
 * krome_load_photoBin_file_2col, with the file read by the caller) */
int krome_spectrum_table(int nbins, const double* Elo, const double* Ehi,
                         int nE, const double* E, const double* J,
                         int loglog, double* phi);
/* bin-mean dust opacity kbin[nbins] [cm^2 g^-1 of dust] from a table
 * kappa(E) at nE ascending energies [eV] (legacy load_opacity_table), or
 * from the embedded Draine R_V = 3.1 absorption opacity (legacy
 * krome_load_average_kabs); -1 also when a bin is outside the table */
int krome_dust_kappa_bins(int nE, const double* E, const double* kappa,
                          int nbins, const double* Elo, const double* Ehi,
                          double* kbin);
int krome_dust_kappa_draine(int nbins, const double* Elo, const double* Ehi,
                            double* kbin);
/* G0 [Draine units] and Av of a field (legacy krome_find_G0_Av): the least
 * squares fit of ln(Phi_b / Phi_Draine,b) = ln G0 - kbin_b d2g 1.8e21 m_p
 * Av / xH over the Draine bins; ln G0 limited to [-7, 7], Av >= 0. kbin:
 * dust opacity per bin (krome_dust_kappa_*); xH: H mass fraction
 * Hnuclei m_p / rho. -1 also with fewer than 2 bins in 5-13.6 eV */
int krome_find_G0_Av(int nbins, const double* Elo, const double* Ehi,
                     const double* phi, const double* kbin, double d2g,
                     double xH, double* G0, double* Av);
/* photo reactions: count, reaction and reactant index, events per absorbed
 * photon (; pdiss=, 1 by default), cross section sigma(E) [cm^2] (NaN for
 * a bad r), and the bin-averaged cross sections sigma[nphoto][nbins] of
 * the bins configured on h (unshielded) */
int krome_nphoto(void);
int krome_photo_reaction(int r);
int krome_photo_reactant(int r);
double krome_photo_pdiss(int r);
/* J_LW rates (k = khat J_LW, the Lyman-Werner band mean intensity of the
 * bins): count, reaction, absorber, pdiss, and the band weights w[nbins]
 * of h's bins (J_LW = sum_b w_b Phi_b). Put bin edges at 11.2 and 13.6 eV
 * (network-format section 5.4). */
int krome_nlw(void);
int krome_lw_reaction(int q);
int krome_lw_reactant(int q);
double krome_lw_pdiss(int q);
int krome_opts_lw_weights(const krome_opts* h, double* w);
double krome_photo_xsec(int r, double E);
int krome_opts_photo_sigma(const krome_opts* h, double* sigma);
/* Habing flux G[ncell] [1.6e-3 erg cm^-2 s^-1, 6-13.6 eV] of Phi */
int krome_get_habing(int ncell, const double* Phi, const krome_opts* h,
                     double* G);
/* optical depth tau[ncell][nbins] with the columns of the generated
 * num2col (codegen radiation.num2col; legacy krome_get_opacity): every photo
 * reaction and J_LW rate contributes sigma_b / pdiss times its absorber
 * column, shielded ones through the curve of growth sigma_b / pdiss
 * int_0^N s dN'. With d2g[ncell] (NULL: no dust) the dust adds kappa_b d2g
 * Sigma_gas, Sigma_gas = N_H rho / n_H (N_H the num2col H, H+, 2 H2
 * column): the dust bins configured on h, else kdust[nbins] [cm^2 per g
 * of DUST], by default (NULL) the Draine R_V = 3.1 opacity. Photon bins
 * required on h. */
int krome_get_opacity(int ncell, const double* n, const double* Tgas,
                      const double* par, int par_stride,
                      const krome_opts* h, const double* d2g,
                      const double* kdust, double* tau);
/* photon removal coefficient kabs[ncell][nbins] [cm^-1] (no length): the
 * photochemistry with its shielding factors and 1/pdiss photons per event
 * (the radiation.evolve sink; a radiation.evolve build has already removed it from
 * the returned field), the dust bins configured on h, and rho_gas kdust[b]
 * when kdust[nbins] (dust opacity per unit GAS mass [cm^2 g^-1], e.g. d2g
 * times krome_dust_kappa_*) is not NULL. Times the cell size it is legacy
 * krome_get_opacity_size(_d2g). */
int krome_get_absorption(int ncell, const double* n, const double* Tgas,
                         const double* par, int par_stride,
                         const krome_opts* h, const double* kdust,
                         double* kabs);

/* ---- species / reaction metadata ----------------------------------------- */
double krome_species_mass(int i);            /* [g]; NaN out of range */
int krome_species_charge(int i);
int krome_nelem(void);
const char* krome_element_name(int e);
int krome_species_natoms(int i, int e);      /* atoms of element e; -1 bad */
const char* krome_reaction_name(int r);      /* equation (+ T range) */

#ifdef __cplusplus
}
#endif
#endif /* KROME_H */
