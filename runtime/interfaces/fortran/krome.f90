! KROME-X Fortran interface: iso_c_binding bindings to the C interface
! (krome.h) plus Fortran-friendly wrappers. The generated module krome_net
! (krome_net.f90) holds the network constants: krome_nspec, krome_npar,
! krome_idx_<species> and krome_par_<param>, all 1-BASED.
!
! Array layouts (no transposes: Fortran column-major == C row-major):
!   n(krome_nspec, ncell)   number densities [cm^-3]
!   par(krome_npar, ncell)  per-cell runtime parameters, or par(krome_npar)
!                           broadcast to every cell
!   Phi(nbins, ncell)       photon flux per bin: held fixed, or
!                           co-evolved and updated in place when the
!                           network was generated with radiation.evolve
!
! Minimal use:
!   use krome
!   type(c_ptr) :: h
!   h = krome_opts_new()
!   call krome_set(h, "evolve_T", 1d0)  ! thermal evolution on
!   call krome_params_default(par)      ! then per-cell edits (parameters
!                                       ! without a default are NaN and
!                                       ! must be set), e.g.
!   par(krome_par_crate, :) = ...
!   ierr = krome_solve(n, Tgas, dt, h, par=par)
!   call krome_opts_free(h)
module krome
  use iso_c_binding
  implicit none
  private
  ! solver return flags (krome.h KROME_FLAG_*: the integrator's failures,
  ! -1 .. -27, and KROME's, <= -100; krome_flag_message_f describes any)
  integer, parameter, public :: krome_flag_too_much_work = -1
  integer, parameter, public :: krome_flag_too_much_acc = -2
  integer, parameter, public :: krome_flag_err_failure = -3
  integer, parameter, public :: krome_flag_conv_failure = -4
  integer, parameter, public :: krome_flag_lsetup_fail = -6
  integer, parameter, public :: krome_flag_lsolve_fail = -7
  integer, parameter, public :: krome_flag_rhsfunc_fail = -8
  integer, parameter, public :: krome_flag_first_rhsfunc_err = -9
  integer, parameter, public :: krome_flag_reptd_rhsfunc_err = -10
  integer, parameter, public :: krome_flag_unrec_rhsfunc_err = -11
  integer, parameter, public :: krome_flag_nls_fail = -16
  integer, parameter, public :: krome_flag_ill_input = -22
  integer, parameter, public :: krome_flag_too_close = -27
  integer, parameter, public :: krome_flag_ok = 0
  integer, parameter, public :: krome_flag_no_dust = -102
  integer, parameter, public :: krome_flag_no_field = -103
  integer, parameter, public :: krome_flag_unsupported = -104
  integer, parameter, public :: krome_flag_param_unset = -105
  integer, parameter, public :: krome_flag_conserve_gas_total = -106
  integer, parameter, public :: krome_flag_surface_bins = -107
  integer, parameter, public :: krome_flag_no_lapack = -108
  integer, parameter, public :: krome_flag_no_kokkos = -109
  integer, parameter, public :: krome_flag_nonpositive_t = -110
  integer, parameter, public :: krome_flag_device_memory = -111
  integer, parameter, public :: krome_flag_conserve_unreached = -112
  integer, parameter, public :: krome_flag_conserve_large = 2
  integer, parameter, public :: krome_flag_host_recovered = 1
  ! the C kinds used in the interface, re-exported for callers
  public :: c_ptr, c_null_ptr, c_int, c_double, c_long, c_loc
  public :: krome_opts_new, krome_opts_free, krome_nspec_c, krome_npar_c, &
       krome_nrea_c
  ! device backend (Kokkos build): 1 / 0; initialize / finalize Kokkos
  ! when the caller does not (-104 without the backend)
  public :: krome_has_device, krome_device_initialize, &
       krome_device_finalize, krome_device_initialize_id
  ! device-resident data (krome.h krome_solve_cells_device; manual 6.8):
  ! every array in the device's memory, passed as a C pointer (CUDA
  ! Fortran: c_loc of a device array; OpenACC: c_loc inside host_data
  ! use_device). Species-major: n(ld, nspec), the default; cell-major:
  ! n(nspec, ncell). The rows (Phi, Td, adust) in a krome_cells_c of device
  ! pointers; the summary on the host
  integer, parameter, public :: krome_layout_species_major = 0
  integer, parameter, public :: krome_layout_cell_major = 1
  public :: krome_cells_c, krome_device_summary, krome_solve_cells_device, &
       krome_solve_cells_x_device, krome_conserve_cells_device, &
       krome_conserve_cells_x_device
  public :: krome_set, krome_get, krome_species_idx, krome_param_idx, &
       krome_species_name_f, krome_params_default, krome_set_photo, &
       krome_version_f, &
       krome_set_dustbins, krome_nbins, krome_solve, krome_solve_x, &
       krome_set_photo_broadening, krome_ndust, krome_set_dust_Td, &
       krome_set_dust_a, krome_get_dust_a, krome_get_dust_counts, &
       krome_get_dust_types, &
       krome_get_dust_Td, &
       krome_set_atol_species, krome_set_rtol_species, &
       krome_get_atol_species, krome_get_rtol_species, &
       krome_set_error_file, krome_flag_message_f
  ! queries and metadata (legacy krome_get_*); 1-based indices
  public :: krome_param_name_f, krome_param_required_by_f, &
       krome_reaction_name_f, krome_element_name_f, krome_heating_name_f, &
       krome_cooling_name_f, krome_nelem_c, krome_nheat_c, krome_ncool_c, &
       krome_species_mass_f, krome_species_charge_f, krome_species_natoms_f
  public :: krome_get_rates, krome_get_fluxes, krome_get_thermal, &
       krome_get_Tdust, &
       krome_get_gamma, krome_get_mu, krome_get_rho, krome_get_ntot, &
       krome_get_Hnuclei, krome_get_charge, krome_get_electrons, &
       krome_get_free_fall_time, krome_get_jeans_length, &
       krome_get_element_totals, krome_get_element_totals_x, &
       krome_conserve, krome_conserve_x, &
       krome_n2x, krome_x2n, krome_equilibrium, &
       krome_get_qeff, krome_get_free_fall_time_rho, krome_num2col, &
       krome_ncoolant_c, krome_coolant_name_f, krome_coolant_species_f, &
       krome_coolant_nlev_f, krome_coolant_levels, &
       krome_get_coolant_cooling, krome_get_coolant_populations, &
       krome_get_photo_heating, krome_get_ratio_flux_draine
  ! legacy conveniences (krome_user), one cell n(nspec)
  public :: krome_names_header, krome_print_best_flux, krome_explore_flux, &
       krome_print_info, krome_scale_Z, krome_print_coolant_populations, &
       krome_dump_cooling, krome_dump_rates, krome_print_params
  ! legacy krome_scale_Z: one cell n(nspec) or cells n(nspec, ncell)
  interface krome_scale_Z
    module procedure scale_z_1_, scale_z_2_
  end interface krome_scale_Z
  ! radiation (legacy krome_set_photoBin_*, krome_get_opacity*)
  public :: krome_bins_lin, krome_bins_log, krome_J_to_phi, krome_phi_to_J, &
       krome_spectrum_blackbody, krome_blackbody_range, &
       krome_spectrum_draine, krome_spectrum_j21, krome_spectrum_hm2012, &
       krome_spectrum_table, krome_dust_kappa_bins, krome_dust_kappa_draine, &
       krome_find_G0_Av, krome_nphoto_c, krome_photo_reaction_f, &
       krome_photo_reactant_f, krome_photo_xsec_f, krome_photo_sigma, &
       krome_get_habing, krome_get_opacity, krome_get_absorption

  ! the optional per-cell arrays of a multi-cell call (krome.h krome_cells)
  type, bind(C) :: krome_cells_c
    type(c_ptr) :: phi = c_null_ptr
    type(c_ptr) :: td = c_null_ptr
    type(c_ptr) :: adust = c_null_ptr
  end type krome_cells_c
  ! the counts of a device-resident call (krome.h krome_device_summary)
  type, bind(C) :: krome_device_summary
    integer(c_long) :: device = 0  ! solved on the device
    integer(c_long) :: host = 0    ! failed there, solved on the host
    integer(c_long) :: failed = 0  ! failed on both
  end type krome_device_summary

  interface
    function krome_has_device() bind(C, name="krome_has_device") result(r)
      import :: c_int
      integer(c_int) :: r
    end function krome_has_device
    function krome_device_initialize() bind(C, name="krome_device_initialize") &
        result(r)
      import :: c_int
      integer(c_int) :: r
    end function krome_device_initialize
    function krome_device_finalize() bind(C, name="krome_device_finalize") &
        result(r)
      import :: c_int
      integer(c_int) :: r
    end function krome_device_finalize
    function krome_device_initialize_id(device_id) &
        bind(C, name="krome_device_initialize_id") result(r)
      import :: c_int
      integer(c_int), value :: device_id
      integer(c_int) :: r
    end function krome_device_initialize_id
    function krome_solve_cells_device(ncell, n, layout, ld, tgas, dt, par, &
        par_stride, ex, flags, summary, h) &
        bind(C, name="krome_solve_cells_device") result(r)
      import :: c_ptr, c_int, c_long, krome_cells_c, krome_device_summary
      integer(c_int), value :: ncell, layout, par_stride
      integer(c_long), value :: ld
      type(c_ptr), value :: n, tgas, dt, par, flags, h
      type(krome_cells_c), intent(in) :: ex
      type(krome_device_summary), intent(out) :: summary
      integer(c_int) :: r
    end function krome_solve_cells_device
    function krome_solve_cells_x_device(ncell, x, rho, layout, ld, tgas, dt, &
        par, par_stride, ex, flags, summary, h) &
        bind(C, name="krome_solve_cells_x_device") result(r)
      import :: c_ptr, c_int, c_long, krome_cells_c, krome_device_summary
      integer(c_int), value :: ncell, layout, par_stride
      integer(c_long), value :: ld
      type(c_ptr), value :: x, rho, tgas, dt, par, flags, h
      type(krome_cells_c), intent(in) :: ex
      type(krome_device_summary), intent(out) :: summary
      integer(c_int) :: r
    end function krome_solve_cells_x_device
    function krome_conserve_cells_device(ncell, n, layout, ld, ref, &
        ref_stride, par, par_stride, flags) &
        bind(C, name="krome_conserve_cells_device") result(r)
      import :: c_ptr, c_int, c_long
      integer(c_int), value :: ncell, layout, ref_stride, par_stride
      integer(c_long), value :: ld
      type(c_ptr), value :: n, ref, par, flags
      integer(c_int) :: r
    end function krome_conserve_cells_device
    function krome_conserve_cells_x_device(ncell, x, layout, ld, ref, &
        ref_stride, par, par_stride, flags) &
        bind(C, name="krome_conserve_cells_x_device") result(r)
      import :: c_ptr, c_int, c_long
      integer(c_int), value :: ncell, layout, ref_stride, par_stride
      integer(c_long), value :: ld
      type(c_ptr), value :: x, ref, par, flags
      integer(c_int) :: r
    end function krome_conserve_cells_x_device
    function krome_nspec_c() bind(C, name="krome_nspec") result(r)
      import :: c_int
      integer(c_int) :: r
    end function krome_nspec_c
    function krome_nrea_c() bind(C, name="krome_nrea") result(r)
      import :: c_int
      integer(c_int) :: r
    end function krome_nrea_c
    function krome_npar_c() bind(C, name="krome_npar") result(r)
      import :: c_int
      integer(c_int) :: r
    end function krome_npar_c
    function c_version() bind(C, name="krome_version") result(r)
      import :: c_ptr
      type(c_ptr) :: r
    end function c_version
    function c_species_name(i) bind(C, name="krome_species_name") result(r)
      import :: c_int, c_ptr
      integer(c_int), value :: i
      type(c_ptr) :: r
    end function c_species_name
    function c_species_index(name) bind(C, name="krome_species_index") &
         result(r)
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int) :: r
    end function c_species_index
    function c_param_index(name) bind(C, name="krome_param_index") result(r)
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int) :: r
    end function c_param_index
    subroutine c_params_default(par) bind(C, name="krome_params_default")
      import :: c_double
      real(c_double) :: par(*)
    end subroutine c_params_default
    function krome_opts_new() bind(C, name="krome_opts_new") result(h)
      import :: c_ptr
      type(c_ptr) :: h
    end function krome_opts_new
    subroutine krome_opts_free(h) bind(C, name="krome_opts_free")
      import :: c_ptr
      type(c_ptr), value :: h
    end subroutine krome_opts_free
    function c_opts_set(h, key, v) bind(C, name="krome_opts_set") result(r)
      import :: c_ptr, c_char, c_double, c_int
      type(c_ptr), value :: h
      character(kind=c_char), intent(in) :: key(*)
      real(c_double), value :: v
      integer(c_int) :: r
    end function c_opts_set
    function c_opts_get(h, key, v) bind(C, name="krome_opts_get") result(r)
      import :: c_ptr, c_char, c_double, c_int
      type(c_ptr), value :: h
      character(kind=c_char), intent(in) :: key(*)
      real(c_double) :: v
      integer(c_int) :: r
    end function c_opts_get
    function c_set_photo(h, nb, elo, ehi) bind(C, name="krome_opts_set_photo") &
         result(r)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: h
      integer(c_int), value :: nb
      real(c_double), intent(in) :: elo(*), ehi(*)
      integer(c_int) :: r
    end function c_set_photo
    function c_set_photo_broadening(h, T, v) &
         bind(C, name="krome_opts_set_photo_broadening") result(r)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: h
      real(c_double), value :: T, v
      integer(c_int) :: r
    end function c_set_photo_broadening
    function c_nbins(h) bind(C, name="krome_opts_nbins") result(r)
      import :: c_ptr, c_int
      type(c_ptr), value :: h
      integer(c_int) :: r
    end function c_nbins
    function c_set_dustbins(h, nd, si, cc, alow, aup, phi) &
         bind(C, name="krome_opts_set_dustbins") result(r)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: h
      integer(c_int), value :: nd, si, cc
      real(c_double), value :: alow, aup, phi
      integer(c_int) :: r
    end function c_set_dustbins
    function c_solve_cells(ncell, n, tgas, dt, par, stride, ex, flags, h) &
         bind(C, name="krome_solve_cells") result(r)
      import :: c_ptr, c_int, c_double, krome_cells_c
      integer(c_int), value :: ncell, stride
      real(c_double) :: n(*), tgas(*)
      real(c_double), intent(in) :: dt(*)
      type(krome_cells_c), intent(in) :: ex
      type(c_ptr), value :: par, flags, h
      integer(c_int) :: r
    end function c_solve_cells
    function c_solve_cells_x(ncell, x, rho, tgas, dt, par, stride, ex, &
         flags, h) bind(C, name="krome_solve_cells_x") result(r)
      import :: c_ptr, c_int, c_double, krome_cells_c
      integer(c_int), value :: ncell, stride
      real(c_double) :: x(*), tgas(*)
      real(c_double), intent(in) :: rho(*), dt(*)
      type(krome_cells_c), intent(in) :: ex
      type(c_ptr), value :: par, flags, h
      integer(c_int) :: r
    end function c_solve_cells_x
    function c_ndust(h) bind(C, name="krome_opts_ndust") result(r)
      import :: c_ptr, c_int
      type(c_ptr), value :: h
      integer(c_int) :: r
    end function c_ndust
    function c_set_dust_td(h, td) bind(C, name="krome_opts_set_dust_Td") &
         result(r)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: h
      real(c_double), intent(in) :: td(*)
      integer(c_int) :: r
    end function c_set_dust_td
    function c_get_dust_td(h, td) bind(C, name="krome_opts_get_dust_Td") &
         result(r)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: h
      real(c_double), intent(out) :: td(*)
      integer(c_int) :: r
    end function c_get_dust_td
    function c_set_atol_sp(h, v) bind(C, name="krome_opts_set_atol_species") &
         result(r)
      import :: c_ptr, c_int
      type(c_ptr), value :: h, v
      integer(c_int) :: r
    end function c_set_atol_sp
    function c_set_rtol_sp(h, v) bind(C, name="krome_opts_set_rtol_species") &
         result(r)
      import :: c_ptr, c_int
      type(c_ptr), value :: h, v
      integer(c_int) :: r
    end function c_set_rtol_sp
    function c_get_atol_sp(h, v) bind(C, name="krome_opts_get_atol_species") &
         result(r)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: h
      real(c_double), intent(out) :: v(*)
      integer(c_int) :: r
    end function c_get_atol_sp
    function c_get_rtol_sp(h, v) bind(C, name="krome_opts_get_rtol_species") &
         result(r)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: h
      real(c_double), intent(out) :: v(*)
      integer(c_int) :: r
    end function c_get_rtol_sp
    function c_set_dust_a(h, a) bind(C, name="krome_opts_set_dust_a") &
         result(r)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: h
      real(c_double), intent(in) :: a(*)
      integer(c_int) :: r
    end function c_set_dust_a
    function c_get_dust_a(h, a) bind(C, name="krome_opts_get_dust_a") &
         result(r)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: h
      real(c_double), intent(out) :: a(*)
      integer(c_int) :: r
    end function c_get_dust_a
    function c_get_dust_counts(h, rho, d2g, nd) &
         bind(C, name="krome_opts_get_dust_counts") result(r)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: h
      real(c_double), value :: rho, d2g
      real(c_double), intent(out) :: nd(*)
      integer(c_int) :: r
    end function c_get_dust_counts
    function c_get_dust_types(h, isc) bind(C, name="krome_opts_get_dust_types") &
         result(r)
      import :: c_ptr, c_int
      type(c_ptr), value :: h
      integer(c_int), intent(out) :: isc(*)
      integer(c_int) :: r
    end function c_get_dust_types
    function c_set_error_file(h, path) bind(C, name="krome_opts_set_error_file") &
         result(r)
      import :: c_ptr, c_int, c_char
      type(c_ptr), value :: h
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: r
    end function c_set_error_file
    function c_flag_message(flag) bind(C, name="krome_flag_message") result(r)
      import :: c_int, c_ptr
      integer(c_int), value :: flag
      type(c_ptr) :: r
    end function c_flag_message
    ! names by 0-based index (C strings)
    function c_param_name(i) bind(C, name="krome_param_name") result(r)
      import :: c_int, c_ptr
      integer(c_int), value :: i
      type(c_ptr) :: r
    end function c_param_name
    function c_param_required_by(i) bind(C, name="krome_param_required_by") &
         result(r)
      import :: c_int, c_ptr
      integer(c_int), value :: i
      type(c_ptr) :: r
    end function c_param_required_by
    function c_reaction_name(i) bind(C, name="krome_reaction_name") result(r)
      import :: c_int, c_ptr
      integer(c_int), value :: i
      type(c_ptr) :: r
    end function c_reaction_name
    function c_element_name(i) bind(C, name="krome_element_name") result(r)
      import :: c_int, c_ptr
      integer(c_int), value :: i
      type(c_ptr) :: r
    end function c_element_name
    function c_heating_name(i) bind(C, name="krome_heating_name") result(r)
      import :: c_int, c_ptr
      integer(c_int), value :: i
      type(c_ptr) :: r
    end function c_heating_name
    function c_cooling_name(i) bind(C, name="krome_cooling_name") result(r)
      import :: c_int, c_ptr
      integer(c_int), value :: i
      type(c_ptr) :: r
    end function c_cooling_name
    function krome_nelem_c() bind(C, name="krome_nelem") result(r)
      import :: c_int
      integer(c_int) :: r
    end function krome_nelem_c
    function krome_nheat_c() bind(C, name="krome_nheat") result(r)
      import :: c_int
      integer(c_int) :: r
    end function krome_nheat_c
    function krome_ncool_c() bind(C, name="krome_ncool") result(r)
      import :: c_int
      integer(c_int) :: r
    end function krome_ncool_c
    function c_species_mass(i) bind(C, name="krome_species_mass") result(r)
      import :: c_int, c_double
      integer(c_int), value :: i
      real(c_double) :: r
    end function c_species_mass
    function c_species_charge(i) bind(C, name="krome_species_charge") &
         result(r)
      import :: c_int
      integer(c_int), value :: i
      integer(c_int) :: r
    end function c_species_charge
    function c_species_natoms(i, e) bind(C, name="krome_species_natoms") &
         result(r)
      import :: c_int
      integer(c_int), value :: i, e
      integer(c_int) :: r
    end function c_species_natoms
    ! queries: (ncell, n, Tgas, par, stride, Phi, h, out...)
    function c_get_rates(ncell, n, tgas, par, stride, phi, h, k) &
         bind(C, name="krome_get_rates") result(r)
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par, phi, h
      real(c_double) :: k(*)
      integer(c_int) :: r
    end function c_get_rates
    function c_get_tdust(ncell, n, tgas, par, stride, phi, h, td) &
         bind(C, name="krome_get_Tdust") result(r)
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par, phi, h
      real(c_double) :: td(*)
      integer(c_int) :: r
    end function c_get_tdust
    function c_get_fluxes(ncell, n, tgas, par, stride, phi, h, k) &
         bind(C, name="krome_get_fluxes") result(r)
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par, phi, h
      real(c_double) :: k(*)
      integer(c_int) :: r
    end function c_get_fluxes
    function c_get_thermal(ncell, n, tgas, par, stride, phi, h, heat, cool, &
         heats, cools, dtdt) bind(C, name="krome_get_thermal") result(r)
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par, phi, h, heat, cool, heats, cools, dtdt
      integer(c_int) :: r
    end function c_get_thermal
    function c_get_gamma(ncell, n, tgas, par, stride, g) &
         bind(C, name="krome_get_gamma") result(r)
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par
      real(c_double) :: g(*)
      integer(c_int) :: r
    end function c_get_gamma
    ! cell-wise scalars (ncell, n, out)
    function c_get_qeff(ncell, n, tgas, par, stride, q) &
         bind(C, name="krome_get_qeff") result(r)
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par
      real(c_double) :: q(*)
      integer(c_int) :: r
    end function c_get_qeff
    function c_get_mu(ncell, n, out) bind(C, name="krome_get_mu") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: n(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_mu
    function c_get_rho(ncell, n, out) bind(C, name="krome_get_rho") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: n(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_rho
    function c_get_ntot(ncell, n, out) bind(C, name="krome_get_ntot") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: n(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_ntot
    function c_get_hnuclei(ncell, n, out) bind(C, name="krome_get_Hnuclei") &
         result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: n(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_hnuclei
    function c_get_charge(ncell, n, out) bind(C, name="krome_get_charge") &
         result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: n(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_charge
    function c_get_electrons(ncell, n, out) &
         bind(C, name="krome_get_electrons") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: n(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_electrons
    function c_get_tff_rho(ncell, rho, out) &
         bind(C, name="krome_get_free_fall_time_rho") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: rho(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_tff_rho
    function c_num2col(ncell, num, n, tgas, par, stride, out) &
         bind(C, name="krome_num2col") result(r)
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: num(*), n(*), tgas(*)
      type(c_ptr), value :: par
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_num2col
    function c_ncoolant() bind(C, name="krome_ncoolant") result(r)
      import :: c_int
      integer(c_int) :: r
    end function c_ncoolant
    function c_coolant_name(i) bind(C, name="krome_coolant_name") result(r)
      import :: c_int, c_ptr
      integer(c_int), value :: i
      type(c_ptr) :: r
    end function c_coolant_name
    function c_coolant_species(i) bind(C, name="krome_coolant_species") &
         result(r)
      import :: c_int
      integer(c_int), value :: i
      integer(c_int) :: r
    end function c_coolant_species
    function c_coolant_nlev(i) bind(C, name="krome_coolant_nlev") result(r)
      import :: c_int
      integer(c_int), value :: i
      integer(c_int) :: r
    end function c_coolant_nlev
    function c_coolant_levels(i, e, g) bind(C, name="krome_coolant_levels") &
         result(r)
      import :: c_int, c_double
      integer(c_int), value :: i
      real(c_double) :: e(*), g(*)
      integer(c_int) :: r
    end function c_coolant_levels
    function c_get_coolant_cooling(ncell, n, tgas, par, stride, phi, h, &
         cool) bind(C, name="krome_get_coolant_cooling") result(r)
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par, phi, h
      real(c_double) :: cool(*)
      integer(c_int) :: r
    end function c_get_coolant_cooling
    function c_get_coolant_pops(ncell, n, tgas, par, stride, phi, h, i, &
         pop) bind(C, name="krome_get_coolant_populations") result(r)
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ncell, stride, i
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par, phi, h
      real(c_double) :: pop(*)
      integer(c_int) :: r
    end function c_get_coolant_pops
    function c_get_photo_heating(ncell, n, tgas, par, stride, phi, h, &
         heat) bind(C, name="krome_get_photo_heating") result(r)
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par, phi, h
      real(c_double) :: heat(*)
      integer(c_int) :: r
    end function c_get_photo_heating
    function c_get_ratio_draine(ncell, phi, h, ratio) &
         bind(C, name="krome_get_ratio_flux_draine") result(r)
      import :: c_ptr, c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: phi(*)
      type(c_ptr), value :: h
      real(c_double) :: ratio(*)
      integer(c_int) :: r
    end function c_get_ratio_draine
    function c_get_tff(ncell, n, out) &
         bind(C, name="krome_get_free_fall_time") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: n(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_tff
    function c_get_element_totals(ncell, n, out) &
         bind(C, name="krome_get_element_totals") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: n(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_element_totals
    function c_get_element_totals_x(ncell, x, out) &
         bind(C, name="krome_get_element_totals_x") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: x(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_element_totals_x
    function c_conserve_cells(ncell, n, ref, ref_stride, par, par_stride, &
         flags) bind(C, name="krome_conserve_cells") result(r)
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: ncell, ref_stride, par_stride
      real(c_double) :: n(*)
      real(c_double), intent(in) :: ref(*)
      type(c_ptr), value :: par, flags
      integer(c_int) :: r
    end function c_conserve_cells
    function c_conserve_cells_x(ncell, x, ref, ref_stride, par, &
         par_stride, flags) bind(C, name="krome_conserve_cells_x") result(r)
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: ncell, ref_stride, par_stride
      real(c_double) :: x(*)
      real(c_double), intent(in) :: ref(*)
      type(c_ptr), value :: par, flags
      integer(c_int) :: r
    end function c_conserve_cells_x
    function c_scale_z(ncell, n, z) bind(C, name="krome_scale_Z") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double) :: n(*)
      real(c_double), value :: z
      integer(c_int) :: r
    end function c_scale_z
    function c_n2x(ncell, n, out) bind(C, name="krome_n2x") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: n(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_n2x
    function c_get_jeans(ncell, n, tgas, out) &
         bind(C, name="krome_get_jeans_length") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: n(*), tgas(*)
      real(c_double) :: out(*)
      integer(c_int) :: r
    end function c_get_jeans
    function c_x2n(ncell, x, rho, n) bind(C, name="krome_x2n") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: x(*), rho(*)
      real(c_double) :: n(*)
      integer(c_int) :: r
    end function c_x2n
    function c_equilibrium(ncell, n, tgas, par, stride, ex, flags, h) &
         bind(C, name="krome_equilibrium_cells") result(r)
      import :: c_ptr, c_int, c_double, krome_cells_c
      integer(c_int), value :: ncell, stride
      real(c_double) :: n(*), tgas(*)
      type(krome_cells_c), intent(in) :: ex
      type(c_ptr), value :: par, flags, h
      integer(c_int) :: r
    end function c_equilibrium
    ! radiation
    function c_bins_lin(nb, lo, up, elo, ehi) bind(C, name="krome_bins_lin") &
         result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb
      real(c_double), value :: lo, up
      real(c_double) :: elo(*), ehi(*)
      integer(c_int) :: r
    end function c_bins_lin
    function c_bins_log(nb, lo, up, elo, ehi) bind(C, name="krome_bins_log") &
         result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb
      real(c_double), value :: lo, up
      real(c_double) :: elo(*), ehi(*)
      integer(c_int) :: r
    end function c_bins_log
    function c_J_to_phi(nb, elo, ehi, j, phi) bind(C, name="krome_J_to_phi") &
         result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb
      real(c_double), intent(in) :: elo(*), ehi(*), j(*)
      real(c_double) :: phi(*)
      integer(c_int) :: r
    end function c_J_to_phi
    function c_phi_to_J(nb, elo, ehi, phi, j) bind(C, name="krome_phi_to_J") &
         result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb
      real(c_double), intent(in) :: elo(*), ehi(*), phi(*)
      real(c_double) :: j(*)
      integer(c_int) :: r
    end function c_phi_to_J
    function c_spectrum_bb(nb, elo, ehi, tbb, phi) &
         bind(C, name="krome_spectrum_blackbody") result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb
      real(c_double), intent(in) :: elo(*), ehi(*)
      real(c_double), value :: tbb
      real(c_double) :: phi(*)
      integer(c_int) :: r
    end function c_spectrum_bb
    function c_bb_range(tbb, lo, up) bind(C, name="krome_blackbody_range") &
         result(r)
      import :: c_int, c_double
      real(c_double), value :: tbb
      real(c_double) :: lo, up
      integer(c_int) :: r
    end function c_bb_range
    function c_spectrum_draine(nb, elo, ehi, phi) &
         bind(C, name="krome_spectrum_draine") result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb
      real(c_double), intent(in) :: elo(*), ehi(*)
      real(c_double) :: phi(*)
      integer(c_int) :: r
    end function c_spectrum_draine
    function c_spectrum_j21(nb, elo, ehi, j21, alpha, phi) &
         bind(C, name="krome_spectrum_j21") result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb
      real(c_double), intent(in) :: elo(*), ehi(*)
      real(c_double), value :: j21, alpha
      real(c_double) :: phi(*)
      integer(c_int) :: r
    end function c_spectrum_j21
    function c_spectrum_hm12(nb, elo, ehi, z, phi) &
         bind(C, name="krome_spectrum_hm2012") result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb
      real(c_double), intent(in) :: elo(*), ehi(*)
      real(c_double), value :: z
      real(c_double) :: phi(*)
      integer(c_int) :: r
    end function c_spectrum_hm12
    function c_spectrum_table(nb, elo, ehi, ne, e, j, loglog, phi) &
         bind(C, name="krome_spectrum_table") result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb, ne, loglog
      real(c_double), intent(in) :: elo(*), ehi(*), e(*), j(*)
      real(c_double) :: phi(*)
      integer(c_int) :: r
    end function c_spectrum_table
    function c_dust_kappa_bins(ne, e, k, nb, elo, ehi, kb) &
         bind(C, name="krome_dust_kappa_bins") result(r)
      import :: c_int, c_double
      integer(c_int), value :: ne, nb
      real(c_double), intent(in) :: e(*), k(*), elo(*), ehi(*)
      real(c_double) :: kb(*)
      integer(c_int) :: r
    end function c_dust_kappa_bins
    function c_dust_kappa_draine(nb, elo, ehi, kb) &
         bind(C, name="krome_dust_kappa_draine") result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb
      real(c_double), intent(in) :: elo(*), ehi(*)
      real(c_double) :: kb(*)
      integer(c_int) :: r
    end function c_dust_kappa_draine
    function c_find_G0_Av(nb, elo, ehi, phi, kb, d2g, xh, g0, av) &
         bind(C, name="krome_find_G0_Av") result(r)
      import :: c_int, c_double
      integer(c_int), value :: nb
      real(c_double), intent(in) :: elo(*), ehi(*), phi(*), kb(*)
      real(c_double), value :: d2g, xh
      real(c_double) :: g0, av
      integer(c_int) :: r
    end function c_find_G0_Av
    function krome_nphoto_c() bind(C, name="krome_nphoto") result(r)
      import :: c_int
      integer(c_int) :: r
    end function krome_nphoto_c
    function c_photo_reaction(i) bind(C, name="krome_photo_reaction") &
         result(r)
      import :: c_int
      integer(c_int), value :: i
      integer(c_int) :: r
    end function c_photo_reaction
    function c_photo_reactant(i) bind(C, name="krome_photo_reactant") &
         result(r)
      import :: c_int
      integer(c_int), value :: i
      integer(c_int) :: r
    end function c_photo_reactant
    function c_photo_xsec(i, e) bind(C, name="krome_photo_xsec") result(r)
      import :: c_int, c_double
      integer(c_int), value :: i
      real(c_double), value :: e
      real(c_double) :: r
    end function c_photo_xsec
    function c_photo_sigma(h, sig) bind(C, name="krome_opts_photo_sigma") &
         result(r)
      import :: c_int, c_double, c_ptr
      type(c_ptr), value :: h
      real(c_double) :: sig(*)
      integer(c_int) :: r
    end function c_photo_sigma
    function c_get_habing(ncell, phi, h, g) bind(C, name="krome_get_habing") &
         result(r)
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: ncell
      real(c_double), intent(in) :: phi(*)
      type(c_ptr), value :: h
      real(c_double) :: g(*)
      integer(c_int) :: r
    end function c_get_habing
    function c_get_opacity(ncell, n, tgas, par, stride, h, d2g, kdust, tau) &
         bind(C, name="krome_get_opacity") result(r)
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par, h, d2g, kdust
      real(c_double) :: tau(*)
      integer(c_int) :: r
    end function c_get_opacity
    function c_get_absorption(ncell, n, tgas, par, stride, h, kdust, kabs) &
         bind(C, name="krome_get_absorption") result(r)
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: ncell, stride
      real(c_double), intent(in) :: n(*), tgas(*)
      type(c_ptr), value :: par, h, kdust
      real(c_double) :: kabs(*)
      integer(c_int) :: r
    end function c_get_absorption
    function c_strlen(s) bind(C, name="strlen") result(r)
      import :: c_ptr, c_size_t
      type(c_ptr), value :: s
      integer(c_size_t) :: r
    end function c_strlen
  end interface

contains

  ! NUL-terminated copy of a (trimmed) Fortran string
  pure function cstr(s) result(r)
    character(len=*), intent(in) :: s
    character(kind=c_char, len=len_trim(s)+1) :: r
    r = trim(s)//c_null_char
  end function cstr

  ! set an option by name (see krome.h for the keys); stops on an unknown
  ! key unless ierr is present
  subroutine krome_set(h, key, v, ierr)
    type(c_ptr), intent(in) :: h
    character(len=*), intent(in) :: key
    real(c_double), intent(in) :: v
    integer, intent(out), optional :: ierr
    integer :: r
    r = c_opts_set(h, cstr(key), v)
    if (present(ierr)) then
      ierr = r
    else if (r /= 0) then
      print *, "krome_set: unknown or invalid option ", trim(key)
      stop 1
    end if
  end subroutine krome_set

  function krome_get(h, key) result(v)
    type(c_ptr), intent(in) :: h
    character(len=*), intent(in) :: key
    real(c_double) :: v
    if (c_opts_get(h, cstr(key), v) /= 0) v = -huge(1d0)
  end function krome_get

  ! 1-based species / parameter index by name (0 if unknown)
  function krome_species_idx(name) result(i)
    character(len=*), intent(in) :: name
    integer :: i
    i = c_species_index(cstr(name)) + 1
  end function krome_species_idx

  function krome_param_idx(name) result(i)
    character(len=*), intent(in) :: name
    integer :: i
    i = c_param_index(cstr(name)) + 1
  end function krome_param_idx

  ! Fortran copy of a C string (blank if NULL)
  function fstr(p) result(s)
    type(c_ptr), intent(in) :: p
    character(len=128) :: s
    character(kind=c_char), pointer :: ch(:)
    integer :: k, l
    s = ""
    if (.not. c_associated(p)) return
    l = int(c_strlen(p))
    call c_f_pointer(p, ch, [l])
    do k = 1, min(l, len(s))
      s(k:k) = ch(k)
    end do
  end function fstr

  ! the KROME-X version (e.g. 2.0.0-dev+g04bb53f)
  function krome_version_f() result(s)
    character(len=64) :: s
    s = trim(fstr(c_version()))
  end function krome_version_f

  ! names by 1-based index
  function krome_species_name_f(i) result(s)
    integer, intent(in) :: i
    character(len=32) :: s
    s = trim(fstr(c_species_name(int(i - 1, c_int))))
  end function krome_species_name_f

  function krome_param_name_f(i) result(s)
    integer, intent(in) :: i
    character(len=32) :: s
    s = trim(fstr(c_param_name(int(i - 1, c_int))))
  end function krome_param_name_f

  ! what requires a parameter without a default (blank: it has one)
  function krome_param_required_by_f(i) result(s)
    integer, intent(in) :: i
    character(len=128) :: s
    s = trim(fstr(c_param_required_by(int(i - 1, c_int))))
  end function krome_param_required_by_f

  function krome_reaction_name_f(i) result(s)
    integer, intent(in) :: i
    character(len=128) :: s
    s = trim(fstr(c_reaction_name(int(i - 1, c_int))))
  end function krome_reaction_name_f

  function krome_element_name_f(i) result(s)
    integer, intent(in) :: i
    character(len=16) :: s
    s = trim(fstr(c_element_name(int(i - 1, c_int))))
  end function krome_element_name_f

  function krome_heating_name_f(i) result(s)
    integer, intent(in) :: i
    character(len=32) :: s
    s = trim(fstr(c_heating_name(int(i - 1, c_int))))
  end function krome_heating_name_f

  function krome_cooling_name_f(i) result(s)
    integer, intent(in) :: i
    character(len=32) :: s
    s = trim(fstr(c_cooling_name(int(i - 1, c_int))))
  end function krome_cooling_name_f

  ! species mass [g], charge, atoms of element e (1-based indices)
  function krome_species_mass_f(i) result(m)
    integer, intent(in) :: i
    real(c_double) :: m
    m = c_species_mass(int(i - 1, c_int))
  end function krome_species_mass_f

  function krome_species_charge_f(i) result(q)
    integer, intent(in) :: i
    integer :: q
    q = c_species_charge(int(i - 1, c_int))
  end function krome_species_charge_f

  function krome_species_natoms_f(i, e) result(a)
    integer, intent(in) :: i, e
    integer :: a
    a = c_species_natoms(int(i - 1, c_int), int(e - 1, c_int))
  end function krome_species_natoms_f

  ! fill one parameter row (or every column of par(npar, ncell)) with the
  ! network defaults
  subroutine krome_params_default(par)
    real(c_double), intent(out) :: par(:, :)
    integer :: c
    do c = 1, size(par, 2)
      call c_params_default(par(:, c))
    end do
  end subroutine krome_params_default

  subroutine krome_set_photo(h, elo, ehi)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in) :: elo(:), ehi(:)
    if (c_set_photo(h, int(size(elo), c_int), elo, ehi) /= 0) then
      print *, "krome_set_photo: invalid photon bins"
      stop 1
    end if
  end subroutine krome_set_photo

  ! Doppler broadening of the tabulated cross sections, fixed at bin setup
  ! (legacy -useBroadening): Tgas [K], vturb [cm/s]; rebuilds set bins
  subroutine krome_set_photo_broadening(h, Tgas, vturb, ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in) :: Tgas, vturb
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = c_set_photo_broadening(h, Tgas, vturb)
    if (present(ierr)) ierr = r
  end subroutine krome_set_photo_broadening

  function krome_nbins(h) result(nb)
    type(c_ptr), intent(in) :: h
    integer :: nb
    nb = c_nbins(h)
  end function krome_nbins

  subroutine krome_set_dustbins(h, nd, use_si, use_c, alow, aup, phi)
    type(c_ptr), intent(in) :: h
    integer, intent(in) :: nd
    logical, intent(in) :: use_si, use_c
    real(c_double), intent(in) :: alow, aup, phi
    integer(c_int) :: si, cc
    si = merge(1, 0, use_si)
    cc = merge(1, 0, use_c)
    if (c_set_dustbins(h, int(nd, c_int), si, cc, alow, aup, phi) /= 0) &
         then
      print *, "krome_set_dustbins: invalid dust population"
      stop 1
    end if
  end subroutine krome_set_dustbins

  ! integrate size(Tgas) cells. par: (npar) broadcast or (npar, ncell) per
  ! cell (pass as a rank-2 array; a single column broadcasts); Phi:
  ! (nbins, ncell); flags: per-cell solver flags. Returns 0 on success.
  function krome_solve(n, Tgas, dt, h, par, phi, flags, td, adust) &
       result(ierr)
    real(c_double), intent(inout), target :: n(:, :)
    real(c_double), intent(inout), target :: Tgas(:)
    real(c_double), intent(in) :: dt(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional, target :: par(:, :)
    real(c_double), intent(inout), optional, target :: phi(:, :)
    integer(c_int), intent(out), optional, target :: flags(:)
    ! per-cell dust temperatures td(ndust, ncell) [K]: read and updated
    real(c_double), intent(inout), optional, target :: td(:, :)
    ! per-cell grain radii adust(ndust, ncell) [cm] (dust.evolve)
    real(c_double), intent(inout), optional, target :: adust(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi, pflags, ptd, pad
    integer(c_int) :: stride
    ppar = c_null_ptr
    pphi = c_null_ptr
    pflags = c_null_ptr
    stride = 0
    if (present(par)) then
      ppar = c_loc(par)
      if (size(par, 2) > 1) stride = int(size(par, 1), c_int)
    end if
    if (present(phi)) pphi = c_loc(phi)
    if (present(flags)) pflags = c_loc(flags)
    ptd = c_null_ptr
    if (present(td)) ptd = c_loc(td)
    pad = c_null_ptr
    if (present(adust)) pad = c_loc(adust)
    ierr = c_solve_cells(int(size(Tgas), c_int), n, Tgas, dt, ppar, stride, &
         krome_cells_c(pphi, ptd, pad), pflags, h)
  end function krome_solve

  ! krome_solve with mass fractions x(nspec, ncell) and the gas density
  ! rho(ncell) [g cm^-3] (legacy -useX): n = rho x / m in, x = m n / rho
  ! out (the caller's rho; renormalize_x rescales to the input sum) for the
  ! cells that succeed
  function krome_solve_x(x, rho, Tgas, dt, h, par, phi, flags, td, adust) &
       result(ierr)
    real(c_double), intent(inout), target :: x(:, :)
    real(c_double), intent(in) :: rho(:)
    real(c_double), intent(inout), target :: Tgas(:)
    real(c_double), intent(in) :: dt(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional, target :: par(:, :)
    real(c_double), intent(inout), optional, target :: phi(:, :)
    integer(c_int), intent(out), optional, target :: flags(:)
    real(c_double), intent(inout), optional, target :: td(:, :)
    real(c_double), intent(inout), optional, target :: adust(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi, pflags, ptd, pad
    integer(c_int) :: stride
    ppar = c_null_ptr
    pphi = c_null_ptr
    pflags = c_null_ptr
    stride = 0
    if (present(par)) then
      ppar = c_loc(par)
      if (size(par, 2) > 1) stride = int(size(par, 1), c_int)
    end if
    if (present(phi)) pphi = c_loc(phi)
    if (present(flags)) pflags = c_loc(flags)
    ptd = c_null_ptr
    if (present(td)) ptd = c_loc(td)
    pad = c_null_ptr
    if (present(adust)) pad = c_loc(adust)
    ierr = c_solve_cells_x(int(size(Tgas), c_int), x, rho, Tgas, dt, ppar, &
         stride, krome_cells_c(pphi, ptd, pad), pflags, h)
  end function krome_solve_x

  ! the error-report file (default "KROME_ERROR_REPORT"; in an MPI run
  ! give each rank its own name)
  ! a one-line description of a solver flag (0, a KROME flag -100..-111,
  ! or a negative CVODE flag)
  function krome_flag_message_f(flag) result(s)
    integer, intent(in) :: flag
    character(len=128) :: s
    s = trim(fstr(c_flag_message(int(flag, c_int))))
  end function krome_flag_message_f

  ! the dust bins' temperatures td(ndust) [K]: every cell's starting Td
  ! without a per-cell td, and the Td of the queries (C++ DustBins::Td)
  function krome_ndust(h) result(nd)
    type(c_ptr), intent(in) :: h
    integer :: nd
    nd = c_ndust(h)
  end function krome_ndust

  subroutine krome_set_dust_Td(h, td, ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in) :: td(:)
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = -1
    if (size(td) == c_ndust(h)) r = c_set_dust_td(h, td)
    if (present(ierr)) ierr = r
  end subroutine krome_set_dust_Td

  subroutine krome_get_dust_Td(h, td, ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: td(:)
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = -1
    if (size(td) == c_ndust(h)) r = c_get_dust_td(h, td)
    if (present(ierr)) ierr = r
  end subroutine krome_get_dust_Td

  ! the bins' grain radii a(ndust) [cm] (dust.evolve: every cell's
  ! starting radii without a per-cell adust), the grain counts nd(ndust)
  ! [cm^-3] of a gas density rho [g cm^-3] and d2g, the bins' grain types
  ! (1 carbon, 0 silicate)
  ! per-species tolerances (legacy -customATOL / -customRTOL): tol(nspec)
  ! replaces the scalar atol_n / rtol for the species, a negative entry
  ! takes the scalar; without tol, the scalar for every species. The get
  ! variants return the tolerance each species gets. ierr -1: wrong size,
  ! a NaN, or a zero atol
  subroutine krome_set_atol_species(h, tol, ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional :: tol(:)
    integer, intent(out), optional :: ierr
    real(c_double), allocatable, target :: buf(:)
    integer(c_int) :: r
    r = -1
    if (.not. present(tol)) then
      r = c_set_atol_sp(h, c_null_ptr)
    else if (size(tol) == krome_nspec_c()) then
      buf = tol
      r = c_set_atol_sp(h, c_loc(buf))
    end if
    if (present(ierr)) ierr = r
  end subroutine krome_set_atol_species

  subroutine krome_set_rtol_species(h, tol, ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional :: tol(:)
    integer, intent(out), optional :: ierr
    real(c_double), allocatable, target :: buf(:)
    integer(c_int) :: r
    r = -1
    if (.not. present(tol)) then
      r = c_set_rtol_sp(h, c_null_ptr)
    else if (size(tol) == krome_nspec_c()) then
      buf = tol
      r = c_set_rtol_sp(h, c_loc(buf))
    end if
    if (present(ierr)) ierr = r
  end subroutine krome_set_rtol_species

  subroutine krome_get_atol_species(h, tol, ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: tol(:)
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = -1
    if (size(tol) == krome_nspec_c()) r = c_get_atol_sp(h, tol)
    if (present(ierr)) ierr = r
  end subroutine krome_get_atol_species

  subroutine krome_get_rtol_species(h, tol, ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: tol(:)
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = -1
    if (size(tol) == krome_nspec_c()) r = c_get_rtol_sp(h, tol)
    if (present(ierr)) ierr = r
  end subroutine krome_get_rtol_species

  subroutine krome_set_dust_a(h, a, ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in) :: a(:)
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = -1
    if (size(a) == c_ndust(h)) r = c_set_dust_a(h, a)
    if (present(ierr)) ierr = r
  end subroutine krome_set_dust_a

  subroutine krome_get_dust_a(h, a, ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: a(:)
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = -1
    if (size(a) == c_ndust(h)) r = c_get_dust_a(h, a)
    if (present(ierr)) ierr = r
  end subroutine krome_get_dust_a

  subroutine krome_get_dust_counts(h, rho, d2g, nd, ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in) :: rho, d2g
    real(c_double), intent(out) :: nd(:)
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = -1
    if (size(nd) == c_ndust(h)) r = c_get_dust_counts(h, rho, d2g, nd)
    if (present(ierr)) ierr = r
  end subroutine krome_get_dust_counts

  subroutine krome_get_dust_types(h, isc, ierr)
    type(c_ptr), intent(in) :: h
    integer(c_int), intent(out) :: isc(:)
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = -1
    if (size(isc) == c_ndust(h)) r = c_get_dust_types(h, isc)
    if (present(ierr)) ierr = r
  end subroutine krome_get_dust_types

  subroutine krome_set_error_file(h, path, ierr)
    type(c_ptr), intent(in) :: h
    character(len=*), intent(in) :: path
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = c_set_error_file(h, trim(path) // c_null_char)
    if (present(ierr)) ierr = r
  end subroutine krome_set_error_file

  ! ---- queries: same conventions as krome_solve (par broadcast or per
  ! cell, Phi(nbins, ncell) with photon bins); return 0 on success
  subroutine par_phi_(par, phi, ppar, stride, pphi)
    real(c_double), intent(in), optional, target :: par(:, :), phi(:, :)
    type(c_ptr), intent(out) :: ppar, pphi
    integer(c_int), intent(out) :: stride
    ppar = c_null_ptr
    pphi = c_null_ptr
    stride = 0
    if (present(par)) then
      ppar = c_loc(par)
      if (size(par, 2) > 1) stride = int(size(par, 1), c_int)
    end if
    if (present(phi)) pphi = c_loc(phi)
  end subroutine par_phi_

  ! rate coefficients k(nrea, ncell)
  function krome_get_rates(n, Tgas, h, k, par, phi) result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: k(:, :)
    real(c_double), intent(in), optional, target :: par(:, :), phi(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi
    integer(c_int) :: stride
    call par_phi_(par, phi, ppar, stride, pphi)
    ierr = c_get_rates(int(size(Tgas), c_int), n, Tgas, ppar, stride, pphi, &
         h, k)
  end function krome_get_rates

  ! reaction fluxes k * prod(n), fl(nrea, ncell) [cm^-3 s^-1]
  function krome_get_fluxes(n, Tgas, h, fl, par, phi) result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: fl(:, :)
    real(c_double), intent(in), optional, target :: par(:, :), phi(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi
    integer(c_int) :: stride
    call par_phi_(par, phi, ppar, stride, pphi)
    ierr = c_get_fluxes(int(size(Tgas), c_int), n, Tgas, ppar, stride, &
         pphi, h, fl)
  end function krome_get_fluxes

  ! the dust temperature the rates see, td(ncell) [K] (legacy
  ! krome_get_table_Tdust); -1 when the network has no dust temperature
  function krome_get_Tdust(n, Tgas, h, td, par, phi) result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: td(:)
    real(c_double), intent(in), optional, target :: par(:, :), phi(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi
    integer(c_int) :: stride
    call par_phi_(par, phi, ppar, stride, pphi)
    ierr = c_get_tdust(int(size(Tgas), c_int), n, Tgas, ppar, stride, &
         pphi, h, td)
  end function krome_get_Tdust

  ! thermal balance [erg cm^-3 s^-1]: totals heat(ncell), cool(ncell),
  ! terms heats(nheat, ncell), cools(ncool, ncell) (krome_heating_name_f /
  ! krome_cooling_name_f), dTdt(ncell) [K/s]; all outputs optional
  function krome_get_thermal(n, Tgas, h, heat, cool, heats, cools, dTdt, &
       par, phi) result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out), optional, target :: heat(:), cool(:), &
         dTdt(:), heats(:, :), cools(:, :)
    real(c_double), intent(in), optional, target :: par(:, :), phi(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi, ph, pc, phs, pcs, pdt
    integer(c_int) :: stride
    call par_phi_(par, phi, ppar, stride, pphi)
    ph = c_null_ptr; pc = c_null_ptr; phs = c_null_ptr
    pcs = c_null_ptr; pdt = c_null_ptr
    if (present(heat)) ph = c_loc(heat)
    if (present(cool)) pc = c_loc(cool)
    if (present(heats)) phs = c_loc(heats)
    if (present(cools)) pcs = c_loc(cools)
    if (present(dTdt)) pdt = c_loc(dTdt)
    ierr = c_get_thermal(int(size(Tgas), c_int), n, Tgas, ppar, stride, &
         pphi, h, ph, pc, phs, pcs, pdt)
  end function krome_get_thermal

  function krome_get_gamma(n, Tgas, g, par) result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    real(c_double), intent(out) :: g(:)
    real(c_double), intent(in), optional, target :: par(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi
    integer(c_int) :: stride
    call par_phi_(par=par, ppar=ppar, stride=stride, pphi=pphi)
    ierr = c_get_gamma(int(size(Tgas), c_int), n, Tgas, ppar, stride, g)
  end function krome_get_gamma

  ! energy released per reaction event q(nrea, ncell) [MeV] (; qeff=)
  function krome_get_qeff(n, Tgas, q, par) result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    real(c_double), intent(out) :: q(:, :)
    real(c_double), intent(in), optional, target :: par(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi
    integer(c_int) :: stride
    call par_phi_(par=par, ppar=ppar, stride=stride, pphi=pphi)
    ierr = c_get_qeff(int(size(Tgas), c_int), n, Tgas, ppar, stride, q)
  end function krome_get_qeff

  ! per-cell scalars from n(nspec, ncell) into out(ncell)
  function krome_get_mu(n, out) result(ierr)
    real(c_double), intent(in) :: n(:, :)
    real(c_double), intent(out) :: out(:)
    integer :: ierr
    ierr = c_get_mu(int(size(n, 2), c_int), n, out)
  end function krome_get_mu

  function krome_get_rho(n, out) result(ierr)
    real(c_double), intent(in) :: n(:, :)
    real(c_double), intent(out) :: out(:)
    integer :: ierr
    ierr = c_get_rho(int(size(n, 2), c_int), n, out)
  end function krome_get_rho

  function krome_get_ntot(n, out) result(ierr)
    real(c_double), intent(in) :: n(:, :)
    real(c_double), intent(out) :: out(:)
    integer :: ierr
    ierr = c_get_ntot(int(size(n, 2), c_int), n, out)
  end function krome_get_ntot

  function krome_get_Hnuclei(n, out) result(ierr)
    real(c_double), intent(in) :: n(:, :)
    real(c_double), intent(out) :: out(:)
    integer :: ierr
    ierr = c_get_hnuclei(int(size(n, 2), c_int), n, out)
  end function krome_get_Hnuclei

  function krome_get_charge(n, out) result(ierr)
    real(c_double), intent(in) :: n(:, :)
    real(c_double), intent(out) :: out(:)
    integer :: ierr
    ierr = c_get_charge(int(size(n, 2), c_int), n, out)
  end function krome_get_charge

  function krome_get_electrons(n, out) result(ierr)
    real(c_double), intent(in) :: n(:, :)
    real(c_double), intent(out) :: out(:)
    integer :: ierr
    ierr = c_get_electrons(int(size(n, 2), c_int), n, out)
  end function krome_get_electrons

  function krome_get_free_fall_time(n, out) result(ierr)
    real(c_double), intent(in) :: n(:, :)
    real(c_double), intent(out) :: out(:)
    integer :: ierr
    ierr = c_get_tff(int(size(n, 2), c_int), n, out)
  end function krome_get_free_fall_time

  ! the free-fall time [s] from the mass density rho(ncell) [g cm^-3]
  function krome_get_free_fall_time_rho(rho, out) result(ierr)
    real(c_double), intent(in) :: rho(:)
    real(c_double), intent(out) :: out(:)
    integer :: ierr
    ierr = c_get_tff_rho(int(size(rho), c_int), rho, out)
  end function krome_get_free_fall_time_rho

  ! the column density [cm^-2] of num(ncell) [cm^-3] in each cell
  ! (radiation.num2col; legacy krome_num2col)
  function krome_num2col(num, n, Tgas, out, par) result(ierr)
    real(c_double), intent(in) :: num(:), n(:, :), Tgas(:)
    real(c_double), intent(out) :: out(:)
    real(c_double), intent(in), optional, target :: par(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi
    integer(c_int) :: stride
    call par_phi_(par, ppar=ppar, stride=stride, pphi=pphi)
    ierr = c_num2col(int(size(Tgas), c_int), num, n, Tgas, ppar, stride, out)
  end function krome_num2col

  ! the line coolants (the coolZ.dat metals, then the thermal.coolant_files ones),
  ! 1-based: their number, legacy name, species, levels (energies [K],
  ! degeneracies)
  function krome_ncoolant_c() result(r)
    integer :: r
    r = c_ncoolant()
  end function krome_ncoolant_c

  function krome_coolant_name_f(i) result(s)
    integer, intent(in) :: i
    character(len=32) :: s
    s = trim(fstr(c_coolant_name(int(i - 1, c_int))))
  end function krome_coolant_name_f

  function krome_coolant_species_f(i) result(r)
    integer, intent(in) :: i
    integer :: r
    r = c_coolant_species(int(i - 1, c_int)) + 1
  end function krome_coolant_species_f

  function krome_coolant_nlev_f(i) result(r)
    integer, intent(in) :: i
    integer :: r
    r = c_coolant_nlev(int(i - 1, c_int))
  end function krome_coolant_nlev_f

  subroutine krome_coolant_levels(i, E, g, ierr)
    integer, intent(in) :: i
    real(c_double), intent(out) :: E(:), g(:)
    integer, intent(out), optional :: ierr
    integer(c_int) :: r
    r = -1
    if (size(E) == krome_coolant_nlev_f(i) .and. size(g) == size(E)) &
         r = c_coolant_levels(int(i - 1, c_int), E, g)
    if (present(ierr)) ierr = r
  end subroutine krome_coolant_levels

  ! each coolant's line cooling cool(ncoolant, ncell) [erg cm^-3 s^-1]
  function krome_get_coolant_cooling(n, Tgas, h, cool, par, phi) &
       result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: cool(:, :)
    real(c_double), intent(in), optional, target :: par(:, :), phi(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi
    integer(c_int) :: stride
    call par_phi_(par, phi, ppar, stride, pphi)
    ierr = c_get_coolant_cooling(int(size(Tgas), c_int), n, Tgas, ppar, &
         stride, pphi, h, cool)
  end function krome_get_coolant_cooling

  ! the level populations of coolant i (1-based), pop(nlev_i, ncell)
  ! [cm^-3] (legacy krome_popcool_dump)
  function krome_get_coolant_populations(n, Tgas, h, i, pop, par, phi) &
       result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    type(c_ptr), intent(in) :: h
    integer, intent(in) :: i
    real(c_double), intent(out) :: pop(:, :)
    real(c_double), intent(in), optional, target :: par(:, :), phi(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi
    integer(c_int) :: stride
    ierr = -1
    if (size(pop, 1) /= krome_coolant_nlev_f(i)) return
    call par_phi_(par, phi, ppar, stride, pphi)
    ierr = c_get_coolant_pops(int(size(Tgas), c_int), n, Tgas, ppar, &
         stride, pphi, h, int(i - 1, c_int), pop)
  end function krome_get_coolant_populations

  ! photoheating per photo reaction heat(nphoto, ncell) [erg cm^-3 s^-1]
  ! (legacy krome_get_photoBin_heats, per absorber)
  function krome_get_photo_heating(n, Tgas, h, heat, par, phi) result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: heat(:, :)
    real(c_double), intent(in), optional, target :: par(:, :), phi(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi
    integer(c_int) :: stride
    call par_phi_(par, phi, ppar, stride, pphi)
    ierr = c_get_photo_heating(int(size(Tgas), c_int), n, Tgas, ppar, &
         stride, pphi, h, heat)
  end function krome_get_photo_heating

  ! the field over the Draine ISRF, ratio(ncell) of phi(nbins, ncell)
  ! (legacy krome_get_ratioFluxDraine)
  function krome_get_ratio_flux_draine(phi, h, ratio) result(ierr)
    real(c_double), intent(in) :: phi(:, :)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: ratio(:)
    integer :: ierr
    ierr = c_get_ratio_draine(int(size(phi, 2), c_int), phi, h, ratio)
  end function krome_get_ratio_flux_draine

  function krome_get_jeans_length(n, Tgas, out) result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    real(c_double), intent(out) :: out(:)
    integer :: ierr
    ierr = c_get_jeans(int(size(n, 2), c_int), n, Tgas, out)
  end function krome_get_jeans_length

  ! element totals tot(nelem, ncell) [atoms cm^-3]; _x: per gram of gas
  ! [atoms g^-1] from the mass fractions x(nspec, ncell) (legacy
  ! krome_conserveLinGetRef_x)
  function krome_get_element_totals(n, tot) result(ierr)
    real(c_double), intent(in) :: n(:, :)
    real(c_double), intent(out) :: tot(:, :)
    integer :: ierr
    ierr = c_get_element_totals(int(size(n, 2), c_int), n, tot)
  end function krome_get_element_totals
  function krome_get_element_totals_x(x, tot) result(ierr)
    real(c_double), intent(in) :: x(:, :)
    real(c_double), intent(out) :: tot(:, :)
    integer :: ierr
    ierr = c_get_element_totals_x(int(size(x, 2), c_int), x, tot)
  end function krome_get_element_totals_x

  ! consistent advection (legacy krome_conserveLin_x; manual 6.7): the
  ! cells n(nspec, ncell) [or x, mass fractions] moved back to the
  ! reference element totals ref(nelem, ncell) or ref(nelem, 1) (one row
  ! for all cells) and to zero net charge; par(npar, ncell) or par(npar,
  ! 1), optional (the defaults), for the freeze directives; flags(ncell),
  ! optional: 0, krome_flag_conserve_unreached or the warning
  ! krome_flag_conserve_large. Returns the lowest failing cell's flag, or 0
  function krome_conserve(n, ref, par, flags) result(ierr)
    real(c_double), intent(inout) :: n(:, :)
    real(c_double), intent(in) :: ref(:, :)
    real(c_double), intent(in), target, optional :: par(:, :)
    integer(c_int), intent(out), target, optional :: flags(:)
    integer :: ierr
    ierr = conserve_f_(n, ref, .false., par, flags)
  end function krome_conserve
  function krome_conserve_x(x, ref, par, flags) result(ierr)
    real(c_double), intent(inout) :: x(:, :)
    real(c_double), intent(in) :: ref(:, :)
    real(c_double), intent(in), target, optional :: par(:, :)
    integer(c_int), intent(out), target, optional :: flags(:)
    integer :: ierr
    ierr = conserve_f_(x, ref, .true., par, flags)
  end function krome_conserve_x
  function conserve_f_(v, ref, isx, par, flags) result(ierr)
    real(c_double), intent(inout) :: v(:, :)
    real(c_double), intent(in) :: ref(:, :)
    logical, intent(in) :: isx
    real(c_double), intent(in), target, optional :: par(:, :)
    integer(c_int), intent(out), target, optional :: flags(:)
    integer :: ierr
    integer(c_int) :: rs, ps
    type(c_ptr) :: pp, fp
    rs = 0
    if (size(ref, 2) > 1) rs = int(size(ref, 1), c_int)
    pp = c_null_ptr
    ps = 0
    if (present(par)) then
      pp = c_loc(par)
      if (size(par, 2) > 1) ps = int(size(par, 1), c_int)
    end if
    fp = c_null_ptr
    if (present(flags)) fp = c_loc(flags)
    if (isx) then
      ierr = c_conserve_cells_x(int(size(v, 2), c_int), v, ref, rs, pp, &
                                ps, fp)
    else
      ierr = c_conserve_cells(int(size(v, 2), c_int), v, ref, rs, pp, ps, &
                              fp)
    end if
  end function conserve_f_

  ! number densities <-> mass fractions (rho(ncell) [g cm^-3])
  function krome_n2x(n, x) result(ierr)
    real(c_double), intent(in) :: n(:, :)
    real(c_double), intent(out) :: x(:, :)
    integer :: ierr
    ierr = c_n2x(int(size(n, 2), c_int), n, x)
  end function krome_n2x

  function krome_x2n(x, rho, n) result(ierr)
    real(c_double), intent(in) :: x(:, :), rho(:)
    real(c_double), intent(out) :: n(:, :)
    integer :: ierr
    ierr = c_x2n(int(size(x, 2), c_int), x, rho, n)
  end function krome_x2n

  ! ---- legacy conveniences (krome_user): one cell n(nspec); par(npar)
  ! and phi(nbins) optional as in the queries; h may be c_null_ptr

  ! the species names separated by spaces (legacy krome_get_names_header)
  function krome_names_header() result(s)
    character(len=:), allocatable :: s
    integer :: i
    s = ""
    do i = 1, krome_nspec_c()
      if (i > 1) s = s // " "
      s = s // trim(krome_species_name_f(i))
    end do
  end function krome_names_header

  ! one cell's reaction fluxes fl(nrea) [cm^-3 s^-1]
  subroutine cell_fluxes_(n, Tgas, h, par, phi, fl, ierr)
    real(c_double), intent(in) :: n(:), Tgas
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional :: par(:), phi(:)
    real(c_double), intent(out) :: fl(:)
    integer, intent(out) :: ierr
    real(c_double) :: nn(size(n), 1), T(1), f2(size(fl), 1)
    real(c_double), allocatable :: p2(:, :), ph2(:, :)
    nn(:, 1) = n
    T(1) = Tgas
    ! an unallocated array is an absent optional argument
    if (present(par)) p2 = reshape(par, [size(par), 1])
    if (present(phi)) ph2 = reshape(phi, [size(phi), 1])
    ierr = krome_get_fluxes(nn, T, h, f2, par=p2, phi=ph2)
    fl = f2(:, 1)
  end subroutine cell_fluxes_

  ! the nbest largest reaction fluxes [cm^-3 s^-1], index and equation
  ! (legacy krome_print_best_flux); ierr: 0 or a query flag
  subroutine krome_print_best_flux(n, Tgas, nbest, h, par, phi, ierr)
    real(c_double), intent(in) :: n(:), Tgas
    integer, intent(in) :: nbest
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional :: par(:), phi(:)
    integer, intent(out), optional :: ierr
    real(c_double), allocatable :: fl(:)
    logical, allocatable :: used(:)
    integer :: i, j, r
    allocate(fl(krome_nrea_c()), used(krome_nrea_c()))
    call cell_fluxes_(n, Tgas, h, par, phi, fl, r)
    if (present(ierr)) ierr = r
    if (r /= 0) return
    used = .false.
    print *, "***************"
    do i = 1, min(nbest, size(fl))
      j = maxloc(fl, 1, mask=.not. used)
      used(j) = .true.
      print '(I8,a1,a50,E17.8)', j, " ", krome_reaction_name_f(j), fl(j)
    end do
  end subroutine krome_print_best_flux

  ! every reaction's flux, absolute and relative to the largest and to the
  ! sum, one line per reaction then a blank line, to unit, labelled by xvar
  ! (legacy krome_explore_flux)
  subroutine krome_explore_flux(n, Tgas, unit, xvar, h, par, phi, ierr)
    real(c_double), intent(in) :: n(:), Tgas, xvar
    integer, intent(in) :: unit
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional :: par(:), phi(:)
    integer, intent(out), optional :: ierr
    real(c_double), allocatable :: fl(:)
    real(c_double) :: fmax, fsum
    integer :: i, r
    allocate(fl(krome_nrea_c()))
    call cell_fluxes_(n, Tgas, h, par, phi, fl, r)
    if (present(ierr)) ierr = r
    if (r /= 0) return
    fmax = maxval(fl)
    fsum = sum(fl)
    do i = 1, size(fl)
      write(unit, '(I8,5E17.8e3,a3,a50)') i, xvar, Tgas, fl(i), &
           fl(i) / fmax, fl(i) / fsum, " ", krome_reaction_name_f(i)
    end do
    write(unit, *)
  end subroutine krome_explore_flux

  ! every coolant's level populations to unit, one row per level: name,
  ! level (1-based), Tgas, population, the coolant's total; a blank line
  ! at the end (legacy krome_popcool_dump)
  subroutine krome_print_coolant_populations(n, Tgas, unit, h, par, phi, &
       ierr)
    real(c_double), intent(in) :: n(:), Tgas
    integer, intent(in) :: unit
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional :: par(:), phi(:)
    integer, intent(out), optional :: ierr
    real(c_double) :: nn(size(n), 1), T(1)
    real(c_double), allocatable :: p2(:, :), ph2(:, :), pop(:, :)
    integer :: i, l, r
    nn(:, 1) = n
    T(1) = Tgas
    if (present(par)) p2 = reshape(par, [size(par), 1])
    if (present(phi)) ph2 = reshape(phi, [size(phi), 1])
    r = 0
    do i = 1, krome_ncoolant_c()
      if (allocated(pop)) deallocate(pop)
      allocate(pop(krome_coolant_nlev_f(i), 1))
      r = krome_get_coolant_populations(nn, T, h, i, pop, par=p2, phi=ph2)
      if (r /= 0) exit
      do l = 1, size(pop, 1)
        write(unit, '(a8,I5,3E17.8e3)') trim(krome_coolant_name_f(i)), l, &
             Tgas, pop(l, 1), sum(pop(:, 1))
      end do
    end do
    if (r == 0) write(unit, *)
    if (present(ierr)) ierr = r
  end subroutine krome_print_coolant_populations

  ! one row to unit: Tgas, the total cooling, every cooling term (legacy
  ! krome_dump_cooling)
  subroutine krome_dump_cooling(n, Tgas, unit, h, par, phi, ierr)
    real(c_double), intent(in) :: n(:), Tgas
    integer, intent(in) :: unit
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional :: par(:), phi(:)
    integer, intent(out), optional :: ierr
    real(c_double) :: nn(size(n), 1), T(1), cool(1)
    real(c_double), allocatable :: cools(:, :)
    real(c_double), allocatable :: p2(:, :), ph2(:, :)
    integer :: r
    allocate(cools(max(krome_ncool_c(), 1), 1))
    nn(:, 1) = n
    T(1) = Tgas
    if (present(par)) p2 = reshape(par, [size(par), 1])
    if (present(phi)) ph2 = reshape(phi, [size(phi), 1])
    r = krome_get_thermal(nn, T, h, par=p2, phi=ph2, cool=cool, cools=cools)
    if (r == 0) write(unit, '(99E14.5e3)') Tgas, cool(1), &
         cools(1:krome_ncool_c(), 1)
    if (present(ierr)) ierr = r
  end subroutine krome_dump_cooling

  ! the rate coefficients on npts log temperatures Tmin..Tmax at n =
  ! 1e-40, rows Tgas, reaction, k to unit (legacy krome_dump_rates)
  subroutine krome_dump_rates(Tmin, Tmax, npts, unit, h, par, phi, ierr)
    real(c_double), intent(in) :: Tmin, Tmax
    integer, intent(in) :: npts, unit
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional :: par(:), phi(:)
    integer, intent(out), optional :: ierr
    real(c_double) :: T(1)
    real(c_double), allocatable :: nn(:, :), k(:, :)
    real(c_double), allocatable :: p2(:, :), ph2(:, :)
    integer :: i, j, r
    allocate(nn(krome_nspec_c(), 1), k(max(krome_nrea_c(), 1), 1))
    r = -1
    if (Tmin > 0d0 .and. Tmax >= Tmin .and. npts >= 1) r = 0
    nn = 1d-40
    if (present(par)) p2 = reshape(par, [size(par), 1])
    if (present(phi)) ph2 = reshape(phi, [size(phi), 1])
    do i = 1, npts
      if (r /= 0) exit
      T(1) = 1d1**log10(Tmin)
      if (npts > 1) T(1) = 1d1**(log10(Tmin) + (i - 1) &
           * (log10(Tmax) - log10(Tmin)) / (npts - 1))
      r = krome_get_rates(nn, T, h, k, par=p2, phi=ph2)
      if (r /= 0) exit
      do j = 1, krome_nrea_c()
        write(unit, '(E17.8e3,I8,E17.8e3)') T(1), j, k(j, 1)
      end do
      write(unit, *)
    end do
    if (present(ierr)) ierr = r
  end subroutine krome_dump_rates

  ! every parameter of the row par (absent: the defaults), 'name: value'
  ! (legacy krome_print_phys_variables)
  subroutine krome_print_params(par)
    real(c_double), intent(in), optional :: par(:)
    real(c_double), allocatable :: p(:, :)
    integer :: i
    allocate(p(krome_npar_c(), 1))
    call krome_params_default(p)
    if (present(par)) p(:, 1) = par
    do i = 1, size(p, 1)
      print '(a,": ",ES24.16)', trim(krome_param_name_f(i)), p(i, 1)
    end do
  end subroutine krome_print_params

  ! the species table: index, name, mass [g], charge, abundance, their sum
  ! and Tgas (legacy krome_get_info)
  subroutine krome_print_info(n, Tgas)
    real(c_double), intent(in) :: n(:), Tgas
    integer :: i
    print '(a4,a10,a11,a5,a11)', "#", "Name", "m (g)", "Chrg", "x"
    do i = 1, size(n)
      print '(I4,a10,E11.3,I5,E11.3)', i, " " // krome_species_name_f(i), &
           krome_species_mass_f(i), krome_species_charge_f(i), n(i)
    end do
    print '(a30,E11.3)', " sum", sum(n)
    print '(a14,E11.3)', "Tgas", Tgas
  end subroutine krome_print_info

  ! the neutral atomic metals at the Asplund+09 abundances times 10^Z
  ! (legacy krome_scale_Z; generated with the network)
  subroutine scale_z_1_(n, Z)
    real(c_double), intent(inout) :: n(:)
    real(c_double), intent(in) :: Z
    integer(c_int) :: r
    r = c_scale_z(1_c_int, n, Z)
  end subroutine scale_z_1_

  subroutine scale_z_2_(n, Z)
    real(c_double), intent(inout) :: n(:, :)
    real(c_double), intent(in) :: Z
    integer(c_int) :: r
    r = c_scale_z(int(size(n, 2), c_int), n, Z)
  end subroutine scale_z_2_

  ! integrate every cell to steady state (see krome.h); flags(ncell): 0
  ! converged, 1 eq_tmax reached first, <0 failure (cell untouched)
  function krome_equilibrium(n, Tgas, h, par, phi, flags, td, adust) &
       result(ierr)
    real(c_double), intent(inout), target :: n(:, :), Tgas(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(in), optional, target :: par(:, :)
    real(c_double), intent(inout), optional, target :: phi(:, :)
    integer(c_int), intent(out), optional, target :: flags(:)
    real(c_double), intent(inout), optional, target :: td(:, :)
    real(c_double), intent(inout), optional, target :: adust(:, :)
    integer :: ierr
    type(c_ptr) :: ppar, pphi, pflags, ptd, pad
    integer(c_int) :: stride
    call par_phi_(par, phi, ppar, stride, pphi)
    pflags = c_null_ptr
    if (present(flags)) pflags = c_loc(flags)
    ptd = c_null_ptr
    if (present(td)) ptd = c_loc(td)
    pad = c_null_ptr
    if (present(adust)) pad = c_loc(adust)
    ierr = c_equilibrium(int(size(Tgas), c_int), n, Tgas, ppar, stride, &
         krome_cells_c(pphi, ptd, pad), pflags, h)
  end function krome_equilibrium

  ! ---- radiation: Phi(nbins[, ncell]) photon flux per bin; every spectrum
  ! is the bin integral 4 pi/h int J dlnE (J [eV s^-1 cm^-2 Hz^-1 sr^-1]);
  ! functions return 0, or -1 for invalid arguments (see krome.h)
  function krome_bins_lin(lower, upper, elo, ehi) result(ierr)
    real(c_double), intent(in) :: lower, upper
    real(c_double), intent(out) :: elo(:), ehi(:)
    integer :: ierr
    ierr = c_bins_lin(int(size(elo), c_int), lower, upper, elo, ehi)
  end function krome_bins_lin

  function krome_bins_log(lower, upper, elo, ehi) result(ierr)
    real(c_double), intent(in) :: lower, upper
    real(c_double), intent(out) :: elo(:), ehi(:)
    integer :: ierr
    ierr = c_bins_log(int(size(elo), c_int), lower, upper, elo, ehi)
  end function krome_bins_log

  function krome_J_to_phi(elo, ehi, J, phi) result(ierr)
    real(c_double), intent(in) :: elo(:), ehi(:), J(:)
    real(c_double), intent(out) :: phi(:)
    integer :: ierr
    ierr = c_J_to_phi(int(size(elo), c_int), elo, ehi, J, phi)
  end function krome_J_to_phi

  function krome_phi_to_J(elo, ehi, phi, J) result(ierr)
    real(c_double), intent(in) :: elo(:), ehi(:), phi(:)
    real(c_double), intent(out) :: J(:)
    integer :: ierr
    ierr = c_phi_to_J(int(size(elo), c_int), elo, ehi, phi, J)
  end function krome_phi_to_J

  function krome_spectrum_blackbody(elo, ehi, Tbb, phi) result(ierr)
    real(c_double), intent(in) :: elo(:), ehi(:), Tbb
    real(c_double), intent(out) :: phi(:)
    integer :: ierr
    ierr = c_spectrum_bb(int(size(elo), c_int), elo, ehi, Tbb, phi)
  end function krome_spectrum_blackbody

  function krome_blackbody_range(Tbb, lower, upper) result(ierr)
    real(c_double), intent(in) :: Tbb
    real(c_double), intent(out) :: lower, upper
    integer :: ierr
    ierr = c_bb_range(Tbb, lower, upper)
  end function krome_blackbody_range

  function krome_spectrum_draine(elo, ehi, phi) result(ierr)
    real(c_double), intent(in) :: elo(:), ehi(:)
    real(c_double), intent(out) :: phi(:)
    integer :: ierr
    ierr = c_spectrum_draine(int(size(elo), c_int), elo, ehi, phi)
  end function krome_spectrum_draine

  ! J21 and alpha default to the legacy 1, 1
  function krome_spectrum_j21(elo, ehi, phi, J21, alpha) result(ierr)
    real(c_double), intent(in) :: elo(:), ehi(:)
    real(c_double), intent(out) :: phi(:)
    real(c_double), intent(in), optional :: J21, alpha
    integer :: ierr
    real(c_double) :: j, a
    j = 1d0
    a = 1d0
    if (present(J21)) j = J21
    if (present(alpha)) a = alpha
    ierr = c_spectrum_j21(int(size(elo), c_int), elo, ehi, j, a, phi)
  end function krome_spectrum_j21

  function krome_spectrum_hm2012(elo, ehi, z, phi) result(ierr)
    real(c_double), intent(in) :: elo(:), ehi(:), z
    real(c_double), intent(out) :: phi(:)
    integer :: ierr
    ierr = c_spectrum_hm12(int(size(elo), c_int), elo, ehi, z, phi)
  end function krome_spectrum_hm2012

  ! tabulated J(E) at ascending E [eV]; loglog interpolation when true
  function krome_spectrum_table(elo, ehi, E, J, phi, loglog) result(ierr)
    real(c_double), intent(in) :: elo(:), ehi(:), E(:), J(:)
    real(c_double), intent(out) :: phi(:)
    logical, intent(in), optional :: loglog
    integer :: ierr
    integer(c_int) :: ll
    ll = 0
    if (present(loglog)) ll = merge(1, 0, loglog)
    ierr = c_spectrum_table(int(size(elo), c_int), elo, ehi, &
         int(size(E), c_int), E, J, ll, phi)
  end function krome_spectrum_table

  ! bin-mean dust opacity kbin(nbins) [cm^2 g^-1] from kappa(E), E [eV]
  ! ascending, or from the embedded Draine R_V = 3.1 table
  function krome_dust_kappa_bins(E, kappa, elo, ehi, kbin) result(ierr)
    real(c_double), intent(in) :: E(:), kappa(:), elo(:), ehi(:)
    real(c_double), intent(out) :: kbin(:)
    integer :: ierr
    ierr = c_dust_kappa_bins(int(size(E), c_int), E, kappa, &
         int(size(elo), c_int), elo, ehi, kbin)
  end function krome_dust_kappa_bins

  function krome_dust_kappa_draine(elo, ehi, kbin) result(ierr)
    real(c_double), intent(in) :: elo(:), ehi(:)
    real(c_double), intent(out) :: kbin(:)
    integer :: ierr
    ierr = c_dust_kappa_draine(int(size(elo), c_int), elo, ehi, kbin)
  end function krome_dust_kappa_draine

  ! G0 and Av of a field (legacy krome_find_G0_Av); xH = Hnuclei m_p / rho
  function krome_find_G0_Av(elo, ehi, phi, kbin, d2g, xH, G0, Av) &
       result(ierr)
    real(c_double), intent(in) :: elo(:), ehi(:), phi(:), kbin(:), d2g, xH
    real(c_double), intent(out) :: G0, Av
    integer :: ierr
    ierr = c_find_G0_Av(int(size(elo), c_int), elo, ehi, phi, kbin, d2g, xH, &
         G0, Av)
  end function krome_find_G0_Av

  ! photo reaction i (1-based): its reaction and reactant (1-based), and
  ! cross section [cm^2] at E [eV]
  function krome_photo_reaction_f(i) result(r)
    integer, intent(in) :: i
    integer :: r
    r = c_photo_reaction(int(i - 1, c_int)) + 1
  end function krome_photo_reaction_f

  function krome_photo_reactant_f(i) result(r)
    integer, intent(in) :: i
    integer :: r
    r = c_photo_reactant(int(i - 1, c_int)) + 1
  end function krome_photo_reactant_f

  function krome_photo_xsec_f(i, E) result(s)
    integer, intent(in) :: i
    real(c_double), intent(in) :: E
    real(c_double) :: s
    s = c_photo_xsec(int(i - 1, c_int), E)
  end function krome_photo_xsec_f

  ! bin-averaged photo cross sections sigma(nbins, nphoto) of h's bins
  function krome_photo_sigma(h, sigma) result(ierr)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: sigma(:, :)
    integer :: ierr
    ierr = c_photo_sigma(h, sigma)
  end function krome_photo_sigma

  ! Habing flux G(ncell) of phi(nbins, ncell)
  function krome_get_habing(phi, h, G) result(ierr)
    real(c_double), intent(in) :: phi(:, :)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: G(:)
    integer :: ierr
    ierr = c_get_habing(int(size(phi, 2), c_int), phi, h, G)
  end function krome_get_habing

  ! optical depth tau(nbins, ncell) with the num2col columns (legacy
  ! krome_get_opacity); dust with d2g(ncell): the dust bins of h, else
  ! kdust(nbins) [cm^2 per g of dust], default Draine R_V = 3.1
  function krome_get_opacity(n, Tgas, h, tau, par, d2g, kdust) result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: tau(:, :)
    real(c_double), intent(in), optional, target :: par(:, :), d2g(:), &
         kdust(:)
    integer :: ierr
    type(c_ptr) :: ppar, pphi, pd, pk
    integer(c_int) :: stride
    call par_phi_(par=par, ppar=ppar, stride=stride, pphi=pphi)
    pd = c_null_ptr
    pk = c_null_ptr
    if (present(d2g)) pd = c_loc(d2g)
    if (present(kdust)) pk = c_loc(kdust)
    ierr = c_get_opacity(int(size(Tgas), c_int), n, Tgas, ppar, stride, h, &
         pd, pk, tau)
  end function krome_get_opacity

  ! absorption coefficient kabs(nbins, ncell) [cm^-1]; kdust(nbins): dust
  ! opacity per unit GAS mass [cm^2 g^-1] (optional)
  function krome_get_absorption(n, Tgas, h, kabs, par, kdust) result(ierr)
    real(c_double), intent(in) :: n(:, :), Tgas(:)
    type(c_ptr), intent(in) :: h
    real(c_double), intent(out) :: kabs(:, :)
    real(c_double), intent(in), optional, target :: par(:, :), kdust(:)
    integer :: ierr
    type(c_ptr) :: ppar, pphi, pk
    integer(c_int) :: stride
    call par_phi_(par=par, ppar=ppar, stride=stride, pphi=pphi)
    pk = c_null_ptr
    if (present(kdust)) pk = c_loc(kdust)
    ierr = c_get_absorption(int(size(Tgas), c_int), n, Tgas, ppar, stride, &
         h, pk, kabs)
  end function krome_get_absorption

end module krome

! physical constants (kromex/constants.py: CODATA 2018, IAU; the legacy
! module name krome_constants), the numbers the solver uses:
!   use krome_constants
module krome_constants
  implicit none
  public
  real(8), parameter :: boltzmann_eV = 8.617333262145179d-05  ! eV / K; exact
  real(8), parameter :: boltzmann_J = 1.380649d-23  ! J / K; exact
  real(8), parameter :: boltzmann_erg = 1.380649d-16  ! erg / K; exact
  real(8), parameter :: planck_eV = 4.135667696923859d-15  ! eV s; exact
  real(8), parameter :: planck_J = 6.62607015d-34  ! J s; exact
  real(8), parameter :: planck_erg = 6.62607015d-27  ! erg s; exact
  real(8), parameter :: hbar_erg = 1.0545718176461565d-27  ! erg s; exact
  real(8), parameter :: hc_k = 1.4387768775039336d0  ! cm K; second radiation constant: E/k [K] = hc_k x E [cm^-1]
  real(8), parameter :: gravity = 6.6743d-08  ! cm^3 g^-1 s^-2; CODATA 2018
  real(8), parameter :: e_mass = 9.1093837015d-28  ! g; CODATA 2018
  real(8), parameter :: p_mass = 1.67262192369d-24  ! g; CODATA 2018
  real(8), parameter :: n_mass = 1.67492749804d-24  ! g; CODATA 2018
  real(8), parameter :: amu = 1.6605390666d-24  ! g; CODATA 2018
  real(8), parameter :: clight = 29979245800.0d0  ! cm s^-1; exact
  real(8), parameter :: pi = 3.141592653589793d0
  real(8), parameter :: eV_to_erg = 1.602176634d-12  ! eV -> erg; exact
  real(8), parameter :: ry_to_eV = 13.605693122994d0  ! Rydberg -> eV; CODATA 2018
  real(8), parameter :: ry_to_erg = 2.1798723611035473d-11  ! Rydberg -> erg; CODATA 2018
  real(8), parameter :: seconds_per_year = 31557600.0d0  ! yr -> s; Julian year, 365.25 d
  real(8), parameter :: km_to_cm = 100000.0d0  ! km -> cm
  real(8), parameter :: pc_to_cm = 3.085677581491368d+18  ! pc -> cm; IAU 2015
  real(8), parameter :: Mpc_to_cm = 3.085677581491368d+24  ! Mpc -> cm; IAU 2015
  real(8), parameter :: stefboltz_erg = 5.6703744191844286d-05  ! erg s^-1 cm^-2 K^-4; exact
  real(8), parameter :: N_avogadro = 6.02214076d+23  ! mol^-1; exact
  real(8), parameter :: Rgas_J = 8.31446261815324d0  ! J K^-1 mol^-1; exact
  real(8), parameter :: Rgas_kJ = 0.00831446261815324d0  ! kJ K^-1 mol^-1; exact
end module krome_constants
