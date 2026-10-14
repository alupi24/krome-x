!################################################################
! The Fortran-interface test: the one-zone primordial collapse of
! tests/collapse through the Fortran binding (same loop, parameters and
! output as test_cpp.cpp, which uses the C++ API; check_diffs.py compares
! the two outputs)
!################################################################
program test_krome
  use krome            ! the binding: krome_solve, krome_opts_*, queries,
                       ! krome_names_header, krome_print_best_flux, ...
  use krome_net        ! krome_nspec, krome_idx_*, krome_par_*
  use krome_test_utils ! krome_free_fall_time, krome_check
  implicit none
  integer, parameter :: rstep = 500000
  integer :: i, unit, out
  real(8) :: dtH, tff, dd, dd1, ntot
  real(8) :: x(krome_nspec, 1), Tgas(1), dt(1), par(krome_npar, 1)
  type(c_ptr) :: h

  !INITIAL CONDITIONS
  ntot = 1d0               !total density, cm-3
  Tgas = 1d2               !temperature, K

  !solver options: thermal evolution on
  h = krome_opts_new()
  call krome_set(h, "evolve_T", 1d0)

  !per-cell parameters: the network defaults, then the legacy default
  ! Tcmb = 2.73 K (the legacy driver's krome_redshift = 15 is an unused
  ! user variable: legacy ran at z = 0); tff at every step
  call krome_params_default(par)
  par(krome_par_Tcmb, 1) = 2.73d0

  !species default, cm-3
  x(:, 1) = 1d-40
  x(krome_idx_H, 1)   = ntot
  x(krome_idx_H2, 1)  = 1d-6 * ntot
  x(krome_idx_E, 1)   = 1d-4 * ntot
  x(krome_idx_Hj, 1)  = 1d-4 * ntot
  x(krome_idx_He, 1)  = 0.0775d0 * ntot

  dd = ntot
  open(newunit=unit, file="explore.dat", status="replace")
  open(newunit=out, file="abundances.dat", status="replace")
  write(out, '(a)') "#ntot Tgas " // trim(krome_names_header())
  print '(a5,2a11)', "step", "n(cm-3)", "Tgas(K)"

  do i = 1, rstep
     dd1 = dd
     !free-fall time, s: the compression timescale of COMPRESS heating
     tff = krome_free_fall_time(x(:, 1))
     par(krome_par_tff, 1) = tff
     dtH = 0.01d0 * tff
     dd = dd + dd / tff * dtH
     x(:, 1) = x(:, 1) * (dd / dd1)
     dt = dtH
     if (dd > 1d18) exit
     call krome_check(krome_solve(x, Tgas, dt, h, par=par), "solve")
     write(out, '(99E17.8e3)') dd, Tgas(1), x(:, 1) / dd
     if (mod(i, 100) == 0) then
        print '(I5,99E11.3)', i, dd, Tgas(1)
        call krome_print_best_flux(x(:, 1), Tgas(1), 5, h, par(:, 1))
        call krome_explore_flux(x(:, 1), Tgas(1), unit, dd, h, par(:, 1))
     end if
  end do
  close(unit)
  close(out)
  call krome_opts_free(h)
  print *, "Now run ./test_cpp, then python3 check_diffs.py"

end program test_krome
