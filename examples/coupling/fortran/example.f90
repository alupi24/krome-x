!################################################################
! KROME-X from Fortran (the binding over the C interface): the cell of
! the quick start (README), primordial gas at nH = 1e4 cm-3 and 100 K,
! followed in time to 1e6 yr; the C, C++ and Python examples solve the
! same cell (README.md)
!################################################################
program example
  use krome            ! the binding: krome_solve, krome_opts_*, queries
  use krome_net        ! krome_nspec, krome_idx_*, krome_par_*
  implicit none
  real(8), parameter :: spy = 3.15d7   ! seconds per year
  real(8), parameter :: nH = 1d4       ! cm-3
  integer :: out, i, ierr
  real(8) :: t, x(krome_nspec, 1), Tgas(1), dt(1), par(krome_npar, 1)
  type(c_ptr) :: h

  !solver options: Tgas evolved with the chemistry
  h = krome_opts_new()
  call krome_set(h, "evolve_T", 1d0)

  !per-cell parameters: the network defaults, and the free-fall time of
  ! the COMPRESS heating
  call krome_params_default(par)
  par(krome_par_tff, 1) = 1d14

  !initial densities, cm-3
  x(:, 1) = 0d0
  x(krome_idx_H, 1)  = nH
  x(krome_idx_Hj, 1) = 1d0
  x(krome_idx_E, 1)  = 1d0
  x(krome_idx_He, 1) = 775d0
  Tgas = 1d2

  open(newunit=out, file="abundances.dat", status="replace")
  write(out, '(a)') "#t_yr Tgas " // trim(krome_names_header())
  print '(a5,2a11)', "step", "t(yr)", "Tgas(K)"

  !steps from 1 yr, 10% longer each, to 1e6 yr
  t = 0d0
  dt = spy
  i = 0
  do while (t < 1d6 * spy)
     ierr = krome_solve(x, Tgas, dt, h, par=par)
     if (ierr /= 0) then
        print *, "KROME failed: ", trim(krome_flag_message_f(ierr))
        stop 1
     end if
     t = t + dt(1)
     dt = dt * 1.1d0
     i = i + 1
     write(out, '(99E17.8e3)') t / spy, Tgas(1), x(:, 1) / nH
     if (mod(i, 20) == 0) print '(I5,99E11.3)', i, t / spy, Tgas(1)
  end do
  close(out)
  call krome_opts_free(h)
  print *, "wrote abundances.dat"

end program example
