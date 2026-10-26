"""KROME-X Python interface (optional binding, codegen project.bindings = ['python']).

A ctypes binding of the C interface (krome.h) in the shared library
libkrome.so / libkrome.dylib that `make shared` builds in the generated
project directory. Needs numpy; nothing else (no compilation step of its
own).

    import krome
    net = krome.load()                 # libkrome next to this file, or
                                       # $KROME_LIB, or load("path/lib")
    n = net.densities({"H": 1e4, "H+": 1e-4, "e-": 1e-4})
    opts = net.options(evolve_T=True, rtol=1e-8)
    res = net.solve(n, Tgas=100.0, dt=3.15e13, par={"crate": 1e-16},
                    opts=opts)          # net.params lists the parameters;
                                        # net.required those without a
                                        # default, which must be given
    print(res.Tgas, res.n[net.idx("H2")])

Array layouts are those of the C interface, with numpy's row-major
order: n[ncell, nspec], par[ncell, npar], Phi[ncell, nbins]. A single
cell may be passed as 1-D arrays and scalars, and comes back the same
way. solve() returns NEW arrays; the inputs are not modified.
"""

import ctypes
import os
import sys
from pathlib import Path

import numpy as np

__all__ = ["load", "Network", "Options", "Result", "Thermal",
           "KromeError"]

_D = ctypes.POINTER(ctypes.c_double)
_I = ctypes.POINTER(ctypes.c_int)
# query functions: (ncell, n, Tgas, par, par_stride, Phi, opts, ...)
_QARGS = [ctypes.c_int, _D, _D, _D, ctypes.c_int, _D, ctypes.c_void_p]
_CELLWISE = {"mu": "krome_get_mu", "rho": "krome_get_rho",
             "ntot": "krome_get_ntot", "hnuclei": "krome_get_Hnuclei",
             "charge": "krome_get_charge",
             "electrons": "krome_get_electrons",
             "free_fall_time": "krome_get_free_fall_time",
             "element_totals": "krome_get_element_totals",
             "n2x": "krome_n2x"}

# option values that are enumerations in the C interface
_ENUMS = {"jac": ("auto", "analytic", "fd"),
          "linsol": ("auto", "dense", "lapack", "sparse"),
          "path": ("auto", "percell", "device")}
_BOOLS = ("evolve_T", "clip_negative", "error_report", "renormalize_x",
          "device_cache")
_INTS = ("max_steps", "min_device", "device_memory", "min_interleave",
         "max_restarts", "max_err_fails", "max_conv_fails", "conserve")
_KEYS = ("evolve_T", "jac", "linsol", "rtol", "atol_n", "atol_T", "atol_u",
         "max_steps", "max_err_fails", "max_conv_fails", "max_restarts",
         "clip_negative", "error_report", "renormalize_x", "conserve",
         "path", "device_cache", "device_memory", "min_device",
         "min_interleave", "eq_t0", "eq_tmax", "eq_factor", "eq_rtol",
         "eq_atol")


class _Constants:
    """physical constants (kromex/constants.py: CODATA 2018, IAU), the
    numbers the solver uses (krome.constants.p_mass, ...)"""

    boltzmann_eV = 8.617333262145179e-05  # eV / K; exact
    boltzmann_J = 1.380649e-23  # J / K; exact
    boltzmann_erg = 1.380649e-16  # erg / K; exact
    planck_eV = 4.135667696923859e-15  # eV s; exact
    planck_J = 6.62607015e-34  # J s; exact
    planck_erg = 6.62607015e-27  # erg s; exact
    hbar_erg = 1.0545718176461565e-27  # erg s; exact
    hc_k = 1.4387768775039336  # cm K; second radiation constant: E/k [K] = hc_k x E [cm^-1]
    gravity = 6.6743e-08  # cm^3 g^-1 s^-2; CODATA 2018
    e_mass = 9.1093837015e-28  # g; CODATA 2018
    p_mass = 1.67262192369e-24  # g; CODATA 2018
    n_mass = 1.67492749804e-24  # g; CODATA 2018
    amu = 1.6605390666e-24  # g; CODATA 2018
    clight = 29979245800.0  # cm s^-1; exact
    pi = 3.141592653589793
    eV_to_erg = 1.602176634e-12  # eV -> erg; exact
    ry_to_eV = 13.605693122994  # Rydberg -> eV; CODATA 2018
    ry_to_erg = 2.1798723611035473e-11  # Rydberg -> erg; CODATA 2018
    seconds_per_year = 31557600.0  # yr -> s; Julian year, 365.25 d
    km_to_cm = 100000.0  # km -> cm
    pc_to_cm = 3.085677581491368e+18  # pc -> cm; IAU 2015
    Mpc_to_cm = 3.085677581491368e+24  # Mpc -> cm; IAU 2015
    stefboltz_erg = 5.6703744191844286e-05  # erg s^-1 cm^-2 K^-4; exact
    N_avogadro = 6.02214076e+23  # mol^-1; exact
    Rgas_J = 8.31446261815324  # J K^-1 mol^-1; exact
    Rgas_kJ = 0.00831446261815324  # kJ K^-1 mol^-1; exact

    def __repr__(self):
        return "<krome.constants: %s>" % ", ".join(
            k for k in vars(type(self)) if not k.startswith("_"))


constants = _Constants()


class _Cells(ctypes.Structure):
    """krome.h krome_cells: the optional per-cell arrays"""
    _fields_ = [("Phi", _D), ("Td", _D), ("adust", _D)]


def _cells(Phi, Td, adust):
    """a krome_cells for the given arrays (None: not used), or NULL"""
    if Phi is None and Td is None and adust is None:
        return None
    return ctypes.pointer(_Cells(_dptr(Phi), _dptr(Td), _dptr(adust)))


class KromeError(RuntimeError):
    """a rejected argument or option, or a failed solve."""


def _default_library():
    env = os.environ.get("KROME_LIB")
    if env:
        return Path(env)
    here = Path(__file__).resolve().parent
    ext = ".dylib" if sys.platform == "darwin" else ".so"
    return here / ("libkrome" + ext)


def load(path=None):
    """load the KROME library of a generated project (default: the
    libkrome next to this module, or $KROME_LIB)."""
    return Network(path)


def _dptr(a):
    return None if a is None else a.ctypes.data_as(_D)


class Network:
    """one generated network, compiled into a shared libkrome."""

    def __init__(self, path=None):
        path = Path(path) if path is not None else _default_library()
        if not path.exists():
            raise KromeError("KROME library %s not found: run `make shared` "
                             "in the project directory (or set KROME_LIB)"
                             % path)
        self.path = path
        lib = self._lib = ctypes.CDLL(str(path))
        proto = {
            "krome_network_name": (ctypes.c_char_p, []),
            "krome_nspec": (ctypes.c_int, []),
            "krome_version": (ctypes.c_char_p, []),
            "krome_nrea": (ctypes.c_int, []),
            "krome_npar": (ctypes.c_int, []),
            "krome_species_name": (ctypes.c_char_p, [ctypes.c_int]),
            "krome_param_name": (ctypes.c_char_p, [ctypes.c_int]),
            "krome_param_required_by": (ctypes.c_char_p, [ctypes.c_int]),
            "krome_params_default": (None, [_D]),
            "krome_opts_new": (ctypes.c_void_p, []),
            "krome_opts_free": (None, [ctypes.c_void_p]),
            "krome_opts_set": (ctypes.c_int, [ctypes.c_void_p,
                                              ctypes.c_char_p,
                                              ctypes.c_double]),
            "krome_opts_get": (ctypes.c_int, [ctypes.c_void_p,
                                              ctypes.c_char_p, _D]),
            "krome_opts_set_photo": (ctypes.c_int, [ctypes.c_void_p,
                                                    ctypes.c_int, _D, _D]),
            "krome_opts_nbins": (ctypes.c_int, [ctypes.c_void_p]),
            "krome_opts_set_photo_broadening": (
                ctypes.c_int, [ctypes.c_void_p, ctypes.c_double,
                               ctypes.c_double]),
            "krome_opts_set_dustbins": (
                ctypes.c_int, [ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                               ctypes.c_int, ctypes.c_double,
                               ctypes.c_double, ctypes.c_double]),
            "krome_solve_cells": (
                ctypes.c_int, [ctypes.c_int, _D, _D, _D, _D, ctypes.c_int,
                               ctypes.POINTER(_Cells), _I,
                               ctypes.c_void_p]),
            "krome_solve_cells_x": (
                ctypes.c_int, [ctypes.c_int, _D, _D, _D, _D, _D,
                               ctypes.c_int, ctypes.POINTER(_Cells), _I,
                               ctypes.c_void_p]),
            "krome_opts_set_atol_species": (ctypes.c_int,
                                            [ctypes.c_void_p, _D]),
            "krome_opts_set_rtol_species": (ctypes.c_int,
                                            [ctypes.c_void_p, _D]),
            "krome_opts_get_atol_species": (ctypes.c_int,
                                            [ctypes.c_void_p, _D]),
            "krome_opts_get_rtol_species": (ctypes.c_int,
                                            [ctypes.c_void_p, _D]),
            "krome_opts_set_dust_a": (ctypes.c_int, [ctypes.c_void_p, _D]),
            "krome_opts_get_dust_a": (ctypes.c_int, [ctypes.c_void_p, _D]),
            "krome_opts_get_dust_counts": (
                ctypes.c_int, [ctypes.c_void_p, ctypes.c_double,
                               ctypes.c_double, _D]),
            "krome_opts_get_dust_types": (ctypes.c_int, [ctypes.c_void_p,
                                                         _I]),
            "krome_opts_ndust": (ctypes.c_int, [ctypes.c_void_p]),
            "krome_opts_set_dust_Td": (ctypes.c_int, [ctypes.c_void_p, _D]),
            "krome_opts_get_dust_Td": (ctypes.c_int, [ctypes.c_void_p, _D]),
            "krome_flag_message": (ctypes.c_char_p, [ctypes.c_int]),
            "krome_opts_set_error_file": (
                ctypes.c_int, [ctypes.c_void_p, ctypes.c_char_p]),
            "krome_equilibrium_cells": (
                ctypes.c_int, [ctypes.c_int, _D, _D, _D, ctypes.c_int,
                               ctypes.POINTER(_Cells), _I,
                               ctypes.c_void_p]),
            # queries
            "krome_get_rates": (ctypes.c_int, _QARGS + [_D]),
            "krome_get_coolant_cooling": (ctypes.c_int, _QARGS + [_D]),
            "krome_get_coolant_populations": (ctypes.c_int,
                                              _QARGS + [ctypes.c_int, _D]),
            "krome_get_photo_heating": (ctypes.c_int, _QARGS + [_D]),
            "krome_get_ratio_flux_draine": (
                ctypes.c_int, [ctypes.c_int, _D, ctypes.c_void_p, _D]),
            "krome_get_free_fall_time_rho": (ctypes.c_int,
                                             [ctypes.c_int, _D, _D]),
            "krome_num2col": (ctypes.c_int, [ctypes.c_int, _D, _D, _D, _D,
                                             ctypes.c_int, _D]),
            "krome_ncoolant": (ctypes.c_int, []),
            "krome_coolant_name": (ctypes.c_char_p, [ctypes.c_int]),
            "krome_coolant_species": (ctypes.c_int, [ctypes.c_int]),
            "krome_coolant_nlev": (ctypes.c_int, [ctypes.c_int]),
            "krome_coolant_levels": (ctypes.c_int, [ctypes.c_int, _D, _D]),
            "krome_get_fluxes": (ctypes.c_int, _QARGS + [_D]),
            "krome_get_Tdust": (ctypes.c_int, _QARGS + [_D]),
            "krome_get_thermal": (ctypes.c_int, _QARGS + [_D] * 5),
            "krome_get_gamma": (ctypes.c_int, [ctypes.c_int, _D, _D, _D,
                                               ctypes.c_int, _D]),
            "krome_get_qeff": (ctypes.c_int, [ctypes.c_int, _D, _D, _D,
                                              ctypes.c_int, _D]),
            "krome_get_jeans_length": (ctypes.c_int, [ctypes.c_int, _D, _D,
                                                      _D]),
            "krome_x2n": (ctypes.c_int, [ctypes.c_int, _D, _D, _D]),
            "krome_scale_Z": (ctypes.c_int, [ctypes.c_int, _D,
                                             ctypes.c_double]),
            "krome_nheat": (ctypes.c_int, []),
            "krome_ncool": (ctypes.c_int, []),
            "krome_heating_name": (ctypes.c_char_p, [ctypes.c_int]),
            "krome_cooling_name": (ctypes.c_char_p, [ctypes.c_int]),
            "krome_species_mass": (ctypes.c_double, [ctypes.c_int]),
            "krome_species_charge": (ctypes.c_int, [ctypes.c_int]),
            "krome_nelem": (ctypes.c_int, []),
            "krome_element_name": (ctypes.c_char_p, [ctypes.c_int]),
            "krome_species_natoms": (ctypes.c_int, [ctypes.c_int,
                                                    ctypes.c_int]),
            "krome_reaction_name": (ctypes.c_char_p, [ctypes.c_int]),
            # radiation
            "krome_bins_lin": (ctypes.c_int, [ctypes.c_int, ctypes.c_double,
                                              ctypes.c_double, _D, _D]),
            "krome_bins_log": (ctypes.c_int, [ctypes.c_int, ctypes.c_double,
                                              ctypes.c_double, _D, _D]),
            "krome_J_to_phi": (ctypes.c_int, [ctypes.c_int, _D, _D, _D, _D]),
            "krome_phi_to_J": (ctypes.c_int, [ctypes.c_int, _D, _D, _D, _D]),
            "krome_spectrum_blackbody": (ctypes.c_int, [ctypes.c_int, _D, _D,
                                                        ctypes.c_double, _D]),
            "krome_blackbody_range": (ctypes.c_int, [ctypes.c_double, _D,
                                                     _D]),
            "krome_spectrum_draine": (ctypes.c_int, [ctypes.c_int, _D, _D,
                                                     _D]),
            "krome_spectrum_j21": (ctypes.c_int, [ctypes.c_int, _D, _D,
                                                  ctypes.c_double,
                                                  ctypes.c_double, _D]),
            "krome_spectrum_hm2012": (ctypes.c_int, [ctypes.c_int, _D, _D,
                                                     ctypes.c_double, _D]),
            "krome_spectrum_table": (ctypes.c_int, [ctypes.c_int, _D, _D,
                                                    ctypes.c_int, _D, _D,
                                                    ctypes.c_int, _D]),
            "krome_dust_kappa_bins": (ctypes.c_int, [ctypes.c_int, _D, _D,
                                                     ctypes.c_int, _D, _D,
                                                     _D]),
            "krome_dust_kappa_draine": (ctypes.c_int, [ctypes.c_int, _D, _D,
                                                       _D]),
            "krome_find_G0_Av": (ctypes.c_int, [ctypes.c_int, _D, _D, _D, _D,
                                                ctypes.c_double,
                                                ctypes.c_double, _D, _D]),
            "krome_nphoto": (ctypes.c_int, []),
            "krome_photo_reaction": (ctypes.c_int, [ctypes.c_int]),
            "krome_photo_reactant": (ctypes.c_int, [ctypes.c_int]),
            "krome_photo_xsec": (ctypes.c_double, [ctypes.c_int,
                                                   ctypes.c_double]),
            "krome_opts_photo_sigma": (ctypes.c_int, [ctypes.c_void_p, _D]),
            "krome_get_habing": (ctypes.c_int, [ctypes.c_int, _D,
                                                ctypes.c_void_p, _D]),
            "krome_get_opacity": (ctypes.c_int, [ctypes.c_int, _D, _D, _D,
                                                 ctypes.c_int,
                                                 ctypes.c_void_p, _D, _D,
                                                 _D]),
            "krome_photo_pdiss": (ctypes.c_double, [ctypes.c_int]),
            "krome_nlw": (ctypes.c_int, []),
            "krome_has_device": (ctypes.c_int, []),
            "krome_device_initialize": (ctypes.c_int, []),
            "krome_device_finalize": (ctypes.c_int, []),
            "krome_lw_reaction": (ctypes.c_int, [ctypes.c_int]),
            "krome_lw_reactant": (ctypes.c_int, [ctypes.c_int]),
            "krome_lw_pdiss": (ctypes.c_double, [ctypes.c_int]),
            "krome_opts_lw_weights": (ctypes.c_int, [ctypes.c_void_p, _D]),
            "krome_get_element_totals_x": (ctypes.c_int, [ctypes.c_int, _D,
                                                          _D]),
            "krome_conserve_cells": (ctypes.c_int, [ctypes.c_int, _D, _D,
                                                    ctypes.c_int, _D,
                                                    ctypes.c_int, _I]),
            "krome_conserve_cells_x": (ctypes.c_int, [ctypes.c_int, _D, _D,
                                                      ctypes.c_int, _D,
                                                      ctypes.c_int, _I]),
            "krome_get_absorption": (ctypes.c_int, [ctypes.c_int, _D, _D, _D,
                                                    ctypes.c_int,
                                                    ctypes.c_void_p, _D,
                                                    _D]),
        }
        for name in _CELLWISE.values():
            proto[name] = (ctypes.c_int, [ctypes.c_int, _D, _D])
        for name, (res, args) in proto.items():
            fn = getattr(lib, name)
            fn.restype, fn.argtypes = res, args
        self.name = lib.krome_network_name().decode()
        self.version = lib.krome_version().decode()  # KROME-X version
        self.nspec = lib.krome_nspec()
        self.nrea = lib.krome_nrea()
        self.npar = lib.krome_npar()
        self.species = tuple(lib.krome_species_name(i).decode()
                             for i in range(self.nspec))
        self.params = tuple(lib.krome_param_name(i).decode()
                            for i in range(self.npar))
        # parameters without a default: {name: what requires it}
        self.required = {
            s: lib.krome_param_required_by(i).decode()
            for i, s in enumerate(self.params)
            if lib.krome_param_required_by(i)}
        # species, element and reaction metadata
        self.masses = np.array([lib.krome_species_mass(i)
                                for i in range(self.nspec)])  # [g]
        self.charges = np.array([lib.krome_species_charge(i)
                                 for i in range(self.nspec)])
        self.elements = tuple(lib.krome_element_name(e).decode()
                              for e in range(lib.krome_nelem()))
        self.composition = np.array(
            [[lib.krome_species_natoms(i, e)
              for e in range(len(self.elements))]
             for i in range(self.nspec)], dtype=int).reshape(
                 self.nspec, len(self.elements))
        self.reactions = tuple(lib.krome_reaction_name(r).decode()
                               for r in range(self.nrea))
        self.heating_terms = tuple(lib.krome_heating_name(i).decode()
                                   for i in range(lib.krome_nheat()))
        self.cooling_terms = tuple(lib.krome_cooling_name(i).decode()
                                   for i in range(lib.krome_ncool()))
        # the line coolants (coolZ.dat metals, then the --cool-file ones):
        # legacy names, species, and levels (energies [K], degeneracies)
        nco = lib.krome_ncoolant()
        self.coolants = tuple(lib.krome_coolant_name(i).decode()
                              for i in range(nco))
        self.coolant_species = tuple(
            self.species[lib.krome_coolant_species(i)] for i in range(nco))
        self.coolant_levels = []
        for i in range(nco):
            E = np.zeros(lib.krome_coolant_nlev(i))
            g = np.zeros_like(E)
            lib.krome_coolant_levels(i, _dptr(E), _dptr(g))
            self.coolant_levels.append((E, g))
        self.coolant_levels = tuple(self.coolant_levels)
        # photo reactions: reaction index and absorbing reactant (species)
        self.photo_reactions = tuple(lib.krome_photo_reaction(r)
                                     for r in range(lib.krome_nphoto()))
        self.photo_reactants = tuple(lib.krome_photo_reactant(r)
                                     for r in range(lib.krome_nphoto()))
        self.photo_pdiss = tuple(lib.krome_photo_pdiss(r)
                                 for r in range(lib.krome_nphoto()))
        # J_LW rates: reaction, absorber, pdiss
        self.lw_reactions = tuple(lib.krome_lw_reaction(q)
                                  for q in range(lib.krome_nlw()))
        self.lw_reactants = tuple(lib.krome_lw_reactant(q)
                                  for q in range(lib.krome_nlw()))
        self.lw_pdiss = tuple(lib.krome_lw_pdiss(q)
                              for q in range(lib.krome_nlw()))
        self._sidx = {s: i for i, s in enumerate(self.species)}
        self._pidx = {s: i for i, s in enumerate(self.params)}
        self._opts = None

    def __repr__(self):
        return ("<krome.Network %r: %d species, %d reactions, %d params "
                "(%s)>" % (self.name, self.nspec, self.nrea, self.npar,
                           self.path))

    # -- description -------------------------------------------------------

    def idx(self, species):
        """0-based index of a species in n[..., nspec]."""
        try:
            return self._sidx[species]
        except KeyError:
            raise KeyError("species %r not in network %r (species: %s)"
                           % (species, self.name, ", ".join(self.species))
                           ) from None

    def par_idx(self, param):
        """0-based index of a runtime parameter in par[..., npar]."""
        try:
            return self._pidx[param]
        except KeyError:
            raise KeyError("parameter %r not in network %r (parameters: "
                           "%s)" % (param, self.name, ", ".join(self.params))
                           ) from None

    def params_default(self):
        """one row of the default parameter values, in `params` order;
        a required parameter (see `required`) is NaN until set."""
        row = np.empty(self.npar)
        self._lib.krome_params_default(_dptr(row))
        return row

    def param_rows(self, ncell=1, **values):
        """par[ncell, npar]: the defaults with the given overrides; each
        value is a scalar (every cell) or one value per cell."""
        rows = np.tile(self.params_default(), (ncell, 1))
        for name, v in values.items():
            rows[:, self.par_idx(name)] = np.broadcast_to(
                np.asarray(v, dtype=float), (ncell,))
        return rows

    def densities(self, values, ncell=None, fill=0.0):
        """number densities [cm^-3] from {species: value}: a vector
        n[nspec], or n[ncell, nspec] when ncell is given (a value may then
        be one per cell)."""
        shape = (self.nspec,) if ncell is None else (ncell, self.nspec)
        n = np.full(shape, float(fill))
        for s, v in values.items():
            n[..., self.idx(s)] = v
        return n

    # -- solver ------------------------------------------------------------

    def options(self, **values):
        """a new solver-options handle, e.g. options(evolve_T=True,
        rtol=1e-8, jac="fd")."""
        return Options(self, **values)

    @property
    def has_device(self):
        """True when the library has the device (Kokkos) backend."""
        return bool(self._lib.krome_has_device())

    def device_initialize(self):
        """initialize Kokkos for the device backend (a no-op when already
        initialized); call device_finalize() before exiting."""
        if self._lib.krome_device_initialize() != 0:
            raise KromeError("no device backend in this library")

    def device_finalize(self):
        """finalize Kokkos (after the last device solve)."""
        if self._lib.krome_device_finalize() != 0:
            raise KromeError("no device backend in this library")

    def solve(self, n, Tgas, dt, par=None, Phi=None, opts=None, check=True,
              Td=None, adust=None):
        """integrate every cell over its own dt [s].

        n     n[nspec] or n[ncell, nspec] [cm^-3]
        Tgas  scalar or one per cell [K] (evolved when opts evolve_T)
        dt    scalar or one per cell [s]
        par   None (defaults), {name: scalar or per-cell values}, one row
              par[npar] for every cell, or par[ncell, npar]
        Phi   photon flux per bin, needed when opts.set_photo() was used:
              Phi[nbins] for every cell or Phi[ncell, nbins]
        opts  an Options handle (default: all defaults)
        check raise KromeError when a cell fails (else inspect
              result.flags)
        Td    per-cell dust temperatures [K] with dust bins on opts:
              Td[ndust] for every cell or Td[ncell, ndust]; each cell
              starts from its row, which receives its final Td (None: the
              bins' Td, opts.set_dust_Td)
        adust per-cell grain radii [cm] (dust.evolve), as Td: adust[ndust]
              or adust[ncell, ndust] (None: the bins' radii, opts.dust_a)

        Returns a Result with new arrays n, Tgas, Phi (with radiation.evolve
        the evolved field), Td and adust (when given), the per-cell solver
        flags and ierr."""
        single, n, ncell = self._cells(n)
        T = self._per_cell(Tgas, ncell)
        dtv = self._per_cell(dt, ncell)
        opts = self._options(opts)
        P, stride = self._params(par, ncell)
        F = self._phi(Phi, ncell, opts)
        D = self._td(Td, ncell, opts)
        A = self._td(adust, ncell, opts, "adust")

        flags = np.zeros(ncell, dtype=np.intc)
        ierr = self._lib.krome_solve_cells(
            ncell, _dptr(n), _dptr(T), _dptr(dtv), _dptr(P), stride,
            _cells(F, D, A), flags.ctypes.data_as(_I), opts._h)
        res = Result(n[0] if single else n, T[0] if single else T,
                     None if F is None else (F[0] if single else F),
                     flags[0] if single else flags, ierr,
                     None if D is None else (D[0] if single else D),
                     None if A is None else (A[0] if single else A))
        if check and ierr != 0:
            bad = np.flatnonzero(np.atleast_1d(flags))
            raise KromeError("krome_solve_cells returned %d: %s (failed "
                             "cells: %s); pass check=False to get the "
                             "result anyway" % (ierr, self.flag_message(ierr),
                                                bad.tolist() or "n/a"))
        return res


    def solve_x(self, x, rho, Tgas, dt, par=None, Phi=None, opts=None,
                check=True, Td=None, adust=None):
        """solve() with mass fractions (legacy -useX): x[nspec] or
        x[ncell, nspec], and the gas density rho [g cm^-3], scalar or one
        per cell. n = rho x / m on the way in, x = m n / rho on the way out
        (the caller's rho; renormalize_x rescales to the input sum) for the
        cells that succeed;
        the solver itself works on number densities. Returns a Result
        whose n holds the new mass fractions."""
        single, x, ncell = self._cells(x)
        rh = self._per_cell(rho, ncell)
        T = self._per_cell(Tgas, ncell)
        dtv = self._per_cell(dt, ncell)
        opts = self._options(opts)
        P, stride = self._params(par, ncell)
        F = self._phi(Phi, ncell, opts)
        D = self._td(Td, ncell, opts)
        A = self._td(adust, ncell, opts, "adust")
        flags = np.zeros(ncell, dtype=np.intc)
        ierr = self._lib.krome_solve_cells_x(
            ncell, _dptr(x), _dptr(rh), _dptr(T), _dptr(dtv), _dptr(P),
            stride, _cells(F, D, A), flags.ctypes.data_as(_I), opts._h)
        res = Result(x[0] if single else x, T[0] if single else T,
                     None if F is None else (F[0] if single else F),
                     flags[0] if single else flags, ierr,
                     None if D is None else (D[0] if single else D),
                     None if A is None else (A[0] if single else A))
        if check and ierr != 0:
            bad = np.flatnonzero(np.atleast_1d(flags))
            raise KromeError("krome_solve_cells_x returned %d: %s (failed "
                             "cells: %s); pass check=False to get the "
                             "result anyway" % (ierr, self.flag_message(ierr),
                                                bad.tolist() or "n/a"))
        return res

    def flag_message(self, flag):
        """a one-line description of a solver flag (0, a KROME flag
        -100..-111, or a negative CVODE flag)"""
        return self._lib.krome_flag_message(int(flag)).decode()

    def equilibrium(self, n, Tgas, par=None, Phi=None, opts=None,
                    check=True, Td=None, adust=None):
        """integrate every cell to steady state (legacy krome_equilibrium;
        isothermal unless opts evolve_T; options eq_t0, eq_tmax, eq_factor,
        eq_rtol, eq_atol). Result.flags: 0 converged, 1 eq_tmax reached
        first (state returned), <0 failure (cell untouched; raises unless
        check=False)."""
        single, n, ncell = self._cells(n)
        T = self._per_cell(Tgas, ncell)
        opts = self._options(opts)
        P, stride = self._params(par, ncell)
        F = self._phi(Phi, ncell, opts)
        D = self._td(Td, ncell, opts)
        A = self._td(adust, ncell, opts, "adust")
        flags = np.zeros(ncell, dtype=np.intc)
        ierr = self._lib.krome_equilibrium_cells(
            ncell, _dptr(n), _dptr(T), _dptr(P), stride, _cells(F, D, A),
            flags.ctypes.data_as(_I), opts._h)
        res = Result(n[0] if single else n, T[0] if single else T,
                     None if F is None else (F[0] if single else F),
                     flags[0] if single else flags, ierr,
                     None if D is None else (D[0] if single else D),
                     None if A is None else (A[0] if single else A))
        if check and (np.atleast_1d(flags) < 0).any():
            bad = np.flatnonzero(np.atleast_1d(flags) < 0)
            raise KromeError("equilibrium failed in cells %s; pass "
                             "check=False to get the result anyway"
                             % bad.tolist())
        return res

    # -- input preparation (the krome_solve_cells conventions) ----------------

    def _cells(self, n):
        n = np.array(n, dtype=float, order="C", copy=True)
        if n.ndim == 0 or n.ndim > 2 or n.shape[-1] != self.nspec:
            raise KromeError("n must have shape (%d,) or (ncell, %d), got "
                             "%s" % (self.nspec, self.nspec, n.shape))
        single = n.ndim == 1
        n = n.reshape(-1, self.nspec)
        return single, n, n.shape[0]

    @staticmethod
    def _per_cell(v, ncell):
        return np.array(np.broadcast_to(np.asarray(v, dtype=float),
                                        (ncell,)), order="C", copy=True)

    def _options(self, opts):
        if opts is not None:
            return opts
        if self._opts is None:
            self._opts = Options(self)
        return self._opts

    def _params(self, par, ncell):
        """(rows or None, stride); raises for unset REQUIRED parameters"""
        stride, P = 0, None
        if isinstance(par, dict):
            P = self.param_rows(ncell, **par)
            stride = self.npar
        elif par is not None:
            P = np.ascontiguousarray(par, dtype=float)
            if P.shape == (self.npar,):
                stride = 0
            elif P.shape == (ncell, self.npar):
                stride = self.npar
            else:
                raise KromeError("par must have shape (%d,) or (%d, %d), got "
                                 "%s" % (self.npar, ncell, self.npar,
                                         P.shape))
        rows = (self.params_default()[None, :] if P is None
                else P.reshape(-1, self.npar))
        unset = [s for s in self.required
                 if np.isnan(rows[:, self._pidx[s]]).any()]
        if unset:
            raise KromeError("required parameter(s) not set: %s" % "; ".join(
                "%s (for %s)" % (s, self.required[s]) for s in unset))
        return P, stride

    @staticmethod
    def _phi(Phi, ncell, opts):
        B = opts.nbins
        if B > 0:
            if Phi is None:
                raise KromeError("photon bins are configured (set_photo): "
                                 "Phi is required")
            return np.array(np.broadcast_to(np.asarray(Phi, dtype=float),
                                            (ncell, B)), order="C", copy=True)
        if Phi is not None:
            raise KromeError("Phi given but no photon bins configured "
                             "(opts.set_photo)")
        return None

    @staticmethod
    def _td(Td, ncell, opts, what="Td"):
        if Td is None:
            return None
        nd = opts.ndust
        if nd == 0:
            raise KromeError("%s given but no dust bins configured "
                             "(opts.set_dustbins)" % what)
        return np.array(np.broadcast_to(np.asarray(Td, dtype=float),
                                        (ncell, nd)), order="C", copy=True)

    def _query(self, fname, width, n, Tgas, par, Phi, opts):
        single, n, ncell = self._cells(n)
        T = self._per_cell(Tgas, ncell)
        opts = self._options(opts)
        P, stride = self._params(par, ncell)
        F = self._phi(Phi, ncell, opts)
        out = np.zeros((ncell, width))
        rc = getattr(self._lib, fname)(ncell, _dptr(n), _dptr(T), _dptr(P),
                                       stride, _dptr(F), opts._h, _dptr(out))
        if rc != 0:
            raise KromeError("%s returned %d" % (fname, rc))
        return out[0] if single else out

    # -- queries (legacy krome_get_*) --------------------------------------

    def rates(self, n, Tgas, par=None, Phi=None, opts=None):
        """rate coefficients k[..., nrea] (labels: `reactions`)"""
        return self._query("krome_get_rates", self.nrea, n, Tgas, par, Phi,
                           opts)

    def fluxes(self, n, Tgas, par=None, Phi=None, opts=None):
        """reaction fluxes k * prod(n) [cm^-3 s^-1], [..., nrea]"""
        return self._query("krome_get_fluxes", self.nrea, n, Tgas, par, Phi,
                           opts)

    def coolant_cooling(self, n, Tgas, par=None, Phi=None, opts=None):
        """each line coolant's cooling [erg cm^-3 s^-1], [..., ncoolant]
        (labels: coolants; legacy krome_cooling<metal>): the Z and
        coolant-file cooling terms are their sums"""
        nc = len(self.coolants)
        out = self._query("krome_get_coolant_cooling", max(nc, 1), n, Tgas,
                          par, Phi, opts)
        return out[..., :nc]

    def coolant_populations(self, n, Tgas, coolant, par=None, Phi=None,
                            opts=None):
        """the level populations [cm^-3] of a coolant (name or index),
        [..., nlev] in the file's level numbering (legacy
        krome_popcool_dump; levels: coolant_levels)"""
        i = coolant if isinstance(coolant, (int, np.integer)) \
            else self.coolants.index(coolant)
        single, n, ncell = self._cells(n)
        T = self._per_cell(Tgas, ncell)
        opts = self._options(opts)
        P, stride = self._params(par, ncell)
        F = self._phi(Phi, ncell, opts)
        out = np.zeros((ncell, len(self.coolant_levels[i][0])))
        rc = self._lib.krome_get_coolant_populations(
            ncell, _dptr(n), _dptr(T), _dptr(P), stride, _dptr(F), opts._h,
            int(i), _dptr(out))
        if rc != 0:
            raise KromeError("krome_get_coolant_populations returned %d"
                             % rc)
        return out[0] if single else out

    def photo_heating(self, n, Tgas, opts, par=None, Phi=None):
        """photoheating per photo reaction [erg cm^-3 s^-1], [..., nphoto]
        (labels: photo_reactions), summing to the PHOTO heating term
        (legacy krome_get_photoBin_heats gave it per absorber)"""
        npho = len(self.photo_reactions)
        out = self._query("krome_get_photo_heating", max(npho, 1), n, Tgas,
                          par, Phi, opts)
        return out[..., :npho]

    def ratio_flux_draine(self, Phi, opts):
        """the field over the Draine ISRF (legacy
        krome_get_ratioFluxDraine); Phi[nbins] or Phi[ncell, nbins]"""
        F = np.ascontiguousarray(Phi, dtype=float)
        single = F.ndim == 1
        F = F.reshape(1, -1) if single else F
        out = np.zeros(F.shape[0])
        rc = self._lib.krome_get_ratio_flux_draine(F.shape[0], _dptr(F),
                                                   opts._h, _dptr(out))
        if rc != 0:
            raise KromeError("krome_get_ratio_flux_draine returned %d" % rc)
        return out[0] if single else out

    def free_fall_time_rho(self, rho):
        """free-fall time [s] from the mass density rho [g cm^-3]"""
        r = np.ascontiguousarray(np.atleast_1d(rho), dtype=float)
        out = np.zeros_like(r)
        self._lib.krome_get_free_fall_time_rho(len(r), _dptr(r), _dptr(out))
        return out[0] if np.ndim(rho) == 0 else out

    def num2col(self, num, n, Tgas, par=None):
        """the column density [cm^-2] of the number density num [cm^-3] in
        the cell (radiation.num2col; legacy krome_num2col)"""
        single, n, ncell = self._cells(n)
        T = self._per_cell(Tgas, ncell)
        x = self._per_cell(num, ncell)
        P, stride = self._params(par, ncell)
        out = np.zeros(ncell)
        rc = self._lib.krome_num2col(ncell, _dptr(x), _dptr(n), _dptr(T),
                                     _dptr(P), stride, _dptr(out))
        if rc != 0:
            raise KromeError("krome_num2col returned %d" % rc)
        return out[0] if single else out

    def Tdust(self, n, Tgas, par=None, Phi=None, opts=None):
        """the dust temperature the rates see [K], one per cell (legacy
        krome_get_table_Tdust): the dust.tables value, else the Tdust
        parameter (the binned-dust mean with dust bins on opts)"""
        out = self._query("krome_get_Tdust", 1, n, Tgas, par, Phi, opts)
        return out[..., 0] if np.ndim(out) else out

    def thermal(self, n, Tgas, par=None, Phi=None, opts=None):
        """thermal balance [erg cm^-3 s^-1]: a Thermal with the totals
        heat and cool, the terms heating / cooling ({name: value}; FLOOR is
        the negative floor correction) and dTdt [K/s]"""
        single, n, ncell = self._cells(n)
        T = self._per_cell(Tgas, ncell)
        opts = self._options(opts)
        P, stride = self._params(par, ncell)
        F = self._phi(Phi, ncell, opts)
        nh, nc = len(self.heating_terms), len(self.cooling_terms)
        heat, cool, dT = (np.zeros(ncell) for _ in range(3))
        heats = np.zeros((ncell, max(nh, 1)))
        cools = np.zeros((ncell, max(nc, 1)))
        rc = self._lib.krome_get_thermal(
            ncell, _dptr(n), _dptr(T), _dptr(P), stride, _dptr(F), opts._h,
            _dptr(heat), _dptr(cool), _dptr(heats), _dptr(cools), _dptr(dT))
        if rc != 0:
            raise KromeError("krome_get_thermal returned %d" % rc)
        pick = (lambda a: a[0]) if single else (lambda a: a)
        return Thermal(pick(heat), pick(cool),
                       {s: pick(heats[:, i])
                        for i, s in enumerate(self.heating_terms)},
                       {s: pick(cools[:, i])
                        for i, s in enumerate(self.cooling_terms)},
                       pick(dT))

    def gamma(self, n, Tgas, par=None):
        """adiabatic index (codegen thermal.gamma)"""
        single, n, ncell = self._cells(n)
        T = self._per_cell(Tgas, ncell)
        P, stride = self._params(par, ncell)
        out = np.zeros(ncell)
        rc = self._lib.krome_get_gamma(ncell, _dptr(n), _dptr(T), _dptr(P),
                                       stride, _dptr(out))
        if rc != 0:
            raise KromeError("krome_get_gamma returned %d" % rc)
        return out[0] if single else out

    def qeff(self, n, Tgas, par=None):
        """energy released per reaction event [..., nrea] [MeV] (; qeff=,
        0 elsewhere); times fluxes() it is the nuclear energy rate"""
        single, n, ncell = self._cells(n)
        T = self._per_cell(Tgas, ncell)
        P, stride = self._params(par, ncell)
        out = np.zeros((ncell, max(self.nrea, 1)))
        rc = self._lib.krome_get_qeff(ncell, _dptr(n), _dptr(T), _dptr(P),
                                      stride, _dptr(out))
        if rc != 0:
            raise KromeError("krome_get_qeff returned %d" % rc)
        out = out[:, :self.nrea]
        return out[0] if single else out

    def _cellwise(self, fname, n):
        single, n, ncell = self._cells(n)
        out = np.zeros(ncell)
        getattr(self._lib, fname)(ncell, _dptr(n), _dptr(out))
        return out[0] if single else out

    def mu(self, n):
        """mean molecular weight rho / (ntot m_p)"""
        return self._cellwise("krome_get_mu", n)

    def rho(self, n):
        """mass density [g cm^-3]"""
        return self._cellwise("krome_get_rho", n)

    def ntot(self, n):
        return self._cellwise("krome_get_ntot", n)

    def hnuclei(self, n):
        """H nuclei density [cm^-3]"""
        return self._cellwise("krome_get_Hnuclei", n)

    def charge(self, n):
        """net charge density [e cm^-3]"""
        return self._cellwise("krome_get_charge", n)

    def electrons(self, n):
        """electron density that makes the cell neutral (legacy
        get_electrons)"""
        return self._cellwise("krome_get_electrons", n)

    def free_fall_time(self, n):
        return self._cellwise("krome_get_free_fall_time", n)

    def jeans_length(self, n, Tgas):
        single, n, ncell = self._cells(n)
        T = self._per_cell(Tgas, ncell)
        out = np.zeros(ncell)
        self._lib.krome_get_jeans_length(ncell, _dptr(n), _dptr(T),
                                         _dptr(out))
        return out[0] if single else out

    def element_totals(self, n):
        """element totals [atoms cm^-3], [..., nelem] (labels: elements)"""
        single, n, ncell = self._cells(n)
        out = np.zeros((ncell, max(len(self.elements), 1)))
        self._lib.krome_get_element_totals(ncell, _dptr(n), _dptr(out))
        out = out[:, :len(self.elements)]
        return out[0] if single else out

    def element_totals_x(self, x):
        """element totals per gram of gas [atoms g^-1] from the mass
        fractions, [..., nelem] (legacy krome_conserveLinGetRef_x)"""
        single, x, ncell = self._cells(x)
        out = np.zeros((ncell, max(len(self.elements), 1)))
        self._lib.krome_get_element_totals_x(ncell, _dptr(x), _dptr(out))
        out = out[:, :len(self.elements)]
        return out[0] if single else out

    def conserve(self, n, ref, par=None, flags=None, check=True):
        """consistent advection (legacy krome_conserveLin_x; manual 6.7):
        the cells n ([nspec] or [ncell, nspec], cm^-3) moved back to the
        reference element totals ref ([nelem]: one row for all, or [ncell,
        nelem]; element_totals of a reference cell) and to zero net charge
        (negatives clipped, the elements by the minimal relative correction
        of the conserve option, then the electrons from neutrality). par:
        the parameters as for solve (for the freeze directives). flags: an
        int array [ncell] receiving each cell's flag (0,
        FLAG_CONSERVE_UNREACHED -112, the warning FLAG_CONSERVE_LARGE 2).
        Returns the corrected densities (the input is not modified); raises
        KromeError when a cell cannot reach its reference, unless
        check=False"""
        return self._conserve(n, ref, par, flags, check,
                              "krome_conserve_cells")

    def conserve_x(self, x, ref, par=None, flags=None, check=True):
        """conserve() for mass fractions, ref per gram (element_totals_x)"""
        return self._conserve(x, ref, par, flags, check,
                              "krome_conserve_cells_x")

    def _conserve(self, v, ref, par, flags, check, fname):
        single, v, ncell = self._cells(v)
        v = v.copy()
        nel = len(self.elements)
        R = np.ascontiguousarray(ref, dtype=float)
        if R.shape == (nel,):
            rs = 0
        elif R.shape == (ncell, nel):
            rs = nel
        else:
            raise KromeError("ref must have shape (%d,) or (%d, %d), got %s"
                             % (nel, ncell, nel, R.shape))
        if nel == 0:
            R = np.zeros(1)
        # the parameters only select the freeze directives: no REQUIRED check
        P, ps = None, 0
        if isinstance(par, dict):
            P, ps = self.param_rows(ncell, **par), self.npar
        elif par is not None:
            P = np.ascontiguousarray(par, dtype=float)
            if P.shape == (ncell, self.npar):
                ps = self.npar
            elif P.shape != (self.npar,):
                raise KromeError("par must have shape (%d,) or (%d, %d), got "
                                 "%s" % (self.npar, ncell, self.npar,
                                         P.shape))
        fl = np.zeros(ncell, dtype=np.intc)
        r = getattr(self._lib, fname)(ncell, _dptr(v), _dptr(R), rs,
                                      None if P is None else _dptr(P), ps,
                                      fl.ctypes.data_as(_I))
        if r == -1:
            raise KromeError("%s: bad arguments" % fname)
        if flags is not None:
            flags[...] = fl.reshape(np.shape(flags))
        if check and r != 0:
            bad = np.flatnonzero(fl < 0)
            raise KromeError("%s: %s (cells %s); pass check=False to get "
                             "the result anyway" % (fname,
                                                    self.flag_message(r),
                                                    bad.tolist()))
        return v[0] if single else v

    def n2x(self, n):
        """mass fractions from number densities"""
        single, n, ncell = self._cells(n)
        out = np.zeros_like(n)
        self._lib.krome_n2x(ncell, _dptr(n), _dptr(out))
        return out[0] if single else out

    def x2n(self, x, rho):
        """number densities from mass fractions and rho [g cm^-3]"""
        single, x, ncell = self._cells(x)
        r = self._per_cell(rho, ncell)
        out = np.zeros_like(x)
        self._lib.krome_x2n(ncell, _dptr(x), _dptr(r), _dptr(out))
        return out[0] if single else out

    # -- legacy conveniences (krome_user) -----------------------------------

    def names_header(self):
        """the species names separated by spaces (legacy
        krome_get_names_header)"""
        return " ".join(self.species)

    def best_fluxes(self, n, Tgas, nbest=5, par=None, Phi=None, opts=None):
        """the nbest largest reaction fluxes of one cell: a list of
        (1-based reaction index, equation, flux [cm^-3 s^-1])"""
        fl = self.fluxes(n, Tgas, par, Phi, opts)
        order = sorted(range(len(fl)), key=lambda r: -fl[r])[:nbest]
        return [(r + 1, self.reactions[r], float(fl[r])) for r in order]

    def print_best_flux(self, n, Tgas, nbest=5, par=None, Phi=None,
                        opts=None, file=None):
        """print best_fluxes as legacy krome_print_best_flux"""
        out = file if file is not None else sys.stdout
        print(" ***************", file=out)
        for r, eq, f in self.best_fluxes(n, Tgas, nbest, par, Phi, opts):
            print("%8d %-50s %17.8E" % (r, eq, f), file=out)

    def explore_flux(self, file, n, Tgas, xvar, par=None, Phi=None,
                     opts=None):
        """every reaction's flux of one cell, absolute and relative to the
        largest and to the sum, one line per reaction then a blank line,
        labelled by xvar (legacy krome_explore_flux); file: an open text
        file or a path (appended)"""
        fl = self.fluxes(n, Tgas, par, Phi, opts)
        fmax, fsum = float(np.max(fl)), float(np.sum(fl))
        lines = ["%8d %17.8E %17.8E %17.8E %17.8E %17.8E   %s"
                 % (r + 1, xvar, Tgas, f, f / fmax, f / fsum,
                    self.reactions[r]) for r, f in enumerate(fl)]
        text = "\n".join(lines) + "\n\n"
        if hasattr(file, "write"):
            file.write(text)
        else:
            with open(file, "a") as f:
                f.write(text)

    @staticmethod
    def _emit(file, text):
        if file is None:
            sys.stdout.write(text)
        elif hasattr(file, "write"):
            file.write(text)
        else:
            with open(file, "a") as f:
                f.write(text)

    def print_coolant_populations(self, file, n, Tgas, par=None, Phi=None,
                                  opts=None):
        """every coolant's level populations of one cell, one row per
        level: name, level (1-based), Tgas, population, the coolant's
        total; a blank line at the end (legacy krome_popcool_dump); file:
        an open text file, a path (appended) or None (stdout)"""
        rows = []
        for i, name in enumerate(self.coolants):
            pop = self.coolant_populations(n, Tgas, i, par, Phi, opts)
            tot = float(np.sum(pop))
            rows += ["%8s%5d%17.8E%17.8E%17.8E" % (name, l + 1, Tgas, v, tot)
                     for l, v in enumerate(pop)]
        self._emit(file, "".join(r + "\n" for r in rows) + "\n")

    def dump_cooling(self, file, n, Tgas, par=None, Phi=None, opts=None):
        """one row: Tgas, the total cooling, every cooling term (legacy
        krome_dump_cooling)"""
        th = self.thermal(n, Tgas, par, Phi, opts)
        vals = [Tgas, th.cool] + [th.cooling[s] for s in self.cooling_terms]
        self._emit(file, "".join("%14.5E" % v for v in vals) + "\n")

    def dump_rates(self, file, Tmin, Tmax, npts, par=None, Phi=None,
                   opts=None):
        """the rate coefficients on npts log temperatures Tmin..Tmax at
        n = 1e-40: rows Tgas, reaction (1-based), k; a blank line after
        each temperature (legacy krome_dump_rates)"""
        n = np.full(self.nspec, 1e-40)
        a, b = np.log10(Tmin), np.log10(Tmax)
        out = []
        for i in range(npts):
            T = 10.0 ** (a + i * (b - a) / (npts - 1) if npts > 1 else a)
            k = self.rates(n, T, par, Phi, opts)
            out += ["%17.8E%8d%17.8E" % (T, r + 1, v) for r, v in
                    enumerate(k)]
            out.append("")
        self._emit(file, "\n".join(out) + "\n")

    def print_params(self, par=None, file=None):
        """every parameter, 'name: value' (legacy
        krome_print_phys_variables); par: None (defaults), {name: value}
        or a row par[npar]"""
        if par is None:
            row = self.params_default()
        elif isinstance(par, dict):
            row = self.param_rows(1, **par)[0]
        else:
            row = np.asarray(par, dtype=float).reshape(self.npar)
        self._emit(file, "".join("%s: %.16E\n" % (s, row[i])
                                 for i, s in enumerate(self.params)))

    def print_info(self, n, Tgas, file=None):
        """the species table: index, name, mass [g], charge, abundance,
        their sum, and Tgas (legacy krome_get_info)"""
        out = file if file is not None else sys.stdout
        n = np.asarray(n, dtype=float)
        print("%4s %-9s%11s%5s%11s" % ("#", "Name", "m (g)", "Chrg", "x"),
              file=out)
        for i, s in enumerate(self.species):
            print("%4d %-9s%11.3E%5d%11.3E" % (i + 1, s, self.masses[i],
                                               self.charges[i], n[i]),
                  file=out)
        print("%30s%11.3E" % (" sum", n.sum()), file=out)
        print("%14s%11.3E" % ("Tgas", Tgas), file=out)

    def scale_Z(self, n, Z):
        """a copy of n with the neutral atomic metals at the Asplund+09
        abundances times 10^Z (legacy krome_scale_Z); n[nspec] or
        n[ncell, nspec]"""
        single, n, ncell = self._cells(n)  # a copy
        self._lib.krome_scale_Z(ncell, _dptr(n), float(Z))
        return n[0] if single else n

    # -- radiation (legacy krome_set_photoBin_*, krome_get_opacity*) ----------
    # The photon field is Phi[..., nbins] (photon flux per bin [cm^-2
    # s^-1]); every spectrum is its bin integral 4 pi/h int J dlnE, with J
    # in eV s^-1 cm^-2 Hz^-1 sr^-1.

    def _bins(self, Elo, Ehi):
        lo = np.ascontiguousarray(Elo, dtype=float)
        hi = np.ascontiguousarray(Ehi, dtype=float)
        if lo.shape != hi.shape or lo.ndim != 1:
            raise KromeError("Elo and Ehi must be 1-D arrays of one shape")
        return lo, hi

    def _spec(self, fname, Elo, Ehi, *extra):
        lo, hi = self._bins(Elo, Ehi)
        out = np.zeros(lo.size)
        rc = getattr(self._lib, fname)(lo.size, _dptr(lo), _dptr(hi), *extra,
                                       _dptr(out))
        if rc != 0:
            raise KromeError("%s returned %d" % (fname, rc))
        return out

    def bins_lin(self, nbins, lower, upper):
        """(Elo, Ehi) of nbins linear bins in [lower, upper] eV"""
        lo, hi = np.zeros(nbins), np.zeros(nbins)
        if self._lib.krome_bins_lin(nbins, lower, upper, _dptr(lo),
                                    _dptr(hi)) != 0:
            raise KromeError("invalid binning")
        return lo, hi

    def bins_log(self, nbins, lower, upper):
        """(Elo, Ehi) of nbins logarithmic bins in [lower, upper] eV"""
        lo, hi = np.zeros(nbins), np.zeros(nbins)
        if self._lib.krome_bins_log(nbins, lower, upper, _dptr(lo),
                                    _dptr(hi)) != 0:
            raise KromeError("invalid binning")
        return lo, hi

    def J_to_phi(self, Elo, Ehi, J):
        """Phi of a J constant across each bin (legacy photoBinJ)"""
        j = np.ascontiguousarray(np.broadcast_to(np.asarray(J, dtype=float),
                                                 np.shape(Elo)))
        return self._spec("krome_J_to_phi", Elo, Ehi, _dptr(j))

    def phi_to_J(self, Elo, Ehi, Phi):
        """the J constant across each bin that carries Phi"""
        f = np.ascontiguousarray(Phi, dtype=float)
        return self._spec("krome_phi_to_J", Elo, Ehi, _dptr(f))

    def spectrum_blackbody(self, Elo, Ehi, Tbb):
        """undiluted blackbody at Tbb [K]"""
        return self._spec("krome_spectrum_blackbody", Elo, Ehi, float(Tbb))

    def blackbody_range(self, Tbb):
        """(lower, upper) [eV] where B_nu(Tbb) > 1e-6 (legacy BBlog_auto)"""
        lo, hi = ctypes.c_double(), ctypes.c_double()
        if self._lib.krome_blackbody_range(float(Tbb), ctypes.byref(lo),
                                           ctypes.byref(hi)) != 0:
            raise KromeError("no blackbody range for Tbb = %g" % Tbb)
        return lo.value, hi.value

    def spectrum_draine(self, Elo, Ehi):
        """Draine 1978 interstellar radiation field (5-13.6 eV)"""
        return self._spec("krome_spectrum_draine", Elo, Ehi)

    def spectrum_j21(self, Elo, Ehi, J21=1.0, alpha=1.0):
        """J = J21 1e-21 erg s^-1 cm^-2 Hz^-1 sr^-1 (13.6 eV / E)^alpha"""
        return self._spec("krome_spectrum_j21", Elo, Ehi, float(J21),
                          float(alpha))

    def spectrum_hm2012(self, Elo, Ehi, z):
        """Haardt & Madau 2012 UV background at redshift z (0-15.66)"""
        return self._spec("krome_spectrum_hm2012", Elo, Ehi, float(z))

    def spectrum_table(self, Elo, Ehi, E, J, loglog=False):
        """tabulated J(E) (E [eV]), linear or log-log interpolation, zero
        outside the table"""
        e = np.ascontiguousarray(E, dtype=float)
        j = np.ascontiguousarray(J, dtype=float)
        o = np.argsort(e)
        e, j = np.ascontiguousarray(e[o]), np.ascontiguousarray(j[o])
        return self._spec("krome_spectrum_table", Elo, Ehi, e.size, _dptr(e),
                          _dptr(j), int(bool(loglog)))

    def spectrum_file(self, Elo, Ehi, fname, loglog=False):
        """spectrum_table from a two-column file: E [eV], J (legacy
        krome_load_photoBin_file_2col; '#' comments)"""
        d = np.loadtxt(fname, comments="#", ndmin=2)
        return self.spectrum_table(Elo, Ehi, d[:, 0], d[:, 1], loglog)

    def dust_kappa(self, Elo, Ehi, E=None, kappa=None):
        """bin-mean dust opacity [cm^2 g^-1 of dust] from kappa(E) (E
        [eV]), default the embedded Draine R_V = 3.1 absorption opacity"""
        if E is None:
            return self._spec("krome_dust_kappa_draine", Elo, Ehi)
        lo, hi = self._bins(Elo, Ehi)
        e = np.ascontiguousarray(E, dtype=float)
        k = np.ascontiguousarray(kappa, dtype=float)
        o = np.argsort(e)
        e, k = np.ascontiguousarray(e[o]), np.ascontiguousarray(k[o])
        out = np.zeros(lo.size)
        if self._lib.krome_dust_kappa_bins(e.size, _dptr(e), _dptr(k),
                                           lo.size, _dptr(lo), _dptr(hi),
                                           _dptr(out)) != 0:
            raise KromeError("dust opacity table does not cover the bins")
        return out

    def dust_kappa_file(self, Elo, Ehi, fname, unit="eV"):
        """dust_kappa from a two-column file (legacy load_opacity_table):
        energy or wavelength (unit eV, keV, micron, cm, Angstrom), kappa
        [cm^2 g^-1]"""
        d = np.loadtxt(fname, comments="#", ndmin=2)
        x, k = d[:, 0], d[:, 1]
        hc = constants.planck_eV * constants.clight  # eV cm
        conv = {"ev": lambda v: v, "kev": lambda v: 1e3 * v,
                "micron": lambda v: hc / (v * 1e-4),
                "cm": lambda v: hc / v,
                "angstrom": lambda v: hc / (v * 1e-8)}
        try:
            E = conv[unit.lower()](x)
        except KeyError:
            raise KromeError("unknown unit %r (%s)" % (unit, ", ".join(conv))
                             ) from None
        return self.dust_kappa(Elo, Ehi, E, k)

    def find_G0_Av(self, Elo, Ehi, Phi, kappa, d2g, xH=None, n=None):
        """(G0, Av) of a field (legacy krome_find_G0_Av): kappa the
        per-bin dust opacity (dust_kappa), xH the H mass fraction, or
        computed from the densities n"""
        if xH is None:
            if n is None:
                raise KromeError("give xH or n")
            xH = float(self.hnuclei(n) * constants.p_mass / self.rho(n))
        lo, hi = self._bins(Elo, Ehi)
        f = np.ascontiguousarray(Phi, dtype=float)
        k = np.ascontiguousarray(kappa, dtype=float)
        G0, Av = ctypes.c_double(), ctypes.c_double()
        if self._lib.krome_find_G0_Av(lo.size, _dptr(lo), _dptr(hi), _dptr(f),
                                      _dptr(k), float(d2g), float(xH),
                                      ctypes.byref(G0), ctypes.byref(Av)) != 0:
            raise KromeError("find_G0_Av: invalid input (need >= 2 bins in "
                             "5-13.6 eV)")
        return G0.value, Av.value

    def photo_xsec(self, r, E):
        """cross section [cm^2] of photo reaction r (0..nphoto-1) at E [eV]"""
        E = np.asarray(E, dtype=float)
        out = np.array([self._lib.krome_photo_xsec(r, float(e))
                        for e in E.ravel()]).reshape(E.shape)
        return out if out.ndim else float(out)

    def habing(self, Phi, opts):
        """Habing flux G of Phi[..., nbins] (bins of opts)"""
        B = opts.nbins
        F = np.ascontiguousarray(np.asarray(Phi, dtype=float))
        single = F.ndim == 1
        F = F.reshape(-1, B)
        out = np.zeros(F.shape[0])
        if self._lib.krome_get_habing(F.shape[0], _dptr(F), opts._h,
                                      _dptr(out)) != 0:
            raise KromeError("habing: photon bins not configured")
        return out[0] if single else out

    def _radq(self, fname, n, Tgas, par, opts, *extra):
        single, n, ncell = self._cells(n)
        T = self._per_cell(Tgas, ncell)
        P, stride = self._params(par, ncell)
        out = np.zeros((ncell, opts.nbins))
        rc = getattr(self._lib, fname)(ncell, _dptr(n), _dptr(T), _dptr(P),
                                       stride, opts._h, *extra, _dptr(out))
        if rc != 0:
            raise KromeError("%s returned %d" % (fname, rc))
        return out[0] if single else out

    def opacity(self, n, Tgas, opts, par=None, d2g=None, kdust=None):
        """optical depth [..., nbins] with the num2col columns (legacy
        krome_get_opacity): photo and J_LW channels (shielded ones through
        the curve of growth), plus dust when d2g (scalar or per cell) is
        given: the dust bins of opts, else kdust [cm^2 per g of dust],
        default the Draine R_V = 3.1 opacity"""
        ncell = self._cells(n)[2]
        dg = None if d2g is None else self._per_cell(d2g, ncell)
        kd = None if kdust is None else np.ascontiguousarray(kdust,
                                                            dtype=float)
        return self._radq("krome_get_opacity", n, Tgas, par, opts,
                          _dptr(dg), _dptr(kd))

    def absorption(self, n, Tgas, opts, par=None, kdust=None):
        """photon removal coefficient [..., nbins] [cm^-1] (no length):
        the photochemistry with its shielding and 1/pdiss photons per
        event, the dust bins of opts, and rho kdust (kdust: dust opacity
        per unit GAS mass per bin, e.g. d2g * dust_kappa(...))"""
        kd = None if kdust is None else np.ascontiguousarray(kdust,
                                                            dtype=float)
        return self._radq("krome_get_absorption", n, Tgas, par, opts,
                          _dptr(kd))


class Thermal:
    """thermal() output: totals heat, cool [erg cm^-3 s^-1], the terms
    heating / cooling ({name: value}) and dTdt [K s^-1]."""

    def __init__(self, heat, cool, heating, cooling, dTdt):
        self.heat, self.cool = heat, cool
        self.heating, self.cooling, self.dTdt = heating, cooling, dTdt

    def __repr__(self):
        return "<krome.Thermal heat=%s cool=%s>" % (self.heat, self.cool)


class Options:
    """solver options (a krome_opts handle; docs/manual.md section 5).
    Keys: evolve_T, clip_negative, error_report (bool); jac ("auto" /
    "analytic" / "fd"), linsol ("auto" / "dense" / "lapack" / "sparse"), path
    ("auto" / "percell" / "device"); rtol, atol_n, atol_T, atol_u;
    max_steps, max_err_fails, max_conv_fails, max_restarts; conserve (0
    off, 1 charge, 2 elements + charge, 3 electrons for neutrality);
    device_cache, device_memory, min_device,
    min_interleave;
    eq_t0, eq_tmax, eq_factor, eq_rtol, eq_atol (the
    equilibrium solver); error_file (a string: the error-report file)."""

    def __init__(self, net, **values):
        self._net = net  # keeps the library loaded while the handle lives
        self._h = net._lib.krome_opts_new()
        for k, v in values.items():
            self[k] = v

    def __setitem__(self, key, value):
        if key == "error_file":  # the error-report file (a string)
            if self._net._lib.krome_opts_set_error_file(
                    self._h, str(value).encode()) != 0:
                raise KromeError("invalid error_file %r" % (value,))
            return
        v = value
        if key in _ENUMS and isinstance(value, str):
            try:
                v = _ENUMS[key].index(value.lower())
            except ValueError:
                raise KromeError("option %s: %r is not one of %s"
                                 % (key, value, _ENUMS[key])) from None
        if key not in _KEYS:
            raise KromeError("unknown option %r; options: %s"
                             % (key, ", ".join(_KEYS)))
        if self._net._lib.krome_opts_set(self._h, key.encode(),
                                         float(v)) != 0:
            raise KromeError("invalid value for option %s: %r"
                             % (key, value))

    def __getitem__(self, key):
        out = ctypes.c_double()
        if self._net._lib.krome_opts_get(self._h, key.encode(),
                                         ctypes.byref(out)) != 0:
            raise KeyError("unknown option %r" % key)
        v = out.value
        if key in _BOOLS:
            return v != 0.0
        if key in _ENUMS:
            return _ENUMS[key][int(v)]
        if key in _INTS:
            return int(v)
        return v

    def set(self, **values):
        for k, v in values.items():
            self[k] = v
        return self

    def set_photo(self, Elo, Ehi):
        """photon bins [eV edges]; Phi[..., nbins] is then required by
        solve()."""
        lo = np.ascontiguousarray(Elo, dtype=float)
        hi = np.ascontiguousarray(Ehi, dtype=float)
        if lo.shape != hi.shape or lo.ndim != 1:
            raise KromeError("Elo and Ehi must be 1-D arrays of one shape")
        if self._net._lib.krome_opts_set_photo(self._h, lo.size, _dptr(lo),
                                               _dptr(hi)) != 0:
            raise KromeError("invalid photon bins (need Ehi > Elo)")
        return self

    def set_photo_broadening(self, Tgas, vturb=0.0):
        """Doppler broadening of the tabulated cross sections, applied when
        the bins are built (legacy -useBroadening): Tgas [K], vturb
        [cm/s]; rebuilds the bins already set."""
        if self._net._lib.krome_opts_set_photo_broadening(
                self._h, float(Tgas), float(vturb)) != 0:
            raise KromeError("broadening needs Tgas >= 0 and vturb >= 0")
        return self

    @property
    def nbins(self):
        return self._net._lib.krome_opts_nbins(self._h)

    def lw_weights(self):
        """J_LW band weights [nbins]: J_LW = sum_b w_b Phi_b"""
        out = np.zeros(self.nbins)
        if self._net._lib.krome_opts_lw_weights(self._h, _dptr(out)) != 0:
            raise KromeError("photon bins not configured")
        return out

    def photo_sigma(self):
        """bin-averaged photo cross sections [nphoto, nbins] [cm^2]"""
        out = np.zeros((len(self._net.photo_reactions), self.nbins))
        if out.size and self._net._lib.krome_opts_photo_sigma(
                self._h, _dptr(out)) != 0:
            raise KromeError("photon bins not configured")
        return out

    def set_dustbins(self, nd=10, use_Si=True, use_C=True, alow=5e-7,
                     aup=2.5e-5, phi=-3.5):
        """binned dust population (MRN defaults as in the C++ DustBins)."""
        if self._net._lib.krome_opts_set_dustbins(
                self._h, nd, int(use_Si), int(use_C), alow, aup, phi) != 0:
            raise KromeError("invalid dust bins")
        return self

    @property
    def ndust(self):
        """number of dust bins (nd x grain types; 0 without bins)"""
        return self._net._lib.krome_opts_ndust(self._h)

    def set_dust_Td(self, Td):
        """the dust bins' temperatures [K] (ndust values): every cell's
        starting Td when solve() gets no Td, and the queries' Td"""
        t = np.ascontiguousarray(Td, dtype=float)
        if t.shape != (self.ndust,) or self._net._lib.krome_opts_set_dust_Td(
                self._h, _dptr(t)) != 0:
            raise KromeError("set_dust_Td needs %d positive values (dust "
                             "bins configured)" % self.ndust)
        return self

    @property
    def dust_Td(self):
        """the dust bins' temperatures [K]"""
        t = np.zeros(self.ndust)
        if self.ndust == 0 or self._net._lib.krome_opts_get_dust_Td(
                self._h, _dptr(t)) != 0:
            raise KromeError("no dust bins configured")
        return t

    def _set_species_tol(self, what, tol):
        net = self._net
        fn = getattr(net._lib, "krome_opts_set_%s_species" % what)
        if tol is None:
            fn(self._h, None)
            return self
        if isinstance(tol, dict):  # {species: value}, the rest the scalar
            t = -np.ones(net.nspec)
            for s, v in tol.items():
                t[net.idx(s)] = v
        else:
            t = np.ascontiguousarray(tol, dtype=float)
        if t.shape != (net.nspec,) or fn(self._h, _dptr(t)) != 0:
            raise KromeError("set_%s_species needs %d values (a dict "
                             "{species: value}), no NaN%s"
                             % (what, net.nspec,
                                ", no zero" if what == "atol" else ""))
        return self

    def set_atol_species(self, atol):
        """per-species absolute tolerances [cm^-3] (legacy -customATOL):
        nspec values or a dict {species: atol}, replacing atol_n for those
        species (a negative value, or a species not in the dict, takes
        atol_n); None: atol_n for every species"""
        return self._set_species_tol("atol", atol)

    def set_rtol_species(self, rtol):
        """per-species relative tolerances (legacy -customRTOL), as
        set_atol_species, replacing rtol"""
        return self._set_species_tol("rtol", rtol)

    @property
    def atol_species(self):
        """the absolute tolerance each species gets (nspec values)"""
        t = np.zeros(self._net.nspec)
        self._net._lib.krome_opts_get_atol_species(self._h, _dptr(t))
        return t

    @property
    def rtol_species(self):
        """the relative tolerance each species gets (nspec values)"""
        t = np.zeros(self._net.nspec)
        self._net._lib.krome_opts_get_rtol_species(self._h, _dptr(t))
        return t

    def set_dust_a(self, a):
        """the dust bins' grain radii [cm] (ndust values, dust.evolve):
        every cell's starting radii when solve() gets no adust"""
        t = np.ascontiguousarray(a, dtype=float)
        if t.shape != (self.ndust,) or self._net._lib.krome_opts_set_dust_a(
                self._h, _dptr(t)) != 0:
            raise KromeError("set_dust_a needs %d non-negative values (dust "
                             "bins configured)" % self.ndust)
        return self

    @property
    def dust_a(self):
        """the dust bins' grain radii [cm] (the reference sizes after
        set_dustbins)"""
        t = np.zeros(self.ndust)
        if self.ndust == 0 or self._net._lib.krome_opts_get_dust_a(
                self._h, _dptr(t)) != 0:
            raise KromeError("no dust bins configured")
        return t

    def dust_counts(self, rho_gas, d2g):
        """the grain counts per bin [cm^-3] of a gas density rho_gas
        [g cm^-3] and dust-to-gas parameter d2g"""
        t = np.zeros(self.ndust)
        if self.ndust == 0 or self._net._lib.krome_opts_get_dust_counts(
                self._h, float(rho_gas), float(d2g), _dptr(t)) != 0:
            raise KromeError("no dust bins configured")
        return t

    @property
    def dust_types(self):
        """the grain type of each bin: 1 carbon, 0 silicate"""
        t = np.zeros(self.ndust, dtype=np.intc)
        if self.ndust == 0 or self._net._lib.krome_opts_get_dust_types(
                self._h, t.ctypes.data_as(_I)) != 0:
            raise KromeError("no dust bins configured")
        return t

    def close(self):
        if getattr(self, "_h", None):
            self._net._lib.krome_opts_free(self._h)
            self._h = None

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class Result:
    """solve() output: n, Tgas, Phi (None without photon bins), Td and adust
    (None unless given), per-cell solver flags (0 = success) and ierr (0,
    or the last nonzero flag)."""

    def __init__(self, n, Tgas, Phi, flags, ierr, Td=None, adust=None):
        self.n, self.Tgas, self.Phi = n, Tgas, Phi
        self.flags, self.ierr = flags, ierr
        self.Td = Td  # per-cell dust temperatures, when given
        self.adust = adust  # per-cell grain radii (--evolve-dust), when given

    @property
    def ok(self):
        return self.ierr == 0

    def __repr__(self):
        return "<krome.Result ierr=%d, %d cell(s)>" % (
            self.ierr, np.atleast_1d(self.flags).size)
