# Data sources

KROME-X ships tables and databases compiled by other people. This file
lists where each one comes from and what to cite. The GPL-3.0 (`LICENSE`)
covers the KROME-X code; **the data files in this folder are not covered
by it**: they are distributed under the terms of their original
providers. Where a provider states terms, they are given below; the
others ask for citation. When you publish results that depend on a
dataset, please cite its original source as well as KROME-X.

Files marked *(legacy)* are inherited unchanged from legacy KROME
(Grassi et al. 2014). Files marked *verbatim* are unmodified copies of
the provider's files; the others are reformatted or derived from the
sources named.

## Thermochemistry

| file | content | source |
|---|---|---|
| `thermo_burcat.toml` | NASA 7-coefficient polynomials and formation enthalpies of the gas-phase species of the shipped networks (selected records, reformatted) | E. Goos, A. Burcat & B. Ruscic, *Extended Third Millennium Ideal Gas and Condensed Phase Thermochemical Database for Combustion with Updates from Active Thermochemical Tables*, THERM.DAT of 25 May 2021, https://respecth.elte.hu/burcat/; A. Burcat & B. Ruscic 2005, ANL-05/20 and TAE 960 |
| `thermoNIST.dat` *(legacy)* | Shomate fits, the fallback when a species has no Burcat record | NIST Chemistry WebBook; M. W. Chase, *NIST-JANAF Thermochemical Tables*, 4th ed. (1998) |
| `thermochemistry/*.gfe` *(legacy)* | Gibbs free energy against temperature of Al, Mg, Si and Ti oxide clusters and small molecules (nucleation networks) | not recorded in legacy KROME |

**Terms of the Burcat database**, as stated by its authors: the data are
free of charge on condition that proper quotation is given to the
database and its authors; it is strictly forbidden to include the
database, as is or parts of it, in any commercial software, firmware or
hardware without permission from the authors (Elke.Goos@dlr.de). These
terms apply to `thermo_burcat.toml` and to the coefficients KROME-X
writes from it into a generated project: KROME-X prints them when a
project uses the file (reverse reactions, the `DH` heating and cooling,
the formation enthalpies of `kromex-doc`) and writes them next to the
coefficients in the generated code.

## Reaction rates and atomic data

| file | content | source |
|---|---|---|
| `database/recombination_badnell.dat` | radiative and dielectronic recombination rate fits for whole isoelectronic sequences | N. R. Badnell 2006, ApJS 167, 334 (RR); N. R. Badnell et al. 2003, A&A 406, 1151 (DR); tables RR 2023-05-11 and DR 2023-05-12 from http://amdpp.phys.strath.ac.uk/tamoc/ |
| `database/sources/badnell_RR_clist_K.txt`, `badnell_DR_clist_K.txt` (*verbatim*) | the original coefficient tables behind the file above | © 2023 N. R. Badnell, Department of Physics, University of Strathclyde; free to use with citation of the papers above |
| `database/radiative_recombination.dat`, `database/sources/verner_radiative_recombination.txt` *(legacy)* | radiative recombination fits | D. A. Verner's atomic data, www.pa.uky.edu/~verner/atom.html |
| `database/photoioniziation.dat` *(legacy)* | photoionization cross-section fits | Verner et al. 1996 |
| `database/collisional_ionization.dat` *(legacy)* | collisional ionization fits | Voronov 1997, ADNDT 65, 1, via Verner's `cfit.f` |
| `database/KIDAchemistry.dat` *(legacy)* | gas-phase reactions in legacy KROME's database format | the KIDA database (Wakelam et al. 2012, 2015) |
| `database/hydroHe.dat` *(legacy)* | H and He reactions | not recorded in legacy KROME |
| `database/Alchemistry.dat` *(legacy)* | aluminium chemistry | Starik et al. 2014, Combust. Flame 161, 1659 |
| `database/surface_*.dat` *(legacy)* | grain-surface reactions, adsorption, desorption and chemisorption | Hollenbach & McKee 1979; Cazaux et al. 2004, 2010; Iqbal et al. 2010; Hocuk et al. 2014 (cited block by block in the files) |
| `surface_chemisorption_rates.dat` *(legacy)* | tabulated H chemisorption rates | not recorded in legacy KROME |
| `Ebare_ice.dat` *(legacy)* | binding energies on bare grains and water ice | arXiv:1408.5029, Table 2 |
| `Ebare_ice_Semenov.dat` *(legacy)* | binding energies for the Semenov et al. 2010 benchmark | the OSU network |
| `secion_xray_ratios.dat` | electron-impact ionization ratios of the metals, for the X-ray secondary ionization | Lennon et al. 1988, via Meijerink & Spaans 2005, eq. D.13; Quadri et al. 2025 |
| `asplund.dat` *(legacy)* | solar abundances | Asplund et al. 2009, Table 1 |

## Photochemistry and radiation

| file | content | source |
|---|---|---|
| `database/leiden_xsecs/*.dat` (*verbatim*) | photodissociation and photoionization cross sections | the Leiden database for photodissociation and photoionization of astrophysically relevant molecules, http://home.strw.leidenuniv.nl/~ewine/photo/ (A. Heays; Heays et al. 2017), with the references of each file |
| `database/swri_xsecs/*.dat` *(legacy, verbatim)* | photo cross sections | PHIDRATES, the SwRI photo rate database (Huebner & Mukherjee 2015, https://phidrates.space.swri.edu/), with the references of each file |
| `H2pdB.dat` *(legacy)* | H2 photodissociation line data for the photon bins | not recorded in legacy KROME |
| `HM2012.dat` *(legacy)* | the UV/X-ray background of quasars and galaxies against redshift (re-arranged) | Haardt & Madau 2012, ApJ 746, 125 (the file's header gives the volume as 764) |
| `shield_CO_unigrid.dat` | CO self- and mutual shielding | Visser, van Dishoeck & Black 2009, A&A 503, 323, Table 6, on the 0.2 dex grid of the paper's web page |
| `kabs_draine_Rv31.dat` *(legacy)* | dust absorption opacity per unit dust mass against wavelength | from the Lyman limit to 1 cm, the K_abs of B. T. Draine's `kext_albedo_WD_MW_3.1_60` (https://www.astro.princeton.edu/~draine/dust/dustmix.html, the pre-2003 tables: within 1% after resampling): the Milky Way R_V = 3.1 dust of Weingartner & Draine 2001 (b_C = 60 ppm) with the optical properties of Li & Draine 2001, Draine & Lee 1984 and Laor & Draine 1993; at shorter wavelengths (EUV and X-rays, to 1e-4 micron) the table differs from that file by up to tens of percent near the absorption edges and extends beyond it, and its source is not recorded |

## Cooling

| file | content | source |
|---|---|---|
| `coolZ.dat` | metal fine-structure coolants: levels, radiative and collisional rates | per block, in the file: LAMDA (Schöier et al. 2005; van der Tak et al. 2020); Barinovs et al. 2005; Lique et al. 2018; Klos et al. 2021; Yan & Babb 2023; CHIANTI 11.0.2 (Dere et al. 1997; Dufresne et al. 2024); NIST ASD (Kramida et al.); Tayal 2007, 2008; Zhang & Pradhan 1995; Hollenbach & McKee 1989; Maio et al. 2007 |
| `cooling/coolFeII_16.dat` | 16-level Fe+ coolant | CHIANTI 11.0.2; NIST ASD; Zhang & Pradhan 1995; Hollenbach & McKee 1989 |
| `cooling/coolO2.dat` | O2 coolant | LAMDA (O2: JPL and CDMS, 2010 update), converted by `lamda2krome` |
| `cooling/coolChianti.dat` *(legacy)* | Ne ion coolants | CHIANTI, converted by legacy KROME's `tools/chianti.py` |
| `cooling/sources/catom.dat`, `c+.dat`, `oatom@lique.dat` (*verbatim*) | the LAMDA files behind `coolZ.dat` | LAMDA (Schöier et al. 2005; van der Tak et al. 2020), with the references of each collision partner in the files |
| `cooling/sources/chianti/*` (*verbatim excerpts*) | CHIANTI 11.0.2 level, collision and radiative data of Fe II, O II, Si I and Si II (the kept levels, with the original reference footers) | CHIANTI (Dere et al. 1997; Dufresne et al. 2024) |
| `cooling/sources/nist_asd_*.txt` | NIST line lists (wavelength columns dropped) | A. Kramida, Yu. Ralchenko, J. Reader & NIST ASD Team, NIST Atomic Spectra Database, retrieved 2026-10-02 |
| `coolCO.dat`, `coolH2O.dat`, `coolHCN.dat`, `coolOH.dat` *(legacy)* | molecular cooling tables | Neufeld & Kaufman 1993; Omukai et al. 2010; data kindly provided by K. Omukai |
| `escape_H2.dat` *(legacy)* | H2 line escape fractions | courtesy of K. Omukai (2014) |
| `coolZ_CIE.tar.xz` *(legacy)* | metal cooling and heating in ionization equilibrium (re-arranged) | Shen et al. 2013, ApJ 765, 89 (CLOUDY, with and without the Haardt & Madau 2012 background); tables kindly provided by S. Shen |

## Dust

| file | content | source |
|---|---|---|
| `dust_optics/sources/callindex.out_*` (*verbatim*) | dielectric functions of astronomical silicate and graphite | B. T. Draine, https://www.astro.princeton.edu/~draine/dust/dust.diel.html; Draine 2003b, ApJ 598, 1026 |
| `optC.dat`, `optSi.dat` | absorption and scattering efficiencies and asymmetry parameter, computed with Mie theory from the dielectric functions above | Draine 2003b; Draine & Lee 1984; Mie series after Bohren & Huffman 1983 |
| `dust_tables/*` *(legacy)* | precomputed dust temperature, H2 formation and dust cooling tables (fields HM2012 and MC, and the 3D tables with Av) | Grassi et al. 2017 |

## Other

| file | content | source |
|---|---|---|
| `partition/partCO.dat`, `partH2even.dat`, `partH2odd.dat` *(legacy)* | rotational partition functions of CO and ortho/para H2 | A. Popovas, PhD thesis (2014); Popovas & Jørgensen 2016 |

## Data inside the code

Some literature coefficients are written directly in the source code,
with their references next to them: the RECFAST 1.5.2 atomic data and
fudge factors (`runtime/core/special_rates.hpp`: Seager, Sasselov &
Scott 1999, 2000; Wong, Moss & Scott 2008; Rubiño-Martín et al. 2010),
the Weingartner & Draine 2001 photoelectric heating fits
(`kromex/codegen.py`), and the dust optical tables and radiation fields
compiled from the files above (`runtime/core/data_qabs.hpp`,
`runtime/core/data_radiation.hpp`).
