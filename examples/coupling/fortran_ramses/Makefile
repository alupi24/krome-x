# KROME v2 coupling example (Fortran, RAMSES-like). KROME = the project
# generated with ../coupling.toml (built there: make lib fortran)
KROME    ?= ../../../build_coupling
ifeq ($(origin FC),default)
FC = gfortran
endif
FFLAGS   ?= -O2
ifeq ($(shell uname),Darwin)
CXXLIB   ?= -lc++
else
CXXLIB   ?= -lstdc++
endif
LIBS =

coupling: coupling.f90 $(KROME)/libkrome.a $(KROME)/krome.mod
	$(FC) $(FFLAGS) -I$(KROME) coupling.f90 $(KROME)/krome.o \
	  $(KROME)/krome_net.o $(KROME)/libkrome.a $(LIBS) $(CXXLIB) -o $@

clean:
	rm -f coupling result.dat
