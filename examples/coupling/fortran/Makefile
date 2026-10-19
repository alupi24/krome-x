# KROME-X from Fortran (the binding). KROME = the project generated with
# ../options.toml (built there: make lib fortran)
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

example: example.f90 $(KROME)/libkrome.a $(KROME)/krome.mod
	$(FC) $(FFLAGS) -I$(KROME) example.f90 $(KROME)/krome.o \
	  $(KROME)/krome_net.o $(KROME)/libkrome.a $(CXXLIB) -o $@

clean:
	rm -f example abundances.dat
