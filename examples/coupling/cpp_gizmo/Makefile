# KROME v2 coupling example (C++, GIZMO-like, the header-only C++ API).
# KROME = the project generated with ../coupling.toml (built there: make
# lib, for net_<name>_data.o). OpenMP over the particles: make
# OPENMP=-fopenmp (Apple clang needs libomp)
KROME    ?= ../../../build_coupling
ifeq ($(origin CXX),default)
CXX = c++
endif
CXXFLAGS ?= -O2 -std=c++17
OPENMP   ?=
LIBS =

coupling: coupling.cpp $(KROME)/libkrome.a
	$(CXX) $(CXXFLAGS) $(OPENMP) -I$(KROME) \
	  coupling.cpp $(wildcard $(KROME)/net_*_data.o) $(LIBS) -o $@

clean:
	rm -f coupling result.dat
