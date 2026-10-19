# KROME-X from C++ (the header-only API). KROME = the project generated
# with ../options.toml (built there: make lib, for net_<name>_data.o)
KROME    ?= ../../../build_coupling
ifeq ($(origin CXX),default)
CXX = c++
endif
CXXFLAGS ?= -O2 -std=c++17

example: example.cpp $(KROME)/libkrome.a
	$(CXX) $(CXXFLAGS) -I$(KROME) example.cpp \
	  $(wildcard $(KROME)/net_*_data.o) -o $@

clean:
	rm -f example abundances.dat
