# KROME v2 coupling example (C, ENZO-like). KROME = the project generated
# with ../coupling.toml (built there: make lib)
KROME    ?= ../../../build_coupling
ifeq ($(origin CC),default)
CC = cc
endif
CFLAGS   ?= -O2
ifeq ($(shell uname),Darwin)
CXXLIB   ?= -lc++
else
CXXLIB   ?= -lstdc++
endif
LIBS =

coupling: coupling.c $(KROME)/libkrome.a
	$(CC) $(CFLAGS) -I$(KROME) coupling.c $(KROME)/libkrome.a $(LIBS) \
	  $(CXXLIB) -lm -o $@

clean:
	rm -f coupling result.dat
