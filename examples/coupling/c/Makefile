# KROME-X from C. KROME = the project generated with ../options.toml
# (built there: make lib)
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

example: example.c $(KROME)/libkrome.a
	$(CC) $(CFLAGS) -I$(KROME) example.c $(KROME)/libkrome.a $(CXXLIB) \
	  -lm -o $@

clean:
	rm -f example abundances.dat
