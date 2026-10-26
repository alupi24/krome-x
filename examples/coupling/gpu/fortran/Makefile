# KROME-X from a GPU-resident CUDA Fortran code (NVIDIA HPC SDK). KROME =
# the project generated with ../../options.toml and built with Kokkos for
# the GPU, its Fortran modules by nvfortran (../README.md):
#   make -C $(KROME) lib fortran KOKKOS=... CXX=.../nvcc_wrapper FC=nvfortran
KROME      ?= ../../../../build_coupling_gpu
KOKKOS     ?= $(error set KOKKOS to the Kokkos installation the project was built with)
KOKKOS_LIB ?= $(firstword $(wildcard $(KOKKOS)/lib64) $(KOKKOS)/lib)
NVFORTRAN  ?= nvfortran
# the GPU: cc80 = A100 (cc90 = H100, ...)
CUFFLAGS   ?= -O2 -cuda -gpu=cc80
# the C++ library, and GNU OpenMP when Kokkos runs OpenMP on the host
# (empty OMPLIB otherwise); load the GCC that built Kokkos last, so that
# its libstdc++ comes first
OMPLIB     ?= -lgomp
LIBS       ?= -L$(KOKKOS_LIB) -lkokkoscore -Wl,-rpath,$(KOKKOS_LIB) \
              -lstdc++ -ldl $(OMPLIB)

example: example.cuf $(KROME)/libkrome.a $(KROME)/krome.mod
	$(NVFORTRAN) $(CUFFLAGS) -I$(KROME) example.cuf $(KROME)/krome.o \
	  $(KROME)/krome_net.o $(KROME)/libkrome.a $(LIBS) -o $@

clean:
	rm -f example abundances.dat *.mod
