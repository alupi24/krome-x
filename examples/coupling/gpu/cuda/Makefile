# KROME-X from a GPU-resident CUDA C++ code. KROME = the project
# generated with ../../options.toml and built with Kokkos for the GPU
# (../README.md):
#   make -C $(KROME) lib KOKKOS=... CXX=.../nvcc_wrapper
KROME      ?= ../../../../build_coupling_gpu
KOKKOS     ?= $(error set KOKKOS to the Kokkos installation the project was built with)
KOKKOS_LIB ?= $(firstword $(wildcard $(KOKKOS)/lib64) $(KOKKOS)/lib)
NVCC       ?= nvcc
# the GPU: sm_80 = A100 (sm_90 = H100, ...)
NVCCFLAGS  ?= -O2 -std=c++17 -arch=sm_80
# GNU OpenMP when Kokkos runs OpenMP on the host (empty OMPLIB otherwise)
OMPLIB     ?= -lgomp
LIBS       ?= -L$(KOKKOS_LIB) -lkokkoscore -Xlinker -rpath=$(KOKKOS_LIB) \
              -ldl $(OMPLIB)

example: example.cu $(KROME)/libkrome.a
	$(NVCC) $(NVCCFLAGS) -I$(KROME) example.cu $(KROME)/libkrome.a \
	  $(LIBS) -o $@

clean:
	rm -f example abundances.dat
