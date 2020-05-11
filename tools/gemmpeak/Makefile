CUDA_INSTALL_PATH ?= /opt/cuda
INCLUDES  := -I. -I$(CUDA_INSTALL_PATH)/include
LIB       := -L$(CUDA_INSTALL_PATH)/lib64 -lcublas -lcudart -lcuda

MKLDIR=-L${MKLROOT}/lib/intel64
#MKLDIR=-L/nics/e/sw/local/kfs/opt/intel/composer_xe_2011_sp1.11.339/mkl/lib/intel64
#SMKL=$(MKLDIR)  -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_sequential -lmkl_intel_lp64 -lmkl_core
#MKL>=2016
SMKL=$(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_sequential -lmkl_gf_lp64 -lmkl_core -Wl,--end-group -Wl,-Bdynamic -ldl -lm

#
CXX       := g++ 
NVCC      := nvcc -cubin
#NVCC      := nvcc -cubin
LINK      := g++ -fPIC

all : cu-gemmpeak mt-gemmpeak

mt-gemmpeak : mt-gemmpeak.o
	$(LINK) -o $@ $< -lpthread -lrt $(SMKL)

cu-gemmpeak : cu-gemmpeak.o sgemmN-70.cubin #sgemmN-13.cubin sgemmN-11.cubin
	$(LINK) $(LIB) -o $@ $<

sgemmN-70.cubin : sgemmN.cu
	$(NVCC) -arch sm_70 $< -o $@

sgemmN-13.cubin : sgemmN.cu
	$(NVCC) -m 32 -maxrregcount 32 -arch sm_13 $< -o $@

sgemmN-11.cubin : sgemmN.cu
	$(NVCC) -m 32 -maxrregcount 32 -arch sm_11 $< -o $@

clean :
	rm -f *.o *.cubin cu-gemmpeak mt-gemmpeak
