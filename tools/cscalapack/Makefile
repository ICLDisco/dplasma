#export OMPI_CC=/opt/intel/bin/icc
#composer_xe_2013.3.163/bin/intel64/icc
#export OMPI_FC=/opt/intel/bin/ifort
#composer_xe_2013.3.163/bin/intel64/ifort
MPICC=mpicc
MPIF77=mpif77
CFLAGS+=
#-Wall -g
#-ftrapv
CC=$(MPICC)
LD=mpif77
#LD=mpif77 -nofor-main
MKLDIR=-L${MKLROOT}/lib/intel64
#MKLDIR=-L/nics/e/sw/local/kfs/opt/intel/composer_xe_2011_sp1.11.339/mkl/lib/intel64
#SMKL=$(MKLDIR)  -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_sequential -lmkl_intel_lp64 -lmkl_core
#MKL>=2016
SMKL=$(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
-lmkl_sequential -lmkl_gf_lp64 -lmkl_core -Wl,--end-group -Wl,-Bdynamic -ldl -lm
# Old MKL < 11.xx
#TMKL=-L/opt/mkl/lib/em64t -lmkl_scalapack_lp64 -lmkl_lapack -lmkl_blacs_openmpi_lp64 -lmkl -lmkl_gf_lp64 -lmkl_core -lguide -lpthread
# New MKL > 11.xx
#TMKL=$(MKLDIR)  -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lmkl_blacs_openmpi_lp64 -lm -fopenmp
#MKL>=2016
TMKL=$(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_core \
-lmkl_gnu_thread -lmkl_blacs_openmpi_lp64 -Wl,--end-group -Wl,-Bdynamic -ldl -lm -fopenmp
#TMKL=$(MKLDIR)  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_openmpi_lp64 -lm -openmp

TLIBS=$(TMKL)
SLIBS=$(SMKL)

ORIG_SMKL= $(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
-lmkl_sequential -lmkl_gf_lp64 -lmkl_core -Wl,--end-group -Wl,-Bdynamic -ldl
ORIG_TMKL= $(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_core \
-lmkl_gnu_thread -lmkl_blacs_openmpi_lp64 -Wl,--end-group -Wl,-Bdynamic -lm -fopenmp -ldl

ORIG_SLIBS=$(ORIG_SMKL)
ORIG_TLIBS=$(ORIG_TMKL)

WRAP_SLIBS=-L../../tests/wrappers -Wl,-rpath,../../tests/wrappers -lwrapper $(ORIG_SLIBS)
WRAP_TLIBS=-L../../tests/wrappers -Wl,-rpath,../../tests/wrappers -lwrapper $(ORIG_TLIBS)

TARGETS= \
pdpotrf_orig tpdpotrf_orig \
pdpotrf tpdpotrf \
pdgeqrf_orig tpdgeqrf_orig \
pdgeqrf tpdgeqrf \
pdgetrf_orig tpdgetrf_orig \
pdgetrf tpdgetrf \
pdgemm_orig tpdgemm_orig \
pdgemm tpdgemm \
pdtrmm_orig tpdtrmm_orig \
pdtrmm tpdtrmm

all: $(TARGETS)

common_orig.o: common.c
	$(CC) ${CFLAGS} -o $@ -c $<

p%_orig.o: p%.c
	$(CC) ${CFLAGS} -o $@ -c $<

p%_orig: p%_orig.o common_orig.o
	$(LD) $(LDFLAGS)	-o $@ $^ $(ORIG_SLIBS)

t%_orig: %_orig.o common_orig.o
	$(LD) $(LDFLAGS)	-o $@ $^ $(ORIG_TLIBS)


#WRAPPED VERSION
common.o: common.c
	$(CC) ${CFLAGS} -DDPLASMA_WRAPPER_ON -o $@ -c $<

p%.o: p%.c
	$(CC) ${CFLAGS} -DDPLASMA_WRAPPER_ON -c $<

p%: p%.o common.o
	$(LD) $(LDFLAGS)	-o $@ $^ $(WRAP_SLIBS)

t%: %.o common.o
	$(LD) $(LDFLAGS)	-o $@ $^ $(WRAP_TLIBS)

#t%: %.o common.o
#	$(LD) $(LDFLAGS)	-o $@ $^ $(TLIBS)

#p%: p%.o common.o
#	$(LD) $(LDFLAGS)	-o $@ $^ $(SLIBS)

clean:
	$(RM) *.o $(TARGETS)

